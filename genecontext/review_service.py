"""Loopback report delivery plus bounded asynchronous realignment of saved models.

The HTTP process only validates requests and calls the existing submit.py.
MAFFT and the portable engine always run in its shared resource slice.
"""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit, parse_qs, unquote
import argparse
import hashlib
import hmac
import json
import os
import re
import secrets
import sqlite3
import subprocess
import threading
import time

try:
    from .review_service_worker import sha, write_json, verified_result
except ImportError:
    from review_service_worker import sha, write_json, verified_result


class ReviewService:
    def __init__(self, config, submit=None):
        self.config = json.loads(Path(config).read_text())
        self.state = Path(self.config['state'])
        (self.state / 'jobs').mkdir(parents=True, exist_ok=True)
        self.lock = threading.RLock()
        secret = self.state / 'server.secret'
        if not secret.exists():
            with secret.open('x') as stream:
                stream.write(secrets.token_hex(32))
            secret.chmod(0o600)
        self.secret = secret.read_bytes()
        self.reports = {r['id']: r for r in self.config['reports']}
        self.paths = {p: r for r in self.config['reports'] for p in r['paths']}
        self.submit = submit or self.submit_job

    def session(self, params):
        path = params.get('report', [''])[0]
        report = self.paths.get(path)
        workspace = params.get('workspace', [''])[0]
        if not report or not re.fullmatch(r'[a-zA-Z0-9_-]{16,80}', workspace):
            raise ValueError('Unknown report or invalid workspace')
        if params.get('pool', [''])[0] != report['pool_sha256']:
            raise ValueError('This HTML has a different model pool. Refresh the report.')
        token = hmac.new(self.secret, (report['id'] + workspace + report['dataset_sha256']).encode(), hashlib.sha256).hexdigest()
        return report, workspace, token

    def view(self, job_id):
        if not re.fullmatch(r'[a-f0-9]{24}', job_id):
            raise ValueError('Invalid job ID')
        root = self.state / 'jobs' / job_id
        job = json.loads((root / 'job.json').read_text())
        state = json.loads((root / 'status.json').read_text()) if (root / 'status.json').exists() else {'status': 'queued', 'step': 'Waiting for submission'}
        receipt = root / 'receipt.json'
        if receipt.exists():
            run = json.loads(receipt.read_text())
            completed = Path(run['run_directory']) / 'completed.json'
            if completed.exists():
                done = json.loads(completed.read_text())
                if done.get('returncode') != 0:
                    state = {'status': 'failed', 'step': state['step'] if state.get('status') == 'failed' else 'Worker failed; inspect its receipt'}
                    write_json(root / 'status.json', state)
            elif state.get('status') == 'done':
                state = dict(state, status='running', step='Finalizing result receipt')
        elif state.get('status') == 'done':
            state = dict(state, status='running', step='Finalizing submission receipt')
        if state['status'] == 'done':
            if sha(root / 'result.json') != state['result_sha256']:
                raise ValueError('Saved result checksum failed')
            result = verified_result(state['verified_directory'])
            if result != json.loads((root / 'result.json').read_text()):
                raise ValueError('Saved result differs from sealed result')
            state = dict(state, result=result)
        return dict(state, job_id=job_id, review_id=job['review_id'], group=job['request']['group'], ids=job['request']['ids'], report=job['report'])

    def workspace(self, report, workspace, group):
        key = hashlib.sha256((report['id'] + workspace + group).encode()).hexdigest()
        return self.state / ('workspace-' + key + '.json')

    def status(self, report, workspace, group):
        response = {'ready': True, 'mode': 'hosted', 'review_id': report['review_id'], 'stage': report['stage']}
        latest = self.workspace(report, workspace, group)
        if latest.exists():
            response['job'] = self.view(json.loads(latest.read_text())['job_id'])
        return response

    def align(self, report, workspace, body):
        ids, group = body.get('ids'), body.get('group')
        if not isinstance(ids, list) or not 2 <= len(ids) <= 2000 or not all(isinstance(v, str) for v in ids) or len(set(ids)) != len(ids):
            raise ValueError('Select 2–2000 distinct source sequences')
        if sha(report['dataset']) != report['dataset_sha256']:
            raise ValueError('Pinned source dataset changed')
        dataset = json.loads(Path(report['dataset']).read_text())
        by_id = {s['id']: s for s in dataset['sequences']}
        if any(sid not in by_id or by_id[sid]['group'] != group for sid in ids):
            raise ValueError('Unknown model or different gene group')
        decisions = body.get('decisions')
        if decisions is not None and (not isinstance(decisions, dict) or decisions.get('review_id') != report['review_id']):
            raise ValueError('Decisions belong to a different checkpoint')
        request = {'ids': sorted(ids), 'group': group}
        if decisions is not None:
            request['decisions'] = decisions
        signature = json.dumps([report['id'], report['dataset_sha256'], request, self.config['engine_fingerprint']], sort_keys=True, separators=(',', ':')).encode()
        job_id = hashlib.sha256(signature).hexdigest()[:24]
        root = self.state / 'jobs' / job_id
        with self.lock:
            attempt = 0
            while root.exists() and self.view(root.name)['status'] == 'failed':
                attempt += 1
                job_id = hashlib.sha256(signature + str(attempt).encode()).hexdigest()[:24]
                root = self.state / 'jobs' / job_id
            active = 0
            for path in (self.state / 'jobs').glob('*/status.json'):
                if self.view(path.parent.name).get('status') in ['queued', 'running']:
                    active += 1
            if not root.exists():
                if active >= self.config.get('max_active', 2):
                    raise ValueError('Two alignments are running. Wait for one to finish, then press Recalculate.')
                root.mkdir()
                job = dict(report=report['id'], review_id=report['review_id'], dataset=report['dataset'],
                           dataset_sha256=report['dataset_sha256'], runner=self.config['runner'],
                           runner_sha256=self.config['runner_sha256'], request=request, created=time.time())
                write_json(root / 'job.json', job)
                write_json(root / 'status.json', {'status': 'queued', 'step': 'Submitting bounded MAFFT job'})
                threading.Thread(target=self.submit, args=(root,), daemon=True).start()
            write_json(self.workspace(report, workspace, group), {'job_id': job_id})
        return self.view(job_id)

    def submit_job(self, root):
        try:
            result = subprocess.run(['python3', self.config['submit'], '--name', 'review-mafft-' + root.name[:10], '--threads', '2',
                                     '--', 'python3', self.config['worker'], str(root)], check=True, capture_output=True, text=True, timeout=60)
            receipt = json.loads(result.stdout)
            if receipt.get('submission_returncode') != 0:
                raise ValueError('Submission rejected')
            write_json(root / 'receipt.json', receipt)
        except Exception as error:
            write_json(root / 'status.json', {'status': 'failed', 'step': 'Could not submit: ' + str(error)})


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, service, **kwargs):
        self.service = service
        super().__init__(*args, **kwargs)

    def api(self):
        parsed = urlsplit(self.path)
        if not parsed.path.startswith('/api/review/'):
            return False
        try:
            host = self.headers.get('Host', '')
            if host not in self.service.config['allowed_hosts']:
                raise ValueError('Invalid loopback host')
            origin = self.headers.get('Origin', '')
            if origin and origin != 'http://' + host:
                raise ValueError('Different origin')
            report, workspace, token = self.service.session(parse_qs(parsed.query))
            if self.command == 'POST' and not hmac.compare_digest(self.headers.get('X-Review-Token', ''), token):
                raise ValueError('Session token required')
            action = parsed.path.removeprefix('/api/review/')
            params = parse_qs(parsed.query)
            if self.command == 'GET' and action == 'status':
                value = self.service.status(report, workspace, params.get('group', [''])[0])
                value['token'] = token
            elif self.command == 'POST' and action in ['align', 'export']:
                size = int(self.headers.get('Content-Length', '0'))
                if not 0 < size <= 1_048_576:
                    raise ValueError('Request too large')
                body = json.loads(self.rfile.read(size))
                if not isinstance(body, dict):
                    raise ValueError('Expected a model-selection object')
                if action == 'align':
                    value = self.service.align(report, workspace, body)
                else:
                    path = self.service.workspace(report, workspace, body.get('group', ''))
                    job_id = json.loads(path.read_text())['job_id']
                    state = self.service.view(job_id)
                    if state['status'] != 'done' or state['result']['id'] != body.get('ID'):
                        raise ValueError('No verified result for this workspace')
                    archive = self.service.state / 'jobs' / job_id / 'alignment.zip'
                    if sha(archive) != state['archive_sha256']:
                        raise ValueError('Export checksum failed')
                    self.send_response(200)
                    self.send_header('Content-Type', 'application/zip')
                    self.send_header('Content-Length', str(archive.stat().st_size))
                    self.send_header('Cache-Control', 'no-store')
                    self.end_headers()
                    with archive.open('rb') as stream:
                        self.copyfile(stream, self.wfile)
                    return True
            else:
                self.send_error(404)
                return True
            raw = json.dumps(value).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(raw)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.end_headers()
            self.wfile.write(raw)
        except (ValueError, KeyError, OSError, TypeError) as error:
            raw = json.dumps({'error': str(error)}).encode()
            self.send_response(400)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(raw)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(raw)
        return True

    def do_POST(self):
        if not self.api():
            self.send_error(404)

    def do_GET(self):
        if self.api():
            return
        path = urlsplit(self.path).path
        match = re.fullmatch(r'/(?:original/)?(MYC|KLF4|KLF5|LIN28A|SOX3|SOX2|PRDM14|NANOS3)/loci/(locus_[a-f0-9]+)\.json', path)
        if match:
            root = Path(self.directory)
            dbpath = (root / 'original' if path.startswith('/original/') else root) / match[1] / 'loci.sqlite'
            if not dbpath.is_file():
                dbpath = root / 'original' / match[1] / 'loci.sqlite'
            if not dbpath.is_file():
                self.send_error(404)
                return
            with sqlite3.connect('file:' + str(dbpath) + '?mode=ro', uri=True) as database:
                row = database.execute('SELECT payload,sha256 FROM loci WHERE id=?', (match[2],)).fetchone()
            if row is None:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Encoding', 'gzip')
            self.send_header('Content-Length', str(len(row[0])))
            self.send_header('X-Payload-SHA256', row[1])
            self.end_headers()
            self.wfile.write(row[0])
            return
        super().do_GET()

    def send_head(self):
        path = Path(self.translate_path(urlsplit(self.path).path))
        compressed = Path(str(path) + '.gz')
        if path.suffix == '.html' and compressed.is_file() and compressed.stat().st_mtime >= path.stat().st_mtime:
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Encoding', 'gzip')
            self.send_header('Content-Length', str(compressed.stat().st_size))
            self.send_header('Cache-Control', 'no-cache')
            self.end_headers()
            return compressed.open('rb')
        return super().send_head()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--port', type=int, default=18794)
    args = parser.parse_args()
    service = ReviewService(args.config)
    handler = partial(Handler, service=service, directory=args.root)
    ThreadingHTTPServer(('127.0.0.1', args.port), handler).serve_forever()


if __name__ == '__main__':
    main()
