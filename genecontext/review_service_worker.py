"""Run the existing portable review engine inside a submit.py receipt."""
from pathlib import Path
import hashlib
import json
import os
import re
import selectors
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    pending = path.with_name(path.name + '.pending')
    pending.write_text(json.dumps(value, ensure_ascii=True, indent=2) + '\n')
    os.replace(pending, path)


def verified_result(root):
    root = Path(root)
    seals = json.loads((root / 'checksums.json').read_text())
    for name, digest in seals.items():
        if Path(name).name != name or sha(root / name) != digest:
            raise ValueError('Result seal failed: ' + name)
    return json.loads((root / 'result.json').read_text())


def run(job):
    job = Path(job).resolve()
    config = json.loads((job / 'job.json').read_text())
    process = None
    origin = token = None
    write_json(job / 'status.json', {'status': 'running', 'step': 'Starting MAFFT runner'})
    try:
        source = Path(config['dataset'])
        if sha(source) != config['dataset_sha256'] or sha(config['runner']) != config['runner_sha256']:
            raise ValueError('Pinned dataset or engine changed')
        work = job / 'engine'
        work.mkdir(exist_ok=True)
        data = work / 'review_data.json'
        if not data.exists():
            data.symlink_to(source)
        process = subprocess.Popen([config['runner'], '--root', str(work), '--no-browser'],
                                   stdout=subprocess.PIPE, stderr=(job / 'runner.stderr.log').open('ab'), text=True)
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        if not selector.select(20):
            raise RuntimeError('Runner did not start')
        line = process.stdout.readline()
        match = re.fullmatch(r'Open: (http://127\.0\.0\.1:\d+)/\?review_token=([a-f0-9]+)\s*', line)
        if not match:
            raise RuntimeError('Runner startup failed')
        origin, token = match.groups()

        def api(path, body):
            request = urllib.request.Request(origin + '/api/' + path, data=json.dumps(body).encode(),
                                             headers={'X-Review-Token': token, 'Content-Type': 'application/json'}, method='POST')
            return urllib.request.urlopen(request, timeout=3600)

        write_json(job / 'status.json', {'status': 'running', 'step': 'MAFFT protein alignment and codon backtranslation'})
        with api('align', config['request']) as response:
            result = json.load(response)
        root = work / 'review_results' / result['id']
        saved = verified_result(root)
        if saved != result or sorted(result['ids']) != sorted(config['request']['ids']):
            raise ValueError('Returned result does not match requested models')
        dataset = json.loads(source.read_text())
        sequences = {s['id']: s for s in dataset['sequences']}
        if result['review_id'] != dataset['review_id']:
            raise ValueError('Wrong source review')
        for sid in result['ids']:
            if result['protein'][sid].replace('-', '') != sequences[sid]['protein']:
                raise ValueError('Source protein changed: ' + sid)
        if sha(source) != config['dataset_sha256']:
            raise ValueError('Source changed during alignment')
        with api('export', {'ID': result['id']}) as response:
            archive = response.read()
        (job / 'alignment.zip.pending').write_bytes(archive)
        with zipfile.ZipFile(job / 'alignment.zip.pending') as archive:
            if archive.testzip() is not None:
                raise ValueError('Invalid result ZIP')
            for name, digest in json.loads(archive.read('checksums.json')).items():
                if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                    raise ValueError('Export checksum failed')
        os.replace(job / 'alignment.zip.pending', job / 'alignment.zip')
        write_json(job / 'result.json', result)
        write_json(job / 'status.json', {'status': 'done', 'step': 'Verified and saved', 'result_id': result['id'],
                                        'result_sha256': sha(job / 'result.json'), 'archive_sha256': sha(job / 'alignment.zip'),
                                        'verified_directory': str(root), 'finished': time.time()})
    except Exception as error:
        write_json(job / 'status.json', {'status': 'failed', 'step': str(error)})
        raise
    finally:
        if origin and token:
            try:
                request = urllib.request.Request(origin + '/api/shutdown', data=b'{}', method='POST',
                                                 headers={'X-Review-Token': token, 'Content-Type': 'application/json'})
                urllib.request.urlopen(request, timeout=5).close()
            except Exception:
                pass
        if process:
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=30)


if __name__ == '__main__':
    run(sys.argv[1])
