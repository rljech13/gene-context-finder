"""One lazy CPU process; requests still contain exactly one independent window."""
import json
import subprocess
import time


class WorkerSession:
    def __init__(self, python, worker, args, root, runner):
        self.argv = [str(python), str(worker), '--serve', '--threads', str(args.threads),
                     '--batch-size', str(args.tiberius_batch_size)]
        self.root = root; self.runner = runner; self.process = None; self.log = None

    def predict(self, model, fasta, output, stdout, stderr):
        if self.process is None:
            folder = self.root/'context_prediction'; folder.mkdir(parents=True, exist_ok=True)
            self.log = (folder/'session.stderr.log').open('a')
            self.process = subprocess.Popen(self.argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=self.log, text=True, bufsize=1)
        request = dict(model=str(model), fasta=str(fasta), out=str(output),
                       stdout=str(stdout), stderr=str(stderr))
        entry = dict(argv=self.argv, request=request, output=str(stdout), stderr=str(stderr), returncode=None)
        self.runner.commands.append(entry); started = time.perf_counter()
        try:
            self.process.stdin.write(json.dumps(request)+'\n'); self.process.stdin.flush()
            response = self.process.stdout.readline()
            if not response:
                raise RuntimeError('Tiberius worker stopped; see context_prediction/session.stderr.log')
            result = json.loads(response)
            entry['returncode'] = 0 if result.get('ok') else 1
            if not result.get('ok'):
                # A failed TensorFlow call may leave process state inconsistent.
                self.close()
                raise RuntimeError('Tiberius worker failed: '+result.get('error', 'unknown error'))
        except (BrokenPipeError, json.JSONDecodeError) as error:
            self.close()
            raise RuntimeError('Invalid Tiberius worker response') from error
        except RuntimeError:
            self.close()
            raise
        finally:
            entry['elapsed_seconds'] = time.perf_counter()-started

    def close(self):
        if self.process is not None:
            self.process.stdin.close()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.terminate()
                try:
                    self.process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.process.kill(); self.process.wait()
            self.process.stdout.close(); self.log.close()
            self.process = None
