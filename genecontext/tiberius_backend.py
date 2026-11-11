"""Pinned shared weights and an optional, isolated native-CPU runtime."""
from datetime import datetime, timezone
import json
import hashlib
import time
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile

from .download import transfer
from .sequence import read_fasta, read_gff, reverse_complement, sha256, translate, write_fasta


DEFAULT_PREDICTOR_CACHE = str(Path.home() / ".cache" / "gene-context" / "tiberius")
TIBERIUS_COMMIT = "b3ea48f263a075260dbf1e91a412e82ccdc96d3c"
WORKER = Path(__file__).with_name("tiberius_worker.py")
MODELS = {
    "vertebrates": dict(url="https://bioinf.uni-greifswald.de/bioinf/tiberius/models/vertebrates_weights.tar.gz",
        md5="62ddf2a876a0218756ebf0aea1f15963", sha256="cffe4571ba515adbfe8d68b9fec58d2682320011836c4f167466cba5917c0149",
        files={"model_config.json": "01cfb23a6c086b695ea9e32f0599393600ef6724dfe3968156243c9476e8657b",
               "model_layers.json": "e2a7035266f47803f2a196a63377b53f3e73f1f3d06ebfd6536259ef67dbf8fe",
               "weights.h5": "f80833ebc615c534e897a45a824dca1e148c4834b48720cada7001225c1a1545"}),
    "mammalia": dict(url="https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_nosm_weights_v2.tar.gz",
        md5="0bfd6ee0aee4647915f1761931e557cf", sha256="f57551eb5d6ee04a011d7b5c828e6283ded1e04d26135048185c9b27b1948053",
        files={"model_config.json": "c11be795ddc4c9fd9b87636e6ff483b52cbbeec917d2093a43e1c6aba97cdea4",
               "model_layers.json": "fcdb621baaa6364021e5fb3d239e66908d627299428e50ddd4e95854ad6885b2",
               "weights.h5": "974c18dc5e9e5ec053670f44b57e60290760d4fff97305dcc26a8f08775985bd"})}


def model_name(value):
    aliases = {"mammalia_nosofttmasking_v2": "mammalia", "mammalia": "mammalia",
               "vertebrata": "vertebrates", "vertebrates": "vertebrates"}
    return aliases.get(value.lower(), value)


def validate_model(path, expected=None):
    from .search_cache import fingerprint
    path = Path(path).expanduser().resolve()
    if not (path / "model_config.json").is_file():
        raise ValueError(f"No Tiberius model_config.json in {path}")
    weights = [p for p in (path / "weights.h5", path / "model.weights.h5") if p.is_file()]
    if len(weights) != 1:
        raise ValueError(f"Expected one Tiberius weight file in {path}")
    paths = [path / "model_config.json", *weights]
    if (path / "model_layers.json").is_file():paths.append(path / "model_layers.json")
    snapshots = {p.name:list(fingerprint(p)) for p in paths}
    config = json.loads((path / "model_config.json").read_text())
    if not isinstance(config, dict) or config.get("clamsa") or config.get("hmm"):
        raise ValueError("Use a Tiberius ab initio model without ClaMSA or an embedded HMM")
    if config.get("inp_size") == 6 or config.get("softmasking"):
        raise ValueError("This CPU backend requires an unmasked model; softmasked models are not silently reinterpreted")
    files = {p.name:sha256(p) for p in paths}
    if any(list(fingerprint(p)) != snapshots[p.name] for p in paths):
        raise ValueError('Tiberius weights changed during verification')
    if expected and files != expected:
        raise ValueError(f"Tiberius model fails pinned file checksums; preserved unchanged: {path}")
    return dict(path=str(path), files_sha256=files, files_fingerprint=snapshots)


def ensure_model(name, cache, offline=False):
    name = model_name(name)
    if name not in MODELS:
        return dict(name="custom", **validate_model(name))
    spec = MODELS[name]
    root = Path(cache).expanduser().resolve() / "models" / name
    model = root / "model"
    if model.exists():
        return dict(name=name, url=spec["url"], archive_sha256=spec["sha256"], **validate_model(model, spec["files"]))
    if offline:
        raise ValueError(f"Tiberius model {name} is not cached. Run gene-context setup-predictor first, or provide a local model directory.")
    root.mkdir(parents=True, exist_ok=True)
    lock = root / "prepare.lock"
    try:
        with lock.open("x") as handle:
            handle.write(datetime.now(timezone.utc).isoformat())
    except FileExistsError as error:
        raise RuntimeError(f"Tiberius model preparation is already running: {lock}") from error
    try:
        archive = root / "archive.tar.gz"
        print(f"[predictor] Preparing shared {name} model: {spec['url']}", file=sys.stderr, flush=True)
        transfer(spec["url"], archive, spec["md5"])
        if sha256(archive) != spec["sha256"]:
            raise ValueError(f"Tiberius archive differs from the pinned SHA-256: {archive}")
        with tempfile.TemporaryDirectory(prefix="extract-", dir=root) as temporary:
            staging = Path(temporary) / "model"; staging.mkdir()
            seen = set()
            with tarfile.open(archive, "r:gz") as packed:
                for member in packed:
                    if member.isdir():
                        continue
                    base = Path(member.name).name
                    member_path = PurePosixPath(member.name)
                    if (not member.isfile() or member_path.is_absolute() or ".." in member_path.parts or
                            base not in spec["files"] or base in seen or member.size > 512*1024**2):
                        raise ValueError("Unexpected Tiberius archive member")
                    seen.add(base)
                    with packed.extractfile(member) as inp, open(staging / base, "wb") as out:
                        shutil.copyfileobj(inp, out)
            validate_model(staging, spec["files"])
            staging.rename(model)
        info = dict(name=name, url=spec["url"], archive_sha256=spec["sha256"],
                    **validate_model(model, spec["files"]))
        (root / "receipt.json").write_text(json.dumps(info, indent=2) + "\n")
        return info
    finally:
        lock.unlink()


def runtime_python(args):
    if args.tiberius_python:
        value = shutil.which(str(args.tiberius_python)) or str(Path(args.tiberius_python).expanduser())
        # Resolving a venv executable symlink would select its base Python and
        # lose the isolated environment's installed packages.
        return str(Path(value).absolute()) if Path(value).is_file() else value
    return str(Path(args.predictor_cache).expanduser().resolve() / "runtime" / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python"))


def runtime_info(python):
    try:
        result = subprocess.run([python, str(WORKER), "--probe"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError("Tiberius CPU runtime is unavailable. Run gene-context setup-predictor, or set --tiberius-python to a prepared Python 3.12 environment.") from error
    if result.returncode:
        raise ValueError("Incomplete Tiberius CPU runtime; run gene-context setup-predictor. " + result.stderr[-1500:])
    info = json.loads(result.stdout)
    if info["packages"]["tiberius"] != "2.0.9":
        raise ValueError("This backend requires the tested Tiberius 2.0.9 runtime")
    return info


def setup_predictor(args):
    root = Path(args.predictor_cache).expanduser().resolve(); root.mkdir(parents=True, exist_ok=True)
    runtime = root / "runtime"
    python = str(runtime / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python"))
    if not args.models_only:
        lock = root / "setup.lock"
        try:
            handle = lock.open("x")
        except FileExistsError as error:
            raise RuntimeError(f"Predictor setup is already running: {lock}") from error
        try:
            with handle:
                handle.write(datetime.now(timezone.utc).isoformat())
            if (root / "runtime_receipt.json").exists():
                info = runtime_info(python)
            else:
                candidates = [args.python] if args.python else ([sys.executable] if (3, 12) <= sys.version_info[:2] <= (3, 13) else ["python3.12", "python3.13"])
                selected = next((shutil.which(p) for p in candidates if p and shutil.which(p)), None)
                if not selected:
                    raise ValueError("Install Python 3.12 and run setup-predictor --python /path/to/python3.12")
                version = subprocess.run([selected, "-c", "import sys; print(sys.version_info.minor if sys.version_info.major == 3 else 0)"], capture_output=True, text=True, check=True)
                if version.stdout.strip() not in {"12", "13"}:
                    raise ValueError("Tiberius CPU setup requires Python 3.12 or 3.13")
                log = root / "setup.log"
                commands = [[selected, "-m", "venv", str(runtime)],
                            [python, "-m", "pip", "install", "--disable-pip-version-check", "--no-input", "--no-cache-dir",
                             "-r", str(Path(__file__).with_name("tiberius_cpu_requirements.txt")),
                             "tiberius @ git+https://github.com/Gaius-Augustus/Tiberius.git@" + TIBERIUS_COMMIT]]
                print(f"[predictor] Installing isolated CPU runtime; log: {log}", file=sys.stderr, flush=True)
                with open(log, "w") as output:
                    for command in commands:
                        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
                        if result.returncode:
                            raise RuntimeError(f"CPU runtime setup failed; see {log}. Re-running retries the incomplete setup.")
                info = runtime_info(python)
                info.update(upstream_commit=TIBERIUS_COMMIT, requirements_sha256=sha256(Path(__file__).with_name("tiberius_cpu_requirements.txt")))
                (root / "runtime_receipt.json").write_text(json.dumps(info, indent=2) + "\n")
        finally:
            lock.unlink()
    models = [ensure_model(name, root) for name in args.models]
    receipt = dict(cache=str(root), python=None if args.models_only else python, models=models)
    if args.with_reference:
        from .context_reference import swissprot
        _, receipt["reference"] = swissprot(args.reference_cache)
    print(json.dumps(receipt, indent=2))
    return 0


def validate_prediction(folder, windows):
    """Check fresh and cached exports against the exact native window DNA."""
    output = folder / "predictions.gff3"
    proof = next((line.split("=", 1)[1] for line in (folder / "stdout.log").read_text().splitlines()
                  if line.startswith("GENE_CONTEXT_CPU=")), None)
    cpu = json.loads(proof) if proof else {}
    if (not cpu.get("visible_devices") or cpu.get("input_windows") != 1 or
            any("GPU" in device for device in cpu["visible_devices"])):
        raise ValueError("Tiberius did not confirm isolated CPU-only execution")
    models, _ = read_gff(output, source="tiberius")
    if any(m.contig not in windows for m in models):
        raise ValueError("Tiberius returned an unknown window identifier")
    # Empty predictions still require complete exports, not a truncated worker.
    def exported(suffix):
        path = Path(str(output)+suffix)
        return read_fasta(path) if path.stat().st_size else {}
    exported_cds = exported('.cds.fna')
    exported_proteins = exported('.proteins.faa')
    coding = {key.split("|")[1]: seq for key, seq in exported_cds.items()}
    proteins = {key.split("|")[1]: seq for key, seq in exported_proteins.items()}
    if (len(coding) != len(exported_cds) or len(proteins) != len(exported_proteins) or
            set(coding) != {m.id for m in models} or set(proteins) != set(coding)):
        raise ValueError("Tiberius GFF3/CDS/protein identifiers disagree")
    for model in models:
        if not model.cds or any(not 0 <= a < b <= len(windows[model.contig]) for a,b,_ in model.cds):
            raise ValueError("Tiberius CDS outside its input window")
        sequence = "".join(windows[model.contig][a:b] if model.strand == "+" else
            reverse_complement(windows[model.contig][a:b]) for a,b,_ in model.cds).upper()
        sequence = sequence[model.cds[0][2] or 0:]
        if sequence != coding[model.id] or translate(sequence).rstrip("*") != proteins[model.id].rstrip("*"):
            raise ValueError("Tiberius output fails native CDS/protein roundtrip")
    if not (folder / 'predictions.gff3.raw.gff3').is_file():
        raise ValueError("Tiberius pre-filter evidence is missing")
    return models, cpu


def predict_batches(args, root, rows, windows, runner, errors):
    from .artifact_cache import cached
    from .tiberius_session import WorkerSession
    from .search_cache import fingerprint
    python = runtime_python(args)
    runtime = runtime_info(python)
    mode = getattr(args, 'tiberius_worker_mode', 'isolated')
    provenance = dict(predictor="tiberius", cpu_only=True, runtime=runtime, batches=[],
        worker_mode=mode, window_isolation="independent upstream sequence-length clamping and fresh per-window HMM",
        method="DNA-only prediction before protein identification; no annotation/protein hints")
    implementation = {name:sha256(Path(__file__).with_name(name)) for name in
                      ['tiberius_backend.py','tiberius_worker.py','sequence.py','tiberius_session.py']}
    result = {}; weights = {}
    session = WorkerSession(python, WORKER, args, root, runner)
    try:
        for index, row in enumerate(rows, 1):
            name, locus = row['context_predictor_model'], row['locus_id']
            if row['window_bp'] < 1000:
                row['context_status'] = 'prediction_window_too_short'
                continue
            folder = root / 'context_prediction' / f'tiberius_{index:03d}'
            folder.mkdir(parents=True, exist_ok=True)
            batch = dict(model=name, loci=[locus], status='running')
            provenance['batches'].append(batch)
            started = time.perf_counter()
            commands_start = len(runner.commands)
            try:
                if name not in weights:
                    weights[name] = ensure_model(name, args.predictor_cache, args.offline)
                info = weights[name]; batch['weights'] = info
                def check_weights():
                    if any(list(fingerprint(Path(info['path'])/file)) != before
                           for file,before in info.get('files_fingerprint', {}).items()):
                        raise ValueError('Tiberius weights changed after verification')
                check_weights()
                own = {locus: windows[locus]}
                identity = dict(schema='gcf_tiberius_prediction_v1',
                    window_id=locus, dna_sha256=hashlib.sha256(windows[locus].encode()).hexdigest(),
                    weights_sha256=info['files_sha256'], runtime=runtime, implementation_sha256=implementation,
                    threads=args.threads, batch_size=args.tiberius_batch_size, worker_mode=mode,
                    requested_seq_len=400050, no_softmasking=True)

                def build(staging):
                    fasta = staging/'windows.fna'; output = staging/'predictions.gff3'
                    write_fasta(fasta, own.items())
                    print(f'[predictor] Tiberius CPU: {locus}, model {name}', file=sys.stderr, flush=True)
                    before = len(runner.commands)
                    try:
                        if mode == 'reuse':
                            session.predict(info['path'], fasta, output, staging/'stdout.log', staging/'stderr.log')
                        else:
                            runner.run([python, WORKER, '--threads', args.threads, '--batch-size', args.tiberius_batch_size,
                                        '--model', info['path'], '--fasta', fasta, '--out', output],
                                       staging/'stdout.log', staging/'stderr.log')
                        validate_prediction(staging, own)
                        check_weights()
                    except BaseException:
                        # Failed evidence is useful for diagnosis, but is never
                        # published as a reusable scientific prediction.
                        for file in staging.iterdir():
                            if file.is_file():shutil.copyfile(file, folder/file.name)
                        raise
                    # Keep actual invocation and timing with the cached evidence.
                    (staging/'execution.json').write_text(json.dumps(runner.commands[before:], indent=2)+'\n')

                stored, receipt = cached(Path(args.predictor_cache)/'predictions', identity, build)
                check_weights()
                for file in stored.iterdir():
                    if file.name != 'cache.json':
                        shutil.copyfile(file, folder/file.name)
                models, cpu = validate_prediction(folder, own)
                timing = next((line.split('=', 1)[1] for line in (folder/'stdout.log').read_text().splitlines()
                               if line.startswith('GENE_CONTEXT_TIMING=')), '{}')
                result[locus] = models
                batch.update(status='completed', coding_models=len(models), cpu_execution=cpu,
                    cache=receipt, execution_timing=json.loads(timing),
                    raw_gff3=str((folder/'predictions.gff3').relative_to(root)))
            except (OSError, ValueError, RuntimeError, IndexError) as error:
                batch.update(status='error', error=str(error))
                row['context_status'] = 'context_inference_error'
                errors.append(dict(stage='context_prediction', locus_id=locus, error=str(error)))
            for command in runner.commands[commands_start:]:
                for field in ['output', 'stderr']:
                    if command.get(field):command[field] = str(folder/Path(command[field]).name)
            batch['elapsed_seconds'] = time.perf_counter()-started
    finally:
        session.close()
    return result, provenance
