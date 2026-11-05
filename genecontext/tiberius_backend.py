"""Pinned shared weights and an optional, isolated native-CPU runtime."""
from collections import defaultdict
from datetime import datetime, timezone
import json
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
    path = Path(path).expanduser().resolve()
    if not (path / "model_config.json").is_file():
        raise ValueError(f"No Tiberius model_config.json in {path}")
    weights = [p for p in (path / "weights.h5", path / "model.weights.h5") if p.is_file()]
    if len(weights) != 1:
        raise ValueError(f"Expected one Tiberius weight file in {path}")
    config = json.loads((path / "model_config.json").read_text())
    if not isinstance(config, dict) or config.get("clamsa") or config.get("hmm"):
        raise ValueError("Use a Tiberius ab initio model without ClaMSA or an embedded HMM")
    if config.get("inp_size") == 6 or config.get("softmasking"):
        raise ValueError("This CPU backend requires an unmasked model; softmasked models are not silently reinterpreted")
    files = {p.name: sha256(p) for p in [path / "model_config.json", *weights]}
    if (path / "model_layers.json").is_file():
        files["model_layers.json"] = sha256(path / "model_layers.json")
    if expected and files != expected:
        raise ValueError(f"Tiberius model fails pinned file checksums; preserved unchanged: {path}")
    return dict(path=str(path), files_sha256=files)


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


def predict_batches(args, root, rows, windows, runner, errors):
    python = runtime_python(args)
    provenance = dict(predictor="tiberius", cpu_only=True, runtime=runtime_info(python), batches=[],
                      method="DNA-only prediction before protein identification; no annotation/protein hints")
    groups = defaultdict(list)
    for row in rows:
        groups[row["context_predictor_model"]].append(row)
    result = {}
    for index, (name, loci) in enumerate(sorted(groups.items()), 1):
        for row in loci:
            if row["window_bp"] < 1000:
                row["context_status"] = "prediction_window_too_short"
        loci = [row for row in loci if row["window_bp"] >= 1000]
        if not loci:
            continue
        folder = root / "context_prediction" / f"tiberius_{index:03d}"
        folder.mkdir(parents=True)
        batch = dict(model=name, loci=[r["locus_id"] for r in loci], status="running")
        provenance["batches"].append(batch)
        try:
            info = ensure_model(name, args.predictor_cache, args.offline)
            batch["weights"] = info
            fasta = folder / "windows.fna"
            write_fasta(fasta, [(r["locus_id"], windows[r["locus_id"]]) for r in loci])
            output = folder / "predictions.gff3"
            print(f"[predictor] Tiberius CPU: {len(loci)} windows, model {name}", file=sys.stderr, flush=True)
            runner.run([python, WORKER, "--threads", args.threads, "--batch-size", args.tiberius_batch_size,
                        "--model", info["path"], "--fasta", fasta, "--out", output], folder / "stdout.log", folder / "stderr.log")
            proof = next((line.split("=", 1)[1] for line in (folder / "stdout.log").read_text().splitlines()
                          if line.startswith("GENE_CONTEXT_CPU=")), None)
            batch["cpu_execution"] = json.loads(proof) if proof else {}
            if not batch["cpu_execution"].get("visible_devices") or any("GPU" in device for device in batch["cpu_execution"]["visible_devices"]):
                raise ValueError("Tiberius did not confirm CPU-only execution")
            models, _ = read_gff(output, source="tiberius")
            ids = {r["locus_id"] for r in loci}
            if any(m.contig not in ids for m in models):
                raise ValueError("Tiberius returned an unknown window identifier")
            exported_cds = read_fasta(str(output) + ".cds.fna") if models else {}
            exported_proteins = read_fasta(str(output) + ".proteins.faa") if models else {}
            coding = {key.split("|")[1]: seq for key, seq in exported_cds.items()}
            proteins = {key.split("|")[1]: seq for key, seq in exported_proteins.items()}
            if set(coding) != {m.id for m in models} or set(proteins) != set(coding):
                raise ValueError("Tiberius GFF3/CDS/protein identifiers disagree")
            for model in models:
                sequence = "".join(windows[model.contig][a:b] if model.strand == "+" else
                    reverse_complement(windows[model.contig][a:b])
                    for a, b, _ in model.cds).upper()
                phase = model.cds[0][2] or 0
                sequence = sequence[phase:]
                if sequence != coding[model.id] or translate(sequence).rstrip("*") != proteins[model.id].rstrip("*"):
                    raise ValueError("Tiberius output fails native CDS/protein roundtrip")
            for row in loci:
                result[row["locus_id"]] = [m for m in models if m.contig == row["locus_id"]]
            batch.update(status="completed", coding_models=len(models), raw_gff3=str(output.relative_to(root)))
        except (OSError, ValueError, RuntimeError, IndexError) as error:
            batch.update(status="error", error=str(error))
            for row in loci:
                row["context_status"] = "context_inference_error"
                errors.append(dict(stage="context_prediction", locus_id=row["locus_id"], error=str(error)))
    return result, provenance
