"""Explicit stage-2 restart identity, independent extraction checkpoints and lock."""
import json
from pathlib import Path
import shutil
from uuid import uuid4

from .early_synteny import token
from .sequence import sha256
from .stage1 import new_output
from .search_cache import fingerprint


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.'+uuid4().hex+'.tmp')
    try:
        temporary.write_text(json.dumps(value, indent=2)+'\n')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def prepare(args, plan, by_genome, protein_db, versions, predictor):
    root = Path(args.out).resolve()
    selected = sorted({r['genome_id'] for r in plan['regions']})
    inputs = {}; snapshots = {}
    for gid in selected:
        inputs[gid] = {}; snapshots[gid] = {}
        for field in ['fasta', 'gff3']:
            if by_genome[gid].get(field):
                path = by_genome[gid][field]; before = fingerprint(path)
                inputs[gid][field] = sha256(path)
                if fingerprint(path) != before:
                    raise ValueError('Stage 2 input changed while hashing: '+str(path))
                snapshots[gid][field] = before
    # Retain exact method identity. A different implementation starts a new run;
    # its prediction cache can still be reused when that method's key matches.
    implementation = {p.name:sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
    parameters = {k:v for k,v in vars(args).items() if k not in {'resume','out'} and not callable(v)}
    identity = dict(schema='gcf_stage2_resume_v1', plan=plan, parameters=parameters,
                    input_manifest_sha256=sha256(args.genomes), inputs=inputs,
                    protein_db=protein_db, tools=versions, predictor=predictor,
                    implementation_sha256=implementation)
    record_path = root/'.stage2'/'run.json'
    resumed = getattr(args, 'resume', False)
    if resumed:
        if not record_path.is_file():
            raise ValueError('--resume requires an output created by resumable stage2')
        record = json.loads(record_path.read_text())
        if record.get('identity') != identity:
            raise ValueError('Stage 2 resume inputs, decisions, parameters or implementation changed; use a new --out')
        if (root/'seal.json').exists():
            from .stage2_choices import verify_checkpoint
            _, data = verify_checkpoint(root)
            return root, record, inputs, snapshots, data
        # Never discard a prior attempt. Only these program-owned scratch
        # directories are moved; review choices remain in their original place.
        previous = root/'.stage2'/'attempts'/uuid4().hex
        for name in ['context_prediction','identification','target_support','neighbor_support']:
            if (root/name).exists():
                previous.mkdir(parents=True, exist_ok=True)
                shutil.move(str(root/name), str(previous/name))
    else:
        root = new_output(root)
        record = dict(identity=identity, review_id=str(uuid4()))
        atomic_json(record_path, record)
    return root, record, inputs, snapshots, None


def lock_path(output):
    root = Path(output).resolve()
    return root.parent/'_genecontext_cache'/'run_locks'/(token(str(root))+'.lock')


def protein_step(root, name, identity, runner, compute):
    from .artifact_cache import cached
    destination = root/name
    commands_start = len(runner.commands)
    def build(folder):
        try:
            value = compute(folder/'work')
            atomic_json(folder/'result.json', value)
            atomic_json(folder/'execution.json', runner.commands[commands_start:])
        except BaseException:
            if (folder/'work').exists():
                shutil.copytree(folder/'work', destination, dirs_exist_ok=True)
            raise
    try:
        stored, receipt = cached(root/'.stage2'/'protein_steps', dict(step=name, **identity), build)
        if (stored/'work').exists():
            shutil.copytree(stored/'work', destination)
        receipt['execution_receipt'] = str((stored/'execution.json').relative_to(root))
        return json.loads((stored/'result.json').read_text()), receipt
    finally:
        for command in runner.commands[commands_start:]:
            for field in ['output', 'stderr']:
                if command.get(field):
                    command[field] = str(destination/Path(command[field]).name)
