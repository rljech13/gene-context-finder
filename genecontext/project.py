"""Study configuration and an append-only ledger above the native STOP stages.

Imported checkpoints are referenced in place. Only this project's own .gcf
directory is writable; neither attaching a result nor importing a decision runs
scientific tools. A decision is always bound to one exact checkpoint.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import socket
from uuid import uuid4

from .sequence import sha256

SCHEMA = 'gene_context_project_v1'
STAGES = ('stage1', 'stage2', 'stage3', 'stage4')
PATH_OPTIONS = {'datasets', 'source_datasets', 'neighbor_cache', 'protein_db',
                'predictor_cache', 'prediction_models', 'tiberius_python',
                'reference_panel', 'structures', 'structure_cache', 'iqtree', 'usalign'}
RESERVED = {'command', 'results', 'out', 'decisions', 'selection', 'draft',
            'validate_only', 'resume', 'gene', 'query', 'query_type', 'genomes',
            'query_metadata', 'help', 'version'}


def now():
    return datetime.now(timezone.utc).isoformat()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.'+path.name+'.'+uuid4().hex)
    try:
        with temporary.open('x') as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2)
            handle.write('\n'); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', value) or value in {'.', '..'}:
        raise ValueError('Use a unique target ID of 1–80 letters, digits, dots, underscores or hyphens')
    return value


def path_from(base, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Expected a nonempty file path')
    return (base/Path(value).expanduser()).resolve()


class Project:
    def __init__(self, path):
        path = Path(path).expanduser().resolve()
        self.path = path/'project.json' if path.is_dir() else path
        self.root = self.path.parent
        self.data = json.loads(self.path.read_text())
        if not isinstance(self.data, dict) or self.data.get('schema') != SCHEMA:
            raise ValueError('Expected '+SCHEMA)
        unknown = set(self.data)-{'schema', 'name', 'genomes', 'targets', 'options', 'execution'}
        if unknown:
            raise ValueError('Unknown project fields: '+', '.join(sorted(unknown)))
        if not isinstance(self.data.get('name'), str) or not self.data['name'].strip():
            raise ValueError('Project needs a name')
        self.genomes = path_from(self.root, self.data.get('genomes'))
        targets = self.data.get('targets')
        if not isinstance(targets, list) or not targets:
            raise ValueError('Project needs at least one target')
        self.targets = {}
        for target in targets:
            if not isinstance(target, dict):
                raise ValueError('Each target must be an object')
            unknown = set(target)-{'id', 'gene', 'query', 'query_type', 'query_metadata', 'options'}
            if unknown:
                raise ValueError('Unknown target fields: '+', '.join(sorted(unknown)))
            tid = identifier(target.get('id'))
            if tid in self.targets:
                raise ValueError('Duplicate target ID: '+tid)
            if not isinstance(target.get('gene'), str) or not target['gene'].strip():
                raise ValueError('Target needs a gene/search-group name: '+tid)
            if target.get('query_type') not in {'protein', 'cds', 'dna'}:
                raise ValueError('Target query_type must be protein, cds or dna: '+tid)
            path_from(self.root, target.get('query'))
            self.targets[tid] = target
            self.check_options(target.get('options', {}))
        self.check_options(self.data.get('options', {}))
        execution = self.data.get('execution', {})
        if not isinstance(execution, dict) or set(execution)-{'stage2'}:
            raise ValueError('Only execution.stage2 is configurable')
        profile = execution.get('stage2', {'mode': 'plan'})
        if not isinstance(profile, dict) or profile.get('mode') not in {'plan', 'local', 'queue'}:
            raise ValueError('execution.stage2.mode must be plan, local or queue')
        if set(profile)-{'mode', 'entrypoint', 'root', 'source', 'python'}:
            raise ValueError('Unknown stage2 execution field')
        if profile['mode'] == 'queue':
            for key in ['entrypoint', 'root', 'source']:
                path_from(self.root, profile.get(key))
        self.profile = profile
        self.store = self.root/'.gcf'
        if self.store.is_symlink():
            raise ValueError('The managed .gcf directory must not be a symlink')

    @staticmethod
    def check_options(options):
        if not isinstance(options, dict) or set(options)-set(STAGES):
            raise ValueError('Options must map stage1–stage4 to option objects')
        for stage, values in options.items():
            if not isinstance(values, dict) or any(not isinstance(k, str) or k.startswith('-') or '-' in k for k in values):
                raise ValueError('Options use CLI destinations with underscores: '+stage)
            if set(values) & RESERVED:
                raise ValueError('Project owns input/output/decision options: '+', '.join(sorted(set(values) & RESERVED)))

    def target(self, tid):
        if tid not in self.targets:
            raise ValueError('Unknown target: '+tid)
        return self.targets[tid]

    def options(self, tid, stage):
        return dict(self.data.get('options', {}).get(stage, {}),
                    **self.target(tid).get('options', {}).get(stage, {}))

    def load(self):
        path = self.store/'state.json'
        state = json.loads(path.read_text()) if path.exists() else dict(
            schema='gene_context_project_state_v1', runs=[], decisions=[], current={}, selected={})
        if state.get('schema') != 'gene_context_project_state_v1':
            raise ValueError('Unknown project ledger schema')
        return state

    def save(self, state):
        atomic(self.store/'state.json', state)

    @contextmanager
    def locked(self):
        # OS-released lock; no stale sentinel after a killed process. Compute
        # execution is supported on Linux/macOS; the portable viewer also has Windows builds.
        import fcntl
        self.store.mkdir(parents=True, exist_ok=True)
        with (self.store/'project.lock').open('a') as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError('Another project operation is running; status remains readable') from None
            try:
                yield self.load()
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def current(self, state, tid, stage):
        rid = state['current'].get(tid, {}).get(stage)
        return next((r for r in state['runs'] if r['id'] == rid), None)

    def decision(self, state, run):
        digest = state['selected'].get(run['id'])
        entry = next((d for d in state['decisions'] if d['run_id'] == run['id'] and d['sha256'] == digest), None)
        if entry and sha256(self.root/entry['object']) != digest:
            raise ValueError('Stored decision changed: '+entry['object'])
        return entry

    def snapshot_decision(self, source):
        raw = Path(source).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        extension = '.json' if raw.lstrip().startswith(b'{') else '.tsv'
        path = self.store/'decisions'/digest[:2]/(digest+extension)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            if path.read_bytes() != raw:
                raise ValueError('Decision object integrity failure')
        else:
            with path.open('xb') as handle:
                handle.write(raw); handle.flush(); os.fsync(handle.fileno())
        return path, digest


def checkpoint(stage, folder, full=True):
    """Verify scientific seals, never alter a checkpoint or its review state."""
    root = Path(folder).resolve()
    seal = json.loads((root/'seal.json').read_text())
    if not isinstance(seal, dict) or not seal or 'manifest.json' not in seal or ('evidence.json' if stage in STAGES[:2] else 'finalization.decisions.json') not in seal:
        raise ValueError('Incomplete checkpoint seal: '+stage)
    for name, digest in seal.items():
        path = (root/name).resolve()
        if root not in path.parents or not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest):
            raise ValueError('Unsafe checkpoint seal entry: '+name)
        if full and (not path.is_file() or sha256(path) != digest):
            raise ValueError('Checkpoint evidence changed: '+name)
    manifest = json.loads((root/'manifest.json').read_text())
    if stage in STAGES[:2]:
        data = json.loads((root/'evidence.json').read_text())
        schema = {'stage1': 'blast_loci_v1', 'stage2': 'stage2_review_v1'}[stage]
        if data.get('schema') != schema:
            raise ValueError('Wrong checkpoint stage: '+stage)
        required = {'evidence.json', 'manifest.json'}
        required |= {'queries.fasta', 'loci.tsv'} if stage == 'stage1' else {'windows.fna', 'proteins.faa', 'cds.fna'}
        if not required <= set(seal):
            raise ValueError('Missing scientific inputs in '+stage+' seal')
        review_id = data['review_id']
        gene = data.get('config', {}).get('gene') if stage == 'stage1' else None
        errors = data.get('errors', [])
        if stage == 'stage1':
            errors = errors + [g for g in data.get('genomes', []) if g.get('status') in {'error', 'failed'}]
    else:
        expected = {'stage3': 'stage3_comparison', 'stage4': 'stage4_prototype'}[stage]
        if manifest.get('stage') != expected:
            raise ValueError('Wrong checkpoint stage: '+stage)
        if stage == 'stage3':
            from .stage4 import checked_source
            if full:
                checked_source(root)
        data = json.loads((root/'finalization.decisions.json').read_text())
        review_id = data['review_id']; gene = None
        errors = manifest.get('errors', [])
    status = 'completed_with_errors' if errors or 'error' in str(manifest.get('status', '')) or 'error' in str(data.get('status', '')) else 'awaiting_review'
    if stage == 'stage4' and status != 'completed_with_errors':
        status = 'preview' if manifest.get('preview') else 'completed'
    return dict(path=str(root), seal_sha256=sha256(root/'seal.json'),
                manifest_sha256=sha256(root/'manifest.json'), review_id=review_id,
                evidence_sha256=sha256(root/'evidence.json') if stage in STAGES[:2] else None,
                source_review_id=manifest.get('source_review_id') or data.get('source_review_id') or (review_id if stage == 'stage4' else None),
                gene=gene, status=status)


def verify_record(run, full=False):
    info = checkpoint(run['stage'], run['output'], full)
    if info != run['checkpoint']:
        raise ValueError('Attached checkpoint identity changed: '+run['id'])
    return info


def validate_decision(project, tid, run, path):
    stage = run['stage']; source = Path(run['output'])
    verify_record(run, full=True)
    if run['checkpoint']['status'] == 'completed_with_errors':
        raise ValueError('Resolve checkpoint errors before recording a handoff')
    if stage == 'stage1':
        from .stage1 import verify, validate_decisions
        _, data = verify(source)
        if path.suffix.lower() == '.json':
            from .stage2 import selected_regions
            from .project_run import native_args
            args = native_args(project, tid, 'stage2', source=source, decision=path, out=project.store/'unused', validation=True)
            regions, _, _ = selected_regions(data, args)
            return dict(review_id=data['review_id'], ready=True, selected=len(regions))
        from .comparison import table
        result = validate_decisions(data, table(path))
        if not result['valid']:
            raise ValueError('\n'.join(result['errors']))
        return dict(review_id=data['review_id'], ready=bool(result['regions']), selected=len(result['regions']))
    if stage == 'stage2':
        from .stage2_choices import verify_checkpoint, validate_choices, issue_text
        from .stage3 import prepare_inputs
        _, data = verify_checkpoint(source); decision = json.loads(path.read_text())
        issues = validate_choices(data, decision)
        if issues:
            raise ValueError('\n'.join(issue_text(data, issue) for issue in issues))
        selected = prepare_inputs(data, decision, source)
        return dict(review_id=data['review_id'], ready=bool(selected), selected=len(selected))
    if stage == 'stage3':
        from .stage3_finalization import validate_choices
        evidence = json.loads((source/'finalization.evidence.json').read_text())
        state = validate_choices(json.loads(path.read_text()), evidence['spec'])
        return dict(review_id=state['review_id'], ready=all(c['reviewed'] for c in state['choices'].values()),
                    selected=len(state['choices']))
    raise ValueError('STOP4 is an export; review decisions belong to STOP1–3')


def attach(project, tid, stage, source, new_attempt=False):
    project.target(tid)
    with project.locked() as state:
        if project.store.is_relative_to(Path(source).resolve()):
            raise ValueError('Keep project metadata outside the imported scientific checkpoint')
        info = checkpoint(stage, source)
        current = project.current(state, tid, stage)
        if current and current.get('checkpoint') == info:
            return current
        if current and not new_attempt:
            raise ValueError('A current result already exists; --new-attempt retains it in history')
        target = project.target(tid)
        if stage == 'stage1' and info['gene'] != target['gene']:
            raise ValueError('STOP1 search-group name differs from target gene; describe the shared panel explicitly')
        parent = project.current(state, tid, STAGES[STAGES.index(stage)-1]) if stage != 'stage1' else None
        if parent:
            verify_record(parent, True)
            if info['source_review_id'] != parent['checkpoint']['review_id']:
                raise ValueError('Checkpoint belongs to another upstream review')
            if stage == 'stage2':
                data = json.loads((Path(source)/'evidence.json').read_text())
                if data.get('source_evidence_sha256') != parent['checkpoint']['evidence_sha256']:
                    raise ValueError('STOP2 source evidence differs from attached STOP1')
        run = dict(id=uuid4().hex, target=tid, stage=stage, output=info['path'],
                   imported=True, created=now(), status=info['status'], checkpoint=info,
                   source_run_id=parent['id'] if parent else None, decision_sha256=None)
        # Imported output has historical decisions of its own. Do not claim it
        # consumed the project's currently selected export without checking bytes.
        if parent:
            d = project.decision(state, parent)
            snapshot = Path(source)/('decisions.snapshot.json' if stage in {'stage3', 'stage4'} else 'selection.input.json')
            if d and snapshot.exists() and sha256(snapshot) == d['sha256']:
                run['decision_sha256'] = d['sha256']
        state['runs'].append(run); state['current'].setdefault(tid, {})[stage] = run['id']
        project.save(state)
        return run


def record_decision(project, tid, stage, source):
    with project.locked() as state:
        run = project.current(state, tid, stage)
        if not run:
            raise ValueError('Attach or run '+stage+' first')
        path = Path(source).expanduser().resolve()
        before = sha256(path)
        info = validate_decision(project, tid, run, path)
        if sha256(path) != before:
            raise ValueError('Decision changed during validation; export again')
        obj, digest = project.snapshot_decision(path)
        if digest != before:
            raise ValueError('Decision changed while saving')
        entry = dict(run_id=run['id'], stage=stage, target=tid, sha256=digest,
                     object=str(obj.relative_to(project.root)), original=str(path), created=now(), **info)
        if not any(d['run_id'] == run['id'] and d['sha256'] == digest for d in state['decisions']):
            state['decisions'].append(entry)
        state['selected'][run['id']] = digest
        project.save(state)
        return entry


def init(args):
    root = Path(args.out).expanduser().resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError('Initialize a new/empty project directory; attach existing results afterwards')
    genomes = Path(args.genomes).expanduser().resolve()
    if not genomes.is_file():
        raise ValueError('Genome manifest does not exist')
    target_file = Path(args.targets).expanduser().resolve()
    with target_file.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    targets = []
    for row in rows:
        if not any(row.values()):
            continue
        target = {k: row.get(k, '') for k in ['id', 'gene', 'query_type']}
        for key in ['query', 'query_metadata']:
            if row.get(key):
                target[key] = os.path.relpath(path_from(target_file.parent, row[key]), root)
        targets.append(target)
    data = dict(schema=SCHEMA, name=args.name, genomes=os.path.relpath(genomes, root), targets=targets,
                options={}, execution=dict(stage2=dict(mode=args.stage2_executor)))
    # Validate using a temporary config outside the user's results, then publish.
    root.mkdir(parents=True, exist_ok=True)
    config = root/'project.json'
    try:
        atomic(config, data); Project(config)
    except BaseException:
        config.unlink(missing_ok=True)
        raise
    return dict(project=str(config), targets=len(targets), status='configured', calculations_started=False)
