"""Explicit one-stage execution; no automatic crossing of a manual STOP."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
from uuid import uuid4

from .project import (Project, STAGES, PATH_OPTIONS, atomic, attach, checkpoint,
                      init, now, path_from, record_decision, signature, verify_record)
from .sequence import sha256


def options_argv(project, tid, stage):
    from .cli import parser
    command = next(a for a in parser()._actions if isinstance(a, argparse._SubParsersAction)).choices[stage]
    actions = {a.dest: a for a in command._actions}
    result = []
    for key, value in project.options(tid, stage).items():
        action = actions.get(key)
        if action is None or not action.option_strings:
            raise ValueError('Unknown '+stage+' option: '+key)
        flag = next((v for v in action.option_strings if v.startswith('--')), action.option_strings[0])
        if action.nargs == 0:
            if type(value) is not bool:
                raise ValueError('Boolean option expected: '+key)
            if value:
                result.append(flag)
        else:
            if value is None or isinstance(value, dict) or type(value) is bool:
                raise ValueError('Option needs a value: '+key)
            values = value if isinstance(value, list) else [value]
            if not values or any(not isinstance(v, (str, int, float)) or type(v) is bool for v in values):
                raise ValueError('Invalid option value: '+key)
            if key in PATH_OPTIONS:
                values = [str(path_from(project.root, v)) for v in values]
            if key == 'tiberius_model' and (str(value).startswith(('.', '/', '~')) or '/' in str(value)):
                values = [str(path_from(project.root, value))]
            result += [flag, *map(str, values)]
    return result


def native_args(project, tid, stage, *, out, source=None, decision=None, validation=False):
    from .cli import parser
    target = project.target(tid)
    argv = [stage, '--out', str(out)]
    if stage == 'stage1':
        argv += ['--gene', target['gene'], '--query', str(path_from(project.root, target['query'])),
                 '--query-type', target['query_type'], '--genomes', str(project.genomes)]
        if target.get('query_metadata'):
            argv += ['--query-metadata', str(path_from(project.root, target['query_metadata']))]
        if 'neighbor_cache' not in project.options(tid, stage):
            argv += ['--neighbor-cache', str(project.store/'cache'/'search')]
    else:
        argv += ['--results', str(source)]
        if decision:
            argv += ['--draft' if stage == 'stage2' and Path(decision).suffix == '.json' else '--decisions', str(decision)]
        if stage == 'stage2':
            argv += ['--genomes', str(project.genomes)]
            if 'protein_db' not in project.options(tid, stage):
                if not validation:
                    raise ValueError('Configure options.stage2.protein_db before stage2')
                argv += ['--protein-db', str(project.store/'not-used-during-decision-validation')]
    argv += options_argv(project, tid, stage)
    try:
        return parser().parse_args(argv)
    except SystemExit as error:
        raise ValueError('Invalid '+stage+' configuration') from error


def argv_for(args):
    """Round-trip the existing parser's destinations, keeping native defaults."""
    # Only actions defined by the native parser are serialized.
    from .cli import parser
    command = next(a for a in parser()._actions if isinstance(a, argparse._SubParsersAction)).choices[args.command]
    result = [args.command]
    for action in command._actions:
        if action.dest == 'help' or not action.option_strings:
            continue
        value = getattr(args, action.dest, None)
        if value is None:
            continue
        flag = next((f for f in action.option_strings if f.startswith('--')), action.option_strings[0])
        if action.nargs == 0:
            if value:
                result.append(flag)
        else:
            result += [flag, *map(str, value if isinstance(value, list) else [value])]
    return result


def input_files(project, tid, args):
    """Hash compact inputs and detect dataset changes without rehashing genomes.

    Native stages own exact source hashes. Stat identities here prevent a saved
    managed result being reused after an input has changed at the same path.
    Later comparisons need only their sealed upstream sequences, not the genomes.
    """
    target = project.target(tid)
    paths = {project.genomes} if args.command in {'stage1', 'stage2'} else set()
    if args.command == 'stage1':
        paths.add(path_from(project.root, target['query']))
    if args.command == 'stage1' and target.get('query_metadata'):
        paths.add(path_from(project.root, target['query_metadata']))
    for key in PATH_OPTIONS-{'neighbor_cache', 'predictor_cache', 'protein_db', 'tiberius_python', 'iqtree', 'usalign'}:
        value = getattr(args, key, None)
        if value:
            paths.add(Path(value))
    for path in paths:
        if not path.is_file():
            raise ValueError('Input file is missing: '+str(path))
    fingerprints = {}
    if args.command in {'stage1', 'stage2'}:
        from .early_synteny import datasets
        manifests = {project.genomes}
        for key in ['datasets', 'source_datasets']:
            if getattr(args, key, None):
                manifests.add(Path(getattr(args, key)))
        for manifest in manifests:
            for row in datasets(manifest):
                for key in ['fasta', 'gff3', 'protein', 'cds']:
                    if row.get(key):
                        path = Path(row[key]); stat = path.stat()
                        fingerprints[str(path)] = [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]
    return dict(files={str(p): sha256(p) for p in sorted(paths)}, datasets=fingerprints)


def status_rows(project, state=None):
    state = state or project.load(); rows = []
    for tid in project.targets:
        previous_ok = True
        for stage in STAGES:
            run = project.current(state, tid, stage)
            row = dict(target=tid, gene=project.target(tid)['gene'], stage=stage,
                       status='ready' if stage == 'stage1' else 'blocked', report=None, detail='')
            if run:
                row.update(run_id=run['id'], status=run['status'], report=str(Path(run['output'])/'report.html'))
                if run.get('checkpoint'):
                    try:
                        verify_record(run)
                        decision = project.decision(state, run)
                        row['decision_sha256'] = decision['sha256'] if decision else None
                        if stage != 'stage4' and decision:
                            row['status'] = 'reviewed' if decision['ready'] else 'review_incomplete' if decision['selected'] else 'closed_no_selection'
                        if stage != 'stage1' and run.get('source_run_id'):
                            parent = project.current(state, tid, STAGES[STAGES.index(stage)-1])
                            d = project.decision(state, parent) if parent else None
                            if not previous_ok or not parent or parent['id'] != run['source_run_id'] or (d and d['sha256'] != run.get('decision_sha256')):
                                row['status'] = 'stale'
                                row['detail'] = 'Upstream result or selected decisions differ; this historical run is retained.'
                    except (OSError, ValueError, KeyError) as error:
                        row.update(status='integrity_error', detail=str(error))
                elif row['status'] == 'running':
                    row['status'] = 'running_or_interrupted'
                    row['detail'] = 'No verified completion yet. Use sync after completion; no automatic retry.'
                previous_ok = row['status'] in {'reviewed', 'awaiting_review', 'completed', 'closed_no_selection'}
            else:
                previous_ok = False
            rows.append(row)
    return rows


def plan(project, tid, stage, new_attempt=False, state=None):
    state = state or project.load(); project.target(tid)
    current = project.current(state, tid, stage)
    parent = decision = None
    if stage != 'stage1':
        parent = project.current(state, tid, STAGES[STAGES.index(stage)-1])
        if not parent or not parent.get('checkpoint'):
            raise ValueError('Complete or attach '+STAGES[STAGES.index(stage)-1]+' first')
        verify_record(parent, True)
        parent_status = next(r for r in status_rows(project, state) if r['target'] == tid and r['stage'] == parent['stage'])
        if parent_status['status'] in {'stale', 'integrity_error', 'completed_with_errors'}:
            raise ValueError('Upstream checkpoint is '+parent_status['status'])
        decision = project.decision(state, parent)
        if not decision or not decision['ready']:
            raise ValueError('Import complete explicit '+parent['stage']+' decisions before continuing')
        # Revalidate against native scientific rules at every handoff.
        from .project import validate_decision
        validate_decision(project, tid, parent, project.root/decision['object'])
    output = project.store/'runs'/tid/stage/uuid4().hex
    args = native_args(project, tid, stage, out=output,
                       source=parent['output'] if parent else None,
                       decision=project.root/decision['object'] if decision else None)
    inputs = input_files(project, tid, args)
    config = dict(target=project.target(tid), options=project.options(tid, stage), inputs=inputs,
                  source_run_id=parent['id'] if parent else None,
                  source_seal_sha256=parent['checkpoint']['seal_sha256'] if parent else None,
                  decision_sha256=decision['sha256'] if decision else None)
    key = signature(config)
    if current and not new_attempt:
        if current.get('config_sha256') == key and current.get('checkpoint') and current['status'] not in {'failed', 'completed_with_errors'}:
            verify_record(current, True)
            return dict(action='reuse', target=tid, stage=stage, output=current['output'], run_id=current['id'])
        raise ValueError('A result/attempt already exists. Use --new-attempt to keep it and create another; imported results are never rerun automatically.')
    if stage == 'stage2':
        from .stage1 import verify
        from .stage2 import planned_regions, window_plan
        source, data = verify(parent['output'])
        regions, selection, _, intervals = planned_regions(source, data, args, planning=True)
        windows = window_plan(source, data, args, regions, selection, intervals)
        if not windows['ready_for_prediction']:
            raise ValueError('STOP2 window plan needs review:\n'+'\n'.join(b['message'] for b in windows['execution_blockers'] if not b['waived']))
    else:
        windows = None
    return dict(action='run', target=tid, stage=stage, output=str(output),
                argv=[sys.executable, '-m', 'genecontext', *argv_for(args)],
                config_sha256=key, inputs=inputs, source_run_id=config['source_run_id'],
                decision_sha256=config['decision_sha256'], windows=windows,
                executor=project.profile['mode'] if stage == 'stage2' else 'local',
                next_stop=stage, automatic_next_stage=False)


def queue_spec(project, prepared):
    profile = project.profile
    source = path_from(project.root, profile['source'])
    if not (source/'genecontext').is_dir():
        raise ValueError('Queue source must be a deployed immutable source checkout')
    active = Path(__file__).parent
    for path in active.iterdir():
        if path.is_file() and path.suffix in {'.py', '.txt'}:
            deployed = source/'genecontext'/path.name
            if not deployed.is_file() or sha256(deployed) != sha256(path):
                raise ValueError('Deploy the same reviewed code to a new immutable source before enqueueing; keep live pinned sources unchanged: '+path.name)
    root = path_from(project.root, profile['root'])
    if not (root/'state.json').is_file():
        raise ValueError('Use the existing initialized queue; project commands never create a dispatcher')
    files = {str(p.relative_to(source)): sha256(p) for p in (source/'genecontext').rglob('*')
             if p.is_file() and p.suffix in {'.py', '.txt', '.js', '.json'} and '__pycache__' not in p.parts}
    command = list(prepared['argv'])
    command[0] = str(path_from(project.root, profile['python'])) if profile.get('python') else sys.executable
    windows = prepared['windows']; ids = [r['locus_id'] for r in windows['regions']]
    spec = dict(panel=prepared['target'], name=project.data['name']+' / '+prepared['target'],
                job_id=prepared['config_sha256'][:24], loci=len(ids), locus_ids=ids,
                source_review_id=windows['source_review_id'], source=str(source), source_files=files,
                command=command)
    return spec


def run_stage(project, tid, stage, new_attempt=False):
    with project.locked() as state:
        prepared = plan(project, tid, stage, new_attempt, state)
        if prepared['action'] == 'reuse':
            return prepared
        if prepared['executor'] == 'plan':
            raise ValueError('STOP2 is plan-only. Configure the existing shared queue, or explicitly choose local execution on a separate machine.')
        code = {str(p.relative_to(Path(__file__).parent)): sha256(p)
                for p in Path(__file__).parent.rglob('*') if p.is_file() and p.suffix in {'.py', '.txt'} and '__pycache__' not in p.parts}
        record = dict(prepared, id=uuid4().hex, created=now(), status='running',
                      host=socket.gethostname(), pid=os.getpid(), code_sha256=code, imported=False)
        record.pop('windows', None)
        if prepared['executor'] == 'queue':
            spec = queue_spec(project, prepared)
            spec_path = project.store/'queue_specs'/(record['id']+'.json')
            atomic(spec_path, spec)
            record['queue_spec'] = str(spec_path); record['queue_root'] = str(path_from(project.root, project.profile['root']))
            record['status'] = 'submitting'
        state['runs'].append(record); state['current'].setdefault(tid, {})[stage] = record['id']
        project.save(state)
        try:
            if prepared['executor'] == 'queue':
                entrypoint = path_from(project.root, project.profile['entrypoint'])
                receipt = subprocess.run([spec['command'][0], str(entrypoint), 'enqueue', '--root', record['queue_root'], '--spec', str(spec_path)],
                                         check=True, capture_output=True, text=True, timeout=60, cwd=spec['source'],
                                         env=dict(os.environ, PYTHONPATH=spec['source']))
                record['queue_receipt'] = json.loads(receipt.stdout); record['status'] = 'queued'
                # The queue may deduplicate to an earlier output. Follow the
                # canonical ticket rather than inventing another computation.
                queue = json.loads((Path(record['queue_root'])/'state.json').read_text())
                job = next(j for j in queue['jobs'] if j['id'] == record['queue_receipt']['id'])
                record['output'] = job['output']
            else:
                log = project.store/'logs'/(record['id']+'.log'); log.parent.mkdir(parents=True, exist_ok=True)
                record['log'] = str(log)
                env = dict(os.environ, PYTHONPATH=str(Path(__file__).parent.parent), PYTHONUNBUFFERED='1')
                with log.open('x') as handle:
                    completed = subprocess.run(prepared['argv'], cwd=project.root, env=env, stdout=handle, stderr=subprocess.STDOUT)
                record['returncode'] = completed.returncode
                info = checkpoint(stage, record['output']) if (Path(record['output'])/'seal.json').exists() else None
                if info:
                    record['checkpoint'] = info
                record['status'] = info['status'] if info and completed.returncode == 0 else 'failed'
                if completed.returncode:
                    record['error'] = 'Native stage failed; inspect '+str(log)
        except BaseException as error:
            record.update(status='failed', error=str(error))
            raise
        finally:
            record['updated'] = now(); project.save(state)
        return record


def sync(project):
    changed = []
    with project.locked() as state:
        for run in state['runs']:
            if run.get('checkpoint') or run.get('status') not in {'running', 'queued', 'submitting'}:
                continue
            if run.get('queue_root'):
                queue = json.loads((Path(run['queue_root'])/'state.json').read_text())
                ticket = run.get('queue_receipt', {}).get('id')
                if not ticket:
                    continue
                job = next(j for j in queue['jobs'] if j['id'] == ticket)
                if job['status'] not in {'awaiting_stop2_review', 'completed_with_errors', 'completed', 'failed', 'needs_context_review'}:
                    continue
            if (Path(run['output'])/'seal.json').exists():
                info = checkpoint(run['stage'], run['output'])
                run.update(checkpoint=info, status=info['status'], updated=now()); changed.append(run['id'])
        project.save(state)
    return dict(verified=changed, calculations_started=False)


def register_cli(commands):
    group = commands.add_parser('project', help='Configure a study, attach existing results and run one reviewed STOP at a time')
    sub = group.add_subparsers(dest='project_action', required=True)
    p = sub.add_parser('init', help='Create a study configuration; no calculations')
    p.add_argument('--out', required=True); p.add_argument('--name', required=True)
    p.add_argument('--genomes', required=True); p.add_argument('--targets', required=True, help='TSV: id, gene, query, query_type, optional query_metadata')
    p.add_argument('--stage2-executor', choices=['plan', 'local'], default='plan')
    for name in ['validate', 'status', 'plan', 'run', 'attach', 'decisions', 'sync', 'report']:
        p = sub.add_parser(name); p.add_argument('--project', required=True, help='project.json or its directory')
        if name in {'plan', 'run', 'attach', 'decisions'}:
            p.add_argument('--target', required=True, help='Target ID; run/plan also accept all')
            p.add_argument('--stage', choices=STAGES if name != 'decisions' else STAGES[:3], required=True)
        if name in {'plan', 'run', 'attach'}:
            p.add_argument('--new-attempt', action='store_true', help='Retain the current attempt and explicitly create/select another')
        if name == 'attach':
            p.add_argument('--results', required=True)
        if name == 'decisions':
            p.add_argument('--file', required=True)
        if name == 'status':
            p.add_argument('--json', action='store_true')


def main(args):
    if args.project_action == 'init':
        result = init(args)
    else:
        project = Project(args.project); action = args.project_action
        if action == 'status':
            rows = status_rows(project)
            if not args.json:
                for row in rows:
                    print(f"{row['target']}\t{row['stage']}\t{row['status']}\t{row.get('detail', '')}")
                return 0
            result = rows
        elif action == 'validate':
            for tid in project.targets:
                # Check configured CLI options even before any results exist.
                for stage in STAGES:
                    options_argv(project, tid, stage)
                    native_args(project, tid, stage, out=project.store/'unused',
                                source=project.store/'unused-source', decision=project.store/'unused.tsv' if stage != 'stage1' else None,
                                validation=True)
                args1 = native_args(project, tid, 'stage1', out=project.store/'unused')
                input_files(project, tid, args1)
                from .stage1_proteins import discovery_rows
                from .stage1 import panel
                discovery_rows(args1); panel(args1)
            for run in project.load()['runs']:
                if run.get('checkpoint'):
                    verify_record(run, full=True)
            errors = [r for r in status_rows(project) if r['status'] == 'integrity_error']
            if errors:
                raise ValueError(json.dumps(errors))
            result = dict(valid=True, targets=len(project.targets), calculations_started=False)
        elif action in {'plan', 'run'}:
            tids = list(project.targets) if args.target == 'all' else [args.target]
            result = [(plan if action == 'plan' else run_stage)(project, tid, args.stage, args.new_attempt) for tid in tids]
        elif action == 'attach':
            result = attach(project, args.target, args.stage, args.results, args.new_attempt)
        elif action == 'decisions':
            result = record_decision(project, args.target, args.stage, args.file)
        elif action == 'sync':
            result = sync(project)
        elif action == 'report':
            from .project_report import write_report
            result = write_report(project)
        else:
            raise ValueError('Unknown project action')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 1 if isinstance(result, list) and any(r.get('status') == 'failed' for r in result) else 0
