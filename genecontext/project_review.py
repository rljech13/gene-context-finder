"""Sequential review workspaces above immutable scientific checkpoints.

A review is an editable presentation, never an annotation run or inherited
approval. A handoff is accepted only against its exact workspace and native
model/context evidence. Historical reviews and exports remain in the ledger.
"""
import json
from pathlib import Path
from uuid import uuid4

from .project import atomic, now, verify_record
from .sequence import sha256


def current_review(project, state, tid, full=False):
    rid=state.get('selected_reviews',{}).get(tid)
    review=next((r for r in state.get('reviews',[]) if r['id']==rid),None)
    if not review:return None
    parent=project.current(state,tid,'stage2')
    if not parent or parent['id']!=review['source_run_id']:
        raise ValueError('The review belongs to an earlier prepared model set')
    verify_record(parent,full)
    if parent['checkpoint']['seal_sha256']!=review['source_seal_sha256']:
        raise ValueError('Prepared model source changed')
    root=project.root/review['directory']
    if sha256(root/'seal.json')!=review['seal_sha256']:
        raise ValueError('Review seal changed')
    seal=json.loads((root/'seal.json').read_text())
    if full:
        for name,digest in seal.items():
            path=(root/name).resolve()
            if not path.is_relative_to(root.resolve()) or sha256(path)!=digest:
                raise ValueError('Review presentation changed: '+name)
    return review


def prepare(project,tid,new_attempt=False):
    project.target(tid)
    with project.locked() as state:
        parent=project.current(state,tid,'stage2')
        if not parent:raise ValueError('Complete or attach stage2 model preparation first')
        verify_record(parent,True)
        if parent['checkpoint']['status']=='completed_with_errors':
            raise ValueError('Resolve model-preparation errors before review')
        if state.get('selected_reviews',{}).get(tid) and not new_attempt:
            review=current_review(project,state,tid,True)
            return dict(action='reuse',**review)
        from .alignment_review import build
        identity=uuid4().hex;out=project.store/'reviews'/tid/identity
        root,pool=build(parent['output'],out,legacy=True,title=project.target(tid)['gene'])
        entry=dict(id=identity,target=tid,created=now(),directory=str(root.relative_to(project.root)),
            source_run_id=parent['id'],source_seal_sha256=parent['checkpoint']['seal_sha256'],
            source_review_id=pool['source_review_id'],alignment_review_id=pool['review_id'],
            seal_sha256=sha256(root/'seal.json'),report=str(root/'report.html'),calculations_started=False)
        state.setdefault('reviews',[]).append(entry);state.setdefault('selected_reviews',{})[tid]=identity
        project.save(state)
        return dict(action='created',**entry)


def verify_handoff(project,state,tid,path):
    from .stage4_handoff import read_handoff
    from .stage2_choices import verify_checkpoint,validate_choices,issue_text,constraints
    review=current_review(project,state,tid,True)
    if not review:raise ValueError('Create a project review workspace first')
    before=sha256(path);data,members=read_handoff(path)
    if (data['source_review_id'],data['alignment_review_id'])!=(review['source_review_id'],review['alignment_review_id']):
        raise ValueError('Handoff belongs to another review workspace')
    pool=json.loads((project.root/review['directory']/'review_data.json').read_text())
    exact={s['id']:s for s in pool['sequences']}
    for mid,member in members.items():
        if mid not in exact or member!=exact[mid]:
            raise ValueError('Handoff member differs from the prepared workspace: '+mid)
    parent=project.current(state,tid,'stage2');_,native=verify_checkpoint(parent['output'])
    issues=validate_choices(native,data['stage2'])
    if issues:raise ValueError('\n'.join(issue_text(native,issue) for issue in issues))
    regions={r['region_id']:r for r in native['regions']}
    for region in data['stage3']['synteny']['regions']:
        if region!=regions.get(region['region_id']):raise ValueError('Synteny region differs from prepared evidence')
    retained={r['region_id'] for r in data['stage3']['synteny']['regions']}
    rules=constraints(native)
    expected={m['model_id']:dict(model_id=m['model_id'],region_id=m['region_id'],
        physical_id=rules[m['model_id']]['physical_id'],source=m['source'],contig=m['contig'],
        start=m['start'],end=m['end'],strand=m['strand']) for m in native['models'] if m['region_id'] in retained}
    observed={m['model_id']:{k:m.get(k) for k in expected[next(iter(expected))]} for m in data['stage3']['synteny'].get('neighborhood',[])}
    if len(observed)!=len(data['stage3']['synteny'].get('neighborhood',[])) or observed!=expected:
        raise ValueError('Synteny neighborhood differs from prepared evidence')
    pairs=native.get('neighbor_pairs',[])
    if set(data['stage3']['synteny']['reviewed_pairs'])!=set(data['stage2'].get('pairs',{})):
        raise ValueError('Reviewed synteny correspondences differ from saved decisions')
    for pair in data['stage3']['synteny']['reviewed_pairs'].values():
        if pair not in pairs:raise ValueError('Synteny correspondence differs from prepared evidence')
    if sha256(path)!=before:raise ValueError('Handoff changed during validation')
    return review,before


def record_handoff(project,tid,path):
    path=Path(path).expanduser().resolve()
    with project.locked() as state:
        review,digest=verify_handoff(project,state,tid,path)
        obj,saved=project.snapshot_decision(path)
        if saved!=digest:raise ValueError('Handoff changed while saving')
        entry=dict(target=tid,review_id=review['id'],sha256=digest,object=str(obj.relative_to(project.root)),created=now())
        if not any(h['sha256']==digest and h['review_id']==review['id'] for h in state.get('handoffs',[])):
            state.setdefault('handoffs',[]).append(entry)
        state.setdefault('selected_handoffs',{})[tid]=dict(sha256=digest,review_id=review['id'])
        project.save(state);return entry


def handoff(project,state,tid,full=False):
    selection=state.get('selected_handoffs',{}).get(tid)
    if not selection:return None
    review=current_review(project,state,tid,full)
    if not review or selection['review_id']!=review['id']:raise ValueError('Final handoff belongs to an earlier review')
    entry=next(h for h in state['handoffs'] if h['sha256']==selection['sha256'] and h['review_id']==review['id'])
    path=project.root/entry['object']
    if sha256(path)!=entry['sha256']:raise ValueError('Saved final handoff changed')
    if full:verify_handoff(project,state,tid,path)
    return entry


def plan_final(project,tid,new_attempt=False,state=None):
    from .project import signature
    from .project_run import argv_for,options_argv,input_files
    from .cli import parser
    import sys
    state=state or project.load();entry=handoff(project,state,tid,True)
    if not entry:raise ValueError('Import the explicit final STOP3 export with project handoff first')
    review=current_review(project,state,tid,True)
    args=parser().parse_args(['stage4','--handoff',str(project.root/entry['object']),'--out',str(project.store/'runs'/tid/'stage4'/uuid4().hex),*options_argv(project,tid,'stage4')])
    inputs=dict(input_files(project,tid,args),handoff_sha256=entry['sha256'])
    key=signature(dict(target=project.target(tid),handoff=entry['sha256'],review=review['id'],options=project.options(tid,'stage4'),inputs=inputs))
    current=project.current(state,tid,'stage4')
    if current and not new_attempt:
        if current.get('config_sha256')==key and current.get('checkpoint'):
            verify_record(current,True);return dict(action='reuse',target=tid,stage='stage4',output=current['output'],run_id=current['id'])
        raise ValueError('A final result/attempt exists; --new-attempt retains it in history')
    return dict(action='run',target=tid,stage='stage4',output=args.out,argv=[sys.executable,'-m','genecontext',*argv_for(args)],config_sha256=key,
        source_run_id=review['source_run_id'],decision_sha256=entry['sha256'],inputs=inputs,windows=None,executor='local',next_stop='stage4',automatic_next_stage=False)


def status(project,state,rows):
    for tid in project.targets:
        selected={r['stage']:r for r in rows if r['target']==tid}
        try:
            review=current_review(project,state,tid)
            final=handoff(project,state,tid) if review else None
            if review:
                selected['stage2']['report']=review['report']
                selected['stage3'].update(status='reviewed' if final else 'awaiting_synteny_review',report=review['report']+'?view=synteny')
            else:selected['stage3']['status']='blocked'
            run=project.current(state,tid,'stage4')
            if run and run.get('checkpoint'):
                if not final or run.get('decision_sha256')!=final['sha256']:
                    selected['stage4'].update(status='stale',detail='Final handoff changed; previous result retained.')
            elif final:selected['stage4']['status']='ready'
        except (ValueError,OSError,KeyError) as error:
            selected['stage3'].update(status='stale',detail=str(error));selected['stage4']['status']='stale' if project.current(state,tid,'stage4') else 'blocked'
    return rows
