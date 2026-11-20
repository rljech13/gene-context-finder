"""Prepare every selected protein once, in bounded Server upload batches."""
from collections import defaultdict, deque
import hashlib
import html
import json
from pathlib import Path
import re

from .af3_requests import prepare as prepare_requests
from .sequence import sha256, write_tsv
from .stage1 import dump, new_output, now
from .stage4_structure import ca_atoms


def verified_result(path):
    path=Path(path).resolve();receipt=json.loads((path/'receipt.json').read_text())
    if receipt.get('status')!='available' or receipt.get('provider')!='alphafold_server':
        raise ValueError('Expected a completed AlphaFold Server result')
    seal=json.loads((path/'seal.json').read_text())
    for name,digest in seal.items():
        file=(path/name).resolve()
        if path not in file.parents or sha256(file)!=digest:raise ValueError('Reused Server result changed: '+name)
    if sha256(path/'source.zip')!=receipt['source_archive_sha256']:raise ValueError('Source archive changed')
    sequence=''.join(a['aa'] for a in ca_atoms(path/'model.pdb'))
    if hashlib.sha256(sequence.encode()).hexdigest()!=receipt['protein_sha256']:raise ValueError('Reused protein differs')
    return receipt,dict(result_directory=str(path),receipt_sha256=sha256(path/'receipt.json'),
        catalogue_sha256=sha256(path/'catalogue.json'),job_source_url=receipt.get('job_source_url',''))


def prepare(sequences,out,metadata=None,completed=(),batch_size=30,first_batch_size=None,approval=None):
    from .structure_approval import verify
    approval=verify(approval,sequences)
    if not 1<=batch_size<=100 or first_batch_size is not None and not 1<=first_batch_size<=100:
        raise ValueError('Server upload batch sizes must be 1–100 jobs')
    metadata=metadata or {};labels={}
    for sid,sequence in sequences.items():
        m=metadata.get(sid,{})
        prefix='_'.join(str(m.get(k,'')) for k in ['gene','species']).strip('_') or sid
        prefix=re.sub('[^A-Za-z0-9_.-]','_',prefix)[:70]
        labels[sid]='GCF_'+prefix+'_'+hashlib.sha256(sequence.encode()).hexdigest()[:12]
    # Reuse is by full sequence, not accession, species, or proposed gene name.
    ready={}
    for folder in completed:
        receipt,reference=verified_result(folder)
        ready[receipt['protein_sha256']]=reference
    root=new_output(out)
    rows=prepare_requests(sequences,root/'requests',labels=labels,approval=approval)
    requests=json.loads((root/'requests/requests.json').read_text())['requests']
    jobs=[];queues=defaultdict(deque)
    for request in requests:
        aliases=request['source_sequence_ids'];digest=request['protein_sha256']
        job=dict(protein_sha256=digest,job_name=request['job_name'],source_sequence_ids=aliases,
            metadata=[dict(metadata.get(sid,{}),sequence_id=sid) for sid in aliases],
            server_input_file='requests/'+request['server_input_file'],
            status='available' if digest in ready else 'awaiting_submission')
        if digest in ready:job['result']=ready[digest]
        jobs.append(job)
        if digest not in ready:queues[str(metadata.get(aliases[0],{}).get('gene','all'))].append(job)
    ordered=[]
    while any(queues.values()):
        for queue in queues.values():
            if queue:ordered.append(queue.popleft())
    (root/'batches').mkdir();uploads=[];offset=0
    while offset<len(ordered):
        size=(first_batch_size or batch_size) if not uploads else batch_size
        part=ordered[offset:offset+size];name='batches/upload_%03d.server.json'%(len(uploads)+1)
        payload=[]
        for job in part:
            payload.extend(json.loads((root/job['server_input_file']).read_text()));job['upload_file']=name
        dump(root/name,payload);uploads.append(dict(file=name,jobs=len(part),sha256=sha256(root/name)));offset+=size
    unsupported=[r for r in rows if r['status']=='unsupported_sequence']
    plan=dict(schema='gcf_af3_campaign_v1',created=now(),provider='alphafold_server',model='AlphaFold 3',
        display_policy='one_best_structure_per_protein',server_samples_per_job=5,seed=42,jobs=jobs,
        unsupported=unsupported,uploads=uploads,source_code_sha256=sha256(Path(__file__)),
        counts=dict(selected_sequences=len(sequences),unique_requested=len(jobs),
            available=sum(j['status']=='available' for j in jobs),awaiting_submission=len(ordered),unsupported=len(unsupported)),
        note='Upload batches are not automatic daily schedules. Read the current Server quota before submitting; each job consumes one quota unit.')
    dump(root/'campaign.json',plan)
    export=[]
    for job in jobs:
        for m in job['metadata']:
            export.append(dict(sequence_id=m['sequence_id'],gene=m.get('gene',''),species=m.get('species',''),
                original_model=m.get('original_model',''),protein_sha256=job['protein_sha256'],
                job_name=job['job_name'],status=job['status'],upload_file=job.get('upload_file','')))
    write_tsv(root/'jobs.tsv',export,list(export[0]) if export else ['sequence_id','status'])
    links=''.join('<li><a href="'+u['file']+'">'+html.escape(u['file'])+'</a> — '+str(u['jobs'])+' jobs</li>' for u in uploads)
    (root/'index.html').write_text('<meta charset="utf-8"><title>AlphaFold Server campaign</title><h1>All selected proteins</h1><pre>'+html.escape(json.dumps(plan['counts'],indent=2))+'</pre><p>One best structure is displayed per protein. Identical proteins share a job; completed exact proteins are not resubmitted.</p><ol>'+links+'</ol><p>Upload JSON batches in AlphaFold Server, then run the saved drafts within the account quota. This preparation has submitted no jobs.</p><a href="jobs.tsv">All models and job names</a>')
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    return plan
