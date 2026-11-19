"""Prepare exact-sequence AF3 requests without provisioning or spending credits."""
import hashlib
import json
from pathlib import Path
import re

from .stage1 import dump, new_output, now
from .sequence import sha256


def prepare(sequences, out, limit=None, labels=None, approval=None):
    from .structure_approval import verify
    approval=verify(approval,sequences)
    if limit is not None and (type(limit) is not int or limit<1):
        raise ValueError('AF3 protein limit must be a positive integer or omitted for all')
    root=new_output(out);rows=[];requests=[];seen={}
    dump(root/'final-reference-approval.json',approval)
    for sid,sequence in sequences.items():
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}',sid):
            raise ValueError('Unsafe structure sequence ID')
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        record=seen.get(digest)
        if record is None:
            if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
                record=dict(status='unsupported_sequence',reason='AF3 pilot accepts canonical protein monomers')
            elif not 4<=len(sequence)<=5000:
                record=dict(status='unsupported_sequence',reason='AlphaFold Server monomer length must be 4–5000 residues')
            elif limit is not None and len(requests)>=limit:
                record=dict(status='not_requested',reason='Pilot unique-protein limit reached')
            else:
                name=(labels or {}).get(sid) or 'gcf_'+digest[:24]
                if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z0-9_.-]{1,96}',name):
                    raise ValueError('AF3 job label must be 1–96 ASCII letters, numbers, dots, dashes or underscores')
                job=dict(name=name,dialect='alphafold3',version=4,modelSeeds=[42],
                         sequences=[dict(protein=dict(id='A',sequence=sequence))])
                filename=digest+'.af3.json';dump(root/filename,job)
                server_file=digest+'.server.json'
                dump(root/server_file,[dict(name=job['name'],modelSeeds=['42'],dialect='alphafoldserver',version=1,
                    sequences=[dict(proteinChain=dict(sequence=sequence,count=1))])])
                # A Cloud endpoint's output_dir and credentials belong to its
                # configured account. Do not invent an endpoint or bucket here.
                record=dict(status='awaiting_prediction',reason='Upload the Server JSON and run the draft in AlphaFold Server',
                            job_name=name,provider='alphafold_server',model='AlphaFold 3',adapter_version=1,
                            input_file=filename,input_sha256=sha256(root/filename),protein_sha256=digest,
                            server_input_file=server_file,server_input_sha256=sha256(root/server_file))
                requests.append(dict(record,source_sequence_ids=[]))
            seen[digest]=record
        if 'input_file' in record:
            next(r for r in requests if r['protein_sha256']==digest)['source_sequence_ids'].append(sid)
        rows.append(dict(sequence_id=sid,protein_sha256=digest,**{k:v for k,v in record.items() if k!='protein_sha256'}))
    batch=[]
    for request in requests:
        batch.extend(json.loads((root/request['server_input_file']).read_text()))
    uploads=[]
    if len(batch)<=100:
        dump(root/'upload.server.json',batch);uploads.append('upload.server.json')
    else:
        (root/'batches').mkdir()
        for offset in range(0,len(batch),100):
            filename='batches/upload_%03d.server.json'%(offset//100+1)
            dump(root/filename,batch[offset:offset+100]);uploads.append(filename)
    (root/'README.txt').write_text('AlphaFold Server\n\n1. Sign in at https://alphafoldserver.com/\n2. Use Upload JSON with these files: '+', '.join(uploads)+'\n3. Open each draft, check its name/sequence and seed, and submit within the account quota.\n4. The Server calculates five samples per job. Display one best structure per protein; retain the original ZIP as evidence.\n5. Import each completed ZIP into a new directory with gene-context af3-import.\n\nNo predictions have been submitted by this preparation command.\n')
    dump(root/'requests.json',dict(schema='gcf_af3_requests_v1',created=now(),requests=requests,rows=rows,
        cloud_envelope=dict(instances='Use one prepared AF3 JSON object per instance',parameters=dict(output_dir='Configure the owning project Cloud Storage output URI')),
        prediction_submitted=False,provider='alphafold_server',transport='manual_json_upload',
        final_reference_approval_sha256=sha256(root/'final-reference-approval.json'),
        model='AlphaFold 3',adapter_version=1,upload_files=uploads,display_policy='one_best_structure_per_protein',extensions={},source_code_sha256=sha256(Path(__file__)),
        documentation='https://docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphafold-3'))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p.name!='seal.json'})
    return rows
