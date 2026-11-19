"""Prepare exact-sequence AF3 requests without provisioning or spending credits."""
import hashlib
from pathlib import Path
import re

from .stage1 import dump, new_output, now
from .sequence import sha256


def prepare(sequences, out, limit=3):
    if not 1<=limit<=50:
        raise ValueError('AF3 pilot limit must be between 1 and 50 unique proteins')
    root=new_output(out);rows=[];requests=[];seen={}
    for sid,sequence in sequences.items():
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}',sid):
            raise ValueError('Unsafe structure sequence ID')
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        record=seen.get(digest)
        if record is None:
            if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
                record=dict(status='unsupported_sequence',reason='AF3 pilot accepts canonical protein monomers')
            elif len(seen)>=limit:
                record=dict(status='not_requested',reason='Pilot unique-protein limit reached')
            else:
                job=dict(name='gcf_'+digest[:24],dialect='alphafold3',version=4,modelSeeds=[42],
                         sequences=[dict(protein=dict(id='A',sequence=sequence))])
                filename=digest+'.af3.json';dump(root/filename,job)
                server_file=digest+'.server.json'
                dump(root/server_file,[dict(name=job['name'],modelSeeds=[42],dialect='alphafoldserver',version=1,
                    sequences=[dict(proteinChain=dict(sequence=sequence,count=1))])])
                # A Cloud endpoint's output_dir and credentials belong to its
                # configured account. Do not invent an endpoint or bucket here.
                record=dict(status='awaiting_prediction',reason='AF3 request prepared; endpoint/account not configured',
                            input_file=filename,input_sha256=sha256(root/filename),protein_sha256=digest,
                            server_input_file=server_file,server_input_sha256=sha256(root/server_file))
                requests.append(dict(record,source_sequence_ids=[]))
            seen[digest]=record
        if 'input_file' in record:
            next(r for r in requests if r['protein_sha256']==digest)['source_sequence_ids'].append(sid)
        rows.append(dict(sequence_id=sid,protein_sha256=digest,**{k:v for k,v in record.items() if k!='protein_sha256'}))
    dump(root/'requests.json',dict(schema='gcf_af3_requests_v1',created=now(),requests=requests,rows=rows,
        cloud_envelope=dict(instances='Use one prepared AF3 JSON object per instance',parameters=dict(output_dir='Configure the owning project Cloud Storage output URI')),
        prediction_submitted=False,source_code_sha256=sha256(Path(__file__)),
        documentation='https://docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphafold-3'))
    dump(root/'seal.json',{str(p.relative_to(root)):sha256(p) for p in root.iterdir() if p.is_file() and p.name!='seal.json'})
    return rows
