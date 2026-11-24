"""Only explicitly reviewed final Stage 3 representatives may leave for folding."""
import hashlib
import json
from pathlib import Path

from .sequence import sha256
from .stage3_finalization import validate_choices,reference_approved

SCHEMA='gcf_final_reference_approval_v1'


def from_validated(source,evidence,state,proteins):
    state=validate_choices(state,evidence['spec']);approved={}
    for bucket in evidence['spec']['buckets']:
        choice=state['choices'][bucket['key']];mid=choice['primary']
        if mid and reference_approved(choice,bucket,evidence['spec'],hashlib.sha256(proteins[mid].encode()).hexdigest()):
            approved[mid]=hashlib.sha256(proteins[mid].encode()).hexdigest()
    return dict(schema=SCHEMA,stage=3,policy='explicitly_reviewed_final_primary_reference',
        source_seal_sha256=sha256(Path(source)/'seal.json'),spec=evidence['spec'],decisions=state,proteins=approved)


def verify(approval,sequences):
    if isinstance(approval,(str,Path)):
        approval=json.loads(Path(approval).read_text())
    if isinstance(approval,dict) and approval.get('schema')=='gcf_synteny_reference_approval_v1':
        from .synteny_reference import verify as verify_synteny
        return verify_synteny(approval,sequences)
    if not isinstance(approval,dict) or approval.get('schema')!=SCHEMA or approval.get('stage')!=3:
        raise ValueError('Folding export requires explicit final reference-model approval after Stage 3; saved candidate selections and alignments are insufficient')
    state=validate_choices(approval['decisions'],approval['spec'])
    allowed={state['choices'][b['key']]['primary'] for b in approval['spec']['buckets']
        if reference_approved(state['choices'][b['key']],b,approval['spec'],
            approval.get('proteins',{}).get(state['choices'][b['key']]['primary']))}
    if not sequences or set(sequences)-allowed:
        raise ValueError('Only explicitly reviewed final Stage 3 primary references may be sent for folding')
    if not approval.get('source_seal_sha256') or any(hashlib.sha256(seq.encode()).hexdigest()!=approval.get('proteins',{}).get(sid) for sid,seq in sequences.items()):
        raise ValueError('Final reference approval does not match these exact protein sequences')
    return approval


def issue(source,out,decisions=None):
    from .stage4 import checked_source,load_models
    source,_=checked_source(source)
    evidence=json.loads((source/'finalization.evidence.json').read_text())
    state=json.loads(Path(decisions or source/'finalization.decisions.json').read_text())
    comparison=json.loads((source/'comparison.json').read_text())
    _,proteins,_=load_models(source,comparison,evidence)
    approval=from_validated(source,evidence,state,proteins)
    if not approval['proteins']:raise ValueError('No explicitly approved final Stage 3 reference models; no folding export is permitted')
    target=Path(out)
    if target.exists():raise ValueError('Use a new final approval file; previous decisions are preserved')
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(approval,indent=2)+'\n')
    return approval
