"""Retain sealed STOP-1 context independently of the chosen window policy."""
import json
from pathlib import Path
from .sequence import sha256


def collect(data, prior, prior_sha256):
    if prior.get('review_id') != data['source_review_id']:
        raise ValueError('Prior synteny belongs to another STOP-1 review')
    if prior_sha256 != data['source_evidence_sha256']:
        raise ValueError('Prior synteny source fingerprint differs')
    loci = {(l['genome_id'], l['locus_id']): l for l in prior['loci']}
    return dict(schema='stage2_prior_context_v1', source_review_id=prior['review_id'],
        source_evidence_sha256=prior_sha256,
        references=prior.get('early_synteny', {}).get('references', []),
        regions={r['region_id']: loci.get((r['genome_id'], r['locus_id']), {}).get('early_synteny', {})
                 for r in data['regions']})


def load(root, data):
    context = data.get('prior_synteny_context')
    sidecar = Path(root) / 'prior_synteny_context.json'
    if sidecar.exists():
        saved = json.loads(sidecar.read_text())
        if saved.get('stage2_evidence_sha256') != sha256(Path(root) / 'evidence.json'):
            raise ValueError('Prior synteny supplement belongs to another checkpoint')
        context = saved['context']
    if context is None and data.get('source_checkpoint'):
        source = Path(data['source_checkpoint']) / 'evidence.json'
        if source.is_file():
            context = collect(data, json.loads(source.read_text()), sha256(source))
    if context is not None and (context.get('schema') != 'stage2_prior_context_v1' or
            context.get('source_review_id') != data['source_review_id'] or
            context.get('source_evidence_sha256') != data['source_evidence_sha256']):
        raise ValueError('Prior synteny provenance differs')
    return context
