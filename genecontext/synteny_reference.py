"""Add final approval to the separate synteny stop without editing source data."""
from pathlib import Path
import json

MARKER = '// gcf-synteny-reference-approval-v1'


def attach(html):
    """Small additive hook for alignment-first report publishers, idempotent."""
    if MARKER in html:
        return html
    if 'gcf-separate-synteny-stop-v2' not in html:
        return html
    script = Path(__file__).with_name('synteny_reference_approval.txt').read_text()
    if html.count('</body>') != 1:
        raise ValueError('Expected one report body')
    return html.replace('</body>', '<script>'+script+'</script></body>')


def verify(approval, sequences):
    """Require explicit post-synteny approval of exact selected native proteins."""
    if approval.get('stage') != 3:
        raise ValueError('Expected final Stage 3 synteny approval')
    handoff = approval.get('handoff', {})
    decisions = handoff.get('decisions', {})
    source = approval.get('source_review_id')
    review = approval.get('alignment_review_id')
    if not source or not review or handoff.get('review_id') != source or decisions.get('review_id') != source or decisions.get('schema') != 'stage2_draft_v1' or handoff.get('alignment_review_id') != review:
        raise ValueError('Synteny approval belongs to another confirmed model set')
    canonical = json.dumps(dict(handoff=handoff,synteny_pairs=approval.get('synteny_pairs',{})), sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    import hashlib
    seen_physical = set()
    for mid, sequence in sequences.items():
        a = approval.get('references', {}).get(mid, {})
        choice = decisions.get('models', {}).get(mid, {})
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        region = decisions.get('regions', {}).get(a.get('region_id'), {})
        physical = a.get('physical_id')
        if not physical or physical in seen_physical:
            raise ValueError('Only one final reference per physical locus')
        seen_physical.add(physical)
        if (a.get('schema') != 'synteny_reference_approval_v1' or a.get('stage') != 3
                or a.get('confirmed') is not True or a.get('synteny_reviewed') is not True
                or a.get('model_id') != mid or a.get('source_review_id') != source
                or a.get('alignment_review_id') != review or a.get('context_sha256') != hashlib.sha256(canonical.encode()).hexdigest()
                or region.get('saved') is not True or region.get('decision') != 'inspect'
                or choice.get('decision') != 'inspect' or choice.get('role') not in ('target', 'alignment_candidate')
                or a.get('protein_sha256') != digest or approval.get('proteins', {}).get(mid) != digest):
            raise ValueError('Explicit post-synteny approval does not match this model/protein')
    if not sequences:
        raise ValueError('No explicitly approved post-synteny proteins')
    return approval
