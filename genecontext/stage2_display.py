"""Read-only display fragments and concise automatic model comments."""
from collections import defaultdict
import hashlib
from pathlib import Path

from .sequence import attributes, union_length
from .stage2_quality import shared_blocks

POLICY = dict(size_ratio=2.0, min_span_difference_bp=10000, min_protein_difference_aa=50)


def _raw_models(path):
    """GFF identifiers may repeat on different prediction-window contigs."""
    models = {}
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        z = line.split('\t')
        if len(z) != 9:
            continue
        attrs = attributes(z[8])
        if z[2] in {'mRNA', 'transcript'}:
            models[z[0], attrs['ID']] = dict(original_model=attrs['ID'],
                original_gene=attrs.get('Parent', attrs['ID']), prediction_id=z[0],
                strand=z[6], start=int(z[3])-1, end=int(z[4]), cds_blocks=[])
        elif z[2] == 'CDS':
            for parent in attrs.get('Parent', '').split(','):
                if (z[0], parent) in models:
                    models[z[0], parent]['cds_blocks'].append([int(z[3])-1, int(z[4]), int(z[7])])
    return [m for m in models.values() if m['cds_blocks']]


def filtered_fragments(root, data):
    """Only removed predictions overlapping an existing target CDS are drawn.

    These records never enter DATA.models, model choices or sequence exports.
    No raw prediction is inferred when pre-filter output was not saved.
    """
    root = Path(root)
    intervals = {r['region_id']: r for r in data.get('prediction_intervals', [])}
    regions = {r['region_id']: r for r in data['regions']}
    result = []; seen = set()
    for path in sorted(root.glob('context_prediction/*/predictions.gff3.raw.gff3')):
        final = path.with_name('predictions.gff3')
        if not final.exists():
            continue
        signature = lambda m: (m['prediction_id'], m['strand'], tuple(sorted(tuple(b) for b in m['cds_blocks'])))
        kept = {signature(m) for m in _raw_models(final)}
        for raw in _raw_models(path):
            if signature(raw) in kept:
                continue
            interval = intervals.get(raw['prediction_id'])
            if not interval:
                continue
            offset = interval['window_start']
            blocks = [[a+offset, b+offset, p] for a,b,p in raw['cds_blocks']]
            for rid in interval['region_ids']:
                r = regions[rid]
                candidates = [m for m in data['models'] if m['region_id']==rid and
                    m.get('overlaps_anchor') and m['strand']==raw['strand']]
                ranked = sorted([(union_length(shared_blocks(blocks, m['cds_blocks'])), m['model_id'], m)
                                 for m in candidates], reverse=True, key=lambda v: v[:2])
                if not ranked or not ranked[0][0]:
                    continue
                overlap, _, target = ranked[0]
                key = (rid, signature(raw))
                if key in seen:
                    continue
                seen.add(key)
                nt = sum(b-a for a,b,p in blocks)
                reason = 'Below the predictor 200-nt coding-length filter' if nt < 200 else 'Removed by predictor postprocessing'
                identifier = 'filtered_' + hashlib.sha256(repr(key).encode()).hexdigest()[:24]
                result.append(dict(model_id=identifier, region_id=rid,
                    gene_group=target['gene_group'], genome_id=r['genome_id'],
                    source='tiberius', contig=r['contig'], strand=raw['strand'],
                    original_gene=raw['original_gene'], original_model=raw['original_model'],
                    start=raw['start']+offset, end=raw['end']+offset,
                    cds_blocks=blocks, protein_length=max(0,nt//3-1), cds_length=nt,
                    target_cds_overlap=overlap, annotation_name='', target_support=[],
                    identification=dict(proposed_name='',near_best=[]), overlaps_anchor=True,
                    diagnostic_fragment=True, reviewable=False, filter_reason=reason,
                    provenance=dict(source=str(path.relative_to(root)),
                        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        prediction_interval=interval['region_id'])))
    return result


def automatic_comments(data, fragments=()):
    """Size differences are comments, never automatic exclusions or approvals."""
    grouped = defaultdict(list); result = defaultdict(list)
    for m in data['models']:
        grouped[m['region_id'],m['gene_group']].append(m)
    for (rid, group), members in grouped.items():
        sources = {source: [m for m in members if m['source']==source]
                   for source in ['annotation','tiberius']}
        if not all(sources.values()):
            continue
        findings = []
        for kind, get, minimum, unit in [
            ('Genomic span',lambda m:max(b[1] for b in m['cds_blocks'])-min(b[0] for b in m['cds_blocks']),POLICY['min_span_difference_bp'],'kb'),
            ('Protein length',lambda m:m.get('protein_length',0),POLICY['min_protein_difference_aa'],'aa')]:
            values={source:max(get(m) for m in models) for source,models in sources.items()}
            low,high=min(values.values()),max(values.values())
            if low<=0 or high/low<POLICY['size_ratio'] or high-low<minimum:
                continue
            fmt=lambda n:f'{n/1000:.1f} kb' if unit=='kb' else f'{n} aa'
            findings.append(f"{kind} differs: annotation {fmt(values['annotation'])}; Tiberius {fmt(values['tiberius'])} ({high/low:.1f}× between source maxima).")
        if findings:
            name=next((m.get('annotation_name') for m in sources['annotation'] if m.get('annotation_name')), 'Gene')
            result[rid].append(name+': '+' '.join(findings)+' This does not establish a model error or change the selected model.')
    for m in fragments:
        result[m['region_id']].append(f"Tiberius also predicted a {m['cds_length']}-nt fragment overlapping {m['target_cds_overlap']}-nt of target CDS; it was removed during predictor postprocessing. Display only; not a selectable model.")
    return dict(result)
