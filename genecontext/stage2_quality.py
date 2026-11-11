"""Model recovery checks, separate from ORF integrity and biological approval."""
from collections import defaultdict
import copy
import json
from pathlib import Path
import re

from .sequence import attributes, open_text, sha256, union_length


FLAG_KEYS = {'exception', 'partial', 'pseudo', 'pseudogene', 'transl_except'}


def annotation_flags(attrs):
    return sorted({f'{k}={v}' for k, v in attrs.items() if k in FLAG_KEYS})


def discrepancy_flags(flags):
    return [f for f in flags if ('discrepancy' in f.lower() or
            f.startswith(('transl_except=', 'pseudogene=')) or
            f.lower() in {'partial=true', 'pseudo=true'} or
            'ribosomal slippage' in f.lower())]


def shared_blocks(a, b):
    return [(max(x[0], y[0]), min(x[1], y[1]))
            for x in a for y in b if max(x[0], y[0]) < min(x[1], y[1])]


def model_quality(data):
    models = data['models']; groups = defaultdict(list); result = {}
    for m in models:
        groups[m['genome_id'], m['region_id'], m['gene_group']].append(m)
        coverage = [h['coverage_subject'] for h in m.get('target_support', [])
                    if 'coverage_subject' in h]
        result[m['model_id']] = dict(query_coverage=max(coverage) if coverage else None,
            annotation_flags=m.get('annotation_flags', []), warnings=[], notes=[],
            split_models=[], annotation_cds_fraction=None, annotation_exact=False,
            reconstruction=m.get('reconstruction'))
    for members in groups.values():
        annotated = [m for m in members if m['source'] == 'annotation']
        predicted = [m for m in members if m['source'] == 'tiberius']
        for m in members:
            q = result[m['model_id']]; blocks = m.get('cds_blocks', [])
            exact = [a for a in annotated if blocks and a.get('cds_blocks') == blocks
                     and a.get('cds_sha256') == m.get('cds_sha256')]
            q['annotation_exact'] = bool(exact)
            comparators = exact or sorted(annotated, key=lambda a:
                -union_length(shared_blocks(blocks, a.get('cds_blocks', []))))[:1]
            flags = m.get('annotation_flags', []) if m['source'] == 'annotation' else sorted({
                flag for a in comparators for flag in a.get('annotation_flags', [])})
            q['comparison_annotation_flags'] = flags
            if discrepancy_flags(flags):
                q['warnings'].append('Source annotation has a sequence discrepancy or exceptional CDS')
            if m['source'] == 'tiberius' and comparators and blocks:
                a = comparators[0]; size = union_length([b[:2] for b in a.get('cds_blocks', [])])
                fraction = union_length(shared_blocks(blocks, a.get('cds_blocks', []))) / size if size else None
                q['annotation_cds_fraction'] = fraction
                if exact:
                    q['notes'].append('Exact CDS match to one annotation model')
                elif fraction is not None:
                    q['notes'].append(f'Different CDS structure; {fraction:.0%} of comparator CDS bases retained')
                    if fraction < .8:
                        q['warnings'].append('Prediction covers less than 80% of the comparator annotation CDS')
            if (m['source'] == 'tiberius' and not exact and q['query_coverage'] is not None
                    and q['query_coverage'] < .5):
                q['warnings'].append('Best target-query coverage is below 50%; full-gene recovery is unverified')
        # Flag a split only when separate ORFs jointly recover most of ONE
        # annotated CDS. Overlapping isoforms and repeated window observations
        # must not create false splits.
        split = set()
        for a in annotated:
            ab = a.get('cds_blocks', []); size = union_length([b[:2] for b in ab])
            if not size:
                continue
            for i, x in enumerate(predicted):
                xb = x.get('cds_blocks', []); xs = shared_blocks(xb, ab)
                if not xs:
                    continue
                for y in predicted[i+1:]:
                    if (x.get('original_gene') == y.get('original_gene') or
                            x.get('prediction_model_id', x['model_id']) == y.get('prediction_model_id', y['model_id'])):
                        continue
                    yb = y.get('cds_blocks', []); ys = shared_blocks(yb, ab)
                    if not ys:
                        continue
                    lengths = [union_length([b[:2] for b in bs]) for bs in [xb, yb]]
                    if (min(lengths) and union_length(shared_blocks(xb, yb)) / min(lengths) < .2
                            and max(union_length(xs), union_length(ys)) / size < .8
                            and union_length(xs + ys) / size >= .8):
                        split.update([x['model_id'], y['model_id']])
        if split:
            recovered=next((a for a in annotated if (a.get('reconstruction') or {}).get('resolves_boundary_conflict')),None)
            for m in members:
                q = result[m['model_id']]; q['split_models'] = sorted(split)
                if recovered:
                    q['full_reconstruction_model']=recovered['model_id']
                if m['model_id'] in split:
                    q['warnings'].insert(0, 'Fragmentary prediction relative to the reconstructed full CDS' if recovered else 'Possible split prediction: separate ORFs recover one annotated CDS')
                elif m['source'] == 'annotation':
                    rec=m.get('reconstruction') or {}
                    if rec.get('resolves_boundary_conflict'):
                        q['notes'].insert(0, f"Full CDS reconstructed by exact splice alignment of {rec['accession']}; supplied RefSeq reports RNA-seq support for all introns. Separate Tiberius ORFs are fragmentary alternatives.")
                    else:
                        q['notes'].insert(0, f'One annotated transcript conflicts with {len(split)} separate Tiberius ORFs')
                        q['warnings'].insert(0, 'Conflicting gene boundaries: annotation joins CDS segments assigned to separate predicted ORFs; exon connectivity needs review')
    return result


def collect_annotation_quality(data, manifest, evidence_sha256):
    from .early_synteny import datasets
    wanted = defaultdict(lambda: defaultdict(set))
    for m in data['models']:
        if m['source'] == 'annotation':
            wanted[m['genome_id']][m['original_model']].add(m['model_id'])
            wanted[m['genome_id']][m['original_gene']].add(m['model_id'])
    expected = {r['genome_id']: r.get('annotation_sha256') for r in data['provenance']['genomes']}
    sources = []; flags = {m['model_id']: set() for m in data['models'] if m['source'] == 'annotation'}
    for entry in datasets(manifest):
        gid = entry['genome_id']
        if gid not in wanted:
            continue
        path = Path(entry['gff3']); digest = sha256(path)
        if not expected.get(gid) or digest != expected[gid]:
            raise ValueError('Annotation fingerprint differs for ' + gid)
        pattern = re.compile('|'.join(re.escape(k) + r'(?:;|,|$)' for k in wanted[gid]))
        with open_text(path) as handle:
            for line in handle:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#') or not pattern.search(line):
                    continue
                fields = line.rstrip().split('\t')
                if len(fields) != 9:
                    continue
                attrs = attributes(fields[8]); ids = [attrs.get('ID', '')] + attrs.get('Parent', '').split(',')
                for key in ids:
                    for mid in wanted[gid].get(key, []):
                        flags[mid].update(annotation_flags(attrs))
        sources.append(dict(genome_id=gid, annotation_sha256=digest))
    missing = set(wanted) - {s['genome_id'] for s in sources}
    if missing:
        raise ValueError('Missing annotation inputs: ' + ', '.join(sorted(missing)))
    return dict(schema='stage2_annotation_quality_v1', source_evidence_sha256=evidence_sha256,
                sources=sources, models={mid: sorted(values) for mid, values in flags.items()})


def with_annotation_quality(root, data, evidence_name='evidence.json', supplement_name='annotation_quality.json'):
    path = Path(root) / supplement_name
    if not path.exists():
        return data
    context = json.loads(path.read_text())
    if (context.get('schema') != 'stage2_annotation_quality_v1' or
            context.get('source_evidence_sha256') != sha256(Path(root) / evidence_name)):
        raise ValueError('Annotation quality metadata belongs to another checkpoint')
    expected = {g['genome_id']: g.get('annotation_sha256') for g in data['provenance']['genomes']}
    for source in context['sources']:
        if expected.get(source['genome_id']) != source['annotation_sha256']:
            raise ValueError('Annotation quality source fingerprint differs')
    lookup = {m['model_id']: m for m in data['models'] if m['source'] == 'annotation'}
    if set(context['models']) - set(lookup):
        raise ValueError('Unknown annotation model in quality metadata')
    data = copy.deepcopy(data)
    for m in data['models']:
        if m['model_id'] in context['models']:
            values = context['models'][m['model_id']]
            if not isinstance(values, list) or any(not isinstance(v, str) for v in values):
                raise ValueError('Invalid annotation quality flags')
            m['annotation_flags'] = values
    return data


def run(args):
    from .stage2_choices import verify_checkpoint
    root, data = verify_checkpoint(args.results)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
    context = collect_annotation_quality(data, args.genomes, sha256(root / 'evidence.json'))
    (out / 'annotation_quality.json').write_text(json.dumps(context, indent=2) + '\n')
    print(str(out / 'annotation_quality.json'))
    return 0
