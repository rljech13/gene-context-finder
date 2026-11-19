"""Fresh gene-count neighborhoods from the matching complete GFF, before DNA work."""
from collections import defaultdict
from pathlib import Path

from .early_synteny import GeneIndex, datasets, parse_annotation
from .sequence import sha256


def annotation_neighborhoods(data, regions, manifest, requested=10, minimum=5):
    if not manifest:
        raise ValueError('Annotation-neighbor planning requires --genomes with matching GFF3 files')
    if not 1 <= minimum <= requested:
        raise ValueError('Neighbor minimum must be positive and no greater than the requested count')
    inputs = {g['genome_id']: g for g in datasets(manifest)}
    genomes = {g['genome_id']: g for g in data['genomes']}
    loci = {(l['genome_id'], l['locus_id']): l for l in data['loci']}
    grouped = defaultdict(list)
    for r in regions:
        grouped[r['genome_id']].append(r)
    excluded = defaultdict(set)
    for ref in data.get('early_synteny', {}).get('references', []):
        if ref.get('gene_id'):
            excluded[ref['genome_id']].add(ref['gene_id'])
    contexts = {}
    for gid, own in grouped.items():
        entry = inputs.get(gid)
        if entry is None:
            raise ValueError('Selected genome absent from --genomes: ' + gid)
        path = entry.get('gff3')
        if not path:
            for r in own:
                contexts[gid, r['locus_id']] = dict(genome_id=gid, locus_id=r['locus_id'],
                    status='annotation_unavailable', neighbors=[], display_only_genes=[],
                    neighbors_requested=requested, minimum_neighbors=minimum,
                    limits=dict(left='annotation_unavailable', right='annotation_unavailable'))
            continue
        path = Path(path); before = path.stat(); digest = sha256(path)
        expected = genomes[gid].get('annotation_inputs', {}).get('gff3', {}).get('sha256')
        if expected and digest != expected:
            raise ValueError('Annotation differs from stage 1 dataset: ' + gid)
        genes, models, _ = parse_annotation(path)
        del models
        index = GeneIndex(dict(genes=list(genes.values())))
        for r in own:
            locus = loci[gid, r['locus_id']]
            a = min(locus['start'], r['anchor_start'])
            b = max(locus['end'], r['anchor_end'])
            context = index.neighbors(r['contig'], a, b, r['strand'], requested, 0,
                exclude=excluded[gid], anchor_gene=locus.get('annotation_gene_id'))
            # Preserve the POU contract: existing displayed preliminary genes
            # remain complete when refreshing a previously sparse neighborhood.
            old = locus.get('early_synteny') or {}
            old_genes = old.get('neighbors', []) + old.get('display_only_genes', [])
            retained = []
            for prior in old_genes:
                gene = genes.get(prior.get('gene_id'))
                if not gene or gene['contig'] != r['contig']:
                    raise ValueError('Preliminary gene absent from matching full annotation: ' + gid)
                if (gene['start'],gene['end']) != (prior['start'],prior['end']):
                    raise ValueError('Preliminary/full annotation gene bounds differ: ' + gid)
                retained.append(dict(gene))
            context['retained_preliminary_genes'] = retained
            if retained:
                context['window_start'] = min(context['window_start'], *(g['start'] for g in retained))
                context['window_end'] = max(context['window_end'], *(g['end'] for g in retained))
            for gene in context['neighbors'] + context['display_only_genes']:
                if not 0 <= gene['start'] < gene['end'] <= locus['contig_length']:
                    raise ValueError('Annotation outside the saved contig bounds: ' + gid)
            context.update(genome_id=gid, locus_id=r['locus_id'], status='full_annotation_count_first',
                gff3_sha256=digest, annotation_identity='exact_stage1_sha256' if expected else 'new_annotation_fingerprint',
                neighbors_requested=requested, minimum_neighbors=minimum,
                distance_cap=0, annotation_coding_genes=sum(g['coding'] for g in genes.values()))
            contexts[gid, r['locus_id']] = context
        if (path.stat().st_size, path.stat().st_mtime_ns) != (before.st_size, before.st_mtime_ns):
            raise ValueError('GFF3 changed during neighborhood planning: ' + gid)
    return contexts


def coverage(context, start, end):
    requested = context['neighbors_requested']; minimum = context['minimum_neighbors']
    included = [g for g in context['neighbors'] if start <= g['start'] < g['end'] <= end]
    sides = {side: sum(g['side'] == side for g in included) for side in ['left', 'right']}
    missing = [g['gene_id'] for g in context['neighbors'] if g not in included]
    ready = not missing and all(n >= minimum for n in sides.values())
    return dict(requested_per_side=requested, minimum_per_side=minimum, captured=sides,
        fully_captured_neighbor_ids=[g['gene_id'] for g in included],
        neighbors_outside_reviewed_window=missing, limits=context.get('limits', {}),
        source_status=context['status'],
        status='sufficient_context' if ready else context['status'] if context['status']=='annotation_unavailable' else 'insufficient_context',
        ready_for_prediction=ready,
        note='Distinct coding genes, not transcript count; overlapping/query-source copies do not fill the quota. RNA context is retained separately.')
