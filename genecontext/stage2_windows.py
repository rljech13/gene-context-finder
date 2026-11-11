"""Resolve reannotation bounds from the sealed preliminary gene neighborhood.

Coordinates are native, zero-based and half-open. This plans DNA extraction;
it never relabels an existing prediction as covering a larger interval.
"""
from copy import deepcopy


def resolve_windows(data, regions, decision_rows, policy='pre-synteny'):
    if policy not in {'pre-synteny', 'fixed'}:
        raise ValueError('Unknown stage 2 window policy: '+str(policy))
    loci = {(l['genome_id'], l['locus_id']): l for l in data['loci']}
    decisions = {(r['genome_id'], r['locus_id']): r for r in decision_rows}
    result = []
    for original in regions:
        r = deepcopy(original)
        key = r['genome_id'], r['locus_id']
        locus = loci[key]
        context = locus.get('early_synteny') or {}
        row = decisions[key]
        manual = any(str(row.get(k) or '').strip() for k in ['window_start', 'window_end'])
        detail = dict(policy=policy, previous_start=r['window_start'],
                      previous_end=r['window_end'], neighbor_genes=[], display_only_genes=[])
        if manual:
            detail.update(source='manual', reason='Explicit reviewed window retained')
        elif policy == 'fixed':
            detail.update(source='fixed', reason='Explicit fixed-flank policy')
        elif not context.get('neighbors'):
            detail.update(source='fixed_fallback', reason='No preliminary neighbor coordinates',
                          preliminary_status=context.get('status', 'not_available'))
        else:
            for name, expected in [('genome_id', key[0]), ('locus_id', key[1])]:
                if context.get(name) != expected:
                    raise ValueError('Preliminary context identity differs: '+r['region_id'])
            bounds = [(r['anchor_start'], r['anchor_end']), (locus['start'], locus['end'])]
            a, b = context.get('window_start'), context.get('window_end')
            if type(a) is not int or type(b) is not int or not 0 <= a < b <= locus['contig_length']:
                raise ValueError('Invalid preliminary window: '+r['region_id'])
            bounds.append((a, b))
            # Display-only genes can extend beyond the saved independent-neighbor
            # envelope; include their full span too, rather than cutting them off.
            for field in ['neighbors', 'display_only_genes']:
                genes = context.get(field, [])
                if not isinstance(genes, list):
                    raise ValueError('Invalid preliminary gene list: '+r['region_id'])
                for gene in genes:
                    a, b = gene.get('start'), gene.get('end')
                    if (gene.get('contig') != r['contig'] or type(a) is not int or type(b) is not int
                            or not 0 <= a < b <= locus['contig_length']):
                        raise ValueError('Invalid preliminary gene coordinates: '+r['region_id'])
                    bounds.append((a, b))
                    detail['neighbor_genes' if field == 'neighbors' else 'display_only_genes'].append(
                        dict(gene_id=gene['gene_id'], name=gene.get('name', ''), start=a, end=b,
                             strand=gene.get('strand', ''), side=gene.get('side', '')))
            r['window_start'] = min(a for a, _ in bounds)
            r['window_end'] = max(b for _, b in bounds)
            detail.update(source='pre-synteny', reason='Complete preliminary neighborhood envelope',
                          preliminary_status=context.get('status', ''))
        if not 0 <= r['window_start'] <= r['anchor_start'] < r['anchor_end'] <= r['window_end'] <= locus['contig_length']:
            raise ValueError('Resolved window does not contain the reviewed anchor: '+r['region_id'])
        r.update(window_source=detail['source'], window_resolution=detail)
        result.append(r)
    return result
