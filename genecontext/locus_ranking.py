"""Deterministic multi-query review ranking; never alters discovery evidence."""
from collections import defaultdict

DEFINITION = {
    'version': 'weighted_source_support_v1',
    'chain_weight': 'bit_score * query_coverage * aligned_identity',
    'review_score': 'maximum source-group weight + mean weight of all supplied source groups',
    'source_group': 'source species; species sharing an identical query sequence are merged; unknown sources use unique query IDs',
    'missing_support': 0,
    'tie_break': 'smaller best HSP E-value, larger raw bit score, contig, start, locus ID',
    'note': 'Heuristic ordering within a genome/search group, not a probability or orthology decision. Synteny and display filters do not affect this rank.',
}


def rank_loci(data):
    """Best bit-score chain/query, best weighted query/source, then max + mean."""
    parent = {}
    labels = {}
    query_keys = {}
    sequence_keys = defaultdict(set)

    def find(key):
        parent.setdefault(key, key)
        if parent[key] != key:
            parent[key] = find(parent[key])
        return parent[key]

    for q in data['queries']:
        species = sorted({m.get('species', '').strip() for m in q['members']} - {''})
        keys = ['species:'+s for s in species] or ['query:'+q['id']]
        for key in keys:
            labels[key] = key.split(':', 1)[1]
            find(key)
        for key in keys[1:]:
            parent[find(key)] = find(keys[0])
        query_keys[q['id']] = keys[0]
        sequence_keys[q['sequence']].update(keys)
    for keys in sequence_keys.values():
        keys = sorted(keys)
        for key in keys[1:]:
            parent[find(key)] = find(keys[0])
    source_labels = defaultdict(set)
    for key in parent:
        source_labels[find(key)].add(labels[key])
    source_ids = sorted(source_labels)
    source_names = {key: ' + '.join(sorted(source_labels[key])) for key in source_ids}
    qsource = {qid: find(key) for qid, key in query_keys.items()}
    qmap = {q['id']: q for q in data['queries']}
    metrics = {}
    groups = defaultdict(list)
    for locus in data['loci']:
        best = {}
        for c in locus['chains']:
            if c.get('ambiguous_support'):
                continue
            old = best.get(c['query'])
            if old is None or (c['bitscore'], c['coverage'], c['identity'], c['chain_id']) > (old['bitscore'], old['coverage'], old['identity'], old['chain_id']):
                best[c['query']] = c
        contributions = {}
        for qid, c in best.items():
            weight = c['bitscore'] * c['coverage'] * c['identity']
            source = qsource[qid]
            item = dict(source=source_names[source], query=qid,
                        query_names=[m['query'] for m in qmap[qid]['members']],
                        chain_id=c['chain_id'], weight=weight, bitscore=c['bitscore'],
                        coverage=c['coverage'], identity=c['identity'], evalue=c['evalue'])
            if source not in contributions or (weight, c['bitscore'], qid) > (contributions[source]['weight'], contributions[source]['bitscore'], contributions[source]['query']):
                contributions[source] = item
        values = [contributions[key]['weight'] if key in contributions else 0 for key in source_ids]
        strongest = max(values, default=0)
        mean = sum(values)/len(values) if values else 0
        entry = dict(genome_id=locus['genome_id'], locus_id=locus['locus_id'],
                     discovery_rank=locus['rank_within_genome'], score=strongest+mean,
                     strongest_source_weight=strongest, mean_source_weight=mean,
                     represented_sources=sum(v > 0 for v in values), supplied_sources=len(values),
                     best_evalue=min((c['evalue'] for c in best.values()), default=float('inf')),
                     best_raw_bitscore=max((c['bitscore'] for c in best.values()), default=0),
                     contributions=[contributions[key] if key in contributions else dict(source=source_names[key],query='',query_names=[],chain_id='',weight=0,bitscore=0,coverage=0,identity=0,evalue=None) for key in source_ids])
        metrics[locus['locus_id']] = entry
        groups[locus['genome_id'], locus.get('gene', data['config']['gene'])].append(locus)
    for group in groups.values():
        group.sort(key=lambda l: (-metrics[l['locus_id']]['score'], metrics[l['locus_id']]['best_evalue'], -metrics[l['locus_id']]['best_raw_bitscore'], l['contig'], l['start'], l['locus_id']))
        for rank, locus in enumerate(group, 1):
            metrics[locus['locus_id']]['rank'] = rank
    # JSON has no infinity: unavailable evidence is represented explicitly.
    for entry in metrics.values():
        if entry['best_evalue'] == float('inf'):
            entry['best_evalue'] = None
    return dict(definition=DEFINITION, source_groups=list(source_names.values()), loci=metrics)
