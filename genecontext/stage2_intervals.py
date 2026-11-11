"""Bounded shared prediction windows and projection onto reviewed loci."""
from collections import defaultdict
from dataclasses import replace
import hashlib
import json


def canonical_windows(regions, maximum):
    buckets = defaultdict(list)
    for r in regions:
        buckets[r['genome_id'], r['contig']].append(r)
    result = []
    for (gid, contig), rows in sorted(buckets.items()):
        clusters = []
        for r in sorted(rows, key=lambda x: (x['window_start'], x['window_end'], x['region_id'])):
            if clusters and r['window_start'] < clusters[-1]['window_end']:
                c = clusters[-1]
                end = max(c['window_end'], r['window_end'])
                if end - c['window_start'] > maximum:
                    raise ValueError('Overlapping prediction windows exceed --max-window-bases: '
                                     + gid + ' ' + contig + '; review their boundaries')
                c['window_end'] = end
                c['region_ids'].append(r['region_id'])
            else:
                if r['window_end'] - r['window_start'] > maximum:
                    raise ValueError('Prediction window exceeds --max-window-bases')
                clusters.append(dict(genome_id=gid, contig=contig,
                    window_start=r['window_start'], window_end=r['window_end'],
                    region_ids=[r['region_id']], gene=r.get('gene', 'target')))
        for c in clusters:
            key = [gid, contig, c['window_start'], c['window_end']]
            cid = 'prediction_' + hashlib.sha256(json.dumps(key).encode()).hexdigest()[:24]
            c.update(region_id=cid, locus_id=cid, strand='+',
                     anchor_start=c['window_start'], anchor_end=c['window_end'])
            result.append(c)
    return result


def project_model(model, canonical, region):
    """Keep full CDS, including bases outside a smaller original display window."""
    start = canonical['window_start'] + model.start
    end = canonical['window_start'] + model.end
    if max(start, region['window_start']) >= min(end, region['window_end']):
        return None
    delta = canonical['window_start'] - region['window_start']
    return replace(model, contig=region['region_id'],
                   cds=[(a + delta, b + delta, phase) for a, b, phase in model.cds],
                   stops=[(a + delta, b + delta) for a, b in model.stops])
