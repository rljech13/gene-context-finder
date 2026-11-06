"""Turn local BLAST matches into conservative, auditable candidate loci.

Coordinates are zero-based half-open. Chaining follows query order in the
transcriptional direction; genomic proximity alone never merges two chains.
"""
from collections import defaultdict
import hashlib
import math

from .sequence import union_length


def intersection_length(left, right):
    points = []
    for a, b in left:
        for c, d in right:
            if max(a, c) < min(b, d):
                points.append((max(a, c), min(b, d)))
    return union_length(points)


def oriented(h):
    return (h['start'], h['end']) if h['strand'] == '+' else (-h['end'], -h['start'])


def compatible(a, b, max_intron, max_overlap):
    """Allow small boundary overlaps, but never a query-coordinate reset."""
    x, y = oriented(a), oriented(b)
    if b['qstart'] <= a['qstart'] or b['qend'] <= a['qend'] or y[0] <= x[0] or y[1] <= x[1]:
        return False
    qover = max(0, a['qend'] - b['qstart'])
    if qover > max_overlap * min(a['qend']-a['qstart'], b['qend']-b['qstart']):
        return False
    gover = max(0, x[1] - y[0])
    if gover > max_overlap * min(x[1]-x[0], y[1]-y[0]):
        return False
    return y[0] - x[1] <= max_intron


def summarise(hsps):
    """Query union is exact; identity/score count overlapping query bases once."""
    covered, matches, weight, score = [], 0.0, 0, 0.0
    for h in sorted(hsps, key=lambda h: (-h['bitscore'], h['hsp_id'])):
        interval = [(h['qstart'], h['qend'])]
        fresh = union_length(interval) - intersection_length(interval, covered)
        fraction = fresh / (h['qend']-h['qstart'])
        effective = h['alignment_length'] * fraction
        matches += h['identity'] * effective
        weight += effective
        score += h['bitscore'] * fraction
        covered.extend(interval)
    return dict(coverage=union_length(covered)/hsps[0]['qlen'],
                identity=matches/weight if weight else 0,
                bitscore=score, evalue=min(h['evalue'] for h in hsps))


def chain_hsps(hsps, *, max_intron=200000, max_overlap=0.15, max_locus_span=2000000):
    """Greedy extraction of best collinear DAG paths; every HSP is retained.

    This proposes loci, not validated transcripts. Alternative paths are not
    exhaustively enumerated; conflicting residual chains are retained and flagged.
    """
    groups = defaultdict(list)
    for h in hsps:
        groups[(h['query'], h['contig'], h['strand'])].append(h)
    chains = []
    for key, values in sorted(groups.items()):
        remaining = sorted(values, key=lambda h: (oriented(h)[0], h['qstart'], h['hsp_id']))
        while remaining:
            scores, paths = [], []
            for i, h in enumerate(remaining):
                best, path = h['bitscore'], [i]
                for j in range(i-1, -1, -1):
                    previous = remaining[j]
                    if oriented(h)[0]-oriented(previous)[1] > max_intron:
                        continue  # A long preceding HSP can still overlap this one.
                    if not compatible(previous, h, max_intron, max_overlap):
                        continue
                    first = remaining[paths[j][0]]
                    if oriented(h)[1]-oriented(first)[0] > max_locus_span:
                        continue
                    overlap = max(0, previous['qend']-h['qstart'])
                    gain = h['bitscore'] * (1-overlap/(h['qend']-h['qstart']))
                    gap = max(0, oriented(h)[0]-oriented(previous)[1])
                    candidate = scores[j] + gain - math.log1p(gap)/10
                    if candidate > best:
                        best, path = candidate, paths[j]+[i]
                scores.append(best); paths.append(path)
            winning = max(range(len(remaining)), key=lambda i: scores[i])
            indices = set(paths[winning])
            own = [h for i, h in enumerate(remaining) if i in indices]
            coordinates = ';'.join(h['hsp_id'] for h in own)
            cid = 'chain_' + hashlib.sha256(coordinates.encode()).hexdigest()[:16]
            chains.append(dict(chain_id=cid, query=key[0], contig=key[1], strand=key[2],
                               start=min(h['start'] for h in own), end=max(h['end'] for h in own),
                               blocks=[(h['start'], h['end']) for h in own], hsps=own,
                               **summarise(own)))
            remaining = [h for i, h in enumerate(remaining) if i not in indices]
    return chains


def same_locus(a, b, shared_fraction=0.5):
    if (a['contig'], a['strand']) != (b['contig'], b['strand']):
        return False
    shared = intersection_length(a['blocks'], b['blocks'])
    smaller = min(union_length(a['blocks']), union_length(b['blocks']))
    return bool(smaller and shared / smaller >= shared_fraction)


def group_loci(chains, shared_fraction=0.5):
    """Complete-link clustering avoids a single broad chain bridging two copies."""
    clusters = []
    for chain in sorted(chains, key=lambda c: (-c['bitscore'], c['contig'], c['start'], c['chain_id'])):
        matching = [g for g in clusters if all(same_locus(chain, c, shared_fraction) for c in g)]
        if matching:
            best = max(matching, key=lambda g: max(intersection_length(chain['blocks'], c['blocks']) for c in g))
            best.append(chain)
        else:
            clusters.append([chain])
    return sorted(clusters, key=lambda g: (g[0]['contig'], min(c['start'] for c in g), g[0]['strand']))
