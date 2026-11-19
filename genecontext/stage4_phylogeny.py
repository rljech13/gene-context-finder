"""Small per-group ML trees. Root assumptions never become identity evidence."""
from collections import defaultdict
from html import escape
import hashlib
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

from .sequence import read_fasta, write_fasta

ALPHABET = set('ACDEFGHIKLMNPQRSTVWY')


def unique_alignment(alignment, labels, occupancy=.5):
    """Collapse full-protein duplicates; retain all provenance as leaf aliases."""
    if not alignment or len({len(s) for s in alignment.values()}) != 1:
        raise ValueError('A nonempty, equal-length protein alignment is required')
    aliases = defaultdict(list)
    for sid, seq in alignment.items():
        aliases[hashlib.sha256(seq.replace('-', '').encode()).hexdigest()].append(sid)
    sequences, leaves = {}, []
    for n, (digest, members) in enumerate(aliases.items(), 1):
        tid = 'T%05d' % n
        sequences[tid] = alignment[members[0]]
        leaves.append(dict(tree_id=tid, protein_sha256=digest, sequence_ids=members,
                           labels=[labels[s] for s in members]))
    width = len(next(iter(sequences.values())))
    keep = [i for i in range(width) if sum(s[i] in ALPHABET for s in sequences.values()) / len(sequences) >= occupancy]
    masked = {k: ''.join(s[i] if s[i] in ALPHABET else '-' for i in keep) for k, s in sequences.items()}
    # A fragment must not become an all-gap taxon after masking.
    excluded = [k for k, s in masked.items() if not set(s) & ALPHABET]
    masked = {k: s for k, s in masked.items() if k not in excluded}
    for leaf in leaves:
        leaf['included'] = leaf['tree_id'] not in excluded
    informative = sum(len({s[i] for s in masked.values()} - {'-'}) > 1 for i in range(len(keep)))
    return masked, leaves, dict(original_columns=width, retained_columns=len(keep),
        retained_columns_1based=[i+1 for i in keep], variable_columns=informative,
        occupancy=occupancy, excluded_all_gap_taxa=excluded)


def parse_newick(text):
    tokens = re.findall(r"\(|\)|,|:|;|[^\s(),:;]+", re.sub(r'\[[^\]]*\]', '', text))
    pos = 0
    def node():
        nonlocal pos
        item = dict(children=[], name='', length=0.)
        if pos < len(tokens) and tokens[pos] == '(':
            pos += 1
            item['children'].append(node())
            while pos < len(tokens) and tokens[pos] == ',':
                pos += 1; item['children'].append(node())
            if pos >= len(tokens) or tokens[pos] != ')': raise ValueError('Malformed Newick')
            pos += 1
        if pos < len(tokens) and tokens[pos] not in [':', ',', ')', ';']:
            item['name'] = tokens[pos]; pos += 1
        if pos < len(tokens) and tokens[pos] == ':':
            pos += 1; item['length'] = float(tokens[pos]); pos += 1
            if not math.isfinite(item['length']) or item['length'] < 0: raise ValueError('Invalid tree branch length')
        return item
    tree = node()
    if tokens[pos:] != [';']: raise ValueError('Malformed Newick terminator')
    return tree


def tip_names(node):
    return sum((tip_names(c) for c in node['children']), []) if node['children'] else [node['name']]


def tree_svg(newick, leaves):
    tree = parse_newick(newick)
    by_id = {l['tree_id']: l for l in leaves if l['included']}
    if set(tip_names(tree)) != set(by_id) or len(tip_names(tree)) != len(by_id):
        raise ValueError('Tree leaves differ from the exact deduplicated alignment')
    tips = []
    def place(n, distance=0):
        n['x'] = distance
        if not n['children']:
            n['y'] = 24 + len(tips)*25; tips.append(n)
        else:
            for c in n['children']: place(c, distance+c['length'])
            n['y'] = sum(c['y'] for c in n['children'])/len(n['children'])
    place(tree)
    maximum = max(n['x'] for n in tips)
    scale = 510/max(maximum, .001)
    height = 65+len(tips)*25
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1120 {height}" role="img" aria-label="Protein gene tree"><rect width="100%" height="100%" fill="white"/><g font-family="system-ui" font-size="12" fill="#27464a">']
    def draw(n):
        x = 18+n['x']*scale; y = n['y']
        if n['children']:
            ys = [c['y'] for c in n['children']]
            out.append(f'<path d="M{x},{min(ys)} V{max(ys)}" stroke="#60868b" fill="none"/>')
            for c in n['children']:
                out.append(f'<path d="M{x},{c["y"]} H{18+c["x"]*scale}" stroke="#60868b" fill="none"/>'); draw(c)
            if n['name']: out.append(f'<text x="{x+3}" y="{y-5}" font-size="9">{escape(n["name"])}</text>')
        else:
            leaf = by_id[n['name']]; label = leaf['labels'][0]
            suffix = f' (+{len(leaf["labels"])-1} identical)' if len(leaf['labels']) > 1 else ''
            out.append(f'<text x="{x+7}" y="{y+4}"><title>{escape("; ".join(leaf["labels"]))}</title>{escape(label+suffix)}</text>')
    draw(tree)
    bar = .1 if maximum >= .1 else .01
    out.append(f'<path d="M18,{height-24} h{bar*scale}" stroke="#27464a"/><text x="18" y="{height-7}">{bar:g} substitutions / site</text></g></svg>')
    return ''.join(out)


def build_tree(folder, alignment, labels, runner, threads=2, bootstrap=1000,
               occupancy=.5, outgroups=(), root_basis='', engine=None):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    masked, leaves, mask = unique_alignment(alignment, labels, occupancy)
    write_fasta(folder/'tree.aligned.faa', masked.items())
    result = dict(status='not_computed', leaves=leaves, mask=mask, rooting='unrooted',
        root_basis='', outgroup_sequence_ids=list(outgroups), bootstrap_replicates=bootstrap,
        method='Protein maximum likelihood; ModelFinder BIC; UFBoot; seed 42', commands=[])
    if outgroups:
        if not root_basis.strip(): raise ValueError('Outgroup rooting requires an explicit biological basis')
        unknown = set(outgroups)-set(alignment)
        if unknown: raise ValueError('Unknown outgroup sequence: '+', '.join(sorted(unknown)))
        outer = [l['tree_id'] for l in leaves if l['included'] and set(l['sequence_ids']) & set(outgroups)]
        if not outer or len(outer) == len(masked): raise ValueError('Outgroups must separate from the ingroup')
        if any(set(l['sequence_ids']) & set(outgroups) and set(l['sequence_ids'])-set(outgroups) for l in leaves):
            raise ValueError('An identical protein cannot be both ingroup and outgroup')
        result.update(rooting='outgroup', root_basis=root_basis, outgroup_tree_ids=outer)
    if len(masked) < 4 or mask['retained_columns'] < 30 or mask['variable_columns'] < 5:
        result['reason'] = 'Insufficient taxa or variable aligned sites for this prototype (4 taxa, 30 retained and 5 variable columns required).'
        return result
    engine = engine or next((shutil.which(n) for n in ['iqtree3','iqtree2','iqtree'] if shutil.which(n)), None)
    if not engine:
        result['reason'] = 'IQ-TREE is not installed; alignment and leaf provenance are exported.'
        return result
    argv = [engine, '-s', folder/'tree.aligned.faa', '-st', 'AA', '-m', 'MFP', '-B', str(bootstrap),
            '-bnni', '-seed', '42', '-T', str(threads), '--prefix', folder/'ml']
    if outgroups: argv += ['-o', ','.join(outer)]
    started = time.perf_counter()
    runner.run(argv, folder/'stdout.log', folder/'stderr.log')
    newick = (folder/'ml.treefile').read_text()
    (folder/'tree.svg').write_text(tree_svg(newick, leaves))
    result.update(status='computed', elapsed_seconds=time.perf_counter()-started,
                  tool_path=str(engine), tool_version=subprocess.run([engine,'--version'],capture_output=True,text=True).stdout.splitlines()[:3])
    # A tree containing no named competitors cannot distinguish unsampled paralogs.
    result['identity_use'] = 'Context only; branch support does not establish a target name or orthology.'
    if not outgroups: result['root_basis'] = 'No validated outgroup profile supplied; display orientation is arbitrary.'
    return result
