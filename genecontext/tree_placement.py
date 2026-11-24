"""Automatic, root-invariant placement against an independent species scaffold.

T/15 is a descriptive support-weighted quartet agreement, not an orthology
probability. No leaf, clade, branch or outgroup is chosen to assign this score.
This module is stdlib-only so the exact same bytes can ship in a review kit.
"""
from collections import Counter, defaultdict
import csv
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import random
import re

METHOD = 'species_quartet_placement_v1'
MAX_QUARTETS = 2000
LIMITATION = ('Agreement with the supplied species scaffold, not a probability of correct annotation or orthology. '
              'Gene duplication/loss, incomplete lineage sorting, alignment error and genuine gene-history discordance can lower agreement. '
              'Quartets share data and are not independent replicates. Inspect coverage and support separately.')


def digest(value):
    raw = value if isinstance(value, bytes) else value.encode()
    return hashlib.sha256(raw).hexdigest()


def taxon_name(value):
    return ' '.join(str(value or '').replace('_', ' ').split())


def read_lineages(path=None):
    path = Path(path) if path else Path(__file__).with_name('tree_placement_taxonomy.tsv')
    with path.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {}
    for row in rows:
        name = taxon_name(row.get('taxon') or row.get('species'))
        chain = tuple(x.strip() for x in row['lineage'].split('/') if x.strip())
        if not name or not chain or name in result:
            raise ValueError('Lineages need unique taxa and nonempty root-to-tip paths')
        result[name] = chain
    return result


class Tree:
    """Canonical unrooted edge splits; degree-two rooting is immaterial."""
    def __init__(self, text, *, support_field=None, normalize=False):
        if len(text) > 4 * 1024 * 1024:
            raise ValueError('Tree exceeds 4 MiB')
        tokens = re.findall(r"'(?:[^']|'')*'|[(),:;]|[^\s(),:;]+", re.sub(r'\[[^\]]*\]', '', text))
        self.pos = 0
        names, edges = [], []
        def take():
            if self.pos >= len(tokens):
                raise ValueError('Truncated Newick')
            value = tokens[self.pos]; self.pos += 1
            return value
        def branch(depth=0):
            if depth > 500:
                raise ValueError('Tree nesting exceeds 500')
            children = []
            if self.pos < len(tokens) and tokens[self.pos] == '(':
                take(); children.append(branch(depth + 1))
                while self.pos < len(tokens) and tokens[self.pos] == ',':
                    take(); children.append(branch(depth + 1))
                if take() != ')' or len(children) < 2:
                    raise ValueError('Invalid Newick branching')
            label = ''
            if self.pos < len(tokens) and tokens[self.pos] not in {':', ',', ')', ';'}:
                label = take()
                if label.startswith("'"):
                    label = label[1:-1].replace("''", "'")
            length = None
            if self.pos < len(tokens) and tokens[self.pos] == ':':
                take()
                try: length = float(take())
                except ValueError as error: raise ValueError('Invalid branch length') from error
                if not math.isfinite(length) or length < 0:
                    raise ValueError('Invalid branch length')
            if not children:
                label = taxon_name(label) if normalize else label
                if not label or label in names:
                    raise ValueError('Tree needs unique, nonempty tip names')
                names.append(label)
                return {label}, label, length
            tips = set().union(*(c[0] for c in children))
            for child_tips, child_label, child_length in children:
                edges.append((child_tips, child_label if len(child_tips) > 1 else '', child_length))
            return tips, label, length
        all_tips, _, _ = branch()
        if take() != ';' or self.pos != len(tokens):
            raise ValueError('Expected one terminated Newick tree')
        self.tips = tuple(sorted(all_tips)); self.bits = {n: 1 << i for i, n in enumerate(self.tips)}
        complete = (1 << len(self.tips)) - 1
        buckets = defaultdict(list)
        for side, label, length in edges:
            mask = sum(self.bits[n] for n in side); other = complete ^ mask
            if mask.bit_count() < 2 or other.bit_count() < 2:
                continue
            key = min(mask, other)
            support = None
            try:
                fields = label.split('/')
                if len(fields) == 1 or support_field is not None:
                    support = float(fields[0 if support_field is None else support_field])
                    if not math.isfinite(support) or not 0 <= support <= 100:
                        support = None
            except (ValueError, IndexError):
                pass
            buckets[key].append((support, length))
        self.edges = []
        for mask, values in buckets.items():
            # Newick may split a single unrooted edge at its display root.
            known = {s for s, _ in values if s is not None}
            support = next(iter(known)) if len(known) == 1 else None
            lengths = [length for _, length in values]
            if all(length is not None for length in lengths) and sum(lengths) == 0:
                continue
            self.edges.append((mask, support))
        self.cache = {}

    def quartet(self, tips, *, support=False):
        key = (tuple(tips), support)
        if key in self.cache:
            return self.cache[key]
        if len(set(tips)) != 4:
            return None, None, 'identical_sequence_tie'
        if any(t not in self.bits for t in tips):
            return None, None, 'taxon_not_in_tree'
        bits = [self.bits[t] for t in tips]; whole = sum(bits); resolutions = set(); supports = []
        for mask, value in self.edges:
            part = mask & whole
            if part.bit_count() != 2:
                continue
            own = part if part & bits[0] else whole ^ part
            resolutions.add(next(i for i in (1, 2, 3) if own & bits[i]))
            supports.append(value)
        if len(resolutions) != 1:
            result = (None, None, 'unresolved_topology')
        elif support and any(v is None for v in supports):
            result = (next(iter(resolutions)), None, 'support_unavailable')
        else:
            result = (next(iter(resolutions)), min(supports) if support else None, '')
        if len(self.cache) < 200000:
            self.cache[key] = result
        return result


def lineage_newick(lineages):
    root = {}
    for name, chain in sorted(lineages.items()):
        node = root
        for part in chain:
            node = node.setdefault(part, {})
        node['@' + name] = None
    def quote(name): return "'" + name.replace("'", "''") + "'"
    def emit(node):
        parts = [quote(key[1:]) if key.startswith('@') else emit(value) for key, value in sorted(node.items())]
        return parts[0] if len(parts) == 1 else '(' + ','.join(parts) + ')'
    return emit(root) + ';'


def triples(taxa, focal):
    """Exhaustive for ordinary panels; deterministic bounded sampling for big ones."""
    total = math.comb(len(taxa), 3)
    if total <= MAX_QUARTETS:
        return list(combinations(taxa, 3)), total
    rng = random.Random(METHOD + '|' + focal + '|' + json.dumps(taxa))
    ranks = sorted(rng.sample(range(total), MAX_QUARTETS)); rows = []
    n = len(taxa)
    for rank in ranks:
        i = 0
        while rank >= math.comb(n - i - 1, 2):
            rank -= math.comb(n - i - 1, 2); i += 1
        j = i + 1
        while rank >= n - j - 1:
            rank -= n - j - 1; j += 1
        rows.append((taxa[i], taxa[j], taxa[j + 1 + rank]))
    return rows, total


def assess(newick, models, *, species_tree=None, lineages=None, support_kind='reported_branch_support', support_field=None, details=False):
    """Assess each exact model automatically against one independently supplied scaffold.

    models: model_id, tree_id, taxon, optional protein_sha256 and primary (peer
    eligibility). References are not peers. More than one distinct peer protein
    for a species makes that peer species ambiguous; no best copy is selected.
    """
    if not isinstance(models, (list, tuple)) or any(not isinstance(m, dict) for m in models):
        raise ValueError('Models must be a list of records')
    if support_field is not None and (type(support_field) is not int or support_field < 0):
        raise ValueError('Support field must be a nonnegative integer')
    if lineages is not None and not isinstance(lineages, dict):
        raise ValueError('Lineages must map taxa to root-to-tip lists')
    if species_tree is not None and lineages is not None:
        raise ValueError('Supply a species tree OR taxonomic lineages')
    if lineages is None and species_tree is None:
        lineages = read_lineages()
    normalized = {}
    for key, chain in (lineages or {}).items():
        name = taxon_name(key)
        if not name or name in normalized or not isinstance(chain, (list, tuple)) or not chain or any(not isinstance(x, str) or not x.strip() or x.startswith('@') for x in chain):
            raise ValueError('Lineages need unique normalized taxa and nonempty lists of clade names')
        normalized[name] = tuple(chain)
    reference_text = species_tree if species_tree is not None else lineage_newick(normalized)
    expected = Tree(reference_text, normalize=True); observed = Tree(newick, support_field=support_field)
    records = []; seen = set(); hashes = defaultdict(set)
    for item in models:
        m = dict(item); mid = m.get('model_id') or m.get('id')
        if not isinstance(mid, str) or not mid or mid in seen:
            raise ValueError('Models need distinct nonempty model_id values')
        seen.add(mid); m['model_id'] = mid; m['taxon'] = taxon_name(m.get('taxon') or m.get('species'))
        leaf = m.get('tree_id')
        if leaf and m.get('protein_sha256'):
            hashes[leaf].add(m['protein_sha256'])
        records.append(m)
    if any(len(values) > 1 for values in hashes.values()):
        raise ValueError('One collapsed tree leaf has conflicting protein identities')
    peers = defaultdict(set)
    for m in records:
        if not m.get('reference') and m.get('primary', True) and m['taxon'] in expected.bits and m.get('tree_id') in observed.bits:
            peers[m['taxon']].add(m['tree_id'])
    ambiguous = sorted(t for t, ids in peers.items() if len(ids) > 1)
    unique = {t: next(iter(ids)) for t, ids in peers.items() if len(ids) == 1}
    results = {}
    for m in records:
        mid, taxon, leaf = m['model_id'], m['taxon'], m.get('tree_id')
        ambiguous_other = [t for t in ambiguous if t != taxon]
        result = dict(model_id=mid, taxon=taxon, tree_id=leaf, protein_sha256=m.get('protein_sha256'),
                      method=METHOD, maximum=15, points=None, status='not_assessable', reason='',
                      tree_sha256=digest(newick), reference_sha256=digest(reference_text),
                      support_kind=support_kind, support_field=support_field, counts={}, examples=[])
        results[mid] = result
        if m.get('reference'):
            result['reason'] = 'reference_sequence_not_scored'; continue
        if leaf not in observed.bits:
            result['reason'] = 'model_not_in_tree'; continue
        if not taxon or taxon not in expected.bits:
            result['reason'] = 'species_not_in_reference_scaffold'; continue
        others = sorted(t for t in unique if t != taxon)
        if len(others) < 3:
            result['reason'] = 'fewer_than_three_unambiguous_other_species'; continue
        tests, total = triples(others, taxon); counts = Counter(); weighted = 0.; support_sum = 0.
        examples = []; all_details = []
        for triple in tests:
            taxa = (taxon,) + triple
            want, _, problem = expected.quartet(taxa)
            if problem:
                counts['reference_unresolved'] += 1; continue
            counts['reference_resolved'] += 1
            tips = (leaf,) + tuple(unique[t] for t in triple)
            actual, strength, problem = observed.quartet(tips, support=True)
            if problem:
                counts[problem] += 1; continue
            counts['assessed'] += 1; support_sum += strength
            match = actual == want; counts['concordant' if match else 'discordant'] += 1
            if match: weighted += strength / 100
            row = dict(taxa=list(taxa), expected_partner=taxa[want], observed_partner=taxa[actual], support=strength, concordant=match)
            if details: all_details.append(row)
            if len(examples) < 12: examples.append(row)
        result.update(counts=dict(counts), possible_quartets=total, sampled_quartets=len(tests),
                      sampling='exhaustive' if total == len(tests) else 'deterministic_2000',
                      ambiguous_peer_species=ambiguous_other, examples=examples)
        if details: result['quartets'] = all_details
        assessed = counts['assessed']; eligible = counts['reference_resolved']
        result['coverage'] = assessed / eligible if eligible else None
        if not assessed:
            result['reason'] = 'no_resolved_supported_species_quartets'; continue
        result.update(points=15 * weighted / assessed, weighted_agreement=weighted / assessed,
                      topology_agreement=counts['concordant'] / assessed,
                      mean_support=support_sum / assessed,
                      status='assessed' if assessed == eligible and total == len(tests) and not ambiguous_other else 'partial',
                      reason='')
    for result in results.values():
        result['proposed_points'] = result['points']
        result['included_in_score'] = False
    return dict(schema='gcf_tree_placement_v1', method=METHOD, maximum=15, requires_stage4_tree_review=True,
                formula='15 * sum(minimum central-path support / 100 for concordant assessed quartets) / number of assessed quartets',
                reference_kind='supplied_species_tree' if species_tree is not None else 'taxonomic_scaffold',
                reference_newick=reference_text, reference_sha256=digest(reference_text), tree_sha256=digest(newick),
                limitation=LIMITATION, models=results)


def records_from_members(members, phylogeny, decisions=None):
    alias = {mid: leaf['tree_id'] for leaf in phylogeny.get('leaves', []) if leaf.get('included', True) for mid in leaf['sequence_ids']}
    result = []
    for mid, item in members.items():
        meta = item.get('metadata', {})
        taxon = item.get('species') or meta.get('species') or item.get('label', '').split(' · ')[0]
        if not taxon:
            taxon = re.split(r'_GC[AF]?_', meta.get('genome_id', ''))[0]
        role = (decisions or {}).get(mid, {}).get('role')
        if decisions is not None and role not in {'primary_target', 'special_case'} and not item.get('reference'):
            continue
        result.append(dict(model_id=mid, tree_id=alias.get(mid), taxon=taxon_name(taxon),
                           protein_sha256=meta.get('protein_sha256') or digest(item.get('protein', '')),
                           reference=bool(item.get('reference')), primary=role in {None, 'primary_target'}))
    return result


def reference_options(args):
    tree = getattr(args, 'species_tree', None)
    lineage = getattr(args, 'species_lineages', None)
    if tree and lineage:
        raise ValueError('Use --species-tree OR --species-lineages')
    return {'species_tree': Path(tree).read_text()} if tree else {'lineages': read_lineages(lineage)}


def attach_stage4(root, trees, rows, args):
    """Attach the common assessment to either Stage4 route, without rebuilding trees."""
    root = Path(root); options = reference_options(args)
    for tree in trees:
        group_rows = [r for r in rows if r['group'] == tree['group']]
        alias = {mid: leaf['tree_id'] for leaf in tree.get('leaves', []) if leaf.get('included', True) for mid in leaf['sequence_ids']}
        records = [dict(model_id=r['model_id'], tree_id=alias.get(r['model_id']), taxon=r['species'],
                        protein_sha256=r['protein_sha256'], primary=r['classification'] == 'primary_target') for r in group_rows]
        if tree['status'] == 'computed':
            folder = root / tree['folder']
            newick = tree.get('unrooted_newick') or (folder / 'ml.treefile').read_text()
            result = assess(newick, records, support_kind='UFBoot', **options)
            path = folder / 'placement.json'; path.write_text(json.dumps(result, indent=2)+'\n')
            tree['placement'] = dict(method=METHOD, file=str(path.relative_to(root)),
                                     reference_sha256=result['reference_sha256'], tree_sha256=result['tree_sha256'])
        else:
            result = {'models': {r['model_id']: dict(model_id=r['model_id'], method=METHOD, status='not_assessable',
                      maximum=15, points=None, reason='tree_not_computed', counts={}) for r in group_rows}}
        for row in group_rows:
            value = result['models'][row['model_id']]
            row.update(tree_placement=value, tree_proposed_score=value['points'], tree_score=None, tree_reviewed=False,
                       tree_assessment_status='awaiting_tree_review' if value['points'] is not None else value['status'],
                       tree_agreement=value.get('weighted_agreement'), tree_assessed_quartets=value['counts'].get('assessed', 0),
                       tree_quartet_coverage=value.get('coverage'))


    from .tree_review import write_bundle
    panels = {'analysis': dict(label='Gene Context Finder · tree review', seal='', groups={})}
    for tree in trees:
        if tree['status'] != 'computed':
            continue
        folder = root / tree['folder']
        aligned = {}; current = None
        for line in (folder/'proteins.aligned.faa').read_text().splitlines():
            if line.startswith('>'):
                current = line[1:].split()[0]; aligned[current] = ''
            elif current:
                aligned[current] += line.strip()
        review_rows = {}
        for leaf in tree['leaves']:
            present = [aligned[mid] for mid in leaf['sequence_ids'] if mid in aligned]
            if not present or len(set(present)) != 1:
                raise ValueError('Review leaf aliases differ from the saved alignment')
            review_rows[leaf['tree_id']] = present[0]
        full = json.loads((root/tree['placement']['file']).read_text())
        aliases = {mid: leaf['tree_id'] for leaf in tree['leaves'] for mid in leaf['sequence_ids']}
        records = [dict(model_id=r['model_id'], tree_id=aliases.get(r['model_id']), taxon=r['species'],
                        protein_sha256=r['protein_sha256'], primary=r['classification']=='primary_target',
                        label=r['species']+' · '+r['original_model']) for r in rows if r['group']==tree['group']]
        newick = tree.get('unrooted_newick') or (folder/'ml.treefile').read_text()
        group = dict(label=tree['gene'], rows=review_rows, leaves=tree['leaves'], models=records,
                     trees={'saved':newick}, species_tree=full['reference_newick'], reference_kind=full['reference_kind'],
                     mask=[i-1 for i in tree['mask']['retained_columns_1based']] if tree.get('mask',{}).get('retained_columns_1based') else {'mode':'all'},
                     outgroup_ids=tree.get('outgroup_tree_ids',[]), root_basis=tree.get('root_basis',''),
                     method={'support_kind':'UFBoot'})
        group['inputSha256'] = digest(json.dumps(group,sort_keys=True))
        panels['analysis']['groups'][tree['group']] = group
        tree['review_panel'] = 'analysis'
    panels['analysis']['seal'] = digest(json.dumps(panels['analysis']['groups'],sort_keys=True))
    if panels['analysis']['groups']:
        write_bundle(root, panels)


def run_cli(args):
    raw = json.loads(Path(args.models).read_text())
    models = raw['models'] if isinstance(raw, dict) else raw
    if not isinstance(models, list):
        raise ValueError('Models JSON must contain a list of model_id/tree_id/taxon records')
    result = assess(Path(args.tree_file).read_text(), models, support_kind=args.support_kind,
                    support_field=args.support_field, details=args.details, **reference_options(args))
    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=False)
    (root/'placement.json').write_text(json.dumps(result, indent=2)+'\n')
    columns = ['model_id','taxon','tree_id','proposed_points','included_in_score','maximum','status','reason','topology_agreement','weighted_agreement','mean_support','coverage','assessed','concordant','discordant']
    with (root/'placement.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t'); writer.writeheader()
        for row in result['models'].values():
            writer.writerow({key: row.get(key, row['counts'].get(key, 'NA')) if row.get(key) is not None else 'NA' for key in columns})
    manifest = dict(schema='gcf_tree_placement_run_v1', method=METHOD, model_input_sha256=digest(Path(args.models).read_bytes()),
                    tree_sha256=result['tree_sha256'], reference_sha256=result['reference_sha256'],
                    source_code_sha256=digest(Path(__file__).read_bytes()), biological_approval_asserted=False,
                    counts=dict(Counter(row['status'] for row in result['models'].values())))
    (root/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(output=str(root/'placement.json'), method=METHOD, counts=manifest['counts'])))
    return 0
