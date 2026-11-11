"""Shared scientific constraints for STOP-2 presentation and its consumer."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from .sequence import union_length, sha256


def verify_checkpoint(path):
    root = Path(path).resolve()
    seal=json.loads((root / 'seal.json').read_text())
    if {'evidence.json','cds.fna','proteins.faa','windows.fna'}-set(seal):
        raise ValueError('Incomplete stage 2 scientific evidence seal')
    for name, digest in seal.items():
        file = (root / name).resolve()
        if root not in file.parents or not file.is_file() or sha256(file) != digest:
            raise ValueError('Stage 2 evidence changed or is missing: ' + name)
    data = json.loads((root / 'evidence.json').read_text())
    if data.get('schema') != 'stage2_review_v1':
        raise ValueError('Expected a locus-first stage2 checkpoint')
    from .stage2_quality import with_annotation_quality
    from .stage2_reconstruct import load
    return root, load(root, with_annotation_quality(root, data))


def constraints(data):
    """Reconcile observations without merging distinct annotated genes/fusions."""
    models = data['models']; by_group = defaultdict(list)
    for m in models:
        by_group[m['gene_group']].append(m)
    parent = {g: g for g in by_group}
    def find(g):
        while parent[g] != g:
            parent[g] = parent[parent[g]]; g = parent[g]
        return g
    annotations = {g: {m['original_gene'] for m in ms if m['source'] == 'annotation'}
                   for g, ms in by_group.items()}
    buckets = defaultdict(list)
    for g, ms in by_group.items():
        m = ms[0]; buckets[m['genome_id'], m['contig'], m['strand']].append(g)
    for groups in buckets.values():
        for i, ga in enumerate(groups):
            for gb in groups[i+1:]:
                if annotations[ga] and annotations[gb]:
                    continue
                compatible = False
                for a in by_group[ga]:
                    for b in by_group[gb]:
                        if a['region_id'] == b['region_id']:
                            continue
                        if any('overlaps_multiple_annotated_genes' in m.get('issues', []) for m in [a,b]):
                            continue
                        shared = union_length([(max(x[0], y[0]), min(x[1], y[1]))
                            for x in a['cds_blocks'] for y in b['cds_blocks']
                            if max(x[0], y[0]) < min(x[1], y[1])])
                        sizes = [union_length([x[:2] for x in m['cds_blocks']]) for m in [a,b]]
                        same_bounds = (a['start'], a['end']) == (b['start'], b['end'])
                        junctions = lambda m: set(zip([x[1] for x in sorted(m['cds_blocks'])][:-1],
                                                      [x[0] for x in sorted(m['cds_blocks'])][1:]))
                        if sizes and min(sizes) and shared / max(sizes) >= .6 and shared / min(sizes) >= .8 and (same_bounds or len(junctions(a) & junctions(b)) >= 2):
                            compatible = True
                if compatible:
                    ra, rb = find(ga), find(gb)
                    # A chain of predictions must never bridge two annotated genes.
                    anno_a = set().union(*(annotations[g] for g in groups if find(g) == ra))
                    anno_b = set().union(*(annotations[g] for g in groups if find(g) == rb))
                    if anno_a and anno_b and anno_a != anno_b:
                        continue
                    parent[rb] = ra
    merged = defaultdict(list)
    for g in parent: merged[find(g)].append(g)
    physical = {g: 'physical_' + hashlib.sha256(json.dumps(sorted(values)).encode()).hexdigest()[:24]
                for values in merged.values() for g in values}
    expected = {v.get('gene_symbol','').upper() for members in data.get('query_members', {}).values()
                for v in members if v.get('gene_symbol') and not v['gene_symbol'].upper().startswith('LOC')}
    from .stage2_quality import model_quality
    quality = model_quality(data)
    result = {}
    for m in models:
        warnings = []
        if not m['overlaps_anchor']: warnings.append('Outside selected interval')
        if not m.get('target_support'): warnings.append('No retained target-query support')
        ident = m.get('identification', {})
        names = {h.get('gene_symbol','').upper() for h in ident.get('near_best', [])} - {''}
        if expected and names and not names & expected:
            warnings.append('Protein identification points to another gene/family')
        if 'overlaps_multiple_annotated_genes' in m.get('issues', []):
            warnings.append('Prediction overlaps several annotated genes')
        warnings.extend(quality[m['model_id']]['warnings'])
        result[m['model_id']] = dict(physical_id=physical[m['gene_group']], target_warnings=warnings,
                                    quality=quality[m['model_id']])
    return result


def copy_label(model, choice):
    """A local copy identity, never an inferred cross-species orthology label."""
    return choice.get('copy_label', '').strip() or (
        f"{model.get('original_gene', model['model_id'])} · "
        f"{model['contig']}:{model['start']}–{model['end']} ({model['strand']})")


def target_warnings(model, choice, metadata):
    warnings = metadata[model['model_id']]['target_warnings']
    if choice.get('selection_kind') == 'additional_copy':
        # A deliberately added copy has its own model interval as the seed.
        # Query support, family identity and fusion warnings still apply.
        warnings = [w for w in warnings if w != 'Outside selected interval']
    return warnings


def validate_choices(data, state):
    """Return structured issues. No edits or automatic approvals are made."""
    regions = {r['region_id']: r for r in data['regions']}
    models = {m['model_id']: m for m in data['models']}
    issues = []
    def issue(text, rid=None, mid=None):
        issues.append(dict(message=text, region_id=rid, model_id=mid))
    if not isinstance(state, dict):
        issue('Decisions must be a JSON object'); return issues
    if state.get('schema') != 'stage2_draft_v1' or state.get('review_id') != data['review_id']:
        issue('Decisions belong to another stage 2 checkpoint'); return issues
    for section in ['regions', 'models', 'pairs']:
        if not isinstance(state.get(section, {}), dict):
            issue('Invalid decision section: ' + section); return issues
        if any(not isinstance(value, dict) for value in state.get(section, {}).values()):
            issue('Invalid decision record: ' + section); return issues
    if set(state.get('regions', {})) - set(regions): issue('Unknown region in decisions')
    if set(state.get('models', {})) - set(models): issue('Unknown model in decisions')
    for mid,c in state.get('models',{}).items():
        if (not isinstance(c.get('notes'),str) or not isinstance(c.get('role'),str) or
                type(c.get('target_override',False)) is not bool or
                not isinstance(c.get('comparison_group',''),str) or
                not isinstance(c.get('copy_label',''),str) or
                c.get('selection_kind','') not in {'','additional_copy'}):
            issue('Invalid model comment, role, override or comparison label',models.get(mid,{}).get('region_id'),mid)
            return issues
    metadata = constraints(data); selected = []; retained = []
    for rid, r in regions.items():
        choice = state.get('regions', {}).get(rid, {})
        if choice.get('saved') is not True or choice.get('decision') not in {'inspect','exclude','unresolved'}:
            issue('Save a region decision', rid)
        if not isinstance(choice.get('notes'), str) or not choice['notes'].strip():
            issue('Supply a region reason', rid)
        own = [m for m in models.values() if m['region_id'] == rid]
        kept = [m for m in own if state.get('models', {}).get(m['model_id'], {}).get('decision') == 'inspect']
        targets = [m for m in kept if state['models'][m['model_id']].get('role') == 'target']
        if choice.get('decision') == 'inspect' and not targets: issue('Choose a target CDS model', rid)
        if choice.get('decision') == 'exclude' and kept: issue('Clear selected models in an excluded region', rid)
        if choice.get('decision') != 'inspect' and kept: issue('Selected models conflict with excluded/unresolved region decision', rid)
        for m in kept:
            mid = m['model_id']; c = state['models'][mid]
            retained.append(m)
            if c.get('role') not in {'target','neighbor'} or not isinstance(c.get('notes'), str):
                issue('Invalid model role/comment', rid, mid); continue
            if type(c.get('target_override', False)) is not bool:
                issue('Expert override must be an explicit boolean', rid, mid)
            if m['coding_status'] == 'disrupted': issue('Disrupted CDS cannot enter comparison', rid, mid)
            if (m['coding_status'] != 'complete_predicted' or m['window_edge']) and not c['notes'].strip():
                issue('Explain partial/edge model selection', rid, mid)
            if c['role'] == 'target':
                warnings = target_warnings(m, c, metadata)
                if warnings and not (c.get('target_override') is True and c['notes'].strip()):
                    issue('; '.join(warnings) + ': explicit expert override and reason required', rid, mid)
                if any(x in c.get('copy_label','') for x in '\t\r\n'):
                    issue('Copy label must be a single line', rid, mid)
                label = c.get('comparison_group', r.get('gene','target'))
                if not isinstance(label, str) or not label.strip() or any(x in label for x in '\t\r\n'):
                    issue('Supply a valid comparison-group label', rid, mid)
                else:
                    selected.append((m, label, metadata[mid]['physical_id']))
    for mid, c in state.get('models', {}).items():
        if mid not in models: continue
        if not isinstance(c, dict) or c.get('decision') not in {'','inspect','reject','unresolved','not_selected'}:
            issue('Invalid model decision', models[mid]['region_id'], mid)
    by_retained = defaultdict(list)
    for m in retained:by_retained[metadata[m['model_id']]['physical_id']].append(m)
    for ms in by_retained.values():
        if len({m.get('cds_sha256') for m in ms})>1:
            for m in ms:issue('Conflicting coding models selected for one physical gene',m['region_id'],m['model_id'])
        per_region=defaultdict(list)
        for m in ms:per_region[m['region_id']].append(m)
        for own_models in per_region.values():
            if len(own_models)>1:
                for m in own_models:issue('Choose one model per physical gene in this region',m['region_id'],m['model_id'])
    by_physical = defaultdict(list); by_copy_label = defaultdict(list)
    for m, label, physical in selected:
        by_physical[physical].append(m)
        by_copy_label[m['genome_id'], label, copy_label(m, state['models'][m['model_id']])].append(m)
    for ms in by_physical.values():
        if len(ms) > 1:
            for m in ms: issue('One physical gene was selected more than once; keep one representative model/region', m['region_id'], m['model_id'])
    for ms in by_copy_label.values():
        if len({metadata[m['model_id']]['physical_id'] for m in ms}) > 1:
            for m in ms: issue('Use distinct copy labels for different physical copies in this genome/comparison group', m['region_id'], m['model_id'])
    valid = defaultdict(list)
    for p in data['neighbor_pairs']:
        if p.get('kind') == 'neighbor_proposal':
            valid[':'.join(sorted([p['gene_a'],p['gene_b']]))].append(p)
    for key, p in state.get('pairs', {}).items():
        originals = valid.get(key, [])
        if p.get('decision') != 'inspect' or not any(all(p.get(k) == original.get(k) for k in ['gene_a','gene_b','model_a','model_b']) for original in originals):
            issue('Unknown or changed neighboring-gene proposal'); continue
        for mid in [p['model_a'],p['model_b']]:
            rid = models[mid]['region_id']
            if state.get('regions', {}).get(rid, {}).get('decision') != 'inspect':
                issue('A selected neighbor pair belongs to a region that is excluded or unresolved', rid, mid)
    return issues


def issue_text(data, issue):
    regions = {r['region_id']:r for r in data['regions']}
    r = regions.get(issue.get('region_id'))
    if r:
        species = r['genome_id'].split('_GC')[0].replace('_',' ')
        return f"{species} · {r['contig']}:{r['anchor_start']}–{r['anchor_end']}: {issue['message']}"
    return issue['message']
