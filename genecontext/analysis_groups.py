"""Explicit STOP-1 analysis groups; membership never approves a locus."""
import hashlib
import re

DRAFT_SCHEMAS = {'stage1_pair_draft_v1', 'stage1_pair_draft_v2'}


def normalize_groups(data, groups, selected):
    if not isinstance(groups, list):
        raise ValueError('Analysis groups must be a list')
    loci = {l['locus_id']: l['genome_id'] for l in data['loci']}
    ids, names, assigned, result = set(), set(), set(), []
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError('Invalid analysis group')
        gid, name = group.get('group_id'), group.get('name')
        if not isinstance(gid, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', gid):
            raise ValueError('Invalid analysis group ID')
        if (not isinstance(name, str) or not name.strip() or len(name) > 120
                or name != name.strip() or any(ord(c) < 32 or ord(c) == 127 for c in name)):
            raise ValueError('Analysis group names must be nonempty, at most 120 characters, without control characters')
        if gid in ids or name.casefold() in names:
            raise ValueError('Duplicate analysis group ID or name')
        ids.add(gid); names.add(name.casefold())
        members = group.get('members')
        if not isinstance(members, list):
            raise ValueError('Analysis group members must be a list')
        normalized = []
        for member in members:
            if not isinstance(member, dict):
                raise ValueError('Invalid analysis group member')
            key = member.get('genome_id'), member.get('locus_id')
            if not all(isinstance(v, str) and v for v in key):
                raise ValueError('Invalid analysis group member IDs')
            if loci.get(key[1]) != key[0] or key not in selected:
                raise ValueError('Analysis group contains an unknown or unselected locus')
            if key in assigned:
                raise ValueError('A selected locus must belong to at most one analysis group')
            assigned.add(key)
            normalized.append(dict(genome_id=key[0], locus_id=key[1]))
        result.append(dict(group_id=gid, name=name, members=normalized))
    return result


def draft_groups(data, draft):
    selected = {(p.get('genome_id'), lid) for p in draft['pairs']
                if isinstance(p, dict) and isinstance(p.get('genome_id'), str) and isinstance(p.get('selected'), list)
                for lid in p['selected'] if isinstance(lid, str)}
    return normalize_groups(data, draft.get('analysis_groups', []), selected)


def row_groups(data, rows):
    groups, selected = {}, set()
    for row in rows:
        gid, name = row.get('analysis_group_id', '').strip(), row.get('analysis_group', '').strip()
        if not gid and not name:
            continue
        if row.get('decision', '').strip() != 'inspect':
            raise ValueError('Only explicitly inspected loci may be assigned to analysis groups')
        name = name or gid
        gid = gid or 'group_' + hashlib.sha256(name.encode()).hexdigest()[:16]
        group = groups.setdefault(gid, dict(group_id=gid, name=name, members=[]))
        if group['name'] != name:
            raise ValueError('Conflicting names for one analysis group ID')
        member = dict(genome_id=row.get('genome_id'), locus_id=row.get('locus_id'))
        group['members'].append(member)
        selected.add((member['genome_id'], member['locus_id']))
    return normalize_groups(data, list(groups.values()), selected)


def choose_group(groups, requested=None, active=None):
    def match(value):
        matches = [g for g in groups if value in (g['group_id'], g['name'])]
        if len(matches) != 1:
            raise ValueError('Unknown or ambiguous analysis group: ' + str(value))
        if not matches[0]['members']:
            raise ValueError('The chosen analysis group has no selected loci')
        return matches[0]
    # A scoped draft is deliberately bound to its exported group. Use the full
    # review draft when switching groups, instead of silently changing its scope.
    declared = match(active) if active else None
    chosen = match(requested) if requested else declared
    if declared and chosen['group_id'] != declared['group_id']:
        raise ValueError('The requested group differs from this scoped group draft')
    if groups and chosen is None:
        raise ValueError('This review contains analysis groups; choose --analysis-group ID_OR_NAME or use a separately exported group JSON')
    if requested and not groups:
        raise ValueError('The review has no analysis groups')
    return chosen


def project_rows(rows, group):
    """Limit execution scope without inventing biological exclusions."""
    members = {(m['genome_id'], m['locus_id']) for m in group['members']}
    genomes = {g for g, _ in members}
    result = []
    for original in rows:
        if original.get('genome_id') not in genomes:
            continue
        row = dict(original)
        key = row.get('genome_id'), row.get('locus_id')
        if key in members:
            if row.get('decision') != 'inspect':
                raise ValueError('Group membership cannot approve an unselected locus')
            row.update(analysis_group_id=group['group_id'], analysis_group=group['name'])
        else:
            row.update(original_decision=row.get('decision', ''), decision='exclude',
                       reason='outside_analysis_group', selection_status='outside_analysis_group',
                       notes='Outside this analysis group; original decision is retained in the source input.',
                       analysis_group_id='', analysis_group='')
        row['pair_decision'] = 'inspect'
        result.append(row)
    return result


def validate_row_bindings(data, rows):
    known = {(l['genome_id'], l['locus_id']) for l in data['loci']}
    known.update((g['genome_id'], '') for g in data['genomes']
                 if not any(l['genome_id'] == g['genome_id'] for l in data['loci']))
    seen = set()
    for row in rows:
        key = row.get('genome_id'), row.get('locus_id')
        if row.get('review_id') != data['review_id'] or key not in known or key in seen:
            raise ValueError('Wrong checkpoint, unknown or duplicate locus in grouped decisions')
        seen.add(key)


def annotate_selection(regions, group):
    if group is None:
        return regions
    return [dict(r, analysis_group_id=group['group_id'], analysis_group=group['name']) for r in regions]
