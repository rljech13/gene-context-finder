"""Self-contained local tree review; scores are always recomputed by Python core.

This file also ships unchanged as START_REVIEW.py beside tree_placement.py.
"""
import argparse
from copy import deepcopy
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
import os
from pathlib import Path
import platform
import re
import secrets
import shutil
import subprocess
import tarfile
import tempfile
import threading
from urllib.parse import unquote, urlsplit
from urllib.request import urlopen
import webbrowser
import zipfile

try:
    from . import tree_placement as core
except ImportError:
    import tree_placement as core

REVIEW_SCHEMA = 'gcf_tree_review_v2'
OVERLAY_SCHEMA = 'gcf_tree_review_overlay_v2'
STATE_NAME = 'tree_review_v2.json'
MAX_BODY = 12 * 1024 * 1024


def sha(value):
    return hashlib.sha256(value if isinstance(value, bytes) else value.encode()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)


def atomic(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name('.' + path.name + '.' + secrets.token_hex(8))
    try:
        with tmp.open('x', encoding='utf-8') as f:
            f.write(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
            f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def _parse_tree(text):
    """Small display parser; the scoring core independently validates Newick."""
    core.Tree(text)
    tokens = re.findall(r"'(?:[^']|'')*'|[(),:;]|[^\s(),:;]+", re.sub(r'\[[^\]]*\]', '', text))
    pos = 0
    def node():
        nonlocal pos
        children = []
        if tokens[pos] == '(':
            pos += 1; children.append(node())
            while tokens[pos] == ',': pos += 1; children.append(node())
            pos += 1
        label = ''
        if tokens[pos] not in {':', ',', ')', ';'}:
            label = tokens[pos]; pos += 1
            if label.startswith("'"): label = label[1:-1].replace("''", "'")
        length = None
        if tokens[pos] == ':':
            pos += 1; length = float(tokens[pos]); pos += 1
        return {'name': label, 'length': length, 'children': children}
    return node()


def orient_tree(text, outgroup):
    """Root on a real separating edge only; never invent a polyphyletic clade."""
    if not outgroup: return text
    tree = _parse_tree(text); adjacent = {}; names = {}
    def visit(n, parent=None):
        i = len(adjacent); adjacent[i] = {}; names[i] = n['name']
        if parent is not None:
            edge = (n['length'], n['name'] if n['children'] else '')
            adjacent[i][parent] = edge; adjacent[parent][i] = edge
        for c in n['children']: visit(c, i)
        return i
    root = visit(tree)
    if len(adjacent[root]) == 2:
        a, b = adjacent[root]; la, sa = adjacent[a].pop(root); lb, sb = adjacent[b].pop(root)
        edge = (None if la is None or lb is None else la + lb, sa or sb)
        adjacent[a][b] = edge; adjacent[b][a] = edge
        del adjacent[root]
    all_tips = {names[n] for n in adjacent if len(adjacent[n]) == 1}
    wanted = set(outgroup)
    if not wanted < all_tips: raise ValueError('Outgroup must be a proper nonempty subset of tree tips')
    def tips(n, parent):
        children = [c for c in adjacent[n] if c != parent]
        return {names[n]} if not children else set().union(*(tips(c, n) for c in children))
    match = None
    for a in adjacent:
        for b in adjacent[a]:
            if a < b and (tips(a, b) == wanted or tips(b, a) == wanted): match = (a, b)
    if match is None: raise ValueError('Outgroup is not monophyletic on an existing unrooted edge; keep unrooted or change the selection')
    def quoted(label):
        return label if re.fullmatch(r'[A-Za-z0-9_.+-]+', label) else "'" + label.replace("'", "''") + "'"
    def emit(n, parent, incoming_length, incoming_support):
        children = [c for c in adjacent[n] if c != parent]
        body = '(' + ','.join(emit(c, n, *adjacent[n][c]) for c in children) + ')' if children else quoted(names[n])
        if children and incoming_support: body += incoming_support
        if incoming_length is not None: body += ':' + format(incoming_length, '.15g')
        return body
    a, b = match; length, support = adjacent[a][b]
    half = None if length is None else length / 2
    return '(' + emit(a, b, half, support) + ',' + emit(b, a, half, support) + ');'


def _groups(manifest, panel):
    panels = manifest.get('panels', manifest)
    if panel not in panels: raise ValueError('Unknown panel')
    return panels[panel]


def alignment_input(group, selected, mask=None):
    rows = group.get('rows', {})
    if not isinstance(selected, list) or not selected or not all(isinstance(i, str) for i in selected) or len(set(selected)) != len(selected):
        raise ValueError('Select distinct aligned sequence IDs')
    if set(selected) - set(rows): raise ValueError('Unknown aligned sequence ID')
    selected = sorted(selected)
    width = len(rows[selected[0]])
    if not width or any(len(rows[i]) != width for i in selected): raise ValueError('Alignment rows have inconsistent widths')
    mask = mask if mask is not None else group.get('mask', {'mode': 'all'})
    if isinstance(mask, list):
        if not mask or any(type(i) is not int or not 0 <= i < width for i in mask) or mask != sorted(set(mask)):
            raise ValueError('Mask must be sorted distinct zero-based column indices')
        columns = mask
    elif isinstance(mask, dict) and mask.get('mode') in {'all', 'occupancy_50'}:
        columns = list(range(width)) if mask['mode'] == 'all' else [j for j in range(width)
            if sum(rows[i][j].upper() in 'ACDEFGHIKLMNPQRSTVWY' for i in selected) * 2 >= len(selected)]
    else: raise ValueError('Use all columns, occupancy_50, or exact column indices')
    if not columns: raise ValueError('The alignment mask retains no columns')
    masked = {i: ''.join(rows[i][j] for j in columns) for i in selected}
    if any(not any(c.upper() in 'ACDEFGHIKLMNPQRSTVWY' for c in s) for s in masked.values()):
        raise ValueError('A selected sequence has no canonical residues after masking')
    fasta = ''.join('>' + i + '\n' + masked[i] + '\n' for i in selected)
    return selected, mask, columns, fasta


def _tree_source(group, request):
    source = request.get('source')
    if isinstance(source, str): source = {'id': source}
    if not isinstance(source, dict): raise ValueError('Choose a saved or rebuilt tree source')
    sid = source.get('id')
    trees = group.get('trees', {})
    if sid in trees:
        newick = trees[sid]
        if isinstance(newick, dict): newick = newick.get('newick')
        provenance = {'kind': 'saved', 'sourceId': sid}
    elif isinstance(sid, str) and (sid == 'rebuilt' or sid.startswith('rebuilt:')):
        newick = source.get('newick'); provenance = source.get('provenance')
        if not isinstance(provenance, dict) or provenance.get('kind') != 'rebuilt':
            raise ValueError('Rebuilt trees need their alignment and execution provenance')
    else: raise ValueError('Unknown tree source')
    if not isinstance(newick, str): raise ValueError('Tree source has no Newick')
    if source.get('sha256') is not None and source['sha256'] != sha(newick): raise ValueError('Tree checksum does not match')
    core.Tree(newick)
    return sid, newick, provenance


def _validate_group(group, request, binding=None):
    if not isinstance(request, dict): raise ValueError('Invalid tree group review')
    if request.get('inputSha256') != group.get('inputSha256'): raise ValueError('Review belongs to another alignment/model input')
    sid, newick, provenance = _tree_source(group, request)
    tree = core.Tree(newick)
    rows = group.get('rows', {})
    if set(tree.tips) - set(rows): raise ValueError('Tree contains a tip outside the exact aligned sequence pool')
    # A saved tree can omit models; omitted models receive explicit NA from core.
    selection = request.get('selected_ids', provenance.get('selected_ids', sorted(tree.tips)))
    selected, mask, columns, fasta = alignment_input(group, selection, request.get('mask', provenance.get('mask')))
    stale = set(selected) != set(tree.tips)
    if provenance.get('kind') == 'rebuilt':
        if provenance.get('alignment_sha256') != sha(fasta): stale = True
        if provenance.get('selected_ids') != selected or provenance.get('columns') != columns: stale = True
    else:
        baseline_mask = group.get('mask', {'mode': 'all'})
        if canonical(mask) != canonical(baseline_mask): stale = True
    default_tree_method = {'model': 'MFP', 'seed': 42, 'bootstrap': 1000, 'threads': 2}
    source_tree_method = provenance.get('tree_method', group.get('tree_method', default_tree_method))
    tree_method = request.get('tree_method', source_tree_method)
    if not isinstance(tree_method, dict): raise ValueError('Invalid tree construction method')
    if canonical(tree_method) != canonical(source_tree_method): stale = True
    outgroup = request.get('outgroup_ids', [])
    if not isinstance(outgroup, list) or any(not isinstance(i, str) for i in outgroup) or len(outgroup) != len(set(outgroup)):
        raise ValueError('Outgroup must contain distinct sequence IDs')
    aliases = {leaf['tree_id']: set(leaf.get('sequence_ids', [])) | {leaf['tree_id']} for leaf in group.get('leaves', [])}
    mapped = set()
    for name in outgroup:
        matches = [tip for tip in tree.tips if name == tip or name in aliases.get(tip, set())]
        if len(matches) != 1: raise ValueError('Outgroup is not a unique included tree sequence')
        mapped.add(matches[0])
    basis = request.get('root_basis', '')
    comment = request.get('comment', '')
    if not isinstance(basis, str) or len(basis) > 4000 or not isinstance(comment, str) or len(comment) > 32000:
        raise ValueError('Invalid root basis or review comment')
    if mapped and not basis.strip(): raise ValueError('Explain the outgroup basis')
    # Collapsed exact aliases cannot be divided between ingroup and outgroup.
    for tip in mapped:
        original_aliases = aliases.get(tip, set()) - {tip}
        if tip not in outgroup and original_aliases and not original_aliases.issubset(outgroup):
            raise ValueError('Select the collapsed leaf or all its identical aliases as outgroup')
    display = orient_tree(newick, sorted(mapped))
    method = request.get('method', group.get('method', {})) or {}
    if not isinstance(method, dict) or set(method) - {'support_field', 'support_kind'}: raise ValueError('Unknown placement-method options')
    field = method.get('support_field')
    if field is not None and (type(field) is not int or not 0 <= field <= 4): raise ValueError('Invalid support field')
    kind = method.get('support_kind', 'reported_branch_support')
    if not isinstance(kind, str) or not 1 <= len(kind) <= 100: raise ValueError('Invalid support kind')
    reference = group.get('species_tree')
    if not isinstance(reference, str) or not reference.strip(): raise ValueError('An independent species reference tree is required')
    models = deepcopy(group.get('models', []))
    for model in models:
        if model.get('tree_id') not in selected: model['tree_id'] = None
    placement = core.assess(newick, models, species_tree=reference, support_field=field, support_kind=kind)
    signature_data = {'binding': binding or {}, 'tree_sha256': sha(newick), 'inputSha256': group.get('inputSha256'),
        'selected_ids': selected, 'alignment_sha256': sha(fasta), 'alignment_pool_sha256': sha(canonical(rows)),
        'mask': mask, 'columns': columns, 'outgroup_ids': sorted(outgroup), 'resolved_outgroup_ids': sorted(mapped),
        'root_basis': basis.strip(), 'reference_sha256': sha(reference), 'method': method,
        'core_method': core.METHOD, 'core_sha256': sha(Path(core.__file__).read_bytes()), 'tree_method': tree_method,
        'models': models, 'source_provenance': provenance}
    signature = sha(canonical(signature_data))
    approved = request.get('approved') is True and request.get('confirmed_signature') == signature and not stale
    canonical_request = {'inputSha256': group.get('inputSha256'),
        'source': {'id': sid, 'sha256': sha(newick)}, 'selected_ids': selected, 'mask': mask,
        'outgroup_ids': sorted(outgroup), 'root_basis': basis.strip(), 'comment': comment,
        'method': method, 'tree_method': tree_method, 'approved': approved, 'confirmed_signature': signature if approved else None}
    if provenance.get('kind') == 'rebuilt': canonical_request['source'].update(newick=newick, provenance=deepcopy(provenance))
    return {'newick': newick, 'display_newick': display, 'sourceId': sid, 'outgroup_ids': sorted(outgroup),
        'outgroup_tree_ids': sorted(mapped), 'alignment_fasta': fasta,
        'root_basis': basis.strip(), 'comment': comment, 'approved': approved, 'signature': signature,
        'confirmed_signature': signature if approved else None, 'placement': placement,
        'selected_ids': selected, 'mask': mask, 'columns': columns, 'method': method, 'tree_method': tree_method,
        'requires_rebuild': stale, 'approval_reason': 'tree_input_changed_rebuild_required' if stale else
            'reviewed_exact_tree_and_outgroup' if approved else 'pending_exact_tree_and_outgroup_review',
        'provenance': dict(provenance, signature_inputs=signature_data), 'request': canonical_request}


def validate_panel(manifest, item):
    """Return a recomputed overlay panel; caller-supplied numeric decisions are ignored."""
    if not isinstance(item, dict) or not isinstance(item.get('panel'), str): raise ValueError('Review needs a panel ID')
    panel = item['panel']; source = _groups(manifest, panel)
    if item.get('seal') != source.get('seal'): raise ValueError('Panel seal differs from its source')
    requests = item.get('groups', {})
    if not isinstance(requests, dict) or set(requests) - set(source.get('groups', {})): raise ValueError('Unknown review group')
    result = {'seal': source.get('seal'), 'groups': {}, 'models': {}}
    for gid, request in requests.items():
        group = _validate_group(source['groups'][gid], request, {'panel': panel, 'seal': source.get('seal'), 'group': gid}); result['groups'][gid] = group
        for mid, placement in group['placement']['models'].items():
            if mid in result['models']: raise ValueError('Model occurs in more than one tree group')
            raw = placement['points']
            result['models'][mid] = {'points': raw if group['approved'] else None, 'proposed_points': raw,
                'placement': placement, 'approved': group['approved'], 'group_id': gid,
                'signature': group['signature'], 'status': placement['status'] if group['approved'] else 'pending_review',
                'maximum': 15, 'reason': placement['reason'] or ('' if group['approved'] else group['approval_reason'])}
    # Never inherit old manual points for a model in an unreviewed group.
    for gid, group in source.get('groups', {}).items():
        for model in group.get('models', []):
            mid = model['model_id']
            if mid not in result['models']:
                result['models'][mid] = {'points': None, 'proposed_points': None, 'placement': None,
                    'approved': False, 'group_id': gid, 'status': 'pending_review', 'maximum': 15, 'reason': 'tree_group_not_reviewed'}
    result['request'] = {'panel': panel, 'seal': source.get('seal'),
        'groups': {gid: group['request'] for gid, group in result['groups'].items()}}
    return result


def prepare_overlay(reviews, manifest):
    if not isinstance(reviews, dict) or reviews.get('schema') != REVIEW_SCHEMA:
        raise ValueError('Expected gcf_tree_review_v2; legacy manual tree scores are not approved automatically')
    panels = reviews.get('panels', {})
    if not isinstance(panels, dict): raise ValueError('Review panels must be an object')
    out = {'schema': OVERLAY_SCHEMA, 'method': core.METHOD, 'panels': {}}
    for key, item in panels.items():
        if not isinstance(item, dict) or item.get('panel', key) != key: raise ValueError('Panel ID mismatch')
        out['panels'][key] = validate_panel(manifest, dict(item, panel=key))
    out['request'] = {'schema': REVIEW_SCHEMA, 'panels': {key: panel['request'] for key, panel in out['panels'].items()}}
    return out


def default_review(manifest):
    panels = {}
    for key, panel in manifest.get('panels', manifest).items():
        groups = {}
        for gid, group in panel.get('groups', {}).items():
            if not group.get('trees'): continue
            sid = next(iter(group['trees'])); text = group['trees'][sid]
            if isinstance(text, dict): text = text['newick']
            groups[gid] = {'inputSha256': group.get('inputSha256'), 'source': {'id': sid, 'sha256': sha(text)},
                'selected_ids': sorted(core.Tree(text).tips), 'mask': group.get('mask', {'mode': 'all'}),
                'outgroup_ids': list(group.get('outgroup_ids', group.get('outgroup_tree_ids', group.get('outgroup_sequence_ids', [])))),
                'root_basis': group.get('root_basis', ''), 'comment': group.get('comment', ''), 'approved': False, 'confirmed_signature': None}
        panels[key] = {'panel': key, 'seal': panel.get('seal'), 'groups': groups}
    return {'schema': REVIEW_SCHEMA, 'panels': panels}


def _iqtree(root, download):
    for name in ['iqtree3', 'iqtree2', 'iqtree']:
        found = shutil.which(name)
        if found: return Path(found)
    system = {'Darwin': 'darwin', 'Linux': 'linux', 'Windows': 'windows'}.get(platform.system())
    machine = platform.machine().lower(); arch = 'arm64' if machine in {'arm64', 'aarch64'} else 'amd64' if machine in {'amd64', 'x86_64'} else machine
    if system == 'windows' and arch == 'arm64': arch = 'amd64'
    manifest = json.loads((Path(root) / 'iqtree_manifest.json').read_text()); config = manifest.get(str(system) + '-' + arch)
    if not config: raise ValueError('Install a suitable native IQ-TREE for this platform')
    toolroot = Path(root) / '.review' / 'tools' / config['SHA256']; binary = toolroot / ('iqtree.exe' if system == 'windows' else 'iqtree')
    receipt = toolroot / 'receipt.json'
    if binary.is_file() and receipt.is_file() and sha(binary.read_bytes()) == json.loads(receipt.read_text()).get('binary_sha256'): return binary
    if not download: raise ValueError('IQ-TREE is missing; explicitly enable its pinned download or install it on PATH')
    toolroot.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryFile() as archive:
        digest = hashlib.sha256(); total = 0
        with urlopen(config['URL'], timeout=60) as response:
            while True:
                block = response.read(1024 * 1024)
                if not block: break
                total += len(block)
                if total > 512 * 1024 * 1024: raise ValueError('IQ-TREE download exceeds limit')
                digest.update(block); archive.write(block)
        if digest.hexdigest() != config['SHA256']: raise ValueError('IQ-TREE archive checksum failed')
        archive.seek(0)
        if config['Archive'].endswith('.zip'):
            with zipfile.ZipFile(archive) as z:
                info = z.getinfo(config['Entry'])
                if info.file_size > 256 * 1024 * 1024: raise ValueError('IQ-TREE binary exceeds limit')
                data = z.read(info)
        else:
            with tarfile.open(fileobj=archive, mode='r:*') as z:
                info = z.getmember(config['Entry'])
                if not info.isfile() or info.size > 256 * 1024 * 1024: raise ValueError('Invalid IQ-TREE binary entry')
                data = z.extractfile(info).read()
        temp = toolroot / ('pending-' + secrets.token_hex(6)); temp.write_bytes(data); temp.chmod(0o700); os.replace(temp, binary)
        atomic(receipt, {'archive_sha256': digest.hexdigest(), 'binary_sha256': sha(data), 'url': config['URL']})
    return binary


def rebuild(root, manifest, item, *, runner=None):
    panel = _groups(manifest, item.get('panel'))
    if item.get('seal') != panel.get('seal'): raise ValueError('Panel seal differs')
    gid = item.get('group'); group = panel.get('groups', {}).get(gid)
    if group is None or item.get('inputSha256') != group.get('inputSha256'): raise ValueError('Unknown group/input')
    selected, mask, columns, fasta = alignment_input(group, item.get('selected_ids'), item.get('mask'))
    if len(selected) < 4: raise ValueError('IQ-TREE review requires at least four selected leaves')
    settings = item.get('tree_method', {})
    model = settings.get('model', 'MFP')
    if model not in {'MFP', 'LG', 'WAG', 'JTT'}: raise ValueError('Unsupported substitution-model setting')
    seed = settings.get('seed', 42); bootstrap = settings.get('bootstrap', 1000); threads = settings.get('threads', 2)
    if any(type(v) is not int for v in [seed, bootstrap, threads]) or not 0 <= seed <= 2147483647 or not 1000 <= bootstrap <= 5000 or not 1 <= threads <= 8:
        raise ValueError('Invalid tree seed, bootstrap count or thread count')
    settings = dict(model=model, seed=seed, bootstrap=bootstrap, threads=threads)
    key = sha(canonical({'panel': item['panel'], 'group': gid, 'input': group['inputSha256'], 'fasta': sha(fasta), 'settings': settings}))
    job = Path(root) / '.review' / 'trees' / key; job.mkdir(parents=True, exist_ok=True)
    (job / 'alignment.faa').write_text(fasta)
    engine = _iqtree(root, item.get('allow_download') is True)
    argv = [str(engine), '-s', 'alignment.faa', '-st', 'AA', '-keep-ident', '-m', model, '-B', str(bootstrap),
        '-bnni', '-seed', str(seed), '-T', str(threads), '--prefix', 'ml', '-redo']
    result = (runner or subprocess.run)(argv, cwd=job, capture_output=True, text=True, timeout=1800)
    (job / 'stdout.log').write_text(result.stdout); (job / 'stderr.log').write_text(result.stderr)
    if result.returncode: raise ValueError('IQ-TREE failed; logs remain in the package .review/trees directory')
    newick = (job / 'ml.treefile').read_text(); observed = core.Tree(newick)
    if set(observed.tips) != set(selected): raise ValueError('IQ-TREE output differs from selected leaf set')
    version = subprocess.run([str(engine), '--version'], capture_output=True, text=True, timeout=15).stdout.splitlines()[:3]
    provenance = {'kind': 'rebuilt', 'alignment_sha256': sha(fasta), 'selected_ids': selected, 'mask': mask,
        'columns': columns, 'tree_method': settings, 'engine_version': version, 'command': ['iqtree'] + argv[1:],
        'tree_sha256': sha(newick), 'stdout_sha256': sha(result.stdout), 'stderr_sha256': sha(result.stderr)}
    source = {'id': 'rebuilt:' + key, 'sha256': sha(newick), 'newick': newick, 'provenance': provenance}
    atomic(job / 'result.json', source)
    return {'source': source, 'selected_ids': selected, 'mask': mask, 'approved': False, 'confirmed_signature': None}


class ReviewState:
    def __init__(self, root):
        self.root = Path(root).resolve(); self.lock = threading.RLock(); self.token = secrets.token_hex(32)
        self.manifest = json.loads((self.root / 'tree_inputs.json').read_text())
        self.file = self.root / '.review' / STATE_NAME
        initial = self.root / 'tree_initial_reviews.json'
        if not self.file.exists() and initial.exists():
            # An exported session is a seed only after exact-input/core revalidation.
            verified = prepare_overlay(json.loads(initial.read_text()), self.manifest)
            atomic(self.file, verified['request'])

    def reviews(self):
        result = default_review(self.manifest)
        if not self.file.exists(): return result
        stored = json.loads(self.file.read_text())
        if stored.get('schema') != REVIEW_SCHEMA: raise ValueError('Saved reviews have an incompatible schema')
        # Preserve the old file/history, while new inputs start pending instead
        # of inheriting approval or hiding newly added groups in the UI.
        for key, panel in stored.get('panels', {}).items():
            current = result['panels'].get(key)
            if current is None or panel.get('seal') != current['seal']: continue
            for gid, group in panel.get('groups', {}).items():
                expected = _groups(self.manifest, key)['groups'].get(gid)
                if expected and group.get('inputSha256') == expected.get('inputSha256'):
                    current['groups'][gid] = group
        return result

    def save(self, data):
        with self.lock:
            if data.get('schema') == REVIEW_SCHEMA:
                incoming = data
            else:
                key = data.get('panel'); incoming = {'schema': REVIEW_SCHEMA, 'panels': {key: data}}
            # Validation always invokes core; client placement/points/decisions never enter state.
            calculated = prepare_overlay(incoming, self.manifest)
            old = self.reviews()
            for key, item in incoming['panels'].items():
                clean = {'panel': key, 'seal': item['seal'], 'groups': {}}
                for gid, group in calculated['panels'][key]['groups'].items():
                    original = item['groups'][gid]
                    clean['groups'][gid] = {k: deepcopy(original[k]) for k in ['inputSha256', 'source']}
                    clean['groups'][gid].update({k: group[k] for k in ['outgroup_ids', 'root_basis', 'comment', 'approved',
                        'confirmed_signature', 'selected_ids', 'mask', 'method', 'tree_method']})
                if key not in old['panels'] or old['panels'][key].get('seal') != clean['seal']:
                    old['panels'][key] = clean
                else:
                    old['panels'][key]['groups'].update(clean['groups'])
            # Preserve history without overwriting v1 or any scientific source.
            atomic(self.root / '.review' / 'history' / (sha(canonical(old)) + '.json'), old)
            atomic(self.file, old)
            return {'reviews': old, 'overlay': prepare_overlay(old, self.manifest)}


def create_server(root, port=0):
    state = ReviewState(root)
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args): pass

        def send_json(self, code, value):
            data = json.dumps(value, ensure_ascii=False, allow_nan=False).encode()
            self.send_response(code); self.send_header('Content-Type', 'application/json; charset=utf-8')
            self.send_header('Content-Length', str(len(data))); self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff'); self.end_headers(); self.wfile.write(data)

        def guard(self, write=False):
            host = self.headers.get('Host', '')
            try:
                parsed = urlsplit('http://' + host)
                valid = parsed.hostname in {'127.0.0.1', 'localhost'} and parsed.port == self.server.server_port and not parsed.username and parsed.path == ''
            except ValueError: valid = False
            if not valid: raise ValueError('Loopback Host required')
            for name in ['Origin', 'Referer']:
                raw = self.headers.get(name)
                if raw:
                    origin = urlsplit(raw)
                    if origin.scheme != 'http' or origin.netloc != host: raise ValueError('Same-origin request required')
            if self.headers.get('Sec-Fetch-Site') in {'cross-site', 'same-site'}: raise ValueError('Same-origin request required')
            if write and not secrets.compare_digest(self.headers.get('X-GCF-Token', ''), state.token): raise ValueError('Session token required')

        def do_GET(self):
            try:
                self.guard(); path = unquote(urlsplit(self.path).path)
                if path == '/api/manifest': return self.send_json(200, {'manifest': state.manifest, 'token': state.token})
                if path == '/api/state':
                    with state.lock:
                        saved = state.reviews()
                        try: overlay = prepare_overlay(saved, state.manifest); warnings = []
                        except ValueError as e: overlay = None; warnings = [str(e)]
                    return self.send_json(200, {'reviews': saved, 'overlay': overlay, 'warnings': warnings})
                if path == '/': path = '/TREE_REVIEW.html'
                name = path[1:]; relative = Path(name)
                protected = {'tree_inputs.json', 'tree_initial_reviews.json', 'iqtree_manifest.json'}
                if not name or any(part.startswith('.') for part in relative.parts) or relative.name in protected:
                    return self.send_json(404, {'error': 'Not found'})
                import mimetypes
                public_suffixes = {'.html', '.htm', '.js', '.mjs', '.css', '.svg', '.png', '.jpg', '.jpeg', '.webp', '.gif', '.ico',
                    '.tsv', '.csv', '.json', '.txt', '.fna', '.faa', '.fa', '.fasta', '.gff', '.gff3', '.pdb', '.cif', '.mmcif',
                    '.nwk', '.newick', '.treefile', '.nex', '.nexus', '.pdf', '.zip'}
                if relative.suffix.lower() not in public_suffixes: return self.send_json(404, {'error': 'Not found'})
                file = (state.root / relative).resolve()
                if not file.is_relative_to(state.root) or not file.is_file(): return self.send_json(404, {'error': 'Not found'})
                content_type = mimetypes.guess_type(name)[0] or 'application/octet-stream'
                data = file.read_bytes(); self.send_response(200); self.send_header('Content-Type', content_type); self.send_header('Content-Length', str(len(data)))
                self.send_header('Cache-Control', 'no-store'); self.send_header('X-Content-Type-Options', 'nosniff')
                script_policy = "'self'" if relative.name == 'TREE_REVIEW.html' else "'self' 'unsafe-inline'"
                self.send_header('Content-Security-Policy', "default-src 'self' data: blob:; script-src " + script_policy + "; style-src 'self' 'unsafe-inline'; object-src 'none'; frame-ancestors 'self'; base-uri 'self'")
                self.end_headers(); self.wfile.write(data)
            except (ValueError, FileNotFoundError) as e: self.send_json(403, {'error': str(e)})

        def do_POST(self):
            try:
                self.guard(True)
                if self.headers.get_content_type() != 'application/json': raise ValueError('JSON required')
                length = int(self.headers.get('Content-Length', '-1'))
                if not 0 <= length <= MAX_BODY: raise ValueError('Request exceeds allowed size')
                data = json.loads(self.rfile.read(length)); path = urlsplit(self.path).path
                if not isinstance(data, dict): raise ValueError('Request must be an object')
                if path in {'/api/assess', '/api/preview'}:
                    return self.send_json(200, validate_panel(state.manifest, data))
                if path in {'/api/save', '/api/export'}: return self.send_json(200, state.save(data))
                if path == '/api/rebuild':
                    with state.lock: result = rebuild(state.root, state.manifest, data)
                    return self.send_json(200, result)
                self.send_json(404, {'error': 'Unknown operation'})
            except (ValueError, TypeError, KeyError, OSError, subprocess.SubprocessError) as e:
                self.send_json(400, {'error': str(e)})
    server = ThreadingHTTPServer(('127.0.0.1', port), Handler); server.review_state = state
    server.daemon_threads = True
    return server


def serve(root, port=0, no_browser=False):
    server = create_server(root, port)
    url = 'http://127.0.0.1:' + str(server.server_port) + '/TREE_REVIEW.html'
    print('Tree review: ' + url, flush=True)
    if not no_browser: webbrowser.open(url)
    try: server.serve_forever()
    except KeyboardInterrupt: pass
    finally: server.server_close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--root', default=str(Path(__file__).parent))
    parser.add_argument('--port', type=int, default=0); parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args(); serve(args.root, args.port, args.no_browser)
