"""Generate a portable, Python-standard-library automatic tree-review bundle.

Public API: write_bundle(root, panels, initial_reviews=None),
validate_panel(manifest, item), prepare_overlay(reviews, manifest), serve(...).
Review items and overlays carry canonical `request` fields for revalidation.
"""
import csv
from io import StringIO
import json
from pathlib import Path
import os
import secrets

from .tree_review_service import validate_panel, prepare_overlay, serve, atomic, sha


IQTREE_DOWNLOADS = {
    'darwin-arm64': ('iqtree-3.0.1-macOS-arm.zip', 'iqtree-3.0.1-macOS-arm/bin/iqtree3', '0c9d39915ef773c06d261af738499a71b706068ff07707b1f0c3c260dfaea664'),
    'darwin-amd64': ('iqtree-3.0.1-macOS-intel.zip', 'iqtree-3.0.1-macOS-intel/bin/iqtree3', '366f65557f21662cfa91a96b2d208cfb53efcd6f17c005b341cd6c764c524d13'),
    'linux-amd64': ('iqtree-3.0.1-Linux-intel.tar.gz', 'iqtree-3.0.1-Linux-intel/bin/iqtree3', '2fe1f38daba07a7fb00b696572291e9a72a4e18918619c57c70f63857096bdbc'),
    'linux-arm64': ('iqtree-3.0.1-Linux-arm.tar.gz', 'iqtree-3.0.1-Linux-arm/bin/iqtree3', '8a2c72418356b6e7a43462322fc33d7b68a6a485662ef92b6fa0b512d723db19'),
    'windows-amd64': ('iqtree-3.0.1-Windows.zip', 'iqtree-3.0.1-Windows/bin/iqtree3.exe', '16cd257d58e4677a1f50657f443f6b2490916746a23c8451ebadeaf9ac1e5d43'),
}


def _write(path, data):
    data = data.encode() if isinstance(data, str) else data
    temp = path.with_name('.' + path.name + '.' + secrets.token_hex(8))
    try:
        with temp.open('xb') as f: f.write(data)
        os.replace(temp, path)
    finally: temp.unlink(missing_ok=True)


def write_bundle(root, panels, initial_reviews=None):
    """Write generated review assets only; scientific sources and .review stay intact.

    panels are caller-provided arbitrary groups/genes, never a study-specific list.
    Supply group.mask (zero-based columns or a mask descriptor), outgroup_ids,
    and root_basis when retaining Stage3 construction choices. They cannot be
    recovered reliably from a Newick display root alone.
    Optional v2 initial reviews are recomputed before being retained as a seed.
    """
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).parent
    manifest = {'schema': 'gcf_tree_review_inputs_v2', 'panels': json.loads(json.dumps(panels, allow_nan=False)),
        'core_method': 'species_quartet_placement_v1', 'core_sha256': sha((source / 'tree_placement.py').read_bytes())}
    for panel in manifest['panels'].values():
        if not isinstance(panel.get('groups'), dict): raise ValueError('Panel needs tree groups')
        for group in panel['groups'].values():
            if not isinstance(group.get('inputSha256'), str) or not group['inputSha256']: raise ValueError('Group needs stable inputSha256')
    # Validate definitions without approving or beginning any scientific calculation.
    from .tree_review_service import default_review
    prepare_overlay(default_review(manifest), manifest)
    atomic(root / 'tree_inputs.json', manifest)
    for original, target in [('tree_review_service.py', 'START_REVIEW.py'), ('tree_placement.py', 'tree_placement.py'),
                             ('tree_review.html', 'TREE_REVIEW.html'), ('tree_review_view.js', 'tree_review_view.js')]:
        _write(root / target, (source / original).read_bytes())
    profile = source / 'tree_placement_taxonomy.tsv'
    if profile.exists(): _write(root / profile.name, profile.read_bytes())
    buffer = StringIO(); writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
    writer.writerow(['panel', 'group', 'model_id', 'tree_id', 'taxon', 'protein_sha256', 'primary', 'species_tree_sha256'])
    for key, panel in manifest['panels'].items():
        for gid, group in panel['groups'].items():
            for model in group.get('models', []):
                writer.writerow([key, gid, model['model_id'], model.get('tree_id', ''), model.get('taxon', ''),
                    model.get('protein_sha256', ''), model.get('primary', True), sha(group.get('species_tree', ''))])
    _write(root / 'tree_profile.tsv', buffer.getvalue())
    downloads = {key: {'URL': 'https://github.com/iqtree/iqtree3/releases/download/v3.0.1/' + archive,
        'SHA256': digest, 'Archive': archive, 'Entry': entry} for key, (archive, entry, digest) in IQTREE_DOWNLOADS.items()}
    atomic(root / 'iqtree_manifest.json', downloads)
    if initial_reviews is not None:
        verified = prepare_overlay(initial_reviews, manifest)
        atomic(root / 'tree_initial_reviews.json', verified['request'])
    return manifest
