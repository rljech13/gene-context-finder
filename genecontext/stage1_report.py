"""Offline, locus-first review of existing BLAST evidence and early synteny."""
import json
from pathlib import Path

TEMPLATE = Path(__file__).with_name("stage1_review.txt").read_text().replace(
    "__CONTEXT_VIEW__", Path(__file__).with_name("stage1_context_view.txt").read_text())

def write_report(root, data):
    payload = json.dumps(data).replace('<', chr(92)+'u003c').replace('&', chr(92)+'u0026')
    (root/'report.html').write_text(TEMPLATE.replace('__DATA__', payload))


def refresh(args):
    """Make a presentation revision while retaining the exact evidence/decision IDs."""
    import shutil
    from . import __version__
    from .sequence import sha256
    from .stage1 import verify, new_output, dump, now
    source, data = verify(args.results)
    target = Path(args.out).resolve()
    if target == source or source in target.parents or target in source.parents:
        raise ValueError('Report revision must be separate from the original checkpoint')
    root = new_output(target)
    original = json.loads((source/'seal.json').read_text())
    # Copy only checkpoint artifacts, never whole-genome scratch databases/workers.
    for name in original:
        destination = root/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, destination)
    if (source/'loci.decisions.tsv').exists():
        shutil.copyfile(source/'loci.decisions.tsv', root/'loci.decisions.tsv')
    manifest = json.loads((root/'manifest.json').read_text())
    manifest['presentation'] = dict(version=__version__, created=now(), source_checkpoint=str(source),
        source_evidence_sha256=sha256(source/'evidence.json'), retained_review_id=data['review_id'],
        renderer_sha256={name:sha256(Path(__file__).with_name(name)) for name in
                         ['stage1_report.py','stage1_review.txt','stage1_context_view.txt']},
        note='Presentation revision only; all search evidence, filters, coordinates and review IDs retained.')
    dump(root/'manifest.json', manifest)
    write_report(root, data)
    dump(root/'seal.json', {name:sha256(root/name) for name in original})
    print('Updated report: '+str(root/'report.html'))
    return 0
