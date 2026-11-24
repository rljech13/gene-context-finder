"""Run against an installed wheel, from outside the source checkout.

Synthetic inputs exercise report generation and packaging, not biological claims.
The helper fixtures import production modules from the installed environment.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

from test_project_review import workspace, exported
from test_project import hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--runtime', type=Path)
    args = parser.parse_args()
    import genecontext
    checkout = Path(__file__).resolve().parents[1]
    if Path(genecontext.__file__).resolve().is_relative_to(checkout):
        raise SystemExit('Run from a clean installed environment outside the checkout')
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=False)
    project, source, choices, review = workspace(root)
    with project.locked():
        try:
            with project.locked():
                raise AssertionError('A second writer acquired the project lock')
        except ValueError as error:
            assert 'Another project operation' in str(error)
    before = hashes(source)
    review_root = project.root/review['directory']
    page = (review_root/'report.html').read_text()
    for marker in ['gcf-readable-model-labels-v1', 'gcf-alignment-exon-track-v1',
                   'function analysisComment(', 'syntenyReview',
                   'gcf-review-storage-v2', 'gcfReviewStorageLoader']:
        assert marker in page, marker
    assert 'DNMT3B' not in page
    export = root/'handoff.json'
    data = exported(project, review, choices, export)
    for analysis in data['analyses'].values():
        analysis['comments'] = {'alignment': 'Synthetic alignment note', 'tree': 'Synthetic tree note'}
    from test_project_review import save_export
    save_export(export, data)
    def cli(*argv):
        result = subprocess.run([sys.executable, '-m', 'genecontext', *map(str, argv)],
                                capture_output=True, text=True, check=True)
        return result.stdout
    cli('project', 'handoff', '--project', project.root, '--target', 'custom', '--file', export)
    cli('project', 'run', '--project', project.root, '--target', 'custom', '--stage', 'stage4')
    from genecontext.project_run import plan
    assert plan(project, 'custom', 'stage4')['action'] == 'reuse'
    result = project.current(project.load(), 'custom', 'stage4')
    final = Path(result['output'])
    evidence = json.loads((final/'evidence.json').read_text())
    assert evidence['trees'][0]['review_comments']['tree'] == 'Synthetic tree note'
    package = None
    if args.runtime:
        import os
        env = dict(os.environ, GCF_REVIEW_BIN_DIR=str(args.runtime.resolve()))
        package = root/'portable'
        subprocess.run([sys.executable, '-m', 'genecontext', 'review-package', '--results',
                        str(review_root), '--out', str(package)], env=env, check=True)
        with zipfile.ZipFile(package/'review_app.zip') as archive:
            assert archive.testzip() is None
            packed = archive.read('Gene_context_review/report.html').decode()
            assert 'gcf-review-storage-v2' in packed
            assert 'gcf-readable-model-labels-v1' in packed
    assert hashes(source) == before
    summary = dict(status='passed', synthetic_only=True, installed_module=str(genecontext.__file__),
                   report=str(review_root/'report.html'), final_report=str(final/'report.html'),
                   portable=str(package) if package else None, native_sources_unchanged=True,
                   report_sha256=hashlib.sha256(page.encode()).hexdigest())
    (root/'verification.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
