"""A compact offline overview; it cannot start jobs or approve models."""
from html import escape
from pathlib import Path

from .project import STAGES
from .project_run import status_rows


def write_report(project):
    rows = status_rows(project)
    body = []
    for tid, target in project.targets.items():
        cells = []
        for stage in STAGES:
            row = next(r for r in rows if r['target'] == tid and r['stage'] == stage)
            label = escape(row['status'].replace('_', ' '))
            report = Path(row['report']) if row.get('report') else None
            if report and report.is_file():
                label = '<a href="'+escape(report.as_uri(), quote=True)+'">'+label+'</a>'
            cells.append('<td>'+label+('<small>'+escape(row['detail'])+'</small>' if row['detail'] else '')+'</td>')
        body.append('<tr><th>'+escape(target['gene'])+'<small>'+escape(tid)+'</small></th>'+''.join(cells)+'</tr>')
    html = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gene Context study</title><style>body{font:15px/1.5 system-ui;color:#233a3e;background:#f4f7f6;max-width:1200px;margin:40px auto;padding:20px}table{border-collapse:collapse;width:100%;background:white}th,td{text-align:left;padding:14px;border-bottom:1px solid #d9e3e0;vertical-align:top}small{display:block;color:#52656a;margin-top:5px}a{color:#125d64}code{background:#e4edeb;padding:2px 5px}.scroll{overflow:auto}</style>
<h1>__NAME__</h1><p>Search loci → review models → compare sequences → final dossier.</p>
<p>Each STOP requires an explicit review. Opening this overview starts no calculations and approves no models. Historical attempts and byte-exact decisions remain in the project ledger.</p>
<div class="scroll"><table><thead><tr><th>Gene / group</th><th>STOP 1 · loci</th><th>STOP 2 · models</th><th>STOP 3 · comparison</th><th>STOP 4 · dossier</th></tr></thead><tbody>__ROWS__</tbody></table></div>
<p>Refresh with <code>gene-context project report --project /path/to/project.json</code>. Server checkpoints should be viewed on their owning host or through its existing viewer.</p></html>'''.replace('__NAME__', escape(project.data['name'])).replace('__ROWS__', ''.join(body))
    path = project.root/'project-report.html'
    if path.exists() and '<title>Gene Context study</title>' not in path.read_text():
        raise ValueError('Refusing to replace an unrelated project-report.html')
    path.write_text(html)
    return dict(report=str(path), calculations_started=False)
