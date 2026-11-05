"""Portable local HTML overview; no web services or third-party assets."""
from html import escape


def write_report(root, queries, summary, models, regions, errors):
    names = {q["id"]: q["name"] for q in queries}
    document = ['<!doctype html><html lang="en"><meta charset="utf-8">',
                '<title>Gene and context search</title>',
                '<style>body{font:15px system-ui;margin:32px auto;max-width:1400px;padding:0 24px;color:#172033}'
                'table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #ddd;'
                'padding:7px;text-align:left}th{background:#f1f5f9}a{color:#1856b8}input{padding:10px;width:50%}'
                '.note{background:#f1f5f9;padding:12px}code{font-size:13px}h2{margin-top:32px}</style>',
                '<h1>Sequence-driven gene and context search</h1>',
                f'<p>{len(summary)} genome/query observations · {len(models)} candidate transcript models · '
                f'{len(regions)} DNA regions · {len(errors)} errors</p>',
                '<p class="note">Candidates are sequence-supported models. Transcript counts include isoforms. '
                'Local alignments can cover only a conserved domain and require family/orthology review. No-hit is not a biological absence call; sequence similarity does not establish orthology or expression. '
                'Coordinates are native, zero-based and half-open. Blue context arrows mark recovered models.</p>',
                '<p><a href="manifest.json">Run provenance</a> · <a href="models.tsv">Model table</a> · '
                '<a href="genome_query_summary.tsv">Genome summary</a> · '
                '<a href="dna_regions.tsv">DNA regions</a> · <a href="neighbors.tsv">Neighbours</a></p>',
                '<input id="filter" type="search" placeholder="Filter genome, query or model…" aria-label="Filter rows">',
                '<h2>All genomes and queries</h2><table><thead><tr><th>Genome</th><th>Query</th>'
                '<th>Status</th><th>Complete</th><th>Partial</th><th>Disrupted</th><th>DNA regions</th></tr></thead><tbody>']
    for row in summary:
        values = [row["genome_id"], names[row["query"]], row["status"], row["complete"],
                  row["partial"], row["disrupted"], row["dna_regions"]]
        title = escape(row.get("error", ""), quote=True)
        document.append(f'<tr class="data" title="{title}">' + ''.join(f'<td>{escape(str(v))}</td>' for v in values) + '</tr>')
    document.append('</tbody></table><h2>Query exports</h2><ul>')
    for query in queries:
        qid = query["id"]
        document.append(f'<li>{escape(query["name"])} ({qid}): '
                        f'<a href="sequences/{qid}/cds.fna">CDS</a> · '
                        f'<a href="sequences/{qid}/protein.faa">Protein</a> · '
                        f'<a href="sequences/{qid}/upstream.fna">Upstream</a> · '
                        f'<a href="sequences/{qid}/downstream.fna">Downstream</a> · '
                        f'<a href="sequences/{qid}/dna_regions.fna">DNA hits</a></li>')
    document.append('</ul><h2>Recovered coding models</h2><table><thead><tr><th>Model / context</th>'
                    '<th>Query</th><th>Location</th><th>Strand</th><th>Integrity</th><th>Source</th>'
                    '<th>Protein aa</th><th>Query coverage</th><th>Homology support</th></tr></thead><tbody>')
    for row in models:
        mid, gid = row["model_id"], row["genome_id"]
        link = f'<a href="genomes/{gid}/{mid}.context.svg">{escape(mid)}</a>'
        values = [names[row["query"]], f'{row["contig"]}:{row["start"]}-{row["end"]}', row["strand"],
                  row["status"] + (": " + row["issues"] if row["issues"] else ""), row["source"],
                  row["protein_aa"], f'{row["coverage"]:.1%}', row.get("homology_support", "")]
        document.append('<tr class="data"><td>' + link + '</td>' + ''.join(
            f'<td>{escape(str(v))}</td>' for v in values) + '</tr>')
    document.append('</tbody></table>')
    if errors:
        document.append('<h2>Errors</h2><ul>' + ''.join(
            f'<li>{escape(str(error))}</li>' for error in errors) + '</ul>')
    document.append('<script>document.getElementById("filter").addEventListener("input",function(){'
                    'const term=this.value.toLowerCase();document.querySelectorAll("tr.data").forEach('
                    'row=>row.hidden=!row.textContent.toLowerCase().includes(term));});</script></html>')
    (root / "report.html").write_text("\n".join(document))
