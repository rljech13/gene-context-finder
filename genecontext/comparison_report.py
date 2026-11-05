"""Portable comparative HTML/SVG reports and protein/codon summaries."""
from collections import Counter, defaultdict
import hashlib
from html import escape

from .alignment import backtranslate, pairwise_stats
from .sequence import read_fasta, translate, write_fasta, write_tsv


PAIR_COLUMNS = ["sequence_a", "sequence_b", "paired_positions", "identical_positions", "identity", "coverage_a", "coverage_b"]


def protein_comparison(args, source, root, groups, queries, selected, runner, errors):
    summaries = []
    query_by_id = {q["id"]: q for q in queries}
    for group in groups:
        folder = root / "queries" / group["id"]
        folder.mkdir(parents=True)
        proteins, coding = {}, {}
        for row in selected:
            allowed = {"complete_predicted", "partial"} if args.include_partial else {"complete_predicted"}
            if row["group"] != group["id"] or row["status"] not in allowed:
                continue
            try:
                pid = row["model_id"]
                proteins[row["locus_id"]] = read_fasta(source / "sequences" / row["query"] / "protein.faa")[pid]
                cds = read_fasta(source / "sequences" / row["query"] / "cds.fna")[pid]
                coding[row["locus_id"]] = cds[:-3] if translate(cds).endswith("*") else cds
            except (OSError, ValueError, KeyError) as error:
                errors.append(dict(stage="protein_input", locus_id=row["locus_id"], error=str(error)))
        candidate_ids = set(proteins)
        for qid in group["queries"]:
            query = query_by_id[qid]
            if query["protein"]:
                proteins["reference_"+qid] = query["protein"]
                if query["coding"]:
                    cds = query["coding"]
                    coding["reference_"+qid] = cds[:-3] if translate(cds).endswith("*") else cds
        if not proteins:
            summaries.append(dict(group=group["id"], status="no_protein_models", candidate_proteins=0))
            continue
        try:
            if sum(map(len, proteins.values())) > getattr(args, "max_alignment_bases", 2000000):
                raise ValueError("Protein alignment exceeds --max-alignment-bases")
            path, aligned = folder / "proteins.input.faa", folder / "proteins.aligned.faa"
            write_fasta(path, proteins.items())
            if len(proteins) == 1:
                write_fasta(aligned, proteins.items())
            else:
                runner.run(["mafft", "--thread", args.threads, "--auto", "--amino", path], aligned)
            records = read_fasta(aligned)
            if set(records) != set(proteins) or any(seq.replace("-", "") != proteins[name] for name, seq in records.items()):
                raise ValueError("Protein alignment changed input content")
            pairwise = pairwise_stats(records, set("ACDEFGHIKLMNPQRSTVWY"))
            write_tsv(folder / "protein_pairs.tsv", pairwise, PAIR_COLUMNS)
            codons = backtranslate({k: v for k, v in records.items() if k in coding}, coding)
            write_fasta(folder / "cds.codon.aligned.fna", codons.items())
            write_tsv(folder / "cds_pairs.tsv", pairwise_stats(codons, set("ACGT")), PAIR_COLUMNS)
            profiles = []
            for position in range(len(next(iter(records.values())))):
                counts = Counter(records[k][position] for k in candidate_ids if records[k][position] in "ACDEFGHIKLMNPQRSTVWY")
                n = sum(counts.values())
                consensus, count = counts.most_common(1)[0] if counts else ("", 0)
                profiles.append(dict(alignment_position=position+1, observed_loci=n, eligible_loci=len(candidate_ids),
                                     coverage=n/len(candidate_ids) if candidate_ids else "", consensus=consensus,
                                     consensus_fraction=count/n if n >= 2 else ""))
            write_tsv(folder / "protein_profile.tsv", profiles, ["alignment_position", "observed_loci", "eligible_loci", "coverage", "consensus", "consensus_fraction"])
            _protein_svg(folder / "protein_conservation.svg", profiles)
            summaries.append(dict(group=group["id"], status="aligned" if candidate_ids else "query_references_only", candidate_proteins=len(candidate_ids), pairs=pairwise))
        except (OSError, ValueError, RuntimeError) as error:
            errors.append(dict(stage="protein", group=group["id"], error=str(error)))
            summaries.append(dict(group=group["id"], status="error", error=str(error), candidate_proteins=len(candidate_ids)))
    return summaries


def _svg(width, height, title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            '<style>text{font-family:system-ui,sans-serif;fill:#172033;font-size:12px}</style>',
            f'<text x="20" y="25" style="font-size:17px">{escape(title)}</text>']


def _color(group):
    return "#94a3b8" if not group else "#" + hashlib.sha256(group.encode()).hexdigest()[:6]


def _protein_svg(path, rows):
    shapes = _svg(1500, 240, "Protein conservation across recovered loci (query references excluded)")
    shapes.append('<path d="M50 180H1460M50 50V180" stroke="#64748b" fill="none"/>')
    bin_size = max(1, (len(rows)+1199)//1200)
    coverage_points = []
    for i in range(0, len(rows), bin_size):
        part = rows[i:i+bin_size]
        values = [r["consensus_fraction"] for r in part if r["consensus_fraction"] != ""]
        x, w = 50+1410*i/max(1, len(rows)), max(1, 1410*len(part)/max(1, len(rows)))
        height = 120*sum(values)/len(values) if values else 0
        shapes.append(f'<rect x="{x:.2f}" y="{180-height:.2f}" width="{w:.2f}" height="{height:.2f}" fill="#2563eb"/>')
        coverage = [r["coverage"] for r in part if r["coverage"] != ""]
        if coverage:
            coverage_points.append(f'{x:.2f},{180-120*sum(coverage)/len(coverage):.2f}')
    shapes.append('<polyline points="'+' '.join(coverage_points)+'" fill="none" stroke="#ea580c" stroke-width="2"/>')
    shapes += ['<text x="50" y="210">Alignment column →; blue: consensus fraction; orange: locus coverage. Blank consensus: fewer than two observed loci.</text>', '</svg>']
    path.write_text("\n".join(shapes))


def synteny_svg(path, title, loci, context):
    height = max(180, 105+len(loci)*100)
    shapes = _svg(1600, height, title+" — local gene context")
    if not loci:
        shapes += ['<text x="20" y="70">No selected locus with an accessible genome.</text>', '</svg>']
        path.write_text("\n".join(shapes)); return
    left = max(r["target_relative_start"] for r in loci)
    right = max(r["window_bp"]-r["target_relative_start"] for r in loci)
    scale = 1270/max(1, left+right)
    by_locus = defaultdict(list)
    for feature in context:
        by_locus[feature["locus_id"]].append(feature)
    previous = {}
    for i, row in enumerate(loci):
        y = 75+i*100
        x = lambda position: 300+(position-row["target_relative_start"]+left)*scale
        shapes += [f'<text x="20" y="{y}">{escape(row["genome_id"][:38])}</text>',
                   f'<text x="20" y="{y+18}" style="font-size:10px">{escape(row["locus_id"]+" · "+row["selection_status"])}</text>',
                   f'<path d="M{x(0):.2f} {y}H{x(row["window_bp"]):.2f}" stroke="#cbd5e1"/>']
        features = sorted(by_locus[row["locus_id"]], key=lambda f: f["relative_start"])
        counts = Counter(f.get("homology_group") for f in features)
        current = {}
        label_ends = [-1.0, -1.0]
        for feature in features:
            fam = feature.get("homology_group", "")
            start, end = max(0, feature["relative_start"]), min(row["window_bp"], feature["relative_end"])
            if start >= end:
                continue
            a, b = x(start), x(end)
            center = (a+b)/2
            if fam and counts[fam] == 1:
                current[fam] = (center, y)
                if fam in previous:
                    px, py = previous[fam]
                    shapes.append(f'<path d="M{px:.2f} {py+6}L{center:.2f} {y-6}" stroke="{_color(fam)}" stroke-width="3" opacity="0.25"/>')
            tip = min(7, (b-a)/2)
            points = (f'{a:.2f},{y-6} {b-tip:.2f},{y-6} {b:.2f},{y} {b-tip:.2f},{y+6} {a:.2f},{y+6}'
                      if feature["relative_strand"] == "+" else
                      f'{b:.2f},{y-6} {a+tip:.2f},{y-6} {a:.2f},{y} {a+tip:.2f},{y+6} {b:.2f},{y+6}')
            label = f'{feature["name"]} | {fam or "no protein anchor"} | {feature["evidence"]} | {feature["contig"]}:{feature["start"]}-{feature["end"]} ({feature["strand"]})'
            dash = ' stroke="#172033" stroke-dasharray="2 2"' if feature["evidence"] != "provided_gff3" else ''
            shapes.append(f'<polygon points="{points}" fill="{_color(fam)}"{dash}><title>{escape(label)}</title></polygon>')
            name = feature["name"][:20]
            label_start, label_end = center-len(name)*2.5, center+len(name)*2.5
            lane = next((j for j, end in enumerate(label_ends) if label_start > end+5), None)
            if lane is not None:
                shapes.append(f'<text x="{center:.2f}" y="{y+24+lane*12}" text-anchor="middle" style="font-size:9px">{escape(name)}</text>')
                label_ends[lane] = label_end
        a, b = x(row["target_relative_start"]), x(row["target_relative_end"])
        shapes.append(f'<rect x="{a:.2f}" y="{y-10}" width="{max(3,b-a):.2f}" height="20" fill="#dc2626"><title>{escape(row["model_id"]+" | "+row["anchor_kind"])}</title></rect>')
        shapes.append(f'<text x="{(a+b)/2:.2f}" y="{y-16}" text-anchor="middle" style="font-size:10px;fill:#b91c1c">{escape(row["target_name"][:30])}</text>')
        shapes.append(f'<text x="300" y="{y+48}" style="font-size:10px">{escape(row["contig"])}:{row["window_start"]}-{row["window_end"]} ({row["strand"]}); {escape(row["context_status"])}</text>')
        previous = current
    for tick in range(5):
        position = tick*(left+right)/4
        shapes.append(f'<text x="{300+position*scale:.2f}" y="{height-43}" text-anchor="middle" style="font-size:10px">{(position-left)/1000:.1f} kb</text>')
    shapes += [f'<text x="300" y="{height-27}" style="font-size:10px">Distance from target start; rows oriented with target. Native coordinates are shown per row.</text>',
               f'<text x="20" y="{height-10}">Red: target. Colours: protein homology groups, not proven orthology. Dashed: inferred gene.</text>', '</svg>']
    path.write_text("\n".join(shapes))


def dna_svg(path, title, loci, blocks, reference_id):
    reference = next((r for r in loci if r["locus_id"] == reference_id), None)
    shapes = _svg(1600, max(150, 100+len(loci)*55), title+" — local DNA blocks on reference coordinates")
    if reference is None:
        shapes += ['<text x="20" y="70">No reference locus is available in the requested reference genome.</text>', '</svg>']
        path.write_text("\n".join(shapes)); return
    scale = 1270/max(1, reference["window_bp"])
    for i, row in enumerate(loci):
        y = 65+i*55
        shapes += [f'<text x="20" y="{y+12}">{escape(row["genome_id"][:38])}</text>',
                   f'<rect x="300" y="{y}" width="1270" height="18" fill="#e2e8f0"/>']
        for block in blocks:
            if block["status"] != "accepted" or {block["locus_a"], block["locus_b"]} != {row["locus_id"], reference_id}:
                continue
            side = "a" if block["locus_a"] == reference_id else "b"
            start, end = block["start_"+side], block["end_"+side]
            color = "#9333ea" if block["relative_strand"] == "-" else "#16a34a"
            label = f'identity {block["identity"]:.3f}; {block["paired_bases"]} paired bases; relative strand {block["relative_strand"]}'
            shapes.append(f'<rect x="{300+start*scale:.2f}" y="{y}" width="{max(1,(end-start)*scale):.2f}" height="18" fill="{color}" opacity="{max(0.25,block["identity"]):.3f}"><title>{escape(label)}</title></rect>')
        if row["locus_id"] == reference_id:
            shapes.append(f'<text x="310" y="{y+14}" style="font-size:11px">Reference: {escape(reference_id)}</text>')
    height = max(150,100+len(loci)*55)
    for position in (reference["target_relative_start"], reference["target_relative_end"]):
        x = 300+position*scale
        shapes.append(f'<path d="M{x:.2f} 45V{height-65}" stroke="#dc2626" stroke-dasharray="3 3" opacity="0.6"/>')
    shapes.append(f'<text x="{300+reference["target_relative_start"]*scale:.2f}" y="43" style="fill:#dc2626;font-size:10px">Target boundaries</text>')
    for tick in range(5):
        position = tick*reference["window_bp"]/4
        shapes.append(f'<text x="{300+position*scale:.2f}" y="{height-43}" text-anchor="middle" style="font-size:10px">{(position-reference["target_relative_start"])/1000:.1f} kb</text>')
    shapes.append(f'<text x="300" y="{height-27}" style="font-size:10px">Reference distance from target start; target orientation.</text>')
    shapes += [f'<text x="20" y="{max(150,100+len(loci)*55)-10}">Green: accepted forward blocks; purple: reverse blocks; gray: no accepted alignment. Gray does not mean gene loss.</text>', '</svg>']
    path.write_text("\n".join(shapes))


def html_table(rows, columns):
    out = ['<table><thead><tr>'+''.join(f'<th>{escape(c.replace("_", " "))}</th>' for c in columns)+'</tr></thead><tbody>']
    for row in rows:
        values = []
        for key in columns:
            value = row.get(key, "")
            values.append(f'{value:.4f}' if isinstance(value, float) else str(value))
        out.append('<tr>'+''.join(f'<td>{escape(v)}</td>' for v in values)+'</tr>')
    out.append('</tbody></table>')
    return '\n'.join(out)


STYLE = '<style>body{font:15px system-ui;margin:30px auto;max-width:1600px;padding:0 20px;color:#172033}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{padding:8px;border:1px solid #ddd;text-align:left}th{background:#f1f5f9}img{width:100%;height:auto;border:1px solid #ddd}a{color:#1856b8}.note{padding:12px;background:#f1f5f9}h2{margin-top:30px}</style>'


def write_comparison_reports(root, groups, original, selected, observations, context, synteny, dna, blocks, protein, errors, args):
    index = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Comparative gene contexts</title>', STYLE,
             '<h1>Comparative gene contexts</h1>',
             f'<p>{len(groups)} target groups · {len(original["genomes"])} input genomes · ±{args.window:,} bp windows</p>',
             '<p class="note">Ranked loci are candidates. Multiple loci, isoforms and reference-supported predictions are retained in the tables. '
             'Local protein matches and conserved gene order support a comparison but do not establish individual-copy orthology. '
             'Missing annotation or alignment is not evidence of gene absence or regulatory activity.</p>',
             (f'<p class="note">Context windows use {"Tiberius on CPU" if args.predictor == "tiberius" else "AUGUSTUS"} for DNA-only gene prediction, followed by protein identification. Short-window skips and failures are shown in the locus status. '
              'Unidentified predictions remain in the context. Supplied GFF3 is a separate comparison/reference track.</p>'
              '<p><a href="context_predictions.tsv">All window predictions</a> · <a href="context_prediction.json">Predictor models and CPU provenance</a> · '
              '<a href="context_prediction_cds.fna">Predicted CDS</a> · '
              '<a href="context_prediction_proteins.faa">Predicted proteins</a> · '
              '<a href="provided_neighbors.tsv">Provided annotation neighbours</a></p>' if args.context_mode == "predict" else ''),
             '<p><a href="summary.tsv">All genomes/targets</a> · <a href="candidates.tsv">All candidates and ranks</a> · '
             '<a href="loci.tsv">Selected windows</a> · <a href="neighbors.tsv">Neighbours</a> · '
             '<a href="neighbor_proteins.tsv">Neighbour protein sources</a> · '
             '<a href="context_reference.json">Context reference database</a> · '
             '<a href="synteny.tsv">Synteny</a> · <a href="dna_pairs.tsv">DNA comparisons</a> · '
             '<a href="dna_categories.tsv">DNA by region</a> · <a href="manifest.json">Provenance</a></p><ul>']
    loci_by_id = {r["locus_id"]: r for r in selected if not r.get("comparison_error")}
    def with_genomes(rows):
        return [dict(r, genome_a=loci_by_id[r["locus_a"]]["genome_id"], genome_b=loci_by_id[r["locus_b"]]["genome_id"]) for r in rows]
    for group in groups:
        folder = root / "queries" / group["id"]
        folder.mkdir(parents=True, exist_ok=True)
        loci = [r for r in selected if r["group"] == group["id"] and not r.get("comparison_error")]
        requested = next((r for r in loci if r["genome_id"] == args.reference_genome), None) if args.reference_genome else (loci[0] if loci else None)
        reference_id = requested["locus_id"] if requested else ""
        group["reference_locus"] = reference_id
        synteny_svg(folder / "synteny.svg", group["name"], loci, context)
        dna_svg(folder / "dna_blocks.svg", group["name"], loci, blocks, reference_id)
        index.append(f'<li><a href="queries/{group["id"]}/report.html">{escape(group["name"])}</a> ({group["id"]})</li>')
        page = ['<!doctype html><html lang="en"><meta charset="utf-8">',
                f'<title>{escape(group["name"])}</title>', STYLE, f'<h1>{escape(group["name"])}</h1>',
                '<p><a href="../../report.html">All targets</a> · <a href="../../candidates.tsv">All candidate loci</a></p>',
                '<h2>Loci and evidence</h2>', html_table([r for r in observations if r["group"] == group["id"]],
                    ["genome_id", "status", "candidate_loci", "target_name", "model_status", "coverage", "identity", "model_id", "error"]),
                '<p class="note">One representative is shown per genome. It is chosen by query coverage × identity, then coding integrity. '
                'Ranked alternatives remain in candidates.tsv. This ranking is not an orthology assignment. '
                'Gene boundaries use GFF3 where the selected annotation supplies them; otherwise coding-model/DNA-hit boundaries are used. '
                'Windows are clipped at contig ends. Distances and coordinates are available in loci.tsv.</p>',
                '<h2>Local synteny</h2><img src="synteny.svg" alt="Gene order and orientation across input genomes">',
                html_table(loci, ["genome_id", "annotation_status", "context_status", "context_predictor", "context_predictor_model"]),
                html_table([r for r in context if r["group"] == group["id"]],
                           ["genome_id", "name", "gene_id", "homology_group", "start", "end", "strand", "evidence", "identification_status", "reference_identity", "reference_coverage"]),
                '<p>Neighbour families use reciprocal best protein matches within the collected context/reference catalog. '
                'They are provisional homology groups. At least two unique anchors are required to compare order; repeated families are excluded from that order score. '
                'Unannotated windows can contain additional unpredicted genes. Inferred names are borrowed from supporting proteins, '
                'not confirmed gene symbols. Automatic Swiss-Prot references are attributed to the UniProt Consortium (CC BY 4.0); '
                '<a href="../../context_reference.json">release, checksums and reference sources</a>.</p>',
                html_table(with_genomes([r for r in synteny if r["group"] == group["id"]]),
                           ["genome_a", "genome_b", "unique_anchors", "order_concordance", "relative_orientation_concordance", "status"]),
                '<h2>Genomic DNA conservation</h2><img src="dna_blocks.svg" alt="Local aligned DNA blocks">',
                '<p>All selected locus pairs are compared. The plot uses the stated reference; tables include every pair. '
                'Higher-scoring non-overlapping blocks contribute once. Identity excludes gaps and ambiguous bases; coverage counts paired canonical bases. '
                'Separate coding, between-CDS intron, target noncoding and flank results are in dna_categories.tsv. '
                'Input lowercase masking is preserved and the engines apply their own seed/low-complexity rules; no repeat library is downloaded or inferred. '
                'Repetitive matches and alignment failure need review before claiming conserved noncoding elements.</p>',
                html_table(with_genomes([r for r in dna if r["group"] == group["id"]]),
                           ["genome_a", "genome_b", "status", "accepted_blocks", "identity", "coverage_a", "coverage_b", "error"]),
                '<h2>Amino-acid conservation</h2>']
        psummary = next(r for r in protein if r["group"] == group["id"])
        page.append(html_table([psummary], ["status", "candidate_proteins", "error"]))
        if psummary["status"] in {"aligned", "query_references_only"}:
            page += ['<p><a href="proteins.aligned.faa">Protein alignment</a> · <a href="protein_pairs.tsv">Protein pairs</a> · '
                     '<a href="protein_profile.tsv">Protein profile</a> · <a href="cds.codon.aligned.fna">Codon alignment</a> · '
                     '<a href="cds_pairs.tsv">CDS pairs</a></p>',
                     '<img src="protein_conservation.svg" alt="Amino-acid conservation along the alignment">',
                     ('<p>Complete and clean partial coding models enter this alignment. ' if args.include_partial else
                      '<p>Only complete coding models enter this alignment. ') + 'Query references are included in the alignment, '
                     'but excluded from the conservation profile. Single-observation columns have no conservation score.</p>']
            labels = {r["locus_id"]: r["genome_id"] for r in loci}
            labels.update({"reference_"+q["id"]: "query: "+q["name"] for q in original["queries"]})
            page.append(html_table([dict(r, sequence_a=labels.get(r["sequence_a"],r["sequence_a"]),
                                         sequence_b=labels.get(r["sequence_b"],r["sequence_b"])) for r in psummary["pairs"]], PAIR_COLUMNS))
        page.append('</html>')
        (folder / "report.html").write_text('\n'.join(page))
    index += ['</ul>', html_table(observations, ["group_name", "genome_id", "status", "candidate_loci", "target_name", "model_status", "model_id", "error"])]
    combined_errors = original.get("errors", []) + errors
    if combined_errors:
        index += ['<h2>Search/comparison errors</h2><pre>'+escape(str(combined_errors))+'</pre>']
    index.append('</html>')
    (root / "report.html").write_text('\n'.join(index))
