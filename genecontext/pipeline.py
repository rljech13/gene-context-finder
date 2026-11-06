from __future__ import annotations

import csv
import html
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .alignment import align_sequences
from .report import write_report
from .sequence import (Intervals, coding_sequence, read_fasta, read_gff,
                       sha256, translate, union_length, write_fasta, write_tsv)
from .tools import Genome, Runner, tool_versions
from .protein_discovery import SCORE_FIELDS, translated_search, score_candidates


def load_genomes(args):
    rows = []
    if args.genomes:
        manifest = Path(args.genomes).resolve()
        with open(manifest, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"genome_id", "fasta"}.issubset(reader.fieldnames or []):
                raise ValueError("Genome TSV requires genome_id and fasta columns")
            for row in reader:
                if not any(row.values()):
                    continue
                for key in ("fasta", "gff3", "protein", "cds", "annotation_provenance"):
                    value = row.get(key, "") or ""
                    row[key] = str((manifest.parent / value).resolve()) if value else ""
                if row.get("tiberius_model"):
                    from .tiberius_backend import MODELS, model_name
                    value = model_name(row["tiberius_model"].strip())
                    row["tiberius_model"] = value if value in MODELS else str((manifest.parent / Path(value).expanduser()).resolve())
                rows.append(row)
    else:
        for path in sorted(Path(args.genome_dir).resolve().rglob("*")):
            if not path.is_file() or not re.search(r"\.(fa|fna|fasta)(\.gz)?$", path.name, re.I):
                continue
            stem = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", path.name, flags=re.I)
            gid = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)
            rows.append(dict(genome_id=gid, fasta=str(path), gff3=""))
    seen = set()
    for row in rows:
        gid = row["genome_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", gid):
            raise ValueError(f"Unsafe/empty genome_id {gid!r}; use letters, digits, . _ -")
        if gid in seen:
            raise ValueError(f"Duplicate genome_id {gid!r}; use an explicit unique-ID TSV")
        seen.add(gid)
        row.setdefault("gff3", "")
    if not rows:
        raise ValueError("No input genomes found")
    return rows


def load_queries(args):
    queries = []
    for i, (original, seq) in enumerate(read_fasta(args.query).items(), 1):
        coding, protein = "", ""
        if args.query_type == "protein":
            protein = seq.removesuffix("*")
            if not set(protein) <= set("ACDEFGHIKLMNPQRSTVWYX"):
                raise ValueError(f"{original}: invalid protein alphabet/internal stop")
        else:
            seq = seq.replace("U", "T")
            if not set(seq) <= set("ACGTRYMKSWBDHVN"):
                raise ValueError(f"{original}: invalid nucleotide alphabet")
            translated = translate(seq)
            inferred = (seq.startswith("ATG") and len(seq) % 3 == 0 and
                        translated.endswith("*") and "*" not in translated[:-1]
                        and set(seq) <= set("ACGT") and len(translated) >= 31)
            if args.query_type == "cds" or inferred:
                if len(seq) % 3 or "*" in translated.removesuffix("*") or "X" in translated:
                    raise ValueError(f"{original}: CDS must have unambiguous triplets, no internal stop")
                coding, protein = seq, translated.removesuffix("*")
        if protein and len(protein) < 30:
            raise ValueError(f"{original}: protein/CDS queries must encode >=30 amino acids")
        queries.append(dict(id=f"q{i:04d}", name=original, sequence=seq,
                            coding=coding, protein=protein, type=args.query_type,
                            inferred_cds=(args.query_type == "dna" and bool(coding))))
    return queries


def model_index(models):
    return Intervals(models, lambda m: m.contig, lambda m: m.start, lambda m: m.end)


def coding_intervals(model):
    """Normalise separate stop features for duplicate comparison, not extraction."""
    intervals = sorted([(s, e) for s, e, _ in model.cds] + model.stops)
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    # Preserve a different initial reading phase as a distinct partial model.
    first = sorted(model.cds, reverse=model.strand == "-")[0][2]
    return tuple(merged), first or 0


def validate_annotations(models, features, genome):
    for model in models:
        length = genome.length(model.contig)
        if any(s < 0 or e > length for s, e, _ in model.cds):
            raise ValueError(f"CDS coordinates outside source contig: {model.id}")
    for feature in features:
        if feature["end"] > genome.length(feature["contig"]):
            raise ValueError(f"Gene coordinates outside source contig: {feature['id']}")


def protein_search(args, runner, work, genome, queries, annotations):
    eligible = {q["id"]: q for q in queries if q["protein"]}
    if not eligible:
        return []
    proteins = work / "query.protein.faa"
    write_fasta(proteins, [(qid, q["protein"]) for qid, q in eligible.items()])
    index = work / "genome.mpi"
    runner.run(["miniprot", "-t", args.threads, "-d", index, genome.path], work / "index.stdout")
    gff = work / "miniprot.gff3"
    runner.run(["miniprot", "-t", args.threads, "-G", args.max_intron,
                "-N", args.max_hits, "-p", args.secondary_ratio, "--gff",
                "--outn", args.max_hits, "--outs", args.secondary_ratio,
                "--outc", args.min_coverage, index, proteins], gff)
    models, _ = read_gff(gff, source="miniprot")
    index = model_index(annotations)
    candidates, annotate = [], {}
    for model in models:
        target = model.attrs.get("Target", "").split()
        if len(target) < 3 or target[0] not in eligible:
            continue
        qid, begin, end = target[:3]
        coverage = (int(end)-int(begin)+1) / len(eligible[qid]["protein"])
        identity = float(model.attrs.get("Identity", 0))
        if coverage < args.min_coverage or identity < args.min_identity:
            continue
        evidence = dict(coverage=coverage, identity=identity, evidence="protein_to_genome")
        candidates.append((qid, model, evidence))
        for annotation in index.overlap(model.contig, model.start, model.end):
            if annotation.strand == model.strand and any(
                max(s, a) < min(e, b) for s, e, _ in model.cds for a, b, _ in annotation.cds):
                annotate[annotation.id] = annotation
    # Test overlapping annotated transcripts against the query instead of
    # attributing a long overlapping neighbour solely by envelope intersection.
    annotation_records, name_map = [], {}
    for i, model in enumerate(annotate.values(), 1):
        _, protein, status, _ = coding_sequence(model, genome)
        if protein and status != "disrupted":
            aid = f"a{i:06d}"
            annotation_records.append((aid, protein))
            name_map[aid] = model
    if annotation_records:
        subject = work / "annotation_candidates.faa"
        write_fasta(subject, annotation_records)
        hits = work / "annotation.blastp.tsv"
        runner.run(["blastp", "-query", proteins, "-subject", subject,
                    "-evalue", args.protein_evalue, "-max_hsps", "1",
                    "-outfmt", "6 qseqid sseqid pident qstart qend bitscore evalue"], hits)
        with open(hits) as handle:
            for line in handle:
                qid, aid, identity, start, end, score, evalue = line.split()
                coverage = (int(end)-int(start)+1)/len(eligible[qid]["protein"])
                identity = float(identity)/100
                if coverage >= args.min_coverage and identity >= args.min_identity:
                    candidates.append((qid, name_map[aid], dict(
                        coverage=coverage, identity=identity, bitscore=float(score),
                        evalue=float(evalue), evidence="annotation_and_blastp")))
    return candidates


def nucleotide_search(args, runner, work, genome, queries, annotations):
    if args.query_type == "protein":
        return [], []
    query_map = {q["id"]: q for q in queries}
    fasta = work / "query.dna.fna"
    write_fasta(fasta, [(q["id"], q["sequence"]) for q in queries])
    database = work / "blast_genome"
    # Do not parse/rename user contig IDs or impose BLAST's 50-character local-ID
    # limit. sseqid preserves the original FASTA first token in this database.
    if not Path(str(database)+".nin").exists():
        runner.run(["makeblastdb", "-in", genome.path, "-dbtype", "nucl",
                    "-out", database], work / "makeblastdb.stdout")
    hits = work / "blastn.tsv"
    # Increase subject cap to include every contig, not BLAST's default top 500.
    runner.run(["blastn", "-task", "blastn", "-query", fasta, "-db", database,
                "-num_threads", args.threads, "-evalue", args.dna_evalue,
                "-word_size", args.word_size, "-max_target_seqs", max(1, len(genome.index)),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore evalue"], hits)
    annotation_index, matches, hsps = model_index(annotations), defaultdict(list), []
    with open(hits) as handle:
        for line in handle:
            qid, contig, pid, length, qs, qe, ss, se, bits, evalue = line.split()
            if contig not in genome.index:
                raise ValueError(f"BLAST subject ID {contig!r} does not match the original FASTA IDs")
            if float(pid)/100 < args.dna_identity or int(length) < args.min_dna_hsp:
                continue
            qs, qe, ss, se = map(int, (qs, qe, ss, se))
            strand = "+" if (qe-qs)*(se-ss) >= 0 else "-"
            start, end = min(ss, se)-1, max(ss, se)
            hsp = dict(query=qid, contig=contig, start=start, end=end, strand=strand,
                       qstart=min(qs, qe)-1, qend=max(qs, qe), identity=float(pid)/100,
                       bitscore=float(bits), evalue=float(evalue))
            hsps.append(hsp)
            for model in annotation_index.overlap(contig, start, end):
                # Query genomic DNA can be supplied on either orientation. Both
                # genomic gene strands are valid; sequence evidence is retained.
                overlap = sum(max(0, min(end, b)-max(start, a)) for a, b, _ in model.cds)
                if overlap >= args.min_cds_overlap:
                    matches[(qid, model.id)].append(hsp)
    by_id, candidates = {m.id: m for m in annotations}, []
    for (qid, mid), evidence in matches.items():
        coverage = union_length([(h["qstart"], h["qend"]) for h in evidence])/len(query_map[qid]["sequence"])
        if coverage >= args.min_dna_coverage:
            candidates.append((qid, by_id[mid], dict(coverage=coverage,
                identity=max(h["identity"] for h in evidence),
                evidence="blastn_exon_overlap", coverage_kind="query_HSP_span_union")))
    return candidates, hsps


def merge_dna_regions(hsps, distance):
    groups, regions = defaultdict(list), []
    for hsp in hsps:
        groups[(hsp["query"], hsp["contig"], hsp["strand"])].append(hsp)
    for (qid, contig, strand), values in sorted(groups.items()):
        current = None
        for hsp in sorted(values, key=lambda h: h["start"]):
            if current is None or hsp["start"] > current["end"] + distance:
                current = dict(query=qid, contig=contig, strand=strand, start=hsp["start"],
                               end=hsp["end"], hsp_count=0)
                regions.append(current)
            current["end"] = max(current["end"], hsp["end"])
            current["hsp_count"] += 1
    return regions


def context_features(models, features):
    known = {f["id"] for f in features}
    derived = {}
    for model in models:
        key = (model.gene, model.contig, model.strand)
        if model.gene in known:
            continue
        if key not in derived:
            derived[key] = dict(id=model.gene, name=model.gene, contig=model.contig,
                                start=model.start, end=model.end, strand=model.strand)
        derived[key]["start"] = min(derived[key]["start"], model.start)
        derived[key]["end"] = max(derived[key]["end"], model.end)
    return features + list(derived.values())


def extract_model(args, genome, model, mid, qid, gid, evidence):
    cds, protein, status, issues = coding_sequence(model, genome)
    start = min([model.start] + [s for s, _ in model.stops])
    end = max([model.end] + [e for _, e in model.stops])
    size = genome.length(model.contig)
    low, high = max(0, start-args.flank), min(size, end+args.flank)
    left, right = (low, start), (end, high)
    up, down = (left, right) if model.strand == "+" else (right, left)
    row = dict(model_id=mid, genome_id=gid, query=qid, source=model.source,
               source_model=model.id, source_gene=model.gene, contig=model.contig,
               start=start, end=end, strand=model.strand, status=status,
               issues=";".join(issues), coverage=evidence["coverage"], identity=evidence["identity"],
               evidence=evidence["evidence"], cds_bp=len(cds), protein_aa=len(protein),
               evidence_json=json.dumps(evidence, separators=(",", ":")),
               cds_exons=json.dumps(model.cds, separators=(",", ":")),
               upstream_start=up[0], upstream_end=up[1], downstream_start=down[0],
               downstream_end=down[1], upstream_complete=up[1]-up[0] == args.flank,
               downstream_complete=down[1]-down[0] == args.flank,
               locus_start=low, locus_end=high)
    row.update({key: evidence.get(key, "") for key in SCORE_FIELDS})
    sequences = dict(cds=cds, protein=protein,
                     locus=genome.fetch(model.contig, low, high, model.strand),
                     upstream=genome.fetch(model.contig, *up, model.strand),
                     downstream=genome.fetch(model.contig, *down, model.strand))
    return row, sequences


def write_models_gff(path, records):
    with open(path, "w") as handle:
        handle.write("##gff-version 3\n")
        for row, model in records:
            mid = row["model_id"]
            handle.write(f"{model.contig}\tgenecontext\tmRNA\t{row['start']+1}\t{row['end']}\t.\t"
                         f"{model.strand}\t.\tID={mid};status={row['status']}\n")
            for start, end, phase in model.cds:
                handle.write(f"{model.contig}\t{model.source}\tCDS\t{start+1}\t{end}\t.\t"
                             f"{model.strand}\t{phase if phase is not None else '.'}\tParent={mid}\n")
            for start, end in model.stops:
                handle.write(f"{model.contig}\t{model.source}\tstop_codon\t{start+1}\t{end}\t.\t"
                             f"{model.strand}\t0\tParent={mid}\n")


def context_svg(path, row, neighbors):
    start, end = row["context_start"], row["context_end"]
    width, margin, scale = 1200, 35, 1130/max(1, end-start)
    text = html.escape(f"{row['model_id']} | {row['contig']}:{start}-{end} | native genomic coordinates")
    shapes = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="180" viewBox="0 0 {width} 180">',
              '<rect width="100%" height="100%" fill="white"/>',
              f'<text x="35" y="25" font-family="sans-serif" font-size="13">{text}</text>',
              '<path d="M35 85H1165" stroke="#999"/>']
    # Always draw the recovered model, even for an unannotated assembly.
    all_features = neighbors + [dict(start=row["start"], end=row["end"], strand=row["strand"],
                                    name=row["model_id"], id=row["model_id"], target=True)]
    for i, feature in enumerate(all_features):
        a, b = margin+(max(start, feature["start"])-start)*scale, margin+(min(end, feature["end"])-start)*scale
        y = 65 if feature.get("target") else 105 + (i % 2)*18
        tip = min(7, (b-a)/2)
        points = (f"{a},{y} {b-tip},{y} {b},{y+5} {b-tip},{y+10} {a},{y+10}"
                  if feature["strand"] != "-" else
                  f"{b},{y} {a+tip},{y} {a},{y+5} {a+tip},{y+10} {b},{y+10}")
        label = html.escape(str(feature.get("name", feature.get("id", ""))))
        color = "#2563eb" if feature.get("target") else "#64748b"
        shapes.append(f'<polygon points="{points}" fill="{color}"><title>{label}</title></polygon>')
        if feature.get("target"):
            shapes.append(f'<text x="{a:.2f}" y="57" font-family="sans-serif" font-size="11">{label}</text>')
    shapes.append('<text x="35" y="167" font-family="sans-serif" font-size="11">Blue: recovered model; gray: annotated neighbours. Gene envelopes shown; positions on source genome.</text></svg>')
    path.write_text("\n".join(shapes))


MODEL_COLUMNS = ["model_id", "genome_id", "query", "source", "source_model", "source_gene", "contig",
                 "start", "end", "strand", "status", "issues", "coverage", "identity", "evidence",
                 "evidence_json", "cds_bp", "protein_aa", "cds_exons", "locus_start", "locus_end",
                 "upstream_start", "upstream_end", "downstream_start", "downstream_end",
                 "upstream_complete", "downstream_complete", "context_start", "context_end", *SCORE_FIELDS]


def run(args):
    genomes, queries = load_genomes(args), load_queries(args)
    required = {"samtools"}
    if any(q["protein"] for q in queries):
        required.update({"miniprot", "blastp"})
        if args.protein_search == "tblastn":
            required.update({"tblastn", "makeblastdb"})
    if args.query_type != "protein":
        required.update({"blastn", "makeblastdb"})
    if not args.skip_alignments:
        required.add("mafft")
    versions = tool_versions(sorted(required))
    missing = [name for name, info in versions.items() if not info["path"]]
    if missing:
        raise ValueError("Missing tools: " + ", ".join(missing) + "; see environment.yml")
    root = Path(args.out).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    runner = Runner(root)
    manifest = dict(version=__version__, started=datetime.now(timezone.utc).isoformat(),
                    parameters=vars(args), tools=versions, query_sha256=sha256(args.query),
                    source_code_sha256={p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))},
                    queries=queries, genomes=[], coordinate_system="0-based half-open; GFF3 exports 1-based inclusive",
                    interpretation="Sequence-supported candidates; no automatic orthology, expression or gene-loss calls")
    rows, genome_summary, regions, context_rows = [], [], [], []
    sequences = {q["id"]: {key: [] for key in ("cds", "protein", "locus", "upstream", "downstream")}
                 for q in queries}
    for kinds in sequences.values():
        for key in ("dna_regions", "dna_region_upstream", "dna_region_downstream"):
            kinds[key] = []
    usable = {q["id"]: set() for q in queries}
    errors = []
    for entry in genomes:
        gid, genome = entry["genome_id"], None
        work = root / "genomes" / gid
        work.mkdir(parents=True)
        info = dict(entry)
        manifest["genomes"].append(info)
        print(f"[{gid}] searching all queries", file=sys.stderr, flush=True)
        try:
            info["fasta_sha256"] = sha256(entry["fasta"])
            if entry["gff3"]:
                info["gff3_sha256"] = sha256(entry["gff3"])
            genome = Genome(entry["fasta"], work / "cache", runner)
            annotations, features = read_gff(entry["gff3"]) if entry["gff3"] else ([], [])
            validate_annotations(annotations, features, genome)
            if args.protein_search == "tblastn":
                candidates, seed_regions = translated_search(args, runner, work, genome, queries, annotations)
            else:
                candidates = protein_search(args, runner, work, genome, queries, annotations)
                seed_regions = []
            dna_candidates, hsps = nucleotide_search(args, runner, work, genome, queries, annotations)
            candidates.extend(dna_candidates)
            candidates = score_candidates(args, runner, work, genome, queries, candidates)
            local_regions = merge_dna_regions(hsps, args.dna_merge_distance) + seed_regions
            for region in local_regions:
                region.setdefault("source", "blastn")
            raw_sequences = []
            for i, region in enumerate(local_regions, 1):
                region.update(genome_id=gid, region_id=f"{gid}.r{i:06d}")
                region.update(locus_start=max(0, region["start"]-args.flank),
                              locus_end=min(genome.length(region["contig"]), region["end"]+args.flank))
                left = (region["locus_start"], region["start"])
                right = (region["end"], region["locus_end"])
                up, down = (left, right) if region["strand"] == "+" else (right, left)
                region.update(upstream_start=up[0], upstream_end=up[1], downstream_start=down[0],
                              downstream_end=down[1], upstream_complete=up[1]-up[0] == args.flank,
                              downstream_complete=down[1]-down[0] == args.flank)
                raw_sequences.append((region["query"], region["region_id"], dict(
                    dna_regions=genome.fetch(region["contig"], region["start"], region["end"], region["strand"]),
                    dna_region_upstream=genome.fetch(region["contig"], *up, region["strand"]),
                    dna_region_downstream=genome.fetch(region["contig"], *down, region["strand"]))))
            write_fasta(work / "dna_regions.fna", [(r["region_id"], genome.fetch(
                r["contig"], r["start"], r["end"], r["strand"])) for r in local_regions])
            write_fasta(work / "dna_regions.with_flanks.fna", [(r["region_id"], genome.fetch(
                r["contig"], r["locus_start"], r["locus_end"], r["strand"])) for r in local_regions])
            # Exact exon-coordinate duplicates are the same model; transcript
            # isoforms and different loci remain separate. Prefer annotation.
            unique = {}
            for qid, model, evidence in sorted(candidates, key=lambda x: x[1].source != "annotation"):
                key = (qid, model.contig, model.strand, coding_intervals(model))
                support = dict(source=model.source, model=model.id, **evidence)
                if key not in unique:
                    evidence = dict(evidence, supporting_models=[])
                    unique[key] = (qid, model, evidence)
                unique[key][2]["supporting_models"].append(support)
            feature_index = Intervals(context_features(annotations, features), lambda f: f["contig"],
                                      lambda f: f["start"], lambda f: f["end"])
            local_rows, local_sequences, local_context, gff_records = [], [], [], []
            for i, (qid, model, evidence) in enumerate(sorted(unique.values(), key=lambda x:
                    (x[0], x[1].contig, x[1].start, x[1].end, x[1].id)), 1):
                mid = f"{gid}.{qid}.m{i:06d}"
                row, extracted = extract_model(args, genome, model, mid, qid, gid, evidence)
                row.update(context_start=max(0, row["start"]-args.context),
                           context_end=min(genome.length(model.contig), row["end"]+args.context))
                neighbors = feature_index.overlap(model.contig, row["context_start"], row["context_end"])
                for feature in sorted(neighbors, key=lambda f: (f["start"], f["end"], f["id"])):
                    local_context.append(dict(model_id=mid, genome_id=gid, **feature))
                context_svg(work / f"{mid}.context.svg", row, neighbors)
                local_rows.append(row)
                local_sequences.append((qid, mid, extracted, row))
                gff_records.append((row, model))
            write_models_gff(work / "models.gff3", gff_records)
            # Commit exports only after this genome completed successfully.
            rows.extend(local_rows)
            regions.extend(local_regions)
            context_rows.extend(local_context)
            for qid, rid, extracted in raw_sequences:
                for kind, seq in extracted.items():
                    sequences[qid][kind].append((rid, seq))
            for qid, mid, extracted, row in local_sequences:
                for kind, seq in extracted.items():
                    sequences[qid][kind].append((mid, seq))
                if row["status"] == "complete_predicted" or (args.include_partial and row["status"] == "partial"):
                    usable[qid].add(mid)
            for query in queries:
                found = [r for r in local_rows if r["query"] == query["id"]]
                raw = [r for r in local_regions if r["query"] == query["id"]]
                status = ("models_found" if found else "protein_seeds_no_cds_model" if any(r.get("source") == "tblastn" for r in raw) else "dna_hits_no_cds_model" if raw else "no_hit")
                genome_summary.append(dict(genome_id=gid, query=query["id"], query_name=query["name"],
                    status=status, models=len(found), complete=sum(r["status"] == "complete_predicted" for r in found),
                    partial=sum(r["status"] == "partial" for r in found),
                    disrupted=sum(r["status"] == "disrupted" for r in found), dna_regions=len(raw), error=""))
            info["status"] = "completed"
        except (OSError, ValueError, RuntimeError) as error:
            info.update(status="error", error=str(error))
            errors.append(dict(genome_id=gid, error=str(error)))
            print(f"[{gid}] ERROR: {error}; continuing with remaining genomes", file=sys.stderr)
            for query in queries:
                genome_summary.append(dict(genome_id=gid, query=query["id"], query_name=query["name"],
                    status="error", models=0, complete=0, partial=0, disrupted=0, dna_regions=0, error=str(error)))
        finally:
            if genome:
                genome.close()
    write_tsv(root / "models.tsv", rows, MODEL_COLUMNS)
    write_tsv(root / "genome_query_summary.tsv", genome_summary,
              ["genome_id", "query", "query_name", "status", "models", "complete", "partial", "disrupted", "dna_regions", "error"])
    write_tsv(root / "dna_regions.tsv", regions,
              ["region_id", "genome_id", "query", "contig", "start", "end", "strand", "hsp_count", "locus_start", "locus_end",
               "upstream_start", "upstream_end", "downstream_start", "downstream_end", "upstream_complete", "downstream_complete",
               "source", "seed_bitscore", "seed_query_coverage", "seed_identity", "ranking_score", "homology_support"])
    write_tsv(root / "neighbors.tsv", context_rows,
              ["model_id", "genome_id", "id", "name", "contig", "start", "end", "strand"])
    for qid, kinds in sequences.items():
        directory = root / "sequences" / qid
        directory.mkdir(parents=True)
        for kind, records in kinds.items():
            write_fasta(directory / f"{kind}.{'faa' if kind == 'protein' else 'fna'}", records)
    alignment_summary = []
    if not args.skip_alignments:
        try:
            alignment_summary = align_sequences(args, runner, root, queries, sequences, usable)
        except (OSError, RuntimeError, ValueError) as error:
            errors.append(dict(stage="alignment", error=str(error)))
    manifest.update(finished=datetime.now(timezone.utc).isoformat(),
                    status="completed_with_errors" if errors else "completed",
                    commands=runner.commands, errors=errors, alignment_summary=alignment_summary)
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    write_report(root, queries, genome_summary, rows, regions, errors)
    (root / "README.txt").write_text(
        "Start with genome_query_summary.tsv. Every input genome/query has a row.\n"
        "models.tsv contains sequence-supported candidate transcripts, not proven orthologs or copy counts.\n"
        "All TSV coordinates are native 0-based half-open; GFF3 is 1-based inclusive.\n"
        "CDS, proteins and oriented flanks are in sequences/; models.gff3 and context SVGs in genomes/.\n"
        "Only complete predicted CDS enter alignments unless --include-partial was requested.\n"
        "No-hit does not establish biological absence; absent annotations do not establish gene-free flanks.\n"
        "Noncoding DNA hits without an annotation/translated coding query remain DNA regions, not invented CDS.\n"
        "Protein-guided CDS alignments exclude terminal stops; CDS exports retain them when resolved.\n"
        "Flank similarity is exploratory; a global alignment does not establish regulatory function.\n")
    print(f"Results: {root}; {len(rows)} candidate models, {len(errors)} errors", file=sys.stderr)
    return 1 if errors else 0
