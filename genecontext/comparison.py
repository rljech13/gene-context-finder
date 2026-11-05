"""Comparison of user-supplied loci. No species, gene or database is built in."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

from .sequence import (coding_sequence, read_fasta, read_gff, sha256, union_length,
                       write_fasta, write_tsv)
from .tools import Genome, Runner, tool_versions


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_groups(queries, mapping=None):
    if not mapping:
        return [dict(id=q["id"], name=q["name"], queries=[q["id"]]) for q in queries]
    aliases = {q["id"]: q["id"] for q in queries}
    for q in queries:
        if q["name"] in aliases and aliases[q["name"]] != q["id"]:
            raise ValueError("Ambiguous query ID/name in group map")
        aliases[q["name"]] = q["id"]
    assigned, grouped = {}, {}
    for row in table(mapping):
        if not row.get("query") or not row.get("group") or row["query"] not in aliases:
            raise ValueError("Query-group TSV requires known query IDs/names and a nonempty group")
        qid, name = aliases[row["query"]], row["group"]
        if qid in assigned:
            raise ValueError(f"Duplicate group assignment: {qid}")
        assigned[qid] = name
        grouped.setdefault(name, []).append(qid)
    if set(assigned) != {q["id"] for q in queries}:
        raise ValueError("Query-group TSV must assign every input query exactly once")
    return [dict(id=f"g{i:04d}", name=name, queries=qids)
            for i, (name, qids) in enumerate(grouped.items(), 1)]


def _intervals(row):
    return [(int(a), int(b)) for a, b, *_ in json.loads(row.get("cds_exons") or "[]")]


def same_locus(a, b):
    if (a["contig"], a["strand"]) != (b["contig"], b["strand"]):
        return False
    if a.get("source") == b.get("source") == "annotation" and a.get("source_gene") != b.get("source_gene"):
        return False
    if a.get("source") == b.get("source") == "annotation" and a.get("source_gene"):
        return True
    ia, ib = _intervals(a), _intervals(b)
    if not ia or not ib:
        return (int(a["start"]), int(a["end"])) == (int(b["start"]), int(b["end"]))
    overlap = union_length([(max(x, u), min(y, v)) for x, y in ia for u, v in ib
                            if max(x, u) < min(y, v)])
    return overlap / min(union_length(ia), union_length(ib)) >= 0.8


def rank_key(row):
    status = {"complete_predicted": 0, "partial": 1, "disrupted": 2, "dna_region": 3}
    return (-float(row.get("coverage") or 0) * float(row.get("identity") or 0),
            status.get(row["status"], 4),
            -float(row.get("coverage") or 0), -int(row.get("protein_aa") or 0), row["model_id"])


def select_candidates(groups, genomes, models, regions, overrides=None):
    selected, candidates, observations = [], [], []
    override = {}
    group_aliases = {g["id"]: g["id"] for g in groups}
    for group in groups:
        if group["name"] in group_aliases and group_aliases[group["name"]] != group["id"]:
            raise ValueError("Ambiguous group ID/name in selection TSV")
        group_aliases[group["name"]] = group["id"]
    for row in table(overrides) if overrides else []:
        if row.get("group") not in group_aliases or not row.get("genome_id") or not row.get("model_id"):
            raise ValueError("Selection TSV requires group, genome_id, model_id")
        key = (group_aliases[row["group"]], row["genome_id"])
        if key in override:
            raise ValueError("Duplicate locus selection")
        override[key] = row["model_id"]
    used = set()
    for group in groups:
        for genome in genomes:
            gid, group_id = genome["genome_id"], group["id"]
            rows = [dict(r) for r in models if r["query"] in group["queries"] and r["genome_id"] == gid]
            # DNA regions are fallback observations, never fabricated coding models.
            if not rows:
                rows = [dict(r, model_id=r["region_id"], status="dna_region", source="blastn",
                             source_gene="", source_model=r["region_id"], coverage=0, identity=0,
                             protein_aa=0, cds_exons="[]") for r in regions
                        if r["query"] in group["queries"] and r["genome_id"] == gid]
            clusters = []
            for row in sorted(rows, key=rank_key):
                cluster = next((c for c in clusters if same_locus(row, c[0])), None)
                if cluster is None:
                    clusters.append([row])
                else:
                    cluster.append(row)
            representatives = [c[0] for c in clusters]
            choice = representatives[0] if representatives else None
            key = (group_id, gid)
            if key in override:
                choice = next((r for r in rows if r["model_id"] == override[key]), None)
                if choice is None:
                    raise ValueError(f"Selected model does not belong to {group_id}/{gid}: {override[key]}")
                used.add(key)
            status = ("user_selected" if key in override else "single_candidate_locus" if len(clusters) == 1
                      else "ranked_multiple_loci" if clusters else
                      "source_search_error" if genome.get("status") == "error" else "not_detected")
            observation = dict(group=group_id, group_name=group["name"], genome_id=gid,
                               status=status, candidate_loci=len(clusters),
                               model_id=choice["model_id"] if choice else "", error=genome.get("error", ""))
            observations.append(observation)
            for rank, cluster in enumerate(clusters, 1):
                for row in cluster:
                    candidates.append(dict(group=group_id, genome_id=gid, locus_rank=rank,
                                           model_id=row["model_id"], representative=row is cluster[0],
                                           selected=bool(choice and row["model_id"] == choice["model_id"]),
                                           contig=row["contig"], start=row["start"], end=row["end"],
                                           strand=row["strand"], status=row["status"],
                                           source=row.get("source", ""), source_gene=row.get("source_gene", ""),
                                           source_model=row.get("source_model", ""),
                                           query=row["query"], coverage=row.get("coverage", ""),
                                           identity=row.get("identity", "")))
            if choice:
                selected.append(dict(choice, group=group_id, selection_status=status,
                                     locus_id=f"L{len(selected)+1:06d}", candidate_loci=len(clusters)))
    if used != set(override):
        raise ValueError("Selection TSV contains an unknown genome or unused selection")
    return selected, candidates, observations


def orient_interval(locus, start, end):
    if locus["strand"] == "-":
        return locus["window_end"] - end, locus["window_end"] - start
    return start - locus["window_start"], end - locus["window_start"]


def category(locus, position):
    start, end = locus["target_relative_start"], locus["target_relative_end"]
    if position < start:
        return "upstream"
    if position >= end:
        return "downstream"
    if any(a <= position < b for a, b in locus["relative_cds"]):
        return "coding"
    if locus["relative_cds"] and min(a for a, _ in locus["relative_cds"]) <= position < max(b for _, b in locus["relative_cds"]):
        return "intron_between_cds"
    return "target_dna" if locus["status"] == "dna_region" else "noncoding_target"


def collect_windows(args, source, root, original, selected, runner, errors):
    from .pipeline import context_features, validate_annotations
    by_genome, context, proteins, metadata, windows = defaultdict(list), [], {}, {}, {}
    protein_ids = {}
    inputs = {g["genome_id"]: g for g in original["genomes"]}
    for row in selected:
        by_genome[row["genome_id"]].append(row)
    for gid, loci in by_genome.items():
        genome = None
        try:
            entry = inputs[gid]
            # Reuse the saved plain FASTA without creating another decompressed copy.
            cached = source / "genomes" / gid / "cache" / "genome.fa"
            fasta = cached if cached.exists() else Path(entry["fasta"])
            genome = Genome(fasta, root / "work" / gid, runner)
            gff = entry.get("gff3")
            models, features = read_gff(gff) if gff else ([], [])
            if gff and entry.get("gff3_sha256") and sha256(gff) != entry["gff3_sha256"]:
                raise ValueError("GFF3 changed since the source search")
            validate_annotations(models, features, genome)
            features = context_features(models, features)
            feature_by_id = {f["id"]: f for f in features}
            models_by_gene = defaultdict(list)
            for model in models:
                models_by_gene[model.gene].append(model)
            for row in loci:
                try:
                    gene = feature_by_id.get(row.get("source_gene", "")) if row.get("source") == "annotation" else None
                    start, end = (gene["start"], gene["end"]) if gene else (int(row["start"]), int(row["end"]))
                    low, high = max(0, start-args.window), min(genome.length(row["contig"]), end+args.window)
                    if high-low > args.max_window_bases:
                        raise ValueError(f"Window exceeds --max-window-bases ({high-low} bp)")
                    row.update(window_start=low, window_end=high, target_start=start, target_end=end,
                               anchor_kind="annotated_gene" if gene else "dna_region" if row["status"] == "dna_region" else "coding_model",
                               target_name=gene["name"] if gene else "DNA hit" if row["status"] == "dna_region" else "predicted coding model",
                               left_clipped=low > start-args.window, right_clipped=high < end+args.window,
                               annotation_status="provided_gff3" if gff else "unannotated",
                               relative_cds=[], context_status="provided_gff3" if gff else "no_annotation_evidence")
                    row["target_relative_start"], row["target_relative_end"] = orient_interval(row, start, end)
                    row["relative_cds"] = [orient_interval(row, a, b) for a, b in _intervals(row)]
                    sequence = genome.fetch(row["contig"], low, high, row["strand"], preserve_case=True)
                    path = root / "windows" / (row["locus_id"] + ".fna")
                    write_fasta(path, [(row["locus_id"], sequence)])
                    row.update(window_fasta=str(path.relative_to(root)), window_sha256=sha256(path),
                               softmasked_bases=sum(c.islower() for c in sequence), window_bp=len(sequence))
                    # Check the target against the original extraction before comparing surroundings.
                    if row["status"] != "dna_region":
                        from .sequence import Model
                        model = Model(row["model_id"], row["contig"], row["strand"],
                                      cds=[tuple(v) for v in json.loads(row["cds_exons"])])
                        # The original export can contain a separate terminal stop.
                        current, _, _, _ = coding_sequence(model, genome)
                        expected_path = source / "sequences" / row["query"] / "cds.fna"
                        expected = read_fasta(expected_path)[row["model_id"]]
                        if current != expected:
                            # A separate stop can extend the CDS envelope by up to three bases.
                            length = len(expected)-len(current)
                            extra = ""
                            if 0 < length <= 3:
                                terminal = min(a for a, _ in _intervals(row)) if row["strand"] == "-" else max(b for _, b in _intervals(row))
                                extra = (genome.fetch(row["contig"], terminal-length, terminal, "-") if row["strand"] == "-"
                                         else genome.fetch(row["contig"], terminal, terminal+length, "+"))
                            if current + extra != expected:
                                raise ValueError("Source target CDS no longer matches saved extraction")
                    else:
                        expected = read_fasta(source / "sequences" / row["query"] / "dna_regions.fna")[row["model_id"]]
                        if genome.fetch(row["contig"], int(row["start"]), int(row["end"]), row["strand"]) != expected:
                            raise ValueError("Source DNA hit no longer matches saved extraction")
                    windows[row["locus_id"]] = sequence
                    neighbours = [f for f in features if f["contig"] == row["contig"] and f["start"] < high and f["end"] > low]
                    for feature in neighbours:
                        # Target isoforms are not neighbouring-gene anchors.
                        target_cds = _intervals(row)
                        feature_cds = [(a, b) for m in models_by_gene[feature["id"]] for a, b, _ in m.cds]
                        shared = union_length([(max(a, c), min(b, d)) for a, b in target_cds for c, d in feature_cds
                                               if max(a, c) < min(b, d)])
                        target_like = (feature["strand"] == row["strand"] and target_cds and
                                       shared >= 0.8 * union_length(target_cds))
                        if feature["id"] == row.get("source_gene") or target_like:
                            continue
                        key = (gid, feature["contig"], feature["id"])
                        protein_id = protein_ids.get(key)
                        if protein_id is None:
                            protein_id = f"P{len(metadata)+1:07d}"
                            options = []
                            for model in models_by_gene[feature["id"]]:
                                _, protein, status, _ = coding_sequence(model, genome)
                                if status != "disrupted" and len(protein) >= 30 and "*" not in protein:
                                    options.append((status != "complete_predicted", -len(protein), model.id, protein, status))
                            metadata[protein_id] = dict(key=key, genome_id=gid, gene_id=feature["id"],
                                                       name=feature["name"], evidence="provided_gff3")
                            protein_ids[key] = protein_id
                            if options:
                                _, _, transcript, protein, coding_status = min(options)
                                proteins[protein_id] = protein
                                metadata[protein_id].update(transcript_id=transcript, coding_status=coding_status)
                        a, b = orient_interval(row, feature["start"], feature["end"])
                        context.append(dict(locus_id=row["locus_id"], group=row["group"], genome_id=gid,
                                            protein_id=protein_id, gene_id=feature["id"], name=feature["name"],
                                            contig=feature["contig"], start=feature["start"], end=feature["end"],
                                            strand=feature["strand"], relative_start=a, relative_end=b,
                                            relative_strand="+" if feature["strand"] == row["strand"] else "-",
                                            evidence="provided_gff3", clipped=feature["start"] < low or feature["end"] > high))
                except (OSError, ValueError, RuntimeError) as error:
                    row["comparison_error"] = str(error)
                    errors.append(dict(stage="window", genome_id=gid, model_id=row["model_id"], error=str(error)))
        except (OSError, ValueError, RuntimeError) as error:
            for row in loci:
                row["comparison_error"] = str(error)
            errors.append(dict(stage="context", genome_id=gid, error=str(error)))
        finally:
            if genome:
                genome.close()
    return context, proteins, metadata, windows


def compare(args):
    from . import __version__
    from .context_compare import annotate_missing, homology_groups, synteny_pairs, prune_reference_catalog
    from .dna_compare import compare_dna
    from .comparison_report import write_comparison_reports, protein_comparison
    source, root = Path(args.results).resolve(), Path(args.out).resolve()
    if source == root or source in root.parents or root in source.parents:
        raise ValueError("Comparison output must be separate from source search results")
    original = json.loads((source / "manifest.json").read_text())
    if args.context_mode == "predict":
        from .context_predict import validate_prediction_inputs
        validate_prediction_inputs(args, original["genomes"])
    groups = make_groups(original["queries"], args.query_groups)
    if args.reference_genome and args.reference_genome not in {g["genome_id"] for g in original["genomes"]}:
        raise ValueError("Unknown --reference-genome")
    versions = tool_versions(["samtools", "blastp", "makeblastdb", "mafft", args.dna_engine])
    if any(not item["path"] for item in versions.values()):
        raise ValueError("Missing comparison tools: " + ", ".join(k for k, v in versions.items() if not v["path"]))
    if root.exists() and any(root.iterdir()):
        raise ValueError("Comparison output must be new or empty")
    selected, candidates, observations = select_candidates(groups, original["genomes"], table(source / "models.tsv"),
                                                         table(source / "dna_regions.tsv"), args.select_loci)
    root.mkdir(parents=True, exist_ok=True)
    (root / "windows").mkdir()
    runner, errors = Runner(root), []
    manifest = dict(version=__version__, started=datetime.now(timezone.utc).isoformat(), stage="comparative_analysis",
                    source_results=str(source), source_manifest_sha256=sha256(source / "manifest.json"),
                    source_models_sha256=sha256(source / "models.tsv"), source_run_status=original.get("status"),
                    source_run_errors=original.get("errors", []), parameters=vars(args), tools=versions,
                    groups=groups, source_code_sha256={p.name: sha256(p) for p in Path(__file__).parent.glob('*.py')},
                    extra_inputs_sha256={k: sha256(v) for k, v in vars(args).items()
                                         if k in {"context_proteins", "query_groups", "select_loci"} and v})
    manifest["tools"].update(tool_versions(["blastx", "miniprot", "augustus"]))
    try:
        context, proteins, metadata, windows = collect_windows(args, source, root, original, selected, runner, errors)
        valid = [r for r in selected if not r.get("comparison_error")]
        try:
            if args.context_mode == "predict":
                from .context_predict import predict_contexts
                predict_contexts(args, root, original["genomes"], valid, context, proteins, metadata, runner, errors)
            else:
                annotate_missing(args, root, valid, context, proteins, metadata, windows, runner, errors)
        except (OSError, ValueError, RuntimeError) as error:
            errors.append(dict(stage="context_inference", error=str(error)))
            for row in valid:
                if args.context_mode == "predict" or row["annotation_status"] == "unannotated":
                    row["context_status"] = "context_inference_error"
            reference_receipt = root / "context_reference.json"
            if reference_receipt.exists():
                info = json.loads(reference_receipt.read_text())
                info.update(status="error", error=str(error))
                reference_receipt.write_text(json.dumps(info, indent=2) + "\n")
        finally:
            # Never send an entire failed/unused search library to all-vs-all BLASTP.
            prune_reference_catalog(proteins, metadata)
        if (root / "context_reference.json").exists():
            manifest["context_reference"] = json.loads((root / "context_reference.json").read_text())
        if (root / "context_prediction.json").exists():
            manifest["context_prediction"] = json.loads((root / "context_prediction.json").read_text())
        try:
            families, matches = homology_groups(args, root, context, proteins, metadata, runner)
        except (OSError, ValueError, RuntimeError) as error:
            errors.append(dict(stage="neighbor_homology", error=str(error)))
            families, matches = {}, []
        for feature in context:
            feature["homology_group"] = families.get(feature["protein_id"], "")
        synteny = synteny_pairs(valid, context)
        dna, blocks, categories = compare_dna(args, root, valid, windows, runner, errors)
        protein = protein_comparison(args, source, root, groups, original["queries"], valid, runner, errors)
        columns = ["group", "group_name", "genome_id", "status", "candidate_loci", "model_id", "target_name", "model_status", "coverage", "identity", "error"]
        for observation in observations:
            chosen = next((r for r in selected if r["group"] == observation["group"] and r["genome_id"] == observation["genome_id"]), None)
            if chosen and chosen.get("comparison_error"):
                observation.update(error=chosen["comparison_error"], status="comparison_error")
            if chosen:
                observation.update(target_name=chosen.get("target_name", ""), model_status=chosen["status"],
                                   coverage=chosen.get("coverage", ""), identity=chosen.get("identity", ""))
        write_tsv(root / "summary.tsv", observations, columns)
        write_tsv(root / "candidates.tsv", candidates, ["group", "genome_id", "locus_rank", "model_id", "representative", "selected", "contig", "start", "end", "strand", "status", "source", "source_gene", "source_model", "query", "coverage", "identity"])
        write_tsv(root / "loci.tsv", selected, ["locus_id", "group", "genome_id", "model_id", "selection_status", "candidate_loci", "target_name", "source", "source_gene", "source_model", "coverage", "identity", "protein_aa", "contig", "strand", "status", "window_start", "window_end", "target_start", "target_end", "anchor_kind", "left_clipped", "right_clipped", "window_bp", "softmasked_bases", "window_fasta", "window_sha256", "annotation_status", "context_status", "context_predictor_model", "comparison_error"])
        write_tsv(root / "neighbors.tsv", context, ["locus_id", "group", "genome_id", "gene_id", "name", "contig", "start", "end", "strand", "relative_start", "relative_end", "relative_strand", "protein_id", "homology_group", "evidence", "clipped", "reference_protein", "reference_identity", "reference_coverage", "prediction_coverage", "coding_status", "identification_status", "provided_gff_overlap"])
        write_tsv(root / "neighbor_proteins.tsv", [dict(protein_id=pid, **item,
                  protein_aa=len(proteins.get(pid, "")), protein_sha256=hashlib.sha256(proteins[pid].encode()).hexdigest() if pid in proteins else "")
                  for pid, item in sorted(metadata.items())],
                  ["protein_id", "genome_id", "gene_id", "name", "transcript_id", "coding_status", "evidence", "reference_protein", "reference_identity", "reference_coverage", "description", "protein_aa", "protein_sha256"])
        write_tsv(root / "neighbor_matches.tsv", matches, ["protein_a", "protein_b", "identity", "coverage_a", "coverage_b", "bitscore", "reciprocal_best_in_local_catalog"])
        write_tsv(root / "synteny.tsv", synteny, ["group", "locus_a", "locus_b", "shared_groups", "unique_anchors", "order_concordance", "relative_orientation_concordance", "status"])
        write_comparison_reports(root, groups, original, selected, observations, context, synteny, dna, blocks, protein, errors, args)
        manifest.update(context_reference_proteins=len(proteins), context_evidence=dict(
            provided=sum(r["evidence"] == "provided_gff3" for r in context),
            inferred=sum(r["evidence"] in {"protein_supported_prediction", "ab_initio_prediction"} for r in context)))
    except (OSError, ValueError, RuntimeError) as error:
        errors.append(dict(stage="comparison", error=str(error)))
        print(f"Comparison error: {error}", file=sys.stderr)
        if not (root / "report.html").exists():
            from html import escape
            (root / "report.html").write_text('<!doctype html><meta charset="utf-8"><h1>Comparison incomplete</h1><pre>'
                + escape(str(errors)) + '</pre><p><a href="manifest.json">Provenance and errors</a></p>')
    manifest.update(finished=datetime.now(timezone.utc).isoformat(), commands=runner.commands, errors=errors,
                    status="completed_with_errors" if errors or original.get("errors") else "completed")
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(f"Comparison: {root}; {len(errors)} new errors", file=sys.stderr)
    return 1 if errors or original.get("errors") else 0


def analyze(args):
    from .pipeline import run
    if args.context_mode == "predict":
        from .pipeline import load_genomes
        from .context_predict import validate_prediction_inputs
        validate_prediction_inputs(args, load_genomes(args))
    if args.word_size < 4 or args.protein_evalue <= 0 or args.dna_evalue <= 0:
        raise ValueError("Word size must be >=4 and E-values must be positive")
    missing = [name for name, info in tool_versions(["samtools", "blastp", "makeblastdb", "mafft", args.dna_engine]).items() if not info["path"]]
    if missing:
        raise ValueError("Missing analysis tools: " + ", ".join(missing))
    root = Path(args.out).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError("Analysis output must be new or empty")
    root.mkdir(parents=True, exist_ok=True)
    search_args = argparse.Namespace(**vars(args))
    search_args.out, search_args.skip_alignments = str(root / "search"), True
    search_args.flank = max(args.flank, args.window)
    search_args.context = max(args.context, args.window)
    search_code = run(search_args)
    compare_args = argparse.Namespace(**vars(args))
    compare_args.results, compare_args.out = str(root / "search"), str(root / "comparison")
    comparison_code = compare(compare_args)
    page = (root / "comparison" / "report.html").read_text()
    page = re.sub(r'href="(?!https?://|#)([^"]+)"', r'href="comparison/\1"', page)
    page = page.replace('<h1>', '<p><a href="search/report.html">All search candidates and search evidence</a></p><h1>', 1)
    (root / "report.html").write_text(page)
    return max(search_code, comparison_code)
