"""Sensitive translated-genome discovery before local coding reconstruction.

Short local matches are explicit candidates, never automatic gene/orthology calls.
"""
from collections import defaultdict
from dataclasses import replace
import json
from pathlib import Path
import sys

from .sequence import coding_sequence, read_gff, union_length, write_fasta


OUTFMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen"
SCORE_FIELDS = ["homology_support", "alignment_query_coverage", "alignment_identity", "protein_bitscore", "ranking_score"]


def hsp_summary(values, query_length):
    """Sum non-overlapping HSPs once, on both query and subject axes."""
    qa, sa, chosen = [], [], []
    for h in sorted(values, key=lambda h: -h["bitscore"]):
        q, s = (h["qstart"], h["qend"]), (h["start"], h["end"])
        if any(max(q[0], a) < min(q[1], b) for a, b in qa) or any(max(s[0], a) < min(s[1], b) for a, b in sa):
            continue
        qa.append(q); sa.append(s); chosen.append(h)
    length = sum(h["aligned_aa"] for h in chosen)
    return dict(coverage=union_length(qa)/query_length,
                identity=sum(h["identity"]*h["aligned_aa"] for h in chosen)/length if length else 0,
                aligned_aa=length, bitscore=sum(h["bitscore"] for h in chosen),
                evalue=min((h["evalue"] for h in chosen), default=1))


def blast_proteins(args, runner, folder, queries, models, genome):
    """Score coding proteins and self references with one consistent BLASTP pass."""
    folder.mkdir(parents=True, exist_ok=True)
    qmap = {q["id"]: q for q in queries if q["protein"]}
    proteins, model_map = [], {}
    for i, model in enumerate(models, 1):
        _, protein, _, _ = coding_sequence(model, genome)
        if not protein:
            continue
        sid = f"p{i:07d}"
        proteins.append((sid, protein.replace("*", "X")))
        model_map[sid] = model
    if not proteins or not qmap:
        return {}, {}
    proteins.extend(("self_"+qid, q["protein"]) for qid, q in qmap.items())
    query, subject, hits = folder/"queries.faa", folder/"models.faa", folder/"blastp.tsv"
    write_fasta(query, [(qid, q["protein"]) for qid, q in qmap.items()]); write_fasta(subject, proteins)
    runner.run(["blastp", "-query", query, "-subject", subject, "-num_threads", args.threads,
                "-evalue", args.protein_evalue, "-max_target_seqs", len(proteins), "-outfmt", OUTFMT], hits)
    grouped = defaultdict(list)
    for line in hits.read_text().splitlines():
        f = line.split("\t")
        if len(f) != 11 or f[0] not in qmap or (f[1] not in model_map and f[1] != "self_"+f[0]):
            if len(f) == 11 and f[1].startswith("self_"):
                continue
            raise ValueError("Unexpected candidate BLASTP output")
        grouped[f[0], f[1]].append(dict(qstart=int(f[4])-1, qend=int(f[5]), start=int(f[6])-1, end=int(f[7]),
                                      identity=float(f[2])/100, aligned_aa=int(f[3]), evalue=float(f[8]), bitscore=float(f[9])))
    selfs = {qid: hsp_summary(grouped[qid, "self_"+qid], len(q["protein"]))["bitscore"] for qid, q in qmap.items()}
    scores = {}
    for (qid, sid), hsps in grouped.items():
        if sid in model_map:
            scores[qid, model_map[sid].id] = hsp_summary(hsps, len(qmap[qid]["protein"]))
    return scores, selfs


def score_candidates(args, runner, work, genome, queries, candidates):
    """Avoid ranking a forced long alignment above stronger local homology."""
    eligible = {q["id"] for q in queries if q["protein"]}
    if not eligible:
        return candidates
    unique = {model.id: model for qid, model, _ in candidates if qid in eligible}
    scores, selfs = blast_proteins(args, runner, work/"candidate_scoring", queries, unique.values(), genome)
    result = []
    for qid, model, evidence in candidates:
        if qid not in eligible:
            result.append((qid, model, evidence)); continue
        score = scores.get((qid, model.id), {})
        coverage, identity = score.get("coverage", 0), score.get("identity", 0)
        support = ("extended_alignment" if coverage >= args.min_coverage and identity >= args.min_identity else
                   "local_alignment" if score else "no_blastp_alignment")
        result.append((qid, model, dict(evidence, homology_support=support,
            alignment_query_coverage=coverage, alignment_identity=identity,
            protein_bitscore=score.get("bitscore", 0),
            ranking_score=score.get("bitscore", 0)/selfs[qid] if selfs.get(qid) else 0)))
    return result


def translated_search(args, runner, work, genome, queries, annotations):
    missing = [q for q in queries if q["protein"]]
    if not missing:
        return [], []
    from .pipeline import model_index, merge_dna_regions
    folder = work/"protein_discovery"; folder.mkdir()
    receipt = dict(method="independent TBLASTn seeds; annotated CDS or dense local miniprot refinement",
                   queries=[q["id"] for q in missing], status="running", annotation_candidates=0, refinement_limit=args.rescue_max_regions)
    path = folder/"discovery.json"
    path.write_text(json.dumps(receipt, indent=2)+"\n")
    print(f"[protein discovery] TBLASTn: {', '.join(receipt['queries'])}", file=sys.stderr, flush=True)
    query = folder/"queries.faa"; write_fasta(query, [(q["id"], q["protein"]) for q in missing])
    database = work/"blast_genome"
    if not Path(str(database)+".nin").exists():
        runner.run(["makeblastdb", "-in", genome.path, "-dbtype", "nucl", "-out", database], folder/"makeblastdb.stdout")
    hits = folder/"tblastn.tsv"
    runner.run(["tblastn", "-query", query, "-db", database, "-num_threads", args.threads,
                "-evalue", args.protein_evalue, "-max_target_seqs", max(1, len(genome.index)), "-outfmt", OUTFMT], hits)
    hsps = []
    for line in hits.read_text().splitlines():
        f = line.split("\t")
        if len(f) != 11 or f[1] not in genome.index:
            raise ValueError("TBLASTn returned an unknown source contig or malformed hit")
        start, end = sorted((int(f[6]), int(f[7])))
        h = dict(query=f[0], contig=f[1], start=start-1, end=end, strand="+" if int(f[7]) >= int(f[6]) else "-",
                 qstart=int(f[4])-1, qend=int(f[5]), aligned_aa=int(f[3]), identity=float(f[2])/100,
                 bitscore=float(f[9]), evalue=float(f[8]), qlen=int(f[10]))
        if h["aligned_aa"] >= args.rescue_min_aa and h["identity"] >= args.min_identity:
            hsps.append(h)
    regions = merge_dna_regions(hsps, args.dna_merge_distance)
    for region in regions:
        own = [h for h in hsps if h["query"] == region["query"] and h["contig"] == region["contig"] and
               h["strand"] == region["strand"] and region["start"] <= h["start"] and h["end"] <= region["end"]]
        region.update(source="tblastn", seed_bitscore=sum(h["bitscore"] for h in own),
                      seed_query_coverage=union_length([(h["qstart"], h["qend"]) for h in own])/own[0]["qlen"],
                      seed_identity=sum(h["identity"]*h["aligned_aa"] for h in own)/sum(h["aligned_aa"] for h in own),
                      ranking_score=sum(h["bitscore"] for h in own)/own[0]["qlen"], homology_support="protein_seed_only")
    index = model_index(annotations)
    supported = defaultdict(dict)
    for h in hsps:
        for model in index.overlap(h["contig"], h["start"], h["end"]):
            overlap = sum(max(0, min(b,h["end"])-max(a,h["start"])) for a,b,_ in model.cds)
            if model.strand == h["strand"] and overlap >= args.min_cds_overlap:
                supported[h["query"]][model.id] = model
    models = {m.id: m for items in supported.values() for m in items.values()}
    scores, _ = blast_proteins(args, runner, folder/"annotation_support", missing, models.values(), genome)
    candidates = []
    for qid, items in supported.items():
        for mid, model in items.items():
            score = scores.get((qid, mid), {})
            if (score.get("coverage", 0) >= args.rescue_min_coverage and score.get("identity", 0) >= args.min_identity
                    and score.get("aligned_aa", 0) >= args.rescue_min_aa):
                candidates.append((qid, model, dict(score, evidence="tblastn_annotation_and_blastp")))
    receipt["annotation_candidates"] = len(candidates)
    # Regions with coding annotation support need no local reconstruction.
    unresolved = [region for region in regions if not any(qid == region["query"] and model.contig == region["contig"] and
                   model.strand == region["strand"] and any(max(a,region["start"]) < min(b,region["end"]) for a,b,_ in model.cds)
                   for qid,model,_ in candidates)]
    refined = []; selected = []
    for q in missing:
        selected.extend(sorted([x for x in unresolved if x["query"] == q["id"]], key=lambda x: -x["seed_bitscore"])[:args.rescue_max_regions])
    for i, region in enumerate(selected, 1):
        lo, hi = max(0,region["start"]-args.rescue_flank), min(genome.length(region["contig"]),region["end"]+args.rescue_flank)
        fasta, output = folder/f"region_{i:04d}.fna", folder/f"region_{i:04d}.gff3"
        write_fasta(fasta, [("window", genome.fetch(region["contig"],lo,hi,preserve_case=True))])
        one_query = folder/f"region_{i:04d}.query.faa"
        write_fasta(one_query, [(q["id"],q["protein"]) for q in missing if q["id"] == region["query"]])
        runner.run(["miniprot", "-k", 3, "-M", 0, "-n", 1, "-l", 3, "-t", args.threads, "-G", args.max_intron,
                    "-N", args.max_hits, "-p", args.secondary_ratio, "--gff", "--outn", args.max_hits,
                    "--outs", args.secondary_ratio, "--outc", args.rescue_min_coverage, fasta, one_query], output)
        local, _ = read_gff(output, source="miniprot_rescue")
        for model in local:
            if model.contig != "window" or model.strand != region["strand"]:
                continue
            native = replace(model, id=f"RES{i:04d}."+model.id, gene=f"RES{i:04d}."+model.gene,
                             contig=region["contig"], cds=[(a+lo,b+lo,p) for a,b,p in model.cds],
                             stops=[(a+lo,b+lo) for a,b in model.stops])
            if sum(max(0,min(b,region["end"])-max(a,region["start"])) for a,b,_ in native.cds) >= args.min_cds_overlap:
                refined.append((region["query"], native, dict(region)))
    refined_models = {m.id:m for _,m,_ in refined}
    scores, _ = blast_proteins(args, runner, folder/"refinement_support", missing, refined_models.values(), genome)
    for qid, model, region in refined:
        score = scores.get((qid,model.id), {})
        if (score.get("coverage",0) >= args.rescue_min_coverage and score.get("identity",0) >= args.min_identity
                and score.get("aligned_aa",0) >= args.rescue_min_aa):
            candidates.append((qid,model,dict(score,evidence="tblastn_and_local_miniprot", seed_start=region["start"],seed_end=region["end"])))
    receipt.update(status="completed", passing_hsps=len(hsps), seed_regions=len(regions),
                   refined_regions=len(selected), unrefined_regions=max(0,len(unresolved)-len(selected)), coding_candidates=len(candidates))
    path.write_text(json.dumps(receipt,indent=2)+"\n")
    return candidates, regions
