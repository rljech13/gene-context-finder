"""Protein-supported context annotation and exploratory local synteny anchors."""
from collections import Counter, defaultdict
from itertools import combinations
from urllib.parse import quote

from .sequence import coding_sequence, read_fasta, read_gff, union_length, write_fasta
from .tools import Genome, tool_versions


def _overlap(a, b):
    return union_length([(max(x, u), min(y, v)) for x, y in a for u, v in b if max(x, u) < min(y, v)])


def annotate_missing(args, root, selected, context, proteins, metadata, windows, runner, errors):
    if args.context_proteins:
        for name, seq in read_fasta(args.context_proteins).items():
            seq = seq.removesuffix("*")
            if len(seq) < 30 or not set(seq) <= set("ACDEFGHIKLMNPQRSTVWYX"):
                raise ValueError(f"Invalid context reference protein: {name}")
            pid = f"R{len(metadata)+1:07d}"
            proteins[pid] = seq
            metadata[pid] = dict(genome_id="external_reference", gene_id=name, name=name, evidence="supplied_reference")
    missing = [r for r in selected if r["annotation_status"] == "unannotated"]
    if not missing or not proteins:
        return
    tools = tool_versions(["blastx", "miniprot"])
    if any(not v["path"] for v in tools.values()):
        raise ValueError("Context inference requires blastx and miniprot")
    work = root / "context_inference"
    work.mkdir()
    reference = work / "reference.faa"
    reference_ids = set(proteins)
    write_fasta(reference, sorted(proteins.items()))
    runner.run(["makeblastdb", "-in", reference, "-dbtype", "prot", "-out", work / "reference_db"], work / "makeblastdb.stdout")
    query = work / "windows.fna"
    write_fasta(query, [(r["locus_id"], windows[r["locus_id"]]) for r in missing])
    screening = work / "blastx.tsv"
    runner.run(["blastx", "-query", query, "-db", work / "reference_db", "-num_threads", args.threads,
                "-evalue", "1e-5", "-max_target_seqs", len(proteins),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"], screening)
    hits = defaultdict(set)
    for line in screening.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 10 or fields[1] not in reference_ids:
            raise ValueError("Unexpected BLASTX screening output")
        if float(fields[2]) >= args.neighbor_identity * 100 and int(fields[3]) >= 20:
            hits[fields[0]].add(fields[1])
    native_records = []
    for row in missing:
        genome = None
        try:
            row["context_status"] = "reference_screen_no_hits"
            ids = hits[row["locus_id"]]
            if not ids:
                continue
            folder = work / row["locus_id"]
            folder.mkdir()
            refs = folder / "hits.faa"
            write_fasta(refs, [(pid, proteins[pid]) for pid in sorted(ids)])
            path = root / row["window_fasta"]
            predictions = folder / "miniprot.gff3"
            runner.run(["miniprot", "-t", args.threads, "-G", max(1000, row["window_bp"]),
                        "--gff", path, refs], predictions)
            models, _ = read_gff(predictions, source="miniprot_context")
            genome = Genome(path, folder / "cache", runner)
            accepted = []
            for model in sorted(models, key=lambda m: (m.start, m.end, m.id)):
                target = model.attrs.get("Target", "").split()
                if len(target) < 3 or target[0] not in ids:
                    continue
                coverage = (int(target[2])-int(target[1])+1) / len(proteins[target[0]])
                if coverage < args.neighbor_coverage or float(model.attrs.get("Identity", 0)) < args.neighbor_identity:
                    continue
                cds, protein, status, _ = coding_sequence(model, genome)
                if status == "disrupted" or len(protein) < 30 or "*" in protein:
                    continue
                intervals = [(a, b) for a, b, _ in model.cds]
                if row["relative_cds"] and _overlap(intervals, row["relative_cds"]) >= 0.8 * min(union_length(intervals), union_length(row["relative_cds"])):
                    continue
                # Multiple reference proteins can support the same inferred gene.
                if any(model.strand == m.strand and _overlap(intervals, [(a, b) for a, b, _ in m.cds]) >=
                       0.8 * min(len(cds), sum(b-a for a, b, _ in m.cds)) for m in accepted):
                    continue
                accepted.append(model)
                model_start = min([model.start] + [a for a, _ in model.stops])
                model_end = max([model.end] + [b for _, b in model.stops])
                if row["strand"] == "+":
                    start, end, strand = row["window_start"]+model_start, row["window_start"]+model_end, model.strand
                    native_cds = [(row["window_start"]+a, row["window_start"]+b, phase) for a, b, phase in model.cds]
                    native_stops = [(row["window_start"]+a, row["window_start"]+b) for a, b in model.stops]
                else:
                    start, end = row["window_end"]-model_end, row["window_end"]-model_start
                    strand = "-" if model.strand == "+" else "+"
                    native_cds = [(row["window_end"]-b, row["window_end"]-a, phase) for a, b, phase in model.cds]
                    native_stops = [(row["window_end"]-b, row["window_end"]-a) for a, b in model.stops]
                key = (row["genome_id"], row["contig"], start, end, strand)
                pid = next((p for p, m in metadata.items() if m.get("prediction_key") == key), None)
                if pid is None:
                    pid = f"N{len(metadata)+1:07d}"
                    ref = metadata[target[0]]
                    metadata[pid] = dict(prediction_key=key, genome_id=row["genome_id"], gene_id=pid,
                                         name=ref["name"], evidence="protein_supported_prediction",
                                         reference_protein=target[0])
                    proteins[pid] = protein
                    native_records.append((row["genome_id"], row["contig"], start, end, strand, pid, ref["name"], status, native_cds, native_stops))
                context.append(dict(locus_id=row["locus_id"], group=row["group"], genome_id=row["genome_id"],
                                    protein_id=pid, gene_id=pid, name=metadata[pid]["name"], contig=row["contig"],
                                    start=start, end=end, strand=strand, relative_start=model_start,
                                    relative_end=model_end, relative_strand=model.strand,
                                    evidence="protein_supported_prediction",
                                    clipped=model.start == 0 or model.end == row["window_bp"] or status == "partial"))
            row["context_status"] = "protein_supported_predictions" if accepted else "no_passing_context_predictions"
        except (OSError, ValueError, RuntimeError) as error:
            row["context_status"] = "context_inference_error"
            errors.append(dict(stage="context_inference", locus_id=row["locus_id"], error=str(error)))
        finally:
            if genome:
                genome.close()
    prediction_dir = root / "context_predictions"
    prediction_dir.mkdir()
    for gid in sorted({r[0] for r in native_records}):
        with open(prediction_dir / (gid + ".gff3"), "w") as handle:
            handle.write("##gff-version 3\n")
            for record_gid, contig, start, end, strand, pid, name, status, cds, stops in native_records:
                if record_gid != gid:
                    continue
                handle.write(f"{contig}\tminiprot_context\tgene\t{start+1}\t{end}\t.\t{strand}\t.\tID={pid};Name={quote(name, safe='')};status={status}\n")
                handle.write(f"{contig}\tminiprot_context\tmRNA\t{start+1}\t{end}\t.\t{strand}\t.\tID={pid}.t1;Parent={pid}\n")
                for a, b, phase in cds:
                    handle.write(f"{contig}\tminiprot_context\tCDS\t{a+1}\t{b}\t.\t{strand}\t{phase if phase is not None else '.'}\tParent={pid}.t1\n")
                for a, b in stops:
                    handle.write(f"{contig}\tminiprot_context\tstop_codon\t{a+1}\t{b}\t.\t{strand}\t0\tParent={pid}.t1\n")


def homology_groups(args, root, context, proteins, metadata, runner):
    if not proteins:
        return {}, []
    work = root / "neighbor_homology"
    work.mkdir()
    path = work / "proteins.faa"
    write_fasta(path, sorted(proteins.items()))
    runner.run(["makeblastdb", "-in", path, "-dbtype", "prot", "-out", work / "db"], work / "makeblastdb.stdout")
    hits = work / "blastp.tsv"
    runner.run(["blastp", "-query", path, "-db", work / "db", "-num_threads", args.threads,
                "-evalue", "1e-5", "-max_target_seqs", len(proteins),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"], hits)
    hsps = defaultdict(list)
    for line in hits.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 12 or fields[0] not in proteins or fields[1] not in proteins:
            raise ValueError("Unexpected neighbour BLASTP output")
        a, b = fields[:2]
        if metadata[a]["genome_id"] != metadata[b]["genome_id"]:
            hsps[a, b].append(fields)
    passing = {}
    for (a, b), values in hsps.items():
        chosen, ia, ib = [], [], []
        for h in sorted(values, key=lambda h: -float(h[9])):
            qa, qb = int(h[4])-1, int(h[5])
            sa, sb = int(h[6])-1, int(h[7])
            if any(max(qa, x) < min(qb, y) for x, y in ia) or any(max(sa, x) < min(sb, y) for x, y in ib):
                continue
            ia.append((qa, qb)); ib.append((sa, sb)); chosen.append(h)
        ca, cb = union_length(ia)/len(proteins[a]), union_length(ib)/len(proteins[b])
        identity = sum(float(h[2])*int(h[3]) for h in chosen)/sum(int(h[3]) for h in chosen)/100
        if min(ca, cb) >= args.neighbor_coverage and identity >= args.neighbor_identity:
            passing[a, b] = dict(identity=identity, coverage_a=ca, coverage_b=cb, bitscore=sum(float(h[9]) for h in chosen))
    best = defaultdict(list)
    for (a, b), metrics in passing.items():
        best[a, metadata[b]["genome_id"]].append((metrics["bitscore"], b))
    unique = {}
    for key, values in best.items():
        values.sort(reverse=True)
        if len(values) == 1 or values[0][0] > values[1][0] + 1e-6:
            unique[key] = values[0][1]
    parent = {p: p for p in proteins}
    def find(p):
        while parent[p] != p:
            parent[p] = parent[parent[p]]
            p = parent[p]
        return p
    matches = []
    for (a, b), metrics in sorted(passing.items()):
        if a >= b:
            continue
        reciprocal = (unique.get((a, metadata[b]["genome_id"])) == b and
                      unique.get((b, metadata[a]["genome_id"])) == a)
        if reciprocal:
            parent[find(b)] = find(a)
        matches.append(dict(protein_a=a, protein_b=b, **metrics,
                            reciprocal_best_in_local_catalog=reciprocal))
    groups, names = {}, {}
    for pid in sorted(proteins):
        representative = find(pid)
        names.setdefault(representative, f"H{len(names)+1:05d}")
        groups[pid] = names[representative]
    return groups, matches


def synteny_pairs(selected, context):
    by_locus = defaultdict(list)
    for feature in context:
        if feature.get("homology_group"):
            by_locus[feature["locus_id"]].append(feature)
    result = []
    for a, b in combinations(selected, 2):
        if a["group"] != b["group"]:
            continue
        fa = sorted(by_locus[a["locus_id"]], key=lambda f: f["relative_start"])
        fb = sorted(by_locus[b["locus_id"]], key=lambda f: f["relative_start"])
        ca, cb = Counter(f["homology_group"] for f in fa), Counter(f["homology_group"] for f in fb)
        common = ca.keys() & cb.keys()
        unique = {h for h in common if ca[h] == cb[h] == 1}
        aa = {f["homology_group"]: (i, f["relative_strand"]) for i, f in enumerate(fa) if f["homology_group"] in unique}
        bb = {f["homology_group"]: (i, f["relative_strand"]) for i, f in enumerate(fb) if f["homology_group"] in unique}
        pairs = list(combinations(sorted(unique), 2))
        order = sum((aa[x][0] < aa[y][0]) == (bb[x][0] < bb[y][0]) for x, y in pairs)/len(pairs) if pairs else ""
        orientation = sum(aa[h][1] == bb[h][1] for h in unique)/len(unique) if unique else ""
        status = ("context_evidence_missing" if any(r["context_status"] in {"no_annotation_evidence", "context_inference_error"} for r in (a, b))
                  else "insufficient_unique_anchors" if len(unique) < 2
                  else "conserved_local_order" if order == orientation == 1
                  else "different_local_order_or_orientation")
        result.append(dict(group=a["group"], locus_a=a["locus_id"], locus_b=b["locus_id"],
                           shared_groups=len(common), unique_anchors=len(unique), order_concordance=order,
                           relative_orientation_concordance=orientation, status=status))
    return result
