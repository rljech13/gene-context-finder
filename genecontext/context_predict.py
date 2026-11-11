"""Predict every window from DNA first; identify the recovered proteins second."""
from collections import defaultdict
import csv
import json
from pathlib import Path
import re
from urllib.parse import quote

from .context_compare import _overlap
from .context_reference import add_references, swissprot
from .sequence import coding_sequence, read_gff, sha256, union_length, write_fasta, write_tsv
from .tools import Genome, tool_versions


def prediction_model(args, entry):
    if args.prediction_models:
        path = Path(args.prediction_models).resolve()
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"genome_id", "model"}.issubset(reader.fieldnames or []):
                raise ValueError("Prediction model TSV requires genome_id and model")
            mapping = {}
            for item in reader:
                gid, value = item["genome_id"].strip(), item["model"].strip()
                if not gid or not value or gid in mapping:
                    raise ValueError("Empty or duplicate prediction model assignment")
                if args.predictor == "tiberius":
                    from .tiberius_backend import MODELS, model_name
                    value = model_name(value)
                    value = value if value in MODELS else str((path.parent / Path(value).expanduser()).resolve())
                mapping[gid] = value
        if entry["genome_id"] in mapping:
            return mapping[entry["genome_id"]]
    if args.predictor == "augustus":
        return entry.get("augustus_species") or args.augustus_species
    from .tiberius_backend import model_name
    return model_name(entry.get("tiberius_model") or args.tiberius_model)


def validate_prediction_inputs(args, genomes):
    if args.predictor == "tiberius":
        from .tiberius_backend import runtime_info, runtime_python, ensure_model
        runtime = runtime_info(runtime_python(args))
        models = {name:ensure_model(name, args.predictor_cache, args.offline)
                  for name in sorted({prediction_model(args, g) for g in genomes})}
        return dict(runtime=runtime, models=models)
    if not tool_versions(["augustus"])["augustus"]["path"]:
        raise ValueError("Window gene prediction requires AUGUSTUS. Install the updated environment, or use --context-mode annotation for the legacy workflow.")
    missing = [g["genome_id"] for g in genomes if not prediction_model(args, g)]
    if missing:
        raise ValueError("Choose an AUGUSTUS trained model with --augustus-species MODEL, or an augustus_species column in genomes.tsv. Missing: " + ", ".join(missing))
    for entry in genomes:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", prediction_model(args, entry)):
            raise ValueError(f"Invalid AUGUSTUS model for {entry['genome_id']}")


def native_interval(row, start, end):
    if row["strand"] == "+":
        return row["window_start"] + start, row["window_start"] + end
    return row["window_end"] - end, row["window_end"] - start


def identify_proteins(args, root, predicted, references, metadata, runner):
    """BLASTP identification after prediction, retaining no-hit and tied models."""
    work = root / "context_identification"
    work.mkdir()
    provenance = dict(mode=args.context_reference, source="local_contexts", local_proteins=len(references),
                      status="preparing", predicted_proteins=len(predicted))
    receipt = root / "context_reference.json"
    receipt.write_text(json.dumps(provenance, indent=2) + "\n")
    if not predicted:
        provenance["status"] = "no_usable_predictions"
        receipt.write_text(json.dumps(provenance, indent=2) + "\n")
        return {}, references
    if args.context_proteins:
        provenance.update(source="supplied_and_local", supplied=add_references(
            args.context_proteins, references, metadata, "supplied_reference"))
    if args.context_reference == "swissprot" or (args.context_reference == "auto" and not references):
        provenance.update(source="automatic_swissprot", status="preparing_swissprot")
        receipt.write_text(json.dumps(provenance, indent=2) + "\n")
        path, info = swissprot(args.reference_cache, allow_download=not args.offline)
        provenance["swissprot"] = dict(info, **add_references(path, references, metadata, "swissprot_reference", skip_invalid=True))
    provenance.update(status="ready" if references else "no_reference_evidence", screening_proteins=len(references))
    receipt.write_text(json.dumps(provenance, indent=2) + "\n")
    if not references:
        return {}, references
    library, queries = work / "references.faa", work / "predicted.faa"
    write_fasta(library, sorted(references.items()))
    write_fasta(queries, sorted(predicted.items()))
    runner.run(["makeblastdb", "-in", library, "-dbtype", "prot", "-out", work / "db"], work / "makeblastdb.stdout")
    matches = work / "blastp.tsv"
    runner.run(["blastp", "-query", queries, "-db", work / "db", "-num_threads", args.threads,
                "-evalue", "1e-5", "-max_target_seqs", len(references),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"], matches)
    hsps = defaultdict(list)
    for line in matches.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 10 or fields[0] not in predicted or fields[1] not in references:
            raise ValueError("Unexpected protein identification BLASTP output")
        hsps[fields[0], fields[1]].append(fields)
    candidates = defaultdict(list)
    for (pid, rid), values in hsps.items():
        qa, sa, selected = [], [], []
        for fields in sorted(values, key=lambda f: -float(f[9])):
            q, s = (int(fields[4])-1, int(fields[5])), (int(fields[6])-1, int(fields[7]))
            if any(max(q[0], x) < min(q[1], y) for x, y in qa) or any(max(s[0], x) < min(s[1], y) for x, y in sa):
                continue
            qa.append(q); sa.append(s); selected.append(fields)
        cov_q, cov_s = union_length(qa) / len(predicted[pid]), union_length(sa) / len(references[rid])
        identity = sum(float(f[2]) * int(f[3]) for f in selected) / sum(int(f[3]) for f in selected) / 100
        if identity >= args.neighbor_identity and min(cov_q, cov_s) >= args.neighbor_coverage:
            candidates[pid].append(dict(reference_protein=rid, reference_identity=identity,
                                        reference_coverage=cov_s, prediction_coverage=cov_q,
                                        bitscore=sum(float(f[9]) for f in selected), name=metadata[rid]["name"]))
    assignments = {}
    for pid, options in candidates.items():
        options.sort(key=lambda r: (-r["bitscore"], -r["reference_identity"], r["reference_protein"]))
        best = options[0]
        ties = [r for r in options if abs(r["bitscore"] - best["bitscore"]) <= 1e-6]
        priority = {"swissprot_reference": 0, "supplied_reference": 1}
        preferred = min(priority.get(metadata[r["reference_protein"]]["evidence"], 2) for r in ties)
        preferred_ties = [r for r in ties if priority.get(metadata[r["reference_protein"]]["evidence"], 2) == preferred]
        best = preferred_ties[0]
        ambiguous = len({r["name"] for r in preferred_ties}) > 1
        assignments[pid] = dict(best, identification_status="ambiguous_best_matches" if ambiguous else "homology_assignment",
                                name="uncharacterized" if ambiguous else best["name"],
                                tied_references=",".join(r["reference_protein"] for r in ties))
    provenance.update(status="completed", assigned_proteins=len(assignments), no_passing_match=len(predicted)-len(assignments))
    receipt.write_text(json.dumps(provenance, indent=2) + "\n")
    return assignments, references


def predict_contexts(args, root, inputs, rows, context, proteins, metadata, windows, runner, errors):
    # GFF3 is an independent comparison/reference track, never prediction hints.
    provided = list(context)
    write_tsv(root / "provided_neighbors.tsv", provided,
              ["locus_id", "genome_id", "gene_id", "name", "contig", "start", "end", "strand", "protein_id"])
    context.clear()
    references = dict(proteins)
    for item in metadata.values():
        item["source_genome_id"] = item["genome_id"]
        item["genome_id"] = "annotation_reference:" + item["genome_id"]
    work, predictions = root / "context_prediction", root / "context_predictions"
    work.mkdir(); predictions.mkdir()
    entries = {g["genome_id"]: g for g in inputs}
    raw_records, cds_records, protein_records, pending = [], [], [], []
    configurations = {}
    native = defaultdict(list)
    for row in rows:
        row["context_status"] = "prediction_pending"
        row["context_predictor"] = args.predictor
        row["context_predictor_model"] = prediction_model(args, entries[row["genome_id"]])
    batch_models, provenance = {}, dict(predictor="augustus", method="DNA-only prediction before protein identification; no annotation/protein hints")
    if args.predictor == "tiberius":
        from .tiberius_backend import predict_batches
        batch_models, provenance = predict_batches(args, root, rows, windows, runner, errors)
    for row in rows:
        genome = None
        if row["context_status"] in {"context_inference_error", "prediction_window_too_short"}:
            continue
        try:
            folder = work / row["locus_id"]
            folder.mkdir()
            if args.predictor == "augustus":
                gff = folder / "augustus.gff3"
                runner.run(["augustus", "--species=" + row["context_predictor_model"], "--gff3=on",
                        "--strand=both", "--genemodel=partial", "--alternatives-from-sampling=false",
                        "--softmasking=on", "--UTR=off", root / row["window_fasta"]], gff)
                config_match = re.search(r"config directory (.+?)/?\s*\.\.\.", gff.read_text())
                if config_match:
                    config = Path(config_match.group(1).strip())
                    species = config / "species" / row["context_predictor_model"]
                    key = str(species)
                    if key not in configurations:
                        configurations[key] = {str(p.relative_to(config)): sha256(p)
                                           for p in sorted(species.rglob("*")) if p.is_file()}
                models, _ = read_gff(gff, source="augustus")
            else:
                models = batch_models.get(row["locus_id"], [])
            genome = Genome(root / row["window_fasta"], folder / "cache", runner)
            for model in models:
                if model.contig != row["locus_id"] or not 0 <= model.start < model.end <= row["window_bp"]:
                    raise ValueError("Predictor returned a model outside its window")
                pid = row["locus_id"] + "." + model.id
                dna, protein, status, issues = coding_sequence(model, genome)
                model_start = min([model.start] + [a for a, _ in model.stops])
                model_end = max([model.end] + [b for _, b in model.stops])
                native_start, native_end = native_interval(row, model_start, model_end)
                strand = model.strand if row["strand"] == "+" else ("-" if model.strand == "+" else "+")
                intervals = [(a, b) for a, b, _ in model.cds]
                target_overlap = (model.strand == "+" and bool(row["relative_cds"]) and
                                  _overlap(intervals, row["relative_cds"]) >= 0.8 *
                                  min(union_length(intervals), union_length(row["relative_cds"])))
                record = dict(protein_id=pid, locus_id=row["locus_id"], genome_id=row["genome_id"],
                              contig=row["contig"], start=native_start, end=native_end, strand=strand,
                              coding_status=status, issues=",".join(issues), target_overlap=target_overlap,
                              predictor=args.predictor, prediction_model=row["context_predictor_model"])
                raw_records.append(record)
                cds_records.append((pid, dna)); protein_records.append((pid, protein))
                native[row["genome_id"]].append((row, model, record))
                if target_overlap:
                    continue
                usable = status != "disrupted" and len(protein) >= 30 and not "*" in protein
                if usable:
                    proteins[pid] = protein
                    metadata[pid] = dict(genome_id=row["genome_id"], gene_id=pid, name="uncharacterized",
                                         evidence="ab_initio_prediction", coding_status=status)
                a, b = native_interval(row, model_start, model_end)
                validation = [p["gene_id"] for p in provided if p["locus_id"] == row["locus_id"] and p["strand"] == strand
                              and max(a, p["start"]) < min(b, p["end"])]
                feature = dict(locus_id=row["locus_id"], group=row["group"], genome_id=row["genome_id"],
                               protein_id=pid if usable else "", gene_id=pid, name="uncharacterized", contig=row["contig"],
                               start=a, end=b, strand=strand, relative_start=model_start, relative_end=model_end,
                               relative_strand=model.strand, evidence="ab_initio_prediction", coding_status=status,
                               identification_status="no_passing_match" if usable else "unusable_protein",
                               clipped=status == "partial", provided_gff_overlap=",".join(sorted(set(validation))))
                context.append(feature)
                pending.append((feature, record))
            row["context_status"] = "predicted_genes" if models else "no_predicted_genes"
        except (OSError, ValueError, RuntimeError) as error:
            row["context_status"] = "context_inference_error"
            errors.append(dict(stage="context_prediction", locus_id=row["locus_id"], error=str(error)))
        finally:
            if genome:
                genome.close()
    write_fasta(root / "context_prediction_cds.fna", cds_records)
    write_fasta(root / "context_prediction_proteins.faa", protein_records)
    (root / "context_prediction.json").write_text(json.dumps(dict(
        **provenance,
        models={r["genome_id"]: r["context_predictor_model"] for r in rows},
        configuration_files_sha256=configurations), indent=2) + "\n")
    # Save predictions before any network lookup; identification failure must
    # not erase the predicted, still unnamed genes from the context diagrams.
    assignments = {}
    try:
        predicted = {f["protein_id"]: proteins[f["protein_id"]] for f, _ in pending if f["protein_id"]}
        assignments, references = identify_proteins(args, root, predicted, references, metadata, runner)
    except (OSError, ValueError, RuntimeError) as error:
        errors.append(dict(stage="context_identification", error=str(error)))
        for feature, _ in pending:
            feature["identification_status"] = "identification_error"
        path = root / "context_reference.json"
        info = json.loads(path.read_text()) if path.exists() else {}
        info.update(status="error", error=str(error))
        path.write_text(json.dumps(info, indent=2) + "\n")
    for feature, record in pending:
        pid = feature["protein_id"]
        assignment = assignments.get(pid, {})
        feature.update(assignment)
        record.update(identification_status=feature["identification_status"], name=feature["name"], **{
            key: value for key, value in assignment.items() if key not in {"identification_status", "name"}})
        if pid:
            metadata[pid].update(assignment)
    # Only identified supporting references join the local neighbour catalogue.
    supporting = {r["reference_protein"] for r in assignments.values()}
    for rid in supporting:
        proteins[rid] = references[rid]
    for rid in list(metadata):
        if rid in references and rid not in supporting:
            metadata.pop(rid)
            proteins.pop(rid, None)
    columns = ["protein_id", "locus_id", "genome_id", "contig", "start", "end", "strand", "predictor", "prediction_model", "coding_status", "issues", "target_overlap",
               "name", "identification_status", "reference_protein", "reference_identity", "reference_coverage", "prediction_coverage", "bitscore", "tied_references"]
    write_tsv(root / "context_predictions.tsv", raw_records, columns)
    for gid, records in native.items():
        with open(predictions / (gid + ".gff3"), "w") as handle:
            handle.write("##gff-version 3\n")
            for row, model, record in records:
                pid = record["protein_id"]
                common = f'{row["contig"]}\t{args.predictor.upper() if args.predictor == "augustus" else "Tiberius"}\t'
                attrs = f'ID={pid};Name={quote(record.get("name", "uncharacterized"), safe="")};coding_status={record["coding_status"]};target_overlap={str(record["target_overlap"]).lower()}'
                handle.write(f'{common}gene\t{record["start"]+1}\t{record["end"]}\t.\t{record["strand"]}\t.\t{attrs}\n')
                handle.write(f'{common}mRNA\t{record["start"]+1}\t{record["end"]}\t.\t{record["strand"]}\t.\tID={pid}.t1;Parent={pid}\n')
                for a, b, phase in model.cds:
                    start, end = native_interval(row, a, b)
                    handle.write(f'{common}CDS\t{start+1}\t{end}\t.\t{record["strand"]}\t{phase if phase is not None else "."}\tParent={pid}.t1\n')
                for a, b in model.stops:
                    start, end = native_interval(row, a, b)
                    handle.write(f'{common}stop_codon\t{start+1}\t{end}\t.\t{record["strand"]}\t0\tParent={pid}.t1\n')
