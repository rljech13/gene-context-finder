from __future__ import annotations

from itertools import combinations
import sys

from .sequence import read_fasta, translate, write_fasta, write_tsv


def backtranslate(alignment, coding):
    result = {}
    for name, protein in alignment.items():
        cds, position, out = coding[name], 0, []
        for aa in protein:
            if aa == "-":
                out.append("---")
            else:
                codon = cds[position:position+3]
                if len(codon) != 3 or translate(codon) != aa:
                    raise ValueError(f"Codon/protein mismatch during backtranslation: {name}")
                out.append(codon)
                position += 3
        if position != len(cds):
            raise ValueError(f"Unused CDS bases during backtranslation: {name}")
        result[name] = "".join(out)
    return result


def pairwise_stats(alignment, alphabet):
    lengths = {len(seq) for seq in alignment.values()}
    if len(lengths) > 1:
        raise ValueError("Alignment rows have different lengths")
    rows = []
    for (a, sa), (b, sb) in combinations(alignment.items(), 2):
        paired = [(x, y) for x, y in zip(sa, sb) if x in alphabet and y in alphabet]
        na, nb = sum(x in alphabet for x in sa), sum(y in alphabet for y in sb)
        matches = sum(x == y for x, y in paired)
        rows.append(dict(sequence_a=a, sequence_b=b, paired_positions=len(paired),
                         identical_positions=matches,
                         identity=(matches/len(paired) if paired else ""),
                         coverage_a=(len(paired)/na if na else ""),
                         coverage_b=(len(paired)/nb if nb else "")))
    return rows


def align_sequences(args, runner, root, queries, sequences, usable):
    """Align all accepted candidates, not an arbitrary representative per genome."""
    summary = []
    for query in queries:
        qid = query["id"]
        output = root / "alignments" / qid
        output.mkdir(parents=True, exist_ok=True)
        records = {key: dict(values) for key, values in sequences[qid].items()}
        # Stop codons are excluded only in the codon alignment, not CDS exports.
        records["cds"] = {name: (cds[:-3] if translate(cds).endswith("*") else cds)
                          for name, cds in records["cds"].items() if name in usable[qid]}
        records["protein"] = {name: seq for name, seq in records["protein"].items()
                              if name in usable[qid]}
        if query["protein"]:
            ref = "reference_" + qid
            records["protein"][ref] = query["protein"]
            if query["coding"]:
                records["cds"][ref] = query["coding"]
                if translate(records["cds"][ref]).endswith("*"):
                    records["cds"][ref] = records["cds"][ref][:-3]
        # Flanks of partial or disrupted models remain exported for inspection,
        # but uncertain CDS boundaries do not define comparable aligned flanks.
        for key in ("upstream", "downstream", "locus"):
            records[key] = {name: seq for name, seq in records[key].items()
                            if name in usable[qid]}
        kinds = ["protein", "upstream", "downstream"]
        if query["type"] == "dna" and not query["protein"]:
            records["dna_regions"]["reference_" + qid] = query["sequence"]
            kinds.extend(["dna_regions", "dna_region_upstream", "dna_region_downstream"])
        if args.align_locus:
            kinds.append("locus")
        for kind in kinds:
            chosen = {name: seq for name, seq in records[kind].items() if seq}
            if not chosen:
                summary.append(dict(query=qid, kind=kind, status="no_eligible_sequences", n=0))
                continue
            if sum(map(len, chosen.values())) > args.max_alignment_bases:
                raise ValueError(f"{qid}/{kind}: exceeds --max-alignment-bases; inspect candidates "
                                 "or raise the limit explicitly")
            source = output / f"{kind}.input.fa"
            target = output / f"{kind}.aligned.fa"
            write_fasta(source, chosen.items())
            if len(chosen) == 1:
                write_fasta(target, chosen.items())
            else:
                options = (["--auto", "--amino"] if kind == "protein" else
                           ["--auto", "--nuc"] if getattr(args, "dna_alignment", "fast") == "auto" else
                           ["--retree", "2", "--maxiterate", "0", "--nuc"])
                runner.run(["mafft", "--thread", args.threads, *options, source], target)
            aligned = read_fasta(target)
            if set(aligned) != set(chosen):
                raise ValueError(f"MAFFT changed sequence identifiers: {target}")
            if any(aligned[key].replace("-", "") != seq for key, seq in chosen.items()):
                raise ValueError(f"MAFFT changed sequence content: {target}")
            alphabet = set("ACGT" if kind != "protein" else "ACDEFGHIKLMNPQRSTVWY")
            write_tsv(output / f"{kind}.pairwise.tsv", pairwise_stats(aligned, alphabet),
                      ["sequence_a", "sequence_b", "paired_positions", "identical_positions",
                       "identity", "coverage_a", "coverage_b"])
            summary.append(dict(query=qid, kind=kind, status="aligned" if len(chosen)>1
                                else "single_sequence", n=len(chosen)))
            if kind == "protein" and records["cds"]:
                # Protein-only query references have no invented reference CDS.
                nt = backtranslate({name: seq for name, seq in aligned.items()
                                    if name in records["cds"]}, records["cds"])
                write_fasta(output / "cds.codon.aligned.fna", nt.items())
                write_tsv(output / "cds.pairwise.tsv", pairwise_stats(nt, set("ACGT")),
                          ["sequence_a", "sequence_b", "paired_positions", "identical_positions",
                           "identity", "coverage_a", "coverage_b"])
    return summary


def realign(args):
    """Re-use saved extraction; write a new alignment run without editing it."""
    import csv
    import json
    from datetime import datetime, timezone
    from pathlib import Path
    from .sequence import sha256
    from .tools import Runner, tool_versions

    source, root = Path(args.results).resolve(), Path(args.out).resolve()
    if source == root or source in root.parents or root in source.parents:
        raise ValueError("Realignment output must be separate from the source result directory")
    parent = source / "manifest.json"
    original = json.loads(parent.read_text())
    with open(source / "models.tsv") as handle:
        models = list(csv.DictReader(handle, delimiter="\t"))
    versions = tool_versions(["mafft"])
    if not versions["mafft"]["path"]:
        raise ValueError("Missing tool: mafft")
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    queries = original["queries"]
    sequences, usable, hashes = {}, {}, {}
    kinds = ("cds", "protein", "locus", "upstream", "downstream", "dna_regions",
             "dna_region_upstream", "dna_region_downstream")
    for query in queries:
        qid = query["id"]
        sequences[qid], usable[qid] = {}, set()
        for kind in kinds:
            path = source / "sequences" / qid / f"{kind}.{'faa' if kind == 'protein' else 'fna'}"
            records = read_fasta(path) if path.exists() and path.stat().st_size else {}
            sequences[qid][kind] = list(records.items())
            if path.exists():
                hashes[str(path.relative_to(source))] = sha256(path)
    for row in models:
        if row["status"] == "complete_predicted" or (args.include_partial and row["status"] == "partial"):
            usable[row["query"]].add(row["model_id"])
    runner = Runner(root)
    manifest = dict(started=datetime.now(timezone.utc).isoformat(), stage="realignment",
                    source_results=str(source), source_manifest_sha256=sha256(parent),
                    source_models_sha256=sha256(source / "models.tsv"), sequence_sha256=hashes,
                    source_run_status=original.get("status"), parameters=vars(args), tools=versions,
                    source_run_errors=original.get("errors", []),
                    source_code_sha256={p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))
                                        if not p.name.startswith('.')})
    errors = []
    try:
        manifest["alignment_summary"] = align_sequences(args, runner, root, queries, sequences, usable)
    except (OSError, RuntimeError, ValueError) as error:
        errors.append(str(error))
        print(f"Alignment error: {error}", file=sys.stderr)
    manifest.update(finished=datetime.now(timezone.utc).isoformat(), commands=runner.commands,
                    errors=errors, status="failed" if errors else "completed")
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (root / "README.txt").write_text(
        "Realignment of saved extraction. Original search/model tables are unchanged.\n"
        f"Source results: {source}\n"
        "Consult both manifests; successful realignment does not erase source search errors.\n")
    print(f"Alignments: {root}; {len(errors)} errors", file=sys.stderr)
    return 1 if errors else 0
