from __future__ import annotations

from itertools import combinations

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
                runner.run(["mafft", "--thread", args.threads, "--auto",
                            "--amino" if kind == "protein" else "--nuc", source], target)
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
