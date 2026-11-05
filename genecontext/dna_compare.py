"""Local genomic alignments with non-overlapping, base-counted summaries."""
from collections import Counter, defaultdict
from itertools import combinations

from .sequence import reverse_complement, write_tsv


def aligned_columns(text_a, text_b, start_a, start_b, step_a=1, step_b=1):
    if len(text_a) != len(text_b):
        raise ValueError("Different alignment row lengths")
    pa, pb = start_a, start_b
    for x, y in zip(text_a.upper(), text_b.upper()):
        ia, ib = (pa if x != "-" else None), (pb if y != "-" else None)
        if x != "-":
            pa += step_a
        if y != "-":
            pb += step_b
        yield ia, ib, x, y


def parse_blast(path):
    blocks = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 8:
            raise ValueError("Unexpected nucleotide BLAST output")
        a, b, c, d = map(int, fields[:4])
        blocks.append(dict(start_a=min(a, b)-1, end_a=max(a, b), start_b=min(c, d)-1, end_b=max(c, d),
                           step_a=1 if b >= a else -1, step_b=1 if d >= c else -1,
                           first_a=a-1, first_b=c-1, text_a=fields[4], text_b=fields[5], score=float(fields[6])))
    return blocks


def parse_lastz(path):
    blocks = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 8:
            raise ValueError("Unexpected LASTZ output")
        score, a, b, c, d, strand, ta, tb = fields
        a, b, c, d = map(int, (a, b, c, d))
        if strand not in {"+", "-"}:
            raise ValueError("Unexpected LASTZ strand")
        blocks.append(dict(start_a=a, end_a=b, start_b=min(c, d), end_b=max(c, d), step_a=1,
                           step_b=1 if strand == "+" else -1, first_a=a,
                           first_b=min(c, d) if strand == "+" else max(c, d)-1,
                           text_a=ta, text_b=tb, score=float(score)))
    return blocks


def summarize_blocks(blocks, seq_a, seq_b, min_block=50, min_identity=0.6):
    accepted_a, accepted_b, covered_a, covered_b = [], [], set(), set()
    output, columns = [], []
    for block in sorted(blocks, key=lambda b: (-b["score"], b["start_a"], b["start_b"])):
        b = dict(block)
        if not (0 <= b["start_a"] < b["end_a"] <= len(seq_a) and 0 <= b["start_b"] < b["end_b"] <= len(seq_b)):
            raise ValueError("Alignment coordinates outside context window")
        rows = list(aligned_columns(b["text_a"], b["text_b"], b["first_a"], b["first_b"], b["step_a"], b["step_b"]))
        for ia, ib, x, y in rows:
            if ia is not None:
                expected = seq_a[ia].upper() if b["step_a"] == 1 else reverse_complement(seq_a[ia].upper())
                if x != expected:
                    raise ValueError("DNA alignment does not roundtrip to first input sequence")
            if ib is not None:
                expected = seq_b[ib].upper() if b["step_b"] == 1 else reverse_complement(seq_b[ib].upper())
                if y != expected:
                    raise ValueError("DNA alignment does not roundtrip to second input sequence")
        if sum(x != "-" for _, _, x, _ in rows) != b["end_a"]-b["start_a"] or sum(y != "-" for _, _, _, y in rows) != b["end_b"]-b["start_b"]:
            raise ValueError("DNA alignment span and text disagree")
        paired = [(ia, ib, x, y) for ia, ib, x, y in rows if x in "ACGT" and y in "ACGT"]
        matches = sum(x == y for _, _, x, y in paired)
        identity = matches/len(paired) if paired else ""
        conflict = (any(max(b["start_a"], s) < min(b["end_a"], e) for s, e in accepted_a) or
                    any(max(b["start_b"], s) < min(b["end_b"], e) for s, e in accepted_b))
        reason = ("short_block" if len(paired) < min_block else "low_identity" if identity < min_identity
                  else "overlap_with_higher_scoring_block" if conflict else "accepted")
        b.update(paired_bases=len(paired), identical_bases=matches, identity=identity, status=reason,
                 relative_strand="+" if b["step_a"] == b["step_b"] else "-")
        if reason == "accepted":
            accepted_a.append((b["start_a"], b["end_a"]))
            accepted_b.append((b["start_b"], b["end_b"]))
            covered_a.update(ia for ia, _, _, _ in paired)
            covered_b.update(ib for _, ib, _, _ in paired)
            columns.extend(paired)
        output.append(b)
    n = len(columns)
    comparable_a = sum(c.upper() in "ACGT" for c in seq_a)
    comparable_b = sum(c.upper() in "ACGT" for c in seq_b)
    summary = dict(status="aligned_blocks" if n else "no_passing_blocks", accepted_blocks=sum(b["status"] == "accepted" for b in output),
                   paired_bases=n, identity=sum(x == y for _, _, x, y in columns)/n if n else "",
                   coverage_a=len(covered_a)/comparable_a if comparable_a else "",
                   coverage_b=len(covered_b)/comparable_b if comparable_b else "",
                   comparable_bases_a=comparable_a, comparable_bases_b=comparable_b)
    return summary, output, columns


def compare_dna(args, root, selected, windows, runner, errors):
    from .comparison import category
    folder = root / "dna"
    folder.mkdir()
    summaries, all_blocks, categories = [], [], []
    for a, b in combinations(selected, 2):
        if a["group"] != b["group"]:
            continue
        stem = a["locus_id"] + "__" + b["locus_id"]
        raw = folder / (stem + ".raw.tsv")
        common = dict(group=a["group"], locus_a=a["locus_id"], locus_b=b["locus_id"], engine=args.dna_engine)
        try:
            pa, pb = root / a["window_fasta"], root / b["window_fasta"]
            if args.dna_engine == "blastn":
                runner.run(["blastn", "-query", pa, "-subject", pb, "-task", "blastn", "-word_size", "11",
                            "-dust", "yes", "-lcase_masking", "-evalue", "1e-5", "-num_threads", args.threads,
                            "-outfmt", "6 qstart qend sstart send qseq sseq bitscore evalue"], raw)
                blocks = parse_blast(raw)
            else:
                runner.run(["lastz", pa, pb, "--strand=both",
                            "--format=general:score,zstart1,end1,zstart2+,end2+,strand2,text1,text2"], raw)
                blocks = parse_lastz(raw)
            sa, sb = windows[a["locus_id"]], windows[b["locus_id"]]
            stats, blocks, columns = summarize_blocks(blocks, sa, sb, args.dna_min_block, args.dna_min_identity)
            summaries.append(dict(common, **stats, error=""))
            all_blocks.extend(dict(common, block=i, **r) for i, r in enumerate(blocks, 1))
            for side, locus, seq, index in (("a", a, sa, 0), ("b", b, sb, 1)):
                totals = Counter(category(locus, i) for i, c in enumerate(seq) if c.upper() in "ACGT")
                counts, matches = Counter(), Counter()
                for column in columns:
                    label = category(locus, column[index])
                    counts[label] += 1
                    matches[label] += column[2] == column[3]
                for label, denominator in sorted(totals.items()):
                    count = counts[label]
                    categories.append(dict(common, side=side, locus_id=locus["locus_id"], category=label,
                                           comparable_bases=denominator, paired_bases=count,
                                           coverage=count/denominator, identity=matches[label]/count if count else ""))
        except (OSError, ValueError, RuntimeError) as error:
            errors.append(dict(stage="dna", **common, error=str(error)))
            summaries.append(dict(common, status="error", error=str(error)))
    write_tsv(root / "dna_pairs.tsv", summaries, ["group", "locus_a", "locus_b", "engine", "status", "accepted_blocks", "paired_bases", "identity", "coverage_a", "coverage_b", "comparable_bases_a", "comparable_bases_b", "error"])
    write_tsv(root / "dna_blocks.tsv", all_blocks, ["group", "locus_a", "locus_b", "engine", "block", "start_a", "end_a", "start_b", "end_b", "relative_strand", "score", "paired_bases", "identical_bases", "identity", "status"])
    write_tsv(root / "dna_categories.tsv", categories, ["group", "locus_a", "locus_b", "engine", "side", "locus_id", "category", "comparable_bases", "paired_bases", "coverage", "identity"])
    return summaries, all_blocks, categories
