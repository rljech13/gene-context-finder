"""Small, dependency-free sequence and annotation primitives.

Internal coordinates are always zero-based, half-open, on the source genome.
GFF3 phases are checked in transcript order; internal exon bases are NEVER
discarded to account for a phase (a codon can span a splice junction).
"""
from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path, preserve_case=False):
    records, name, seq = {}, None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq) if preserve_case else "".join(seq).upper()
                words = line[1:].split()
                if not words:
                    raise ValueError(f"Empty FASTA identifier: {path}")
                name, seq = words[0], []
                if name in records:
                    raise ValueError(f"Duplicate FASTA identifier {name}: {path}")
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                seq.append("".join(line.split()))
    if name is not None:
        records[name] = "".join(seq) if preserve_case else "".join(seq).upper()
    if not records or any(not seq for seq in records.values()):
        raise ValueError(f"Empty FASTA or empty sequence: {path}")
    return records


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            if not seq:
                continue
            handle.write(f">{name}\n")
            for offset in range(0, len(seq), 80):
                handle.write(seq[offset:offset + 80] + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn",
                                       "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"))[::-1]


_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a+b+c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in _BASES for b in _BASES for c in _BASES), _AAS)}


def translate(seq):
    return "".join(CODE.get(seq[i:i+3], "X") for i in range(0, len(seq)-2, 3))


def union_length(intervals):
    total, last = 0, -1
    for start, end in sorted(intervals):
        total += max(0, end - max(start, last))
        last = max(last, end)
    return total


@dataclass
class Model:
    id: str
    contig: str
    strand: str
    cds: list[tuple[int, int, int | None]] = field(default_factory=list)
    stops: list[tuple[int, int]] = field(default_factory=list)
    source: str = "annotation"
    gene: str = ""
    attrs: dict = field(default_factory=dict)

    @property
    def start(self):
        return min(s for s, _, _ in self.cds)

    @property
    def end(self):
        return max(e for _, e, _ in self.cds)


def attributes(text):
    return {unquote(k): unquote(v) for item in text.split(";") if "=" in item
            for k, v in [item.split("=", 1)]}


def read_gff(path, source="annotation"):
    """Read GFF3, including shared CDS parents and CDS parented directly to genes.

    GTF is intentionally rejected: guessing its relationships loses isoforms.
    """
    metadata, parts, features = {}, {}, []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) != 9:
                raise ValueError(f"{path}:{number}: expected nine GFF3 columns")
            contig, _, kind, begin, end, _, strand, phase, raw = row
            attrs = attributes(raw)
            start, end = int(begin)-1, int(end)
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: invalid GFF3 interval")
            fid = attrs.get("ID", "")
            if fid:
                metadata[fid] = (contig, strand, attrs)
            if kind == "gene":
                features.append(dict(contig=contig, start=start, end=end,
                                     strand=strand, id=fid, name=attrs.get("Name", fid)))
            if kind not in {"CDS", "stop_codon"}:
                continue
            parents = attrs.get("Parent", "").split(",")
            if not parents[0] or strand not in {"+", "-"}:
                raise ValueError(f"{path}:{number}: CDS/stop requires Parent and strand")
            p = int(phase) if phase != "." else None
            if p not in {0, 1, 2, None}:
                raise ValueError(f"{path}:{number}: invalid phase")
            for parent in parents:
                model = parts.setdefault(parent, Model(parent, contig, strand, source=source))
                if model.contig != contig or model.strand != strand:
                    raise ValueError(f"{path}: transcript {parent} spans contigs/strands")
                if kind == "CDS":
                    model.cds.append((start, end, p))
                else:
                    model.stops.append((start, end))
    models = []
    for fid, model in parts.items():
        if not model.cds:
            continue
        _, _, attrs = metadata.get(fid, (None, None, {}))
        model.attrs = attrs
        model.gene = attrs.get("Parent", fid)
        model.cds.sort(reverse=model.strand == "-")
        models.append(model)
    return models, features


class Intervals:
    """Sorted interval index, with prefix maxima so nested long features survive."""
    def __init__(self, records, contig, start, end):
        self.groups = {}
        for record in records:
            self.groups.setdefault(contig(record), []).append(record)
        self.data = {}
        for key, values in self.groups.items():
            values.sort(key=start)
            starts, maxima, high = [], [], -1
            for record in values:
                starts.append(start(record))
                high = max(high, end(record))
                maxima.append(high)
            self.data[key] = (values, starts, maxima)
        self.start, self.end = start, end

    def overlap(self, contig, start, end):
        values, starts, maxima = self.data.get(contig, ([], [], []))
        i = bisect.bisect_left(starts, end)-1
        result = []
        while i >= 0 and maxima[i] > start:
            if self.end(values[i]) > start:
                result.append(values[i])
            i -= 1
        return result


def coding_sequence(model, genome):
    """Join source exons, validate phases, and include a separate terminal stop."""
    issues, pieces, used = [], [], 0
    ordered = sorted(model.cds, reverse=model.strand == "-")
    if union_length([(s, e) for s, e, _ in ordered]) != sum(e-s for s, e, _ in ordered):
        issues.append("overlapping_cds_features")
    for i, (start, end, phase) in enumerate(ordered):
        seq = genome.fetch(model.contig, start, end, model.strand)
        if i == 0:
            trim = phase or 0
            if trim:
                issues.append("partial_first_codon")
            seq = seq[trim:]
        elif phase is not None and phase != (3 - used % 3) % 3:
            issues.append("inconsistent_phase")
        pieces.append(seq)
        used += len(seq)
    cds = "".join(pieces)
    # GENCODE-style separate stop_codon features can overlap CDS in other sources.
    stop_extra = []
    for start, end in sorted(model.stops, reverse=model.strand == "-"):
        runs = [(start, end)]
        for cs, ce, _ in model.cds:
            runs = [(a, min(b, cs)) for a, b in runs if a < cs] + [
                (max(a, ce), b) for a, b in runs if b > ce]
        stop_extra.extend(runs)
    if stop_extra:
        stop_extra.sort(reverse=model.strand == "-")
        extra = "".join(genome.fetch(model.contig, s, e, model.strand) for s, e in stop_extra)
        terminal = (min(s for s, _, _ in ordered) if model.strand == "-"
                    else max(e for _, e, _ in ordered))
        adjacent = (max(e for s, e in stop_extra) == terminal if model.strand == "-"
                    else min(s for s, e in stop_extra) == terminal)
        if len(extra) <= 3 and adjacent and translate(cds + extra).endswith("*"):
            cds += extra
        else:
            issues.append("unresolved_separate_stop")
    protein = translate(cds)
    if len(cds) % 3:
        issues.append("non_triplet_cds")
    if any(base not in "ACGT" for base in cds):
        issues.append("ambiguous_cds")
    if "*" in protein.rstrip("*") or protein.endswith("**"):
        issues.append("internal_stop")
    if any(int(model.attrs.get(key, "0")) for key in ("Frameshift", "StopCodon")):
        issues.append("alignment_disruption")
    if not cds.startswith("ATG"):
        issues.append("missing_start")
    if not protein.endswith("*"):
        issues.append("missing_stop")
    bad = set(issues) - {"missing_start", "missing_stop", "partial_first_codon"}
    status = "disrupted" if bad else ("partial" if issues else "complete_predicted")
    return cds, protein.removesuffix("*"), status, sorted(set(issues))


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
