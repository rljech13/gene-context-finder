from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path


class Runner:
    def __init__(self, root):
        self.root = Path(root)
        self.commands = []

    def run(self, args, output, log=None):
        args = [str(arg) for arg in args]
        output = Path(output)
        log = Path(log) if log else output.with_name(output.name + ".log")
        output.parent.mkdir(parents=True, exist_ok=True)
        entry = dict(argv=args, output=str(output), stderr=str(log), returncode=None)
        self.commands.append(entry)
        with open(output, "w") as out, open(log, "w") as err:
            result = subprocess.run(args, stdout=out, stderr=err)
        entry["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError(f"{args[0]} exited {result.returncode}; see {log}")


def tool_versions(names):
    versions = {}
    for name in names:
        resolved = shutil.which(name)
        if not resolved:
            versions[name] = dict(path=None, version="missing")
            continue
        flag = "-version" if name in {"blastn", "blastp", "makeblastdb"} else "--version"
        result = subprocess.run([resolved, flag], capture_output=True, text=True,
                                errors="replace")
        lines = (result.stdout + result.stderr).splitlines()
        versions[name] = dict(path=resolved, version="\n".join(lines[:3]))
    return versions


class Genome:
    """Own a cached FASTA index; never create indices beside the user's inputs."""
    def __init__(self, source, cache, runner):
        self.path = Path(cache) / "genome.fa"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        source = Path(source).resolve()
        with open(source, "rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
        if compressed:
            with gzip.open(source, "rb") as inp, open(self.path, "wb") as out:
                shutil.copyfileobj(inp, out)
        else:
            self.path.symlink_to(source)
        runner.run(["samtools", "faidx", self.path], self.path.parent / "faidx.stdout")
        self.index = {}
        with open(str(self.path) + ".fai") as handle:
            for line in handle:
                name, length, offset, bases, width = line.rstrip().split("\t")[:5]
                if name in self.index:
                    raise ValueError(f"Duplicate contig {name}: {source}")
                self.index[name] = tuple(map(int, (length, offset, bases, width)))
        if not self.index:
            raise ValueError(f"No contigs: {source}")
        # samtools warns and ignores duplicate names; detect this rather than
        # allowing an ambiguous coordinate lookup to silently select a sequence.
        with open(self.path.parent / "faidx.stdout.log") as handle:
            if "duplicate sequence" in handle.read().lower():
                raise ValueError(f"Duplicate contig identifier in {source}")
        self.handle = open(self.path, "rb")
        first = next(iter(self.index))
        sampled = self.fetch(first, 0, min(self.length(first), 4096))
        if not set(sampled) <= set("ACGTRYMKSWBDHVN"):
            self.handle.close()
            raise ValueError(f"Genome FASTA contains non-nucleotide symbols in its first contig: {source}")

    def close(self):
        self.handle.close()

    def length(self, contig):
        if contig not in self.index:
            raise ValueError(f"Contig {contig!r} absent from FASTA (GFF3/FASTA mismatch)")
        return self.index[contig][0]

    def fetch(self, contig, start, end, strand="+"):
        from .sequence import reverse_complement
        length, offset, bases, width = self.index[contig]
        if not (0 <= start <= end <= length):
            raise ValueError(f"Invalid extraction {contig}:{start}-{end} (length {length})")
        if end == start:
            return ""
        begin_byte = offset + (start // bases) * width + start % bases
        end_byte = offset + ((end-1) // bases) * width + (end-1) % bases + 1
        self.handle.seek(begin_byte)
        seq = b"".join(self.handle.read(end_byte-begin_byte).split()).decode().upper()
        if len(seq) != end-start:
            raise ValueError(f"FASTA index/sequence mismatch: {contig}:{start}-{end}")
        return reverse_complement(seq) if strand == "-" else seq
