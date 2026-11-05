from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .download import run_download
from .alignment import realign
from .pipeline import run
from .tools import tool_versions


def fraction(value):
    value = float(value)
    if not 0 < value <= 1:
        raise argparse.ArgumentTypeError("must be >0 and <=1")
    return value


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def nonnegative(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError("must be >=0")
    return value


def parser():
    parser = argparse.ArgumentParser(description="Find genes and their surroundings from DNA or protein sequences")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("doctor", help="Show external executable paths and versions")
    download = commands.add_parser("download", help="Download genome/GFF3 inputs and verify official MD5 checksums")
    download.add_argument("--manifest", required=True, help="TSV: genome_id, genome_url, optional gff_url and reference URLs")
    download.add_argument("--out", required=True, help="Managed download directory; produces genomes.tsv")
    download.add_argument("--only", nargs="+", help="Select genome IDs or exact assembly accessions")
    download.add_argument("--include-reference", action="store_true", help="Also fetch listed reference CDS/proteomes")
    download.add_argument("--dry-run", action="store_true", help="Print plan without network calls or filesystem changes")
    align = commands.add_parser("align", help="Re-align saved extraction without repeating genome searches")
    align.add_argument("--results", required=True, help="Existing run directory with manifest.json, models.tsv and sequences/")
    align.add_argument("--out", required=True, help="Separate new/empty output directory")
    align.add_argument("--threads", type=positive, default=4)
    align.add_argument("--dna-alignment", choices=["fast", "auto"], default="fast")
    align.add_argument("--include-partial", action="store_true")
    align.add_argument("--align-locus", action="store_true")
    align.add_argument("--max-alignment-bases", type=positive, default=2000000)
    run = commands.add_parser("run", help="Search every input genome, extract models/context and align sequences")
    run.add_argument("--query", required=True, help="FASTA with one or more query sequences")
    run.add_argument("--query-type", required=True, choices=["protein", "cds", "dna"])
    inputs = run.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--genomes", help="TSV: genome_id, fasta, optional gff3 (paths relative to TSV)")
    inputs.add_argument("--genome-dir", help="Recursively discover every .fa/.fna/.fasta[.gz]; no annotation assumed")
    run.add_argument("--out", required=True, help="New or empty output directory")
    run.add_argument("--threads", type=positive, default=4)
    run.add_argument("--flank", type=nonnegative, default=5000, help="Oriented upstream/downstream length (bp)")
    run.add_argument("--context", type=nonnegative, default=50000, help="Neighbour-gene search radius (bp)")
    run.add_argument("--max-intron", type=positive, default=200000)
    run.add_argument("--max-hits", type=positive, default=1000, help="Miniprot alignment cap per query/genome")
    run.add_argument("--secondary-ratio", type=fraction, default=0.1)
    run.add_argument("--min-coverage", type=fraction, default=0.5, help="Minimum aligned query fraction")
    run.add_argument("--min-identity", type=fraction, default=0.3, help="Protein identity fraction")
    run.add_argument("--protein-evalue", type=float, default=1e-5)
    run.add_argument("--dna-identity", type=fraction, default=0.7)
    run.add_argument("--min-dna-coverage", type=fraction, default=0.1,
                     help="Minimum union of nucleotide query HSP spans for annotated model association")
    run.add_argument("--dna-evalue", type=float, default=1e-5)
    run.add_argument("--word-size", type=positive, default=11, help="BLASTN word size (>=4)")
    run.add_argument("--min-dna-hsp", type=positive, default=30)
    run.add_argument("--min-cds-overlap", type=positive, default=30)
    run.add_argument("--dna-merge-distance", type=nonnegative, default=10000,
                     help="Merge nearby DNA HSP envelopes for region export only")
    run.add_argument("--include-partial", action="store_true",
                     help="Also align clean but incomplete CDS (never disrupted models)")
    run.add_argument("--skip-alignments", action="store_true")
    run.add_argument("--dna-alignment", choices=["fast", "auto"], default="fast",
                     help="fast: progressive FFT-NS-2 DNA alignment; auto: allow MAFFT iterative refinement")
    run.add_argument("--align-locus", action="store_true", help="Also align whole gene+flank regions")
    run.add_argument("--max-alignment-bases", type=positive, default=2000000)
    return parser


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command == "doctor":
        versions = tool_versions(["samtools", "miniprot", "blastn", "blastp", "makeblastdb", "mafft"])
        print(json.dumps(versions, indent=2))
        return 1 if any(not value["path"] for value in versions.values()) else 0
    if args.command == "download":
        try:
            return run_download(args)
        except (ValueError, OSError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "align":
        try:
            return realign(args)
        except (ValueError, OSError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.word_size < 4 or args.protein_evalue <= 0 or args.dna_evalue <= 0:
        parser().error("word size must be >=4 and E-values must be positive")
    try:
        return run(args)
    except (ValueError, OSError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
