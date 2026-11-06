from __future__ import annotations

import argparse
import json
import subprocess
import sys

from . import __version__
from .download import run_download
from .alignment import realign
from .pipeline import run
from .tools import tool_versions
from .context_reference import DEFAULT_CACHE
from .tiberius_backend import DEFAULT_PREDICTOR_CACHE, MODELS


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


def comparison_options(command):
    command.add_argument("--context-mode", choices=["predict", "annotation"], default="predict",
                         help="predict: predict genes in every window, then identify proteins; annotation: legacy GFF3/reference-supported contexts")
    command.add_argument("--augustus-species", help="AUGUSTUS trained species model; genome TSV augustus_species values override this default")
    predictor_options(command)
    command.add_argument("--offline", action="store_true", help="Require local model/reference caches; do not download missing files")
    command.add_argument("--include-partial", action="store_true",
                         help="Also compare clean partial coding models; disrupted models stay excluded")
    command.add_argument("--window", type=nonnegative, default=50000,
                         help="Context DNA on each side of the gene/model (bp)")
    command.add_argument("--dna-engine", choices=["blastn", "lastz"], default="blastn")
    command.add_argument("--context-proteins", help="Optional naming reference proteins for predicted context genes")
    command.add_argument("--context-reference", choices=["auto", "swissprot", "none"], default="swissprot",
                         help="auto: use Swiss-Prot if no local context proteins exist; swissprot: always add it; none: only local naming references")
    command.add_argument("--reference-cache", default=DEFAULT_CACHE,
                         help="Reusable Swiss-Prot cache; an existing release is pinned and reused offline")
    command.add_argument("--query-groups", help="Optional TSV: query, group; combine references of one target")
    command.add_argument("--select-loci", help="Optional TSV: group, genome_id, model_id; override ranked candidates")
    command.add_argument("--reference-genome", help="Genome ID used for the DNA plot; tables still compare all pairs")
    command.add_argument("--neighbor-coverage", type=fraction, default=0.5)
    command.add_argument("--neighbor-identity", type=fraction, default=0.3)
    command.add_argument("--dna-min-block", type=positive, default=50)
    command.add_argument("--dna-min-identity", type=fraction, default=0.6)
    command.add_argument("--max-window-bases", type=positive, default=2000000)


def predictor_options(command):
    command.add_argument("--predictor", choices=["tiberius", "augustus"], default="tiberius",
                         help="DNA-only window predictor (default: Tiberius on CPU)")
    command.add_argument("--tiberius-model", default="vertebrates",
                         help="vertebrates, mammalia, or an unmasked local model directory; TSV tiberius_model overrides this")
    command.add_argument("--prediction-models", help="Optional TSV: genome_id, model; overrides models for new or saved searches")
    command.add_argument("--predictor-cache", default=DEFAULT_PREDICTOR_CACHE)
    command.add_argument("--tiberius-python", help="Optional Python executable of an existing Tiberius 2.0.9 CPU environment")
    command.add_argument("--tiberius-batch-size", type=positive, default=1,
                         help="Neural inference batch size; 1 limits CPU memory")


def search_options(command, comparative=False):
    command.add_argument("--query", required=True, help="FASTA with one or more query sequences")
    command.add_argument("--query-type", required=True, choices=["protein", "cds", "dna"])
    inputs = command.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--genomes", help="TSV: genome_id, fasta, optional gff3 (paths relative to TSV)")
    inputs.add_argument("--genome-dir", help="Recursively discover every .fa/.fna/.fasta[.gz]; no annotation assumed")
    command.add_argument("--out", required=True, help="New or empty output directory")
    command.add_argument("--threads", type=positive, default=4)
    command.add_argument("--max-intron", type=positive, default=200000)
    command.add_argument("--max-hits", type=positive, default=1000, help="Miniprot alignment cap per query/genome")
    command.add_argument("--secondary-ratio", type=fraction, default=0.1)
    command.add_argument("--min-coverage", type=fraction, default=0.5, help="Fast-route minimum aligned query fraction; local discovery has its own lower gate")
    command.add_argument("--min-identity", type=fraction, default=0.3, help="Protein identity fraction")
    command.add_argument("--protein-search", choices=["tblastn", "miniprot"], default="tblastn",
                         help="tblastn: sensitive genome discovery then local coding reconstruction; miniprot: fast whole-genome mapping")
    command.add_argument("--candidate-min-coverage", dest="rescue_min_coverage", type=fraction, default=0.15,
                         help="Minimum BLASTP query fraction for a coding candidate; not a per-exon seed cutoff")
    command.add_argument("--seed-min-aa", dest="rescue_min_aa", type=positive, default=30)
    command.add_argument("--refine-flank", dest="rescue_flank", type=positive, default=50000,
                         help="DNA on either side of a translated seed for dense local coding reconstruction")
    command.add_argument("--max-refine-regions", dest="rescue_max_regions", type=positive, default=20,
                         help="Maximum unannotated seed regions refined per query/genome; all passing seeds are retained")
    command.add_argument("--protein-evalue", type=float, default=1e-5)
    command.add_argument("--dna-identity", type=fraction, default=0.7)
    command.add_argument("--min-dna-coverage", type=fraction, default=0.1,
                     help="Minimum union of nucleotide query HSP spans for annotated model association")
    command.add_argument("--dna-evalue", type=float, default=1e-5)
    command.add_argument("--word-size", type=positive, default=11, help="BLASTN word size (>=4)")
    command.add_argument("--min-dna-hsp", type=positive, default=30)
    command.add_argument("--min-cds-overlap", type=positive, default=30)
    command.add_argument("--dna-merge-distance", type=nonnegative, default=10000,
                     help="Merge nearby DNA HSP envelopes for region export only")
    command.add_argument("--max-alignment-bases", type=positive, default=2000000)
    if comparative:
        command.set_defaults(flank=0, context=0, skip_alignments=True, include_partial=False,
                             align_locus=False, dna_alignment="fast")
    else:
        command.add_argument("--flank", type=nonnegative, default=5000, help="Oriented flank length (bp)")
        command.add_argument("--context", type=nonnegative, default=50000, help="Neighbour-gene search radius (bp)")
        command.add_argument("--include-partial", action="store_true",
                         help="Also align clean but incomplete CDS (never disrupted models)")
        command.add_argument("--skip-alignments", action="store_true")
        command.add_argument("--dna-alignment", choices=["fast", "auto"], default="fast",
                         help="fast: progressive FFT-NS-2 DNA alignment; auto: allow MAFFT iterative refinement")
        command.add_argument("--align-locus", action="store_true", help="Also align whole gene+flank regions")


def parser():
    parser = argparse.ArgumentParser(description="Find genes and their surroundings from DNA or protein sequences")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor", help="Check sequence tools and the selected predictor")
    predictor_options(doctor)
    setup = commands.add_parser("setup-predictor", help="Install the pinned Tiberius CPU runtime and prepare shared model weights")
    setup.add_argument("--predictor-cache", default=DEFAULT_PREDICTOR_CACHE)
    setup.add_argument("--python", help="Python 3.12/3.13 used to create the isolated runtime")
    setup.add_argument("--models", nargs="+", choices=sorted(MODELS), default=["vertebrates", "mammalia"])
    setup.add_argument("--models-only", action="store_true", help="Prepare weights without installing Python packages")
    setup.add_argument("--with-reference", action="store_true", help="Also prepare Swiss-Prot for subsequent offline identification")
    setup.add_argument("--reference-cache", default=DEFAULT_CACHE)
    reference = commands.add_parser("prepare-context", help="Download and checksum-pin Swiss-Prot for offline context annotation")
    reference.add_argument("--reference-cache", default=DEFAULT_CACHE)
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
    search_options(run)
    analyze = commands.add_parser("analyze", help="Automatic diagnostic search/comparison without manual stops")
    search_options(analyze, comparative=True)
    comparison_options(analyze)
    compare = commands.add_parser("compare", help="Automatic diagnostic comparison of saved search loci")
    compare.add_argument("--results", required=True)
    compare.add_argument("--out", required=True)
    compare.add_argument("--threads", type=positive, default=4)
    compare.add_argument("--max-alignment-bases", type=positive, default=2000000)
    comparison_options(compare)
    search = commands.add_parser("search", help="Stage 1: discover candidates and STOP for manual locus review")
    search_options(search, comparative=True)
    search.add_argument("--query-groups", help="Optional TSV grouping reference queries of one target")
    review = commands.add_parser("review-search", help="Create STOP 1 from a saved search without searching genomes again")
    review.add_argument("--results", required=True)
    review.add_argument("--out", required=True)
    review.add_argument("--query-groups")
    contexts = commands.add_parser("contexts", help="Stage 2: inspect explicitly selected loci, predict genes, and STOP for CDS/neighbour review")
    contexts.add_argument("--review", required=True, help="Stage 1 review directory")
    contexts.add_argument("--search-results", help="Relocated source search; fingerprints must match the checkpoint")
    contexts.add_argument("--decisions", required=True, help="Completed loci.decisions.tsv")
    contexts.add_argument("--out", required=True)
    contexts.add_argument("--threads", type=positive, default=4)
    contexts.add_argument("--max-alignment-bases", type=positive, default=2000000)
    comparison_options(contexts)
    imported = commands.add_parser("review-contexts", help="Create STOP 2 from saved comparisons without repeating prediction")
    imported.add_argument("--results", required=True, help="Existing comparison directory")
    imported.add_argument("--search-results", help="Relocated source search; fingerprints must match")
    imported.add_argument("--out", required=True)
    reviewed = commands.add_parser("compare-reviewed", help="Stage 3: compare approved CDS/neighbours, reusing predictions")
    reviewed.add_argument("--review", required=True, help="Stage 2 review directory")
    reviewed.add_argument("--search-results", help="Relocated source search; preserves the checkpoint review ID")
    reviewed.add_argument("--models", required=True, help="Completed models.decisions.tsv")
    reviewed.add_argument("--neighbors", required=True, help="Completed neighbors.decisions.tsv")
    reviewed.add_argument("--out", required=True)
    reviewed.add_argument("--threads", type=positive, default=4)
    reviewed.add_argument("--dna-engine", choices=["blastn", "lastz"])
    reviewed.add_argument("--dna-min-block", type=positive)
    reviewed.add_argument("--dna-min-identity", type=fraction)
    reviewed.add_argument("--reference-genome")
    reviewed.add_argument("--max-alignment-bases", type=positive)
    return parser


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command in {"search", "review-search", "contexts", "compare-reviewed", "review-contexts"}:
        from .review import search, review_search, contexts, compare_reviewed, review_contexts
        function = {"search": search, "review-search": review_search, "contexts": contexts, "compare-reviewed": compare_reviewed, "review-contexts": review_contexts}[args.command]
        try:
            return function(args)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "doctor":
        versions = tool_versions(["samtools", "miniprot", "blastn", "blastp", "blastx", "tblastn", "makeblastdb", "mafft"] + (["augustus"] if args.predictor == "augustus" else []))
        if args.predictor == "tiberius":
            from .tiberius_backend import runtime_info, runtime_python, ensure_model
            python = runtime_python(args)
            try:
                versions["tiberius"] = dict(path=python, **runtime_info(python), model=ensure_model(args.tiberius_model, args.predictor_cache, offline=True))
            except (ValueError, OSError, RuntimeError) as error:
                versions["tiberius"] = dict(path=None, error=str(error))
        print(json.dumps(versions, indent=2))
        return 1 if any(not value["path"] for value in versions.values()) else 0
    if args.command == "setup-predictor":
        from .tiberius_backend import setup_predictor
        try:
            return setup_predictor(args)
        except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "download":
        try:
            return run_download(args)
        except (ValueError, OSError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "prepare-context":
        from .context_reference import swissprot
        try:
            _, receipt = swissprot(args.reference_cache)
            print(json.dumps(receipt, indent=2))
            return 0
        except (ValueError, OSError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "align":
        try:
            return realign(args)
        except (ValueError, OSError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command in {"compare", "analyze"}:
        from .comparison import compare, analyze
        try:
            return (compare if args.command == "compare" else analyze)(args)
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
