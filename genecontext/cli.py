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


def neighbor_limit(value):
    value = positive(value)
    if value > 10:
        raise argparse.ArgumentTypeError("must be between 1 and 10 genes per side")
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
    command.add_argument("--tiberius-worker-mode", choices=["isolated", "reuse"], default="isolated",
                         help="Reuse neural weights across independent windows, or start one process per window")
    command.add_argument("--tiberius-batch-size", type=positive, default=1,
                         help="Neural inference batch size; 1 limits CPU memory")


def early_synteny_options(command, disable_option=True):
    command.add_argument("--query-metadata-lookup", choices=["auto", "off"], default="auto", help="Resolve missing source GeneIDs from exact versioned RefSeq records; cached, sequence-verified, no query replacement")
    if disable_option:
        command.add_argument("--early-synteny", choices=["auto", "off"], default="auto",
                             help="Automatically check annotated neighbor contexts before STOP 1")
    else:
        command.set_defaults(early_synteny="auto", genomes=None)
    command.add_argument("--datasets", help="Optional complete datasets TSV; default: the search genome list")
    command.add_argument("--source-datasets", help="Optional source-reference datasets TSV; query-to-gene mapping is automatic")
    command.add_argument("--neighbor-cache", help="Reusable context/proteome cache; default: _genecontext_cache beside the output")
    command.add_argument("--synteny-neighbors", type=neighbor_limit, default=5)
    command.add_argument("--synteny-max-neighbors", type=neighbor_limit, default=10,
                         help="Maximum neighboring genes per side (1–10; hard upper bound 10)")
    command.add_argument("--synteny-max-distance", type=positive, default=1000000)


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
    download.add_argument("--annotation-only", action="store_true", help="Fetch source GFF/proteins/CDS without the full genome; suitable for preliminary source contexts")
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
    check = commands.add_parser("check-targets", help="Check saved candidate proteins against a local reference panel; no downloads or decisions")
    check.add_argument("--proteins", required=True, help="Candidate protein FASTA, e.g. model_options.proteins.faa")
    check.add_argument("--references", required=True, help="Local reference protein FASTA, including close paralogs")
    check.add_argument("--reference-labels", help="Optional complete TSV: reference, family; exact reference FASTA IDs")
    check.add_argument("--out", required=True, help="Separate new/empty directory")
    check.add_argument("--threads", type=positive, default=4)
    loci = commands.add_parser("stage1", help="Find annotated protein loci first, check contexts, then STOP for review")
    loci.add_argument("--gene", required=True, help="One target/search-group label, e.g. POU5")
    loci.add_argument("--query", required=True)
    loci.add_argument("--query-type", required=True, choices=["protein", "cds", "dna"])
    loci.add_argument("--query-metadata", help="Complete TSV: query, species, optional group/gene_symbol/transcript")
    loci.add_argument("--genomes", required=True, help="TSV: genome_id, fasta; optional gff3, protein, cds (paths relative to TSV)")
    loci.add_argument("--discovery", choices=["auto", "blastp", "genome"], default="auto",
                      help="auto: BLASTP all annotated products when GFF exists, otherwise TBLASTN; blastp: require annotation; genome: explicit TBLASTN rescue. DNA always uses BLASTN")
    loci.set_defaults(genome_dir=None)
    loci.add_argument("--out", required=True)
    loci.add_argument("--threads", type=positive, default=4)
    loci.add_argument("--jobs", type=positive, default=1, help="Concurrent genomes; --threads is per genome, so CPU demand is jobs × threads")
    loci.add_argument("--search-evalue", type=float, default=1e-3, help="BLAST retrieval ceiling; weaker matches cannot be recovered by refiltering")
    loci.add_argument("--seed-min-length", type=positive, default=20, help="Minimum HSP alignment length: aa for protein/CDS, nt for DNA")
    loci.add_argument("--seed-min-identity", type=fraction, default=0.20)
    loci.add_argument("--max-intron", type=positive, default=200000)
    loci.add_argument("--chain-max-overlap", type=fraction, default=0.15)
    loci.add_argument("--max-locus-span", type=positive, default=2000000)
    loci.add_argument("--locus-shared-fraction", type=fraction, default=0.5, help="Shared aligned genomic bases / smaller chain; not envelope overlap")
    loci.add_argument("--refine", choices=["auto", "miniprot", "none"], default="none", help="Optional local boundary refinement for genomic searches; annotated BLASTP loci use GFF bounds")
    loci.add_argument("--refine-padding", type=nonnegative, default=50000)
    loci.add_argument("--refine-max-models", type=positive, default=10, help="Miniprot alternative-model cap per query/local window; does not cap discovered loci")
    loci.add_argument("--window", type=nonnegative, default=50000, help="Proposed synteny flank length; no context prediction runs here")
    early_synteny_options(loci)
    annotated = commands.add_parser("stage1-context", help="Add automatic early synteny to a saved locus search without repeating BLAST")
    annotated.add_argument("--results", required=True)
    annotated.add_argument("--out", required=True)
    annotated.add_argument("--query-source-fasta", help="Restore explicit source metadata from the unchanged original FASTA in an older checkpoint")
    annotated.add_argument("--threads", type=positive, default=4)
    early_synteny_options(annotated, disable_option=False)
    presentation = commands.add_parser("stage1-report", help="Refresh a saved report without searching, changing evidence or invalidating decision IDs")
    presentation.add_argument("--results", required=True)
    presentation.add_argument("--out", required=True, help="Separate new/empty report directory")
    filtered = commands.add_parser("stage1-filter", help="Change locus review thresholds without another genome search")
    filtered.add_argument("--results", required=True)
    filtered.add_argument("--out", required=True)
    selected = commands.add_parser("stage1-select", help="Validate explicit locus decisions and export chosen synteny windows")
    selected.add_argument("--results", required=True)
    selected.add_argument("--decisions", required=True)
    selected.add_argument("--out")
    selected.add_argument("--validate-only", action="store_true", help="List all decision errors and proposed windows without creating output")
    for command, defaults in [(loci, True), (filtered, False)]:
        command.add_argument("--min-query-coverage", type=fraction, default=0.20 if defaults else None)
        command.add_argument("--min-identity", type=fraction, default=0.25 if defaults else None)
        command.add_argument("--max-evalue", type=float, default=1e-5 if defaults else None)
    db = commands.add_parser("prepare-stage2-db", help="Prepare a pinned local BLASTP catalogue for window protein identification")
    source = db.add_mutually_exclusive_group(required=True)
    source.add_argument("--fasta", help="Local protein FASTA, plain or gzip; UniProt GN/OX metadata retained")
    source.add_argument("--download-vertebrates", action="store_true", help="Explicitly download reviewed UniProt vertebrate proteins and isoforms")
    db.add_argument("--out", required=True)
    s2 = commands.add_parser("stage2", help="Reannotate approved genomic windows and STOP for model/neighbor review")
    s2p = commands.add_parser("stage2-plan", help="Resolve preliminary-neighborhood window bounds without extraction or prediction")
    for command in [s2,s2p]:
        command.add_argument("--results", required=True, help="Sealed locus-first stage 1 checkpoint")
        choice = command.add_mutually_exclusive_group(required=True)
        choice.add_argument("--decisions", help="Completed stage 1 decisions TSV")
        choice.add_argument("--selection", help="Approved output directory from stage1-select")
        choice.add_argument("--draft", help="Portable stage 1 draft; unsaved pairs are rejected unless --committed-only")
        command.add_argument("--committed-only", action="store_true", help="Explicitly process saved draft pairs only; never process unsaved checkboxes")
        command.add_argument("--only-genomes", nargs="+", help="Explicit subset of approved genome IDs, for a pilot")
        command.add_argument("--window-policy", choices=["pre-synteny","fixed"], default="pre-synteny", help="Default: full preliminary gene envelope; fixed retains configured flanks. Explicit manual windows take precedence.")
        command.add_argument("--max-window-bases", type=positive, default=2000000)
        command.add_argument("--out", required=True)
    s2.add_argument("--resume", action="store_true", help="Continue an interrupted stage2 output with identical inputs and parameters")
    s2.add_argument("--genomes", required=True, help="Dataset TSV with local genomic FASTA and optional matching GFF3")
    s2.add_argument("--protein-db", required=True, help="Verified local prepare-stage2-db output; no downloads during analysis")
    s2.add_argument("--threads", type=positive, default=4)
    s2.add_argument("--identification-coverage", type=fraction, default=.5)
    s2.add_argument("--identification-identity", type=fraction, default=.3)
    s2.add_argument("--identification-evalue", type=float, default=1e-5)
    s2.add_argument("--competing-ratio", type=fraction, default=.95)
    predictor_options(s2)
    s2.set_defaults(predictor_cache=None)
    next(a for a in s2._actions if a.dest=="predictor").choices=["tiberius"]
    compact = commands.add_parser("compact-report", help="Create a smaller standalone HTML beside an existing report; retain the original")
    compact.add_argument("--report", required=True, help="Existing report.html")
    compact.add_argument("--out", required=True, help="New HTML filename in the same directory")
    s2r = commands.add_parser("stage2-report", help="Refresh STOP 2 presentation without prediction or protein searches")
    s2r.add_argument("--results", required=True)
    s2r.add_argument("--out", required=True)
    s2r.add_argument("--genomes", help="Add all gene-level annotation features from the original, hash-matched GFF3 files; no prediction/search")
    s2q = commands.add_parser("stage2-quality", help="Collect hash-matched original annotation flags for an existing STOP-2 checkpoint; no prediction")
    s2q.add_argument("--results", required=True)
    s2q.add_argument("--genomes", required=True)
    s2q.add_argument("--out", required=True)
    s2t = commands.add_parser("stage2-reconstruct", help="Reconstruct existing annotation CDS from exact splice-aware reference-transcript mapping")
    s2t.add_argument("--results", required=True)
    s2t.add_argument("--transcripts", required=True, help="Offline GenBank XML transcript records, including CDS coordinates")
    s2t.add_argument("--mapping", required=True, help="TSV with model_id and accession columns")
    s2t.add_argument("--out", required=True)
    s2t.add_argument("--minimap2", default="minimap2")
    s2t.add_argument("--threads", type=positive, default=4)
    s2t.add_argument("--max-intron", type=positive, default=1000000)
    s2t.add_argument("--attach", action="store_true", help="Attach verified supplement and FASTA/GFF exports; keep models, review ID and decisions unchanged")
    s2i = commands.add_parser("stage2-reidentify", help="Repeat broad protein identification from a sealed STOP 2 result, keeping predictions unchanged")
    s2i.add_argument("--results", required=True)
    s2i.add_argument("--protein-db", required=True)
    s2i.add_argument("--out", required=True)
    s2i.add_argument("--threads", type=positive, default=4)
    s3 = commands.add_parser("stage3", help="Compare models explicitly chosen in the new STOP-2 review; export selected sequences")
    s3.add_argument("--results", required=True, help="Sealed locus-first stage2 result")
    s3.add_argument("--decisions", required=True, help="Exported stage2.decisions.json (or a fully saved portable draft)")
    s3.add_argument("--out", required=True)
    s3.add_argument("--validate-only", action="store_true", help="Check decisions and sequence integrity without creating output or running aligners")
    s3.add_argument("--threads", type=positive, default=4)
    s3.add_argument("--include-partial", action="store_true")
    s3.add_argument("--dna-engine", choices=["blastn", "lastz"], default="blastn")
    s3.add_argument("--dna-min-block", type=positive, default=50)
    s3.add_argument("--dna-min-identity", type=fraction, default=.6)
    s3.add_argument("--max-alignment-bases", type=positive, default=2000000)
    s3.add_argument("--control-run", action="store_true", help="Mark all outputs as a software control, not a biological approval")
    s3report=commands.add_parser("stage3-report",help="Refresh a sealed stage-3 report without repeating alignments")
    s3report.add_argument("--results",required=True)
    s3report.add_argument("--out",required=True)
    final=commands.add_parser("stage3-finalize",help="Apply exported primary-target/special-case choices without repeating alignments")
    final.add_argument("--results",required=True)
    final.add_argument("--decisions",required=True,help="Exported stage3.final.decisions.json")
    final.add_argument("--out",required=True)
    return parser


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command == 'compact-report':
        from .compact_report import run as compact
        try:
            return compact(args)
        except (OSError, ValueError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command in {"stage3","stage3-report","stage3-finalize"}:
        from .stage3 import run as stage3
        from .stage3_report import refresh
        try:
            return (stage3 if args.command=="stage3" else refresh)(args)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command in {'stage2-quality','stage2-reconstruct'}:
        if args.command=='stage2-quality':
            from .stage2_quality import run
        else:
            from .stage2_reconstruct import run
        try:
            return run(args)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command in {"stage2", "stage2-plan", "stage2-report", "stage2-reidentify", "prepare-stage2-db"}:
        from .stage2 import run as stage2, reidentify, plan
        from .stage2_db import prepare
        from .stage2_report import refresh
        try:
            return {"stage2":stage2,"stage2-plan":plan,"stage2-report":refresh,"stage2-reidentify":reidentify,"prepare-stage2-db":prepare}[args.command](args)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command in {"stage1", "stage1-filter", "stage1-select", "stage1-context", "stage1-report"}:
        from .stage1 import search, refilter, approve
        from .stage1_context import context
        from .stage1_report import refresh
        if args.command == "stage1" and args.jobs > 1:
            from .stage1_parallel import search
        try:
            return {"stage1": search, "stage1-filter": refilter, "stage1-select": approve, "stage1-context": context, "stage1-report": refresh}[args.command](args)
        except (OSError, ValueError, RuntimeError, KeyError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    if args.command == "check-targets":
        from .target_check import check_targets
        try:
            return check_targets(args)
        except (OSError, ValueError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
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
