# Validation of the generic sequence workflow

Validation below covers the `genecontext` package and its defined controls.
It is not an exhaustive genome-wide sensitivity/specificity study.

## Executables used locally

- Python 3.14.2 (standard library only).
- miniprot 0.18-r281, compiled from
  [81f9b93481cec16622eff339fc998e4f38da344e](https://github.com/lh3/miniprot/tree/81f9b93481cec16622eff339fc998e4f38da344e).
- NCBI BLAST+ 2.17.0+.
- samtools 1.24.
- MAFFT 7.526.

The workflow also targets Python >=3.10. The GitHub workflow passed on
3.10–3.13 for code commit `935de82`; external-tool integration tests skip in CI
when their required binaries are absent. All 28 tests, including actual-tool
integration tests, passed locally and on the Linux host.
A Linux-64 dry-run solve of `environment.yml` succeeded (81 packages). Actual conda environment creation was not required for the local or server runs, which used existing executables.

## Tests

`python -m unittest discover -s tests -v` runs 15 sequence/input unit tests, 7 downloader tests and 6 real-tool
integration tests when every required binary is on PATH. Missing executables
cause explicit integration-test skips, not mocked search results.

Downloader tests use a local HTTP server and cover official-MD5 checking, verified cache reuse, interrupted transfers, checksum failures, preservation of altered files, subset selection and dry-run behaviour. No large remote genome download is needed for those tests. The checksum parser preserves relative paths for nested NCBI entries that share a basename.

Coverage includes source-exon extraction across split codons, minus-strand
models, separate stop features, phase errors, internal stops, clipped flanks,
codon backtranslation, interval lookup, shared CDS parents, duplicate query IDs,
relative manifest paths and recursive FASTA discovery.

Real-tool integration covers two genomic copies on opposite strands including
a two-exon gene, multiple queries, a gzip-compressed negative genome, genomic
DNA queries with and without annotation, coding-DNA/protein routes, exact-model
deduplication, a missing input followed by a valid genome, exact-assembly GFF3
mismatch, and long/pipe-containing contig IDs. Search executables are not mocked.
It also exercises the saved-extraction alignment command and its progressive DNA
preset, including preservation of the original manifest.

## Biological example

The pinned miniprot repository supplies a mouse DPP3 protein and a human genomic
DPP3 test interval. `examples/run_miniprot_demo.sh` downloads those two small
FASTA files and runs the actual workflow. The observed output is one predicted
complete model, 17 CDS exons, 2,214 coding nucleotides including the stop, 737
amino acids, 100% query coverage and miniprot amino-acid identity 0.9282.

The source genomic sequence is an interval, not a whole chromosome. Its FASTA
identifier embeds a chromosome location, but all output coordinates refer to
the supplied FASTA sequence; the workflow does not parse header text as a
coordinate transform. This interval has no available external flanks at the
predicted coding boundaries, correctly recorded as clipped/empty windows.

These checks establish operation on defined controls. They do not establish
sensitivity across every gene family, assembly size, genetic code or evolutionary
distance. Whole-genome discovery remains bounded by the methods, input quality,
search thresholds, maximum intron length and miniprot reporting cap. Individual
copy orthology, expression and CRE activity require additional evidence.

## Linux server checks

On 2026-11-05, all 28 tests in this standalone repository also passed on a Linux compute host with Python 3.10.12, miniprot 0.18-r281, BLAST+ 2.17.0+, samtools 1.13 (htslib 1.21 development build), and MAFFT 7.526. This validates the existing-tool installation route as well as the source package; the conda file specifies newer minimum samtools for fresh installations. A separate five-assembly protein/CDS pilot uses public chicken NANOG and Pou5f3 reference queries. See the [full-genome pilot report](FULL_GENOME_PILOT.md) for measured outcomes, independent source-sequence checks, alignment recovery and limits of a single reference query.
