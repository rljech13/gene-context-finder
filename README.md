# Gene Context Finder

Find candidate gene loci in multiple genomes using an annotation-derived panel of reference sequences.

The current focus is **stage 1**: search all supplied genomes, assemble local matches into loci, automatically check available reference-neighbor contexts, and stop for manual selection. Multiple transcript variants and source species support one target/search group. The report ranks loci, maps matches for every query and supports selection by checkbox.

## Install

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
```

Stage 1 needs BLAST+ and samtools, already listed in `environment.yml`. Miniprot is optional. Genome annotation and annotated proteins/CDS enable the fast early-synteny pass; missing annotation is reported explicitly. This stage does not need Tiberius, model weights or Swiss-Prot.

## Try the POU5 panel

The included [test panel](examples/queries/pou5_test/README.md) contains 16 annotated transcript/CDS variants from human, zebrafish, chicken and Xenopus tropicalis. Source gene symbols and accession versions are retained. The POU5 group includes separately labelled family members and is not a one-to-one orthology assertion.

```bash
gene-context stage1 --gene POU5 \
  --query examples/queries/pou5_test/proteins.faa --query-type protein \
  --query-metadata examples/queries/pou5_test/queries.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv \
  --out results/pou5_stage1 --threads 4
```

This searches six small reference windows from four species, **not complete genomes**. To search the study assemblies, replace `--genomes` with their explicit TSV (`genome_id`, `fasta`; relative or absolute paths, gzip accepted).

With annotated datasets the default search is BLASTP over all available protein isoforms, grouped by their annotated gene loci. `--discovery genome` requests whole-genome TBLASTN for rescue. DNA queries use BLASTN.

Open `results/pou5_stage1/report.html`. Inspect the ranked loci, query tracks and available early-synteny evidence. Select loci, choose a short reason and save each pair decision. Validate and download the decision TSV. The next stage does not start automatically. The small control-window example does not include full annotation/proteome datasets, so its early context is explicitly unavailable; for study assemblies, include `gff3` and `protein` or `cds` in the genome TSV. Additional reference datasets can be supplied with `--source-datasets`.

```bash
gene-context stage1-select --results results/pou5_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv \
  --out results/pou5_selected
```

The output contains selected loci and windows in TSV/BED, a validation summary and the decision snapshot. It does not start gene prediction. The [stage 1 guide](docs/STAGE1.md) explains automatic context preparation, source-query resolution, caching, the review screen and reuse of a saved search through `stage1-context`.

On host2, use the [constrained launcher](docs/LIMITED_EXPERIMENTS.md) for every experiment. Its shared hard limit is 250 GiB RAM with swap disabled; use absolute `/mnt/data` paths. The installed source snapshot must contain v0.6.0 or later.

## Existing full pipeline

The earlier `search → contexts → compare-reviewed` workflow is preserved and documented in [the legacy workflow](docs/LEGACY_WORKFLOW.md). Its model-first checkpoints differ from the new locus-first stage. The new selected TSV/BED is not yet wired into legacy `contexts`; this interface is the next step, after reviewing stage 1 behaviour.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests use local fixtures. Tool-dependent integration tests skip when executables are unavailable. The POU5 controls check annotated reference windows and neighbouring frog copies; they do not establish all-genome discovery accuracy.

Refresh an existing review into the compact layout without rerunning any search:

```bash
gene-context stage1-report --results results/nanog_stage1 --out results/nanog_compact
```

The original evidence and review/decision IDs are retained. The default view shows one candidate locus, a compact support table for every query variant and a small diagram of the candidate and all available query-source neighborhoods; detailed evidence and controls expand on demand. See [metric definitions](docs/STAGE1_METRICS.md).
