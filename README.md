# Gene Context Finder

Find candidate gene loci in multiple genomes using an annotation-derived panel of reference sequences.

The current [work plan](docs/WORK_PLAN.md) maps implemented stages to the remaining work for target-gene conservation and sequence selection across the study species. The active workflow is `stage1 → stage2 → stage3`, with two manual review stops. The [stage 3 guide](docs/STAGE3.md) explains the reviewed-sequence catalogue and conservation comparisons. The [NANOG verification](docs/WORKFLOW_FIXES_2026-11-09.md) records the actual run scope, fixes and remaining work.

The first review stop is **stage 1**: search all supplied genomes, assemble local matches into loci, automatically check available reference-neighbor contexts, and stop for manual selection. Multiple transcript variants and source species support one target/search group. The report ranks loci, maps matches for every query and supports selection by checkbox. Standard FASTA source headers are read automatically for every gene; a common source-dataset inventory can be reused across query panels. Missing source contexts remain explicit, and duplicate self comparisons are hidden by default.

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

## Reannotate selected windows: stage 2

Prepare a broad local protein reference once, then run the selected windows through Tiberius and BLASTP:

```bash
gene-context prepare-stage2-db --download-vertebrates --out /mnt/project/references/vertebrates_reviewed
gene-context setup-predictor --predictor-cache /mnt/project/cache/tiberius --models vertebrates mammalia
gene-context stage2 --results results/pou5_stage1 --selection results/pou5_selected \
  --genomes /mnt/project/genomes.tsv --protein-db /mnt/project/references/vertebrates_reviewed \
  --predictor-cache /mnt/project/cache/tiberius --out results/pou5_stage2 --threads 8
```

Open `results/pou5_stage2/report.html` for **STOP 2**: compare provided annotation and DNA-only predictions, inspect broad-database naming alternatives and target-query support, and choose target models/neighbors. No model is pre-approved. Overlapping selected windows share one bounded prediction interval. The downloadable draft/decisions preserve exact stage 2 IDs. Full details, draft-input rules and host2 memory limits are in [the stage 2 guide](docs/STAGE2.md). Analysis is offline; reference/model preparation is explicit.

After saving the review, download `stage2.decisions.json` and run:

```bash
gene-context stage3 --results results/pou5_stage2 \
  --decisions /path/to/stage2.decisions.json \
  --out results/pou5_comparison --threads 8
```

This reuses the selected sequences, runs protein/codon and local genomic comparisons, and writes a traceable sequence catalogue. It does not repeat search or prediction. Use `--validate-only` first to check the handoff without starting tools. See [stage 3](docs/STAGE3.md) for outputs and interpretation limits.

The report ends with two tables: at most one target representative per gene/genome, and separate cases for additional copies or uncertain assignments. Initial choices are explicit suggestions with source identity, CDS quality and reference identity/coverage evidence. You can change the representative, select none, and record the reason. Download split package exports the current tables and exact CDS/protein FASTAs. To apply exported choices to a separate full report without rerunning alignments:

```bash
gene-context stage3-finalize --results results/pou5_comparison \
  --decisions /path/to/stage3.final.decisions.json \
  --out results/pou5_final
```

This final assignment retains the original STOP-2 selection and all-by-all comparisons; it does not resolve orthology automatically.

## Earlier pipeline

The earlier `search → contexts → compare-reviewed` workflow is preserved and documented in [the legacy workflow](docs/LEGACY_WORKFLOW.md). Its model-first checkpoints differ from the new locus-first stage. Use the new `stage2` command for locus-first selections. Legacy `contexts` and `compare-reviewed` retain their different checkpoint format; they are not consumers of the new STOP-2 decisions.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests use local fixtures. Tool-dependent integration tests skip when executables are unavailable. The POU5 controls check annotated reference windows and neighbouring frog copies; they do not establish all-genome discovery accuracy.

Refresh an existing review into the compact layout without rerunning any search:

```bash
gene-context stage1-report --results results/nanog_stage1 --out results/nanog_compact
```

The original evidence and review/decision IDs are retained. The default view shows the first ten candidate loci in combined review order, each with a compact support table for every query variant. The score combines coverage- and identity-weighted BLAST support across source groups; contributions are available in Details and review_ranking.tsv. Synteny and technical details are collapsed; opening Synteny shows all available query-source neighborhoods. Curated GeneIDs can bind source neighborhoods across isoform differences while retaining sequence-discrepancy flags. See [metric definitions](docs/STAGE1_METRICS.md).
