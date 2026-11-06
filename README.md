# Gene Context Finder

Find candidate gene loci in multiple genomes using an annotation-derived panel of reference sequences.

The current focus is **stage 1**: search all supplied genomes, assemble local matches into loci, and stop for manual selection of the regions to examine for synteny. Multiple transcript variants and source species support one target/search group.

## Install

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
```

Stage 1 needs BLAST+, samtools and miniprot, already listed in `environment.yml`. It does not need Tiberius, model weights, Swiss-Prot or a genome annotation.

## Try the POU5 panel

The included [test panel](examples/queries/pou5_test/README.md) contains 16 annotated transcript/CDS variants from human, zebrafish, chicken and Xenopus tropicalis. Source gene symbols and accession versions are retained. The POU5 group includes separately labelled family members and is not a one-to-one orthology assertion.

```bash
gene-context stage1 --gene POU5 \
  --query examples/queries/pou5_test/proteins.faa --query-type protein \
  --query-metadata examples/queries/pou5_test/queries.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv \
  --refine miniprot --out results/pou5_stage1 --threads 4
```

This searches six small reference windows from four species, **not complete genomes**. To search the study assemblies, replace `--genomes` with their explicit TSV (`genome_id`, `fasta`; relative or absolute paths, gzip accepted).

Open `results/pou5_stage1/report.html`. Inspect the query tracks, change display filters, select loci and correct anchor/window coordinates if needed. Download the decision TSV. Nothing advances automatically.

```bash
gene-context stage1-select --results results/pou5_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv \
  --out results/pou5_selected
```

The output contains selected loci and synteny windows in TSV/BED. It does not start gene prediction. The [stage 1 guide](docs/STAGE1.md) explains assembly rules, persistent filtering without another BLAST run, manual splits and limitations.

On host2, use the [constrained launcher](docs/LIMITED_EXPERIMENTS.md) for every experiment. Its shared hard limit is 250 GiB RAM with swap disabled; use absolute `/mnt/data` paths. The installed source snapshot must contain v0.6.0 or later.

## Existing full pipeline

The earlier `search → contexts → compare-reviewed` workflow is preserved and documented in [the legacy workflow](docs/LEGACY_WORKFLOW.md). Its model-first checkpoints differ from the new locus-first stage. The new selected TSV/BED is not yet wired into legacy `contexts`; this interface is the next step, after reviewing stage 1 behaviour.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests use local fixtures. Tool-dependent integration tests skip when executables are unavailable. The POU5 controls check annotated reference windows and neighbouring frog copies; they do not establish all-genome discovery accuracy.
