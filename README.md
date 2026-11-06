# Gene Context Finder

Find candidate genes and their surroundings in a set of genomes using DNA or protein sequences.

```text
Query FASTA + genome FASTAs + optional matching GFF3
                       ↓
              1. Search all genomes
              STOP: choose loci to inspect
                       ↓
       2. Extract ±50 kb, predict and identify genes
              STOP: choose CDS models and neighbours
                       ↓
       3. Compare the reviewed DNA, CDS and proteins
```

The recommended workflow has **two manual review stops**. Each stage writes a report and exits. Blank decisions never approve a candidate; changing a decision file does not automatically restart a process. See the [staged workflow](docs/STAGED_WORKFLOW.md).

## What to prepare

- **Query FASTA:** one or more proteins, CDS sequences or DNA regions. Use one sequence type per file. Gene names alone are not sequence input.
- **Genome FASTA:** one file per assembly; `.fa`, `.fna`, `.fasta` and gzip files are accepted.
- **GFF3 annotation, if available:** it must match the exact genome assembly. Gzip is accepted.
- **Prediction model:** the shared Tiberius `vertebrates` model is the default; select `mammalia` for mammals. These are two shared models, not one download per species. For other clades, provide a suitable unmasked Tiberius model or use the optional AUGUSTUS backend.

Use already downloaded files directly, with their original filenames.

## Install once

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
gene-context setup-predictor --with-reference
gene-context doctor
```

Setup creates an isolated CPU runtime and downloads both shared models (66.2 MB compressed in total), plus Swiss-Prot for protein identification. No GPU is required. Allow about 2 GB on macOS or 3 GB on Linux for the predictor installation, in addition to genomes and results. See [setup](docs/SETUP_AND_DOWNLOAD.md) for cache locations and offline use.

For experiments on host2, use the [constrained launcher and shared RAM policy](docs/LIMITED_EXPERIMENTS.md). Thread counts alone do not enforce a memory budget.

## 1. List your local genomes

Create `genomes.tsv` with **tab-separated** columns. Replace the example paths with your actual filenames:

```tsv
genome_id	fasta	gff3	tiberius_model
chicken	data/chicken/genome.fna.gz	data/chicken/annotation.gff.gz	vertebrates
human	data/human/genome.fna.gz		mammalia
```

Add one row per assembly. Paths are relative to `genomes.tsv`, or absolute. Leave `gff3` empty when annotation is unavailable. Gzip inputs need no manual decompression.

`tiberius_model` is optional: omitted entries use `--tiberius-model vertebrates`. Select `mammalia` explicitly for mammals. A separate `--prediction-models models.tsv` table can override models, including for saved searches; the [18-species example](examples/vertebrates18.models.tsv) matches the download manifest. No taxonomy is inferred from filenames.

## 2. Search, review, then compare

A chicken Pou5f3 [protein query](examples/queries/chicken_Pou5f3.protein.faa) is included, with [CDS and provenance](examples/queries/README.md). Put your query proteins in `targets.faa`.

```bash
gene-context search --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/01_search --threads 8
```

**Stop 1:** open `results/01_search/review/report.html`. Fill `loci.decisions.tsv`: inspect one or several candidate loci, exclude a pair, or leave it explicitly unresolved. Rankings remain suggestions.

```bash
gene-context contexts --review results/01_search/review \
  --decisions results/01_search/review/loci.decisions.tsv \
  --out results/02_contexts --threads 8
```

This extracts the selected windows and predicts their genes with Tiberius on CPU, including genomes with GFF3. **Stop 2:** open `results/02_contexts/report.html`. Choose the actual target CDS in `models.decisions.tsv` and review the neighbours in `neighbors.decisions.tsv`. Source isoforms, prediction models, exon tracks and sequence exports are retained separately.

For additional evidence before choosing, `check-targets` compares the saved proteins against a small local reference panel including close paralogs. It reports coverage of both proteins and competing matches without downloading data or changing decisions. See the [reference check](docs/TARGET_REFERENCE_CHECK.md).

```bash
gene-context compare-reviewed --review results/02_contexts \
  --models results/02_contexts/models.decisions.tsv \
  --neighbors results/02_contexts/neighbors.decisions.tsv \
  --out results/03_comparison --threads 8
```

The final command reuses the saved predictions. It compares only explicitly kept targets/neighbours and writes `results/03_comparison/report.html`: context maps, synteny, local genomic DNA blocks, protein alignments and codon-guided CDS alignments. Changing the chosen model changes the actual sequences used downstream.

Use `--query-type cds` for spliced coding DNA, or `--query-type dna` for genomic/noncoding regions. Use `--window` on `contexts` to change the default 50 kb flanks. After setup, `--offline` on `contexts` requires cached predictor/reference files. Paths and caches are explained in the [stage-by-stage instructions](docs/STAGED_WORKFLOW.md).

## 3. Outputs to inspect

| Output | Purpose |
|---|---|
| `01_search/search/` | All search candidates, extracted sequences and original evidence |
| `01_search/review/report.html`, `loci.decisions.tsv` | First review: which loci advance |
| `02_contexts/report.html`, `model_options.tsv` | CDS options, exon tracks and provisional surroundings |
| `02_contexts/models.decisions.tsv`, `neighbors.decisions.tsv` | Second review: selected CDS and neighbours |
| `03_comparison/report.html` | Comparative reports for reviewed targets |
| `03_comparison/reviewed_targets.tsv`, `reviewed_targets.gff3` | Exact chosen models and actual flank lengths |
| `03_comparison/reviewed_target_cds.fna`, `reviewed_target_proteins.faa` | Sequences actually used for comparison |
| `03_comparison/approved_models.tsv`, `approved_neighbors.tsv` | Decision snapshots |
| Each stage's `manifest.json` | State, provenance, fingerprints and errors |

Each stage needs a new/empty output directory. Evidence is sealed; edit the decision TSVs, not the saved models/windows. Partial models and shortened flanks require explicit acknowledgement in the decisions. Human selection records a review decision and does not independently prove orthology.

## Reuse existing work

```bash
gene-context review-search --results previous/search --out results/01_review
gene-context review-contexts --results previous/comparison --out results/02_review
```

These generate review checkpoints without repeating search or prediction. Imported automatic selections remain unapproved. See [reuse and decision formats](docs/STAGED_WORKFLOW.md).

The existing `analyze` and `compare` commands still provide uninterrupted diagnostic runs and bypass manual review stops. `run` and `align` remain available for extraction/ordinary multiple alignments. See [automatic comparison](docs/COMPARATIVE_ANALYSIS.md).

**Validation:** the [18-assembly NANOG/POU5 audit](docs/VERTEBRATES18_VALIDATION.md) exposed other-family selections, copy ambiguity and incomplete models. This is why explicit review precedes comparative interpretation. The staged workflow adds review gates; automatic family/copy classification remains future work.

## More detail

- [Manual stages, decision formats and resuming](docs/STAGED_WORKFLOW.md)
- [Setup, local inputs and optional downloads](docs/SETUP_AND_DOWNLOAD.md)
- [Comparative analysis and report interpretation](docs/COMPARATIVE_ANALYSIS.md)
- [Search methods, options and file reference](docs/METHODS.md)
- [Tests and validation](docs/VALIDATION.md)
- [Tiberius integration and end-to-end validation](docs/TIBERIUS_INTEGRATION_VALIDATION.md)
- [Measured Tiberius CPU runtime, RAM and model sizes](docs/TIBERIUS_CPU_BENCHMARK.md)
- [Five-assembly protein/CDS pilot](docs/FULL_GENOME_PILOT.md)

This standalone project grew out of [HIST_Cluster_Finder](https://github.com/rljech13/HIST_Cluster_Finder).
