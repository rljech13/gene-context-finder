# Gene Context Finder

Find candidate genes and their surroundings in a set of genomes using DNA or protein sequences.

```text
Your query FASTA + genome FASTAs + optional matching GFF3 annotations
                              ↓
              Search every query in every genome
                              ↓
         Extract candidate loci and ±50 kb gene contexts
                              ↓
        Predict genes from the DNA of every context window
                              ↓
        Identify predicted proteins against a reference library
                              ↓
       Compare neighbours, genomic DNA and amino acids
                              ↓
                 One comparative report per target
```

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

Setup creates an isolated CPU runtime and downloads both shared models (66.2 MB compressed in total), plus Swiss-Prot for protein identification. No GPU is required. Allow about 2 GB for the predictor installation, in addition to genomes and results. See [setup](docs/SETUP_AND_DOWNLOAD.md) for cache locations and offline use.

## 1. List your local genomes

Create `genomes.tsv` with **tab-separated** columns. Replace the example paths with your actual filenames:

```tsv
genome_id	fasta	gff3	tiberius_model
chicken	data/chicken/genome.fna.gz	data/chicken/annotation.gff.gz	vertebrates
human	data/human/genome.fna.gz		mammalia
```

Add one row per assembly. Paths are relative to `genomes.tsv`, or absolute. Leave `gff3` empty when annotation is unavailable. Gzip inputs need no manual decompression.

`tiberius_model` is optional: omitted entries use `--tiberius-model vertebrates`. Select `mammalia` explicitly for mammals. A separate `--prediction-models models.tsv` table can override models, including for saved searches; the [18-species example](examples/vertebrates18.models.tsv) matches the download manifest. No taxonomy is inferred from filenames.

## 2. Analyze your genes

A ready-to-use chicken Pou5f3 query is included: [protein FASTA](examples/queries/chicken_Pou5f3.protein.faa), with [CDS and provenance](examples/queries/README.md).

Put one or more query proteins in `targets.faa`, then run:

```bash
gene-context analyze --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/analysis --threads 8
```

This searches every query in every genome, checks target coding models, and extracts each selected locus with **50 kb on either side**. It then predicts genes from each window's DNA with Tiberius on CPU, identifies their proteins, compares neighbours, and compares genomic DNA and target proteins. The prediction step runs for **every genome**, including genomes with GFF3. Use `--window` to change the context size. Add `--include-partial` to include clean incomplete target proteins; disrupted target models stay excluded.

**Without GFF3:** the same DNA-only prediction step runs. Genes remain on the diagram even if their proteins cannot be identified. BLASTP names predicted proteins against Swiss-Prot (default), supplied `--context-proteins reference.faa`, and any usable proteins from the provided annotations. GFF3 supplies an independent comparison/reference track, without constraining prediction. Unnamed predictions can still match across genomes by their protein sequences. See [window prediction](docs/WINDOW_PREDICTION.md).

After setup, add `--offline` to require cached models/references and prevent downloads during analysis. Use `--context-reference none` to identify only against your local proteins; prediction still runs.

For spliced coding DNA without introns or UTRs, use `--query-type cds`. For genomic/noncoding DNA, use `--query-type dna`; these matches can remain DNA regions when coding evidence is unavailable.

Use a new or empty output directory for each run. Gzip genomes are decompressed inside that directory, which also holds search indexes; allow additional disk space beyond the compressed downloads.

## 3. Open the results

Start with `results/analysis/report.html`, then open a target's page for its **synteny diagram, DNA blocks and protein comparison**.

| File or folder | What it contains |
|---|---|
| `search/` | Every search candidate, extracted CDS/proteins/flanks and original evidence |
| `comparison/summary.tsv` | Outcome for every genome/target group |
| `comparison/candidates.tsv` | All loci, isoforms, ranks and the displayed representative |
| `comparison/loci.tsv`, `comparison/windows/` | Selected coordinates and actual context DNA |
| `comparison/neighbors.tsv`, `comparison/synteny.tsv` | Neighbour homology groups, order and orientation |
| `comparison/context_predictions/`, `comparison/context_predictions.tsv` | Native-coordinate predicted genes, including unnamed models |
| `comparison/context_prediction_cds.fna`, `comparison/context_prediction_proteins.faa` | CDS/proteins from window predictions |
| `comparison/provided_neighbors.tsv` | Independent neighbours from supplied GFF3 |
| `comparison/context_reference.json` | Identification source, pinned database release and checksums |
| `comparison/dna_pairs.tsv`, `comparison/dna_categories.tsv` | DNA identity and aligned coverage, including region-specific results |
| `comparison/queries/` | Target pages, SVG diagrams, protein/codon alignments and comparison tables |
| `comparison/manifest.json` | Input/window hashes, parameters, versions, commands and errors |

Queries become `q0001`, `q0002`, etc., in FASTA order. Use an optional query-group TSV to combine several references of one gene. Multiple candidate loci remain visible; the displayed representative is ranked by sequence evidence, with coding integrity breaking ties. Neighbour families are provisional protein-homology groups, not proven orthology assignments.

## Reuse an existing search

Build the comparative report without searching the genomes again:

```bash
gene-context compare --results results/search --out results/comparison --threads 8 \
  --prediction-models examples/vertebrates18.models.tsv
```

`analyze` combines search and comparison. The existing `run` and `align` commands remain available for search/extraction or ordinary multiple alignments. See [comparative analysis](docs/COMPARATIVE_ANALYSIS.md) for grouping, manual locus selection, DNA engines and interpretation.

The v0.2 annotation-first context workflow is still available explicitly with `--context-mode annotation`.

## More detail

- [Setup, local inputs and optional downloads](docs/SETUP_AND_DOWNLOAD.md)
- [Comparative analysis and report interpretation](docs/COMPARATIVE_ANALYSIS.md)
- [Search methods, options and file reference](docs/METHODS.md)
- [Tests and validation](docs/VALIDATION.md)
- [Tiberius integration and end-to-end validation](docs/TIBERIUS_INTEGRATION_VALIDATION.md)
- [Measured Tiberius CPU runtime, RAM and model sizes](docs/TIBERIUS_CPU_BENCHMARK.md)
- [Five-assembly protein/CDS pilot](docs/FULL_GENOME_PILOT.md)

This standalone project grew out of [HIST_Cluster_Finder](https://github.com/rljech13/HIST_Cluster_Finder).
