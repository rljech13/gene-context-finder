# Setup and input preparation

The [README](../README.md) shows the usual search workflow. Start with files you or a colleague have already downloaded. The downloader below is optional.

## Install

Run from a terminal on the machine where you want the analysis to run:

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
gene-context setup-predictor --with-reference
gene-context doctor
```

Create the environment from the repository root so `pip: -e .` resolves correctly. Mamba can replace conda. The environment provides Python 3.12, Git, samtools, miniprot, BLAST+ and MAFFT. `setup-predictor` installs the pinned Tiberius 2.0.9 CPU runtime in a separate virtual environment and prepares the two shared models. It avoids TensorFlow CUDA extras; no GPU is required. Existing installations can use `conda env update -f environment.yml` and then run setup.

With existing tools, use Python >=3.10 and install the package with `python -m pip install -e .`. Alternatively, run `python -m genecontext` from this checkout. External tools must be on PATH. Linux/macOS are supported where tool builds are available; use Linux under WSL on Windows. Miniprot needs x86_64 or ARM hardware. See [validation](VALIDATION.md) for tested environments.

## Predictor and naming caches

The default predictor cache is `~/.cache/gene-context/tiberius`: `runtime/` contains Python packages; `models/vertebrates/` and `models/mammalia/` contain verified weights and archives. Both model archives total 66,174,115 bytes; extracted model files total 71,436,162 bytes. Allow about 2 GB on macOS or 3 GB on Linux for the runtime and models, excluding genomes/results. These two models are shared across species. Setup writes `setup.log` and a runtime receipt, and checks the model archives and extracted files against pinned hashes.

Python 3.12 or 3.13 and Git are required for setup. If the driver uses another Python, pass `setup-predictor --python /path/to/python3.12`. The main package itself remains Python >=3.10 with no Python dependencies. The predictor is pinned to upstream commit `b3ea48f263a075260dbf1e91a412e82ccdc96d3c`; its independent dependency lock is shipped with the package. Install a new release in a fresh `--predictor-cache` directory if its runtime requirements change.

Swiss-Prot is separate, under `~/.cache/gene-context/references/swissprot`. `--with-reference` prepares it once; analysis otherwise prepares it on first use. The verified cached release stays pinned. The pipeline builds naming BLAST databases inside each comparison output, not in your input folders.

For custom locations, use matching flags during setup **and analysis**:

```bash
gene-context setup-predictor --predictor-cache /data/gcf-predictor \
  --reference-cache /data/gcf-references --with-reference
gene-context analyze --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/analysis --threads 4 --offline \
  --predictor-cache /data/gcf-predictor --reference-cache /data/gcf-references
```

`--offline` rejects missing caches rather than downloading them. Model/reference caches can be copied; create the Python runtime on the destination machine, since virtual environments are not portable between paths or operating systems. `setup-predictor --models-only` prepares weights without installing Python packages. `--tiberius-python /path/to/prepared/venv/bin/python` selects an existing compatible runtime.

For optional AUGUSTUS: `conda install -c conda-forge -c bioconda augustus`, then select `--predictor augustus --augustus-species MODEL`; per-genome `augustus_species` values override that model. See [window prediction](WINDOW_PREDICTION.md).

## Use files downloaded by a colleague

For each assembly, collect its genome FASTA, its matching GFF3 if available, and the assembly accession/version. Keep each assembly in its own folder; original filenames can be retained. `.gz` files are accepted directly. Convert GTF to GFF3 before using it.

Create `genomes.tsv` with tab-separated columns:

```tsv
genome_id	fasta	gff3	tiberius_model
chicken	data/chicken/genome.fna.gz	data/chicken/genomic.gff.gz	vertebrates
human	data/human/genome.fna.gz		mammalia
```

Use unique genome IDs. Paths are relative to this TSV, or absolute. Empty `gff3` means no annotation. The GFF3 must match the exact assembly and FASTA contig IDs. A `.fai` index does not replace its FASTA. A chromosome or locus excerpt is processed as the supplied sequence, not as a whole genome.

`tiberius_model` selects shared model weights, not a species-specific download. Omitted values use `vertebrates`; choose `mammalia` for mammals. You can instead pass `--prediction-models models.tsv` (`genome_id`, `model`). For custom unmasked Tiberius models, use a directory path relative to the TSV. Every selected window is predicted from DNA before its proteins are identified.

The search reads this local-path table and does not download genome files. Run on the machine that has these paths: files on a colleague's computer must first be copied there or made accessible. Inputs are read-only; working FASTA and indexes are created under the output directory. Default protein naming uses cached Swiss-Prot, downloaded once if missing. `--offline` prevents downloads; `--context-reference none` uses only local naming inputs while retaining predicted unnamed genes.

## Prepare the query FASTA

Use your real protein sequences, CDS sequences or genomic DNA regions. A file can contain multiple records of one type. Give each record a unique FASTA identifier. The pipeline searches each record independently in every genome.

To obtain the chicken NANOG and Pou5f3 example proteins:

```bash
mkdir -p queries
curl -fL 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_001139614.1,NP_001296301.2&rettype=fasta&retmode=text' \
  -o queries/chicken_targets.faa
```

These are [NANOG](https://www.ncbi.nlm.nih.gov/gene/100272166) and [Pou5f3](https://www.ncbi.nlm.nih.gov/gene/427781). Known reference sequences are search examples; a single reference may miss divergent homologs. For their coding DNA:

```bash
curl -fL 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_001146142.2,NM_001309372.2&rettype=fasta_cds_na&retmode=text' \
  -o queries/chicken_targets.cds.fna
```

Check identifiers and CDS/protein agreement. `--query-type cds` expects spliced coding DNA without UTRs, not the entire mRNA. Protein/CDS queries must encode at least 30 amino acids. For genomic or noncoding queries, use `--query-type dna`.

Run with the local genome table:

```bash
gene-context analyze --query queries/chicken_targets.faa --query-type protein \
  --genomes genomes.tsv --out results/analysis --threads 8
```

Open `results/analysis/report.html`. To share the report with its linked sequences and diagrams, copy the full result directory. See [comparative analysis](COMPARATIVE_ANALYSIS.md) for the synteny/DNA/protein report and [methods](METHODS.md) for the underlying search rules.

## Optional: download through the pipeline

The downloader reads a **different table**, containing source URLs rather than local file paths. The supplied [18-assembly download table](../examples/vertebrates18.download.tsv) already contains exact NCBI URLs.

```text
URL table → download → local files + genomes.tsv → run
```

Preview one assembly without downloading:

```bash
gene-context download --manifest examples/vertebrates18.download.tsv \
  --out data/vertebrates18 --only GCF_016699485.2 --dry-run
```

Download that assembly's genome and available GFF3:

```bash
gene-context download --manifest examples/vertebrates18.download.tsv \
  --out data/vertebrates18 --only GCF_016699485.2
```

Omit `--only` for all 18; add `--include-reference` for available protein and CDS files as well. Those extra files are stored as references; the main search reads genome/GFF3 and the separate query FASTA.

The files land in `data/vertebrates18/` on the machine where the command runs. The downloader checks official MD5 sums, supports partial-download resumption where the server allows it, rechecks cached files, and preserves existing files that fail verification. It writes `download_receipt.json` with URLs, byte counts and errors, plus `genomes.tsv` listing verified entries. It continues past failed assemblies; check the receipt before assuming the panel is complete. Re-running a subset regenerates the local-path table for that invocation's selected, verified entries.

Use the generated table in the search:

The downloader supplies genome paths. Use the matching model table to select the two shared models:

```bash
gene-context analyze --query queries/chicken_targets.faa --query-type protein \
  --genomes data/vertebrates18/genomes.tsv --out results/downloaded_analysis \
  --prediction-models examples/vertebrates18.models.tsv --threads 8
```

For a custom download table, required columns are `genome_id` and `genome_url`; optional columns are `accession`, `gff_url`, `protein_url`, `cds_url` and `md5_url`. Sources need a compatible official checksum list. Without an explicit checksum URL, the downloader uses `md5checksums.txt` beside each source file.

### The example panel and storage

The panel was frozen on 2026-11-05: 18 assemblies, with 15 Chromosome, two Complete Genome and one Scaffold assembly. The selected *Lethenteron camtschaticum* and *Pavo cristatus* packages lack GFF3, protein and CDS files.

HTTP file-size checks on that date gave **31.029 GB of gzip inputs** for all available genome/GFF3/protein/CDS files: 18 genomes and 16 files of each other type. Genomes alone total 30.128 GB. GB means 10^9 bytes. This excludes decompressed working genomes, indexes and results. The assemblies contain about 98.4 billion genomic bases; lungfish and axolotl dominate the storage and index requirements.

### Alternative: download from NCBI directly

Manual downloads and the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/) produce usable inputs too:

```bash
datasets download genome accession GCF_016699485.2 \
  --include genome,gff3,cds,protein --filename chicken_dataset.zip
unzip chicken_dataset.zip -d data/chicken_dataset
```

The NCBI CLI is a separate installation; see its [installation guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/). Locate the assembly FASTA and GFF3 in the extracted package and list their actual paths in your local `genomes.tsv`.
