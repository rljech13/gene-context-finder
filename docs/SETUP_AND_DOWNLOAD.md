# Setup, genome downloads, and first run

## 1. Get this implementation and create the environment

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
gene-context doctor
```

Run environment creation from the repository root so the editable Python package
(`pip: -e .`) resolves correctly. `mamba env create -f environment.yml` is an
alternative when mamba is installed. The environment contains Python, samtools,
miniprot, BLAST+ and MAFFT; the application itself has no Python dependencies.
Use `conda run -n gene-context gene-context ...` in batch jobs if activation is
inconvenient. The supported platforms are Linux/macOS with available tool builds;
miniprot needs x86_64 or ARM hardware. Windows users can use Linux under WSL.

A Linux-64 dependency solve was checked for this environment file. The exact
packages used in a real run are recorded in `manifest.json`; version constraints
are not a fully locked environment. When using existing system tools, Python
>=3.10 and all executables listed by `gene-context doctor` must be available on
PATH. No activation script or global environment is modified by the pipeline.

## 2. Download a small subset first

The supplied [18-species manifest](../examples/vertebrates18.download.tsv) lists
exact assembly versions and their HTTPS genome/annotation/reference URLs. It is
an optional example panel; the pipeline accepts your own genomes and sequences.
The panel was frozen on 2026-11-05. It contains 15 Chromosome assemblies, two
Complete Genome assemblies and one Scaffold assembly. The two deliberately
unannotated exceptions are *Lethenteron camtschaticum* and *Pavo cristatus*.
No GFF3 from a different assembly is silently substituted for these species.

Preview chicken without downloading anything:

```bash
gene-context download --manifest examples/vertebrates18.download.tsv \
  --out data/vertebrates18 --only GCF_016699485.2 --dry-run
```

Download and verify chicken genome + matching GFF3:

```bash
gene-context download --manifest examples/vertebrates18.download.tsv \
  --out data/vertebrates18 --only GCF_016699485.2
```

Select more than one assembly using `--only ACCESSION1 ACCESSION2`. Download all
18 by omitting `--only`. Add `--include-reference` if you also want the listed
complete reference CDS and protein FASTA files; they are not needed for the main
protein-to-genome search.

```bash
gene-context download --manifest examples/vertebrates18.download.tsv \
  --out data/vertebrates18 --include-reference
```

The complete panel has approximately 30.1 GB of compressed genome data and
98.4 Gb of genomic sequence, before adding annotations, decompressed working
FASTA, miniprot/BLAST indices and results. The lungfish and axolotl account for
much of this. Start with a subset rather than downloading the full panel onto a
small disk. The main search processes assemblies sequentially with configurable
threads; index memory and maximum intron length matter for giant genomes.

The downloader:

- retrieves official `md5checksums.txt` for each file's source directory;
- downloads to `.part`, resumes when HTTP Range is supported, and publishes a
  final filename only after MD5 agreement;
- reuses existing files only after rechecking their MD5;
- preserves an existing file that disagrees with its source checksum and reports
  an error instead of overwriting it;
- writes `download_receipt.json` with source URLs, checksums, byte counts and
  outcomes, plus `genomes.tsv` containing only fully verified selected entries;
- records failures and continues with the other selected assemblies.

Re-running a subset regenerates `genomes.tsv` for that invocation's selected,
verified genomes. It does not automatically add every older cached download.
The downloaded genomes and reference files are excluded from Git by `.gitignore`.
Genome records are version-pinned; the recorded MD5 also captures the annotation
bytes actually obtained if NCBI updates annotation on the same assembly later.

For your own download panel, create a tab-separated file with `genome_id` and
`genome_url`; optional columns are `accession`, `gff_url`, `protein_url`, `cds_url`
and `md5_url`. The default checksum URL is `md5checksums.txt` beside each source
file. A custom source needs a compatible checksum file. An empty GFF URL means
annotation is unavailable, not an unsuccessful download.

## 3. Supply query sequences

Use a FASTA you already have, or download known reference sequences. For example,
these chicken NANOG and Pou5f3 RefSeq proteins are suitable demonstration queries:

```bash
mkdir -p queries
curl -fL 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_001139614.1,NP_001296301.2&rettype=fasta&retmode=text' \
  -o queries/chicken_targets.faa
```

The records correspond to [NANOG, Gene 100272166](https://www.ncbi.nlm.nih.gov/gene/100272166)
and [Pou5f3, Gene 427781](https://www.ncbi.nlm.nih.gov/gene/427781).
For known spliced coding DNA, obtain their transcript CDS:

```bash
curl -fL 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_001146142.2,NM_001309372.2&rettype=fasta_cds_na&retmode=text' \
  -o queries/chicken_targets.cds.fna
```

Check sequence identifiers and that CDS translation agrees with the intended
protein. `--query-type cds` expects intron-free coding DNA in coding orientation,
not an entire mRNA containing UTRs. Use `--query-type dna` for genomic or noncoding
queries. Arbitrary genomic DNA without a coding annotation is retained as matched
regions and flanks; it does not establish a uniquely defined coding model.

## 4. Run on downloaded or already available genomes

```bash
gene-context run --query queries/chicken_targets.faa --query-type protein \
  --genomes data/vertebrates18/genomes.tsv --out results/protein_run --threads 8

gene-context run --query queries/chicken_targets.cds.fna --query-type cds \
  --genomes data/vertebrates18/genomes.tsv --out results/cds_run --threads 8
```

If genomes are already available, do not download them again. Write a TSV with
`genome_id`, `fasta`, and optional `gff3`, using their actual paths on the machine
where you will run the search. Paths may be absolute or relative to the TSV.
Keep different assembly versions separate. A `.fai` file is only an index and
cannot replace its FASTA; a locus/chromosome excerpt is not a complete genome.

```tsv
genome_id	fasta	gff3
my_assembly	/path/to/genome.fna	/path/to/matching_assembly.gff3.gz
another_assembly	/path/to/another.fna.gz
```

For entirely unannotated FASTA collections, `--genome-dir /path/to/genomes`
recursively uses every `.fa`, `.fna`, `.fasta` and gzipped equivalent. Keep query
FASTA files outside that directory. Prefer the explicit TSV when annotations,
multiple assembly versions or unusual directory layouts are involved.

Open `results/protein_run/report.html` and inspect `genome_query_summary.tsv`.
Every input genome/query has an outcome, including no hits and errors. The output
directory must be new or empty. The original genome/annotation files are not
modified. Copy the full result directory when sharing the HTML report and SVGs.

## Alternative: NCBI Datasets CLI

The [official install guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
and [genome download guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
describe the NCBI package downloader. It is an alternative to the built-in
checksum-verified URL manifest workflow above:

```bash
conda install -c conda-forge ncbi-datasets-cli
datasets download genome accession GCF_016699485.2 \
  --include genome,gff3,cds,protein --filename chicken_dataset.zip
unzip chicken_dataset.zip -d data/chicken_dataset
```

Locate the assembly FASTA and `genomic.gff` in the extracted `ncbi_dataset/data/`
assembly directory, then add their exact paths to an input TSV. For very large
packages follow NCBI's dehydrated-download/rehydration guide instead of treating a
metadata-only package as downloaded sequence data.
