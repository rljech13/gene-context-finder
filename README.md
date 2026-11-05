# Gene Context Finder

Find candidate genes and their surroundings in a set of genomes using DNA or protein sequences.

```text
Your query FASTA + genome FASTAs + optional matching GFF3 annotations
                              ↓
              Search every query in every genome
                              ↓
        Extract candidate CDS, proteins, loci and flanks
                              ↓
             Reports and sequences; optional alignments
```

## What to prepare

- **Query FASTA:** one or more proteins, CDS sequences or DNA regions. Use one sequence type per file. Gene names alone are not sequence input.
- **Genome FASTA:** one file per assembly; `.fa`, `.fna`, `.fasta` and gzip files are accepted.
- **GFF3 annotation, if available:** it must match the exact genome assembly. Gzip is accepted.

Use already downloaded files directly, with their original filenames.

## Install once

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
gene-context doctor
```

## 1. List your local genomes

Create `genomes.tsv` with **tab-separated** columns. Replace the example paths with your actual filenames:

```tsv
genome_id	fasta	gff3
chicken	data/chicken/genome.fna.gz	data/chicken/annotation.gff.gz
lamprey	data/lamprey/genome.fna.gz
```

Add one row per assembly. Paths are relative to `genomes.tsv`, or absolute. Leave `gff3` empty when annotation is unavailable. Gzip inputs need no manual decompression.

## 2. Search using your proteins

Put your real query proteins in `targets.faa`, then run:

```bash
gene-context run --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/search --threads 8 --skip-alignments
```

This searches every query in every genome, checks the recovered coding models, and extracts CDS, proteins and flanks (5 kb per side by default).

**Without GFF3:** protein queries still recover predicted coding models through miniprot. Neighbouring genes are reported from GFF3 when available; the current pipeline does not independently annotate all genes in the surrounding DNA. Flanks are measured from coding-model boundaries, not a measured transcription start site.

For spliced coding DNA without introns or UTRs, use `--query-type cds`. For genomic/noncoding DNA, use `--query-type dna`; these matches can remain DNA regions when coding evidence is unavailable.

Use a new or empty output directory for each run. Gzip genomes are decompressed inside that directory, which also holds search indexes; allow additional disk space beyond the compressed downloads.

## 3. Open the results

Start with `results/search/report.html`.

| File or folder | What it contains |
|---|---|
| `genome_query_summary.tsv` | Outcome for every genome/query pair, including no hits and errors |
| `models.tsv` | Candidate coordinates, evidence and complete/partial/disrupted status |
| `sequences/q0001/` | CDS, proteins, loci and flanks for the first query across genomes |
| `genomes/` | Coding models in GFF3 and context diagrams for each genome |
| `manifest.json` | Query names, input checksums, settings, commands and errors |

Queries become `q0001`, `q0002`, etc., in FASTA order. **Each query has its own results:** several POU5F3 references can recover the same locus in separate groups. Candidates include possible isoforms and paralogs; review them before counting genes or comparing orthologs.

## Optional: align the extracted sequences

Alignments support sequence-conservation comparisons. To add them after extraction:

```bash
gene-context align --results results/search --out results/aligned --threads 8
```

`results/aligned/alignments/` contains per-query protein, codon-CDS and separate upstream/downstream alignments, plus identity/coverage tables. All complete candidates enter by default, so review the candidate set before interpretation.

Omit `--skip-alignments` to search and align in one run. The separate `align` command reuses saved sequences.

## More detail

- [Setup, local inputs and optional downloads](docs/SETUP_AND_DOWNLOAD.md)
- [Search methods, options and file reference](docs/METHODS.md)
- [Tests and validation](docs/VALIDATION.md)
- [Five-assembly protein/CDS pilot](docs/FULL_GENOME_PILOT.md)

This standalone project grew out of [HIST_Cluster_Finder](https://github.com/rljech13/HIST_Cluster_Finder).
