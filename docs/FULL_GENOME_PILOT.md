# Full-genome pilot — 2026-11-05

Subsequent correction: v0.4.1 recovered the previously missed NANOG loci using independent TBLASTn discovery with the same chicken query. See [NANOG diagnosis and validation](NANOG_SEARCH_VALIDATION.md). The historical counts below describe the earlier source version.

Two chicken references were searched against five whole-assembly FASTA inputs with their matching NCBI GFF3 annotations. Each assembly was processed in both protein-input and coding-DNA-input modes, using the default discovery gates, 16 threads, 5-kb oriented flanks and a 50-kb neighbour radius. These are candidate-model counts, including isoforms and other homologous genes, not ortholog or biological-copy counts.

## Queries

- NANOG: [NP_001139614.1](https://www.ncbi.nlm.nih.gov/protein/NP_001139614.1), 309 aa; [NM_001146142.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_001146142.2), CDS 930 nt including stop.
- Pou5f3: [NP_001296301.2](https://www.ncbi.nlm.nih.gov/protein/NP_001296301.2), 389 aa; [NM_001309372.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_001309372.2), CDS 1,170 nt including stop.

The transcript CDS translations were checked against the corresponding RefSeq proteins before execution.

## Assemblies and observed models

Counts below are all retained candidates, with structurally complete predicted models in parentheses. Partial/disrupted models remain in the inventory and do not enter the default alignments.

| Species / assembly | Protein: NANOG | Protein: Pou5f3 | CDS: NANOG | CDS: Pou5f3 |
|---|---:|---:|---:|---:|
| Lepisosteus oculatus / [GCF_040954835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040954835.1/) | 0 (0) | 13 (5) | 0 (0) | 33 (25) |
| Danio rerio / [GCF_052040795.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_052040795.2/) | 0 (0) | 11 (3) | 0 (0) | 106 (98) |
| Monodelphis domestica / [GCF_027887165.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_027887165.2/) | 0 (0) | 12 (5) | 0 (0) | 21 (14) |
| Anolis carolinensis / [GCF_035594765.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_035594765.1/) | 2 (1) | 13 (4) | 2 (1) | 36 (27) |
| Chelonia mydas / [GCF_015237465.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015237465.2/) | 4 (3) | 13 (3) | 4 (3) | 29 (19) |

## Technical outcomes

| Mode | Genome/query rows | Candidate models | Complete predicted | Partial | Disrupted | Wall time | Maximum RSS |
|---|---:|---:|---:|---:|---:|---|---:|
| protein | 10 | 68 | 24 | 22 | 22 | 6:56.09 | 35.8 GiB |
| cds | 10 | 231 | 187 | 22 | 22 | 15:16.95 | 35.8 GiB |

The protein run finished with exit status zero. Both discovery/extraction passes recorded an outcome for every genome/query without an assembly-stage error. The initial CDS pass reached MAFFT auto refinement for 183 Pou5f3-query flank sequences; that expensive alignment subtree was deliberately stopped, leaving its original manifest with an alignment error and exit status one. Gene searches and extraction were not repeated. The saved extraction was then aligned successfully with the v0.1.1 `align` command and progressive DNA mode. The original and resumed manifests are preserved separately. The initial CDS timing above includes the stopped auto attempt and is not a measurement of an end-to-end fast run. RSS is the GNU time maximum including child processes; timing depends on hardware, filesystem/cache state and candidate-set size.

## Independent checks

The output audit used Biopython 1.88 for translation and samtools faidx on the original input genomes, independently of the pipeline translation and FASTA extraction code. It reconstructed complete CDS from exported exon/stop coordinates and checked oriented flanks and codon alignment roundtrips.

| Mode | CDS/protein translation checks | Complete CDS re-extracted from source | Oriented flanks checked | Codon alignments checked | Failures |
|---|---:|---:|---:|---:|---:|
| protein | 68 | 24 | 48 | 24 | 0 |
| cds | 231 | 187 | 374 | 187 | 0 |

The separate progressive alignment stage finished in 4:09.46, with maximum RSS 8.0 GiB and exit status zero. All 187 complete candidate models were retained; the protein/codon alignments also include the two reference sequences. A separate Biopython audit checked 563 ungapped protein/flank rows and 189 codon translations, with zero failures. Sequence comparisons ignore FASTA letter case: MAFFT emits nucleotide rows in lowercase. This validates content preservation, not biological homology of every aligned flank.

## Interpretation limits

The Pou5f3 query also recovered annotated POU1/POU2/POU3-family genes and POU5F1 models. These are retained as sequence-supported candidates; the pipeline does not rename them Pou5f3 or resolve individual-copy orthology. Multiple models at one locus can be isoforms or alternative predictions.

The permissive nucleotide association gate (10% query-HSP-span union) deliberately retains fragment-supported transcript candidates, including many other homeobox/POU genes. These broad candidates require curation before target-specific conservation inference; `--min-dna-coverage` can be raised for a narrower nucleotide search. Structural CDS completeness is a separate question from homology specificity. The single chicken NANOG protein did not produce a passing model in spotted gar, zebrafish or opossum at the default protein discovery gates. That is a search-sensitivity result, not evidence of gene absence. Additional taxon-appropriate references or a separately validated profile-based discovery strategy may be needed for divergent targets. The table records the coding-DNA route independently rather than extrapolating from the protein route.

Full models, alignments, contexts, input checksums, source-file hashes and logs are retained with the run. Search and extraction source hashes are retained with the initial run; the separate realignment manifest fingerprints the v0.1.1 alignment implementation and its input sequences. The downloader was not needed for this cached-input pilot; its source-checksum parsing is covered by separate tests. Raw genome files and run-specific filesystem paths are not bundled into this repository.

## Reproduce the input set

Use `gene-context download` with `examples/vertebrates18.download.tsv` and select:

```text
GCF_040954835.1 GCF_052040795.2 GCF_027887165.2 GCF_035594765.1 GCF_015237465.2
```

Follow [setup/download instructions](SETUP_AND_DOWNLOAD.md) for the reference sequences, then run once with `--query-type protein` and once with `--query-type cds`, using `--threads 16 --flank 5000 --context 50000` and separate new output directories. Current v0.1.1 defaults to progressive DNA alignment; `gene-context align` can also apply it to an existing extraction without repeating discovery. The initial v0.1.0 timing/auto-refinement attempt is retained above as a diagnostic, not a requirement for reproduction.
