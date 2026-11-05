# Tiberius CPU integration validation — v0.4.0

Date: 2026-11-05. This validates operation, sequence extraction and report generation. It is not a gene-prediction accuracy benchmark, an orthology study, or a run of the complete 18-species panel.

## Tested installation

The actual `setup-predictor --with-reference` command created an isolated Python 3.12.9 runtime on macOS arm64 (Apple M3 Pro, 36 GiB RAM), installed the shipped CPU dependency lock, and prepared both shared model directories. Model archives were reused from the independent benchmark cache, then extracted and verified by the new backend. The previously verified Swiss-Prot cache was reused. No CUDA extras were installed and no GPU was visible during inference.

The predictor is Tiberius 2.0.9, upstream commit `b3ea48f263a075260dbf1e91a412e82ccdc96d3c`, TensorFlow 2.20.0, bricks2marble 0.1.3 and hidten 1.0.0. Runtime package versions, installed Tiberius source hashes, model hashes and CPU-visible devices are included in every prediction manifest. The driver remains dependency-free; the built wheel was checked to include the independent requirements lock and worker.

`doctor` passed with actual samtools 1.24, miniprot 0.18-r281, BLAST+ 2.17.0+, MAFFT 7.526, the managed Tiberius runtime and cached weights. AUGUSTUS remains optional; its existing real-tool tests are retained.

## Automated checks

All **53 tests passed locally**, with no skips, in 282.497 seconds. This includes the optional LASTZ and AUGUSTUS controls and all three real Tiberius integration tests.

The complete local suite is run with real sequence tools and prepared Tiberius weights; integration tests do not mock sequence search or neural inference. Missing executables or runtime caches cause explicit skips in environments without those dependencies. GitHub CI covers the dependency-free suite on Python 3.10–3.13; it does not install TensorFlow or run the neural integration tests.

The seven added tests cover:

- Local HTTP model transfer, pinned archive and extracted-file checksums, offline reuse, preservation of corrupt inputs, and rejection of archive path traversal.
- Model-map overrides, relative custom-model paths and preservation of a virtual-environment Python executable.
- Complete protein-query analysis of two opposite-strand fixture genomes, including an irrelevant GFF3 track, a genome without GFF3, a negative query, shared CPU batching, PTGS2 identification, DNA/protein comparison and native CDS/protein roundtrips.
- DNA-only seeds and different per-genome models producing separate CPU batches.
- Identification failure preserving unnamed predictions and the report, plus explicit skipping of windows shorter than 1,000 bp.

The integration fixture combines the official AUGUSTUS HS04636 human PTGS2 sequence with chicken Pou5f3 CDS and N spacers. It is an artificial computational control, not a biological genomic arrangement. See [fixture provenance](../tests/fixtures/README.md).

## Real-sequence annotation-free pilot

Two chicken protein queries (NANOG and Pou5f3) were searched afresh against five FASTA inputs containing seven real genomic fragments (706,653 bases in total): *Anolis carolinensis*, *Chelonia mydas*, *Lepisosteus oculatus*, *Danio rerio* and *Monodelphis domestica*. These fragments came from the earlier [five-assembly pilot](FULL_GENOME_PILOT.md); native orientations were restored before searching. GFF3 was omitted from every input. The mammalian window used `mammalia`; the other windows used `vertebrates`. The run used 4 CPU threads, `--include-partial`, cached Swiss-Prot and `--offline`.

This deliberately exercises annotation-free operation. It cannot measure discovery sensitivity in complete genomes, and the selected windows/sequence-ranked hits are not all confirmed NANOG/Pou5f3 orthologues.

| Output | Observed |
|---|---:|
| Input genome/target combinations | 10 |
| Search coding candidates | 9 |
| Selected candidate windows | 7 |
| Tiberius CPU model groups | 2 |
| Predicted coding models | 22 |
| Predicted neighbours after target-overlap filtering | 15 |
| Neighbours with tentative protein-homology names | 11 |
| Unnamed neighbours retained | 4 |
| Compared DNA window pairs | 11 |
| Synteny pairs with insufficient unique anchors | 10 |
| Synteny pairs with assessable order/orientation | 1 |
| Pipeline errors | 0 |

Four unnamed neighbours comprise three without a passing identification and one with conflicting top names. They remain in the context output. Three of the ten input genome/target combinations had no detected coding candidate with this query/search configuration; this is not evidence of biological gene absence. Only clean complete/partial target models enter the requested amino-acid/codon comparisons; disrupted target candidates remain visible and excluded from those alignments.

Independent checks use BioPython without importing `genecontext`: raw input DNA → oriented window → native predicted CDS → translated protein, phase continuity, codon-alignment backtranslation, HTML link resolution, SVG XML parsing and exact source-module fingerprints. The report's SVGs were rendered and visually inspected. The [independent audit script](benchmarks/audit_tiberius_integration.py) and [recorded result](benchmarks/tiberius_integration_2026-11-05.json) are included. All 7 windows, 22 raw models, 22 native-coordinate models, 4 codon-alignment sequences, 39 report links, 6 SVGs and 17 exact source-module fingerprints passed. The final run took 0.47 s for fragment search and 149.95 s for comparison; setup, downloads and the earlier full-genome searches are excluded, and other local tests ran concurrently.

## Linux host1 check against full-assembly search results

The same managed setup succeeded on Linux x86_64 with Python 3.12.14 and the same dependency lock. The server has AMD EPYC 7763 processors; inference was limited to four threads and CPU-only visibility. The installed runtime occupied 2,286,232,070 bytes on disk, excluding models/references/genomes, versus about 1.4 GB on the tested Mac. Installation guidance therefore allows about 3 GB on Linux.

`compare` reused the original five-whole-assembly protein-search results; target discovery was **not rerun**. It extracted the selected loci from the original full genome FASTAs, retained the supplied GFF3 as an independent reference track, ran both CPU model groups, identified proteins against Swiss-Prot plus the provided annotations, and generated all comparisons offline. It used the example model map and `--include-partial`.

The comparison ran from 12:01:45.995 to 12:08:07.620 UTC: **381.63 seconds (6 min 22 s)**, including full-assembly indexing/annotation parsing, inference, identification and report generation. The full CLI invocation, including preflight checks, took **386.86 seconds (6 min 27 s)** according to GNU time and exited zero. Setup and the earlier discovery pass are excluded.

| Server result | Observed |
|---|---:|
| Selected windows | 7 |
| Predicted coding models | 24 (4 mammalian, 20 vertebrate) |
| Neighbours after target-overlap filtering | 17 |
| Tentatively identified neighbours | 15 |
| Unnamed neighbours retained | 2 |
| DNA window pairs | 11 |
| Pairs with insufficient unique synteny anchors | 10 |
| Pipeline errors | 0 |

The independent audit passed all 7 window roundtrips, 24 raw model roundtrips, 24 native-coordinate model roundtrips, 4 codon-alignment sequences, 39 report links, 6 SVGs, 17 source-module hashes and both CPU batches. See [server run summary](benchmarks/tiberius_host1_2026-11-05.json) and [server audit](benchmarks/tiberius_host1_2026-11-05.audit.json). Differences from the annotation-free fragment pilot reflect its different target models/window boundaries and reference inputs; the two counts are not a controlled accuracy comparison.

## Limits

The two shared models avoid per-species training/downloads; they do not establish equal accuracy across all vertebrate groups or in nonvertebrates. Boundaries, short windows, repeats and incomplete assemblies can affect predictions. Names are homology-based and homology groups are local, provisional groups. A ±50-kb window often has fewer than two unique anchors, so the report may correctly leave synteny unassessed. Neither sequence similarity nor gene order establishes HH1–3 expression or CRE function.

For standalone window prediction timings and RAM measurements, see the earlier [CPU benchmark](TIBERIUS_CPU_BENCHMARK.md). Its 18-window workload uses repeated windows from five species, not 18 species, and excludes the rest of the pipeline.
