# NANOG missed-hit diagnosis and sensitive discovery — v0.4.1

Date: 2026-11-05. This corrects false negatives in the previous target-discovery stage. The v0.4.0 checks established that sequence extraction and reports operated; they did not establish adequate recovery of these known NANOG loci.

## Why the fast search missed NANOG

The previous protein route used whole-genome miniprot mapping with a 50% query-coverage gate. Annotation proteins were checked only near its reported mappings. A gene missed by that mapper could therefore remain undetected even when the input GFF3 already annotated it.

Direct inspection found NANOG in the exact input annotations of all three missed species:

- *Danio rerio*: GeneID 792333, NM_001098392.1, NC_141044.1:14526662–14531384.
- *Lepisosteus oculatus*: GeneID 107075953, XM_069188084.1 / XM_069188085.1, NC_090699.1:14289174–14294722.
- *Monodelphis domestica*: GeneID 100013584, XM_007503619.3, NC_077234.2:125046152–125051077.

Coordinates above are 1-based inclusive gene envelopes from the supplied assembly-matched GFF3, not inferred from a gene-name web search. Exact assemblies are in the [original full-genome pilot](FULL_GENOME_PILOT.md).

The failure had two components. The fish proteins have strong local homology that covers less than half of the chicken query. In addition, lowering miniprot's output coverage to 10% still returned no mapping for either fish; the opossum's reported short mapping was at a different locus. Even mapping just the known NANOG intervals with the standard settings failed. Dense local seeding/chaining recovered alignments (`-k3 -M0 -n1 -l3`), so loosening only the final coverage filter was insufficient. This combined test does not isolate one individual miniprot heuristic as the sole cause.

Miniprot's default index uses sampled 6-mers in a reduced amino-acid alphabet. TBLASTn searches translated DNA with scored protein words (default word size 3) and can detect a conserved exon/domain independently of a spliced gene alignment. These word sizes are not directly interchangeable sensitivity scores. TBLASTn does not itself produce a validated exon/CDS model. Sources: [miniprot source/manual](https://github.com/lh3/miniprot/tree/81f9b93481cec16622eff339fc998e4f38da344e), [official TBLASTn options](https://www.ncbi.nlm.nih.gov/books/NBK279684/), [TBLASTN algorithm paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1779365/).

For the chicken query versus annotated zebrafish NANOG, BLASTP aligned **68 of 309 query residues (22.0% query coverage)** with **57.35% identity**. Whole-genome TBLASTn ranked the correct zebrafish locus first (E = 5.21e-12; 71.6 bits). That separate TBLASTn run took 1.95 s with 16 threads against an already built nucleotide database; database construction and annotation processing are excluded. This is not a general end-to-end runtime claim.

## Change in the pipeline

The default protein/CDS discovery is now:

```text
query protein → TBLASTn across every contig → local genomic seeds
              → matching GFF3 CDS or dense local miniprot reconstruction
              → independent BLASTP verification and coding-integrity checks
              → whole selected target ±50 kb → DNA-only context prediction
              → neighbour identification, synteny, DNA and protein comparisons
```

Seed HSPs require E <=1e-5, at least 30 aligned amino acids and at least 30% identity. An individual exon seed is not required to cover a fixed fraction of the full query. The 15% query-coverage gate (`--candidate-min-coverage`) applies to the subsequent BLASTP check of a coding candidate. `--seed-min-aa`, `--min-identity` and `--protein-evalue` control the other gates.

Without a usable matching annotation, local miniprot refinement uses the source DNA, both relevant coding structure and native coordinates, with dense settings. A model must overlap its original translated seed on the same strand and pass the protein check. At most 20 unannotated regions per query/genome are refined by default; all passing seeds are retained, including unresolved/unrefined regions. The report distinguishes a protein-supported DNA region, a partial coding model and a complete coding model.

New searches rank coding candidates by non-overlapping BLASTP HSP bitscore normalized to the query's self score. This replaces ranking primarily by the span of a possibly forced protein-to-genome alignment. `local_alignment` flags a limited domain/fragment match separately from CDS completeness. A structurally complete gene is not thereby a confirmed orthologue. Other homeobox/POU proteins, isoforms and alternative predictions remain in the candidate tables; they are not all renamed NANOG/Pou5f3.

`--protein-search miniprot` retains the fast whole-genome mapping option. The default TBLASTn route reads only the supplied genome and optional GFF3; it needs no downloaded target-reference library. Naming predicted neighbours still uses the separately configured protein reference cache.

## Regression controls

The actual zebrafish genomic NANOG interval and matched annotation are included as [test fixtures](../tests/fixtures/README.md). The chicken query has an opaque identifier, annotation NANOG names are erased during the test, and both native strands are exercised. The fast route misses the locus; the default route must recover the exact annotated 384-aa protein. Without GFF3 it must retain an explicitly partial sequence-supported reconstruction. A negative query must stay negative.

A separate 600-aa synthetic gene with twelve short coding exons checks that significant exon seeds shorter than 15% of the query are not discarded. Unit controls also check overlap-safe HSP scoring and ranking stronger local protein evidence above a forced long alignment. These are defined recovery controls, not a genome-wide specificity/sensitivity estimate.

All **58 local tests passed without skips** in 375.961 seconds, including the real sequence tools and CPU Tiberius integration controls. The wheel was built and checked against the source modules; `doctor` also verified the installed TBLASTn executable.

Run these real-tool controls with:

```bash
python -m unittest discover -s tests -p test_protein_discovery.py -v
```

## Full-assembly recheck

The same two original chicken proteins were searched afresh against all five original whole-assembly FASTAs on host1, using their matching GFF3 annotations. No extra fish/mammalian query or gene-name lookup was supplied to discovery. The top NANOG candidate in every species matched its input annotation. Full search, coding exports and the downstream comparison are rechecked for this release; measured results and independent audits accompany the report.

The observed top annotated NANOG candidates in the final whole-genome search were:

| Species | Previous miniprot-only route | BLASTP query coverage | BLASTP identity | Bitscore |
|---|---|---:|---:|---:|
| *Lepisosteus oculatus* | Missed | 48.22% | 34.81% | 84.7 |
| *Danio rerio* | Missed | 22.01% | 57.35% | 91.3 |
| *Monodelphis domestica* | Missed | 75.08% | 39.11% | 140.0 |
| *Anolis carolinensis* | Recovered | 82.20% | 50.18% | 236.0 |
| *Chelonia mydas* | Recovered | 100.00% | 58.26% | 362.0 |

Coverage is the union of query-aligned positions; BLASTP identity uses aligned HSP columns, including gaps. The gar's second annotated isoform is retained separately. All five top models have structurally complete CDS; the two fish matches remain explicitly labelled local homology. These percentages describe the query-to-candidate protein alignments, not overall cross-species genomic conservation.

The final complete host1 run exited zero: **517 candidate transcript models, 10 selected windows, 35 context predictions, 25 predicted neighbours and 20 DNA window pairs**, with no search or comparison errors. All ten selected target coding models are structurally complete. The broad candidate inventory includes other homeobox/POU families; 517 is not a count of NANOG/Pou5f3 copies. Nineteen synteny pairs had fewer than two unique anchors; one had assessable but different local order/orientation.

Discovery took 403.64 seconds and comparison took 431.00 seconds according to their stage manifests. The full command took **838.37 seconds (13 min 58 s)** including preflight checks, using 16 threads on the AMD EPYC 7763 server. GNU time reported maximum process RSS 7,947,572 KiB (about 7.58 GiB); this is not the simultaneous aggregate of parent/child memory. Setup/downloads are excluded. Full-assembly annotation processing and indexes are part of this measurement; the small isolated TBLASTn timing above must not be substituted for it.

The independent BioPython audit passed **all 517 source CDS/protein roundtrips against native whole-genome DNA, 10 window roundtrips, 35 raw and 35 native-coordinate context-model roundtrips, 10 codon-alignment sequences, 39 HTML links, 6 SVGs and all 18 exact source-module fingerprints**. Both prediction batches confirmed CPU-only execution. SVGs were rendered and inspected. The executed source matches the released package byte-for-byte.

Evidence: [full run summary](benchmarks/nanog_v041_2026-11-05.json), [independent audit result](benchmarks/nanog_v041_2026-11-05.audit.json), [independent audit script](benchmarks/audit_nanog_search.py). The script imports BioPython, not `genecontext`, and requires access to the original genome paths/indexes recorded by the run.

The increased candidate inventory is expected: a shared homeodomain can recover other homeobox families. Gene-family identification and individual-copy orthology still require review of alignment scope, competing hits and genomic context. This correction establishes that the known NANOG loci are recovered; it does not make every local hit a NANOG orthologue or establish HH1–3 expression.
