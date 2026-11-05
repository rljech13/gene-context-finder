# Full 18-assembly validation of NANOG and POU5 queries

This is an operational and biological-evidence audit of v0.4.2, not a confirmed-ortholog conservation dataset. The original two chicken proteins were searched across every contig of all 18 pinned assemblies on host1. Automatic target ranking exposes several biologically unsuitable selections. Successful program execution does not make those selections NANOG or POU5F3 orthologs.

The executed source is commit `3853f9851f1386ccedace253d1f5f91f3ef00617`. Its 18 Python modules were checked against the published checkout. The source inputs contain 98,388,872,280 bp according to the pinned assembly inventory. Sixteen assemblies have matching supplied GFF3 files; lamprey and peacock were processed without GFF3. The supplied panel itself contains one scaffold-level assembly, lamprey `GCA_018977245.1`; the other 17 are Chromosome or Complete Genome. This exception must remain visible when describing the dataset.

## What was actually executed

```bash
python -m genecontext analyze \
  --query inputs/queries/targets.protein.faa --query-type protein \
  --genomes inputs/genomes.tsv --out results/vertebrates18 \
  --prediction-models examples/vertebrates18.models.tsv \
  --predictor-cache cache/tiberius --reference-cache cache/references \
  --threads 16 --include-partial --offline \
  --reference-genome Gallus_gallus_GCF_016699485.2
```

The default protein discovery route uses TBLASTn seeds and local protein-supported reconstruction, alongside supplied annotation. Each selected target gets its own ±50 kb context window and a separate CPU Tiberius worker. Both input-annotation and DNA-only prediction evidence are retained. Two cached weight sets are shared: `vertebrates` (75,580,561 bytes) and `mammalia` (62,030,997 bytes), not one downloaded model per species. Whole-genome search and window prediction are separate computational stages. GNU time measured a maximum resident set of 68,292,492 KiB (about 65.1 GiB) for the full 98.4-Gb input run; the first isolated Tiberius worker reported 5.36 GiB. The worker footprint must not be presented as the memory requirement of the complete 18-genome analysis.

This run uses default BLASTn for local genomic-window comparisons, across all species pairs. MAFFT aligns selected proteins; coding DNA is back-translated from the protein alignment without inventing reference CDS for the protein-only queries. Candidate/provisional neighborhood similarity does not establish orthology or conserved regulatory activity.

## Independent checks

The audit scripts do not import `genecontext`. The roundtrip auditor initially treated resolvable ambiguous codons differently from source exports: BioPython resolves CCN as Proline, whereas the pipeline conservatively exports X for any non-ACGT codon. The corrected audit checks that explicit export convention and counts all three affected codons/models. Native DNA, CDS and protein files were unchanged. The initial failure and successful rerun are recorded in the completion provenance.

The checks are implemented in [sequence roundtrips](benchmarks/audit_sequence_roundtrips.py), [window prediction agreement](benchmarks/audit_window_predictions.py), and [independent protein identity](benchmarks/audit_protein_identity.py):

1. Re-extract all exported source CDS from the exact native genome coordinates, check translated proteins, raw/native Tiberius coding models, window orientation/checksums, codon alignments, HTML links, SVG syntax, source fingerprints and CPU provenance.
2. Compare Tiberius CDS with supplied GFF3 on the same strand. Agreement requires at least 80% reference CDS recall and 80% prediction CDS precision. Every overlapping isoform is considered, and only fully contained coding genes enter the agreement fraction. This measures agreement with annotation, not ground-truth sensitivity or specificity.
3. BLASTP selected and alternative target candidates, plus every Tiberius model, against the full cached Swiss-Prot library. Best reference hits and their coverage are supporting evidence, not individual-copy orthology assignments.

See [independent window prediction](PREDICTION_ISOLATION.md) for the fixed batch-composition problem and the 59-test real-tool regression run. See [the earlier NANOG search diagnosis](NANOG_SEARCH_VALIDATION.md) for why a five-genome pilot did not establish complete vertebrate target recovery.

## Findings that prevent an unqualified conservation claim

- The human POU5F3 query selects **POU3F2**, another POU class, although POU5F1 candidates are also found. Its BLASTP score is 229 versus 224 for one complete POU5F1 isoform. Sorting by sequence score alone is not a family discriminator.
- NANOG queries select LBX-like candidates in shark and lungfish, and an engrailed-like candidate in lamprey. These are not confirmed NANOG loci. This does not establish that NANOG is absent from those species.
- The lamprey POU query selects a POU3-like candidate.
- Anole selects an annotated `pou5f1` candidate. A shark POU5 candidate's closest reference uses the name `pou5f1.1`; legacy reference nomenclature and actual paralogy require separate review. Neither reference name alone proves a POU5F3 ortholog or a wrong-family hit.
- The selected axolotl NANOG-like model is only 99 aa. Longer annotated NANOG-like candidates remain in the alternative-locus evidence. Local family support does not establish completeness or the correct copy.

- Tiberius itself recovers about 51% of the selected lungfish POU5 candidate's annotated CDS, and about 79% of the selected opossum/anole NANOG CDS. In axolotl, a target-overlapping Tiberius model is 296 aa while the discovery model used in the comparison is only 99 aa; that model/copy discrepancy also requires review.

The comparative alignments and maps retain the default selections so the failure is inspectable. They are diagnostic outputs. Before using them for biological conservation, target-family/copy validation and selection must precede inclusion in comparative analysis. A missing or rejected candidate must remain unresolved rather than being replaced by a different homeobox/POU family.

## Completed-run measurements

The single end-to-end command returned **0** after **109.2 minutes**, followed by three successful independent audits. Search produced **1,724 coding candidate models**, and all **36 context windows** received independent CPU Tiberius runs. The predictor exported **149 models**. Of **121 fully-contained reference gene observations across windows**, **99** met both 80% CDS-agreement thresholds. These are gene observations within windows, not a genome-wide accuracy estimate.

Reference-identity categories: **5** different reference gene review, **4** local family reference support, **25** extended family reference support, **2** POU5 paralog review. The POU5 reference-name flags include nomenclature requiring review; they are not all proven wrong-family calls.

| Species / exact assembly | NANOG query: closest reference; Tiberius models | POU5F3 query: closest reference; Tiberius models |
|---|---|---|
| Lethenteron camtschaticum · [GCA_018977245.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_018977245.1/) · Scaffold | eng1a; 3 models (different reference gene review) | pou3f2-b; 2 models (different reference gene review) |
| Scyliorhinus canicula · [GCF_902713615.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902713615.2/) · Chromosome | lbx2; 3 models (different reference gene review) | pou5f1.1; 3 models (POU5 paralog review) |
| Lepisosteus oculatus · [GCF_040954835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040954835.1/) · Chromosome | NANOG; 5 models (local family reference support) | POU5F3; 3 models (extended family reference support) |
| Danio rerio · [GCF_052040795.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_052040795.2/) · Complete Genome | NANOG; 4 models (local family reference support) | pou5f3; 3 models (extended family reference support) |
| Protopterus annectens · [GCF_019279795.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019279795.1/) · Chromosome | lbx1; 1 model (different reference gene review) | POU5F3; 2 models (extended family reference support) |
| Ambystoma mexicanum · [GCF_040938575.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040938575.1/) · Chromosome | NANOG; 1 model (local family reference support) | POU5F3; 1 model (extended family reference support) |
| Monodelphis domestica · [GCF_027887165.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_027887165.2/) · Chromosome | NANOG; 3 models (extended family reference support) | POU5F3; 4 models (local family reference support) |
| Homo sapiens · [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) · Chromosome | NANOG; 1 model (extended family reference support) | POU3F2; 1 model (different reference gene review) |
| Anolis carolinensis · [GCF_035594765.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_035594765.1/) · Chromosome | NANOG; 3 models (extended family reference support) | POU5F1; 3 models (POU5 paralog review) |
| Chelonia mydas · [GCF_015237465.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015237465.2/) · Chromosome | NANOG; 5 models (extended family reference support) | POU5F3; 3 models (extended family reference support) |
| Alligator mississippiensis · [GCF_030867095.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030867095.1/) · Chromosome | NANOG; 3 models (extended family reference support) | POU5F3; 3 models (extended family reference support) |
| Gallus gallus · [GCF_016699485.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/) · Chromosome | NANOG; 3 models (extended family reference support) | POU5F3; 11 models (extended family reference support) |
| Dromaius novaehollandiae · [GCF_036370855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036370855.1/) · Chromosome | NANOG; 4 models (extended family reference support) | POU5F3; 3 models (extended family reference support) |
| Anas platyrhynchos · [GCF_047663525.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_047663525.1/) · Complete Genome | NANOG; 4 models (extended family reference support) | POU5F3; 5 models (extended family reference support) |
| Coturnix japonica · [GCF_054131305.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_054131305.1/) · Chromosome | NANOG; 4 models (extended family reference support) | POU5F3; 14 models (extended family reference support) |
| Meleagris gallopavo · [GCF_905368555.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_905368555.1/) · Chromosome | NANOG; 3 models (extended family reference support) | POU5F3; 10 models (extended family reference support) |
| Pavo cristatus · [GCA_045791835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_045791835.1/) · Chromosome | NANOG; 4 models (extended family reference support) | POU5F3; 10 models (extended family reference support) |
| Taeniopygia guttata · [GCF_048771995.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048771995.1/) · Chromosome | NANOG; 5 models (extended family reference support) | POU5F3; 9 models (extended family reference support) |

See the [machine-readable run receipt](benchmarks/vertebrates18_v042_2026-11-05.json), [target evidence table](benchmarks/vertebrates18_v042_target_review.tsv), and [per-window predictor audit](benchmarks/vertebrates18_v042_window_prediction_audit.tsv).

The portable inspection report starts at `validation/target_review.html`; `validation/report.html` contains separate GFF3/Tiberius exon tracks for all 36 windows. The comparative report and coloured protein/CDS alignment viewers explicitly retain and flag automatic selections. The bundle contains raw/native prediction GFF3, CDS, proteins and CPU logs, with per-file checksums. Whole genome files and runtime/reference caches remain external.
