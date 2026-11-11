# Full NANOG stage 2 — 25 saved loci, 13 species

All **25 saved candidate loci** reached STOP 2. Unsaved draft choices were not promoted to approvals. The full extraction/prediction run completed in **24 min 32 s** with exit code 0. Broad protein identification was then repeated in **75 s** after fixing symbol-case handling; genomic sequences, CDS, proteins and prediction models were unchanged.

Outputs contain **260 model alternatives: 158 annotation transcripts and 102 Tiberius predictions**, representing 81 candidate gene groups. Overlapping windows and alternative transcripts repeat physical genes; these are not 260 distinct genes. There are 259 structurally complete CDS models and one partial annotation model (missing start). Structural completeness does not prove function or expression.

Tiberius predicts a coding model over the selected seed in **23/25 loci**. The second chicken and second quail candidates have annotated alternatives but no ab initio model over the seed. Those annotations remain available for review.

The five annotation-free peafowl windows were processed. Two focal predictions have NANOG-family support. The other three have broader-database support for **NKX2-3**, **HOXD3 / hoxd3a**, and **HMX3**. They should be inspected as non-NANOG candidates; seed similarity alone did not establish NANOG identity.

Final naming: **195 homology assignments, 16 ambiguous assignments, 49 models without a passing database match**. The catalogue contains 121,120 reviewed vertebrate UniProt sequences including available isoforms (release 2026_04), not every vertebrate protein. Original symbols/accessions/taxa are retained. Case-only conventions such as AICDA/Aicda no longer create false competing-gene assignments; this resolved 54 model-level false ambiguity flags. Paralog suffixes remain distinct.

Neighbor evidence contains **88 unordered candidate neighboring-gene correspondences** (161 directional records), including 73 unique reciprocal proposals. These are local protein-supported proposals, not confirmed genome-wide orthologs. Target-like matches are separate and are not independent neighboring-gene anchors.

## Per-locus result

Coordinates are zero-based, half-open. A/T = annotation/Tiberius model counts. Seed A/T counts only overlapping models. Query coverage is the fraction of the best-scoring supplied query covered by an overlapping predicted protein; every query remains available in the HTML evidence details.

| Species | Candidate locus | Window A/T | Seed A/T | Predicted reference label(s) | Query covered | Review note |
|---|---|---:|---:|---|---:|---|
| Gallus gallus | `locus_4563bf8945a9762e` | 6/3 | 1/1 | NANOG | 57.5% | Pending manual review |
| Gallus gallus | `locus_f6c19f86fb9951eb` | 6/3 | 2/0 | — | — | No ab initio model over seed; annotated alternatives retained |
| Homo sapiens | `locus_0b4eb3c52dbfe7c8` | 14/1 | 2/1 | NANOG / NANOGP8 | 100.0% | Pending manual review |
| Alligator mississippiensis | `locus_7976ae58fc7c759f` | 3/3 | 1/1 | NANOG | 57.5% | Pending manual review |
| Alligator mississippiensis | `locus_7406913ef1011698` | 3/3 | 1/1 | NANOG | 54.1% | Pending manual review |
| Anas platyrhynchos | `locus_66f47eec132b42cb` | 8/4 | 2/1 | NANOG | 68.8% | Pending manual review |
| Anas platyrhynchos | `locus_3057e38726d20c82` | 7/4 | 1/1 | NANOG | 85.2% | Pending manual review |
| Anolis carolinensis | `locus_79510857289034d3` | 8/3 | 1/1 | NANOG | 56.3% | Pending manual review |
| Chelonia mydas | `locus_01217fa462e7dfb5` | 4/5 | 1/1 | NANOG | 58.4% | Pending manual review |
| Chelonia mydas | `locus_72556e5873ad2c9f` | 4/4 | 1/1 | NANOG | 54.1% | Pending manual review |
| Coturnix japonica | `locus_a7acf33399eed00e` | 13/4 | 1/1 | NANOG | 68.8% | Pending manual review |
| Coturnix japonica | `locus_966b22c7b00b3570` | 13/4 | 1/0 | — | — | No ab initio model over seed; annotated alternatives retained |
| Dromaius novaehollandiae | `locus_871d75a5b027d8be` | 5/4 | 1/1 | NANOG | 68.8% | Pending manual review |
| Dromaius novaehollandiae | `locus_fb29486e933484fa` | 5/3 | 1/1 | NANOG | 54.6% | Pending manual review |
| Lepisosteus oculatus | `locus_ae01e82ea7dd3e1d` | 18/5 | 3/1 | — | 66.4% | Pending manual review |
| Lepisosteus oculatus | `locus_2fb92248f35a489a` | 26/6 | 2/1 | — | 69.8% | Pending manual review |
| Meleagris gallopavo | `locus_4e37f0f96e65bc29` | 4/3 | 1/1 | NANOG | 68.8% | One partial annotation model in window |
| Monodelphis domestica | `locus_0188c55e49ed8541` | 3/3 | 1/1 | NANOG / NANOGP8 | 80.0% | Pending manual review |
| Taeniopygia guttata | `locus_885736efab882bd7` | 4/5 | 1/1 | NANOG | 57.0% | Pending manual review |
| Taeniopygia guttata | `locus_06cfe1513747fb0e` | 4/4 | 1/1 | NANOG | 54.1% | Pending manual review |
| Pavo cristatus | `locus_09b5824e2f013565` | 0/4 | 0/1 | NANOG | 57.5% | Pending manual review |
| Pavo cristatus | `locus_1ac19b6ff832833f` | 0/3 | 0/1 | NANOG | 85.2% | Pending manual review |
| Pavo cristatus | `locus_b5532e29f2aa44d7` | 0/10 | 0/1 | NKX2-3 | 65.9% | Non-NANOG homeobox-family candidate; inspect/reject |
| Pavo cristatus | `locus_bc0a0ae2f253daca` | 0/7 | 0/1 | HOXD3 / Hoxd3 / hoxd3 / hoxd3a | 44.3% | Non-NANOG homeobox-family candidate; inspect/reject |
| Pavo cristatus | `locus_23830a361185651e` | 0/4 | 0/1 | HMX3 | 38.4% | Non-NANOG homeobox-family candidate; inspect/reject |

## Verification and timings

- **117 tests run: 116 passed, one optional LASTZ test skipped.** Actual CPU Tiberius integration tests ran with shared cached weights.
- Bounded stress: 8,192 interval extractions from a 16 MiB gzip assembly; 76,800 isoform alignments collapsed to 48 competing gene-pair proposals; 256 randomized spliced CDS models on both strands; malformed BLAST statistics; corrupted databases; non-boolean draft states; predictor failure and truncated-gzip handoff.
- All eight independent audit blocks passed. BioPython reread the original native genome intervals and independently reconstructed/translated all 260 CDS models. GFF3 coordinates, all query support, raw BLAST statistics, near-best assignments and neighbor competition were checked. All 33 previous pilot models reproduced exactly.
- The first audit incorrectly demanded one particular representative among equal-score isoforms. It was repaired to verify the actual exported model pair and maximum group score; the corrected core audit was rerun. Its independent native roundtrips were then reused only after exact model/FASTA/GFF/input/prediction identity checks for the naming-only replay. The earlier audit log is retained.
- The real 25-locus report passed export/reload and wrong-checkpoint tests in an isolated DOM control. A dense synthetic UI control checked 45 gene groups per region, the hard ten-neighbor-per-side limit, both scales, missing-correspondence gaps and an empty failed-region handoff. No test decisions were written into scientific evidence.

Genome hashing / annotation parsing / interval extraction / annotation CDS: **395.3 s**. Tiberius: **997.8 s**. Original broad BLASTP: **72.0 s**. Final naming replay: **75.0 s**. Independent core audit: **105.8 s**. Final audit: **3.5 s**. Timings overlap with parallel regression work and must not all be summed as elapsed wall time.

All host2 computation ran through the resource launcher under the **aggregate 250 GiB RAM hard cap**, with **swap disabled**. Shared and job cgroups recorded no OOM events. The largest production child RSS was 1179472 KiB; it is not aggregate peak memory. The shared slice also contains a large reclaimable filesystem cache from the ongoing assembly download. No capacity claim is made from the bounded software stress controls.

For the 12 protein-discovery inputs, stage 1 had not recorded a genomic FASTA hash. Stage 2 records fresh genome fingerprints and checks the matching annotation hash and contig lengths. Peafowl additionally matched its exact original stage-1 genome fingerprint. The audit preserves that distinction instead of claiming a stronger historical identity check.

## Files and handoff

Open `report.html` (also served at http://127.0.0.1:8781/report.html). `windows.fna`, `cds.fna`, `proteins.faa`, `models.tsv`, `regions.tsv`, `models.gff3` and native per-genome GFF3 exports are beside it. `evidence.json` is the scientific checkpoint; `diagnostics/benchmark.json`, `diagnostics/audits/*verified.json` and the regression logs contain the checks and exact receipts.

All STOP-2 biological choices remain pending. No genomic alignment, codon alignment, protein MSA, final orthology call or expression analysis was performed. The five unsaved lamprey checkboxes remain outside this saved 25-locus run.
