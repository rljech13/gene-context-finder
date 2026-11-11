# Stage 2 pilot — 9 November 2026

The first connected STOP-2 run used three **saved** stage-1 selections: two overlapping chicken windows and one human window. Unsaved selections were not promoted to approvals. This was a two-species pilot, not an 18-genome run.

| Region | Window bp | Annotated transcript models | Tiberius models | Annotated gene groups without a prediction |
|---|---:|---:|---:|---:|
| Chicken `locus_4563bf8945a9762e` | 104,584 | 6 | 3 | 1 of 4 |
| Chicken `locus_f6c19f86fb9951eb` | 110,302 | 6 | 3 | 1 of 4 |
| Human `locus_0b4eb3c52dbfe7c8` | 109,745 | 14 | 1 | 3 of 4 |

Window lengths above include all selected DNA. The pinned Tiberius worker independently clamps its effective neural sequence length; exact execution settings and raw GFF3 remain in the result evidence.

The run completed in **217.12 seconds**, including real CPU Tiberius inference and BLASTP. It produced 33 model alternatives (26 annotation + 7 prediction), representing eight candidate gene groups across the two assemblies. Overlapping chicken windows repeat the same physical genes. Tiberius reconstructed the annotated NANOG CDS in chicken and human, but did not recover all annotated neighboring genes. Annotation and predictions are therefore retained independently; a missing prediction does not establish gene absence.

BLASTP identification used **121,120 reviewed vertebrate UniProtKB protein sequences including available isoforms**, release `2026_04`. Database preparation took 40.15 seconds. This catalogue is a curated subset, not the complete vertebrate proteome. The CLI also accepts a user-supplied broad protein FASTA.

There were 22 homology assignments, seven ambiguous assignments, and four models with no passing database match. These are counts of model alternatives, not distinct genes. Ambiguous names remain unresolved. The local protein comparison proposed **one unordered neighboring-gene pair** (two directional records). It is a family-level candidate correspondence for manual inspection, not confirmed orthology or a final synteny conclusion. Database naming can disagree with annotation naming; both remain visible in the report.

Host2 execution used eight CPU threads and the shared **250 GiB aggregate hard RAM cap**, with swap disabled. The largest child RSS was 1,449,760 KiB (about 1.38 GiB); that value is **not** the aggregate cgroup memory peak. Original evidence and pending biological decisions were preserved. No DNA alignment, codon alignment, protein MSA, or final conservation analysis was run.

The exact pilot command, source archive hash, input provenance, counts and timings are recorded in [the machine-readable benchmark](benchmarks/stage2_nanog_v090_2026-11-09.json). The pilot used an immutable source snapshot. Subsequent report/provenance guards were tested locally and did not replace the recorded prediction output.

See [STAGE2.md](STAGE2.md) for input preparation, the commands, and the second manual-review handoff.
