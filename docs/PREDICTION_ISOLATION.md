# Independent window prediction — v0.4.2

Previously, all windows sharing Tiberius weights were supplied to one worker. The pinned upstream implementation clamps its requested inference length to the longest input sequence when every sequence is shorter than that requested length. Consequently, adding a longer unrelated window can change the inference length and neural/HMM chunking of an existing locus.

This was found while checking the five-genome NANOG/Pou5f3 pilot. Successful extraction, translation and predictor exit status had not established adequate exon recovery. An independent CDS-overlap audit found agreement for 25 of 35 fully contained annotated coding genes at a threshold of 80% reference CDS recall and 80% prediction CDS precision. The previous NANOG predictions recovered about 79% of annotated CDS in opossum and anole. These are comparisons to supplied annotation, not ground-truth sensitivity or specificity.

Since v0.4.2, each selected window runs in its own worker. Both shared weight sets are still cached once. Prediction failures are isolated to that window. CPU provenance records input-window count and requested/effective inference length; the upstream short-sequence clamping remains enabled independently for each window. Worker startup costs more than a shared multi-window process.

All **59 tests passed without skips** in 395.279 seconds with the real sequence tools and CPU Tiberius runtime. See the [test receipt](benchmarks/prediction_isolation_v042_2026-11-05.json).

The real Tiberius integration regression predicts the same PTGS2 genomic interval alone and after a longer unrelated window is added, then compares the returned coding structures. Every worker must report exactly one input window. This checks input isolation; it does not establish biological accuracy for arbitrary genes.

Query reports now show the predictor's own model counts, including target-overlapping predictions, and links to raw and native-coordinate GFF3. The red discovery target can come from annotation or protein-supported reconstruction and must not be mistaken for a Tiberius model.

Independent report checks are available as scripts:

- `docs/benchmarks/audit_window_predictions.py COMPARISON --panel PANEL.tsv` compares predicted and supplied-annotation CDS, draws separate exon tracks, and retains every requested genome in its overview. The panel requires `genome_id`, `species`, `accession` and `assembly_level`.
- `docs/benchmarks/audit_protein_identity.py --comparison COMPARISON --reference SWISSPROT.fasta.gz --work AUDIT_CACHE` compares selected/alternative target proteins and every window prediction to the full reference database using BLASTP. Top local matches are supporting evidence, not orthology assignments.

Neither script imports the `genecontext` package. Gene-family identity, individual-copy orthology, missing models and clipped genes require review beyond these checks.
