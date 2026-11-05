# Window-prediction validation (v0.3.0)

The v0.3 test suite has 46 tests: the previous 40 regression controls plus three reference-library controls, one prediction unit control and two real-tool prediction/identification controls. External integration controls explicitly skip when the required tools are absent; a green dependency-free CI run alone does not demonstrate that AUGUSTUS was executed.

All **46 tests passed locally**, including the actual external-tool controls, in 138.019 seconds. The three reference-library tests were rerun successfully after adding the incomplete-cache-receipt assertion. A separate audit importing no `genecontext` modules checked the Swiss-Prot control's two native windows, two CDS extractions, two protein translations, all 15 source-module fingerprints, 18 HTML links and both SVG files; all checks passed.

Local validation used AUGUSTUS 3.5.0, miniprot 0.18-r281, NCBI BLAST+ 2.17.0+, samtools, MAFFT 7.526, and optional LASTZ 1.04.60. The new controls cover:

- The same DNA produces the same context predictions with a supplied GFF3 track and with no GFF3. A fictitious annotated gene does not force AUGUSTUS to predict it.
- Actual multi-exon predictions are mapped to both native genome strands. Re-extracted CDS and translated proteins round-trip to the exported sequences.
- Unnamed predicted genes are retained and can form shared protein homology anchors without a naming database.
- Post-prediction BLASTP distinguishes a missing match from tied matches with conflicting names.
- A local HTTP server exercises real library transfer, official-checksum validation, offline pinned-cache reuse, preservation of corrupt files, and rejection of incomplete cache receipts.
- Reference headers, skipped proteins and removal of unused search-library sequences from local synteny comparisons are checked.

## Actual Swiss-Prot control

The public AUGUSTUS human example [HS04636](https://github.com/Gaius-Augustus/Augustus/blob/v3.5.0/examples/HS04636.fa), 9,453 bp, was analysed in its original and reverse-complement orientations. These are **two representations of one biological sequence**, not two species. Both input GFF3 fields were empty. The query was a 400-bp noncoding seed from the start of that sequence, so the gene in the context was independently predicted rather than supplied as a target protein.

AUGUSTUS with the `human` model recovered a 10-CDS-exon model on each representation. Subsequent BLASTP against checksum-verified Swiss-Prot release **2026_03** identified both proteins as **PTGS2**, matching [P35354](https://www.uniprot.org/uniprotkb/P35354/entry) with 100% identity and complete coverage of the 604-aa reference. The native strands were opposite, while the oriented context and protein sequences agreed.

The identification library contained 569,692 usable proteins; 6,056 short or nonstandard records were skipped and counted. Only the one supporting Swiss-Prot reference and two predicted neighbour proteins entered the local homology catalogue. Both contexts shared one anchor; order comparison correctly remained `insufficient_unique_anchors`, because at least two distinct anchors are required.

The actual downloaded FASTA was 93,801,562 bytes. Its official MD5 was `bc9d398533e6df582b563c6c03093bd0`, and the recorded SHA-256 was `a9c3496aa727e61d25d9e1291c9ba28d2c2ad19c0e26569dc083d53ec4169536`. The official release metadata, model configuration hashes, predictor output and raw BLASTP matches were retained in the run.

This is an execution, coordinate and naming control for the new order of operations. It does not estimate genome-wide sensitivity, validate a predictor model across all vertebrates, or replace the v0.2 five-assembly target-search pilot. No complete 18-genome panel was rerun for this release.
