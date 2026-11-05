# Prediction integration fixture

`HS04636.fna` is the uppercase DNA sequence from the official AUGUSTUS 3.5.0 example `examples/HS04636.fa`:
https://github.com/Gaius-Augustus/Augustus/blob/v3.5.0/examples/HS04636.fa

This small human sequence is used to exercise actual multi-exon prediction, coordinate conversion on both native strands, annotation-independent predictions, and retention of unnamed genes. The test does not claim that AUGUSTUS predictions are a validated biological gene annotation.

`PTGS2.reference.faa` is the 604-aa human PTGS2 protein matching UniProtKB/Swiss-Prot [P35354](https://www.uniprot.org/uniprotkb/P35354/entry). Source attribution: UniProt Consortium, CC BY 4.0. It supplies an independent protein-identification reference for the Tiberius integration controls. The positive target control inserts the included chicken Pou5f3 CDS into a constructed sequence separated from the human example by Ns; it is explicitly a computational integration fixture, not a natural genomic arrangement.

## NANOG search regression (2026-11-05)

`Danio_nanog_locus.fna` is the exact original FASTA interval NC_141044.1:14521832–14535827 (1-based inclusive) from assembly GCF_052040795.2. `Danio_nanog_locus.gff3` is the matching RefSeq nanog annotation, shifted to that interval and renamed to contig `chr_test`; all source CDS phases are preserved. Source transcript: NM_001098392.1; GeneID: 792333. The supplied assembly/GFF3 from the full-genome pilot are the source of these coordinates; no alignment was used to invent the gene model. The expected 384-aa genomic translation has SHA-256 `3c4596ae8dc9a468015fef0d36971a2a4b4a9271b89e546c79763b67f7400d4d`.

`chicken_nanog.faa` is the 309-aa chicken NANOG protein NP_001139614.1 used in the original search, with an opaque test identifier. Tests also erase NANOG annotation names and reverse-complement a copy, to require sequence evidence and correct native coordinate conversion. Both fixtures derive from public NCBI RefSeq records. The baseline fast route misses this real locus; the translated-genome rescue must recover the annotated full model and an explicitly partial coding model without GFF3.
