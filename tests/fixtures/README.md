# Prediction integration fixture

`HS04636.fna` is the uppercase DNA sequence from the official AUGUSTUS 3.5.0 example `examples/HS04636.fa`:
https://github.com/Gaius-Augustus/Augustus/blob/v3.5.0/examples/HS04636.fa

This small human sequence is used to exercise actual multi-exon prediction, coordinate conversion on both native strands, annotation-independent predictions, and retention of unnamed genes. The test does not claim that AUGUSTUS predictions are a validated biological gene annotation.

`PTGS2.reference.faa` is the 604-aa human PTGS2 protein matching UniProtKB/Swiss-Prot [P35354](https://www.uniprot.org/uniprotkb/P35354/entry). Source attribution: UniProt Consortium, CC BY 4.0. It supplies an independent protein-identification reference for the Tiberius integration controls. The positive target control inserts the included chicken Pou5f3 CDS into a constructed sequence separated from the human example by Ns; it is explicitly a computational integration fixture, not a natural genomic arrangement.
