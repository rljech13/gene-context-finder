# Chicken Pou5f3 query sequences

These versioned RefSeq sequences were retrieved and checked for the 2026-11-05 pilot.

- Species: *Gallus gallus*; gene: [Pou5f3, Gene 427781](https://www.ncbi.nlm.nih.gov/gene/427781).
- [chicken_Pou5f3.protein.faa](chicken_Pou5f3.protein.faa): [NP_001296301.2](https://www.ncbi.nlm.nih.gov/protein/NP_001296301.2), 389 amino acids. Use `--query-type protein`.
- [chicken_Pou5f3.cds.fna](chicken_Pou5f3.cds.fna): spliced CDS from [NM_001309372.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_001309372.2), 1,170 nucleotides including the stop codon. Use `--query-type cds`.

The CDS translates exactly to this protein plus its terminal stop. These are the same sequence contents used in the full-genome pilot, separated into single-query files. Sequence SHA-256 values (uppercase letters, no whitespace or header):

```text
protein  b0f214cdb0dfd533de7c0d8106c6bfdf166144ee0ecce49684037953d449c30a
CDS      d97d75123dd9f7900018eb98963539d74539f2b0f0149a730892b0e672530275
```

From the repository root, with your own local genome table:

```bash
gene-context run --query examples/queries/chicken_Pou5f3.protein.faa \
  --query-type protein --genomes genomes.tsv --out results/chicken_Pou5f3 \
  --threads 8 --skip-alignments
```
