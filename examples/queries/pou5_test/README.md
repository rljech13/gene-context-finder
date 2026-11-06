# POU5 test reference panel

This is a test panel for multi-reference locus discovery. It contains **16 RefSeq transcript/CDS variants, 15 distinct protein sequences, four species and six named source genes**. The source gene set was selected explicitly; this is not a claim that every isoform has pluripotency activity or early embryonic expression.

| Species | NCBI gene | Included transcript variants |
|---|---|---:|
| Homo sapiens | [POU5F1, 5460](https://www.ncbi.nlm.nih.gov/gene/5460) | 5 |
| Danio rerio | [pou5f3, 30333](https://www.ncbi.nlm.nih.gov/gene/30333) | 7 |
| Gallus gallus | [Pou5f3, 427781](https://www.ncbi.nlm.nih.gov/gene/427781) | 1 |
| Xenopus tropicalis | [pou5f3, 100498076](https://www.ncbi.nlm.nih.gov/gene/100498076), [pou5f3.2, 100170181](https://www.ncbi.nlm.nih.gov/gene/100170181), [pou5f3.3, 100170593](https://www.ncbi.nlm.nih.gov/gene/100170593) | 3, one per gene |

The three frog entries are **different genes**, not splice isoforms. Human POU5F1 and nonmammalian POU5F3 are labelled separately. `POU5` is a search-panel label, not a one-to-one orthology assignment. Xenopus tropicalis was chosen as the diploid frog reference for this pilot.

## Files

- `proteins.faa`: annotated translations, one record per transcript variant; recommended query for the first test.
- `cds.fna`: matching spliced CDS, including annotated terminal stops where present.
- `queries.tsv`: exact FASTA IDs, species, official source symbols, transcript/protein accession versions, RefSeq prediction status, sequence hashes and direct NCBI links.
- `reference_loci.tsv`: NCBI Gene primary-assembly gene intervals; zero-based, half-open. These gene intervals can include UTRs and are not identical to CDS bounds.
- `provenance.json` and `sources/`: retrieval URLs, timestamps, hashes and original NCBI XML/FASTA responses.
- `controls.genomes.tsv`, `controls/`: **six small reference windows from four species**, each source gene with 10 kb padding. These are smoke-test inputs, not complete genomes. Frog windows overlap and intentionally contain neighbouring POU5 copies.

All NM_/XM_ transcript products listed in the fetched Gene records were considered, across their reported assemblies. This is a snapshot of **RefSeq**, not an exhaustive union of Ensembl, historical annotations and experimental transcripts. Predicted XM_/XP_ products remain explicitly labelled. No listed variant was excluded in this snapshot.

Each CDS was extracted from its annotated transcript feature and checked against its annotated protein translation. An annotated alternative initiation residue is permitted only at the first amino acid and is recorded in `initiation_exception`.

To refresh into a new directory, from the repository root:

```bash
python scripts/build_pou5_test_query.py --out data/pou5_reference_refresh --controls
```

The builder retains cached responses in its output directory. Use a new directory to retrieve a fresh snapshot; do not overwrite the shipped panel to reproduce this benchmark.
