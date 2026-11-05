# Window gene prediction (v0.3)

The default workflow is:

```text
query → target search in every genome → target ±50 kb
      → predict genes from each window's DNA
      → extract predicted CDS and proteins
      → identify predicted proteins
      → compare neighbour families, order and orientation
```

Gene prediction and gene identification are separate operations. An unknown protein match does **not** delete a predicted gene. The window DNA and target CDS/protein comparisons are also retained independently of neighbour identification.

## Predictor

AUGUSTUS predicts splice-aware protein-coding gene structures from the DNA of **every** selected window, regardless of GFF3 availability. The same settings are used across windows, with an explicitly selected trained species model for each genome. Predictions use both strands, allow partial models at window boundaries, preserve lowercase softmasking and disable UTR prediction and alternative sampling. No annotation, reference protein or synteny hints are supplied to AUGUSTUS.

List installed models with `augustus --species=help`. Select a model using `--augustus-species MODEL`, or override it per genome in the input TSV:

```tsv
genome_id	fasta	gff3	augustus_species
bird	bird.fna.gz	bird.gff3.gz	chicken
human	human.fna.gz		human
```

These are predictor model identifiers, not a universal assignment rule for all birds or mammals. For a species without a trained model, choose and validate a suitable related model or supply a trained model through AUGUSTUS configuration. A gene-prediction model is required; the pipeline does not silently assume that one model fits all input organisms. AUGUSTUS targets eukaryotic gene structures. Input generality does not guarantee prediction accuracy at every evolutionary distance.

The window contains the whole selected target envelope plus 50 kb at each side, clipped at contig ends. A gene crossing that boundary can be partial or missed. Lowercase input is treated as softmasking; an entirely lowercase genome can suppress prediction. Increasing the window or checking input masking can help investigate such cases. The window does not force a fixed number of genes to be predicted.

## Identification after prediction

Predicted neighbour proteins are searched with BLASTP against proteins extracted from provided neighbouring annotations and any supplied `--context-proteins` FASTA. With `--context-reference auto` (default), a lack of usable local reference proteins triggers Swiss-Prot preparation. Use `--context-reference swissprot` to add that database even when local references exist, or `--context-reference none` to use only local inputs.

Swiss-Prot comes from the [official UniProt release directory](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/). The downloader reads `RELEASE.metalink`, verifies its MD5 and expected file size, and records the release, source URL and SHA-256. The cached release is pinned: later runs verify and reuse it without contacting UniProt. Choose a new `--reference-cache` directory to obtain a newer release. The default cache is `~/.cache/gene-context/references`. Short (<30-aa) or nonstandard reference proteins are skipped with counts in the receipt; supplied reference files are validated strictly. Swiss-Prot attribution: UniProt Consortium, [CC BY 4.0](https://www.uniprot.org/help/license).

Prepare the database on an internet-connected machine before an offline analysis:

```bash
gene-context prepare-context --reference-cache /data/gene-context-references
# Copy that entire cache to the analysis machine if necessary.
gene-context analyze --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/analysis \
  --reference-cache /data/gene-context-references --threads 8
```

Naming requires ≥30% identity and ≥50% coverage of both predicted and reference proteins by default. Non-overlapping HSPs contribute once. The highest summed bitscore selects the supporting match; tied top scores with different names are marked ambiguous and retain the label `uncharacterized`. Names from UniProt `GN=` headers or reference identifiers are tentative homology-based assignments, not confirmed species-specific gene symbols or functions. All raw BLASTP matches remain available.

A predicted gene lacking a passing match remains `uncharacterized`. Its usable protein can still form a local homology group with unnamed predictions in other genomes, so lack of a name alone does not prevent synteny comparisons. Unusable/disrupted predictions remain displayed without a protein anchor. Prediction or identification failure is reported separately from gene absence.

## Annotation comparison and outputs

Provided annotations are retained as a separate track in `provided_neighbors.tsv`. `provided_gff_overlap` records same-strand overlap with a predicted neighbour; overlap alone is not validation of its exact exon structure. Provided annotations can supply naming reference proteins but do not constrain predictions.

- `context_prediction/<locus>/augustus.gff3`: original predictor output in oriented-window coordinates, including parameters and the tool's configuration location.
- `context_predictions/<genome>.gff3`: predictions converted back to native genome coordinates, with unique IDs, CDS phases and stops. The files include target-overlapping predictions, flagged `target_overlap=true`; these are excluded from neighbouring-gene anchors.
- `context_predictions.tsv`: all coding predictions, integrity and identification evidence. Target-overlapping predictions are flagged and do not receive neighbour names.
- `context_prediction_cds.fna`, `context_prediction_proteins.faa`: sequence extraction from all coding predictions.
- `neighbors.tsv`: predicted neighbours, including unidentified genes, strands, clipping, protein support, homology groups and overlap with supplied annotation.
- `context_identification/`: the exact search library, predicted queries and raw BLASTP evidence.
- `context_reference.json`: reference provenance, database release/checksums and identification outcome.
- `context_prediction.json`: per-genome predictor models and SHA-256 of model configuration files when the configuration directory is reported by AUGUSTUS.

The complete search library is used only for naming. Only references supporting assignments enter the much smaller local synteny catalogue. Genomic DNA comparisons still compare the actual windows; target amino-acid/CDS comparisons use the original target search models. Ab initio predictions can split, merge or miss genes, so these reports are exploratory local comparisons rather than a complete validated genome annotation.

## Previous behaviour

`--context-mode annotation` preserves the v0.2 route: use supplied GFF3 neighbours, and infer missing annotations by BLASTX screening followed by miniprot protein-supported reconstruction. That is a different method from the new default of predicting genes from DNA before identifying them.

Method: [official AUGUSTUS prediction tutorial](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/tutorial/prediction.html).
