# Window gene prediction (v0.4)

```text
query → find target in every genome → whole target ±50 kb
      → predict genes from the DNA of every window
      → extract predicted CDS and proteins
      → identify predicted proteins
      → compare neighbour families, order and orientation
```

Prediction and identification are separate. A protein without a confident name stays visible as `uncharacterized` and can still match unnamed predictions in other genomes. Window DNA comparisons and target CDS/protein comparisons remain independent of neighbour naming.

## Default: Tiberius on CPU

Prepare the isolated runtime, two shared models and naming database once:

```bash
gene-context setup-predictor --with-reference
gene-context doctor
```

The backend runs pinned Tiberius 2.0.9 on CPU, with GPUs hidden explicitly. Both strands are predicted without GFF3, protein or synteny hints. It uses unmasked DNA (`--no_softmasking`), requests a 400,050-base inference length, and uses neural batch size 1 to limit memory. Upstream Tiberius reduces the inference length for shorter input sequences; the effective value is recorded in CPU provenance. Original window FASTAs retain their original letter case for DNA comparison. No model training or RNA-seq is required.

Since v0.4.2, each window runs in a separate worker, sequentially. This prevents a longer unrelated window from changing another window's upstream length clamping and neural/HMM chunking. The two model weight sets remain shared cached files, not a separate download per species. Worker startup adds overhead. `--threads` limits the TensorFlow intra-operation pool; the inter-operation pool is 1. `--tiberius-batch-size` controls neural chunks within that one window and can be increased at the cost of memory. Provenance records input-window count, requested/effective inference length, visible devices, package versions, source-file hashes and weight checksums. Query reports show separate counts and GFF3 links for the predictor's own models; a discovery target from GFF3 is not itself a Tiberius prediction.

Model choices:

| Value | Use |
|---|---|
| `vertebrates` | Default shared vertebrate model |
| `mammalia` | Shared mammalian unmasked model |
| Local directory | Compatible custom unmasked ab initio model |

There is no model download per species. The two built-in archives total 66.2 MB; extracted files total 71.4 MB. The isolated Python runtime is much larger; allow about 2 GB on macOS or 3 GB on Linux for setup, excluding genomes and results. See [measured CPU benchmark](TIBERIUS_CPU_BENCHMARK.md) and [setup](SETUP_AND_DOWNLOAD.md).

Use an optional manifest column:

```tsv
genome_id	fasta	gff3	tiberius_model
bird	bird.fna.gz	bird.gff3.gz	vertebrates
mammal	mammal.fna.gz		mammalia
```

Omitted values use `--tiberius-model vertebrates`. No taxonomy is inferred. Alternatively, `--prediction-models models.tsv` overrides by `genome_id`, using columns `genome_id` and `model`; custom paths are relative to that table. This also works with saved searches. Extra rows in a model map can describe genomes outside the selected subset. The [18-species map](../examples/vertebrates18.models.tsv) matches the example downloader IDs.

The interface accepts arbitrary supplied genomes, but built-in model accuracy is not universal across all organisms. For nonvertebrates, select and validate a suitable custom model or an AUGUSTUS model. The Tiberius adapter rejects models requiring softmasking, ClaMSA or embedded HMM inputs rather than silently running them incorrectly.

The window includes the whole selected target envelope plus the requested flanks, clipped at contig ends. Boundary genes can be partial or missed. Windows shorter than 1,000 bp are explicitly marked `prediction_window_too_short`; no prediction there does not mean biological absence. A window has no fixed number of expected genes.

## Identification after prediction

BLASTP searches predicted neighbour proteins against Swiss-Prot by default, plus proteins extracted from provided neighbouring annotations and any `--context-proteins reference.faa`. `--context-reference auto` adds Swiss-Prot only when usable local references are absent; `--context-reference none` uses only local naming references. These switches do not disable gene prediction.

Swiss-Prot comes from the [official UniProt release directory](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/). Its download uses `RELEASE.metalink`, verifies the official MD5 and expected byte count, and records the release, URL and SHA-256. An existing cache is verified and reused without contacting UniProt; use a new `--reference-cache` to obtain a newer release. Add `--offline` to analysis to prohibit downloads and require local caches. Short (<30-aa) or nonstandard Swiss-Prot proteins are skipped with counts; supplied reference proteins are validated strictly. Attribution: UniProt Consortium, [CC BY 4.0](https://www.uniprot.org/help/license).

Naming requires at least 30% identity and 50% coverage of both predicted and reference proteins by default. Non-overlapping HSPs contribute once. Highest summed bitscore selects the match. At exactly tied bitscores, curated Swiss-Prot references take priority over supplied references, then annotation-derived references. Conflicting names among equally preferred top matches stay `uncharacterized` with `ambiguous_best_matches`; all tied IDs remain recorded. Names are tentative protein-homology assignments, not confirmed species-specific gene symbols or functions. Raw BLASTP evidence remains available.

Unusable or disrupted predictions stay displayed without a protein anchor. A prediction or identification error is reported separately from a no-hit/no-prediction result. If identification fails, already predicted genes, sequence exports and DNA/protein comparisons remain available, and the command exits with an error status.

## Independent annotation track and outputs

GFF3 is optional and never constrains window predictions. `provided_neighbors.tsv` preserves its independent neighbours. `provided_gff_overlap` records same-strand overlap with predicted neighbours; overlap alone does not validate exact exon structures.

| Output within `comparison/` | Meaning |
|---|---|
| `context_prediction/tiberius_*/` | Batched oriented-window FASTA, raw GFF3/CDS/protein exports and logs |
| `context_predictions/<genome>.gff3` | Predictions in native input-genome coordinates, with unique IDs and CDS phases |
| `context_predictions.tsv` | All coding predictions, integrity and identification evidence |
| `context_prediction_cds.fna`, `context_prediction_proteins.faa` | CDS/proteins extracted from predicted coding models |
| `neighbors.tsv` | Predicted neighbours, including unnamed genes, strands, clipping and homology groups |
| `context_identification/` | Exact naming library, predicted query proteins and BLASTP matches |
| `context_reference.json` | Naming source, database release/checksums and identification outcome |
| `context_prediction.json` | Predictor, per-genome models, CPU proof, runtime and weight fingerprints |

Predictions strongly overlapping the original target CDS are retained and flagged `target_overlap=true`, but excluded from neighbour anchors. Target DNA-only regions without a coding model have no target CDS overlap test; inspect overlapping predictions there. Target amino-acid/codon comparisons use the original sequence-supported target models; the pipeline does not replace them silently with ab initio predictions. Only references supporting naming assignments join the smaller synteny catalogue; the full Swiss-Prot library is used only for identification.

## Optional AUGUSTUS and legacy mode

Install AUGUSTUS separately, then use `--predictor augustus --augustus-species MODEL` or per-genome `augustus_species` entries. List installed models with `augustus --species=help`. This backend predicts both strands with partial models allowed, preserves lowercase softmasking, disables UTR prediction and alternative sampling, and supplies no external hints. Raw output is `context_prediction/<locus>/augustus.gff3`; model configuration hashes are saved when the configuration location is reported. See the [official prediction tutorial](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/tutorial/prediction.html).

`--context-mode annotation` preserves the v0.2 route: provided GFF3 neighbours, or BLASTX screening followed by miniprot reconstruction for missing annotations. This is a separate method from predicting genes from DNA before identifying them.

Ab initio models can split, merge or miss genes. These are exploratory local comparisons, not a complete validated genome annotation or proof of orthology, expression or regulatory activity.
