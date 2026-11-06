# A workflow with two manual review stops

The recommended v0.5 workflow has three execution stages. Each stage ends and writes its output to disk. Nothing waits at a terminal prompt, and no process automatically resumes after a review file is edited.

| Stage | Automatic work | Required human decision |
|---|---|---|
| 1. Search | Search every input genome; export candidates, scores, coding status, CDS/proteins and a review report. | Select loci to inspect, explicitly exclude a pair, or mark it unresolved. |
| 2. Contexts | Extract the selected windows, predict genes from DNA, identify proteins and draw provisional contexts/CDS structures. | Choose the target coding model and decide which neighbours enter comparison. |
| 3. Compare reviewed inputs | Reuse saved windows/predictions; compare reviewed neighbours, genomic DNA, proteins and codon-guided CDS. | Inspect the resulting alignments and evidence before scientific interpretation. |

A successful exit from stages 1/2 means that a checkpoint was written. Its manifest status is `awaiting_review`, not biological validation. An empty decision never approves a candidate.

## Stage 1 — search, then stop

```bash
gene-context search --query targets.faa --query-type protein \
  --genomes genomes.tsv --out results/01_search --threads 8
```

Open `results/01_search/review/report.html`. Edit `results/01_search/review/loci.decisions.tsv` with a text editor or a TSV-compatible spreadsheet editor, preserving tab separators and exact IDs. Keep one row for every target-group/genome pair.

| Field | Value |
|---|---|
| `review_id`, `group`, `genome_id` | Keep the exported identifiers. `review_id` binds the decisions to this checkpoint. |
| `decision` | `inspect`, `exclude`, or `unresolved`. |
| `model_ids` | For `inspect`, copy one or more exact IDs from `candidates.tsv`, separated by semicolons. Empty for other decisions. |
| `notes` | Required for `exclude`/`unresolved`; record why the pair is not advancing. |
| `suggested_model_id` | A ranking suggestion only. Copying it to `model_ids` is an explicit choice; it is never accepted automatically. |

Several distinct loci can advance for one pair. Pick one seed per locus; all matching source isoforms are offered again at STOP 2. Rankings, source gene names and short domain matches do not establish the correct gene family or copy. A pair without candidates still needs `exclude` or `unresolved`; the latter does not mean biological absence.

## Stage 2 — predict contexts, then stop

```bash
gene-context contexts --review results/01_search/review \
  --decisions results/01_search/review/loci.decisions.tsv \
  --out results/02_contexts --threads 8
```

The default window is the selected gene/model plus 50 kb on either side. Matching GFF3 supplies an independent reference; the default DNA-only Tiberius prediction runs for every selected window. The usual predictor/reference options remain available, including `--prediction-models`, `--tiberius-python`, `--offline`, and `--context-reference`.

Open `results/02_contexts/report.html`. It links:

- Source isoforms and predicted target-model options, lengths, completeness, strand, coordinates and actual flank lengths.
- Separate per-locus CDS tracks: blue source models and green predictions.
- Provisional context diagrams and neighbour evidence.
- Exact CDS/protein FASTA for coding options, plus raw/native prediction GFF3 and logs in the saved directories.

No target protein/CDS alignments or genomic DNA comparisons are run by this stage.

An optional [offline reference check](TARGET_REFERENCE_CHECK.md) can compare `model_options.proteins.faa` against your target/paralog panel before decisions. It runs separately and leaves the sealed checkpoint and all decisions unchanged.

Complete **two** files:

**`models.decisions.tsv`** has one row per original target-group/genome pair. Set `decision=keep` and copy exactly one `option_id`, or set `exclude`/`unresolved` with notes. Choosing a prediction replaces the discovery model's CDS/protein in the final comparison. A `dna_only` option retains DNA/context analysis without inventing a coding sequence. Opposite-strand and disrupted coding models cannot be promoted as the target protein; revisit the locus choice or use DNA-only analysis.

A partial coding model requires `allow_partial=yes` and notes. If a prepared window cannot supply the requested flanks, require `allow_short_flanks=yes` and notes. Selecting a new model preserves the prepared DNA window; shortage is shown rather than silently claiming a fresh ±50 kb extraction. Re-run the contexts stage with a larger `--window` when more DNA/prediction coverage is needed.

**`neighbors.decisions.tsv`** lists every predicted gene, including models previously overlapping the discovery target. For neighbours of the final kept targets, choose `keep`, `unnamed`, or `exclude`. A name override and an exclusion require notes. `unnamed` retains coordinates/protein homology without endorsing the borrowed gene name. Rows belonging to unchosen windows or overlapping the final target CDS may remain blank; the target-overlap role is recomputed from the chosen model. Gene-level renaming does not edit sequences or prove orthology.

If every pair is explicitly excluded/unresolved, the next command writes a terminal decision report with no comparable targets and runs no prediction/alignment tools.

## Stage 3 — compare the reviewed inputs

```bash
gene-context compare-reviewed --review results/02_contexts \
  --models results/02_contexts/models.decisions.tsv \
  --neighbors results/02_contexts/neighbors.decisions.tsv \
  --out results/03_comparison --threads 8
```

Open `results/03_comparison/report.html`. This command does not search genomes or execute Tiberius again. It checks the chosen CDS/protein against its exact window DNA, then runs the neighbour comparisons, local genomic alignments, MAFFT protein alignment and codon back-translation. A revised stage-2 decision can be compared in a **new output directory**, reusing the same stage-2 checkpoint. Changing a stage-1 locus decision requires a new contexts output; an existing checkpoint keeps the original decision snapshot.

The result includes:

- `approved_models.tsv`, `approved_neighbors.tsv`: exact decision snapshots.
- `reviewed_targets.tsv`: discovery model, final model, origin and actual flank lengths.
- `reviewed_targets.gff3`, `reviewed_target_cds.fna`, `reviewed_target_proteins.faa`: the actual selected coding structures/sequences. DNA-only targets use a GFF3 `region` feature with no invented gene/CDS.
- `summary.tsv`: every requested pair, including excluded and unresolved pairs.
- Context diagrams, synteny tables, DNA blocks, protein/codon alignments and their coverage/identity tables.
- `manifest.json`: checkpoint/decision fingerprints and newly executed commands. `review.prediction_reused=true` records the resume path.

Original prediction GFF3 remains raw evidence. The final prediction table recomputes target-overlap flags and retains the original flags separately. Reviewed target GFF3 is the explicit record of the selected model.

## Reuse existing results

Create STOP 1 from a saved search:

```bash
gene-context review-search --results previous/search --out results/01_review
```

Create STOP 2 from a completed comparative run:

```bash
gene-context review-contexts --results previous/comparison --out results/02_review
```

For a relocated portable report, pass `--search-results PATH_TO_SAVED_SEARCH`. Its manifest/models fingerprints must match the original comparison. The imported selections remain unapproved: every pair still needs a stage-2 decision. This import covers only the previously prepared windows. To inspect another locus, use STOP 1 and run contexts for that locus.

## Review on one computer, execute on another

Copy the **whole checkpoint directory**, including `seal.json`, rather than just its tables. Its review ID stays the same. Edit the decision files on the reviewing computer. When resuming elsewhere, `contexts` and `compare-reviewed` accept `--search-results PATH_TO_THE_SAME_SEARCH`. Every recorded source fingerprint must match; this changes the source location without regenerating the checkpoint or silently reusing decisions against different evidence.

Prediction of newly selected loci still requires the original genome files. Comparison from an existing context checkpoint uses its saved windows/predictions and needs the matching saved search exports, BLAST+ and MAFFT; it does not need the complete genomes or Tiberius runtime.

## Evidence and decision handling

Each decision row contains the unique `review_id` of its checkpoint. Decisions from a different or regenerated checkpoint are rejected even if model IDs happen to match. Evidence files are sealed by SHA-256. Source manifest, candidate tables and source sequence exports are also checked. Editing a decision TSV is expected; editing a model, window, state file or other sealed evidence invalidates the checkpoint. Create a fresh stage if the evidence needs to change. Unknown IDs, duplicate/missing pair rows and blank decisions stop execution before downstream computation. The decision files are review data, not shell commands.

Existing `analyze`/`compare` commands remain available for uninterrupted diagnostic runs. They bypass these manual gates and should be used explicitly for that purpose. The staged workflow does not yet replace family/copy validation with an automatic classifier: human review supplies that decision.
