# Stage 2: selected windows → gene models → STOP 2

This release handles **protein-coding gene models**. Noncoding-RNA prediction is not implemented. This command accepts the locus-first `stage1` checkpoint. It does **not** use the older model-first `contexts` input format. It stops before alignments; [stage 3](STAGE3.md) consumes its reviewed choices for conservation comparisons and selected-sequence exports.

## Prepare the protein reference once

Use a local catalogue built from a broad protein FASTA (e.g. a versioned vertebrate RefSeq export), or explicitly download reviewed vertebrate UniProtKB proteins and their available isoforms:

```bash
python -m genecontext prepare-stage2-db \
  --download-vertebrates --out /mnt/project/references/vertebrates_reviewed

# Alternatively:
python -m genecontext prepare-stage2-db \
  --fasta /mnt/project/references/vertebrates.faa.gz \
  --out /mnt/project/references/custom_db
```

The UniProt query is `(taxonomy_id:7742) AND (reviewed:true)` with isoforms enabled. This is the reviewed vertebrate subset, **not all vertebrate proteins**. A supplied FASTA is not assumed to be vertebrate-only. UniProt GN/OX fields, original accessions and descriptions are retained in SQLite. Accessions without gene names remain accessions. Duplicate IDs are rejected; different accessions with identical sequences remain distinct. The receipt records release/query/source SHA-256, preparation commands and hashes of every database file. A prepared catalogue is immutable; choose a new directory to update it. Analysis validates it and performs no network lookups.

Reference symbol case alone does not create a competing gene (`NKX2-3` / `Nkx2-3`). Original spellings, accessions and taxa remain in every near-best hit; paralog suffixes and distinct accessions stay separate. A proposed label uses the best-scoring reference spelling and is not an official naming or orthology assertion for the target species.

## Prepare Tiberius once

The tested CPU runtime and `vertebrates`/`mammalia` weights are shared between windows/species. No model per species is required. Configure an explicit non-home cache:

```bash
python -m genecontext setup-predictor \
  --predictor-cache /mnt/project/cache/tiberius \
  --models vertebrates mammalia
```

## Start from saved selections

```bash
python -m genecontext stage2 \
  --results /mnt/project/nanog_stage1 \
  --selection /mnt/project/nanog_selected \
  --genomes /mnt/project/genomes.tsv \
  --protein-db /mnt/project/references/vertebrates_reviewed \
  --predictor-cache /mnt/project/cache/tiberius \
  --out /mnt/project/nanog_stage2 --threads 8
```

`--selection` is a `stage1-select` output. Its decision snapshot, source evidence fingerprint and review ID must match the checkpoint. Alternatively use `--decisions loci.decisions.tsv`, or `--draft loci.review.draft.json`. A draft with any unsaved pair is refused unless `--committed-only` is explicitly set. With that flag, unsaved pairs and their checkboxes are omitted and listed in provenance. They never become biological approvals. `--only-genomes ID ...` explicitly narrows an approved set for a pilot.

The genomic FASTA is required even if stage 1 searched only annotation proteins. Dataset TSV fields are `genome_id,fasta` and optional matching `gff3,protein,cds,species,assembly_accession,annotation_provenance`. See stage 1 for path resolution. A stage 1 FASTA hash, when present, must match; contig lengths and annotation hashes are checked. For protein-only discoveries without a genomic hash, stage 2 records a **new** genome fingerprint and identifies that weaker cross-stage identity check explicitly. The caller must supply the matching assembly.

## Continue an interrupted run

Repeat the original `stage2` command with `--resume`, including the same `--out`.
The run checks source evidence, selections, genomic/annotation content hashes,
protein database, runtime/weights, tools, parameters and implementation before
reusing an extraction checkpoint or completed protein-search step. Every cache
artifact is verified by SHA-256. Changed inputs require a new output directory.
An already sealed run returns its existing status without rewriting the report,
review ID or decisions. `completed_with_errors` is also sealed: use a new output
to retry those errors; successful prediction entries remain reusable.

Extraction is saved per genome. Tiberius results are shared under
`--predictor-cache/predictions/`, keyed by exact DNA and window identifier,
verified weight content, runtime/source fingerprints, threads, batch size and
worker mode. Cache hits still undergo the native CDS/protein round trip. A
failed or incomplete inference is never published. Concurrent requests lock the
same cache key and publish complete results atomically. Interrupted scratch
attempts remain under `.stage2/attempts/`; they never become human approvals.

`--tiberius-worker-mode reuse` is optional. One process retains a single neural
network; an unchanged model avoids loading the weights again. Each request has
one window, an independently clamped sequence length, a fresh `PredictionGTF`
and fresh HMM. Switching models clears the retained network. CPU proof and
model-loading/prediction timings are recorded per window. The compatibility
mode `isolated` remains the default. A synthetic host2 control compares exact
GFF3, pre-filter GFF3, CDS and protein exports, including short → long → short
windows and model switching; this is not a full NANOG performance benchmark.

`compact-report --report RESULT/report.html --out RESULT/report.compact.html`
creates an optional smaller offline transport beside the original. It restores
the exact application script at load time; scientific evidence, stylesheet,
review IDs, labels, zoom and manual decisions are not rewritten. This measures
file-size reduction, not lower browser memory or visual acceptance.

## Computation and outputs

1. Validate selections before starting. Preserve the exact region/anchor/window IDs and coordinates, including split anchors and custom windows.
2. Stream the selected intervals from plain or ordinary gzip FASTA. Do not decompress or index an entire genome into scratch. All selected windows from one genome are extracted together, along with full CDS segments needed by overlapping annotation models.
3. Keep every overlapping annotated coding transcript. Models extending outside the display window retain their full CDS and an edge flag; they are not silently clipped into apparently complete genes.
4. Merge overlapping selected windows into bounded canonical prediction intervals on each assembly/contig, preserving each original display window and anchor. Refuse a transitive merge exceeding `--max-window-bases`. Run DNA-only Tiberius once on each canonical interval and project complete models back to the overlapping display windows. Independent upstream sequence-length clamping is retained. CPU execution and CDS/protein/native-coordinate round trips are checked by the existing pinned backend. No protein/annotation hints are supplied to the predictor.
5. Identify usable proteins with BLASTP against the prepared broad database. Defaults: E ≤ 1e-5, identity ≥ .30, coverage of **both** proteins ≥ .50. Preserve all retrieved passing hits and alternatives scoring within .95 of the best. Competing names yield `ambiguous_family_assignment`, not a forced name. Identification is homology evidence, not orthology or expression proof.
6. Separately compare each model against **all** stage 1 target queries. Retain per-query support; it is not reused as independent gene-name confirmation.
7. Compare local model proteins across genomes. Deduplicate transcript/overlapping-annotation alternatives into candidate gene groups before proposing neighboring-gene correspondences. A prediction overlapping several annotation genes remains unresolved and does not merge them. Near-best alternatives and non-reciprocal assignments remain ambiguous. Reciprocal proposals within the selected catalogue do not establish genome-wide orthology.
8. Stop at `report.html`. Errors preserve available evidence, are shown in the report, and return a nonzero analysis exit code. An existing annotation may still be available after a predictor failure; the failure remains explicit.

Main files:

- `selection.snapshot.tsv`, `selection.input.json` or `.tsv`: normalized decision rows (including manual window overrides) and the exact source decision file.
- `report.html`: independent annotation/prediction tracks, models, target-query support, broad-database alternatives and pairwise neighborhood proposals.
- `regions.tsv`, `windows.fna`: selected native regions and genomic DNA.
- `models.tsv`, `cds.fna`, `proteins.faa`: all model alternatives, quality and sequences.
- `models.gff3`: combined GFF3 with genome-namespaced sequence IDs. `annotations.tsv` maps genomes to individual GFF3 files under `annotations/`, with native contig IDs.
- `evidence.json`: authoritative models, all near-best assignments, selection provenance, parameters and timings; raw searches remain under the search subdirectories.
- `models.decisions.tsv`, `regions.decisions.tsv`, `neighbors.decisions.tsv`: initially pending/empty; excluded from the evidence seal so they remain editable.

`stage2-report --results RESULT --out NEW_PRESENTATION` refreshes the offline report without prediction or BLAST. For an older coding-only checkpoint, add `--genomes /mnt/project/genomes.tsv` to attach all overlapping gene-level features from the original GFF3 files. Each GFF must match the original annotation SHA-256. The sealed `annotation_context.json` supplement is bound to the unchanged evidence hash; CDS models, scientific evidence and review IDs remain unchanged. New stage-2 runs retain these features automatically. It preserves the stage 2 review ID and pending/manual TSVs. Portable report drafts include that exact review ID; another checkpoint's draft is rejected.

To repeat broad-database naming after updating the catalogue or identifier rules, use `stage2-reidentify --results RESULT --protein-db DATABASE --out NEW_RESULT --threads 8`. It verifies the sealed parent, preserves genomic windows, prediction models, target-query support and neighbor searches, and repeats only broad BLASTP identification. The parent stays unchanged. The new evidence gets a new review ID and pending decisions, so previous model approvals cannot silently carry over to revised names.

## Second manual review

**Compare models**, next to a locus's model selector, opens one compact comparison above the decision table. It shows source, protein length and CDS blocks for the selector's candidates on one shared linear genomic ruler. Introns retain their native lengths. Exact CDS alternatives collapse only when their physical gene, native coordinates, strand, CDS/stop blocks and CDS sequence fingerprint agree; all IDs remain in the tooltip and original selector. **Select** changes the proposed model and preserves pending review text; the gene decision still needs **Save**. **Close** hides the comparison. The panel adds no predictions, alignments or automatic approvals.

The read-only **STOP 1 comment** block preserves the original selection note and reason. It appears below the selected-locus controls in each comparison viewer and beside each main locus in **Genes for the next stage**. A shared genome-selection note applies to every locus selected by that STOP-1 decision; it is not an independent gene-specific comment. The editable STOP-2 review comment remains separate and starts from the STOP-2 draft only. Viewing an inherited comment does not save or approve a new decision.

Use **+ Add comparison** to stack additional independent viewers below the first. Each viewer has its own selected locus, compared loci, source filter, view, zoom and scroll position. **Compare with → Choose loci…** opens checkboxes for a custom subset, so different copies can be compared in separate viewers. **Remove comparison** removes only that extra viewer. Viewer settings are saved in this browser separately from review decisions. All viewers share the same evidence and the single **Genes for the next stage** decision table; clicking a gene opens its shared review cards. Adding, changing or removing a viewer does not select models for the next stage.

Locus selectors, custom-subset checkboxes, comparison headings and diagram rows show the focal gene name from the proposed or manually selected main model. A KLF2 candidate found in a KLF4 query panel therefore remains labelled KLF2. Existing provisional names remain provisional; a locus without a target model says so. Names are display labels, not new orthology assignments.

The main diagram follows target orientation and defaults to annotation + Tiberius. Gene clusters focuses on the shown neighborhood with its own exact kb ruler per row; source filters retain this same focused range. Coding blocks show source-wise CDS coverage and thin lines show genomic spans. Overlapping genes occupy separate lanes. CDS model details keeps a shared bp scale; gene-order view remains schematic. Choose five or ten neighbors per side, all loci, one comparison locus, or a custom subset. Source-specific links use saved local BLASTP hits for the displayed source; this does not rerun BLAST or infer orthology. Noncoding GFF context is hollow and read-only. Curves between rows mean saved protein correspondence to the focal locus, not a new direct pairwise alignment. The full DNA-window coordinates, orientation and gene counts remain available. Hidden groups and partial models do not establish gene absence. All alternatives remain in their physical-gene card and the advanced library; display controls do not approve models or change matching thresholds.

The upper panel is for visual inspection; it keeps its own scroll area so the decision table stays close below. **Genes for the next stage** lists every STOP-1 locus at once. Each starts with one proposed main CDS, not all genes in its window. STOP 1 selected genomic intervals, not final CDS isoforms. The proposal is an overlapping, same-strand model ranked by target warnings, coding completeness, window-edge status and saved query bit score, with annotation as a tie breaker. An existing explicit model choice takes priority. A proposal does not approve or export itself.

For each row choose **Use for next stage**, **Unresolved** or **Exclude**, then **Save**. Use **View** to inspect that locus above; changing the viewed locus does not change decisions in other rows. Complete unflagged models need no typed explanation: Save records the visual-review outcome. Partial, edge or unsupported target models still require an actual comment and, where indicated, an explicit expert exception. The model selector offers the exact source/CDS alternatives; evidence and comparison-group options stay collapsed. Editing a saved row makes it unsaved and blocks export until it is saved again. Pending comments and form selections survive navigation and draft reload.

A protein-database symbol difference alone does not require an expert exception when the model matches both a named target query and a retained database protein at exactly 100% identity and full coverage in both directions. Evidence / options explains this case; source names remain unchanged. Partial matches, high but non-identical similarity, and unlinked query IDs do not clear family warnings. Other model-quality and locus checks still apply; identical proteins alone do not establish copy orthology.

**+ Gene** adds one coding gene from that row's window to the table. Added genes default to **Synteny anchor only**, and are not included until their own decision is saved. Choose **Paralog / target copy** when adding another target copy: it joins the same comparison group and uses its own model interval as the seed. Target-query and family checks still apply. **Also compare another gene** uses a separate comparison-group label and retains the explicit exception checks for unsupported or off-anchor targets. Added-gene purpose is visible directly in its table row. Copy labels under Evidence / options identify each physical locus independently of the comparison group. Adding a gene does not replace the main CDS. The all-model library and original detailed controls remain under advanced review. Excluding or leaving a primary locus unresolved clears its retained model and pair selections. Existing drafts keep their exact decisions; the optional `review_plan` field remembers presentation choices but stage 3 consumes only the authoritative saved `models`/`regions`/`pairs` sections. Disrupted CDS cannot enter comparison. Partial/edge model use requires an explanation where flagged. Unsupported or off-anchor target choices require an explicit expert override and reason. One physical gene has one selected representative; distinct physical copies in a species may share one comparison group, with separate copy labels (gene ID and coordinates by default). The collapsed pair-review table retains the original all-model proposals and their export IDs; source-specific diagram links are presentation evidence, not new review approvals. Pair follow-up requests are separate from target-model decisions. Nothing is selected automatically, and an unselected model is recorded as **not selected**, not as individually reviewed/rejected.

The primary validated export is one `stage2.decisions.json`; optional TSVs are available separately. The report checks target/query/anchor conflicts, explicit expert overrides, duplicate physical targets, comparison labels, partial-model explanations and neighbor-pair region states. [Stage 3](STAGE3.md) repeats these checks on the backend and consumes the JSON. The older `compare-reviewed` command retains its legacy contract.

## Resource policy and scope

Tiberius uses one canonical prediction interval per worker, batch size one by default, and `--threads` controls CPU use. A thread count is **not** a hard memory bound. On host2 launch **all** experiments, database preparation and predictions through `gene_context_experiments/submit.py`, under the existing aggregate 250 GiB cgroup and zero swap. Keep source snapshots, databases, temporary files and results on `/mnt`. Do not change the `current` symlink merely to run a pilot: use the requested source snapshot as the command's working directory and PYTHONPATH.

This initial implementation does not auto-expand edge windows or add a miniprot rescue. Edge models, missing predictions and disagreements remain visible for manual review. It does not claim to recover every gene or every splice isoform ab initio. After STOP 2, `stage3` performs protein/codon and genomic comparisons and draws the explicitly chosen neighbor context. RNA evidence and biological orthology decisions remain separate from these software measurements.

The [first real STOP-2 pilot](STAGE2_PILOT.md) documents measured runtime, database scope and annotation/prediction disagreements on three saved NANOG regions.


## Genomic ruler in STOP 2

CDS model details uses one shared ruler based on the saved review windows and the viewer width. Switching locus, comparison, annotation/Tiberius source or neighbor count does not change the bp scale at a given zoom. The selected windows keep their exact native endpoints; longer windows occupy proportionally more space. Long names are shortened inside the overview window when necessary, with full names retained on hover and in the model cards. Label packing never stretches genomic distances. Zoom preserves the genomic viewing position; species captions remain visible during horizontal scrolling.

The report opens in **Gene clusters (kb)**. Each row focuses on displayed neighboring genes, omitting unused DNA-window flanks only from the overview. Each row retains its own exact linear kb ruler and native orientation. Overlapping genes and their labels occupy separate lanes; row height grows as needed. Annotation and Tiberius have separate source tracks. Thin lines show genomic spans / introns; solid or outlined blocks show the union of CDS coverage from that gene's models in that source. This is a gene overview, not an asserted transcript. Badges give source model counts; individual models remain in the inline comparison. Correspondence curves connect gene centres and do not represent a DNA alignment. Pixel lengths across rows are comparable only through the respective rulers. No coordinate, sequence, review ID or decision changes.

One compact track combines annotation and Tiberius under a single gene label. Identical full native-coordinate CDS structures and sequence fingerprints within that gene group are drawn once with an **A+T** badge. Distinct alternatives retain separate thin strokes with **A** (annotation) or **T** (Tiberius) badges. Mere overlap, matching symbols or identical proteins do not merge structures; differences outside the clipped window also remain distinct. Every original model stays separately reviewable in its card. Filled blocks are CDS segments at exact genomic positions and widths, with connectors over intervening genomic sequence. Tiny CDS blocks are not widened; use Zoom 2×/4× and horizontal scrolling. Gene-level noncoding GFF features are hollow spans, not invented coding/exon models. Names and LOC identifiers are printed separately in full; collision packing moves labels/models vertically without altering x coordinates. Hover a model for its original transcript ID and native coordinates; click a selected-locus gene to open its existing review card.

Saved protein correspondences are drawn in the gaps between panels at the displayed gene-centre coordinates. Their hover text identifies the underlying model pair. These curves are candidate correspondences, not genomic DNA alignments or proven orthology. As before, no more than ten neighboring gene groups per side are drawn; hidden counts and all model-review cards remain available. **Gene order (schematic)** remains an explicit optional view with equal spacing. Refreshing this presentation does not search or predict again and retains the checkpoint review ID and decisions.


A source filter that hides available models now produces an explicit notice and a **Show both sources** button; missing annotation must not look like a failed target search. Tiberius-only overview arrows use a stronger outline and light fill. A report URL may include `?locus=REGION_ID` to open a specific review locus. This changes only navigation and preserves all existing review decisions; unknown locus IDs fall back to the normal initial selection.

## Specific copy names

Reports keep a specific locus label separate from protein-database similarity and the selected comparison group. Chicken LOC101748883 is displayed as **NANOG-like**, retaining the original identifier. Unnamed homologues remain family candidates; a shared best hit does not name a particular paralog. `gene_identities.tsv` records labels, source IDs, naming evidence and discrepancies. See [copy-identity rules](COPY_IDENTITIES.md).

### Model integrity and recovery warnings

`ORF intact` means that the extracted DNA has an intact coding frame with start
and stop codons. It does not establish that a full gene was recovered.
STOP 2 separately checks target-query coverage, retained annotation CDS bases,
possible splitting into separate predicted ORFs, and source GFF exceptions.
These warnings are also enforced by the stage 3 decision validator. A target
with warnings requires an explicit expert override and a written reason;
no candidate is silently discarded or approved.

The current review thresholds flag predictions retaining less than 80% of a
comparator annotation CDS, and predictions covering less than 50% of their best
query unless they exactly match an annotation CDS. These are review triggers,
not orthology or absence thresholds. Separate ORFs are flagged as a possible
split when they jointly retain at least 80% of one annotated CDS while each
retains less than 80%; overlapping isoforms and repeated observations are excluded.

For existing checkpoints, recover GFF quality metadata with:

```bash
gene-context stage2-quality --results results/pou5_stage2 \
  --genomes genomes.tsv --out results/pou5_quality
```

The small `annotation_quality.json` is bound to the evidence and annotation
SHA-256 fingerprints. Place it beside `evidence.json` in the checkpoint being
reviewed, retain its fingerprint in the presentation receipt/seal, and refresh
the report. New stage 2 runs retain these flags directly. This operation changes
quality metadata, not predicted exon coordinates or scientific review choices.

A split or short prediction needs model reconstruction with independent evidence
(e.g. protein-guided alignment or transcript support) or selection of a supported
annotation alternative. Do not concatenate predicted proteins to repair a split.

### Reannotation bounds from preliminary synteny

New runs default to `--window-policy annotation-neighbors`. For each reviewed
locus, stage 2 reopens the exact matching full GFF3 and selects ten distinct
coding genes on each side, without a distance cutoff. Sparse older STOP-1
contexts do not limit the new search. Multiple transcripts of one gene count
once; overlapping genes and query-source copies remain display-only. The full
spans of the refreshed neighbors, the target anchor and already displayed
preliminary genes define the DNA interval. Noncoding genes remain visible
context and do not fill the coding-neighbor quota.

Before prediction, the complete window must capture at least five coding genes
per side (the requested default is ten). A missing annotation, short contig or
manually narrowed window remains explicit and blocks automatic prediction.
`--allow-incomplete-neighborhood` is an explicit exception, recorded through the
command and coverage receipt; it never waives the DNA size limit. Manual window
bounds retain priority and are not silently enlarged.

`stage2-plan` reads annotations only and requires `--genomes` for this default
policy. It writes the complete candidate intervals even above the execution
budget, with `ready_for_prediction=false` and exact blockers. No sequence
extraction, Tiberius or BLAST starts during planning. A window above
`--max-window-bases` is retained for review, never clipped to look acceptable.

Use explicit `--window-policy pre-synteny` to reproduce the historical POU
policy: reuse the full saved STOP-1 envelope. This legacy mode does not refresh
missing neighbor genes. `--window-policy fixed` retains configured flanks.
Historical preliminary searches started with five coding neighbors per side;
current count-first searches request ten. Fewer may be available. Isoforms do not
count as separate neighbors. Bounds are native, zero-based, half-open and do not depend
on the strand or the display zoom. No extra 50 kb flanks are added to this envelope.

In the explicit historical `pre-synteny` mode, when no preliminary
neighbors are available the configured flank window is retained as a recorded
`fixed_fallback`, rather than reducing the extraction to the target alone.
`--window-policy fixed` explicitly restores the previous flank policy. Invalid
or mismatched preliminary coordinates are errors, not an excuse to silently use
another locus. The default 2 Mb cap also applies to merged prediction intervals;
oversized windows require explicit boundary review or a higher cap, never an
automatic cut through an outer gene.

Inspect the actual planned bounds before running the predictor:

```bash
gene-context stage2-plan --results results/stage1 --selection results/selected \
  --genomes genomes.tsv --out results/stage2_window_plan
```

`windows.tsv` lists the coordinates, length, number of preliminary neighbors and
boundary source. `window_plan.json` retains the included gene IDs and coordinates,
the original window, source checkpoint fingerprint and merged prediction intervals.
The plan performs no DNA extraction or prediction. The same planner runs inside
`stage2`; its resolved bounds drive DNA extraction, annotation collection and
Tiberius. The original review decisions remain unchanged.

Changing these bounds requires a fresh stage 2 output and a new STOP-2 review ID;
expanding an old report does not reannotate the extra DNA. Stage 1 discovery does
not need to be repeated. Old STOP-2 decisions are not automatically transferred.
The new interval can include additional genes between the preliminary anchors;
Tiberius is not required to reproduce every supplied annotation model.

For historical checkpoints and the explicit fixed policy, the reannotation window
is the selected gene interval plus its configured flanks.
Its size is not the CDS length: long introns can make a short protein span hundreds
of kilobases. The overview now uses pale genomic spans with filled CDS blocks for
spliced models, and reports span/CDS lengths for strongly intron-rich targets.
Zoom 0.25×/0.5× reduces the shared scale for all rows equally.

In those historical/fixed checkpoints, STOP-1 context can extend beyond this window. It is carried
forward with the source review ID and evidence fingerprint, and shown separately
under “Broader annotation context from STOP 1”. Its matches remain prior evidence;
they are not Tiberius predictions in the current window or new review decisions.
If the local window contains no neighboring genes, local synteny is explicitly
unassessable. Target-protein similarity ribbons alone do not establish synteny.

Conflicts between a single annotated transcript and separate predicted ORFs now
require explicit review for **both** sources. The annotation is not exempt from
this check. A conflict can reflect an annotation misjoin, a split prediction, or
insufficient evidence; the check does not decide which source is correct.

Coding-model highlights cover CDS segments, rather than filling the full intron
span. Spliced models use a thin backbone and individual CDS blocks in the cluster
view; native genomic coordinates and the complete extraction window remain exact.
A large annotated genomic span must not be cropped solely because it looks long.

### Reconstruct a full model from a supplied transcript

For a conflict between an annotation transcript and split DNA-only predictions,
project a full reference RNA onto the saved genomic window:

```bash
gene-context stage2-reconstruct --results results/pou5_stage2 \
  --transcripts references.gb.xml --mapping transcript_models.tsv \
  --out results/transcript_reconstruction --threads 4 --attach
```

The TSV has two columns, `model_id` and the exact transcript `accession` (with its
version). GenBank XML supplies the full RNA and transcript CDS coordinates. There
are no downloads during analysis. `minimap2` must be installed; `--minimap2` can
name its executable and `--max-intron` sets the mapping limit (default 1 Mb).

This command requires an exact end-to-end, unambiguous splice alignment, projects
the CDS on either strand, checks an intact ORF, and compares every CDS boundary
and sequence to the existing annotation model. Different candidate structures
are rejected for automatic attachment. It exports native GFF3, CDS and protein
FASTA, plus mapping evidence. It never concatenates raw predicted ORFs.

`--attach` adds small sealed support/sequence files to the existing checkpoint;
model IDs, sequences, coordinates, windows and manual choices remain unchanged.
A reference-derived transcript alone is a consistency check. A source record
reporting RNA-seq support for **all** introns can resolve the structural split
warning for that matching annotation model; raw predictions remain flagged.
The interface distinguishes this reference-backed full model from Tiberius.
RNA-seq support is a source statement, not a new RNA-seq analysis or proof of
orthology. The user still chooses whether the model enters stage 3.

Attachment also refreshes the report in place using saved searches; it does not
rerun gene discovery or Tiberius. Repeating the same transcript/model inputs is
idempotent with respect to model sequences and decisions. A different attachment
requires a separate reviewed revision.

### Compact display and automatic size comments

Show the supplied annotation and Tiberius predictions. When pre-filter GFF output is available, removed predictions overlapping target CDS are also drawn as read-only dashed fragments. They are not selectable models and do not enter CDS/protein exports. Miniprot query projections are not drawn in the STOP-2 overview.

Large differences are recorded in automatic comments, separate from the reviewer's notes and decisions. Default: at least 2× between the per-source maxima, with an absolute difference of at least 10 kb for genomic CDS-model span or 50 aa for protein length. Span includes introns; it is never substituted for CDS/protein length. These thresholds flag a difference, not a failed gene or an automatic rejection. Automatic comments are available under Evidence / options, exported in the decisions JSON, and carried into the next-stage decision audit. No extra warning badges or rows are added to the diagram.

New Tiberius runs retain `predictions.gff3.raw.gff3` before the predictor's own coding-length and in-frame-stop filters. Old runs without that file do not imply an absence of filtered predictions.

### NANOG-style overview

The default genomic-distance overview uses filled annotation arrows, outlined Tiberius arrows and correspondence ribbons. Short gene symbols are drawn in two compact label rows. Full descriptions, source badges and naming alternatives remain in tooltips and model cards. Exon structures remain in CDS model details. This presentation keeps native window bounds, filtered display-only fragments, collapsed automatic comments, and manual selections unchanged.

At the default **Fit** zoom, each cluster window fills the same available row width. The longest reviewed window no longer shrinks the other cluster rows. Higher zoom levels enlarge all rows together and allow horizontal panning. The overview represents genomic spans including introns; equal protein lengths do not imply equal arrow lengths.

Each cluster row has an independent target-relative kb ruler drawn by exactly the same linear transform as that row’s model arrows. The full DNA window length is shown in the row caption. Absolute assembly coordinates remain in boundary tooltips and locus controls. The span includes introns and is not CDS or protein length; fragments are not merged into a model.

The default cluster overview shows short gene symbols only. Family qualifiers and alternative-name lists are omitted from the diagram; unknown LOC identifiers are explicitly abbreviated with an ellipsis. Complete names and source details remain in arrow tooltips, accessible labels, and model cards. Clicking a comparison-row arrow focuses that locus and opens the corresponding model card. Rows now expand only to separate overlapping genes and labels; native coordinates remain unchanged.

The cluster overview keeps short names and source counts, replaces opaque intronic arrows with thin spans and actual CDS coverage, and retains a separate exact kb ruler per row.

### Row-specific scale in Gene clusters

Gene clusters maps the selected neighborhood linearly to the available row width, using bounds calculated from both sources and their display limits. Source filters do not change that row's scale. Fit uses the actual panel width; zoom and pan remain local to each comparison pane. Long introns cannot mask nested genes because filled glyphs cover CDS only and overlapping groups occupy different lanes. The full DNA-window coordinates remain available in the row tooltip and window details. The CDS model details view retains its shared bp scale. No genomic bounds, models, evidence, review IDs or decisions change.

### Run STOP-1 groups separately

A **Download group JSON** file from STOP 1 can be passed directly as `--draft`. It selects its saved `active_analysis_group`; the normal database, predictor and dataset options still apply. Give each group a separate output directory. The resulting STOP-2 report shows the group name alongside the source gene.

For a full grouped draft, choose a group by stable ID or exact name:

```bash
gene-context stage2-plan --results saved_stage1 --draft all_groups.stage1.draft.json \
  --analysis-group 'Family A' --genomes genomes.tsv --out family_A_plan
```

`stage2` accepts the same `--analysis-group` option, including with a grouped `--decisions` TSV. Omitting it on a full grouped input is an error, preventing accidental joint processing. A scoped group JSON cannot be switched to another group by a conflicting flag; use the full draft to switch. Pending decisions in other groups are outside the run. Unsaved member pairs remain errors unless the existing explicit `--committed-only` option is used. Membership never approves a locus. Rows outside the chosen group are recorded as outside the analysis scope, with their original source decisions retained; they are not biological exclusions.

Discovery and review IDs are reused from the unchanged stage-1 checkpoint. The default full-annotation planner refreshes neighbor coverage; it retains complete earlier context rather than trusting a sparse old shortlist. Prediction windows are planned only for the chosen group, so even overlapping loci from different groups cannot merge into one run. Each resulting STOP-2 report groups and compares that group's loci independently.

### Alternative CDS candidates and interactive alignment

`Compare models` compares isoforms of one physical gene. The shared ruler has
an exon view with marked compressed long introns, and an exact genomic view.
Models are stacked as full-width cards and automatically fit the panel. Length labels wrap below each track at a readable fixed font size; no horizontal scroll is required. The shared ruler and all tracks use the same transform.
`Include in alignment` records `role: alignment_candidate`; it does not approve
that CDS as the target and does not replace an approved main model. Several
alternative CDS models of the same gene may be retained. Save the region and
export the matching decisions JSON; the native protein/CDS QC still applies.

Every new STOP2/STOP3 report includes a portable local alignment app. Extract its
ZIP and run `Start_review.command`, `Start_review.sh` or `Start_review.cmd`.
No Python, Conda or Docker is needed. A suitable installed MAFFT (7.487+) is used
first. Otherwise the pinned official MAFFT package is fetched only after the
explicit Recalculate button. The runtime is shared by concurrent report sessions;
Finish and clean up waits for active jobs, saves results and removes only the
last session's owned temporary runtime. User sequences, choices and results stay.
The included launchers cover macOS, Windows and Linux, x86_64 and ARM64. Automatic
portable MAFFT fallback is available for macOS and Windows and Linux x86_64;
Linux ARM64 needs a native MAFFT already on PATH. Windows ARM64 fallback depends
on the system's x64 compatibility. Native Windows execution requires its own QA.

Compare models labels each model's coding exon intervals E1, E2… and intervening
introns I1, I2… in its own 5′→3′ order. Adjacent/overlapping CDS features are
merged for genomic structure display. Every length chip is an actual genomic bp
length; row totals show covered CDS bases, intron bases and first-to-last-CDS
span. Their invariant is CDS + introns = genomic CDS span. Full mRNA/UTR exon
boundaries are not inferred. Details use 1-based inclusive coordinates.

All rows share one genomic coordinate transform. In the readable view, only
shared noncoding intervals longer than 200 bp are compressed (up to 22 pixels per gap, with at most 35% of the available track reserved for gaps);
exons and other uncompressed DNA use one bp scale. An intron retained as coding
sequence in another model can therefore be only partly compressed. Break marks
follow the exact shared compressed intervals. The exact genomic mode compresses
nothing. Screen compression never changes the displayed bp measurements.

Each independent synteny panel has +/− buttons, a continuous zoom slider and Fit.
Ctrl/Command + wheel zooms around the pointer; ordinary wheel scrolling remains
scrolling. Drag the diagram background to pan. Genomic zoom preserves native
bp lengths on the appropriate ruler, the vertical review position and existing
manual decisions. Schematic gene-order view supports visual zoom as well.

Hosted reports connect to the loopback review runner automatically. Select model
membership and press **Recalculate alignment**. MAFFT runs through the existing
bounded launcher; status, exact source IDs, checked results and exports persist
outside the sealed scientific checkpoint. Reloading a hosted report reconnects
to its running or completed job without resubmitting. Changing membership and
pressing the same button creates a new alignment. Identical in-flight/verified
requests reuse one job. A different browser workspace has a separate latest-result
pointer. Source protein/CDS checks and result seals are required before display.

For a standalone ZIP, start its bundled Start_review launcher once. Recalculate
then uses that local runner, and a completed alignment survives HTML reload in
the same browser. A raw file:// HTML cannot execute MAFFT by itself. Hosted and
portable comparison choices remain separate from biological final approval.


### Inline CDS model review

Compare models expands within one physical-gene card. All alternative structures share a coordinate transform and compact stacked tracks. The physical-gene dropdown is removed; it is determined by the card. Intron lengths, coordinate lists and advanced scale options are collapsed by default. Checked isoforms stay in this same card as alignment candidates; they do not create extra gene cards. Set main changes a proposal, not approval. The common Save records the gene decision and checked candidates with the reviewer comment. Existing partial, disrupted-CDS and expert-override checks still apply. Review IDs, exported decision schema and prior choices are preserved; no calculation starts from opening or saving the comparison.
