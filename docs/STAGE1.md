# Stage 1: discover candidate loci, check available contexts, then review

The primary route is **BLASTP against all annotated protein products**, then exact annotation identifiers locate each gene in the assembly. TBLASTN remains available for missing annotations and rescue.

Input: the curated FASTA variants of **one target gene/search group**, a genome list, and the available source datasets. Output: one offline report where the user selects genomic loci. No candidate is approved automatically.

## One command

```bash
gene-context stage1 --gene NANOG \
  --query NANOG_query.fasta --query-type protein \
  --genomes genomes.tsv --out results/nanog_stage1 \
  --jobs 8 --threads 8
```

`--threads` is per genome; `--jobs 8 --threads 8` can use 64 tool threads during discovery. The annotation-based early-synteny pass runs once in the parent after the genome workers finish. On host2, submit every computation through the constrained launcher described in [LIMITED_EXPERIMENTS.md](LIMITED_EXPERIMENTS.md): all concurrent descendants share 250 GiB RAM and zero swap, and data/cache/output must stay on `/mnt/data`.

The genome TSV requires `genome_id` and `fasta`; useful optional columns are `gff3`, `protein`, `cds`, `species` and `assembly_accession`. Paths are relative to the TSV or absolute; gzip files are accepted. Each row identifies one actual assembly. For protein-first search, `fasta` may be empty or point to a not-yet-downloaded assembly when GFF3 and protein/CDS are supplied and GFF contains contig lengths; genomic DNA is needed later for sequence extraction. `protein` is the annotated protein FASTA; `cds` is the annotated CDS FASTA. These files let discovery and the automatic context pass avoid reading the entire genome. If only GFF and genome are present, annotated coding products are reconstructed once and cached. Annotation is optional for the genomic discovery route. If it is missing, the report explicitly marks the context evidence unavailable.

Use `--datasets datasets.tsv` when the complete annotation/protein inventory is separate from the target genome list. Use `--source-datasets references.tsv` for additional source-reference assemblies that are not study targets. The program resolves query origins from identifiers, exact sequences and supplied curated GeneIDs in the provided data. A GeneID can bind the source neighborhood across annotation isoform differences; sequence conflicts are retained and do not become sequence-verified bindings. It preserves ambiguous matches and reports unresolved origins. No manually constructed source-locus list is required. Standard FASTA headers are read automatically: `[organism=Species name]` / `[Species name]`, explicit `GeneID`, and versioned RefSeq protein/transcript accessions. The full source header is retained. An exact versioned accession mapping uniquely to one annotated source locus can bind its neighborhood across sequence differences; such a binding is explicitly context-only, never sequence-verified. Ambiguous accession mappings stay unresolved. No source gene is guessed from its symbol or a best BLAST hit. `--query-metadata queries.tsv` can supply or override source species, gene symbols and accession provenance; if supplied, it must describe every input FASTA ID exactly once.

## What runs

1. The default `--discovery auto` uses BLASTP for protein/CDS queries when a GFF is supplied. **All usable annotated isoforms** enter the search; identical full protein sequences are searched once and keep every transcript/gene binding. CDS queries are translated first. `--discovery blastp` requires annotation, while `--discovery genome` explicitly requests whole-genome TBLASTN. In auto mode a genome without GFF uses TBLASTN. DNA queries always use BLASTN.
2. BLASTP HSPs are chained within each query/subject protein. Exact product → transcript → gene relationships map them to genomic loci; all matched isoforms of one annotated gene share one locus. Separate annotated genes stay separate, even if their proteins are identical. Genomic searches instead chain compatible query/genome HSPs and group overlapping aligned blocks. No gene-name lookup substitutes for sequence search.
3. The original discovery order uses the strongest single-query chain score within each genome. Overlapping query positions count once. The review adds a combined order from coverage- and identity-weighted support across source groups; see **Combined review order** below. Original scores, candidates and raw BLAST evidence remain available.
4. Available reference and candidate neighborhoods are checked automatically. Annotation and proteome indexes are cached. Reference neighbor proteins are matched against target proteomes in batches; coordinates, order, side and direction provide preliminary context evidence. Target genes are not counted as their own neighbor evidence.
5. The report stops for human selection. No Tiberius prediction, target-protein MSA or comparative CDS alignment runs here.

A BLASTP no-hit or weak result does not silently trigger a potentially expensive whole-genome search. Inspect the result and repeat the relevant genome with `--discovery genome` if annotation rescue is needed. Invalid/mismatched annotation is an explicit error. Missing annotation products and identical proteins at multiple gene loci are reported.

The new default is `--refine none`. Annotated candidates use GFF gene boundaries; local miniprot refinement applies only to genomic searches. Optional `--refine miniprot` performs local spliced boundary refinement; these predictions stay in technical details and do not establish gene identity. Older imported checkpoints may already contain those predictions.

Early contexts start with five genes on each side, expanding to at most ten per side when cross-species preliminary support is insufficient, within 1 Mb. The hard maximum is ten genes per side; larger `--synteny-neighbors` or `--synteny-max-neighbors` values are rejected. Same-assembly controls remain visible but do not establish cross-species support or stop expansion. With controls alone, the status is `within_assembly_controls_only`. These are computational limits, not biological thresholds. The actual inspected span is shown; it is separate from the proposed ±50 kb window for detailed downstream work. `--neighbor-cache` selects a reusable cache. By default it is `_genecontext_cache/early_synteny` beside the output, not an implicit home-directory download cache. No reference/model downloads are performed by this pass. `--early-synteny off` requests discovery alone.

## Read the report

Open `report.html` and choose a genome. The first ten candidates appear in combined review order; **Show next** reveals the next ten, and selected loci remain visible. The unit of review is one candidate locus. Its header uses an annotation-derived gene name when available; an unassigned locus is labeled **Candidate + ID**, never the search gene's name.

Each card has one compact query table: visible source species/gene and accession, query-coverage percentage with a position bar, aligned identity, and chain score. Every supplied query variant remains represented, including unmatched variants. The displayed product is the highest-scoring non-ambiguous chain passing all current display thresholds; if no chain passes, the strongest retained chain is shown with a warning. **Passing product shown** explains when this differs from the highest raw-score product. The query-accession disclosure exposes the product, best HSP E-value, subject coverage/positions and the alternative raw-score product. Identical sequences share evidence and do not supply independent votes. Ranking still uses the saved evidence and does not change with display filters.

**Synteny** shows one diagram containing the candidate and every available query-source reference locus. Variants from one annotated source gene share a reference row. The identical candidate/source locus is hidden by default and can be shown with Show self comparison. Other copies in the same assembly remain visible. All controls remain in the saved evidence. Source availability counts and unresolved sources are explicit; a self comparison never replaces missing cross-species evidence. The default view shows up to ten saved independent neighbors per side; the controls beside the diagram switch to five, the saved ±50 kb subset, or a genomic-distance scale. Equal spacing means gene order, not genomic distance. Actual displayed candidate coordinates and the later ±50 kb extraction window are labeled separately. Historical checkpoints with more than ten neighbors remain intact, but display is capped at ten.

Solid links connect saved proposed protein matches through a shared candidate neighbor. Dashed links include same-assembly controls; amber links/outlines mark unresolved assignments. These are preliminary correspondences, not established orthology. Gray genes remain visible. Other query-source copies and overlapping annotated genes have purple outlines and are excluded from independent anchors. New context computations retain these in `display_only_genes`; for older checkpoints, known omitted query-source copies can be drawn from saved reference coordinates. This does not invent protein matches or rerun context inference.

**Details** builds its contents only when opened: ranking contributions, annotated CDS/genomic match positions, all chains/isoforms/HSPs, and editable boundaries/raw files. BLASTP alignments use amino-acid axes; CDS coordinates come from annotation. Genomic BLAST shows its DNA match coordinates and does not estimate coverage of an unknown complete gene. Candidate-DNA availability is indicated per locus; mixed BLASTP/TBLASTN runs retain access to the available genomic-search FASTA.

**Query origins and supplied datasets** is the single source-availability table. A query sequence alone does not specify genomic neighbors: source GFF/proteins/CDS must resolve its original locus. Exact sequence bindings are preferred; a curated GeneID can bind an annotation neighborhood when a supplied isoform differs between releases, with an explicit sequence-mismatch note. Symbols alone do not establish the source. Add source-only assemblies with `--source-datasets`. Source-data or resolver changes require `stage1-context`; `stage1-report` preserves scientific evidence.

Optional neighboring-gene/pair requests remain under **Protein match details and optional neighbor follow-up**, persist in draft JSON, and export as `neighbors.to_inspect.tsv`. They do not launch neighbor analysis and are not consumed by `stage1-select`. There is no second copy of the synteny diagram or genome-wide maximum of unrelated candidates' context scores.

Default display gates are coverage 0.20, identity 0.25 and best HSP E-value 1e-5. All must pass in one non-ambiguous chain. Below-threshold candidates remain available unless **Hide below-threshold candidates** is enabled; selected loci stay visible. Filters do not delete evidence. The draft records current filters, and each newly saved pair records its own filter snapshot. Older drafts retain unknown filter history rather than receiving an invented snapshot.

The standalone HTML uses content-based interning of repeated annotation models, HSPs, neighbor records and comparisons. Hydration preserves all report chains/products and exact coordinates; objects with identical identifiers but different contents remain distinct. Full scientific data and receipts stay in `evidence.json`, raw files and `manifest.json`. Closed per-locus details are constructed lazily. See [metric definitions](STAGE1_METRICS.md).

## Choose loci

For each gene/genome pair, check one or more loci and choose **Investigate selected loci**, or choose **Exclude** / **Unresolved** with no locus selections. Select a short reason; a comment is optional except for **Other**. Save the pair decision. Different CDS predictions stay inside the locus; this stop chooses genomic regions, not a definitive transcript.

The report supports a portable draft JSON and automatic browser-local saving. A draft is not an approved decision file. **Validate all pairs** lists every pending or inconsistent pair and previews the selected windows. **Validate and download decisions TSV** exports exact IDs only after all pair decisions are saved and valid. For an inspected pair, unselected loci are excluded from this run with `selection_status=not_selected_for_run`; the export does not claim that each was individually reviewed. Selected loci, explicit pair exclusion and unresolved decisions have distinct status values. The TSV also retains per-pair display-filter snapshots (`unknown_legacy` for old drafts). No original evidence or decision file is overwritten.

Validate or export the selected windows on the command line:

```bash
gene-context stage1-select --results results/nanog_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv --validate-only

gene-context stage1-select --results results/nanog_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv \
  --out results/nanog_selected
```

The output contains `approved_loci.tsv`, `synteny_windows.bed`, `validation.json` and a decision snapshot. It does not start downstream prediction. Anchor/window overrides are available in the locus technical details. The proposed ±50 kb window uses the selected anchor; it is clipped at contig boundaries. BLASTP defaults to the annotated gene bounds. Genomic-search anchors are similarity-region bounds, not guaranteed full-gene bounds. Annotation completeness still requires review.

Legacy TSVs with explanations in `notes` remain accepted. The new report also exports a short `reason` and repeated `pair_decision` metadata. The selector checks both per-locus decisions and pair consistency, reports all errors, and creates no selection output if invalid.

## Reuse a saved search

To add/recompute available early context and regenerate the report without another genome BLAST:

```bash
gene-context stage1-context --results results/nanog_search \
  --datasets datasets.tsv --source-datasets references.tsv \
  --neighbor-cache /mnt/data/project/cache/early_synteny \
  --out results/nanog_review --threads 8
```

To save different display thresholds without another search or context calculation:

```bash
gene-context stage1-filter --results results/nanog_review \
  --out results/nanog_review_cov20 --min-query-coverage .20 \
  --min-identity .25 --max-evalue 1e-5
```

Each command creates a new checkpoint and leaves the source untouched. Searches use a more permissive initial E-value ceiling of 1e-3; HSP seed gates and chaining limits are recorded in the manifest. A match omitted by discovery cannot be recovered by display filtering. Source errors and missing data remain explicit; no-hit results do not establish gene absence.

The new [`stage2`](STAGE2.md) command consumes the locus-first checkpoint and saved selections, reannotates windows and stops for model/neighbor review. The older `contexts` command uses a separate model-first checkpoint format. Final conservation comparisons remain subsequent work.

## Refresh the presentation only

```bash
gene-context stage1-report --results results/nanog_review \
  --out results/nanog_compact
```

This verifies the source checkpoint and copies its sealed artifacts to a separate directory. It preserves the exact scientific evidence, filters, locus IDs and review ID, so existing drafts and decision TSVs remain compatible. Only the report and presentation provenance change; no BLAST or context calculation is run. The source is untouched.

## Minimal card layout

The first ten cards are shown immediately. Each default card contains locus coordinates, one compact row per query (coverage, aligned identity, score), and two collapsed entries: **Synteny** and **Details**. The synteny overview connects shared candidate-neighbor correspondences between adjacent rows: solid for retained protein proposals, dashed when a same-assembly control is involved, amber dotted for unresolved assignments. Hover a gene or link to trace its candidate-linked chain. A missing correspondence breaks the line; a connection between reference rows does not represent a direct reference-reference search. Source species/labels are visible in each query row; the accession opens product details. Sequence evidence, model/coordinate controls and raw data links remain inside **Details**. No synteny diagram or long metric explanation expands automatically.

The v0.7.6 NANOG source-context run verified all five query-source neighborhoods across the 16 available annotated target genomes. Both chicken loci are GeneID-linked with sequence differences retained; human, mouse and zebrafish have exact identifier/sequence bindings. The 2,546 target loci and seven target-discovery exports are unchanged. This was a context recomputation only, under host2's shared 250 GiB / zero-swap limit. See the [actual-run receipt](benchmarks/nanog_sources_v076_2026-11-06.json).

## Combined review order

The report ranks each genome/search group using bit score, query coverage, aligned identity and support across source groups. For the best non-ambiguous bit-score chain for each query, compute `weight = bitscore × coverage × identity`. Keep the largest query weight per source species; species sharing an identical query sequence are merged so duplicate sequences do not add votes. Unknown sources use unique query IDs. Missing groups contribute zero. The combined score is `largest source weight + mean weight across all supplied source groups`. Equal scores use best HSP E-value, then raw bit score and deterministic genomic coordinates.

This is a transparent review heuristic, not a calibrated probability or an orthology decision. Synteny and display filters do not enter the score. **Details → How this locus is ranked** gives all contributions. `review_ranking.tsv` exports the displayed order; `loci.tsv` retains the original bit-score discovery order. `stage1-report` preserves evidence, locus/review IDs and existing decisions while adding this derived ranking to the sealed checkpoint.

## Specific copy names

Reports keep a specific locus label separate from protein-database similarity and the selected comparison group. Chicken LOC101748883 is displayed as **NANOG-like**, retaining the original identifier. Unnamed homologues remain family candidates; a shared best hit does not name a particular paralog. `gene_identities.tsv` records labels, source IDs, naming evidence and discrepancies. See [copy-identity rules](COPY_IDENTITIES.md).

## One source-resolution path for every gene

New stage-1 runs always extract explicit provenance from the supplied FASTA headers. A per-gene query metadata TSV is optional. Use a common `--source-datasets` inventory containing the source assemblies/annotations, shared across genes; the resolver discovers source loci from the actual supplied identifiers and sequences. Only sources needed by the current query panel become reference rows. The pipeline does not silently download missing assemblies or invent their neighborhoods. Missing source contexts are listed in the main notice and the full query-source table.

Source GFF and protein/CDS files can be downloaded once without the complete genome:

```bash
gene-context download --manifest source_urls.tsv --annotation-only \
  --include-reference --out references/source_annotations
```

The annotation-only manifest requires `genome_id` and `gff_url`; it can supply `accession`, `species`, `protein_url`, `cds_url` and `md5_url`. Official MD5 verification and cache reuse are identical to the ordinary download route. The generated `genomes.tsv` retains protein/CDS paths and can be passed to `--source-datasets`. No source full-genome FASTA is needed when GFF and usable protein/CDS are available.

For an older checkpoint that discarded descriptions, restore explicit metadata only from the unchanged original FASTA and rerun context without repeating target discovery:

```bash
gene-context stage1-context --results saved_stage1 --query-source-fasta original_query.fasta \
  --datasets targets.tsv --source-datasets shared_sources.tsv \
  --neighbor-cache cache/early_synteny --out stage1_with_sources
```

The original FASTA hash and exact record identifiers must match the saved search. Target sequences, BLAST scores, chains, locus IDs and coordinates are retained; source metadata and derived source-group ranking can be updated. This creates a new review checkpoint. Older decisions are preserved in their original checkpoint, rather than silently applied to new context.

When an exact RefSeq protein accession is available but the input lacks a GeneID, the standard context pass can retrieve that exact versioned NCBI record. It verifies the complete protein sequence and organism before accepting a unique GeneID or source symbol. An incorrect accession/sequence, ambiguous GeneID, species conflict or network failure remains explicit; no query sequence is replaced. Records and their SHA-256 receipts are cached and copied into the result. This is metadata retrieval, not a genome, model or protein-reference search. `--query-metadata-lookup off` prevents network retrieval and uses only verified existing cache records and supplied metadata. The same rule applies to every gene; old per-gene TSVs are optional provenance overrides, not a prerequisite.

A shared cache can register its source inventory at `NEIGHBOR_CACHE/source_datasets.tsv`. When `--source-datasets` is omitted, the standard stage-1/context pass reads that exact registered manifest automatically and snapshots its hash in the result. An explicit `--source-datasets` always takes precedence. This is a dataset inventory, shared across all genes; it contains no manually chosen source loci. Missing or conflicting datasets remain explicit.
