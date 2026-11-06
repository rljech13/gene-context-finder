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

Use `--datasets datasets.tsv` when the complete annotation/protein inventory is separate from the target genome list. Use `--source-datasets references.tsv` for additional source-reference assemblies that are not study targets. The program resolves query origins from identifiers and exact sequences in the supplied data. It preserves ambiguous matches and reports unresolved origins. No manually constructed source-locus list is required. `--query-metadata queries.tsv` can supply source species, gene symbols and accession provenance; if supplied, it must describe every input FASTA ID exactly once.

## What runs

1. The default `--discovery auto` uses BLASTP for protein/CDS queries when a GFF is supplied. **All usable annotated isoforms** enter the search; identical full protein sequences are searched once and keep every transcript/gene binding. CDS queries are translated first. `--discovery blastp` requires annotation, while `--discovery genome` explicitly requests whole-genome TBLASTN. In auto mode a genome without GFF uses TBLASTN. DNA queries always use BLASTN.
2. BLASTP HSPs are chained within each query/subject protein. Exact product → transcript → gene relationships map them to genomic loci; all matched isoforms of one annotated gene share one locus. Separate annotated genes stay separate, even if their proteins are identical. Genomic searches instead chain compatible query/genome HSPs and group overlapping aligned blocks. No gene-name lookup substitutes for sequence search.
3. Loci are ranked within each genome by the strongest single-query chain score. Overlapping query positions are counted once; scores from separate isoforms are not added. All candidates and the raw BLAST evidence remain available.
4. Available reference and candidate neighborhoods are checked automatically. Annotation and proteome indexes are cached. Reference neighbor proteins are matched against target proteomes in batches; coordinates, order, side and direction provide preliminary context evidence. Target genes are not counted as their own neighbor evidence.
5. The report stops for human selection. No Tiberius prediction, target-protein MSA or comparative CDS alignment runs here.

A BLASTP no-hit or weak result does not silently trigger a potentially expensive whole-genome search. Inspect the result and repeat the relevant genome with `--discovery genome` if annotation rescue is needed. Invalid/mismatched annotation is an explicit error. Missing annotation products and identical proteins at multiple gene loci are reported.

The new default is `--refine none`. Annotated candidates use GFF gene boundaries; local miniprot refinement applies only to genomic searches. Optional `--refine miniprot` performs local spliced boundary refinement; these predictions stay in technical details and do not establish gene identity. Older imported checkpoints may already contain those predictions.

Early contexts start with five genes on each side, expanding to ten/twenty when insufficient preliminary support is found, within 1 Mb. These are configurable limits (`--synteny-neighbors`, `--synteny-max-neighbors`, `--synteny-max-distance`), not biological thresholds. The actual inspected span is shown; it is separate from the proposed ±50 kb window for detailed downstream work. `--neighbor-cache` selects a reusable cache. By default it is `_genecontext_cache/early_synteny` beside the output, not an implicit home-directory download cache. No reference/model downloads are performed by this pass. `--early-synteny off` requests discovery alone.

## Read the report

Open `report.html`. Choose a genome or open the overview. Each locus is one row; the strongest retained candidate is displayed first, with a button to show the next ten. Selected loci remain visible. Nothing is preselected. Each card is one candidate locus and immediately lists **every supplied query variant**, including weak/no-match variants, with per-query covered positions/percentage, aligned identity and chain score. Each row uses that query's strongest product chain, not the locus-wide best query or a sum across variants. Identical sequence variants share evidence; their count does not represent independent support.

**Where this locus matches each query** shows:

- For BLASTP: annotated CDS blocks on a genomic axis, and separate query/subject amino-acid match axes. CDS blocks are annotation features, not direct BLASTP-to-DNA alignments.
- For genomic BLAST: aligned DNA blocks on a shared genomic axis.
- Every supplied query variant, including variants without retained support at this locus.
- Covered query positions as `N / full query length`, coverage, HSP-weighted identity, chain score and best HSP E-value.
- A separate position map along each query and the underlying HSP coordinates/metrics.

These are maps of the existing search matches. They do not require new comparisons between candidate proteins. BLASTP also reports coverage of the annotated subject protein, without assuming the annotation is complete. Genomic BLAST does not estimate coverage of an unknown complete candidate gene. The miniprot models, if present, are optional technical evidence within their locus.

**Neighbor genes** shows the candidate once, followed by **every available query-source reference locus** in the same diagram. Query variants/isoforms from one annotated source gene share one reference neighborhood. Each row follows its own gene order, oriented with the target pointing right. Matching colors propose correspondence to the same candidate neighbor; they do not imply an orthogroup. The default is five saved neighbors per side with equal schematic spacing. Each reference row states how many proposed pairs are visible versus retained in the full saved comparison. Same-assembly/self-control rows are labeled and do not receive cross-species support colors.

The caption also states how many query variants lack a resolved source neighborhood. **How this is checked, missing references and pair details** explains the procedure, lists availability for every query, and exposes optional scope/scale controls and individual-reference inspection. The latter retains protein-match metrics, ambiguity, gene choices and the two-row diagram. The backend already compares every resolved reference against every candidate; the report does not select just one source for computation.

A query sequence alone cannot provide its genomic neighbors. The supplied source GFF/proteins/CDS must resolve its original locus. Missing source data and identifier/sequence mismatches are reported, not replaced by invented references. Source resolution currently requires identifier plus matching sequence or a permitted exact-sequence binding. The saved NANOG pilot therefore has two accepted reference loci (human and zebrafish), two chicken variants with identifier/sequence mismatches, and one mouse variant without a resolved source. This is a source-data limitation, not a claim of missing synteny.

Sequence tables, metric explanations, neighbor selection, matching diagnostics and decision controls are collapsed by default. See [metric definitions and worked examples](STAGE1_METRICS.md). Neighbor choices persist in draft JSON and export separately as `neighbors.to_inspect.tsv`; they request follow-up inspection and are not consumed by `stage1-select` yet. No heuristic score is shown as an orthology verdict. Missing annotation, missing proteins, unresolved source loci and insufficient proposed pairs remain distinct evidence states. Lack of context support never automatically removes a locus.

Saved display thresholds are coverage 0.20, identity 0.25, E-value 1e-5. One chain must pass all three together; it need not be the highest-scoring chain shown in the locus summary. The strongest candidate remains visible even if below these gates. Further alternatives are available through **Show next candidates**. Check **Hide below-threshold candidates** to narrow the view. Selections remain visible even when below those filters. Browser filters do not delete evidence or change the saved checkpoint.

`loci.tsv` contains every ranked candidate; `passing_loci.tsv` is the saved passing subset; `query_support.tsv` contains per-query/product chains, covered positions and separate subject-protein or genomic intervals. BLASTP loci include original gene/model identifiers and all matched isoforms; `raw/<genome>.annotation/` contains product mappings, matched proteins and cache/provenance receipts. Annotated locus DNA extraction is deferred until selection, so `loci.fna` contains only any genomic-search results; it may be empty for a purely BLASTP run. Genome coordinates are zero-based half-open. Query positions in the report's HSP tables are one-based inclusive; TSV interval fields are zero-based half-open. Chain identity is an estimate from local HSP identities, with overlaps prorated. Chain E-value is the best HSP E-value, not a combined gene-level probability. Ranking and context evidence do not prove orthology.

## Choose loci

For each gene/genome pair, check one or more loci and choose **Investigate selected loci**, or choose **Exclude** / **Unresolved** with no locus selections. Select a short reason; a comment is optional except for **Other**. Save the pair decision. Different CDS predictions stay inside the locus; this stop chooses genomic regions, not a definitive transcript.

The report supports a portable draft JSON and automatic browser-local saving. A draft is not an approved decision file. **Validate all pairs** lists every pending or inconsistent pair and previews the selected windows. **Validate and download decisions TSV** exports exact IDs only after all pair decisions are saved and valid. For an inspected pair, unselected loci are explicitly recorded as excluded from this run. No original evidence or decision file is overwritten.

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

The locus-first selected-window export is not yet connected to the legacy `contexts` model-first interface. Detailed prediction and final conservation comparisons remain separate downstream work.

## Refresh the presentation only

```bash
gene-context stage1-report --results results/nanog_review \
  --out results/nanog_compact
```

This verifies the source checkpoint and copies its sealed artifacts to a separate directory. It preserves the exact scientific evidence, filters, locus IDs and review ID, so existing drafts and decision TSVs remain compatible. Only the report and presentation provenance change; no BLAST or context calculation is run. The source is untouched.
