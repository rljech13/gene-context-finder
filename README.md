# Gene Context Finder

Find, review and compare gene loci across user-supplied genomes and reference panels.

The primary workflow is **configure study → STOP1 loci → STOP2 models → STOP3 sequence comparison and final selection → STOP4 dossier**. Each transition requires the exact exported review decisions. Pluripotency is one application; gene names, families, species and references are configuration inputs. The [project guide](docs/PROJECTS.md) covers multi-target studies, existing checkpoint adoption, safe review history, shared-queue execution and optional portable apps. Native single-stage commands below remain available.

After installation, try a project using the bundled small POU5 control windows:

```sh
gene-context project init --name 'POU5 controls' --out results/pou5-study \
  --targets examples/studies/pou5.targets.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv
gene-context project validate --project results/pou5-study
gene-context project run --project results/pou5-study --target all --stage stage1
gene-context project report --project results/pou5-study
```

Open `project-report.html` in the study directory. Review STOP1, import its export with `project decisions`, then explicitly run the next stage. New projects default to planning-only STOP2 until a compute executor is configured. The six control windows are not whole-genome validation. Existing calculations continue on their pinned source; use `project attach` on the owning host after they complete, without copying or editing their data.

The optional [STOP 4 prototype](docs/STAGE4.md) turns a sealed STOP 3 result into
a compact candidate dossier with separate identity/model-quality assessments,
per-group ML gene trees and optional exact-sequence structural comparisons.
Unreviewed final suggestions require explicit `--preview`; no folding jobs or
queues are started automatically.

The first review stop is **stage 1**: search all supplied genomes, assemble local matches into loci, automatically check available reference-neighbor contexts, and stop for manual selection. Multiple transcript variants and source species support one target/search group. The report ranks loci, maps matches for every query and supports selection by checkbox. Standard FASTA source headers are read automatically for every gene; a common source-dataset inventory can be reused across query panels. Missing source contexts remain explicit, and duplicate self comparisons are hidden by default.

## Install

```bash
git clone https://github.com/rljech13/gene-context-finder.git
cd gene-context-finder
conda env create -f environment.yml
conda activate gene-context
```

Stage 1 needs BLAST+ and samtools, already listed in `environment.yml`. Miniprot is optional. Genome annotation and annotated proteins/CDS enable the fast early-synteny pass; missing annotation is reported explicitly. This stage does not need Tiberius, model weights or Swiss-Prot.

## Try the POU5 panel

The included [test panel](examples/queries/pou5_test/README.md) contains 16 annotated transcript/CDS variants from human, zebrafish, chicken and Xenopus tropicalis. Source gene symbols and accession versions are retained. The POU5 group includes separately labelled family members and is not a one-to-one orthology assertion.

```bash
gene-context stage1 --gene POU5 \
  --query examples/queries/pou5_test/proteins.faa --query-type protein \
  --query-metadata examples/queries/pou5_test/queries.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv \
  --out results/pou5_stage1 --threads 4
```

This searches six small reference windows from four species, **not complete genomes**. To search the study assemblies, replace `--genomes` with their explicit TSV (`genome_id`, `fasta`; relative or absolute paths, gzip accepted).

With annotated datasets the default search is BLASTP over all available protein isoforms, grouped by their annotated gene loci. `--discovery genome` requests whole-genome TBLASTN for rescue. DNA queries use BLASTN.

Stage 1 shares prepared search inputs across query panels. Reuse the same
`--neighbor-cache` directory for different genes: annotation-product BLAST databases,
expanded genomic FASTA/faidx and nucleotide BLAST databases are built once per exact
input/tool version. Without this option, serial and parallel runs use the same
default cache beside the result directory. Cache publication is locked and atomic;
all cached output hashes are checked before reuse. Corrupted entries are reported
and preserved. `genomes.summary.tsv` records reuse and database preparation time.
This changes preparation only; search parameters, locus IDs and manual review
stops retain their existing meaning. See [the optimization plan](docs/OPTIMIZATION_PLAN.md).

Open `results/pou5_stage1/report.html`. Inspect the ranked loci, query tracks and available early-synteny evidence. Select loci, choose a short reason and save each pair decision. Validate and download the decision TSV. The next stage does not start automatically. The small control-window example does not include full annotation/proteome datasets, so its early context is explicitly unavailable; for study assemblies, include `gff3` and `protein` or `cds` in the genome TSV. Additional reference datasets can be supplied with `--source-datasets`.

```bash
gene-context stage1-select --results results/pou5_stage1 \
  --decisions /path/to/downloaded/loci.decisions.tsv \
  --out results/pou5_selected
```

The output contains selected loci and windows in TSV/BED, a validation summary and the decision snapshot. It does not start gene prediction. The [stage 1 guide](docs/STAGE1.md) explains automatic context preparation, source-query resolution, caching, the review screen and reuse of a saved search through `stage1-context`.

On host2, use the [constrained launcher](docs/LIMITED_EXPERIMENTS.md) for every experiment. Its shared hard limit is 250 GiB RAM with swap disabled; use absolute `/mnt/data` paths. The installed source snapshot must contain v0.6.0 or later.

## Reannotate selected windows: stage 2

Prepare a broad local protein reference once, then run the selected windows through Tiberius and BLASTP:

The default stage-2 planner first refreshes ten coding neighbors per side from the matching full GFF3, without a distance cutoff. Run `stage2-plan` with `--genomes` to inspect the complete intervals and captured gene counts. Missing/short neighborhoods and oversized windows block prediction; dry planning retains their full spans for review. See [window planning](docs/STAGE2.md#reannotation-bounds-from-preliminary-synteny).

```bash
gene-context prepare-stage2-db --download-vertebrates --out /mnt/project/references/vertebrates_reviewed
gene-context setup-predictor --predictor-cache /mnt/project/cache/tiberius --models vertebrates mammalia
gene-context stage2 --results results/pou5_stage1 --selection results/pou5_selected \
  --genomes /mnt/project/genomes.tsv --protein-db /mnt/project/references/vertebrates_reviewed \
  --predictor-cache /mnt/project/cache/tiberius --out results/pou5_stage2 --threads 8
```

Open `results/pou5_stage2/report.html` for **STOP 2**: compare provided annotation and DNA-only predictions, inspect broad-database naming alternatives and target-query support, and choose target models/neighbors. No model is pre-approved. Overlapping selected windows share one bounded prediction interval. The downloadable draft/decisions preserve exact stage 2 IDs. Full details, draft-input rules and host2 memory limits are in [the stage 2 guide](docs/STAGE2.md). Analysis is offline; reference/model preparation is explicit.

Interrupted stage 2 runs can continue with the **same command plus `--resume`**.
Completed genome extraction, verified per-window predictions and completed protein
searches are reused. Use the same `--predictor-cache` across genes for prediction
reuse. Inputs, decisions, parameters and method identity must match; completed
review choices are retained. `--tiberius-worker-mode reuse` optionally keeps one
neural network in memory while each window gets independent chunking and a fresh
HMM. The default remains one process per uncached window.

For a smaller standalone HTML, keep the original and create a sibling file:

```bash
gene-context compact-report --report results/pou5_stage2/report.html \
  --out results/pou5_stage2/report.compact.html
```

This preserves the application script, coordinates, labels and zoom. It reduces
file size; it does not reduce the hydrated evidence or demonstrate lower browser
memory use. The compact file uses native browser gzip decompression and offers a
link to the original if that capability is unavailable.

At STOP 1, **Groups for stage 2 — families or copies** lets you assign selected loci to named sets across assemblies. Download each group's JSON and submit it as `stage2 --draft` to its own output directory. Full grouped drafts and TSVs require `--analysis-group ID_OR_NAME`. This reuses the saved BLAST evidence and keeps families in separate STOP-2 reports; see [STOP-1 grouping](docs/STAGE1.md#separate-familycopy-groups-after-stop-1).

After saving the review, download `stage2.decisions.json` and run:

```bash
gene-context stage3 --results results/pou5_stage2 \
  --decisions /path/to/stage2.decisions.json \
  --out results/pou5_comparison --threads 8
```

This reuses the selected sequences, runs protein/codon and local genomic comparisons, and writes a traceable sequence catalogue. It does not repeat search or prediction. Use `--validate-only` first to check the handoff without starting tools. See [stage 3](docs/STAGE3.md) for outputs and interpretation limits.

The report ends with two tables: at most one target representative per gene/genome, and separate cases for additional copies or uncertain assignments. Initial choices are explicit suggestions with source identity, CDS quality and reference identity/coverage evidence. You can change the representative, select none, and record the reason. Download split package exports the current tables and exact CDS/protein FASTAs. To apply exported choices to a separate full report without rerunning alignments:

```bash
gene-context stage3-finalize --results results/pou5_comparison \
  --decisions /path/to/stage3.final.decisions.json \
  --out results/pou5_final
```

This final assignment retains the original STOP-2 selection and all-by-all comparisons; it does not resolve orthology automatically.

In the alignment-first route, STOP2 covers models → alignment → tree →
confirmed-set handoff. **STOP3 is the separate synteny/context review.** At its
end, **«Утвердить модель как эталон»** enables preparation of a Server JSON for
that exact protein. The legacy comparison table has no reference-approval
button. Saved choices remain intact; old reports retain their legacy route to
STOP4. [Post-synteny approval](docs/SYNTENY_REFERENCE.md).

## Earlier pipeline

The earlier `search → contexts → compare-reviewed` workflow is preserved and documented in [the legacy workflow](docs/LEGACY_WORKFLOW.md). Its model-first checkpoints differ from the new locus-first stage. Use the new `stage2` command for locus-first selections. Legacy `contexts` and `compare-reviewed` retain their different checkpoint format; they are not consumers of the new STOP-2 decisions.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests use local fixtures. Tool-dependent integration tests skip when executables are unavailable. The POU5 controls check annotated reference windows and neighbouring frog copies; they do not establish all-genome discovery accuracy.

Refresh an existing review into the compact layout without rerunning any search:

```bash
gene-context stage1-report --results results/nanog_stage1 --out results/nanog_compact
```

The original evidence and review/decision IDs are retained. The default view shows the first ten candidate loci in combined review order, each with a compact support table for every query variant. The score combines coverage- and identity-weighted BLAST support across source groups; contributions are available in Details and review_ranking.tsv. Synteny and technical details are collapsed; opening Synteny shows all available query-source neighborhoods. Curated GeneIDs can bind source neighborhoods across isoform differences while retaining sequence-discrepancy flags. See [metric definitions](docs/STAGE1_METRICS.md).

Stage 1 retains at most **20 candidate loci per gene in each assembly** before preliminary synteny and manual review. Each assembly keeps its best loci by the existing weighted BLAST review ranking (bit score, query coverage and aligned identity across supplied source groups); synteny does not change this ordering. This is a review-pool limit, not a biological rejection of omitted hits. Raw search evidence is retained.

Use `--max-loci-per-assembly N` (1–20) with `stage1` for a smaller pool. To restrict a completed checkpoint without another search or synteny calculation:

```sh
gene-context stage1-limit --results original_checkpoint --out review_top20
```

The new checkpoint has its own review ID and blank decisions, retains the original locus IDs and full evidence for each kept locus, and points to the verified original checkpoint for raw BLAST files. Keep that original directory available. Old decisions belong to the original review and are not automatically approved for the restricted pool.
