# One study, arbitrary targets and assemblies

Gene Context Finder is a general locus/model comparison tool. Pluripotency is
one study using that tool. Gene names, species and reference panels are study
inputs, never dispatch rules in the program. Supported inputs are protein, CDS
or genomic-DNA queries and local genome/annotation manifests. Predictor models
must be appropriate for the organisms: the supplied Tiberius vertebrate and
mammal weights do not establish support for every taxon. Custom model paths and
per-genome model mappings use the native stage-2 options.

The primary workflow is:

| STOP | Computation | Human handoff |
| --- | --- | --- |
| 1 | Multi-reference search, loci, preliminary annotated synteny | Select loci, exclude or retain unresolved cases; optional family/copy groups |
| 2 | Reviewed windows, source annotation, independent prediction, protein support | Select exact target/alternative/neighbor models with reasons |
| 3 | Selected protein/CDS alignments, identity and coverage; optional genomic context | Final primary model and separate cases for each gene/genome |
| 4 | Exact sequence dossier, optional gene trees and imported exact structures | Interpret/export; no automatic biological approval |

Project commands call these existing stages. They never cross a manual STOP
automatically. Old stage commands and checkpoint formats remain usable.

## Configure and inspect

Create `targets.tsv` with columns `id`, `gene`, `query`, `query_type`, optional
`query_metadata`. Paths in this TSV are relative to the TSV. `id` is a safe,
unique software identifier; `gene` is the reference search-group name, which
can also name a family. The genome TSV remains the existing format: `genome_id`,
`fasta`, optional matching `gff3`, `protein`, `cds` and model settings. Its paths
are relative to that TSV. No species list is hardcoded.

The bundled POU5 control uses six small genomic windows, not whole genomes:

```sh
gene-context project init --name 'POU5 controls' --out results/pou5-study \
  --targets examples/studies/pou5.targets.tsv \
  --genomes examples/queries/pou5_test/controls.genomes.tsv
gene-context project validate --project results/pou5-study
gene-context project plan --project results/pou5-study --target all --stage stage1
gene-context project run --project results/pou5-study --target all --stage stage1
gene-context project report --project results/pou5-study
```

Open `project-report.html` for the gene × STOP overview and native report links.
`project status` prints the same status rows; `--json` is available. The overview
is an offline snapshot: regenerate it to refresh status. Viewing starts no tools.
`plan`, `status` and `validate` do not create computation outputs. Validation
checks the configured local inputs and fully verifies attached checkpoint seals;
it is not biological validation. Status uses compact identity pins, with full
scientific verification repeated at every handoff and execution.

`project.json` paths are relative to that file, independent of the current shell
directory. Stage options are native CLI destinations written with underscores:

```json
{
  "schema": "gene_context_project_v1",
  "name": "My study",
  "genomes": "../inputs/genomes.tsv",
  "targets": [
    {"id": "family_a", "gene": "FAMILY_A", "query": "../inputs/family_a.faa", "query_type": "protein"},
    {"id": "family_b", "gene": "FAMILY_B", "query": "../inputs/family_b.cds.fna", "query_type": "cds",
     "options": {"stage2": {"analysis_group": "reviewed_copy_b"}}}
  ],
  "options": {
    "stage1": {"threads": 4, "neighbor_cache": "../shared-cache"},
    "stage2": {"threads": 4, "protein_db": "../references/protein-db", "predictor_cache": "../shared-predictor"},
    "stage3": {"threads": 4, "genomic_context": false},
    "stage4": {"tree": "auto", "threads": 2}
  },
  "execution": {"stage2": {"mode": "plan"}}
}
```

Per-target options override project defaults. The project owns result paths and
decision bindings; options cannot redirect `out`, `results`, `decisions` or the
query/genome identity. To keep a shared family STOP1, configure each analysis
target with that family search name, attach the same STOP1, and specify its
reviewed `analysis_group` under stage2. Each downstream group has its own ledger.

## Import a review and run exactly one next stage

After reviewing STOP1, export its TSV or fully saved/scoped group JSON:

```sh
gene-context project decisions --project STUDY --target family_a --stage stage1 --file /path/to/loci.decisions.tsv
gene-context project plan --project STUDY --target family_a --stage stage2
gene-context project run --project STUDY --target family_a --stage stage2
```

Configure the stage2 executor first, as described below. After STOP2:

```sh
gene-context project decisions --project STUDY --target family_a --stage stage2 --file /path/to/stage2.decisions.json
gene-context project run --project STUDY --target family_a --stage stage3
```

After reviewing the final models in STOP3:

```sh
gene-context project decisions --project STUDY --target family_a --stage stage3 --file /path/to/stage3.final.decisions.json
gene-context project run --project STUDY --target family_a --stage stage4
```

An unreviewed STOP3 suggestion is retained as such and blocks the final handoff.
An explicitly excluded set of STOP1 loci is retained without launching STOP2.
Empty or contradictory model decisions, foreign review IDs, sequence mismatches
and changed seals are rejected before scientific execution.

## Keep existing work in place

Initialize a new compact control directory, then attach any completed native
checkpoint. A STOP2 or STOP3 can be attached independently when the historical
upstream directory is not locally available:

```sh
gene-context project attach --project STUDY --target family_a --stage stage2 --results /path/to/existing/STOP2
gene-context project decisions --project STUDY --target family_a --stage stage2 --file /path/to/exact/export.json
```

Attach fully verifies the seal. It references the original directory without
copying, resealing, rerendering, moving or editing it. When the upstream result
is also attached, the review binding is checked. Attaching an existing output
does not endorse its historical choices or assume it consumed today's export.
Reattaching the same exact result is idempotent. A distinct result requires
`--new-attempt`, which keeps the former result in history.

Project metadata lives only in `.gcf/state.json`, logs, managed run directories
and content-addressed decision objects. Every export's original bytes are saved
once by SHA-256. Changing a selected export makes dependent results `stale`;
their directories and old exports remain available. Rerunning a stage requires
`--new-attempt` and allocates a new output directory. An unchanged completed
managed stage is verified and reused. Failed/interrupted attempts are retained
and never retried implicitly. Native `stage2 --resume` remains available for
exact existing resumable runs; the project runner never silently restarts one.

Live jobs continue under their original launcher and pinned source. Wait for
completion and attach their sealed results. Run project commands on the host
owning the data. Do not restore archived trees or make local genome/result
copies merely to register them. Existing review IDs, decision registries, public
aliases, Sheets and the corpus remain authoritative; this ledger does not
rewrite them or automatically import browser localStorage.

## Stage2 execution and the existing shared queue

New projects default to `mode: plan`. This permits configuration and window
planning without starting a predictor. For a separate local compute machine,
explicitly select `--stage2-executor local` during init or set `mode: local`.
Compute execution uses Python plus the native tools on Linux/macOS; the
portable review app is also available for Windows.

On host2 all new stage2 jobs must use the existing persistent queue and resource
policy described in [STAGE2_QUEUE.md](STAGE2_QUEUE.md). Configure its existing
entrypoint, queue root, immutable deployed source and Python:

```json
{"execution":{"stage2":{
  "mode":"queue",
  "entrypoint":"/absolute/path/to/submit_stage2.py",
  "root":"/absolute/path/to/existing/queue",
  "source":"/absolute/path/to/immutable/deployed/source",
  "python":"/absolute/path/to/python3"
}}}
```

These are placeholders, not host2's actual locators. Project execution submits
one native pinned spec to that entrypoint. It never creates another queue,
dispatcher, predictor pool or shared-cache copy. Existing FIFO/resource limits
and deduplication remain in the queue. The deployed source must match the
runner's reviewed code; publish a new immutable source instead of changing the
source of live calculations. The returned ticket and its canonical
output are retained, including when the queue reuses an earlier job. Run
`project sync --project STUDY` after completion to verify and register the sealed
result; sync starts no computations. The downstream stage still needs review.

## Viewing and optional portable packages

Ordinary STOP2/3 rendering writes HTML plus `review_data.json`. It does not need
Go binaries and does not copy six launchers into every scientific result. The
hosted review service and saved scientific comparisons keep their normal roles.
For an independent portable review, create a separate package explicitly:

```sh
gene-context review-package --results /path/to/sealed/STOP3 \
  --out /path/to/new/review-package --platform darwin-arm64
```

`--platform current` is the default; `all` builds a universal distribution.
Other choices are `darwin-amd64`, `linux-amd64`, `linux-arm64`, `windows-amd64`
and `windows-arm64`. MAFFT availability still follows the portable runtime's
documented architecture support. Packaging never approves models, downloads
genomes or starts alignments. Existing scientific results and old portable
packages are left untouched.

## Optional STOP4 structures

Set `options.stage4.structure_provider` to `afdb` for existing exact-sequence
structures, or `af3-prepare` for native AF3 and AlphaFold Server request files.
`structure_limit` defaults to three unique proteins; `structure_cache` can point
to a shared local cache. The default remains offline with no automatic folding
or paid API job. See [the prototype guide](STAGE4.md#experimental-api-preparation-and-retrieval)
and [provider/access comparison](STRUCTURE_SERVICES.md).
