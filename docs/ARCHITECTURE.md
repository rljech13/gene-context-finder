# Architecture and preservation boundaries

The project configuration describes inputs and execution policy. It contains
no special dispatch rules for a particular gene family or species list.

| Layer | Modules | Responsibility |
| --- | --- | --- |
| CLI | `cli.py` | Explicit commands and options; no work on import. |
| Study | `project.py`, `project_run.py`, `project_review.py`, `project_report.py` | Input paths, exact decision objects, review workspaces, execution, history and stale-result checks. |
| Scientific checkpoints | `stage1*`, `stage2*`, `stage4_handoff.py` | Discovery, preparation and final export. Existing scientific engines retain their checkpoint schemas. |
| Current review | `alignment_review.py`, current `*_review*.txt`, `stop4_handoff.txt` | One current model/alignment/tree → synteny → final-outcome flow. |
| Optional local runtime | `portable-review/` | Authenticated localhost app, explicit MAFFT/IQ-TREE/trimAl work, cached results and persistent choices. Built separately. |
| Compatibility | `legacy_review.py`, `vendor/legacy_review/` | Frozen inputs/templates needed to reproduce completed-study migration. They do not define the current normal workflow. |
| Structures | `alphafold_db.py`, `af3_*`, `synteny_reference.py`, `stage4_structure.py` | Exact-sequence retrieval/import and approval-gated preparation; no automatic Server submission. |

New projects declare `workflow: sequential_v2`; configurations without that
field retain `legacy_v1`. A new workflow does not rewrite an existing ledger.
Model preparation remains the native `stage2` engine; numeric legacy command
names and current UI stop labels are deliberately documented separately.

The current generic route is covered by a source-preserving integration test:
attach prepared models → create review → validate exact handoff → run STOP4 →
reuse the same result. A new review or handoff preserves history and invalidates
the previous final result. Native sequence/context validation is repeated when
accepting a handoff and before final execution.

`docs/report_source_snapshot.json` records the source hashes imported from the
implementation used by current reports and the final adapted source hashes.
The working research checkout and its live services are never deployment targets
for repository consolidation. New package builds include the reviewed sources;
existing research jobs keep their own pinned implementation.

Core installation does not include predictor weights, genomes, experimental
results or precompiled portable launchers. Optional launchers are reproducibly
built from `portable-review/`; the build receipt records compiler, source and
binary hashes. Native scientific tools remain explicit external dependencies.
