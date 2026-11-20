# Studies with arbitrary targets and assemblies

A study is a `project.json` file plus its `.gcf/` ledger. Inputs stay read-only.
References can represent a gene, a family or a reviewed copy group. Their
labels do not select different code paths.

```json
{
  "schema": "gene_context_project_v1",
  "workflow": "sequential_v2",
  "name": "My families",
  "genomes": "../inputs/genomes.tsv",
  "targets": [
    {"id": "family_a", "gene": "FAMILY_A", "query": "../inputs/a.faa", "query_type": "protein"},
    {"id": "family_b", "gene": "FAMILY_B", "query": "../inputs/b.cds.fna", "query_type": "cds"}
  ],
  "options": {"stage1": {"threads": 4}},
  "execution": {"stage2": {"mode": "plan"}}
}
```

Paths are relative to `project.json`. Genome-manifest paths are relative to the
TSV. Query types are `protein`, `cds` and `dna`. Options use native CLI destination
names with underscores; unknown options fail validation. Input/output and
handoff paths are owned by the project, never arbitrary option overrides.
Global options can be overridden within a target's `options` object.

## Commands

| Command | Effect |
| --- | --- |
| `project init` | Create a new sequential configuration from target/genome TSV files; start no calculations. |
| `project validate` | Validate configured local inputs/options and attached seals. |
| `project plan --target ID --stage stage1` | Plan discovery without running it. |
| `project run --target ID --stage stage1` | Run discovery, then stop for locus review. |
| `project decisions --target ID --stage stage1 --file EXPORT` | Store exact exported locus decisions after native validation. |
| `project plan/run --target ID --stage stage2` | Plan/run annotation and Tiberius preparation after locus selection. |
| `project attach --target ID --stage stage2 --results CHECKPOINT` | Reference an existing sealed preparation checkpoint in place. |
| `project review --target ID` | Create the current STOP2/STOP3 workspace from prepared models; run no prediction/alignment/folding. |
| `project handoff --target ID --file STOP4-input.json` | Validate and retain the exact final STOP3 export. |
| `project plan/run --target ID --stage stage4` | Build the final dossier from the accepted handoff, reusing its MSA/tree. |
| `project status/report` | Read state or refresh an offline target × STOP overview. |

Every command above except `init` requires `--project PATH`. Plan/run also accept
`--target all`. `project run --stage stage3` is deliberately rejected for new
projects: STOP3 is the synteny page, not the legacy alignment computation.
Preparation's `stage2` CLI name is retained as a stable engine/checkpoint API.

The model/alignment/tree page and separate synteny view preserve source comments
and warnings. Final outcomes and post-synteny reference approval are distinct.
A candidate checkbox never authorizes folding. Export preparation never submits
a Server job. Cached alignments/trees are reused only for the bound exact input.

## Existing data and revised decisions

`attach` references the original sealed checkpoint; `review` writes only a new
workspace under the study's `.gcf/`. Decisions and handoffs are content-addressed
objects. The same valid review/result is reused; `--new-attempt` explicitly
creates another while keeping history. Old review IDs and source choices are
not overwritten or inherited as new approvals.

Handoff validation binds the source and alignment review IDs, all member
sequences/metadata, model decisions, genomic regions and neighbor
correspondences to the prepared workspace. The saved MSA/tree checksums and
sequence identities are checked. A changed review or selected handoff makes
an earlier final result stale. `project handoff` is repeated for the new export;
existing finals are retained and rebuilt only with `--new-attempt`.

Configurations without `workflow` retain [legacy behavior](PROJECTS_LEGACY.md).
Create a new sequential study and attach their existing preparation checkpoints
when adopting the current flow; do not relabel a live legacy ledger in place.

## Execution and dependencies

Preparation is plan-only by default. Set `execution.stage2.mode` to `local` for
explicit local computation, or configure `queue` with the existing deployment's
`entrypoint`, `root`, `source` and optional `python`. A queued source must match
the reviewed package. The tool creates no competing dispatcher. Live server
jobs and their pinned source stay untouched.

The study runner currently targets Linux/macOS; Windows users can run analysis
under WSL. Optional portable review launchers also support Windows. External
search/prediction tools and appropriate organism models are separate dependencies.

`examples/minimal/` contains two unrelated synthetic targets, two assemblies,
matching annotation and both protein/CDS queries. It is a discovery software
control. `tests/test_project_review.py` separately covers the complete current
prepared-models → review → final dossier route, without pretending to run a
predictor on artificial genomes or assert a biological result.
