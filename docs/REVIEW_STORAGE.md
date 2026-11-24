# Review files and project entry point

The current vertebrate project has one compact control directory:
`outputs/gene_context` in the ChatGPT project workspace. It is also linked from
`~/Downloads/Gene_context` and served at
<http://127.0.0.1:8794/project/START_HERE.html>.

This directory is the entry point for review decisions and report links. It is
not a copy of scientific working directories. Genomes, live jobs and sealed
results remain on host2. Historical output folders keep their original paths.

Review identity is `gene/group + review_id + SHA-256`. A download suffix or newer
mtime is insufficient to select a version. Store each byte-distinct decision once
under `decisions/objects`; keep every original locator in `registry/decisions.json`.
`decisions/current` links to the exact versions consumed by verified runs.

To retain a new export without starting a calculation:

```sh
python3 tools/register_review.py /path/to/stage2.decisions.json
```

Run this from the compact control directory. The helper deduplicates exact bytes,
records unmatched review IDs explicitly, and does not silently replace current
decisions. Checkpoint/model validation remains a separate required handoff.

The canonical tracker is the control directory's `tools/collect_tracker.py` on
host2. It applies `registry/deferred_cases.json`; MYCa, MYCb, MYCh and MYCL are
deferred by the user for a separate investigation, not waiting for Stage 2.
The companion sheet planner preserves custom final decisions, reviewer notes,
stable record IDs, formulas and current public Stage 2 viewer links.

Report renderer/presentation updates are independent of scientific checkpoints.
Stage 3 exports the current chosen models and comments inside each individual
gene report; the old global default export remains historical.

Stage 2 uses the registered `presentation_methods/stop2_ui_20261118` renderer.
The canonical project helper `tools/stop2_presentation.py` builds a separate
presentation and portable ZIP from the currently served payload. Its publish
phase requires `presentations/stop2_ui_20261118/ui-validation.json` to match the
exact report SHA for every delivered gene. It updates existing public aliases
under the queue lock, never scientific results or current decision aliases.
Version marker `gcf-stop2-ui-20261118-v1` and exact renderer hashes are recorded
in `presentation.json`; the package manifest verifies every ZIP member.
