# Scientific execution

For scientific corpus work, read `$HOME/corpora/AGENTS.md` and supplement the
existing question card through the corpus writer.

## host2 stage 2

All stage-2 additions share the existing persistent server queue. Read
`docs/STAGE2_QUEUE.md` before submitting. Use the deployed `submit_stage2.py`
queue entry point from the server resource policy; do not start another
per-gene batch pool or direct stage-2 service. Keep its common CPU/RAM/swap
limits, FIFO order, immutable scientific source and manual STOP 2. Preserve
existing running calculations and deduplicate repeated submissions.
