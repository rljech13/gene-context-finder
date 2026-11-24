# Shared host2 stage-2 queue

All additions to stage 2 use the same persistent queue. Do not launch another
`run_batch.py`, per-panel pool, or direct stage-2 service. Prepare and validate
the exact checkpoint/form/window plan, then enqueue its pinned command through
the existing server `submit_stage2.py` entry point. The entry point accepts
`enqueue --root QUEUE_ROOT --spec prepared-job.json`; `status` reads the same
queue. The queue root and deployed implementation are recorded in the server's
resource policy. This applies across genes and chats.

A prepared spec contains `panel`, `name`, `job_id`, `loci`, `locus_ids`,
`source_review_id`, optional `group_id`, `source`, `source_files` (SHA-256), and
`command` (an argument array beginning with Python, `-m genecontext stage2`).
The command includes the verified input draft, matching genomes, model and
protein catalogues, reviewed full-window budget, output directory and at most
32 threads. The queue derives input fingerprints, preserves FIFO order, and
returns an existing ticket for a repeated submission. Never edit sealed
scientific checkpoints or a deployed immutable source to add a job.

The shared dispatcher admits at most eight analyses and 256 threads in total.
This live policy is stored in queue state and resource policy; the conservative
code defaults remain two analyses / 64 threads for a newly initialized queue.
The existing aggregate hard limit is 250 GiB with zero swap. A 200 GiB
`memory.high` threshold slows allocations/reclaims cache before the hard limit.
Admission requires measured nonreclaimable memory plus 32 GiB for the
starting job and 16 GiB growth reserve for each active job to stay within
192 GiB. Starts/resumes are spaced by at least 30 seconds. File cache is reported separately; dirty/writeback pages are included
conservatively in the reserve calculation. At 128 GiB of nonreclaimable memory,
active job trees are paused. An increased OOM-kill counter holds the queue for
operator investigation; it does not silently retry failed scientific work.
These controls reduce OOM risk; they are not a proof of every predictor's peak.

Live legacy calculations are adopted by PID plus Linux process start time.
Paused calculations retain their process state and are resumed instead of
restarted. Previous dispatchers stay frozen until all their adopted children
finish, then are retired without executing their old pending pools. A stopped
or lost process is never replaced silently. New jobs run in separate bounded
services through the standard submission launcher; a dispatcher restart does
not terminate those services. Durable process/exit receipts distinguish an
interrupted dispatcher from a completed scientific command.

Results retain separate groups, checkpoints, model review IDs and ZIPs. Queue
completion does not approve models or start stage 3. The unified dashboard
links the original results and records missing/short predicted neighborhoods.

Presentation release is a separate gate. A scientific job may pin an older
renderer; successful completion and HTTP 200 do not establish that its public
HTML is current. Before delivering a new Stage 2, run the registered presentation
method against its current public payload, retain its review ID and exact
payload/decision bytes, validate the generated report in the UI, then atomically
switch the existing report and ZIP aliases. Record the renderer version and
hashes separately from the scientific source. The existing monitor applies this
gate to Q0029–Q0033 without restarting their pinned calculations. See
`docs/REVIEW_STORAGE.md` for the canonical method and validation record.

Use the existing shared predictor runtime, weights and protein/reference
catalogues. Queue parallelism must not create per-job model, runtime or catalogue
copies. Existing scientific commands pin their shared caches; retain them.
