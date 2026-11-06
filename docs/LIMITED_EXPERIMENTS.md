# Run experiments with a hard shared RAM budget

`--threads` controls tool parallelism; it does not impose a memory limit. For Linux experiments, enforce memory limits on a cgroup containing the complete process tree. A virtual-address-space limit such as `ulimit -v` is not a physical RAM budget.

The host2 installation configured on 2026-11-06 has 256 logical CPUs and 1,081,802,870,784 bytes of physical RAM (about 1,008 GiB). Its shared experiment budget is **250 GiB = 268,435,456,000 bytes**, below one quarter of physical RAM. Swap is disabled for the group.

## Shared budget across simultaneous runs

The user-systemd slice is `gcfexperiments.slice`:

```ini
[Unit]
Description=Gene Context Finder experiments: shared hard memory budget

[Slice]
MemoryAccounting=yes
MemoryMax=250G
MemorySwapMax=0
```

Each experiment service belongs to this slice and also sets `MemoryMax=250G`, `MemorySwapMax=0` and `OOMPolicy=kill`. The slice limits the **sum** of simultaneous jobs and their descendants; separate per-job limits alone would not do that. Memory charged by the kernel, including charged page cache, counts towards the cgroup budget. An OOM event can terminate a job; it is not treated as a successful result.

This requires Linux cgroup v2 with a usable memory controller and a working user-systemd manager. The deployed launcher starts the slice, checks its properties and refuses submission if the expected policy is not active. A wrapper then checks the actual kernel `memory.max` and `memory.swap.max` files from inside the job before starting the requested command. User lingering is enabled so detached jobs survive SSH logout.

## Current host2 installation

```text
/mnt/data/i.rubinshteyn/gene_context_experiments/
  current -> source_v052_abdecd3/   immutable source snapshot
  submit.py                       detached, constrained submission
  run_job.py                      in-cgroup verification and run receipt
  resource_policy.json            exact shared budget and tool paths
  inputs/                         preserved manual-review checkpoints and panel
  cache/                          isolated Linux Tiberius runtime and shared references
  tools/                          pinned miniprot build
  runs/                           submissions, logs, start/completion receipts
  results/                        scientific outputs
```

Run from a shell on host2, choosing a fresh result directory:

```bash
GCF_ROOT=/mnt/data/i.rubinshteyn/gene_context_experiments
python3 "$GCF_ROOT/submit.py" --name panel-repeat --threads 32 -- \
  /usr/bin/python3 -m genecontext check-targets \
  --proteins "$GCF_ROOT/inputs/manual_review_v050_2026-11-06/02_models/model_options.proteins.faa" \
  --references "$GCF_ROOT/inputs/target_reference_check_2026-11-06/reference_panel.faa" \
  --reference-labels "$GCF_ROOT/inputs/target_reference_check_2026-11-06/reference_labels.tsv" \
  --out "$GCF_ROOT/results/panel-repeat" --threads 32
```

Submission prints the service name and log directory. Submission success is not computation success: inspect `completed.json` and the service's `Result`/`ExecMainStatus`. If OOM kills the wrapper, a completion file may be absent; query the service and its journal. The launcher's thread option records the intended count and sets thread environment variables; pass the tool's own `--threads` argument as shown. It does not impose a CPU quota.

All experiments for this installation must use this launcher. Directly launching a program outside the slice does not inherit its limit. Choose an appropriate budget on another machine; the 250 GiB value is specific to this host.

## Validation and data availability

A deliberately oversized 256 MiB allocation inside an isolated 64 MiB test service was killed by the kernel (`Result=oom-kill`, `ExecMainStatus=9`). Its kernel limits were read from inside the test process. This validates enforcement without approaching the production 250 GiB budget.

The saved coding-model reference check for all 18 assemblies has also run under the shared budget. The imported review checkpoints remain unapproved, with their evidence and pending decision files unchanged.

The transferred checkpoints include saved windows and sequence exports, not the complete original genomes. Comparisons from those saved windows can resume with the matching `--search-results` exports. New searches and extraction of other loci require the complete assembly inputs. Existing full-genome candidates on host2 are inventoried separately; they must be checked against assembly/source fingerprints before reuse. The host1 transfer attempt timed out, so no complete 18-genome migration is claimed.
