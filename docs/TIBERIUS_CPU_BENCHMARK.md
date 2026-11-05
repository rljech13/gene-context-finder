# Tiberius CPU resource benchmark — 2026-11-05

**Tiberius successfully ran locally without a GPU.** The two selected pretrained models are shared across the 18-species panel; no model needs to be downloaded per species.

## Models and disk space

The current panel assigns `vertebrates` to 16 non-mammalian vertebrates and `mammalia_nosofttmasking_v2` to *Homo sapiens* and *Monodelphis domestica*. This is a clade-based starting configuration, not proof of accuracy for each species, particularly lamprey, lungfish and axolotl. See the [complete 18-species mapping](benchmarks/tiberius_cpu_2026-11-05/species_model_mapping.tsv).

| Model | Downloaded archive | Extracted model files | Measured download time |
|---|---:|---:|---:|
| Vertebrates | 36,368,198 bytes (36.37 MB) | 39,211,722 bytes (39.21 MB) | 220.59 s |
| Mammalia, no softmasking | 29,805,917 bytes (29.81 MB) | 32,224,440 bytes (32.22 MB) | 180.38 s |
| **Total** | **66,174,115 bytes (66.17 MB)** | **71,436,162 bytes (71.44 MB)** | **3 min 41 s, downloaded concurrently** |

MB/GB here use decimal units. The measured Python environment occupied about **1.37 GB** of logical file content. The full working installation, including the upstream source checkout, both archives and extracted models, occupied about **1.68 GB**, excluding the installer cache. The retained `uv` cache added about **1.35 GB** of logical content; shared/COW filesystem blocks mean logical totals are not necessarily additional physical storage. These figures exclude source genome assemblies and their search indexes.

Source model configurations and weights: [Tiberius model definitions](https://github.com/Gaius-Augustus/Tiberius/tree/b3ea48f263a075260dbf1e91a412e82ccdc96d3c/model_cfg), [vertebrate archive](https://bioinf.uni-greifswald.de/bioinf/tiberius/models/vertebrates_weights.tar.gz), [mammalian archive](https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_nosm_weights_v2.tar.gz). The archives are downloaded once; repeated prediction uses their local extracted files. Locally recorded SHA-256 fingerprints are in the [benchmark data](benchmarks/tiberius_cpu_2026-11-05/results.json).

## CPU measurements

Hardware: Apple M3 Pro, 36 GiB system RAM, macOS arm64. Software: Python 3.12.9, Tiberius **2.0.9**, upstream commit `b3ea48f263a075260dbf1e91a412e82ccdc96d3c`, TensorFlow **2.20.0**, bricks2marble **0.1.3**, hidten **1.0.0**. GPU visibility was explicitly disabled and the visible-device list contained CPU only. TensorFlow intra-op threads were limited to four and inter-op threads to one. These are thread limits, not a claim that all four cores were continuously saturated.

| Run | Input | Wall time, including process/model startup | Peak RSS |
|---|---:|---:|---:|
| Vertebrates, one anole window | 110,928 bp | 26.24 s | 1.506 GB |
| Mammalia, one opossum window | 108,367 bp | 29.43 s | 1.384 GB |
| Vertebrates, 16-window workload | 1,580,609 bp | 234.01 s | 1.653 GB |
| Mammalia, two-window workload | 216,734 bp | 58.71 s | 1.476 GB |
| **18-window workload, two sequential model runs** | **1,797,343 bp** | **292.72 s = 4 min 53 s** | **1.653 GB maximum, not a sum** |

Timing used `/usr/bin/time -l`, an outer wall clock and process-tree memory samples every 250 ms. RAM is process RSS including resident shared libraries, not a recommended total machine-RAM specification. macOS reported a smaller private memory footprint; RSS is retained as the conservative resource measure. Existing system swap was already present and is not attributed to this benchmark.

The single-window and workload runs were separate predictor processes. The model's configured maximum chunk length was 400,050 bp; Tiberius automatically reduced it to a valid length close to the longest input window, as reported in its logs. `--batch_size 1` was used to limit peak memory. The native CPU environment installed in about 49 seconds on this connection. The first-ever TensorFlow import was separately observed at about 48 seconds; subsequent process startup was faster. Download and setup times depend on the connection and package cache.

## What the 18-window workload represents

The input comprises **seven real extracted windows from five species**: anole, green turtle, spotted gar, zebrafish and opossum. They came from the earlier NANOG/Pou5f3 target-search pilot and range from 65,610 to 113,327 bp, totalling 706,653 bp uniquely. To measure a workload with the panel's 16+2 model distribution, the six non-mammalian windows were repeated with unique sequence IDs until there were 16 sequences, and the opossum window was repeated twice. Thus the 18-window timing is a **measured computational workload using repeated real fragments**, not a completed analysis of 18 distinct species or 18 confirmed orthologous loci.

The two workload runs returned 55 + 8 coding models, including repeat observations by design. A separate Biopython audit, importing neither Tiberius nor genecontext, checked the GFF3 coordinates against each exact input, CDS exon phases, CDS FASTA roundtrips and protein translations. All **63 workload models** passed. Including the two single-window controls, all 70 exported models passed. The audit also checked CPU-only execution and downloaded model fingerprints. It does not validate biological gene accuracy.

These measurements cover **window gene prediction and its CDS/protein exports**. They exclude target discovery in complete genomes, BLAST identification, DNA/protein comparisons, repeat masking, gene-prediction accuracy evaluation and model training. The main pipeline still uses AUGUSTUS by default; this was a separate feasibility benchmark for a future Tiberius backend.

## Practical interpretation

A CPU-only workflow is feasible for local context prediction at this scale. Two shared weight sets occupy only 71 MB extracted; the larger disk cost is the TensorFlow environment. On this machine, a workload approximating one ~100-kb context per panel species took about five minutes.

For similarly sized inputs, simple proportional planning estimates are roughly ten minutes for 36 windows and fifty minutes for 180 windows, reusing one process per model. **Those larger workloads were not timed.** Longer loci, predictor reprocessing near chunk boundaries, processor speed and memory pressure can change the result. A benchmark on a machine with 8 GB total RAM, Linux or Windows was not performed.

## Reproduce the native CPU setup

The upstream `.[from_source]` extra currently requests `tensorflow[and-cuda]`, including indirectly through hidten's TensorFlow extra. On this Apple Silicon machine, both documented-extra resolution attempts failed because NVIDIA CUDA wheels are unavailable. A native CPU installation with the same core packages succeeded. The upstream prediction code and weights were not modified.

A standalone benchmark kit contains `cpu_launcher.py`, `run_benchmark.py`, `prepare_models.py`, `audit_outputs.py`, exact input FASTAs, mapping files, the tested dependency lock and raw results. It excludes installed environments, caches, upstream checkout and model binaries.

From a fresh extracted kit directory, with Python 3.12 available:

```bash
python3.12 -m venv venv
git clone https://github.com/Gaius-Augustus/Tiberius.git Tiberius_source
git -C Tiberius_source checkout b3ea48f263a075260dbf1e91a412e82ccdc96d3c
venv/bin/python -m pip install -r requirements.cpu.lock.txt
venv/bin/python -m pip install --no-deps ./Tiberius_source
venv/bin/python prepare_models.py
venv/bin/python run_benchmark.py --name repeat_vertebrates_16 \
  --model vertebrates --input inputs/workload18_vertebrates.fna
venv/bin/python run_benchmark.py --name repeat_mammalia_2 \
  --model mammalia --input inputs/workload18_mammalia.fna
venv/bin/python audit_outputs.py
```

The lock was tested on macOS arm64 and is not claimed to be a universal lock for every platform. Use fresh output names: benchmark runs are preserved and will not be overwritten. The launcher disables GPU and sets CPU thread limits, then calls the unmodified upstream launcher with local model files. Repeated runs do not need to redownload the weights.
