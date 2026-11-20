# Gene Context Finder

A sequence-driven tool for finding and comparing candidate gene models across
assemblies. Genes, families, reference sequences and species are study inputs.
The pluripotency analyses are one application of the tool.

Supply protein, CDS or genomic DNA queries and local genome/annotation files.
The tool keeps exact sequences, evidence, comments and manual decisions tied to
their source checkpoints. It does not infer orthology or functional equivalence
from a name, tree or predicted structure alone.

## Workflow

| Step | What happens |
| --- | --- |
| STOP 1 | Search candidate loci and inspect preliminary context; explicitly select regions. |
| Preparation | Retain supplied annotation and run independent Tiberius prediction in selected windows. |
| STOP 2 | Choose exact models, inspect MAFFT protein/codon alignment, build and inspect the gene tree. |
| STOP 3 | Review genomic context and synteny, record final model outcomes, and export a verified handoff. |
| STOP 4 | Build the final sequence/context/tree dossier from that handoff, with optional exact structures. |

Reference-model approval is a separate action at the end of STOP3. Only that
explicit approval enables preparation of an AlphaFold Server upload for the
exact protein. Nothing is submitted to a folding service automatically.

## Install

Python 3.10+; the core Python package has no mandatory Python dependencies:

```sh
python -m pip install .
gene-context --help
```

BLAST+ and samtools are needed for discovery. Model preparation additionally
needs Tiberius and a local naming database. MAFFT and IQ-TREE are used during
explicit alignment/tree calculations. Existing exact structures can be imported;
`pip install '.[structures]'` adds the optional mmCIF parser.
See [setup and input preparation](docs/SETUP_AND_DOWNLOAD.md).

Gene names and species are unrestricted labels. Prediction quality depends on
appropriate models and data: supplied vertebrate/mammal predictor weights do
not establish support for every taxon. Custom predictor models can be configured.

## Start a study

Create a genome TSV (`genome_id`, `fasta`, optional matching `gff3`, `protein`,
`cds`) and a target TSV (`id`, `gene`, `query`, `query_type`). File paths are
relative to their TSV. Each target can have multiple reference FASTA records.

```sh
gene-context project init --name 'My study' --out work/my-study \
  --genomes inputs/genomes.tsv --targets inputs/targets.tsv
gene-context project validate --project work/my-study
gene-context project plan --project work/my-study --target all --stage stage1
gene-context project run --project work/my-study --target all --stage stage1
```

Review STOP1, import its explicit exported choices with `project decisions`,
then plan/run preparation. New projects default to planning preparation; choose
local execution on your analysis machine or configure an existing managed queue.
Opening a report never starts prediction or crosses a manual STOP.

After preparation, or to reuse an existing sealed model checkpoint:

```sh
gene-context project attach --project work/my-study --target family_a \
  --stage stage2 --results /data/prepared-models
gene-context project review --project work/my-study --target family_a
```

The resulting page contains STOP2 and a separate STOP3 synteny view. Once you
have reviewed and exported `STOP4-input-*.json`:

```sh
gene-context project handoff --project work/my-study --target family_a \
  --file /path/to/STOP4-input.json
gene-context project run --project work/my-study --target family_a --stage stage4
gene-context project report --project work/my-study
```

The handoff must match the exact prepared models, sequences, genomic context,
review workspace, alignment and tree. Changed inputs make downstream results
stale; prior attempts and decision exports stay in history.
[Study configuration and commands](docs/PROJECTS.md).

## Reports and optional portable apps

An ordinary report needs no bundled executable. Build an optional local
launcher for your platform from the same source when local interactive
calculation is needed (Go 1.23+):

```sh
python scripts/build_review.py --out /tmp/gcf-review-runtime
export GCF_REVIEW_BIN_DIR=/tmp/gcf-review-runtime
gene-context review-package --results /path/to/model-review --out work/review-app
```

Extract the resulting app and run `Start_review` for your operating system.
Tools are obtained only on an explicit calculation when missing. Saved results
and choices stay in the extracted folder. `--platform all` is an explicit build
option; six platform executables are not shipped in the Python wheel or source tree.

## Examples, compatibility and development

[Minimal synthetic study](examples/minimal/README.md): two unrelated artificial
families and two artificial assemblies, with protein and CDS input examples.
These are software controls, not biological results.

Existing checkpoints stay at their original locations. Projects created before
the sequential workflow retain their legacy route. Completed reports can be
[rebuilt without repeating scientific calculations](docs/COMPLETED_REVIEW_REBUILD.md).
[Earlier CLI reference](docs/CLI_REFERENCE.md).

```sh
python -m unittest discover -s tests
npm ci && npm test
(cd portable-review && go test ./...)
```

See [architecture](docs/ARCHITECTURE.md) for the boundaries between study state,
scientific engines, review pages and optional runtimes.
