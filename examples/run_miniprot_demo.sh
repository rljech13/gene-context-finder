#!/usr/bin/env bash
# Run from the repository root, after installing the environment.
set -euo pipefail
DEMO_ROOT="${1:?Usage: bash examples/run_miniprot_demo.sh /path/to/new/demo}"
SOURCE="https://raw.githubusercontent.com/lh3/miniprot/81f9b93481cec16622eff339fc998e4f38da344e/test"
mkdir -p "$DEMO_ROOT/inputs/genomes" "$DEMO_ROOT/inputs/queries"
curl -fL "$SOURCE/DPP3-hs.gen.fa.gz" -o "$DEMO_ROOT/inputs/genomes/DPP3-hs.gen.fa.gz"
curl -fL "$SOURCE/DPP3-mm.pep.fa.gz" -o "$DEMO_ROOT/inputs/queries/DPP3-mm.pep.fa.gz"
python -m genecontext run \
  --query "$DEMO_ROOT/inputs/queries/DPP3-mm.pep.fa.gz" \
  --query-type protein --genome-dir "$DEMO_ROOT/inputs/genomes" \
  --out "$DEMO_ROOT/results" --threads 2 --flank 1000 --context 10000
