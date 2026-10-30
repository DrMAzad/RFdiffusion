#!/usr/bin/env bash
# Route A: validate Myc binder designs with ColabFold (needs a GPU).
# Install localcolabfold first: https://github.com/YoshitakaMo/localcolabfold
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SEQS=${1:-"$HERE/../mpnn/seqs"}        # ProteinMPNN output dir
OUT=${2:-"$HERE/af2_out"}

python "$HERE/make_colabfold_fasta.py" "$SEQS" "$HERE/complexes.fasta"

# AlphaFold2-multimer for the binder:target complex. ~seconds-minutes/design on GPU.
colabfold_batch \
  --num-recycle 3 \
  --model-type alphafold2_multimer_v3 \
  --num-models 5 \
  "$HERE/complexes.fasta" "$OUT"

python "$HERE/parse_af2.py" "$OUT" "$HERE/af2_metrics.csv"
echo "Done -> $HERE/af2_metrics.csv  (filter: pae_interaction<10 & binder pLDDT>80)"
