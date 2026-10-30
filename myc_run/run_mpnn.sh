#!/usr/bin/env bash
# ProteinMPNN sequence design on RFdiffusion binder backbones.
# Fixes the Myc target chain (A) and designs the binder chain (B).
set -euo pipefail
source /home/user/se3cpu/bin/activate
MPNN=/home/user/ProteinMPNN
PDBDIR=${1:-/home/user/drmazad/rfdiffusion/outputs/scan}   # dir of RFdiffusion .pdb designs
OUT=${2:-/home/user/myc_run/mpnn}                           # output dir
NSEQ=${NSEQ:-8}

mkdir -p "$OUT"
python "$MPNN/helper_scripts/parse_multiple_chains.py" \
  --input_path="$PDBDIR" --output_path="$OUT/parsed.jsonl"
python "$MPNN/helper_scripts/assign_fixed_chains.py" \
  --input_path="$OUT/parsed.jsonl" --output_path="$OUT/assigned.jsonl" --chain_list "B"
python "$MPNN/protein_mpnn_run.py" \
  --jsonl_path "$OUT/parsed.jsonl" \
  --chain_id_jsonl "$OUT/assigned.jsonl" \
  --out_folder "$OUT" \
  --num_seq_per_target "$NSEQ" --sampling_temp "0.1" --seed 37 --batch_size 1

echo "MPNN_DONE  -> sequences in $OUT/seqs/*.fa"
