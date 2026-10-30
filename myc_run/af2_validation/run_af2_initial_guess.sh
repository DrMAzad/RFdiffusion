#!/usr/bin/env bash
# Route B (gold standard): AF2 initial-guess via Baker-lab dl_binder_design.
# Seeds AF2 with the RFdiffusion backbone -> far more accurate for binders than
# vanilla AF2. Needs a GPU, AF2 params, and the dl_binder_design repo.
#
# One-time setup:
#   git clone https://github.com/nrbennet/dl_binder_design
#   cd dl_binder_design && (create its conda env per its README)
#   # download AF2 model params into af2_initial_guess/model_weights/params
#
# Inputs: the RFdiffusion design PDBs (outputs/scan/*.pdb) with ProteinMPNN sequences
# threaded in. dl_binder_design does the threading itself from your MPNN run.
set -euo pipefail
DLBD=${DLBD:-/path/to/dl_binder_design}          # <-- set this
DESIGNS=${1:-/home/user/drmazad/rfdiffusion/outputs/scan}
WORK=${2:-$PWD/af2ig_work}
mkdir -p "$WORK"

# 1. pack RFdiffusion PDBs into a silent file
"$DLBD/include/silent_tools/silentfrompdbs" "$DESIGNS"/*.pdb > "$WORK/designs.silent"

# 2. ProteinMPNN sequence design on the backbones (binder chain only)
python "$DLBD/mpnn_fr/dl_interface_design.py" \
  -silent "$WORK/designs.silent" -outsilent "$WORK/mpnn.silent"

# 3. AF2 initial-guess -> scorefile with pae_interaction
python "$DLBD/af2_initial_guess/predict.py" \
  -silent "$WORK/mpnn.silent" -outsilent "$WORK/af2.silent" -scorefilename "$WORK/af2_scores.sc"

echo "Done. Filter $WORK/af2_scores.sc on pae_interaction < 10 (and plddt_binder > 80)."
echo "e.g.:  awk 'NR==1||\$<pae_col> < 10' $WORK/af2_scores.sc"
