#!/usr/bin/env bash
# RFdiffusion: shortest-possible linear (non-cyclic) peptide binder to c-Myc,
# anchored at the reactive leucine-zipper cysteine (Cys984 in 1NKP numbering).
#
# Requires: a CUDA GPU, the RFdiffusion SE3nv conda env, and the model weights
# in ./models (download via scripts/download_models.sh in the repo).
# Run from the root of the RFdiffusion checkout.
set -euo pipefail

# --- target ---------------------------------------------------------------
# Use the trimmed C-terminal zipper for speed/focus, or the full bHLHLZ helix
# for more context. Hotspots are the Max-interface knobs around the cysteine.
TARGET=${TARGET:-myc_run/myc_target_zipperC.pdb}   # or myc_target_full.pdb
HOTSPOTS='[A984,A981,A978]'                         # Cys984 = covalent anchor

# contig: keep the whole target chain, chain-break (/0), then an N-residue binder
# NOTE: if you use myc_target_full.pdb, change A955-984 -> A897-984 below.
TARGET_CONTIG='A955-984'

NUM=${NUM:-50}          # designs per length
OUTDIR=${OUTDIR:-outputs/myc_cys984}
mkdir -p "$OUTDIR"

# --- length scan: "as short as possible" ----------------------------------
# RFdiffusion rarely makes a stable interface below ~8-10 res; scan upward and
# keep the shortest length that still passes downstream filters (AF2 pAE/pLDDT).
for L in 6 8 10 12 14; do
  echo ">>> designing ${L}-residue binders to Myc Cys984"
  ./scripts/run_inference.py \
    inference.input_pdb="${TARGET}" \
    inference.output_prefix="${OUTDIR}/len${L}" \
    "contigmap.contigs=[${TARGET_CONTIG}/0 ${L}-${L}]" \
    "ppi.hotspot_res=${HOTSPOTS}" \
    inference.num_designs="${NUM}" \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
done

echo
echo "Done. Backbones + trajectories are in ${OUTDIR}/."
echo "Trajectories (${OUTDIR}/traj/*_pX0_traj.pdb) load straight into the viewer."
