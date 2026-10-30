#!/usr/bin/env bash
set -euo pipefail
cd /home/user/drmazad/rfdiffusion
source /home/user/se3cpu/bin/activate
export DGLBACKEND=pytorch OMP_NUM_THREADS=4
for L in 6 8 10 12 14; do
  echo "===== LENGTH $L ====="
  ./scripts/run_inference.py \
    inference.input_pdb=/home/user/myc_run/myc_target_zipperC.pdb \
    inference.output_prefix=outputs/scan/len${L} \
    "contigmap.contigs=[A955-984/0 ${L}-${L}]" \
    'ppi.hotspot_res=[A984,A981,A978]' \
    inference.num_designs=8 \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
done
echo "ALL_SCAN_DONE"
