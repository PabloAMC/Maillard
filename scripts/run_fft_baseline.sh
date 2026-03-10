#!/bin/bash
# scripts/run_fft_baseline.sh
# Baseline simulation for Ribose/Cys/Leu FFT bottleneck investigation.

mkdir -p results/fft_baseline

export PYTHONPATH=$PYTHONPATH:$(pwd)
.venv/bin/python scripts/run_cantera_kinetics.py \
  --precursors "ribose:0.1,cysteine:0.05,leucine:0.05" \
  --temp 150 \
  --time 3600 \
  --from-smirks \
  --input results/maillard_results.db \
  --track "H2S,H2,furfural,2-furfurylthiol" \
  --verbose-reactions \
  --output results/fft_baseline/ribose_cys_leu \
  --no-gating
