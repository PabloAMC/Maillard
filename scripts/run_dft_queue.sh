#!/usr/bin/env bash
# run_dft_queue.sh — Sequential DFT runner for all computational-gap targets.
# Runs each target one at a time to avoid memory exhaustion.
# Each target gets its own output JSON for fine-grained status tracking.

set -euo pipefail

LOCKFILE="/tmp/maillard_dft_queue.lock"
if ! ( set -o noclobber; echo "$$" > "$LOCKFILE" ) 2> /dev/null; then
    echo "[DFT QUEUE] Error: Another instance of this queue is already running (pid $(cat "$LOCKFILE"))."
    echo "            If this is a mistake, run: rm $LOCKFILE"
    exit 1
fi
trap 'rm -f "$LOCKFILE"; exit $?' INT TERM EXIT

TARGETS=(
  "lysinoalanine_crosslink"
  "aa_ring_open_dicarbonyl"
  "quinone_cys_michael"
  "pe_schiff_base"
  "pe_amadori"
  "asparagine_sugar_explicit_water_cluster"
  "hexanal_radical_quench"
)

# Allow overriding start index: e.g. START_FROM=3 ./scripts/run_dft_queue.sh
START_FROM="${START_FROM:-0}"
SCRIPT="scripts/run_computational_gap_dft.py"
RESULTS_DIR="results/computational_gap_refinement"
CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_validation}"
CONDA_SH="/opt/conda/etc/profile.d/conda.sh"
ENV_NAME="maillard"

run_target() {
  local target="$1"
  local output="${RESULTS_DIR}/${target}_dft_execution.json"
  echo "=================================================="
  echo "[DFT QUEUE] Starting: $target  ($(date -u '+%Y-%m-%dT%H:%M:%SZ'))"
  echo "=================================================="
  python3 -c "
import subprocess, sys
try:
    subprocess.run(sys.argv[1:], timeout=86400, check=True)
except subprocess.TimeoutExpired:
    print('Command timed out after 24h')
    sys.exit(124)
except subprocess.CalledProcessError as e:
    sys.exit(e.returncode)
" docker exec "$CONTAINER_NAME" bash -lc \
    "set -eo pipefail; \
     source '$CONDA_SH'; \
     export MKL_INTERFACE_LAYER=LP64; \
     set +u; conda activate '$ENV_NAME'; set -u; \
     cd /workspace; \
     exec python $SCRIPT --execute --target '$target' --output '$output'"
  echo "[DFT QUEUE] Finished: $target  ($(date -u '+%Y-%m-%dT%H:%M:%SZ'))"
}

idx=0
for target in "${TARGETS[@]}"; do
  if [ "$idx" -lt "$START_FROM" ]; then
    echo "[DFT QUEUE] Skipping $target (START_FROM=$START_FROM)"
    idx=$((idx + 1))
    continue
  fi
  run_target "$target" || {
    echo "[DFT QUEUE] WARNING: $target exited with error. Continuing to next target."
  }
  idx=$((idx + 1))
done

echo "=================================================="
echo "[DFT QUEUE] All targets processed. $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "=================================================="
