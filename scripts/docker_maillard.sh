#!/usr/bin/env bash

set -euo pipefail

CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_validation}"
IMAGE_NAME="${MAILLARD_IMAGE_NAME:-condaforge/miniforge3}"
WORKSPACE_DIR="${MAILLARD_WORKSPACE_DIR:-$PWD}"
WORKSPACE_MOUNT="/workspace"
CONDA_SH="/opt/conda/etc/profile.d/conda.sh"
ENV_NAME="maillard"

usage() {
  cat <<'EOF'
Usage: ./scripts/docker_maillard.sh <command> [args...]

Commands:
  up           Create or start the Docker container.
  bootstrap    Create/update the maillard conda env and apply required patches.
  shell        Open an interactive shell in /workspace with the env activated.
  run CMD...   Run an arbitrary command inside the activated env.
  pytest ...   Run pytest inside the activated env.
  core         Run the core correctness lane.
  scientific   Run the scientific validation lane.
  qm-heavy     Run the QM / external-backend lane.
  hofmann      Generate the Hofmann diagnostic snapshot.
  targets BENCH [TYPE]
               Generate a benchmark target snapshot (default TYPE=desirable; aliases: off_flavour, off-flavour, competing).
  targets-report
               Generate results/validation/benchmark_targets.{md,json}.
  index        Generate results/validation/benchmark_index.{md,json}.
  summary      Generate results/validation/benchmark_summary.{md,json}.
  status       Show container and environment status.
EOF
}

core_lane() {
  run_in_env "python -m pytest tests/unit tests/integration"
}

scientific_lane() {
  run_in_env "python scripts/generate_benchmark_summary.py"
  run_in_env "python scripts/generate_benchmark_index.py"
  run_in_env "python scripts/generate_benchmark_targets.py"
  run_in_env "python -m pytest tests/scientific tests/unit/test_budget_projection.py tests/unit/test_safety_and_flux.py tests/integration/test_recommendation_engine.py"
}

qm_heavy_lane() {
  run_in_env "python -m pytest tests/qm tests/integration/test_cantera_sim.py tests/integration/test_fft_bottleneck.py tests/integration/test_regression.py"
}

hofmann_diagnostic() {
  run_in_env "python scripts/diagnose_benchmark_selectivity.py --lit data/benchmarks/cys_ribose_140C_Hofmann1998.json"
}

targets_snapshot() {
  local benchmark_path="$1"
  local target_type="${2:-desirable}"
  run_in_env "python scripts/diagnose_benchmark_selectivity.py --lit '$benchmark_path' --targets --target-type '$target_type'"
}

container_exists() {
  docker ps -a --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"
}

container_running() {
  docker ps --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"
}

ensure_container() {
  if ! container_exists; then
    docker run -d \
      --platform linux/amd64 \
      --name "$CONTAINER_NAME" \
      -v "$WORKSPACE_DIR:$WORKSPACE_MOUNT" \
      -w "$WORKSPACE_MOUNT" \
      "$IMAGE_NAME" \
      sleep infinity >/dev/null
  elif ! container_running; then
    docker start "$CONTAINER_NAME" >/dev/null
  fi
}

env_exists() {
  docker exec "$CONTAINER_NAME" bash -lc "source '$CONDA_SH' && conda env list | awk '{print \$1}' | grep -qx '$ENV_NAME'"
}

run_in_env() {
  ensure_container
  local command="$1"
  docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; set +u; conda activate '$ENV_NAME'; set -u; cd '$WORKSPACE_MOUNT'; $command"
}

bootstrap_env() {
  ensure_container

  if ! env_exists; then
    docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; conda create -n '$ENV_NAME' python=3.12 -y"
  fi

  run_in_env "conda install -y -c conda-forge jax jaxlib wget xz"
  run_in_env "pip install --index-url https://download.pytorch.org/whl/cpu torch"
  run_in_env "conda env update -n '$ENV_NAME' --file environment.yml"

  run_in_env '
    if ! command -v xtbiff >/dev/null 2>&1; then
      wget -q https://github.com/grimme-lab/xtbiff/releases/download/v1.1/xtbiff.tar.xz
      tar -xf xtbiff.tar.xz
      mv xtbiff "$CONDA_PREFIX/bin/xtbiff"
      chmod +x "$CONDA_PREFIX/bin/xtbiff"
      rm -f xtbiff.tar.xz
    fi
  '

  run_in_env "python - <<'PY'
from pathlib import Path
import os

patches = [
    (
        Path(os.environ['CONDA_PREFIX']) / 'lib/python3.12/site-packages/e3nn/o3/_wigner.py',
        \"torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'))\",
        \"torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'), weights_only=False)\",
    ),
    (
        Path(os.environ['CONDA_PREFIX']) / 'lib/python3.12/site-packages/mace/calculators/mace.py',
        \"torch.load(f=model_path, map_location=device)\",
        \"torch.load(f=model_path, map_location=device, weights_only=False)\",
    ),
]

for path, old, new in patches:
    if not path.exists():
        continue
    text = path.read_text()
    if new in text:
        continue
    if old in text:
        path.write_text(text.replace(old, new))
PY"
}

status() {
  ensure_container
  docker ps -a --filter "name=$CONTAINER_NAME"
  if env_exists; then
    run_in_env "python --version && which crest && which xtbiff"
  else
    echo "Conda environment '$ENV_NAME' has not been created yet. Run: ./scripts/docker_maillard.sh bootstrap"
  fi
}

cmd="${1:-}"
case "$cmd" in
  up)
    ensure_container
    ;;
  bootstrap)
    bootstrap_env
    ;;
  shell)
    ensure_container
    docker exec -it "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; set +u; conda activate '$ENV_NAME'; set -u; cd '$WORKSPACE_MOUNT'; exec bash -i"
    ;;
  run)
    shift
    run_in_env "$*"
    ;;
  pytest)
    shift
    if [ "$#" -eq 0 ]; then
      run_in_env "python -m pytest tests/"
    else
      run_in_env "python -m pytest $*"
    fi
    ;;
  core)
    core_lane
    ;;
  scientific)
    scientific_lane
    ;;
  qm-heavy)
    qm_heavy_lane
    ;;
  hofmann)
    hofmann_diagnostic
    ;;
  targets)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh targets BENCHMARK_JSON [TYPE]" >&2
      exit 1
    fi
    targets_snapshot "$1" "${2:-desirable}"
    ;;
  targets-report)
    run_in_env "python scripts/generate_benchmark_targets.py"
    ;;
  index)
    run_in_env "python scripts/generate_benchmark_index.py"
    ;;
  summary)
    run_in_env "python scripts/generate_benchmark_summary.py"
    ;;
  status)
    status
    ;;
  *)
    usage
    exit 1
    ;;
esac