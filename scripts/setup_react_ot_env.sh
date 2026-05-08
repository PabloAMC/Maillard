#!/usr/bin/env bash
# scripts/setup_react_ot_env.sh
#
# Provisions a *secondary* conda env named `react_ot` inside the
# `maillard_validation` Docker container, isolated from the main `maillard`
# env. React-OT (deepprinciple/react-ot) is installed there with a CPU
# torch / torch_geometric stack so the main env stays unaffected.
#
# Idempotent: re-running only adds what's missing.
#
# Override the defaults with environment variables:
#   MAILLARD_CONTAINER_NAME       (default: maillard_validation)
#   REACT_OT_REPO_URL             (default: https://github.com/deepprinciple/react-ot.git)
#   REACT_OT_REF                  (default: main)
#   REACT_OT_REPO_DIR             (default: /opt/react-ot inside the container)
#   REACT_OT_TORCH_VERSION        (default: 2.2.1)
#   REACT_OT_TORCH_INDEX          (default: https://download.pytorch.org/whl/cpu)
#   REACT_OT_PYG_WHEEL_INDEX      (default: https://data.pyg.org/whl/torch-${REACT_OT_TORCH_VERSION}+cpu.html)

set -euo pipefail

CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_validation}"
CONDA_SH="/opt/conda/etc/profile.d/conda.sh"
SECONDARY_ENV="${MAILLARD_REACT_OT_ENV:-react_ot}"
REACT_OT_REPO_URL="${REACT_OT_REPO_URL:-https://github.com/deepprinciple/react-ot.git}"
REACT_OT_REF="${REACT_OT_REF:-main}"
REACT_OT_REPO_DIR="${REACT_OT_REPO_DIR:-/opt/react-ot}"
TORCH_VERSION="${REACT_OT_TORCH_VERSION:-2.2.1}"
TORCH_INDEX="${REACT_OT_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"
PYG_WHEEL_INDEX="${REACT_OT_PYG_WHEEL_INDEX:-https://data.pyg.org/whl/torch-${TORCH_VERSION}+cpu.html}"

ensure_secondary_env() {
  docker exec "$CONTAINER_NAME" bash -lc "
    set -eo pipefail
    source '$CONDA_SH'
    if ! conda env list | awk '{print \$1}' | grep -qx '$SECONDARY_ENV'; then
      conda create -n '$SECONDARY_ENV' python=3.10 -y
    fi
    conda activate '$SECONDARY_ENV'
    conda install -y -c conda-forge ase=3.23.0 rdkit git wget
    pip install --index-url '$TORCH_INDEX' 'torch==${TORCH_VERSION}'
    pip install \
      'numpy==1.26.4' \
      'pandas==2.2.3' \
      'pymatgen==2024.11.13' \
      'pytorch-lightning==2.4.0' \
      torchdiffeq \
      networkx \
      timm \
      lmdb \
      rich
    pip install torch_geometric
    if ! pip install torch_scatter torch_sparse torch_cluster -f '$PYG_WHEEL_INDEX'; then
      echo 'WARNING: CPU wheels for torch_scatter / torch_sparse / torch_cluster are not'
      echo '         available for this Torch version. React-OT inference will likely'
      echo '         fail; fall back to notebooks/react_ot_colab_gpu.ipynb.'
    fi
  "
}

clone_react_ot_repo() {
  docker exec "$CONTAINER_NAME" bash -lc "
    set -eo pipefail
    if [ ! -d '$REACT_OT_REPO_DIR/.git' ]; then
      git clone '$REACT_OT_REPO_URL' '$REACT_OT_REPO_DIR'
    fi
    cd '$REACT_OT_REPO_DIR'
    git fetch --tags --quiet
    git checkout '$REACT_OT_REF'
  "
}

install_react_ot_package() {
  docker exec "$CONTAINER_NAME" bash -lc "
    set -eo pipefail
    source '$CONDA_SH'
    conda activate '$SECONDARY_ENV'
    cd '$REACT_OT_REPO_DIR'
    pip install --no-deps -e .
  "
}

ensure_secondary_env
clone_react_ot_repo
install_react_ot_package

cat <<EOF
React-OT secondary env is ready.

Verify import:
  docker exec $CONTAINER_NAME bash -lc "source $CONDA_SH; conda activate $SECONDARY_ENV; python -c 'import reactot; print(reactot.__file__)'"

Download the pretrained checkpoint into the repo (CPU smoke needs it):
  docker exec $CONTAINER_NAME bash -lc "wget -q https://zenodo.org/records/13131875/files/sb-pretrained.ckpt -O /workspace/models/external/react_ot/sb-pretrained.ckpt"

Then run the smoke gate:
  ./scripts/docker_maillard.sh react-ot-smoke
EOF
