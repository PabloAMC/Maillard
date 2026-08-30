#!/usr/bin/env bash
# Build the `maillard` conda environment on a Linux CI runner.
#
# This mirrors the sequence in the repo's Dockerfile, which is the only
# environment build in this repo that is known to work on Linux. Two deviations
# from a plain `conda env create -f environment.yml` matter, and both come
# straight from that Dockerfile:
#
#   1. The `pytorch::pytorch=*=*cpu*` line is stripped. Pinning a build string
#      out of the `pytorch` channel is a slow and fragile solve on linux-64; the
#      Dockerfile filters the same line out for the same reason.
#
#   2. CPU torch is installed from PyTorch's CPU wheel index *before* the env
#      file is applied. environment.yml's pip section contains mace-torch, whose
#      dependency resolution otherwise pulls the default CUDA torch off PyPI --
#      several GB of CUDA payload that a CPU runner cannot use and that
#      frequently exhausts the runner's disk.
#
# Usage (inside a conda-activated `maillard` env, e.g. after setup-miniconda):
#     bash scripts/ci/setup_env.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SOURCE_ENV="${ROOT}/environment.yml"
# Written outside the repo tree on purpose: `environment.ci.yml` is not in
# .gitignore, and a stray untracked file in the working tree is a nuisance for
# anyone running this script locally.
FILTERED_ENV="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/environment.ci.yml"
ENV_NAME="${MAILLARD_ENV_NAME:-maillard}"

if [[ ! -f "${SOURCE_ENV}" ]]; then
  echo "setup_env: ${SOURCE_ENV} not found" >&2
  exit 1
fi

echo "==> Filtering environment.yml (dropping pytorch::pytorch build-string pin)"
grep -v 'pytorch::pytorch=' "${SOURCE_ENV}" > "${FILTERED_ENV}"
echo "    removed: $(grep -c 'pytorch::pytorch=' "${SOURCE_ENV}" || true) line(s)"

echo "==> Installing CPU-only torch ahead of the pip section"
python -m pip install --no-cache-dir \
  --index-url https://download.pytorch.org/whl/cpu \
  torch

# No --prune: it can evict the pip-installed CPU torch we just placed, which is
# the whole point of installing it first. The Dockerfile omits it for the same reason.
echo "==> Applying ${FILTERED_ENV} to conda env '${ENV_NAME}'"
conda env update -n "${ENV_NAME}" --file "${FILTERED_ENV}"

# environment.yml pins numpy>=2.0,<2.4 (2026-08-26): NumPy 2.4 changed the
# einsum_path contraction-tuple arity and breaks PySCF 2.11.0, taking down the
# whole QM lane. Assert the pin actually held rather than trusting the solve.
python - <<'PY'
import sys
import numpy

# Parsed by hand so this check does not itself depend on `packaging` being present.
major, minor = (int(part) for part in numpy.__version__.split(".")[:2])
print(f"numpy {numpy.__version__}")
if not ((major, minor) >= (2, 0) and (major, minor) < (2, 4)):
    sys.exit(
        f"FAIL: numpy {numpy.__version__} is outside the environment.yml pin "
        ">=2.0,<2.4. PySCF 2.11.0 breaks on numpy>=2.4."
    )
PY

echo "==> Environment ready"
conda list -n "${ENV_NAME}" | grep -Ei '^(numpy|scipy|cantera|rdkit|pyscf|torch|ase|matplotlib) ' || true
