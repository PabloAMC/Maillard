#!/usr/bin/env bash
# Build the `maillard` conda environment on a Linux CI runner.
#
# 2026-09-01 (cleaning branch, Phase 1a): with the QM/DFT lane gone, environment.yml
# is a plain conda-forge solve. The former torch-first / pytorch-pin-stripping /
# numpy<2.4 dance existed only for mace-torch and PySCF and is no longer needed.
#
# Usage (inside a conda-activated `maillard` env, e.g. after setup-miniconda):
#     bash scripts/ci/setup_env.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SOURCE_ENV="${ROOT}/environment.yml"
ENV_NAME="${MAILLARD_ENV_NAME:-maillard}"

if [[ ! -f "${SOURCE_ENV}" ]]; then
  echo "setup_env: ${SOURCE_ENV} not found" >&2
  exit 1
fi

echo "==> Applying ${SOURCE_ENV} to conda env '${ENV_NAME}'"
conda env update -n "${ENV_NAME}" --file "${SOURCE_ENV}"

echo "==> Environment ready"
conda list -n "${ENV_NAME}" | grep -Ei '^(numpy|scipy|cantera|rdkit|matplotlib|pyyaml) ' || true
