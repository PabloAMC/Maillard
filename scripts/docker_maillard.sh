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
  stability    Run the Tier 0/1 stability gate.
  core         Run the core correctness lane.
  scientific-fast
               Run the fast scientific regression lane selected by pytest markers.
  reporting-fast
               Generate the lightweight matrix reporting bundle for fast product iteration.
  ml-accelerator-quick-tests
               Run the fast no-wet-lab iteration lane for ML accelerator adapters, drift guards, and TS-seed gating.
  kinetics-validation
               Run the slower Cantera/kinetics validation lane selected by pytest markers.
  scientific   Run the scientific validation lane.
  qm-heavy     Run the QM / external-backend lane.
  hofmann      Generate the Hofmann diagnostic snapshot.
  targets BENCH [TYPE]
               Generate a benchmark target snapshot (default TYPE=desirable; aliases: off_flavour, off-flavour, competing).
  targets-report
               Generate results/validation/benchmark_targets.{md,json}.
  matrix-deltas
               Generate results/validation/matrix_benchmark_deltas.{md,json}.
  matrix-assertions
               Generate results/validation/matrix_benchmark_assertions.{md,json}.
  matrix-evidence
               Generate results/validation/matrix_benchmark_evidence.{md,json}.
  matrix-readiness
               Generate results/validation/matrix_promotion_readiness.{md,json}.
  matrix-promotion-contract
               Generate results/validation/matrix_promotion_contract.{md,json}.
  matrix-closure-audit
               Generate results/validation/matrix_observable_closure_audit.{md,json}.
  matrix-primary-benchmark-campaign
               Generate results/validation/matrix_primary_benchmark_campaign.{md,json}.
  extrusion-external-closure
               Generate results/validation/extrusion_external_closure.{md,json}.
  dha-lysinoalanine-external-package
               Generate results/validation/dha_lysinoalanine_external_package.{md,json}.
  hexanal-nonanal-calibration
               Generate results/validation/hexanal_nonanal_calibration_closure.{md,json}.
  experiment-intake-schema
               Generate results/validation/matrix_experiment_intake_schema.{md,json}.
  compare-experiment INTAKE
               Generate support-delta artifacts for a matrix experiment intake payload.
  materialize-experiment-benchmark INTAKE OUTPUT
               Convert a matrix experiment intake payload into a benchmark JSON artifact.
  literature-learning-loop
               Generate results/validation/literature_learning_loop.{md,json}.
  family-promotion-state
               Generate results/validation/family_promotion_state.{md,json}.
  matrix-family-coverage
               Generate results/validation/matrix_family_coverage.{md,json}.
  matrix-family-ranking
               Generate results/validation/matrix_family_priority_ranking.{md,json}.
  matrix-family-next-action
               Generate results/validation/matrix_family_next_action.{md,json}.
  family-sensitivity [TARGET_TAG] [DELTA_KCAL]
               Generate results/validation/family_sensitivity.{md,json}.
  mycoprotein-reference
               Generate results/validation/mycoprotein_reference.{md,json}.
  pea-soy-external-evidence
               Generate results/validation/pea_soy_external_evidence.{md,json}.
  pea-soy-mixed-external-package
               Generate results/validation/pea_soy_mixed_external_package.{md,json}.
  objective-progress
               Generate results/validation/objective_progress.{md,json} and docs/assets/objective_progress.png.
  primary-matrix-external-package
               Generate results/validation/primary_matrix_external_package.{md,json}.
  hexanal-nonanal-resolution
               Generate results/validation/hexanal_nonanal_resolution.{md,json}.
  scope-gap-guard
               Generate results/validation/scope_gap_guard.{md,json}.
  family-barrier-progress
               Generate results/validation/family_barrier_progress.{md,json} and docs/assets/family_barrier_progress.png.
  dft-coverage-map
               Generate data/lit/dft_coverage_map.json and results/validation/dft_coverage_map.{md,json}.
  mechanistic-refinement
               Generate the selective mechanistic refinement artifact bundle including governance.
  geometry-preopt-assessment
               Generate the geometry preoptimization benchmark and assessment artifacts.
  ts-seed-recovery-assessment
               Generate the TS-seed recovery benchmark and assessment artifacts.
  ml-accelerator-assessment|mlp-assessment
               Generate the full offline ML accelerator artifact bundle: reaction benchmark, geometry preoptimization benchmark, TS-seed recovery benchmark, geometry and TS-seed assessments, ML accelerator assessment, external landscape, and adoption notes.
  Legacy aliases: p3-refinement, p4-geometry-assessment, p4-ts-seed-assessment, p4-mlp-assessment, p4-quick-tests
  skip-registry
               Generate results/validation/skip_registry.{md,json}.
  matrix-branch-deltas [BASE_REF]
               Generate results/validation/matrix_branch_delta_report.{md,json} against BASE_REF (default: main).
  kinetics-timeseries ARGS...
               Run scripts/run_cantera_kinetics.py inside Docker to generate CSV/plot time-series outputs.
  coverage-gaps
               Generate results/validation/benchmark_coverage_gaps.{md,json}.
  validation-figures
               Generate validation_overview, family_validation_overview, and validated_envelope artifacts.
  thermo-gating
               Generate results/validation/thermodynamic_gating_audit.{md,json}.
  validated-envelope
               Generate results/validation/validated_envelope.{md,json}.
  explain-formulation NAME [TARGET_TAG] [MINIMIZE_TAG]
               Generate a formulation explainability artifact in results/validation.
  campaign SPEC [OUTPUT_DIR]
               Run a shareable campaign spec and generate campaign artifacts.
  index        Generate results/validation/benchmark_index.{md,json}.
  summary      Generate results/validation/benchmark_summary.{md,json}.
  notebook     Launch a Jupyter notebook server on port 8888.
  status       Show container and environment status.
EOF
}

core_lane() {
  run_in_env "python -m pytest tests/unit tests/integration"
}

stability_lane() {
  run_in_env "python -m pytest \
    tests/unit/test_arrhenius_params.py \
    tests/unit/test_headspace.py \
    tests/unit/test_safety_and_flux.py \
    tests/integration/test_recommendation_engine.py \
    tests/integration/test_cantera_sim.py \
    tests/integration/test_fft_bottleneck.py \
    tests/integration/test_regression.py \
    tests/integration/test_barrier_calibration.py"
}

scientific_lane() {
  run_in_env "python scripts/generators/generate_benchmark_summary.py"
  run_in_env "python scripts/generators/generate_benchmark_index.py"
  run_in_env "python scripts/generators/generate_benchmark_targets.py"
  run_in_env "python scripts/generators/generate_matrix_benchmark_deltas.py"
  run_in_env "python scripts/generators/generate_matrix_benchmark_assertions.py"
  run_in_env "python scripts/generators/generate_matrix_benchmark_evidence.py"
  run_in_env "python scripts/generators/generate_matrix_promotion_readiness.py"
  run_in_env "python scripts/generators/generate_matrix_promotion_contract.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_observable_closure_audit.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_primary_benchmark_campaign.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_extrusion_external_closure.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_hexanal_nonanal_calibration.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_experiment_intake_schema.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_family_coverage.py"
  run_in_env "python scripts/generators/generate_matrix_family_priority_ranking.py"
  run_in_env "python scripts/generators/generate_matrix_family_next_action.py"
  run_in_env "python scripts/generators/generate_mycoprotein_reference.py"
  run_in_env "python scripts/generators/generate_pea_soy_external_evidence.py"
  run_in_env "python scripts/generators/generate_pea_soy_mixed_external_package.py"
  run_in_env "python scripts/generators/generate_extrusion_external_closure.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_dha_lysinoalanine_external_package.py"
  run_in_env "python scripts/generators/generate_objective_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_primary_matrix_external_package.py"
  run_in_env "python scripts/generators/generate_hexanal_nonanal_resolution.py"
  run_in_env "python scripts/generators/generate_scope_gap_guard.py"
  run_in_env "python scripts/generators/generate_family_barrier_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_refinement_campaign.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_mlp_geometry_assessment.py"
  run_in_env "python scripts/generators/generate_mlp_ts_seed_assessment.py"
  run_in_env "python scripts/generators/generate_mlp_assessment.py"
  run_in_env "python scripts/generators/generate_literature_learning_loop.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_promotion_state.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_lane_validation.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_deviation_audit.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_validation_figures.py"
  run_in_env "python scripts/generators/generate_family_validation_figures.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_validated_envelope_report.py"
  run_in_env "python scripts/generators/generate_thermodynamic_gating_audit.py"
  s9_authority_lane
  scientific_fast_lane
}

scientific_fast_lane() {
  run_in_env "python -m pytest -m 'scientific_regression and not slow' tests/scientific tests/unit/test_budget_projection.py tests/unit/test_safety_and_flux.py tests/integration/test_recommendation_engine.py"
}

reporting_fast_lane() {
  run_in_env "python scripts/generators/generate_matrix_target_status.py"
  run_in_env "python scripts/generators/generate_hexanal_nonanal_calibration.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_family_coverage.py"
  run_in_env "python scripts/generators/generate_matrix_family_priority_ranking.py"
  run_in_env "python scripts/generators/generate_matrix_family_next_action.py"
  run_in_env "python scripts/generators/generate_mycoprotein_reference.py"
  run_in_env "python scripts/generators/generate_pea_soy_external_evidence.py"
  run_in_env "python scripts/generators/generate_pea_soy_mixed_external_package.py"
  run_in_env "python scripts/generators/generate_extrusion_external_closure.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_dha_lysinoalanine_external_package.py"
  run_in_env "python scripts/generators/generate_objective_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_primary_matrix_external_package.py"
  run_in_env "python scripts/generators/generate_hexanal_nonanal_resolution.py"
  run_in_env "python scripts/generators/generate_scope_gap_guard.py"
  run_in_env "python scripts/generators/generate_family_barrier_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
}

s9_authority_lane() {
  run_in_env "python -m pytest tests/unit/test_authority_benchmark_data.py tests/benchmarks/test_quasi_harmonic_correction.py tests/benchmarks/test_barrier_benchmarks.py tests/benchmarks/test_irc_validation.py tests/unit/test_dft_refiner_irc_api.py tests/unit/test_skip_policy_registry.py"
  run_in_env "python scripts/generators/generate_skip_registry.py --output-dir results/validation"
}

p4_quick_tests_lane() {
  run_in_env "python -m pytest -m 'not slow' tests/unit/test_mlp_optimizer_unit.py tests/unit/test_mlp_backend_adapters.py tests/unit/test_chemistry_benchmark_validator.py tests/unit/test_ts_seed_benchmark.py tests/unit/test_ts_seed_benchmark_validator.py"
}

kinetics_validation_lane() {
  run_in_env "python -m pytest -m kinetics_validation tests/integration/test_cantera_sim.py tests/integration/test_temp_profiles.py"
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
      -p 8888:8888 \
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
  # We use 'export MKL_INTERFACE_LAYER=LP64' and 'set +u' to bypass common Conda activation script bugs
  docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; export MKL_INTERFACE_LAYER=LP64; set +u; conda activate '$ENV_NAME'; set -u; cd '$WORKSPACE_MOUNT'; $command"
}

shell_join() {
  local quoted=()
  local arg
  for arg in "$@"; do
    quoted+=("$(printf '%q' "$arg")")
  done
  printf '%s' "${quoted[*]}"
}

bootstrap_env() {
  ensure_container

  # Install system LaTeX packages needed for scienceplots LaTeX rendering mode
  docker exec "$CONTAINER_NAME" bash -lc "apt-get update -qq && apt-get install -y texlive texlive-latex-extra texlive-fonts-recommended dvipng cm-super 2>&1 | tail -3" || true

  if ! env_exists; then
    docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; conda create -n '$ENV_NAME' python=3.12 -y"
  fi

  run_in_env "conda install -y -c conda-forge jax jaxlib wget xz jupyter ipywidgets"
  run_in_env "pip install --index-url https://download.pytorch.org/whl/cpu torch==2.5.1"
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
  stability)
    stability_lane
    ;;
  core)
    core_lane
    ;;
  scientific)
    scientific_lane
    ;;
  scientific-fast)
    scientific_fast_lane
    ;;
  reporting-fast)
    reporting_fast_lane
    ;;
  family-sensitivity)
    shift
    run_in_env "python scripts/generators/generate_family_sensitivity.py --output-dir results/validation --target-tag '${1:-meaty}' --delta-kcal '${2:-3.0}'"
    ;;
  ml-accelerator-quick-tests|p4-quick-tests)
    p4_quick_tests_lane
    ;;
  kinetics-validation)
    kinetics_validation_lane
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
    run_in_env "python scripts/generators/generate_benchmark_targets.py"
    ;;
  matrix-deltas)
    run_in_env "python scripts/generators/generate_matrix_benchmark_deltas.py"
    ;;
  matrix-assertions)
    run_in_env "python scripts/generators/generate_matrix_benchmark_assertions.py"
    ;;
  matrix-evidence)
    run_in_env "python scripts/generators/generate_matrix_benchmark_evidence.py"
    ;;
  matrix-readiness)
    run_in_env "python scripts/generators/generate_matrix_promotion_readiness.py"
    ;;
  matrix-promotion-contract)
    run_in_env "python scripts/generators/generate_matrix_promotion_contract.py --output-dir results/validation"
    ;;
  matrix-closure-audit)
    run_in_env "python scripts/generators/generate_matrix_observable_closure_audit.py --output-dir results/validation"
    ;;
  matrix-primary-benchmark-campaign)
    run_in_env "python scripts/generators/generate_matrix_primary_benchmark_campaign.py --output-dir results/validation"
    ;;
  extrusion-external-closure)
    run_in_env "python scripts/generators/generate_extrusion_external_closure.py --output-dir results/validation"
    ;;
  dha-lysinoalanine-external-package)
    run_in_env "python scripts/generators/generate_dha_lysinoalanine_external_package.py"
    ;;
  hexanal-nonanal-calibration)
    run_in_env "python scripts/generators/generate_hexanal_nonanal_calibration.py --output-dir results/validation"
    ;;
  experiment-intake-schema)
    run_in_env "python scripts/generators/generate_matrix_experiment_intake_schema.py --output-dir results/validation"
    ;;
  literature-learning-loop)
    run_in_env "python scripts/generators/generate_literature_learning_loop.py --output-dir results/validation"
    ;;
  family-promotion-state)
    run_in_env "python scripts/generators/generate_family_promotion_state.py --output-dir results/validation"
    ;;
  matrix-family-coverage)
    run_in_env "python scripts/generators/generate_matrix_family_coverage.py"
    ;;
  matrix-family-ranking)
    run_in_env "python scripts/generators/generate_matrix_family_priority_ranking.py"
    ;;
  matrix-family-next-action)
    run_in_env "python scripts/generators/generate_matrix_family_next_action.py"
    ;;
  mycoprotein-reference)
    run_in_env "python scripts/generators/generate_mycoprotein_reference.py"
    ;;
  pea-soy-external-evidence)
    run_in_env "python scripts/generators/generate_pea_soy_external_evidence.py"
    ;;
  pea-soy-mixed-external-package)
    run_in_env "python scripts/generators/generate_pea_soy_mixed_external_package.py"
    ;;
  objective-progress)
    run_in_env "python scripts/generators/generate_objective_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
    ;;
  primary-matrix-external-package)
    run_in_env "python scripts/generators/generate_primary_matrix_external_package.py"
    ;;
  hexanal-nonanal-resolution)
    run_in_env "python scripts/generators/generate_hexanal_nonanal_resolution.py"
    ;;
  scope-gap-guard)
    run_in_env "python scripts/generators/generate_scope_gap_guard.py"
    ;;
  family-barrier-progress)
    run_in_env "python scripts/generators/generate_family_barrier_progress.py --output-dir results/validation --docs-asset-dir docs/assets"
    ;;
  dft-coverage-map)
    run_in_env "python scripts/generators/generate_dft_coverage_map.py"
    ;;
  mechanistic-refinement|p3-refinement)
    run_in_env "python scripts/generators/generate_refinement_campaign.py --output-dir results/validation"
    ;;
  geometry-preopt-assessment|p4-geometry-assessment)
    run_in_env "python scripts/generators/generate_mlp_geometry_assessment.py"
    ;;
  ts-seed-recovery-assessment|p4-ts-seed-assessment)
    run_in_env "python scripts/generators/generate_mlp_ts_seed_assessment.py"
    ;;
  ml-accelerator-assessment|mlp-assessment|p4-mlp-assessment)
    run_in_env "python scripts/generators/generate_mlp_assessment.py"
    ;;
  skip-registry|s9-skip-registry)
    run_in_env "python scripts/generators/generate_skip_registry.py --output-dir results/validation"
    ;;
  compare-experiment)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh compare-experiment INTAKE_FILE" >&2
      exit 1
    fi
    run_in_env "python scripts/generators/compare_matrix_experiment_intake.py --experiment '$1' --output-dir results/validation"
    ;;
  materialize-experiment-benchmark)
    shift
    if [ "$#" -lt 2 ]; then
      echo "Usage: ./scripts/docker_maillard.sh materialize-experiment-benchmark INTAKE_FILE OUTPUT_BENCHMARK" >&2
      exit 1
    fi
    run_in_env "python scripts/generators/materialize_matrix_experiment_benchmark.py $(shell_join "$1" "$2")"
    ;;
  matrix-branch-deltas)
    shift
    run_in_env "python scripts/compare_matrix_benchmark_branches.py --base-ref '${1:-main}'"
    ;;
  kinetics-timeseries)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh kinetics-timeseries --precursors ribose:0.1,cysteine:0.1 [run_cantera_kinetics args...]" >&2
      exit 1
    fi
    run_in_env "python scripts/run_cantera_kinetics.py $(shell_join "$@")"
    ;;
  coverage-gaps)
    run_in_env "python scripts/generators/generate_benchmark_coverage_gaps.py"
    ;;
  validation-figures)
    run_in_env "python scripts/generators/generate_family_lane_validation.py --output-dir results/validation"
    run_in_env "python scripts/generators/generate_family_deviation_audit.py --output-dir results/validation"
    run_in_env "python scripts/generators/generate_validation_figures.py"
    run_in_env "python scripts/generators/generate_family_validation_figures.py --output-dir results/validation --docs-asset-dir docs/assets"
    run_in_env "python scripts/generators/generate_validated_envelope_report.py"
    ;;
  thermo-gating)
    run_in_env "python scripts/generators/generate_thermodynamic_gating_audit.py"
    ;;
  validated-envelope)
    run_in_env "python scripts/generators/generate_validated_envelope_report.py"
    ;;
  explain-formulation)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh explain-formulation NAME [TARGET_TAG] [MINIMIZE_TAG]" >&2
      exit 1
    fi
    run_in_env "python scripts/explain_formulation.py --name '$1' --target-tag '${2:-meaty}' --minimize-tag '${3:-beany}'"
    ;;
  campaign)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh campaign SPEC [OUTPUT_DIR]" >&2
      exit 1
    fi
    if [ "$#" -ge 2 ]; then
      run_in_env "python scripts/run_campaign.py --spec '$1' --output-dir '$2'"
    else
      run_in_env "python scripts/run_campaign.py --spec '$1'"
    fi
    ;;
  index)
    run_in_env "python scripts/generators/generate_benchmark_index.py"
    ;;
  summary)
    run_in_env "python scripts/generators/generate_benchmark_summary.py"
    ;;
  notebook)
    run_in_env "jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root"
    ;;
  status)
    status
    ;;
  *)
    usage
    exit 1
    ;;
esac