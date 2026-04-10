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
    scientific_lane
    ;;
  scientific-fast)
    scientific_fast_lane
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
    run_generator_script generate_benchmark_targets
  matrix-promotion-contract
               Generate results/validation/matrix_promotion_contract.{md,json}.
    run_generator_script generate_matrix_benchmark_deltas
               Generate results/validation/matrix_observable_closure_audit.{md,json}.
  experiment-intake-schema
    run_generator_script generate_matrix_benchmark_assertions
  compare-experiment INTAKE
               Generate support-delta artifacts for a matrix experiment intake payload.
    run_generator_script generate_matrix_benchmark_evidence
               Convert an extrusion external-closure workbook into intake payloads and support-delta artifacts.
  extrusion-follow-on-workbook WORKBOOK
    run_generator_script generate_matrix_promotion_readiness
  extrusion-diagnostic-examples
               Generate and execute synthetic diagnostic examples for both extrusion workbook flows.
    run_generator_script generate_matrix_promotion_contract --output-dir results/validation
               Generate the no-wet-lab computational-gap refinement plan and the xTB/DFT job manifests.
  computational-gap-xtb [TARGET]
    run_generator_script generate_matrix_observable_closure_audit --output-dir results/validation
  computational-gap-dft [TARGET]
               Execute the computational-gap DFT refinement job set (or a single TARGET id).
    run_generator_script generate_matrix_experiment_intake_schema --output-dir results/validation
               Generate the computational-gap DFT ingestion report from the current execution artifacts.
  computational-gap-dft-promote
    run_generator_script generate_literature_learning_loop --output-dir results/validation
  literature-learning-loop
               Generate results/validation/literature_learning_loop.{md,json}.
    run_generator_script generate_family_ingestion_plan --output-dir results/validation
               Generate results/validation/family_ingestion_plan.{md,json}.
  family-promotion-state
    run_generator_script generate_family_promotion_state --output-dir results/validation
  matrix-family-coverage
               Generate results/validation/matrix_family_coverage.{md,json}.
    run_generator_script generate_matrix_family_coverage
               Generate the selective mechanistic refinement governance artifact bundle.
  mlp-assessment
    run_generator_script generate_refinement_governance --output-dir results/validation
  matrix-branch-deltas [BASE_REF]
               Generate results/validation/matrix_branch_delta_report.{md,json} against BASE_REF (default: main).
    run_generator_script generate_mlp_assessment
               Generate results/validation/benchmark_coverage_gaps.{md,json}.
  validation-figures
    shift
  thermo-gating
               Generate results/validation/thermodynamic_gating_audit.{md,json}.
      exit 1
               Generate results/validation/validated_envelope.{md,json}.
  explain-formulation NAME [TARGET_TAG] [MINIMIZE_TAG]
    ;;
  campaign SPEC [OUTPUT_DIR]
               Run a shareable campaign spec and generate campaign artifacts.
    if [ "$#" -lt 1 ]; then
  summary      Generate results/validation/benchmark_summary.{md,json}.
  deep-research-audit
    fi
  notebook     Launch a Jupyter notebook server on port 8888.
  status       Show container and environment status.
  extrusion-follow-on-workbook)
    shift
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
  run_in_env "python scripts/generators/generate_matrix_experiment_intake_schema.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_matrix_family_coverage.py"
  run_in_env "python scripts/generators/generate_refinement_governance.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_mlp_assessment.py"
  run_in_env "python scripts/generators/generate_computational_gap_refinement_plan.py --output-dir results/validation --manifest-dir results/computational_gap_refinement"
  run_in_env "python scripts/generators/ingest_computational_gap_dft_results.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_literature_learning_loop.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_promotion_state.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_lane_validation.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_deviation_audit.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_validation_figures.py --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_family_validation_figures.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_validated_envelope_report.py"
  run_in_env "python scripts/generators/generate_thermodynamic_gating_audit.py"
  scientific_fast_lane
}

scientific_fast_lane() {
  run_in_env "python -m pytest -m scientific_regression tests/scientific tests/unit/test_budget_projection.py tests/unit/test_safety_and_flux.py tests/integration/test_recommendation_engine.py"
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

run_generator_script() {
  local generator_name="$1"
  shift || true
  if [ "$#" -eq 0 ]; then
    run_in_env "python scripts/generators/${generator_name}.py"
  else
    run_in_env "python scripts/generators/${generator_name}.py $*"
  fi
}

validation_figures_lane() {
  local commands=(
    "generate_family_lane_validation --output-dir results/validation"
    "generate_family_deviation_audit --output-dir results/validation"
    "generate_validation_figures --docs-asset-dir docs/assets"
    "generate_family_validation_figures --output-dir results/validation --docs-asset-dir docs/assets"
    "generate_validated_envelope_report"
  )
  local command
  for command in "${commands[@]}"; do
    run_generator_script $command
  done
}

run_generator_alias() {
  case "$1" in
    targets-report) run_generator_script generate_benchmark_targets ;;
    matrix-deltas) run_generator_script generate_matrix_benchmark_deltas ;;
    matrix-assertions) run_generator_script generate_matrix_benchmark_assertions ;;
    matrix-evidence) run_generator_script generate_matrix_benchmark_evidence ;;
    matrix-readiness) run_generator_script generate_matrix_promotion_readiness ;;
    matrix-promotion-contract) run_generator_script generate_matrix_promotion_contract --output-dir results/validation ;;
    matrix-closure-audit) run_generator_script generate_matrix_observable_closure_audit --output-dir results/validation ;;
    experiment-intake-schema) run_generator_script generate_matrix_experiment_intake_schema --output-dir results/validation ;;
    literature-learning-loop) run_generator_script generate_literature_learning_loop --output-dir results/validation ;;
    family-ingestion-plan) run_generator_script generate_family_ingestion_plan --output-dir results/validation ;;
    family-promotion-state) run_generator_script generate_family_promotion_state --output-dir results/validation ;;
    matrix-family-coverage) run_generator_script generate_matrix_family_coverage ;;
    refinement-governance) run_generator_script generate_refinement_governance --output-dir results/validation ;;
    mlp-assessment) run_generator_script generate_mlp_assessment ;;
    coverage-gaps) run_generator_script generate_benchmark_coverage_gaps ;;
    thermo-gating) run_generator_script generate_thermodynamic_gating_audit ;;
    validated-envelope) run_generator_script generate_validated_envelope_report ;;
    index) run_generator_script generate_benchmark_index ;;
    summary) run_generator_script generate_benchmark_summary ;;
  esac
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
  targets-report|matrix-deltas|matrix-assertions|matrix-evidence|matrix-readiness|matrix-promotion-contract|matrix-closure-audit|experiment-intake-schema|literature-learning-loop|family-ingestion-plan|family-promotion-state|matrix-family-coverage|refinement-governance|mlp-assessment|coverage-gaps|thermo-gating|validated-envelope|index|summary)
    run_generator_alias "$cmd"
    ;;
  compare-experiment)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh compare-experiment INTAKE_FILE" >&2
      exit 1
    fi
    run_in_env "python scripts/generators/compare_matrix_experiment_intake.py --experiment '$1' --output-dir results/validation"
    ;;
  extrusion-closure-workbook)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh extrusion-closure-workbook WORKBOOK_FILE" >&2
      exit 1
    fi
    run_in_env "python scripts/generators/process_extrusion_external_closure_workbook.py --workbook '$1' --output-dir results/validation"
    ;;
  extrusion-follow-on-workbook)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh extrusion-follow-on-workbook WORKBOOK_FILE" >&2
      exit 1
    fi
    run_in_env "python scripts/generators/process_extrusion_disulfide_follow_on_workbook.py --workbook '$1' --output-dir results/validation"
    ;;
  extrusion-diagnostic-examples)
    run_generator_script generate_extrusion_diagnostic_examples
    ;;
  computational-gap-refinement-plan)
    run_generator_script generate_computational_gap_refinement_plan --output-dir results/validation --manifest-dir results/computational_gap_refinement
    ;;
  computational-gap-xtb)
    shift
    if [ "$#" -ge 1 ]; then
      run_in_env "python scripts/run_computational_gap_xtb.py --target '$1' --execute"
    else
      run_in_env "python scripts/run_computational_gap_xtb.py --execute"
    fi
    ;;
  computational-gap-dft)
    shift
    if [ "$#" -ge 1 ]; then
      run_in_env "exec python scripts/run_computational_gap_dft.py --target '$1' --execute"
    else
      run_in_env "exec python scripts/run_computational_gap_dft.py --execute"
    fi
    ;;
  computational-gap-dft-ingest)
    run_generator_script ingest_computational_gap_dft_results --output-dir results/validation
    ;;
  computational-gap-dft-promote)
    run_generator_script promote_computational_gap_dft_results --output-dir results/validation
    ;;
  matrix-branch-deltas)
    shift
    run_in_env "python scripts/compare_matrix_benchmark_branches.py --base-ref '${1:-main}'"
    ;;
  validation-figures)
    validation_figures_lane
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
  deep-research-audit)
    run_in_env "python scripts/deep_research_tracker.py"
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