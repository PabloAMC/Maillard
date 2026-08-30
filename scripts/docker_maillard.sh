#!/usr/bin/env bash

set -euo pipefail

CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_validation}"
IMAGE_NAME="${MAILLARD_IMAGE_NAME:-condaforge/miniforge3}"
WORKSPACE_DIR="${MAILLARD_WORKSPACE_DIR:-$PWD}"
WORKSPACE_MOUNT="/workspace"
CONDA_SH="/opt/conda/etc/profile.d/conda.sh"
ENV_NAME="maillard"

# Stage a host file into the workspace so it is visible inside the container.
# Echoes the path (workspace-relative) that the container should use.
stage_into_workspace() {
  local src="$1"
  # Expand a leading ~ if the caller passed a literal tilde.
  case "$src" in
    "~"|"~/"*) src="${HOME}${src:1}" ;;
  esac
  if [ ! -e "$src" ]; then
    echo >&2 "[docker_maillard.sh] ERROR: file not found on host: $src"
    echo >&2 "  Hint: list candidates with:  ls -lt ~/Downloads | grep -i react_ot"
    exit 2
  fi
  local abs_src
  abs_src="$(cd "$(dirname "$src")" && pwd)/$(basename "$src")"
  case "$abs_src" in
    "$WORKSPACE_DIR"/*)
      # Translate host path to the container mount.
      echo "$WORKSPACE_MOUNT${abs_src#$WORKSPACE_DIR}"
      return 0
      ;;
  esac
  local stage_dir="$WORKSPACE_DIR/.cache/react_ot_uploads"
  mkdir -p "$stage_dir"
  local dest="$stage_dir/$(basename "$abs_src")"
  cp "$abs_src" "$dest"
  echo >&2 "[docker_maillard.sh] staged $abs_src -> $dest (mounted as $WORKSPACE_MOUNT/.cache/react_ot_uploads/$(basename "$abs_src"))"
  echo "$WORKSPACE_MOUNT/.cache/react_ot_uploads/$(basename "$abs_src")"
}

usage() {
  cat <<'EOF'
Usage: ./scripts/docker_maillard.sh <command> [args...]

Single Docker entrypoint for the Maillard workbench. All commands run inside
the maillard conda env (Python 3.12) in the validated container.

══ Container lifecycle ══
  up                                 Boot the validated container (one-time).
  bootstrap                          Install the conda env + dependencies.
  shell                              Open an interactive bash shell in the env.
  status                             Show container + env status.
  notebook                           Launch Jupyter on port 8888.

══ Scientist quickstart ══
  help                               Print this help screen.
  quickstart                         End-to-end demo: run a pea-isolate baseline
                                     and a cysteine-enrichment comparison, then
                                     point at the resulting report bundles.

══ Scientist daily loop ══
  run "<command>"                    Run any command inside the env.
  pytest [PYTEST_ARGS...]            Run pytest inside the env (default: tests/).
  ingest --file RESULTS.csv [...]    Preview + (with --confirm) land GC-MS data
                                     through scripts/ingest_results.py.
  campaign SPEC [OUTPUT_DIR]         Run a shareable campaign spec.
  campaign --names "A,B" [args...]   Run a named comparison through run_campaign.
  explain-formulation NAME [TARGET]  Explain a formulation result.
                  [MINIMIZE]
  compare-experiment INTAKE          Support-delta artifacts for a matrix
                                     experiment intake payload.

══ Trust-loop artifacts (regenerate published evidence) ══
  summary                            Benchmark summary (md + json).
  index                              Benchmark index.
  validated-envelope                 Per-compound 90% envelope report.
  validation-figures                 All validation figures + report visual
                                     examples in one bundle.
  coverage-gaps                      Benchmark coverage gaps.
  thermo-gating                      Thermodynamic gating audit.
  literature-learning-loop           Literature learning loop summary.
  external-inventory                 External validation inventory.
  refinement-watchlist               Refinement watchlist.
  objective-progress                 Strategic objective progress.
  scope-guard                        Out-of-scope demotion guard.
  family-priority                    Matrix family priority ranking.
  family-implementation              Which of the 16 lanes actually carry reaction templates
                                     (derived by engine enumeration; backs the README table).
  family-next-action                 Per-family next-action recommendation.
  family-promotion-state             Family promotion state.
  family-ingestion-plan              Family ingestion plan.
  matrix-family-coverage             Matrix × family coverage.
  matrix-deltas                      Matrix benchmark deltas.
  matrix-assertions                  Matrix benchmark assertions.
  matrix-evidence                    Matrix benchmark evidence summary.
  matrix-readiness                   Matrix promotion readiness.
  matrix-promotion-contract          Matrix promotion contract.
  matrix-closure-audit               Matrix observable closure audit.
  matrix-calibration-refit [...]     Re-run the matrix observable recalibration.
  matrix-branch-deltas [BASE_REF]    Matrix branch delta report (default: main).
  refinement-governance              Selective mechanistic refinement governance.
  mlp-assessment                     MLP backend assessment.
  experiment-intake-schema           Matrix experiment intake schema.
  experiment-value-ranking [...]     Rank (benchmark, compound) pairs by VoI;
                                     emits results/validation/experiment_value_ranking.{md,json}.
  next-experiment [--top N] [...]    Top-N experiment recommendations: prints
                                     ranked list and writes per-candidate intake
                                     YAML + protocol Markdown under
                                     data/protocols/ and results/validation/experiment_requests/.
  targets-report                     Benchmark targets report.
  targets BENCHMARK_JSON [TYPE]      Per-benchmark target snapshot.
  hofmann                            Hofmann 1998 selectivity diagnostic.
  deep-research-audit                Deep-research backlog audit.

══ Test lanes ══
  core                               Unit + integration tests.
  scientific                         Full scientific generator + regression lane.
  scientific-fast                    Scientific regression tests only.
  stability                          Targeted kinetics + safety tests.
  kinetics-validation                Cantera + temperature-profile tests.

══ Extrusion observable ingestion ══
  extrusion-closure-workbook FILE    Convert an external-closure workbook.
  extrusion-follow-on-workbook FILE  Convert a 5.8 disulfide follow-on workbook.
  extrusion-diagnostic-examples      Run synthetic diagnostic example bundles.

See README.md and docs/guides/QUICKSTART.md for the recommended scientist
flow. Trust artifacts live under results/validation/.
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
  run_in_env "python scripts/generators/generate_literature_learning_loop.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_promotion_state.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_lane_validation.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_family_deviation_audit.py --output-dir results/validation"
  run_in_env "python scripts/generators/generate_validation_figures.py --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_family_validation_figures.py --output-dir results/validation --docs-asset-dir docs/assets"
  run_in_env "python scripts/generators/generate_report_visual_examples.py --output-dir results/validation/report_visual_examples --docs-asset-dir docs/assets"
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

run_computational_gap_job() {
  local script_name="$1"
  local output_suffix="$2"
  local extra_args="$3"
  local target="${4:-}"
  local prefix="${5:-}"
  local command="python scripts/run_${script_name}.py ${extra_args}"

  if [ -n "$target" ]; then
    command="$command --target '$target' --output 'results/computational_gap_refinement/${target}_${output_suffix}_execution.json'"
  fi

  run_in_env "${prefix}${command}"
}

validation_figures_lane() {
  local commands=(
    "generate_family_lane_validation --output-dir results/validation"
    "generate_family_deviation_audit --output-dir results/validation"
    "generate_validation_figures --docs-asset-dir docs/assets"
    "generate_family_validation_figures --output-dir results/validation --docs-asset-dir docs/assets"
    "generate_report_visual_examples --output-dir results/validation/report_visual_examples --docs-asset-dir docs/assets"
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
    coverage-gaps) run_generator_script generate_benchmark_coverage_gaps ;;
    thermo-gating) run_generator_script generate_thermodynamic_gating_audit ;;
    validated-envelope) run_generator_script generate_validated_envelope_report ;;
    index) run_generator_script generate_benchmark_index ;;
    summary) run_generator_script generate_benchmark_summary ;;
    objective-progress) run_generator_script generate_objective_progress ;;
    external-inventory) run_generator_script generate_external_validation_inventory --output-dir results/validation ;;
    family-priority) run_generator_script generate_matrix_family_priority_ranking ;;
    family-next-action) run_generator_script generate_matrix_family_next_action ;;
    scope-guard) run_generator_script generate_scope_gap_guard ;;
    family-implementation) run_generator_script generate_family_implementation_status ;;
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
      --memory=10g \
      --memory-swap=10g \
      -p 8888:8888 \
      -v maillard_conda:/opt/conda \
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

quickstart_lane() {
  echo "[quickstart] 1/2 baseline run -> results/quickstart/baseline/"
  run_in_env "python scripts/run_pipeline.py \
    --sugars ribose,glucose \
    --amino-acids cysteine,leucine \
    --ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
    --ph 5.5 --temp 105 --time-minutes 45 \
    --protein-type pea_iso \
    --target meaty --minimize beany \
    --report --output-dir results/quickstart/baseline"
  echo "[quickstart] 2/2 baseline vs cysteine-enrichment comparison -> results/quickstart/comparison/"
  run_in_env "python scripts/run_campaign.py \
    --names 'Soy/Pea Base (Untreated),Cysteine Enrichment (Basic)' \
    --ph 5.5 --temp 105 \
    --target-tag meaty --minimize-tag beany \
    --campaign-name 'Quickstart pea-isolate head-to-head' \
    --output-dir results/quickstart/comparison"
  cat <<'MSG'

[quickstart] Done. Inspect:
  - results/quickstart/baseline/report.md                   (per-compound table + confidence + glossary)
  - results/quickstart/baseline/compound_confidence_overlay.png
  - results/quickstart/comparison/comparison.md             (side-by-side + intervention waterfall)
  - results/quickstart/comparison/comparison_intervention_waterfall.png

Next: ./scripts/docker_maillard.sh help        # full command surface
      docs/guides/QUICKSTART.md                # the daily loop, including ingest
MSG
}

cmd="${1:-}"
case "$cmd" in
  help|--help|-h)
    usage
    ;;
  quickstart)
    quickstart_lane
    ;;
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
  targets-report|matrix-deltas|matrix-assertions|matrix-evidence|matrix-readiness|matrix-promotion-contract|matrix-closure-audit|experiment-intake-schema|literature-learning-loop|family-ingestion-plan|family-promotion-state|matrix-family-coverage|coverage-gaps|thermo-gating|validated-envelope|index|summary|objective-progress|external-inventory|family-priority|family-next-action|scope-guard|family-implementation)
    run_generator_alias "$cmd"
    ;;
  ingest)
    shift
    if [ "$#" -lt 1 ]; then
      echo "Usage: ./scripts/docker_maillard.sh ingest --file RESULTS.csv [INGEST_ARGS...]" >&2
      exit 1
    fi
    ingest_cmd="python scripts/ingest_results.py"
    for arg in "$@"; do
      printf -v quoted_arg '%q' "$arg"
      ingest_cmd+=" $quoted_arg"
    done
    run_in_env "$ingest_cmd"
    ;;
  matrix-calibration-refit)
    shift
    refit_cmd="python scripts/generators/generate_matrix_calibration_refit.py"
    for arg in "$@"; do
      printf -v quoted_arg '%q' "$arg"
      refit_cmd+=" $quoted_arg"
    done
    run_in_env "$refit_cmd"
    ;;
  experiment-value-ranking)
    shift
    rank_cmd="python scripts/generators/generate_experiment_value_ranking.py"
    for arg in "$@"; do
      printf -v quoted_arg '%q' "$arg"
      rank_cmd+=" $quoted_arg"
    done
    run_in_env "$rank_cmd"
    ;;
  next-experiment)
    shift
    next_cmd="python scripts/request_experiment.py"
    for arg in "$@"; do
      printf -v quoted_arg '%q' "$arg"
      next_cmd+=" $quoted_arg"
    done
    run_in_env "$next_cmd"
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
      echo "   or: ./scripts/docker_maillard.sh campaign --names \"A,B\" [extra run_campaign args]" >&2
      exit 1
    fi
    if [ "${1#--}" != "$1" ]; then
      cmd="python scripts/run_campaign.py"
      for arg in "$@"; do
        printf -v quoted_arg '%q' "$arg"
        cmd+=" $quoted_arg"
      done
      run_in_env "$cmd"
    elif [ "$#" -ge 2 ]; then
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