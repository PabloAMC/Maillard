#!/usr/bin/env bash

set -euo pipefail

# 2026-09-14: two runtimes behind one set of verbs. The container is optional. Since the QM lane
# left (2026-09-01) environment.yml is a plain conda env, so the same commands run NATIVELY in a
# local conda env called `maillard`. Selection: MAILLARD_RUNTIME=native|docker; unset, native is
# chosen when no `docker` binary is on PATH, docker otherwise. Every verb below is unchanged; only
# `run_in_env`, `up`, `bootstrap`, `shell` and `status` know which backend they are on.
if [ -n "${MAILLARD_RUNTIME:-}" ]; then
  RUNTIME="$MAILLARD_RUNTIME"
elif command -v docker >/dev/null 2>&1; then
  RUNTIME="docker"
else
  RUNTIME="native"
fi
case "$RUNTIME" in native|docker) ;; *) echo >&2 "MAILLARD_RUNTIME must be native or docker, got '$RUNTIME'"; exit 2 ;; esac

IMAGE_NAME="${MAILLARD_IMAGE_NAME:-condaforge/miniforge3}"
# 2026-09-03: the container follows the HOST architecture. Until this date it was pinned to
# linux/amd64, which on Apple Silicon runs under Rosetta -- where forked worker processes serialise
# once numpy is imported (the envelope's six workers gave no speed-up; measurements in
# tasks/data_restructure_plan.md §7). Natively the tiers run 2.4x faster and the envelope ~8x. The
# amd64 environment stays reachable as MAILLARD_PLATFORM=linux/amd64 (container maillard_validation,
# volume maillard_conda); results agree across architectures to 7e-8 relative.
case "$(uname -m)" in
  arm64|aarch64) HOST_PLATFORM="linux/arm64" ;;
  *)             HOST_PLATFORM="linux/amd64" ;;
esac
PLATFORM="${MAILLARD_PLATFORM:-$HOST_PLATFORM}"
if [ "$PLATFORM" = "linux/amd64" ]; then
  CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_validation}"
  CONDA_VOLUME="${MAILLARD_CONDA_VOLUME:-maillard_conda}"
else
  CONTAINER_NAME="${MAILLARD_CONTAINER_NAME:-maillard_arm64}"
  CONDA_VOLUME="${MAILLARD_CONDA_VOLUME:-maillard_conda_arm64}"
fi
JUPYTER_PORT="${MAILLARD_JUPYTER_PORT:-8888}"
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

# NOTE: bash reads this file incrementally. Do not edit it while a `run` is in flight -- on 2026-09-03 an
# edit made a running invocation jump into another alias's case arm.
usage() {
  cat <<'EOF'
Usage: ./scripts/docker_maillard.sh <command> [args...]

Single entrypoint for the Maillard workbench. Every command runs inside the
maillard conda env (Python 3.12), either natively (a local conda env; the
default when no `docker` binary is found, or MAILLARD_RUNTIME=native) or in
the validated container (MAILLARD_RUNTIME=docker, or docker on PATH).
(2026-09-03, retirement step B5: the legacy screening lane's commands --
quickstart pipeline runs, campaigns, ingest, the matrix/family report lanes --
are gone with it. The front door is scripts/maillard.py.)

== Environment lifecycle ==
  up                                 Boot the container (docker runtime; a no-op natively).
  bootstrap                          Create or update the maillard conda env from environment.yml.
  shell                              Open an interactive bash shell in the env.
  status                             Show runtime, env and interpreter.
  notebook                           Launch Jupyter on port 8888.

== Daily loop ==
  run "<command>"                    Run any command inside the env.
  pytest [PYTEST_ARGS...]            Run pytest inside the env (default: tests/).
  quickstart                         Print a two-arm spec and run `maillard compare` on it.
  compare SPEC [args...]             python scripts/maillard.py compare SPEC ...
  predict SPEC [args...]             python scripts/maillard.py predict SPEC ...
  explain COMPOUND                   python scripts/maillard.py explain COMPOUND
  rank [args...]                     python scripts/maillard.py rank ...
  score DOC [args...]                python scripts/maillard.py score DOC ... (your own measurements)

== Evidence artifacts (regenerate what README quotes) ==
  core-scores                        results/validation/core_panel_scores.{json,md}   (~15 s)
  coverage                           per-module line coverage of src/ and scripts/ over both tiers (~35 min)
  results-readme                     results/README.md (the generated map of every tracked artifact)
  wishlist                           results/validation/data_wishlist.{json,md}: what to measure next, what it unlocks
  core-directional                   results/validation/core_directional_scores.*     (~90 s)
  core-envelope [args...]            results/validation/core_prediction_uncertainty.* (~40 min at n=200)
  model-card                         Re-splice the README model card from the artifacts.
  model-card-check                   Fail if the README model card is stale.
  gates                              Run the six CI gates (incl. artifact_freshness_gate, ~1 min).
  data-readme                        Regenerate data/README.md (fails on undescribed files).
  keys                               Rebuild data/keys/compounds.yml and papers.yml.
  experiment-value-ranking [...]     Rank (benchmark, compound) pairs by VoI.
  deep-research-audit                Literature backlog audit.
  family-ingestion-plan | family-priority | family-next-action | matrix-family-coverage | scope-guard
                                     Literature-side registries and reports.

== Test lanes ==
  core                               Unit + scripts tests (fast).
  scientific                         Integration + scientific tests.
  all                                Everything, plus the gates.

See README.md and docs/guides/QUICKSTART.md for the scientist flow.
EOF
}

core_lane() {
  run_in_env "python -m pytest tests/unit -q"
}

scientific_lane() {
  run_in_env "python -m pytest tests/scientific -q"
}

gates_lane() {
  local gate
  for gate in citation_gate data_readonly_gate fit_target_gate holdout_guard schema_gate artifact_freshness_gate; do
    run_in_env "python scripts/ci/${gate}.py"
  done
}

quickstart_lane() {
  echo "[quickstart] 1/2 writing a two-arm spec -> results/quickstart/my_comparison.yml"
  run_in_env "mkdir -p results/quickstart && python scripts/maillard.py compare --template > results/quickstart/my_comparison.yml"
  echo "[quickstart] 2/2 comparing the two arms on the kinetic core"
  run_in_env "python scripts/maillard.py compare results/quickstart/my_comparison.yml --report results/quickstart/compare.html"
  cat <<'MSG'
[quickstart] Done. Edit results/quickstart/my_comparison.yml (precursors in mM, temp_C, time_min, ph, aw)
and re-run:  ./scripts/docker_maillard.sh compare results/quickstart/my_comparison.yml
Open results/quickstart/compare.html for the OAV chart, refusal cards and declared assumptions.
MSG
}

forward_args() {
  local out="$1"; shift
  local arg quoted
  for arg in "$@"; do
    printf -v quoted '%q' "$arg"
    out+=" $quoted"
  done
  echo "$out"
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

# ---- native backend -------------------------------------------------------------------------
native_conda_sh() {
  # The conda.sh of whichever conda is on PATH (miniforge, miniconda, anaconda), or the usual homes.
  local base
  if command -v conda >/dev/null 2>&1; then
    base="$(conda info --base 2>/dev/null || true)"
    if [ -n "$base" ] && [ -f "$base/etc/profile.d/conda.sh" ]; then echo "$base/etc/profile.d/conda.sh"; return 0; fi
  fi
  local cand
  for cand in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/homebrew/Caskroom/miniforge/base /opt/conda /usr/local/miniforge3; do
    if [ -f "$cand/etc/profile.d/conda.sh" ]; then echo "$cand/etc/profile.d/conda.sh"; return 0; fi
  done
  echo >&2 "[docker_maillard.sh] no conda found. Install miniforge (brew install --cask miniforge) or set MAILLARD_RUNTIME=docker."
  exit 2
}

native_env_exists() {
  local sh; sh="$(native_conda_sh)"
  bash -lc "source '$sh' && conda env list | awk '{print \$1}' | grep -qx '$ENV_NAME'"
}

native_run_in_env() {
  local sh; sh="$(native_conda_sh)"
  local command="$1"
  bash -lc "set -eo pipefail; source '$sh'; export MKL_INTERFACE_LAYER=LP64; set +u; conda activate '$ENV_NAME'; set -u; cd '$WORKSPACE_DIR'; $command"
}

native_bootstrap_env() {
  local sh; sh="$(native_conda_sh)"
  if native_env_exists; then
    bash -lc "set -eo pipefail; source '$sh'; conda env update -n '$ENV_NAME' --file '$WORKSPACE_DIR/environment.yml'"
  else
    bash -lc "set -eo pipefail; source '$sh'; conda env create -n '$ENV_NAME' --file '$WORKSPACE_DIR/environment.yml'"
  fi
}

# ---- docker backend -------------------------------------------------------------------------
container_exists() {
  docker ps -a --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"
}

container_running() {
  docker ps --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"
}

ensure_container() {
  if ! container_exists; then
    docker run -d \
      --platform "$PLATFORM" \
      --name "$CONTAINER_NAME" \
      --memory=10g \
      --memory-swap=10g \
      -p "$JUPYTER_PORT:8888" \
      -v "$CONDA_VOLUME:/opt/conda" \
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
  if [ "$RUNTIME" = "native" ]; then native_run_in_env "$1"; return; fi
  ensure_container
  local command="$1"
  # We use 'export MKL_INTERFACE_LAYER=LP64' and 'set +u' to bypass common Conda activation script bugs
  docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; export MKL_INTERFACE_LAYER=LP64; set +u; conda activate '$ENV_NAME'; set -u; cd '$WORKSPACE_MOUNT'; $command"
}

bootstrap_env() {
  if [ "$RUNTIME" = "native" ]; then native_bootstrap_env; return; fi
  ensure_container

  if ! env_exists; then
    docker exec "$CONTAINER_NAME" bash -lc "set -eo pipefail; source '$CONDA_SH'; conda create -n '$ENV_NAME' python=3.12 -y"
  fi

  run_in_env "conda env update -n '$ENV_NAME' --file environment.yml"
}

status() {
  echo "runtime: $RUNTIME"
  if [ "$RUNTIME" = "native" ]; then
    if native_env_exists; then
      run_in_env "python --version && python -c 'import rdkit, numpy, scipy; print(\"rdkit\", rdkit.__version__, \"numpy\", numpy.__version__, \"scipy\", scipy.__version__)'"
    else
      echo "Conda environment '$ENV_NAME' has not been created yet. Run: ./scripts/docker_maillard.sh bootstrap"
    fi
    return
  fi
  ensure_container
  docker ps -a --filter "name=$CONTAINER_NAME"
  if env_exists; then
    run_in_env "python --version && python -c 'import rdkit; print(rdkit.__version__)'"
  else
    echo "Conda environment '$ENV_NAME' has not been created yet. Run: ./scripts/docker_maillard.sh bootstrap"
  fi
}

cmd="${1:-}"
case "$cmd" in
  help|--help|-h)
    usage
    ;;
  up)
    if [ "$RUNTIME" = "native" ]; then echo "runtime: native (nothing to boot)"; else ensure_container; fi
    ;;
  bootstrap)
    bootstrap_env
    ;;
  shell)
    if [ "$RUNTIME" = "native" ]; then
      sh="$(native_conda_sh)"
      exec bash -lc "source '$sh'; set +u; conda activate '$ENV_NAME'; cd '$WORKSPACE_DIR'; exec bash -i"
    fi
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
  quickstart)
    quickstart_lane
    ;;
  compare|predict|explain|rank|rank-experiments|score|wishlist)
    verb="$cmd"; shift
    run_in_env "$(forward_args "python scripts/maillard.py $verb" "$@")"
    ;;
  core)
    core_lane
    ;;
  scientific)
    scientific_lane
    ;;
  all)
    core_lane
    scientific_lane
    gates_lane
    ;;
  coverage)
    # per-module line coverage of src/ and scripts/ over both tiers (wave generators omitted on purpose)
    run_in_env "coverage erase && coverage run -m pytest tests/unit tests/scientific -q -p no:cacheprovider && coverage combine -q && coverage report --sort=cover"
    ;;
  gates)
    gates_lane
    ;;
  core-scores)
    run_generator_script generate_core_panel_scores
    ;;
  core-directional)
    run_generator_script generate_core_directional_scores
    ;;
  core-envelope)
    shift
    run_in_env "$(forward_args "python scripts/generators/generate_core_prediction_uncertainty.py" "$@")"
    ;;
  model-card)
    run_generator_script generate_model_card
    ;;
  model-card-check)
    run_generator_script generate_model_card --check
    ;;
  wishlist)
    run_generator_script generate_data_wishlist
    ;;
  results-readme)
    run_generator_script build_results_readme
    ;;
  data-readme)
    run_generator_script build_data_readme
    ;;
  keys)
    run_generator_script build_compound_registry
    run_generator_script build_paper_registry
    ;;
  experiment-value-ranking)
    shift
    run_in_env "$(forward_args "python scripts/generators/generate_experiment_value_ranking.py" "$@")"
    ;;
  family-ingestion-plan) run_generator_script generate_family_ingestion_plan --output-dir results/validation ;;
  family-priority) run_generator_script generate_matrix_family_priority_ranking ;;
  family-next-action) run_generator_script generate_matrix_family_next_action ;;
  matrix-family-coverage) run_generator_script generate_matrix_family_coverage ;;
  scope-guard) run_generator_script generate_scope_gap_guard ;;
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
