# Maillard — Agent Instructions

## Communication
- **Always reply in English**, even if the user writes in another language. The user prefers English-only assistant responses.
- Be concise. Use markdown links for file references (workspace-relative paths).

## Project Snapshot
Computational screening framework for meat-like Maillard chemistry in plant-based matrices. Combines deterministic kinetic ODEs (SMIRKS-based reaction families) with matrix-aware retention/headspace physics. The former selective-QM (xTB → DFT) lane was removed on 2026-08-30/09-01; barriers are measured literature values or labelled surrogates.

- Mission and architecture: [README.md](README.md); the retirement plan and the improvement backlog:
  [tasks/data_restructure_plan.md](tasks/data_restructure_plan.md)
- Validation contract & benchmark surface: [docs/reference/VALIDATION_CONTRACT.md](docs/reference/VALIDATION_CONTRACT.md)
- Active roadmap & lessons: [tasks/todo.md](tasks/todo.md), [tasks/lessons.md](tasks/lessons.md)

## Layout
- `src/` — runtime package: `src/kinetic_core/` (the ONE engine: lanes, parameters, panel, scoring, envelope, fit-target ledger), `comparative_cli.py` (the front door's verbs), `report_html.py`, `explain_compound.py`, `experiment_value.py` (the `rank` verb), `model_card.py`, the literature-side registries. Import as `from src.<module> import ...`; do not mutate `sys.path`. The legacy SMIRKS lane was deleted 2026-09-03 (retirement step B5); its README is `docs/history/`, its artifacts `results/legacy_lane/`.
- `scripts/` — `maillard.py` (compare / predict / explain / rank), `generators/` (the core artifacts, the fit waves, the literature ledgers, the registries), `ci/` (five gates), `deep_research_tracker.py`.
- `data/` — curated inputs only, READ-ONLY at runtime (`scripts/ci/data_readonly_gate.py`). Map: [data/README.md](data/README.md) (generated). Paths come from `src/data_paths.py`, loads from `src/data_access.py` (missing files raise), names resolve via `data/keys/` (`src/compound_keys.py`, `src/paper_keys.py`), benchmarks validate against `data/schemas/` (`scripts/ci/schema_gate.py`). Never write a `"data/..."` string or a `json.load` of a curated file in code. Restructure record: [tasks/data_restructure_plan.md](tasks/data_restructure_plan.md).
- `results/` — generated artifacts (`results/validation/` is tracked as frozen evidence). Do not hand-edit. Generators write here, never into `data/`.
- `tests/` — `unit/` (fast), `scientific/` (regression), plus `scripts/` integration tests.

## Execution Environment
**ALWAYS run code inside the Docker container with the `maillard` conda env (Python 3.12).** No exceptions:
- Tests, scripts, validation, benchmarks, smoke checks, and one-off `python -c '...'` invocations all go through `./scripts/docker_maillard.sh run "<cmd>"`.
- Host Python is for editing/static analysis only — never `pytest`, `python`, or `pip` directly on the host.
- If the container is not up yet, run `./scripts/docker_maillard.sh up && ./scripts/docker_maillard.sh bootstrap` first.
- Subagents and execution helpers must follow the same rule; wrap every command in `./scripts/docker_maillard.sh run "..."`.

```bash
./scripts/docker_maillard.sh up           # boot container
./scripts/docker_maillard.sh bootstrap    # install deps
./scripts/docker_maillard.sh run "<cmd>"  # arbitrary command in container
./scripts/docker_maillard.sh core-scores  # the core's panel scorecard
./scripts/docker_maillard.sh gates        # the five CI gates
```

## Tests
```bash
# Inside container (preferred):
./scripts/docker_maillard.sh run "pytest tests/unit -q"
./scripts/docker_maillard.sh run "pytest tests/scientific -q"
./scripts/docker_maillard.sh run "pytest tests/ -q"   # full
```
Markers: `regression`, `slow`, `scientific_regression`, `kinetics_validation` (see [pytest.ini](pytest.ini)).

## Project-Specific Conventions
- **Observable-first governance**: never promote a target to a tighter prior without an artifact in `results/validation/` justifying it. No write-back if the result only reproduces existing surrogate uncertainty.
- **Confidence tiers**: records carry `confidence_tier` / `provenance_tier` / `uncertainty_posture` labels. Preserve them through the pipeline; never silently upgrade.
- **No computed barriers.** Every barrier is a measured literature value or an explicitly labelled surrogate; `assert_no_dft_*()` guards run at import in `src/kinetic_core/`.
- **Naming**: avoid opaque shorthand (`p3`, `wave1`, `c4_c5`). Name scripts/artifacts after the scientific job (e.g. `matrix_sigma_residual_derivation`).
- **Calibration single-application rule**: `HeadspaceModel.get_matrix_benchmark_headspace_factor()` already applies the matrix observable factor — never multiply `calibration_observable_factor` again downstream.
- **No synthetic closure**: internally constructed mixed-matrix benchmarks do not count as external promotion evidence.

## Pitfalls (see [tasks/lessons.md](tasks/lessons.md) for full list)
- LaTeX-backed plots: failure must be explicit; no silent fallback.
- Before deleting any script, confirm with the user — scripts are often invoked ad-hoc via `docker_maillard.sh run`.
- A directory under `data/` that `.gitignore` hides is invisible to every audit and gate (this is how `data/qm` shipped 18 unsourced barriers for four months). `data/*` is ignored by default; whitelist explicitly.
- Tests must never write into `data/`; `tests/integration/test_matrix_calibration_loop.py` once left 102 calibration files there. Monkeypatch output dirs to `tmp_path`.
- Regenerating a file that has been hand-corrected since destroys the correction (`conditions.buffer`, Wave W values). Generators of frozen evidence are write-once (`write_holdout_bundles`).

## Workflow Orchestration

### 1. Plan Node Default
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- If something goes sideways, STOP and re-plan immediately – don't keep pushing
- Use plan mode for verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Subagent Strategy
- Use subagents liberally to keep main context window clean
- Offload research, exploration, and parallel analysis to subagents
- For complex problems, throw more compute at it via subagents
- One tack per subagent for focused execution

### 3. Self-Improvement Loop
- After ANY correction from the user: update `tasks/lessons.md` with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

### 4. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between main and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness

### 5. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution"
- Skip this for simple, obvious fixes – don't over-engineer
- Challenge your own work before presenting it

### 6. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding
- Point at logs, errors, failing tests – then resolve them
- Zero context switching required from the user
- Go fix failing CI tests without being told how

## Task Management

1. **Plan First**: Write plan to `tasks/todo.md` with checkable items
2. **Verify Plan**: Check in before starting implementation
3. **Track Progress**: Mark items complete as you go
4. **Explain Changes**: High-level summary at each step
5. **Document Results**: Add review section to `tasks/todo.md`
6. **Capture Lessons**: Update `tasks/lessons.md` after corrections

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Changes should only touch what's necessary. Avoid introducing bugs.
