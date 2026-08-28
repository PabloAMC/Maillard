# Maillard — Agent Instructions

## Communication
- **Always reply in English**, even if the user writes in another language. The user prefers English-only assistant responses.
- Be concise. Use markdown links for file references (workspace-relative paths).

## Project Snapshot
Computational screening framework for meat-like Maillard chemistry in plant-based matrices. Combines deterministic kinetic ODEs (SMIRKS-based reaction families), matrix-aware retention/headspace physics, and selective QM (xTB → DFT via Sella/PySCF) for barrier refinement.

- Mission, architecture, SMIRKS tier conventions: [README.md](README.md) (the "Architecture and
  computational methods" section; `docs/architecture.md` and `docs/reference/SMIRKS_SYSTEM.md`
  were folded into it and deleted on 2026-08-28, Wave S5)
- Validation contract & benchmark surface: [docs/reference/VALIDATION_CONTRACT.md](docs/reference/VALIDATION_CONTRACT.md)
- Selective-QM runbook: [docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md](docs/guides/COMPUTATIONAL_GAP_RUNBOOK.md)
- Active roadmap & lessons: [tasks/todo.md](tasks/todo.md), [tasks/lessons.md](tasks/lessons.md)

## Layout
- `src/` — runtime package (kinetics, retention, headspace, QM refiners, recommend/optimize, reports). Import as `from src.<module> import ...`; do not mutate `sys.path`.
- `scripts/` — CLI entrypoints and one-shot research scripts. Scientist-facing CLIs route through `src/usability_reports.py`.
- `data/` — curated inputs only (priors, geometries, lit, protocols). No generated artifacts; no Python execution code (curated pathways live in [src/curated_pathways.py](src/curated_pathways.py)).
- `results/` — generated artifacts. `.gitignore`d; do not hand-edit.
- `tests/` — `unit/` (fast), `scientific/` (regression), plus `scripts/` integration tests.
- `models/external/` — third-party ML checkpoints with `provenance.json`.

## Execution Environment
**ALWAYS run code inside the Docker container with the `maillard` conda env (Python 3.12).** No exceptions:
- Tests, scripts, validation, QM, benchmarks, smoke checks, and one-off `python -c '...'` invocations all go through `./scripts/docker_maillard.sh run "<cmd>"`.
- Host Python is for editing/static analysis only — never `pytest`, `python`, or `pip` directly on the host.
- If the container is not up yet, run `./scripts/docker_maillard.sh up && ./scripts/docker_maillard.sh bootstrap` first.
- Subagents and execution helpers must follow the same rule; wrap every command in `./scripts/docker_maillard.sh run "..."`.

```bash
./scripts/docker_maillard.sh up           # boot container
./scripts/docker_maillard.sh bootstrap    # install deps
./scripts/docker_maillard.sh run "<cmd>"  # arbitrary command in container
./scripts/docker_maillard.sh summary      # benchmark summary artifact
./scripts/docker_maillard.sh computational-gap-{xtb,dft-preflight,dft,dft-ingest,dft-promote} <target>
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
- **Observable-first governance**: never promote a target to a tighter prior or selective-DFT anchor without an artifact in `results/validation/` justifying it. No write-back if the result only reproduces existing surrogate uncertainty.
- **Confidence tiers**: priors carry `bounded_calibration` / `transferred_literature` / `surrogate_family` / `xtb_derived` labels. Preserve them through the pipeline; never silently upgrade.
- **xTB is a pathfinder, not a barrier authority.** Final barriers come from DFT (r2SCAN/def2-svp + ddCOSMO water via PySCF/Sella) or explicit literature surrogates.
- **TS guess gates** (defense-in-depth before multi-hour DFT): minimum pairwise RMS Δ between reactant/TS ≥ 0.3 Å; xTB path multi-frame files must read the **last** frame; check imag-mode score before Sella.
- **Naming**: avoid opaque shorthand (`p3`, `wave1`, `c4_c5`). Name scripts/artifacts after the scientific job (e.g. `computational_gap_dft_ingestion_report`).
- **Calibration single-application rule**: `HeadspaceModel.get_matrix_benchmark_headspace_factor()` already applies the matrix observable factor — never multiply `calibration_observable_factor` again downstream.
- **No synthetic closure**: internally constructed mixed-matrix benchmarks do not count as external promotion evidence.
- **Multi-fragment geometry**: embed and MMFF-optimize each fragment independently before combining with explicit spacing; use `MMFFAddDistanceConstraint` for pre-docking, never dummy bonds.
- **DFT checkpointing**: pass `checkpoint_dir` to `calculate_barrier()` so optimized geometries survive failed downstream SCFs.

## Pitfalls (see [tasks/lessons.md](tasks/lessons.md) for full list)
- Stale execution JSON: refresh DFT preflight before trusting per-target status.
- e3nn version: pin `e3nn==0.5.1` for MACE 0.3.x — other ML installs silently upgrade it.
- LaTeX-backed plots: failure must be explicit; no silent fallback.
- `Sella` + `JAX` optional imports: catch runtime failures (CPU feature gaps), not just `ImportError`.
- Before deleting any script, confirm with the user — scripts are often invoked ad-hoc via `docker_maillard.sh run`.

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
