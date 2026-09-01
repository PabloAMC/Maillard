# Contributing to Maillard

This file covers project layout, conventions for developers and AI agents, and how to
run tests. For scientist-facing usage, start with [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md).

---

## Project Layout

```
src/            Runtime package — kinetics, retention, headspace,
                recommend/optimize, reports.
                Import as: from src.<module> import ...
                Never mutate sys.path.

scripts/        CLI entrypoints and one-shot research scripts.
                Scientist-facing CLIs route through src/usability_reports.py.

data/           Curated inputs only; read-only at runtime (scripts/ci/data_readonly_gate.py).
                Generated map of every file: data/README.md. Sub-trees: keys/ (compound and paper
                identity, generated from seeds), schemas/ (JSON Schemas enforced by
                scripts/ci/schema_gate.py), species/, lit/, benchmarks/, protocols/.
                Examples live in docs/examples/, intake templates in docs/templates/, test
                fixtures in tests/fixtures/. Ingestion workflow: data/lit/README.md.
                Restructure record: tasks/data_restructure_plan.md.

results/        Generated artifacts. .gitignore'd. Do not hand-edit.

tests/          unit/ (fast), scientific/ (regression), scripts/ (integration).

docs/           Human-facing documentation. (Wave S5, 2026-08-28: this map used to omit
                four real subtrees and name three files that no longer exist. Nine documents
                were folded into their nearest living home and deleted; see AUDIT.md.)
  guides/                  Onboarding: QUICKSTART (incl. the command reference), GLOSSARY.
  reference/               Deep technical: SCIENTIFIC_REFERENCE, VALIDATION_CONTRACT.
  protocols/               Benchmark intake specs and the PPI/SPI primary protocol.
  validation/              Frozen evidence: the directional-accuracy panel and report,
                           the isotope-topology dossier. Read-only records of measurement.
  notebooks/               Executable walkthroughs.
  assets/                  Images and diagrams.
  slr_benchmark_evaluation.md  The systematic literature review two data registries cite.

data/Gemini_Deep_Research/  All literature synthesis reports (Elicit, SLR, kinetics).

tasks/          todo.md (active roadmap), lessons.md (process lessons).
```

---

## Execution Environment

**Always run code inside the Docker container with the `maillard` conda env (Python 3.12).**

Host Python is for editing and static analysis only. Never run `pytest`, `python`, or `pip`
directly on the host.

```bash
./scripts/docker_maillard.sh up                # boot container
./scripts/docker_maillard.sh bootstrap         # install deps (first time)
./scripts/docker_maillard.sh run "<cmd>"       # arbitrary command in container
./scripts/docker_maillard.sh shell             # interactive shell
```

If the container is not running: `./scripts/docker_maillard.sh up && ./scripts/docker_maillard.sh bootstrap`

---

## Running Tests

```bash
./scripts/docker_maillard.sh run "pytest tests/unit -q"
./scripts/docker_maillard.sh run "pytest tests/scientific -q"
./scripts/docker_maillard.sh run "pytest tests/ -q"     # full suite
```

Markers: `regression`, `slow`, `scientific_regression`, `kinetics_validation` (see `pytest.ini`).

---

## Key Conventions

### Observable-first governance
Never promote a target to a tighter prior without a justifying
artifact in `results/validation/`. No write-back if the result only reproduces existing
surrogate uncertainty.

### Confidence tiers
Records carry `confidence_tier` / `provenance_tier` / `uncertainty_posture` labels. Preserve
these through the pipeline; never silently upgrade.

### No computed barriers
Every barrier is a measured literature value or an explicitly labelled surrogate. The
former xTB → DFT refinement lane was removed on 2026-08-30/09-01; `assert_no_dft_*()`
guards in `src/kinetic_core/` run at import.

### Data access
Paths: `from src import data_paths` (one constant per curated file; never a `"data/..."` string).
Loads: `data_access.load_json / load_yaml` (missing or malformed files raise `DataFileError`).
Compound names: `compound_keys.resolve(name)`; papers: `paper_keys.for_doi(doi)`. Adding a
curated file means: a `data_paths` constant, a line in `build_data_readme.py`, and (for
benchmarks) passing `data/schemas/benchmark.schema.json`.

### Naming
Avoid opaque shorthand (`p3`, `wave1`, `c4_c5`). Name scripts/artifacts after the scientific
job (e.g. `matrix_sigma_residual_derivation`).

### Calibration single-application rule
`HeadspaceModel.get_matrix_benchmark_headspace_factor()` already applies the matrix observable
factor — never multiply `calibration_observable_factor` again downstream.

### No synthetic closure
Internally constructed mixed-matrix benchmarks do not count as external promotion evidence.

---

## Common Pitfalls

| Pitfall | Rule |
|---------|------|
| LaTeX-backed plots | Failure must be explicit; no silent fallback |
| Deleting scripts | Confirm with user first — scripts are often invoked ad-hoc |

Full lessons list: [tasks/lessons.md](tasks/lessons.md)

---

## Workflow for Non-Trivial Changes

1. **Plan first** — write plan to `tasks/todo.md` with checkable items
2. **Check in** before starting implementation
3. **Track progress** — mark items complete as you go
4. **Run tests** — never mark complete without proving it works
5. **Capture lessons** — update `tasks/lessons.md` after any correction
