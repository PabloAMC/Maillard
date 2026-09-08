# Contributing to Maillard

This file covers the project layout, the rules for developers and AI agents, and how to run
the tests and gates. For using the tool, start with [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md).
For what the model is and how well it does, start with
[docs/guides/INTRODUCTION.md](docs/guides/INTRODUCTION.md).

---

## Project layout

```
src/            Runtime package. Import as `from src.<module> import ...`; never mutate sys.path.
  kinetic_core/   The one engine: lanes (trunk, sulfur, acrylamide, lipid), their parameters,
                  the declared condition terms, the panel, scoring, the Monte-Carlo envelope,
                  the directional scorer and the fit-target ledger.
  comparative_cli.py, explain_compound.py, experiment_value.py, report_html.py, model_card.py
                  The front door's verbs (compare, predict, explain, rank, score, wishlist).
  data_paths.py, data_access.py, compound_keys.py, paper_keys.py
                  The only way to reach data/: one constant per curated file, loads that raise
                  on a missing or malformed file, names resolved through data/keys/.

scripts/        maillard.py (the front door); docker_maillard.sh (every command runs through it);
                generators/ (the artifacts under results/, the frozen fit waves listed in
                generators/WAVES.md, the figures and the generated documents); ci/ (the six gates).

data/           Curated inputs only, read-only at runtime (scripts/ci/data_readonly_gate.py).
                Map of every file: data/README.md (generated). keys/ (compound and paper identity),
                schemas/ (enforced by scripts/ci/schema_gate.py), species/, lit/ (dossiers, re-typed
                tables, registries), benchmarks/ (the panel; external_validation/ is never fitted),
                articles/ (PDFs on disk, not tracked). Ingestion workflow: data/lit/README.md.

results/        Generated artifacts. results/validation/ is tracked as frozen evidence (scorecard,
                envelope, directional scorecard, every fit report and pre-registration); do not
                hand-edit; the freshness gate checks it against the generators. Map: results/README.md
                (generated). results/legacy_lane/ is the archive of the retired lane.

tests/          unit/ (fast), scientific/ (the headline guards that pin the README's numbers to
                the artifacts, the frozen-wave hashes, the physics regressions), support.py.

docs/           guides/ (INTRODUCTION, REACTION_TREES, SOURCES, QUICKSTART, GLOSSARY),
                USING_THE_TOOL.md, reference/ (validation contract, fit/hold-out declaration),
                protocols/, validation/ (the directional claims panel), assets/ (figures),
                history/ (the retired lane's README and quick start, the August 2026 audit, old
                roadmaps).

tasks/          data_restructure_plan.md (the living record; section 7 is the backlog),
                test_audit.md, lessons.md.
```

---

## Execution environment

**Always run code inside the Docker container with the `maillard` conda env (Python 3.12).**
Host Python is for editing and static analysis only.

```bash
./scripts/docker_maillard.sh up                # boot the container
./scripts/docker_maillard.sh bootstrap         # build the env from environment.yml (first time)
./scripts/docker_maillard.sh run "<cmd>"       # any command in the container
./scripts/docker_maillard.sh shell             # interactive shell
./scripts/docker_maillard.sh gates             # the six gates
```

Dependencies: everything in `environment.yml` and `pyproject.toml` has an importer under
`src/`, `scripts/` or `tests/`. Add a dependency only together with its consumer.

---

## Running tests

```bash
./scripts/docker_maillard.sh run "pytest tests/unit -q"
./scripts/docker_maillard.sh run "pytest tests/scientific -q"
./scripts/docker_maillard.sh run "pytest tests/ -q"     # everything
./scripts/docker_maillard.sh gates
```

Commit before running the gates: the data read-only gate needs a clean tree, and the freshness
gate compares tracked artifacts with what the generators produce.

---

## Rules

### Fit on primary evidence, validate on levels
Rate constants, activation energies, fed-intermediate yields, conversions and within-study
ratios may be fitted. End-of-cook concentrations validate; they are never fitted. The
declaration is `docs/reference/FIT_HOLDOUT_DECLARATION.md`; `scripts/ci/fit_target_gate.py`
and `holdout_guard.py` enforce it.

### A re-calibration is a new wave, never an edit
The fit generators are frozen by hash (`scripts/generators/WAVES.md`). To change one: copy it
under the next wave id, write the pre-registration with the test it must pass, run it, freeze
it, rebuild the manifest, add the line to WAVES.md, and add a row to the history table in the
introduction. A refused wave is kept as a record, not deleted.

### No number in code without a source
Every constant is a measured literature value with its dossier anchor, a value from a frozen
fit report, or a declared assumption with its band. Nothing is computed from quantum chemistry;
`assert_no_dft_*()` guards run at import.

### Data access
Paths through `src/data_paths.py`, loads through `src/data_access.py`, names through
`compound_keys.resolve()` and `paper_keys.for_doi()`. Adding a curated file means a
`data_paths` constant, a line in `build_data_readme.py`, and, for a benchmark, passing
`data/schemas/benchmark.schema.json`. Nothing writes into `data/` at runtime.

### Headline numbers move together
The README's numbers are pinned by `tests/scientific/test_core_headline_guards.py`. A change
that moves a number regenerates the artifact, the model card (`docker_maillard.sh model-card`)
and the README in the same commit.

### Names say what things are
No wave shorthand in user-facing text (the tags live in file names, artifacts and the glossary's
Part 3). Scripts and artifacts are named after the scientific job.

---

## Workflow for a non-trivial change

1. Write the plan into the backlog (`tasks/data_restructure_plan.md`, section 7) with checkable items.
2. Run both test tiers, commit, run the gates.
3. Record the outcome in the same backlog entry; add a lesson to `tasks/lessons.md` after any correction.
