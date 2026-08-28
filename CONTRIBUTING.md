# Contributing to Maillard

This file covers project layout, conventions for developers and AI agents, and how to
run tests. For scientist-facing usage, start with [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md).

---

## Project Layout

```
src/            Runtime package — kinetics, retention, headspace, QM refiners,
                recommend/optimize, reports.
                Import as: from src.<module> import ...
                Never mutate sys.path.

scripts/        CLI entrypoints and one-shot research scripts.
                Scientist-facing CLIs route through src/usability_reports.py.

data/           Curated inputs only (priors, geometries, lit, protocols).
                No generated artifacts. One exception, stated because this line used to
                deny it: data/reactions/curated_pathways.py IS Python under data/, alongside
                src/curated_pathways.py (corrected 2026-08-28, Wave S5).
                For the literature database and ingestion, see [data/lit/README.md](data/lit/README.md).

results/        Generated artifacts. .gitignore'd. Do not hand-edit.

tests/          unit/ (fast), scientific/ (regression), scripts/ (integration).

docs/           Human-facing documentation. (Wave S5, 2026-08-28: this map used to omit
                four real subtrees and name three files that no longer exist. Nine documents
                were folded into their nearest living home and deleted; see AUDIT.md.)
  guides/                  Onboarding: QUICKSTART (incl. the command reference),
                           COMPUTATIONAL_GAP_RUNBOOK, GLOSSARY.
  reference/               Deep technical: SCIENTIFIC_REFERENCE, VALIDATION_CONTRACT.
  protocols/               Benchmark intake specs and the PPI/SPI primary protocol.
  validation/              Frozen evidence: the directional-accuracy panel and report,
                           the isotope-topology dossier. Read-only records of measurement.
  notebooks/               Executable walkthroughs.
  assets/                  Images and diagrams.
  slr_benchmark_evaluation.md  The systematic literature review two data registries cite.

data/Gemini_Deep_Research/  All literature synthesis reports (Elicit, SLR, kinetics).

tasks/          todo.md (active roadmap), lessons.md (process lessons).
models/external/ Third-party ML checkpoints with provenance.json.
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
Never promote a target to a tighter prior or selective-DFT anchor without a justifying
artifact in `results/validation/`. No write-back if the result only reproduces existing
surrogate uncertainty.

### Confidence tiers
Priors carry `bounded_calibration` / `transferred_literature` / `surrogate_family` /
`xtb_derived` labels. Preserve these through the pipeline; never silently upgrade.

### xTB is a pathfinder, not a barrier authority
Final barriers come from DFT (r2SCAN/def2-svp + ddCOSMO water via PySCF/Sella) or explicit
literature surrogates.

### TS guess gates (before multi-hour DFT)
- Minimum pairwise RMS Δ between reactant/TS ≥ 0.3 Å
- xTB path multi-frame files: read the **last** frame
- Check imag-mode score before launching Sella

### Naming
Avoid opaque shorthand (`p3`, `wave1`, `c4_c5`). Name scripts/artifacts after the scientific
job (e.g. `computational_gap_dft_ingestion_report`).

### Calibration single-application rule
`HeadspaceModel.get_matrix_benchmark_headspace_factor()` already applies the matrix observable
factor — never multiply `calibration_observable_factor` again downstream.

### No synthetic closure
Internally constructed mixed-matrix benchmarks do not count as external promotion evidence.

### Multi-fragment geometry
Embed and MMFF-optimize each fragment independently before combining with explicit spacing.
Use `MMFFAddDistanceConstraint` for pre-docking, never dummy bonds.

### DFT checkpointing
Pass `checkpoint_dir` to `calculate_barrier()` so optimised geometries survive failed
downstream SCFs.

---

## Common Pitfalls

| Pitfall | Rule |
|---------|------|
| Stale execution JSON | Refresh DFT preflight before trusting per-target status |
| e3nn version conflict | Pin `e3nn==0.5.1` for MACE 0.3.x |
| LaTeX-backed plots | Failure must be explicit; no silent fallback |
| Sella + JAX optional imports | Catch runtime failures (CPU feature gaps), not just `ImportError` |
| Deleting scripts | Confirm with user first — scripts are often invoked ad-hoc |

Full lessons list: [tasks/lessons.md](tasks/lessons.md)

---

## Workflow for Non-Trivial Changes

1. **Plan first** — write plan to `tasks/todo.md` with checkable items
2. **Check in** before starting implementation
3. **Track progress** — mark items complete as you go
4. **Run tests** — never mark complete without proving it works
5. **Capture lessons** — update `tasks/lessons.md` after any correction
