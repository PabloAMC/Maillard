# Maillard Reaction Framework — Test Architecture

To ensure both rapid development and scientific accuracy, the Maillard test suite is organized into hierarchical tiers based on computational cost and binary dependencies.

## 🟢 Unit Tests (`tests/unit/`)

- **Cost**: < 1 second.
- **Dependencies**: None (Pure Python/RDKit).
- **Scope**: Mathematical logic (`kinetics`, `thermo`), string processing (`precursors`), and core rule engines (`SmirksEngine`).
- **Command**: `pytest tests/unit/`

### 🟢 unit/ (Sub-second)

- `test_rdkit_logic.py`: Consolidates SMARTS/SMIRKS manipulation.
- `test_kinetics_math.py`: Pure mathematical checks for rates/half-lives.
- `test_data_integrity.py`: Verifies atom balance of curated reaction pathways.
- `test_precursor_resolver.py`: Name-to-SMILES fuzzy matching.

## 🔵 Integration Tests (`tests/integration/`) — EMPTY since retirement step B5b (2026-09-03); the two tiers are `tests/unit` and `tests/scientific`

- **Cost**: 1–10 seconds.
- **Dependencies**: Cantera, xTB (if available).
- **Scope**: Multi-step simulations, CLI argument parsing, recommendation engine logic, and database persistence.
- **Command**: `pytest tests/integration/`

## Stability Gate

- **Purpose**: Fast reproducible Tier 0/1 regression gate for the currently validated deterministic contract.
- **Command**: `./scripts/docker_maillard.sh stability`
- **Includes**: Arrhenius normalization, headspace partitioning, concentration-aware ranking, Cantera simulation stability, parser-safe serialization, and barrier-calibration smoke coverage.

### 🟡 integration/ (Feature-level)

- `test_recommendation_engine.py`: Recommender/Inverse design pipeline.
- `test_cantera_sim.py`: Simulation stability and mass balance.
- `test_db_chemistry.py`: Verifies atom balance of reactions in the results database.

## 🔴 QM & Backend Tests (`tests/qm/`)

- **Cost**: High (Seconds to Minutes).
- **Dependencies**: PySCF, Sella, MACE, CREST.
- **Scope**: Geometry optimization, TS search stability, and explicit solvation cluster generation.
- **Command**: `pytest tests/qm/`

> **`tests/benchmarks/` was removed on 2026-08-27 (Wave J2).** Its Phase 3 QM authority-lane
> tests were deleted on 2026-04-21; all that remained was an orphaned `_lane_policy.py`
> helper that nothing imported. The directory collected zero tests, so the skip scan that
> pointed at it had been reporting an empty lane's zeros as a clean lane ever since. If the
> Phase 3 literature-validation lane is revived, it needs new tests, not that scaffolding.

---

## Markers

- `@pytest.mark.slow`: Deselect with `pytest -m "not slow"` for fast local iteration.
- `@pytest.mark.skipif(...)`: Automatically skips tests if an external backend is unavailable or unsupported in the active environment.
- `xfail_strict = true` is set in `pytest.ini` (2026-08-27). An `xfail`-marked test that starts
  passing is a **failure**, not a silent `xpass`: closing a known gap forces you back here to
  promote the test into a real assertion.

## What counts as a legitimate skip

A skip is legitimate only when its condition names a **real external precondition** that the
repository cannot satisfy for itself — a missing binary (xTB, CREST), an uninstalled optional
backend (PySCF, Sella, MACE), a platform limitation, or a gitignored artifact that must be
generated locally. The reason string must name that precondition specifically.

A skip is **not** legitimate when:

- its condition is "the behaviour under test is broken" — that is a self-excusing skip, and it
  can never fail no matter how bad the code gets. Assert the current behaviour instead, or mark
  the test `xfail(strict=True)` so it goes red when the gap closes;
- its condition can never be true (for example, a file that is tracked in git). Assert the
  precondition so its loss is visible;
- it stands in for an unwritten assertion. Write the assertion or delete the test — a green
  test that checks nothing launders confidence.

Optional QM backend skips should probe the actual backend capability, not hard-code stale paths
such as `conda_env/bin/...`.
