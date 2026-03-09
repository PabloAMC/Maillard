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

## 🔵 Integration Tests (`tests/integration/`)
- **Cost**: 1–10 seconds.
- **Dependencies**: Cantera, xTB (if available).
- **Scope**: Multi-step simulations, CLI argument parsing, recommendation engine logic, and database persistence.
- **Command**: `pytest tests/integration/`

### 🟡 integration/ (Feature-level)
- `test_recommendation_engine.py`: Recommender/Inverse design pipeline.
- `test_cantera_sim.py`: Simulation stability and mass balance.
- `test_db_chemistry.py`: Verifies atom balance of reactions in the results database.

## 🔴 QM & Backend Tests (`tests/qm/`)
- **Cost**: High (Seconds to Minutes).
- **Dependencies**: PySCF, Sella, MACE, CREST.
- **Scope**: Geometry optimization, TS search stability, and explicit solvation cluster generation.
- **Command**: `pytest tests/qm/`

## 🔬 Research Benchmarks (`tests/benchmarks/`)
- **Cost**: Very High (HPC-scale).
- **Dependencies**: High-level DFT (PySCF).
- **Scope**: Phase 3 literature validation and cross-functional benchmarking (e.g. wB97M-V vs. Double-Hybrids). Many tests here are forward-looking placeholders for experimental data.
- **Command**: `pytest tests/benchmarks/`

---

## Markers
- `@pytest.mark.slow`: Deselect with `pytest -m "not slow"` for fast local iteration.
- `@pytest.mark.skipif(...)`: Automatically skips tests if external binaries (xTB, MACE) are missing.
