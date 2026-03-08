# Dev Log: Phase 12 — Cantera Microkinetics Integration
Date: 2026-03-08
Status: COMPLETED ✅

## Context
Phase 12 focuses on integrating Cantera for rigorous ODE-based microkinetic simulations. This bridges the gap between static activation barriers (DFT/xTB) and time-resolved concentration profiles, allowing direct comparison with experimental data (GC-MS).

## Implementation Details
1. **`src/cantera_export.py`**:
    - Generates Cantera YAML mechanism files from the `SmirksEngine` reaction network.
    - **Mass Balance Heuristic**: Implemented a "last-product-absort-remainder" algorithm to ensure atom conservation across all reactions (mandatory for Cantera 3.2+), even when individual species SMILES are incomplete or simplified.
    - **Reversible vs Irreversible**: Switched to reversible reactions (`<=>`) with consistent thermo stubs (`h0=0`, `s0=0`) to match the framework's current focus on kinetic barriers over equilibrium constants.
2. **`src/kinetics.py`**:
    - Added `simulate_network_cantera` which wraps the Cantera Python API.
    - Handles batch reactor simulations at constant pressure/temperature.
    - Automates mapping between auto-generated Cantera species names (`S_0`, `S_1`, ...) and the framework's SMILES-based identity.
3. **`scripts/run_cantera_cli.py`**:
    - Provided a CLI interface to run kinetics simulations directly from precursor sets.
    - Supports `--export` for mechanism files and `--results` for concentration CSVs.

## Verification
- **Test Suite**: `tests/test_cantera_integration.py` and `tests/test_cantera_cli.py`.
- **Total Tests**: 23/23 passing.
- **Key Scenarios Validated**:
    - Isothermal equilibration.
    - Time-temperature profiles.
    - Multiple pathways competition (e.g. Ribose+Cys competing pathways).
    - **Barrier Sensitivity**: Confirmed that a 5 kcal/mol ΔG‡ difference produces >100x rate difference at 150°C, correctly captured by the ODE solver before equilibration effects dominate.

## Lessons Learned
- **Cantera Stoichiometry**: Cantera 3.2+ is extremely strict about atom conservation. Even with `check-balance: False`, it often fails to load unbalanced reactions if the phase definition is sensitive. The remainder-based heuristic is the most robust way to handle the "leaky" pathways inherent in early-stage Maillard modeling.
- **Equilibration Masking**: In a system with identical thermo (`h0=0`, `s0=0`), all reactions eventually reach 50/50 equilibrium. To test barrier sensitivity, it is critical to evaluate concentrations at short time scales ($t \ll 1/k_{fast}$) where kinetics (Ea) governs the flux.

## Next Steps
- **Phase 15**: Temperature Ramp Modeling (extending `kinetics.py` to handle $T(t)$ profiles).
- **Phase 17**: GC-MS Comparison Output (formatting Cantera results for analytical overlay).

---

# Dev Log: Phase 20 — Heuristic DB Population (Deprecated)
Date: 2026-03-08
Status: DEPRECATED ❌ → Absorbed into Phase 21

## Decision Record

### What was attempted

Phase 20 sought to make Cantera simulations possible without DFT by pre-populating `ResultsDB` with the 17 literature-calibrated barriers from `barrier_constants.py`. The approach extracted `ElementaryStep` objects from `data/reactions/curated_pathways.py`, inserted them as SMILES-keyed rows with `method="literature_heuristic"`, and extended `ResultsDB.find_barrier()` to fall back to this method if DFT was absent.

### Why it was reverted

During review, a fundamental architectural flaw was identified:

1. **`ResultsDB` is a precision instrument for computed barriers.** It indexes by exact SMILES pairs, meaning it correctly stores "for *this specific DFT calculation* on *this particular transition state*, the barrier was X kcal/mol." Pre-populating it with generic family constants conflates computed and assumed data.

2. **Loss of generality.** Inserting 16 specific molecule pairs works only for those 16 reactions. If a user runs `ribose:0.1,glucose:0.1` the DB returns nothing, despite `barrier_constants.py` already having the correct generic family rate. The pre-population is simultaneously redundant (for covered cases) and incomplete (for any novel combination).

3. **`barrier_constants.py` already solves this problem.** There is no need to copy generic constants into a relational database.

### Correct architecture (implemented in Phase 21)

The simulation script queries `ResultsDB` first (for any real DFT/xTB calculation on that exact SMILES pair), then falls back to `barrier_constants.get_barrier(reaction_family)` inline. The DB grows organically. No migration script, no pre-population.

### Side effects of this exploration (kept)

- `curated_pathways.py` was audited and two structurally unbalanced reactions were corrected (Cysteine Degradation now yields `Acetaldehyde + CO2`, not `Pyruvaldehyde + CO2`; Thiol Addition now includes explicit H₂ reactant). These corrections are chemically correct and align the file with `smirks_engine.py`.
- `src/cantera_export.py` now enforces **strict mass balance** and raises `ValueError` on any unbalanced reaction, replacing an earlier "virtual ghost atom" hack. This is a net improvement to the codebase.
- `src/results_db.py` now includes `"hf"` and `"literature_heuristic"` in the `find_barrier` method preference list. This is harmless and future-proof.

