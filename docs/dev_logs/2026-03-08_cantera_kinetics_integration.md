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
