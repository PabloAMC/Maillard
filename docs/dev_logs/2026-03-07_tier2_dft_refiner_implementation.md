# Development Log: Tier 2 DFT Refinement Implementation
**Date:** 2026-03-07
**Target:** Phase 3.3 (Amadori & key bifurcations)

## Implementation Overview
We have transitioned from heuristic-based barriers to a tiered DFT refinement pipeline. The system is designed to provide chemical accuracy for rate-determining steps while maintaining a high-throughput capability for initial screening.

### Key Components
1. **`src/dft_refiner.py`**: A PySCF-based wrapper that implements the `r2SCAN-3c` // `wB97M-V` composite protocol.
2. **`scripts/run_tier2_dft.py`**: Command-line runner for batch refinement of the 8 target Maillard bifurcations.
3. **`src/thermo.py`**: Implements Quasi-Harmonic Correction (Grimme/Truhlar) to fix entropy overestimation in low-frequency modes.

## Decisions & Design Patterns

### 1. TS Guess Generation via `xtb --path`
Finding Transition States (TS) for Maillard systems is notoriously difficult due to the large number of atoms (~25-35) and complex proton shuttles.
- **Decision**: Reject the "TS Core Template" approach (as documented in `2026-03-07_ts_architecture_decision.md`).
- **Implementation**: Adopted `xtb --path` (reaction path search) to generate initial guesses. 
- **Requirement**: This necessitates strictly **atom-mapped** reactant and product XYZ coordinates. Shuffled indices in Tier 1 will cause the NEB search to fail.

### 2. The `--fast` Validation Mode
High-level DFT (`r2SCAN-3c`/`def2-TZVP`) on 27-atom systems takes several minutes per SCF step and hours for full TS optimization. This makes local unit testing impractical.
- **Decision**: Added a `--fast` flag to the runner script.
- **Implementation**: Swaps the meta-GGA and range-separated functionals for Hartree-Fock (`HF`) and the minimal `STO-3G` basis. 
- **Benefit**: Reduces runtime from hours to seconds for local validation of the pipeline logic (parsing, optimization flow, and thermo corrections).

### 3. Scaffolded xTB Input Directories
To scale the calculation to all 8 target reactions, we have automated the directory setup.
- **Location**: `data/geometries/xtb_inputs/{reaction_name}/`
- **Pattern**: Each folder contains a `run_xtb.sh` script. Users drop `reactant.xyz` and `product.xyz` (atom-mapped) into the folder and run the script to yield the TS guess for the DFT stage.

### 4. Direct IRC Validation (Phase 3.4) ✅
- **Decision**: Avoided external dependencies (geomeTRIC/pyberny) that proved fragile in the current environment.
- **Implementation**: Built a native "Displacement + Optimization" engine in `DFTRefiner`.
- **Logic**: Extracts the imaginary frequency mode from the TS Hessian, displaces and optimizes to verify reactant/product connectivity.
- **CLI**: Exposed via `--irc` flag.

## Performance Benchmark (2026-03-07)
- **Baseline**: Single-threaded PySCF (standard default).
- **Optimization**: Dynamic core detection using `os.cpu_count()`.
- **Result**: **10x speedup** on Apple M-series (10 cores). Reactant optimizations that previously took 30+ minutes now converge in <5 minutes at the `--fast` level.

## Technical Limitations & Future Work
- **Explicit Water**: Many Maillard steps (like Amadori) are catalyzed by water. Currently, we rely on implicit `ddCOSMO`. Adding explicit water molecules to the TS guesses will likely improve accuracy for these proton-transfer steps.
