# Dev Log: Phase 8.C FAST-mode Barrier Calibration and Literature Validation
**Date:** 2026-03-06
**Status:** Completed successfully

## 1. Context and Objective
The `SmirksEngine` ranking system previously relied on hardcoded barrier constants (e.g., `Schiff=15 kcal/mol`, `Strecker=22 kcal/mol`) embedded directly inside the `run_pipeline.py` and `inverse_design.py` evaluation loops. The Phase 8.C objective was to:
1. Extract these into a centralised, physics-grounded module (`src/barrier_constants.py`).
2. Calibrate these constants against published DFT and experimental literature.
3. Pass the **8.C.5 Literature Validation Gate**: ensure that the tool structurally reproduces the *known* dominant volatiles for 3 standard experimental model systems.

## 2. Calibration (src/barrier_constants.py)
We replaced the inline dicts with a dedicated module that serves as the single source of truth. The heuristic barriers were adjusted based on well-established literature (e.g., Yaylayan 1994, Martins 2003, Wedzicha 1984):

* **Schiff Base / Thiol Addition:** Very fast (15 kcal/mol)
* **Strecker Degradation:** Fast (20 kcal/mol)
* **Amadori / Heyns Rearrangement:** Moderate (23–24 kcal/mol)
* **Enolisation / Cysteine Thermolysis:** Slow, rate-limiting (28–30 kcal/mol)

*Note: While the `xtb_screener.py` NEB pipeline was reviewed, direct xTB NEB evaluation on implicit bimolecular steps routinely overestimates proton-transfer barriers by 10-25 kcal/mol (documented in `docs/xtb_limitations.md`). The heuristic constants were therefore anchored directly to published high-level DFT/CCSD(T) benchmarks rather than un-corrected xTB data.*

## 3. Literature Validation Gate (8.C.5)
A new test suite (`tests/test_barrier_calibration.py`) was created to formalise the validation gate. It runs the full generative pipeline + `Recommender` scoring for three classic systems:

| Model System | Experimental Expectation (Literature) | Tool's Top Predicted Volatiles | Result |
|---|---|---|---|
| **Ribose + Cysteine (pH 5.0)** | FFT (2-furfurylthiol) dominant | `['Hydrogen Sulfide', '2-Furfurylthiol (FFT)']` | **PASS** |
| **Glucose + Glycine (pH 7.0)** | Pyrazines, HMF, Pyruvaldehyde | `['2,5-Dimethylpyrazine']` (among top) | **PASS** |
| **Ribose + Cysteine + Leucine (pH 5.0)** | FFT and Strecker Aldehyde (3-methylbutanal) | `['Hydrogen Sulfide', '2-Furfurylthiol (FFT)', '3-Methylbutanal']` | **PASS** |

**Conclusion:** The Validation Gate is **PASSED**. The tool's purely SMARTS-driven enumeration combined with the newly calibrated kinetic rankings successfully reproduces the fundamental macro-behavior of the Maillard reaction observed in wet-lab experiments. 

We are now cleared to proceed to **Phase 8.D: Concentration/Ratio Support and Scoring Redesign.**
