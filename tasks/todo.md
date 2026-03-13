# Maillard Reactant Framework Improvement Plan

## Checklist

- [ ] **Phase 1: Validation Foundation**
    - [ ] Create `data/benchmarks/cys_ribose_150C_Mottram1994.json` (using feedback data) [Fix 1]
    - [ ] Create `data/benchmarks/cys_glucose_150C_Farmer1999.json` (using feedback data) [Fix 1]
    - [ ] Implement `tests/scientific/test_benchmarks.py` with `TOLERANCE` logic [Fix 2]
    - [ ] Wait for user to provide/confirm additional benchmark JSONs
    - [ ] Run benchmarks and document initial accuracy gap
- [ ] **Phase 2: Core Kinetic & Matrix Refinement**
    - [ ] Implement `src/matrix_correction.py` with `ProteinType` Enum and accessibility constants [Fix 7]
    - [ ] Update `src/conditions.py` (Arrhenius, ionization, Labuza water activity) [Fix 4, 5]
    - [ ] Update `src/recommend.py` to use `matrix_correction` and `results_db` cache [Fix 8/6]
    - [ ] Modify `scripts/run_pipeline.py` and `scripts/optimize_formulation.py` for `--protein-type` [Fix 3]
- [ ] **Phase 3: Sensory & Safety Enhancements**
    - [ ] Create `src/safety.py` with Knol 2009 Acrylamide model [Fix 12]
    - [ ] Update `src/sensory.py` with `ODT_BINDING_CORRECTION` and `MASKING_MATRIX` [Fix 9, 10]
    - [ ] Implement `export_qda_profile()` in `src/sensory.py` [Fix 11]
- [ ] **Phase 4: Lipid Oxidation & Crosstalk**
    - [ ] Implement `src/lipid_oxidation.py` with radical chain initiation logic [Fix 13]
    - [ ] Integrate `lipid_oxidation` precursors into `recommend.py` pipeline [Fix 14]
- [ ] **Phase 5: Optimization & Finalization**
    - [ ] Update `src/bayesian_optimizer.py` to use `predict_with_uncertainty` and penalize extrapolation [Fix 15]
    - [ ] Verify all benchmarks again and update `walkthrough.md`
- [ ] **Phase 6: SOTA Temporal FAST Mode**
    - [x] Implement `Arrhenius Integral` logic with Simpson's Rule in `src/recommend.py`
    - [x] Add CSV ramp parsing and interpolation utility
    - [x] Update `Recommender.predict_from_steps` to handle non-isothermal weights
    - [x] Verify against Cantera-based ground truth benchmarks


## Review Section

[To be populated upon completion]
