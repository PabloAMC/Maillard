# Maillard Reactant Framework Improvement Plan

## Checklist

- [x] **Phase 1: Validation Foundation**
    - [x] Create `data/benchmarks/cys_ribose_150C_Mottram1994.json`
    - [x] Create `data/benchmarks/cys_glucose_150C_Farmer1999.json`
    - [x] Implement `tests/scientific/test_benchmarks.py` with `TOLERANCE` logic
    - [ ] Wait for user to provide/confirm additional benchmark JSONs
    - [x] Run benchmarks and document initial accuracy gap
- [x] **Phase 2: Core Kinetic & Matrix Refinement**
    - [x] Implement `src/matrix_correction.py` with `ProteinType` Enum and accessibility constants
    - [x] Update `src/conditions.py` (Arrhenius, ionization, Labuza water activity)
    - [x] Update `src/recommend.py` to use `matrix_correction` and `results_db` cache
    - [x] Modify `scripts/run_pipeline.py` and `scripts/optimize_formulation.py` for `--protein-type`
- [x] **Phase 3: Sensory & Safety Enhancements**
    - [x] Create `src/safety.py` with Knol 2009 Acrylamide model
    - [x] Update `src/sensory.py` with `ODT_BINDING_CORRECTION` and `MASKING_MATRIX`
    - [x] Implement `export_qda_profile()` in `src/sensory.py`
- [x] **Phase 4: Lipid Oxidation & Crosstalk**
    - [x] Implement `src/lipid_oxidation.py` with radical chain initiation logic
    - [x] Integrate `lipid_oxidation` precursors into `recommend.py` pipeline
- [x] **Phase 5: Optimization & Finalization**
    - [x] Update `src/bayesian_optimizer.py` to use `predict_with_uncertainty` and penalize extrapolation
    - [x] Verify all benchmarks again and update `walkthrough.md`
- [ ] **Phase 6: SOTA Temporal FAST Mode**
    - [x] Implement `Arrhenius Integral` logic with Simpson's Rule in `src/recommend.py`
    - [x] Add CSV ramp parsing and interpolation utility
    - [x] Update `Recommender.predict_from_steps` to handle non-isothermal weights
    - [x] Verify against Cantera-based ground truth benchmarks
- [ ] **Phase 7: Scientific Gaps (Identified by Analysis)**
    - [ ] **Gap #1: Predictive Accuracy vs. Experimental Results** (Difficulty: High)
        - *Gap*: Current FAST mode over-predicts absolute yields (MAE > 6000 ppm) due to heuristic barriers and lack of thermodynamic gating.
        - *Impact*: Users cannot trust absolute concentrations for safety or regulatory compliance.
        - *Enhancement*: 
            - [ ] Implement **Thermodynamic Gating** (ΔG consistency) in `KineticsEngine`.
            - [ ] Create `scripts/calibrate_barriers.py` to automate offset fitting against benchmarks.
            - [ ] Implement **Phase Partitioning** (Henry's Law + Matrix Binding) in `src/headspace.py`.
    - [ ] **Polypeptide-Bound Reactivity** (Difficulty: High)
        - *Gap*: `SmirksEngine` models free amino acids; protein isolates are peptide-bound.
        - *Impact*: Overestimates pyrazine yield; underestimates meat-browning color.
        - *Enhancement*: Add rules for N-terminal & Lysine sidechain peptide reactivity to `SmirksEngine`.
    - [ ] **Polysaccharide & Oligomer Complexity** (Difficulty: Medium)
        - *Gap*: Precursors limited to glucose/ribose; ingredients have maltodextrins/starch.
        - *Impact*: Models miss different kinetics and steric profiles of larger sugars.
        - *Enhancement*: Implement "Sugar Polymer Fragmenter" for (α1->4) and (α1->6) linkages.
    - [ ] **Broad-Spectrum Phytochemical Intervention** (Difficulty: Medium)
        - *Gap*: Lacks mechanism for polyphenol scavenging of dicarbonyls.
        - *Impact*: Over-predicts HMF and Acrylamide by failing to model "quenching" capacity.
        - *Enhancement*: Add `Phytochemical_Capture` family for sequestration by phenolic rings.
    - [ ] **Flavor-Texture Bridge** (Difficulty: High)
        - *Gap*: DHA cross-linking is only heuristic; governs final product firmness.
        - *Impact*: Scientists cannot predict if maximizing flavor will ruin texture.
        - *Enhancement*: Couple DHA yields to a rheological approximation layer.
    - [ ] **High-Shear Processing Dynamics (Extrusion)** (Difficulty: Very High)
        - *Gap*: Static {pH, T, t} model misses mechanical radical generation/burial disruption.
        - *Impact*: Inability to capture "flash" flavor generation in high-moisture extrusion.
        - *Enhancement*: Add `shear_rate` parameter to `ReactionConditions`.


## Review Section

[To be populated upon completion]
