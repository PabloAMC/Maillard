# Maillard Reaction Computational Framework — Working Plan
## Status: SOTA Alignment & Production Scaling

---

## ⚡ ACTIVE EXECUTION & SOTA ROADMAP

The core Tier 0/1/2 pipeline is operational and the zero-DFT laptop pipeline is functional. The **FFT pathway bottleneck (Phase 16) has been resolved** via an H₂S-mediated pathway and pH-dependent rate scaling.

> [!IMPORTANT]
> **Strategic Pivot (2026-03-10):** After a deep audit, we determined that multiple high-impact, non-DFT improvements should be prioritized *before* investing in custom DFT dataset generation and ML potential training. Public pre-trained models (MACE-OFF24, AIMNet2) already cover organic sulfur chemistry with ~0.25 kcal/mol accuracy. See `strategic_assessment.md` for full rationale.

> [!IMPORTANT]
> **Environment Note:** Use `.venv` for full test coverage. The `conda_env` is missing `mace-torch` and will skip MLP/MACE related tests (`tests/qm/test_mlp_*`).

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **🔴 1** | **Phase A: Literature Arrhenius Calibration** | 📋 | Replace heuristic barriers with published experimental (A, Ea) data |
| **🔴 2** | **Phase B: Public ML Potential Integration** | ✅ | Use MACE-OFF24/AIMNet2 instead of custom fine-tuning |
| **🔴 3** | **Phase C: Full Sensory Prediction Model** | ✅ | Expand OAV database from 5 → 30+ compounds; add psychophysical mixing |
| **🟠 4** | **Phase D: Headspace & Volatility Model** | ✅ | Convert matrix → air-phase concentrations for real OAV prediction |
| **🟠 5** | **Phase E: Sigmoid pH Model** | ✅ | Replace step-function pH multipliers with Henderson-Hasselbalch curves |
| **🟡 6** | **Phase F: AGE/Safety Scoring in Optimizer** | ✅ | Integrate toxic marker penalties into InverseDesigner |
| **🟡 7** | **Phase G: Concentration-Aware FAST Ranking** | ✅ | Add reactant-concentration weighting to FAST mode |
| **🟠 8** | **Phase H: Bayesian Formulation Optimization** | 📋 | Continuous optimization over (sugar, AA, pH, T, time) space |
| **🟠 9** | **Phase I: Matrix-Effect Corrections** | 📋 | Model protein/lipid volatile binding for realistic release prediction |
| **✅** | **Phase 12: Advanced Kinetic Features** | ✅ | pH, Temperature, and Thermo gating completed & verified |
| **✅** | **Phase 16: FFT Pathway Bottleneck** | ✅ | Investigated and resolved via H2S-mediated steps |
| **🔵** | **Phases 13–15: Custom DFT + MACE Training** | ⏳ | Deferred. Only needed if public ML potentials prove insufficient |
| **🔵** | **Phase 17: Web Dashboard** | ⏳ | Deferred. GUI for food scientists |
| **🔵** | **Phase 18: Experimental Validation Prep** | ⏳ | Deferred until pipeline is fully validated |

---

### [ACTIVE] Phase A: Literature Arrhenius Calibration `[🔴 CRITICAL | Diff: 4/10]`

> **Why:** `barrier_constants.py` uses coarse midpoint estimates from literature ranges (e.g., "12–20 kcal/mol → 15.0"). Published experimental Arrhenius parameters (A, Ea) exist for the major Maillard steps and would replace TST-with-heuristic-barriers with directly measured kinetics. This is a pure data-curation task.
>
> **Key References:** Martins & van Boekel 2003, Brands & van Boekel 2001, Yaylayan 1994, Hofmann & Schieberle 2000, Wedzicha 1984.

- [ ] **A.1 Literature Search & Data Extraction**: Search for published Arrhenius parameters (A, Ea) for each reaction family in `barrier_constants.py`. Record source, temperature range, pH, and substrate.
  - [ ] Schiff base condensation
  - [ ] Amadori/Heyns rearrangement
  - [ ] Strecker degradation (per amino acid)
  - [ ] 1,2- and 2,3-enolisation
  - [ ] Cysteine thermolysis / H₂S release
  - [ ] Retro-aldol fragmentation
- [ ] **A.2 Create `data/lit/arrhenius_params.yml`**: Structured YAML with per-family entries: `{family, A_value, A_unit, Ea_kcal, Ea_kj, temp_range_K, pH, substrate, source_doi}`.
- [ ] **A.3 Update `barrier_constants.py`**: Add a `get_arrhenius_params(family) → (A, Ea)` function that returns the literature values. Keep `get_barrier()` as a fallback for families without Arrhenius data.
- [ ] **A.4 Update `cantera_export.py`**: Use literature (A, Ea) in Cantera YAML instead of fixed `A=1e13` + TST-derived Ea. Modify `add_reaction()` to call `get_arrhenius_params()`.
- [ ] **A.5 Verification**: Run the Ribose+Cys+Leu benchmark system with new parameters. Compare concentration profiles before/after. Update `tests/integration/test_fft_bottleneck.py` thresholds if needed.
- [ ] **A.6 Tests**: Add `tests/unit/test_arrhenius_params.py` verifying data loading, fallback behavior, and that all families in the YAML have valid numeric entries.

---

### [DONE] Phase B: Public ML Potential Integration `[🔴 CRITICAL | Diff: 5/10]`

> **Why:** MACE-OFF24(M) is trained on ωB97M-D3BJ/def2-TZVPPD (essentially the same theory as our Phase 13 protocol) and already covers H, C, N, O, S. AIMNet2 achieves 1–2 kcal/mol RMSE on organic reactions. Using these as a drop-in for xTB barrier estimation gives near-DFT accuracy without any custom training.
>
> **This replaces old Phases 13–15 entirely** unless benchmarks reveal insufficient accuracy for specific reaction families.

- [x] **B.1 Install & Smoke Test**: Install `mace-off` (or `aimnet2`) in `.venv`. Verify a single-point energy calculation on a Maillard reactant (e.g., furfural).
- [x] **B.2 Create `src/mlp_barrier.py`**: New module wrapping `MACE-OFF24` (or AIMNet2) for barrier estimation:
  - `estimate_barrier(reactant_xyz, product_xyz) → float` using single-point energy differences.
  - `optimize_geometry(xyz_str) → xyz_str` for geometry refinement.
  - Fallback to `xtb_screener.py` if the ML potential fails.
- [x] **B.3 Benchmark**: Run against our 7 existing xTB path results (`data/geometries/xtb_inputs/`). Compare barrier estimates (xTB vs MACE-OFF24 vs literature). Document in `results/mlp_benchmark.json`.
- [x] **B.4 Integration**: Add `mace-off` as a method tier in `results_db.py` `find_barrier()`, positioned between `r2SCAN-3c` and `xtb`. Update `DFTRefiner` to call it before falling back to full DFT.
- [x] **B.5 Tests**: Add `tests/unit/test_mlp_barrier.py` with mocked energies to verify the wrapper handles edge cases (failed optimization, NaN energies, atom type errors).

---

### [DONE] Phase C: Full Sensory Prediction Model `[🔴 CRITICAL | Diff: 5/10]`

> **Why:** `sensory.py` hardcodes only 5 compounds. The YAML databases already contain ODT data for 16 desirable + 5 off-flavour + 9 toxic compounds. The current model also ignores aroma synergy (e.g., FFT + methional together create a superadditive meaty perception).

- [x] **C.1 Unified Sensory Database**: Merge data from `desirable_targets.yml`, `off_flavour_targets.yml`, and `toxic_markers.yml` into a single `SensoryDatabase` class that `SensoryPredictor` loads at init.
  - Each compound should have: `name`, `smiles`, `odour_threshold_ppm`, `descriptors[]`, `type` (desirable/off-flavour/toxic), `sensory_category` (from `sensory_tags.yml`).
  - Remove the hardcoded `SENSORY_DB` dict from `sensory.py`.
- [x] **C.2 Psychophysical Mixing Model**: Implement Stevens' power-law OAV addition:
  - `OAV_total(descriptor) = (Σ OAV_i^n)^(1/n)` where `n ≈ 0.5–0.8` captures sub-additive mixing (Feller model).
  - Add known synergy pairs (FFT + methional → "meaty" boost factor of 1.5×) from Hofmann 2000.
- [x] **C.3 Sensory Radar Chart Output**: Add a `to_radar_dict()` method that returns normalized 0–100 scores per descriptor category (meaty, roasted, beany, malty, earthy, sweet, sulfury). This is the data shape needed for Phase 17 (Web Dashboard).
- [x] **C.4 Integration**: Connect `SensoryPredictor` to `InverseDesigner.evaluate_all()` so formulation results include predicted sensory profiles, not just barrier-based scores.
- [x] **C.5 Tests**: Add `tests/unit/test_sensory_full.py`: verify that OAV calculation works for all 30+ compounds, that synergy pairs produce higher scores, and that the radar output sums correctly.

---

### [DONE] Phase D: Headspace & Volatility Partitioning `[🟠 HIGH | Diff: 6/10]`

> **Why:** The tool currently predicts concentrations in the liquid reaction matrix (kmol/m³ from Cantera). Food scientists care about **what the consumer smells**, which is the headspace (air-phase) concentration. These differ by orders of magnitude depending on the compound's vapor pressure and the matrix composition.
>
> **Key References:** Buttery 1969 (air-water partition), Roberts & Acree 1995 (modified Henry's law for food), Guichard 2002 (protein/lipid binding effects).

- [x] **D.1 Create `data/lit/henry_constants.yml`**: Literature air-water partition coefficients (Kaw) for ~30 key Maillard volatiles. Include temperature coefficients (ΔH_sol) for Clausius-Clapeyron extrapolation.
- [x] **D.2 Create `src/headspace.py`** with class `HeadspaceModel`:
  - `predict_headspace(matrix_conc_dict, temperature_k, fat_fraction, protein_fraction) → air_conc_dict`
  - Uses modified Henry's law: `C_air = Kaw(T) × C_matrix × f(fat, protein)`
  - Fat correction: hydrophobic volatiles partition into fat phase → lower headspace. Model via `Kaw_eff = Kaw / (1 + Kfat × fat_fraction)` where Kfat values are from Buttery.
  - Protein correction: polar volatiles bind to protein → lower headspace. Model via empirical `Kprot` values.
- [x] **D.3 Integration with Sensory Model**: Chain `Cantera output → HeadspaceModel → SensoryPredictor`. The OAV should use headspace concentrations, not matrix concentrations.
- [x] **D.4 Add matrix composition to `ReactionConditions`**: Extend `conditions.py` with `fat_fraction` and `protein_fraction` fields (default: 0.0 and 0.15 for typical plant protein matrix).
- [x] **D.5 Tests**: Add `tests/unit/test_headspace.py`: verify Kaw temperature scaling, fat-phase suppression of hydrophobic compounds (hexanal OAV should drop ~10× with 10% fat), and edge cases (0% fat, 0% protein).

---

### [DONE] Phase E: Sigmoid pH Model `[🟠 HIGH | Diff: 3/10]`

> **Why:** `conditions.py` uses hard step functions (`if pH < 6.0: return 5.0`) that create discontinuities. Real pH effects on Maillard kinetics follow sigmoid curves driven by amino group protonation (pKa-dependent). Van Boekel 2006 (*Biotechnology Advances*) provides explicit modeling approaches.

- [x] **E.1 Replace step functions**: In `conditions.py` `get_ph_multiplier()`, replace the `if/else` blocks with smooth sigmoid functions:
  - 1,2-enolisation (furan pathway): `mult = 1 + 4 / (1 + exp(k × (pH - 6.0)))` — peaks below pH 6.
  - 2,3-enolisation (pyrazine pathway): `mult = 1 + 4 / (1 + exp(-k × (pH - 6.0)))` — peaks above pH 6.
  - Schiff base: `mult = 1 + 2 × exp(−0.5 × ((pH − 5.5)/1.0)²)` — Gaussian peak at 5.5.
  - k ≈ 2.0 (steepness parameter), tunable to match van Boekel 2006 data.
- [x] **E.2 More reaction families**: Add pH multipliers for families not yet covered: `amadori`, `retro_aldol`, `cysteine` (H₂S release is faster at acidic pH), `dha` (alkaline-favored).
- [x] **E.3 Update Tests**: Modify `tests/unit/test_advanced_kinetics.py` to verify smooth behavior (no jumps at pH boundaries) and correct directionality (acidic → furans, alkaline → pyrazines).
- [x] **E.4 Cantera Propagation Check**: Verify that the new pH multipliers still propagate correctly through `cantera_export.py` `add_reaction()` (this should work automatically since Phase 16 fixed the disconnect).

---

### [DONE] Phase F: AGE/Safety Scoring in InverseDesigner `[✅ COMPLETE]`

> **Why:** The `toxic_markers.yml` database already contains CML, CEL, HMF, acrylamide, PhIP, and MeIQx with health risk classifications, but `InverseDesigner` doesn't penalize formulations that produce them. Plant-based scientists need to balance flavor *and* safety.

- [x] **F.1 Safety Score Calculation**: In `inverse_design.py`, add a safety scoring method.
- [x] **F.2 Update `FormulationResult`**: Add `safety_score` and `flagged_toxics`.
- [x] **F.3 Pareto Ranking**: Modify `evaluate_all()` to return formulations ranked by a combined score.
- [x] **F.4 Tests**: Added `tests/unit/test_safety_and_flux.py` for verification.

---

### [DONE] Phase G: Concentration-Aware FAST Ranking `[✅ COMPLETE]`

> **Why:** FAST mode (`recommend.py` `predict_from_steps`) ranked pathways purely by energetic span `exp(−Ea/RT)`, ignoring reactant concentrations. A system with 0.01 M cysteine and one with 1.0 M cysteine get identical rankings. This matters critically for plant-based formulation design.

- [x] **G.1 Boltzmann×Concentration Weighting**: In `predict_from_steps()`, modified the relaxation to weight by concentration.
- [x] **G.2 Bimolecular Correction**: Correct reactant concentrations propagation.
- [x] **G.3 Tracking Output**: Added `weighted_flux`.
- [x] **G.4 Tests**: Verified in `tests/unit/test_safety_and_flux.py`.

---

### [DONE] Phase H: Bayesian Formulation Optimization `[✅ COMPLETE]`

> **Why:** The current `InverseDesigner` evaluates a static grid of 14 formulations. Real formulation design requires optimizing over a continuous space: (sugar type & ratio, amino acid type & ratio, pH, temperature, duration, water activity). Bayesian optimization efficiently explores this space.

- [x] **H.1 Create `src/bayesian_optimizer.py`**: Implement a `FormulationOptimizer` class:
  - Uses `optuna` for Bayesian optimization.
  - Search space: sugar concentration (0.01–1.0 M), amino acid concentrations (0.01–1.0 M each), pH (3–9), temperature (100–200°C), time (10–120 min), water activity (0.3–0.95).
  - Objective: maximize `target_score − λ × safety_penalty` (from InverseDesigner).
- [x] **H.2 Acquisition Function**: Use Expected Improvement (EI) to balance exploration/exploitation. Each evaluation runs `SmirksEngine → Recommender → SensoryPredictor` pipeline.
- [x] **H.3 CLI Script**: Create `scripts/optimize_formulation.py` with:
  - `--target-tag meaty` (sensory target)
  - `--minimize-tag beany` (off-flavour to suppress)
  - `--n-iterations 50` (number of BO iterations)
  - `--risk-aversion` (penalty weight for toxic markers)
- [x] **H.4 Output**: Console optimization trajectory and Pareto front top hit.
- [x] **H.5 Tests**: Add `tests/unit/test_bayesian_optimizer.py` using a mock objective function to verify the optimizer converges and respects parameter bounds.

---

### [DONE] Phase I: Matrix-Effect Corrections `[🟠 HIGH | Diff: 6/10]`

> **Why:** In real PBMAs, plant proteins (soy, pea, wheat gluten) bind volatiles differently than meat proteins. Lipid content is also different (plant oils vs animal fat). These matrix effects can change perceived aroma by 2–10× even at identical chemical concentrations.
>
> **Key References:** Guichard 2002, van Ruth 2001, Kinsella 1989 (binding constants for β-lactoglobulin model).

- [x] **I.1 Protein-Binding Model**: Create lookup table of binding constants (Kp) for key volatiles with common plant proteins. Model: `C_free = C_total / (1 + Kp × [protein])`.
  - Data for soy protein isolate, pea protein, wheat gluten (from Guichard 2002).
- [x] **I.2 Lipid-Phase Partitioning**: Extend the headspace model (Phase D) with lipid-water partition coefficients for hydrophobic volatiles.
- [x] **I.3 Matrix Composition in Conditions**: Add `protein_type` field (soy/pea/wheat/generic) to `ReactionConditions`, which selects the appropriate binding constants.
- [x] **I.4 Integration**: Feed matrix-corrected free concentrations to the headspace model, creating the chain: `Cantera total concs → Matrix correction → Headspace model → Sensory predictor`.
- [x] **I.5 Tests**: Verify that soy protein matrix suppresses hexanal OAV (protein binds aldehydes), and that high-fat formulations suppress polar volatile release.

---

### [DONE] Phase 12: Advanced Kinetic Features `[✅ COMPLETE]`

> **Why:** Standard TST barriers are static. Real Maillard chemistry is heavily pH-dependent (e.g., Schiff Base is optimal at pH 5.5) and influenced by solvent. These features improve accuracy without requiring expensive DFT calculations.
- [x] **12.1 Barrier Auto-Scaling**: Implement pH and temperature-dependent activation scaling (Kirkwood-Onsager approximation) in `src/kinetics.py`.
- [x] **12.2 Multi-Reactant Integration**: Expand `CanteraExporter.add_reaction` to handle higher-order > 2 reactant/product rules.
- [x] **12.3 Dynamic Thermo-Gating**: Use `JobackEstimator` to pre-prune highly endergonic (ΔG > 30 kcal/mol) reactions before ODE integration.
- [x] **12.4 Solvent-Dependent Scaling**: Add empirical barrier scaling based on solvent dielectric constants.
- [x] **12.5 Verification**: Added `tests/unit/test_advanced_kinetics.py` to ensure long-term correctness of pH and solvent scaling.

---

### [DEFERRED] Phases 13–15: Custom DFT Dataset & MACE Fine-Tuning `[🔵 POSTPONED]`

> **Why deferred:** Public ML potentials (MACE-OFF24, AIMNet2) already cover organic sulfur at ~0.25 kcal/mol accuracy. Custom fine-tuning should only be pursued if Phase B benchmarks reveal insufficient accuracy for specific Maillard reaction families.

<details>
<summary><b>Phase 13: Tier 2 Batch Execution (Dataset Gen)</b></summary>

- [ ] **13.1** Ensure all input geometries for the top 500 reactions exist in `data/geometries/xtb_inputs/`.
- [ ] **13.2** Run `src/dft_refiner.py` in batch mode to converge 500+ reactions using PySCF/Sella bridge.
- [ ] **13.3** Save outputs to `results/dft_tier2/*.json`.

</details>

<details>
<summary><b>Phase 14: MACE Fine-Tuning for Sulfur Chemistry</b></summary>

- [ ] **14.1** Extend `scripts/generate_mace_training_data.py` for sulfur species augmentation.
- [ ] **14.2** Add `fine_tune()` entrypoint to `src/mlp_optimizer.py`.
- [ ] **14.3** Fine-tune MACE on Phase 13 data.
- [ ] **14.4** Verify RMSD drift < 0.05 Å on sulfur species.

</details>

<details>
<summary><b>Phase 15: Δ-ML Network Scaling</b></summary>

- [ ] **15.1** Create `src/delta_ml.py` implementing KRR on (E_MACE, E_DFT) pairs.
- [ ] **15.2–15.5** Training script, DB integration, tests, DFTRefiner integration.

</details>

---

### [DEFERRED] Phase 17: Web Dashboard for Food Scientists `[🔵 POSTPONED | Diff: 5/10]`

- [ ] **17.1** Create `app/` directory; scaffold a Streamlit multi-page app.
- [ ] **17.2** **Page 1 — Forward Simulation:** Input form (sugars, amino acids, additives, pH, temp, duration). Output: ranked volatile table + AGE risk scores.
- [ ] **17.3** **Page 2 — Inverse Design:** Input desired sensory tag + precursors-to-minimize. Output: ranked formulation grid with safety scores.
- [ ] **17.4** Connect pages to `SmirksEngine` + `InverseDesigner` + `KineticsEngine` backends.

---

### [DEFERRED] Phase 18: Experimental Validation Preparation `[🔵 POSTPONED]`

- [ ] **18.1** Select 5–10 top-ranked novel formulations from computational output.
- [ ] **18.2** Generate lab protocol summary for a food chemistry partner lab.
- [ ] **18.3** Define GC-MS validation criteria (expected retention times from NIST WebBook).
- [ ] **18.4** Document comparison methodology: predicted vs. observed volatiles (MAE, R²).

---


---

<details>
<summary><b>Phase 16: FFT Pathway Bottleneck Resolution (Mar 2026) ✅</b></summary>

> **Goal:** Investigate and resolve the persistent zero-yield bottleneck for **2-furfurylthiol (FFT)**.

- [x] **Root Cause Identification**: Identified structural deadlock (H₂ requirement) and pH-blindness in ODE simulations.
- [x] **H2S-Mediated Mechanism**: Replaced H₂-dependent template with a robust 2-step H₂S-mediated pathway via a thiohemiacetal intermediate.
- [x] **Numerical Stability**: Decomposed trimolecular reactions into bimolecular steps to resolve CVode integration failure.
- [x] **pH-Kinetics Propagation**: Fixed the disconnect between `conditions.py` and `cantera_export.py`.
- [x] **Verification**: FFT yield confirmed > 0 in Ribose+Cys system; pH sensitivity validated against literature.

</details>

<details>
<summary><b>Phases 9-12: High-Fidelity Refinements & Bug Fixes (Mar 2026) ✅</b></summary>
- [x] **Phase 9: Explicit Solvation (CREST/QCG)**: Wrapper created, added to Pipeline, verified freeze_core logic.
- [x] **Phase 10: MLP Geometry Opt (MACE)**: Integrated MACE with chemical identity guard.
- [x] **Phase 11: Sella TS Search**: Custom PySCF↔ASE bridge implemented with open-shell Hessian support.
- [x] **Phase 12: Advanced Kinetic Features**: Developed dynamic thermo-gating, pH/Temperature barrier scaling.
- [x] **Phase 12b: React-TS Diffusion Model Integration**: Successfully wrapped and verified.
- [x] **Bug Fixes**: Cantera Test Fixtures, Recommender Water/H2 implicit availability, Lysine Scope Bug, MLP NameError, PySCF↔ASE conversion logic, and Resource Warn/Tempfiles audit.
</details>

<details>
<summary><b>Phases 15-24: Scaling & Advanced Infrastructure (Mar 2026) ✅</b></summary>
- [x] **Phase 24: NASA Polynomial Thermodynamics**: Joback group additivity for accurate reverse rates.
- [x] **Phase 23: Barrier Source Unification**: ResultsDB querying unified across Cantera/Inverse.
- [x] **Phase 21: SmirksEngine → Cantera Bridge**: True zero-DFT laptop pipeline via exact barrier matching or constant fallback.
- [x] **Phase 18: Regression Gate**: `tests/test_regression.py` implemented.
- [x] **Phase 17: GC-MS Output**: Comparison metrics against lit.
- [x] **Phase 16: Structured Results Database**: Built SQLite DB tracking kinetics output and calculation provenance.
- [x] **Phase 15: Temperature Ramp**: `kinetics.py` properly dynamic.
</details>

<details>
<summary><b>Phases 6-8: Expansion & Real-World Utility ✅</b></summary>
- [x] Expanded formulation grid, lipid-Maillard synergy integration.
- [x] Implemented core PBMA precursors, inverse design mode, lysine budget & DHA tracking.
- [x] Automated Pathway Generation (SMIRKS rules formalized).
</details>

<details>
<summary><b>Phases 0-5: Infrastructure & Screening Prototype ✅</b></summary>
- [x] Base mechanism scaffolding, Multi-Reactant integration, 3-tier computational roadmap formulated.
- [x] xTB Pathway Screening parallelised.
- [x] Generative logic & core knowledge encoded.
</details>
