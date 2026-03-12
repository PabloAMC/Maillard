# Maillard Reaction Computational Framework — Working Plan
## Status: Scientific Accuracy & SOTA Hardening

---

## ⚡ ACTIVE ROADMAP (Post-Audit 2026-03-11)

> [!IMPORTANT]
> **Architectural Review (2026-03-11):** A comprehensive code review identified 15 findings across the pipeline. 10 out of 12 critical tasks have been completed. The remaining tasks are listed below.

> [!IMPORTANT]
> **Environment Note:** We recommend using the unified conda environment (`environment.yml`) for full QM support (`CREST`/`xtb`).

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **🔥 P0** | **19: Lipid Oxidation Radical Network** | `[ ]` | Predict "beany" off-flavors natively |
| **🚀 P1** | **20: Toxicity Decoupling Interventions** | `[ ]` | BO-driven additive suggestions |
| **🔬 P2** | **21: Enzymatic Pre-Processing Layer** | `[ ]` | Initial matrix cleanup simulation |

---

### Phase 19: Lipid Oxidation Radical Network `[🔥 P0 | Scientific Depth]`

> **Goal:** Model the generation of hexanal/nonanal from PUFAs and their crosstalk with Maillard intermediates.

- [ ] **19.1** Define `LipidRadical` family in `src/smirks_engine.py` using β-scission SMIRKS.
- [ ] **19.2** Add lipid-specific barrier constants for hydroperoxide homolysis (lit-derived).
- [ ] **19.3** Implement `Crosstalk` templates (e.g., $H_2S$ + hexanal → alkylthiazoles).
- [ ] **19.4** Add validation test: verify hexanal generation from linoleic acid model.

### Phase 20: Toxicity Decoupling Interventions `[🚀 P1 | Utility]`

> **Goal:** Enhance the Bayesian Optimizer to suggest ingredients that lower toxicity without compromising flavor.

- [ ] **20.1** Add `InterventionParams` (e.g., calcium carbonate, rosemary extract) to `InverseDesigner`.
- [ ] **20.2** Map interventions to specific barrier modifiers in `mlp_barrier.py` or `kinetics.py`.
- [ ] **20.3** Update `FormulationOptimizer` to include these as searchable categorical variables.
- [ ] **20.4** Test: Verify BO finds a "low-acrylamide" variant of a target flavor profile.

### Phase 21: Enzymatic Pre-Processing Layer `[🔬 P2 | Workflow]`

> **Goal:** Simulate yeast/protease action on raw ingredients before the thermal simulation starts.

- [ ] **21.1** Create `src/pre_processor.py` with standard "Fermentation/Hydrolysis" profiles.
- [ ] **21.2** Implement `PreProcessor.apply(initial_concentrations)` to transform feedstocks.
- [ ] **21.3** Integrate into `scripts/optimize_formulation.py` as a pre-step flag.

---

## ✅ Archived: Completed Work (Post-Audit 2026-03-11)

<details>
<summary><b>R.4: NASA7 Polynomial Format Fix ✅</b></summary>

- Modified `get_nasa_coefficients` to produce two-range NASA7 with 14 coefficients.
</details>

<details>
<summary><b>R.8: Bayesian Optimizer Thread Safety ✅</b></summary>

- Ensured fresh state per trial in BO.
</details>

---

### [DEFERRED] Phases 13–15: Custom DFT Dataset & MACE Fine-Tuning `[🔵 POSTPONED]`

> **Why deferred:** Public ML potentials (MACE-OFF24, AIMNet2) already cover organic sulfur at ~0.25 kcal/mol accuracy. Custom fine-tuning should only be pursued if Phase B benchmarks reveal insufficient accuracy for specific Maillard reaction families.

<details>
<summary>Details</summary>

- [ ] **13.1** Ensure all input geometries for the top 500 reactions exist in `data/geometries/xtb_inputs/`.
- [ ] **13.2** Run `src/dft_refiner.py` in batch mode to converge 500+ reactions using PySCF/Sella bridge.
- [ ] **14.1-14.4** MACE Fine-Tuning for Sulfur Chemistry.
- [ ] **15.1-15.5** Δ-ML Network Scaling.

</details>

---

### [DEFERRED] Phase 17: Web Dashboard for Food Scientists `[🔵 POSTPONED]`

- [ ] Streamlit multi-page app with Forward Simulation and Inverse Design pages.

---

### [DEFERRED] Phase 18: Experimental Validation Preparation `[🔵 POSTPONED]`

- [ ] Select top-ranked formulations, generate lab protocols, define GC-MS validation criteria.

---

## ✅ Archived: Completed Work (Post-Audit 2026-03-11)

<details>
<summary><b>R.1: Cantera Condensed-Phase Model ✅</b></summary>

- Evaluated options, implemented ideal-condensed in `cantera_export.py`, verified correctness for liquid food matrices.
</details>

<details>
<summary><b>R.2: MFT Formation Template ✅</b></summary>

- Added exact template for 2-Methyl-3-furanthiol generation from Ribose/1,4-dideoxyosone.
</details>

<details>
<summary><b>R.3: Arrhenius A-Factor Transparency ✅</b></summary>

- Populated missing A-factors, added `source_quality` tracking to prevent silent fallbacks to `NaN`.
</details>

<details>
<summary><b>R.5: Henry Constants Expansion ✅</b></summary>

- Expanded Kaw database to accurately model the vaporization of meaty and roasted target volatiles.
</details>

<details>
<summary><b>R.6: Sensory Synergy Model ✅</b></summary>

- Added cross-enhancement model for synergistic interactions (e.g., FFT + methional = boosted meaty profile).
</details>

<details>
<summary><b>R.7: Sequential Bottleneck Fix ✅</b></summary>

- Replaced simplistic `max(barriers)` with correct cumulative pathway bottleneck metrics.
</details>

<details>
<summary><b>R.9: Environment Unification ✅</b></summary>

- Consolidated the confusing multi-environment setup into a well-documented `conda`/`pip` architecture in `environment.yml`.
</details>

<details>
<summary><b>R.10: Repository Hygiene ✅</b></summary>

- Cleaned up debug files, formalized `.gitignore`, organized agent lesson logs.
</details>

<details>
<summary><b>R.11: Uncertainty Quantification ✅</b></summary>

- Implemented standard deviations on predictions and propagated uncertainty through the Boltzmann scoring.
</details>

<details>
<summary><b>R.12: Generalized Deamination Logic ✅</b></summary>

- Generalized internal and hydrolytic deamination pathways for amino ketones, strictly enforcing atom conservation.
</details>

## ✅ Archived: Completed Work (Phases 0–Q)

<details>
<summary><b>Phase A: Literature Arrhenius Calibration ✅</b></summary>

- [x] A.1–A.6: Extracted (A, Ea) values, populated `arrhenius_params.yml`, updated `barrier_constants.py` with `get_arrhenius_params()`, updated Cantera exporter, verified integration, added tests.
</details>

<details>
<summary><b>Phase B: Public ML Potential Integration (MACE-OFF24) ✅</b></summary>

- [x] B.1–B.5: Installed MACE-OFF24, created `mlp_barrier.py`, benchmarked against xTB, integrated into `results_db.py`, added unit tests.
</details>

<details>
<summary><b>Phase C: Full Sensory Prediction Model ✅</b></summary>

- [x] C.1–C.5: Unified `SensoryDatabase`, Stevens' power-law OAV, radar chart output, integration with `InverseDesigner`, unit tests.
</details>

<details>
<summary><b>Phase D: Headspace & Volatility Partitioning ✅</b></summary>

- [x] D.1–D.5: `henry_constants.yml`, `HeadspaceModel`, Cantera→Headspace→Sensory chain, matrix composition in `ReactionConditions`, tests.
</details>

<details>
<summary><b>Phase E: Sigmoid pH Model ✅</b></summary>

- [x] E.1–E.4: Smooth sigmoid pH multipliers, all reaction families covered, tests, Cantera propagation verified.
</details>

<details>
<summary><b>Phase F: AGE/Safety Scoring ✅</b></summary>

- [x] F.1–F.4: Safety scoring, `FormulationResult` fields, Pareto ranking, tests.
</details>

<details>
<summary><b>Phase G: Concentration-Aware FAST Ranking ✅</b></summary>

- [x] G.1–G.4: Boltzmann×concentration weighting, bimolecular correction, `weighted_flux` tracking, tests.
</details>

<details>
<summary><b>Phase H: Bayesian Formulation Optimization ✅</b></summary>

- [x] H.1–H.5: `FormulationOptimizer` with Optuna, EI acquisition, CLI script, trajectory output, tests.
</details>

<details>
<summary><b>Phase I: Matrix-Effect Corrections ✅</b></summary>

- [x] I.1–I.5: Protein binding model, lipid partitioning, `protein_type` in conditions, matrix-corrected headspace chain, tests.
</details>

<details>
<summary><b>Phases J–Q: Architecture Fixes (P1–P8) ✅</b></summary>

- [x] J.1: Cantera phase model documentation
- [x] K.1: Barrier constants exact-match dictionary
- [x] L.1: Gaussian water activity model
- [x] M.1: Per-compound protein binding
- [x] N.1: Per-amino-acid BO concentration variables
- [x] O.1: Family-dependent Kirkwood-Onsager correction
- [x] P.1: Sensory tag expansion
- [x] Q.1: Temporal FAST mode
</details>

<details>
<summary><b>Phase 12: Advanced Kinetic Features ✅</b></summary>

- [x] 12.1–12.5: pH/temperature barrier scaling, multi-reactant Cantera, Joback thermo-gating, solvent scaling, tests.
</details>

<details>
<summary><b>Phase 16: FFT Pathway Bottleneck Resolution ✅</b></summary>

- [x] Root cause identification, H₂S-mediated mechanism, numerical stability, pH-kinetics propagation, verification.
</details>

<details>
<summary><b>Phases 0–11: Foundation, Screening & Refinement ✅</b></summary>

- [x] Base mechanism scaffolding, xTB screening, generative logic, SMIRKS formalization.
- [x] Expanded formulation grid, inverse design, lysine budget & DHA tracking.
- [x] CREST/QCG solvation wrapper, MACE geometry opt, Sella TS search, React-TS diffusion.
- [x] NASA polynomial thermo, barrier source unification, results DB, temperature ramps, regression gate.
</details>
