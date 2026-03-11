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
| **🔴 P3** | **R.4: NASA7 Polynomial Format** | `[x]` | Fix single-range → two-range format |
| **🟠 P7** | **R.8: BO Thread Safety** | `[x]` | Stop mutating shared InverseDesigner state |

---

### [x] R.4: NASA7 Polynomial Format Fix `[🔴 P3 | Correctness]`

> **Why:** `get_nasa_coefficients` in `thermo.py` produces 7 coefficients with a single temperature range (300–1000 K). Cantera's NASA7 format requires 14 coefficients (two ranges). This may cause silent thermo errors or Cantera parse failures.

- [x] **R.4.1** Modify `get_nasa_coefficients` to produce two-range NASA7 (low: 300–1000 K, high: 1000–3000 K) with 14 coefficients
- [x] **R.4.2** Alternative: switch Joback-estimated species to `ShomatePoly` which is naturally single-range
- [x] **R.4.3** Update `cantera_export.py` to emit correct `temperature-ranges: [300.0, 1000.0, 3000.0]` with two data arrays
- [x] **R.4.4** Add test: load exported YAML in Cantera, verify thermo evaluations at 300 K, 423 K, and 500 K match Joback predictions within 5%

---

### [x] R.8: Bayesian Optimizer Thread Safety `[🟠 P7 | Robustness]`

> **Why:** `bayesian_optimizer.py` L74 mutates `self.designer.grid` on every trial. This is not thread-safe, corrupts state on crashes, and loses the original grid after the first trial.

- [x] **R.8.1** Create a fresh `InverseDesigner` per trial, or pass formulation as parameter
- [x] **R.8.2** Add a `evaluate_single(formulation, conditions)` method to `InverseDesigner` that doesn't require grid mutation
- [x] **R.8.3** Add test: verify optimizer produces consistent results across repeated runs

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
