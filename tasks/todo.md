# Maillard Reaction Computational Framework — Working Plan
## Status: Scientific Accuracy & SOTA Hardening

---

## ⚡ ACTIVE ROADMAP (Post-Audit 2026-03-11)

> [!IMPORTANT]
> **Architectural Review (2026-03-11):** A comprehensive code review identified 15 findings across the pipeline.
> The priorities below focus on **scientific correctness first**, then predictive power, then engineering quality.
> Full analysis: see `implementation_plan.md` in the brain artifacts.

> [!IMPORTANT]
> **Environment Note:** Use `.venv` for full test coverage. The `conda_env` is missing `mace-torch` and will skip MLP/MACE related tests (`tests/qm/test_mlp_*`).

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **🔴 P0** | **R.1: Cantera Condensed-Phase Fix** | `[ ]` | Fix ideal-gas → forward-only or condensed-phase |
| **🔴 P1** | **R.2: MFT Formation Template** | `[x]` | Add the #1 target compound's pathway |
| **🔴 P2** | **R.3: Arrhenius A-Factor Fallbacks** | `[x]` | Fix silent fallback on NaN pre-exponentials |
| **🔴 P3** | **R.4: NASA7 Polynomial Format** | `[ ]` | Fix single-range → two-range format |
| **🟠 P4** | **R.5: Henry Constants Expansion** | `[ ]` | Cover all target/toxic/off-flavour compounds |
| **🟠 P5** | **R.6: Sensory Synergy Model** | `[ ]` | Implement Hofmann synergy pairs |
| **🟠 P6** | **R.7: Sequential Bottleneck Fix** | `[x]` | Replace `max(barriers)` with `1/Σ(1/kᵢ)` |
| **🟠 P7** | **R.8: BO Thread Safety** | `[ ]` | Stop mutating shared InverseDesigner state |
| **🟡 P8** | **R.9: Environment Unification** | `[ ]` | Merge .venv + conda into single env |
| **🟡 P9** | **R.10: Repository Hygiene** | `[x]` | Clean debug files, add lessons.md |
| **🔵 P10** | **R.11: Uncertainty Quantification** | `[x]` | Add error bars to predictions |
| **🔴 P1** | **R.12: Generalized Deamination** | `[/]` | Fix "N-Skip" hack with proper deamination rules |

---

### [ ] R.1: Cantera Condensed-Phase Model `[🔴 P0 | Correctness]`

> **Why:** The framework uses `ideal-gas` thermodynamics for an aqueous food matrix. Reverse rate constants from K_eq include gas-phase entropy terms that don't apply to condensed-phase reactions. Reactions with mole-count changes (retro-aldol, dehydration, Strecker) will have systematically biased equilibria. Ideal-gas densities are ~3 orders of magnitude lower than liquid water.

- [ ] **R.1.1** Evaluate options: (a) `IdealSolidSolution` phase, (b) disable reverse reactions (forward-only Arrhenius), (c) custom Cantera phase class
- [ ] **R.1.2** Implement chosen fix in `cantera_export.py`
- [ ] **R.1.3** Update NASA7 thermo entries to be consistent with the new phase model
- [ ] **R.1.4** Verify: run `ribose + cysteine + leucine` benchmark and check concentration profiles are physically reasonable (no negative concentrations, monotonic intermediate profiles)
- [ ] **R.1.5** Update `tests/integration/test_cantera_sim.py`

---

### [ ] R.2: MFT Formation Template `[🔴 P1 | #1 Target Compound]`

> **Why:** 2-Methyl-3-furanthiol (MFT) is listed as the single most important meaty odorant (threshold 0.0001 µg/kg), yet no template in `smirks_engine.py` generates it. MFT requires the 1,4-dideoxyosone → norfuraneol → MFT pathway, which is distinct from the FFT (furfural + H₂S) pathway already implemented.

- [x] **R.2.1** Add `_mft_formation` template: ribose → 1,4-deoxyosone (via retro-aldol variant) → 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol) → MFT (via H₂S addition)
- [ ] **R.2.2** Add methanethiol formation template (cysteine degradation side-product)
- [ ] **R.2.3** Add DMDS/DMTS formation via methanethiol radical dimerization
- [x] **R.2.4** Verify: run `ribose + cysteine` system, confirm MFT appears in output
- [x] **R.2.5** Add `tests/unit/test_mft_pathway.py`

---

### [ ] R.3: Arrhenius A-Factor Transparency `[🔴 P2 | Credibility]`

> **Why:** 6/11 reaction families in `arrhenius_params.yml` have `A_value: .nan`. The code silently falls back to heuristic barriers, but `todo.md` marks Phase A as ✅ complete—creating a false sense of literature calibration for users.

- [x] **R.3.1** Populate missing A-factors from literature (Martins & van Boekel 2003, Nursten 2005, Hodge 1953 cover most families)
- [x] **R.3.2** Where no literature A exists, use generic TST: `A = kT/h ≈ 6×10¹² s⁻¹` at 150°C as a physically-motivated default
- [x] **R.3.3** Add a `source_quality` field to each Arrhenius entry: `"literature"`, `"estimated"`, or `"TST_default"`
- [ ] **R.3.4** Surface the barrier source in pipeline output so users know which reactions use calibrated vs. estimated kinetics
- [x] **R.3.5** Update `tests/unit/test_arrhenius_params.py` to assert no NaN A-values remain

---

### [ ] R.4: NASA7 Polynomial Format Fix `[🔴 P3 | Correctness]`

> **Why:** `get_nasa_coefficients` in `thermo.py` produces 7 coefficients with a single temperature range (300–1000 K). Cantera's NASA7 format requires 14 coefficients (two ranges). This may cause silent thermo errors or Cantera parse failures.

- [ ] **R.4.1** Modify `get_nasa_coefficients` to produce two-range NASA7 (low: 300–1000 K, high: 1000–3000 K) with 14 coefficients
- [ ] **R.4.2** Alternative: switch Joback-estimated species to `ShomatePoly` which is naturally single-range
- [ ] **R.4.3** Update `cantera_export.py` to emit correct `temperature-ranges: [300.0, 1000.0, 3000.0]` with two data arrays
- [ ] **R.4.4** Add test: load exported YAML in Cantera, verify thermo evaluations at 300 K, 423 K, and 500 K match Joback predictions within 5%

---

### [ ] R.5: Henry Constants Expansion `[🟠 P4 | Headspace Accuracy]`

> **Why:** Only 10 compounds have Henry's Law constants. Compounds not in the database fall back to `Kaw = 0.01`, which is orders of magnitude wrong for most volatiles.

- [ ] **R.5.1** Add Kaw data for all compounds in `desirable_targets.yml` (16 compounds, ~6 missing)
- [ ] **R.5.2** Add Kaw data for all compounds in `off_flavour_targets.yml` and `toxic_markers.yml`
- [ ] **R.5.3** Add Kaw for key intermediates (acetaldehyde, pyruvaldehyde, aminoacetone)
- [ ] **R.5.4** Sources: NIST Chemistry WebBook, Buttery 1969, Roberts & Acree 1995
- [ ] **R.5.5** Update `tests/unit/test_headspace.py` to verify all target compounds have non-default Kaw

---

### [ ] R.6: Sensory Synergy Model `[🟠 P5 | Sensory Accuracy]`

> **Why:** The current `SensoryPredictor` is purely additive. Aroma perception is non-linear: FFT + methional together create "meaty" perception that neither produces alone. Without synergy, the tool cannot distinguish "generic roasted" from "specifically meaty."

- [ ] **R.6.1** Add synergy pair database to `sensory_tags.yml` (Hofmann & Schieberle 2000 provide major pairs)
- [ ] **R.6.2** Implement cross-enhancement model: `OAV_synergy = OAV_base × (1 + Σ boost_factor × OAV_partner)`
- [ ] **R.6.3** Key pairs: FFT + methional (1.5× meaty), MFT + FFT (2× meaty), pyrazine + aldehyde (1.3× roasted)
- [ ] **R.6.4** Add suppression interactions: hexanal suppresses sulfur perception at high concentrations
- [ ] **R.6.5** Tests: verify synergy pairs produce higher "meaty" score than sum of individual contributions

---

### [ ] R.7: Sequential Bottleneck Fix `[🟠 P6 | Ranking Accuracy]`

> **Why:** `recommend.py` uses `path_span = max(barriers)` (energetic span model) which is correct for catalytic cycles but not for linear Maillard cascades. Two consecutive 20 kcal/mol barriers are slower than one 20 kcal/mol barrier, but the current model treats them identically.

- [x] **R.7.1** Replace `max(barriers)` with cumulative bottleneck: `k_eff = 1 / Σ(1/kᵢ)` in `predict_from_steps`
- [x] **R.7.2** Keep energetic span as a secondary metric for comparison/validation
- [x] **R.7.3** Add test: pathway [10, 25, 15] should rank lower than pathway [25] alone

---

### [ ] R.8: Bayesian Optimizer Thread Safety `[🟠 P7 | Robustness]`

> **Why:** `bayesian_optimizer.py` L74 mutates `self.designer.grid` on every trial. This is not thread-safe, corrupts state on crashes, and loses the original grid after the first trial.

- [ ] **R.8.1** Create a fresh `InverseDesigner` per trial, or pass formulation as parameter
- [ ] **R.8.2** Add a `evaluate_single(formulation, conditions)` method to `InverseDesigner` that doesn't require grid mutation
- [ ] **R.8.3** Add test: verify optimizer produces consistent results across repeated runs

---

### [ ] R.9: Environment Unification `[🟡 P8 | Usability]`

> **Why:** Two separate environments (`.venv` and `conda_env`) create confusion. `conda_env` is missing `mace-torch`.

- [ ] **R.9.1** Decide: single conda env with pip packages (recommended for RDKit/Cantera compatibility)
- [ ] **R.9.2** Update `environment.yml` to include all pip dependencies (`mace-torch`, `optuna`, etc.)
- [ ] **R.9.3** Update README installation instructions
- [ ] **R.9.4** Remove `.venv` references or document the migration path

---

### [ ] R.10: Repository Hygiene `[🟡 P9 | Quality]`

- [x] **R.10.1** Add to `.gitignore`: `debug_*`, `tmp_*`, `test_compete.yaml`, `test_failures.log`, `xtb*`, `wbo`, `charges`, `minimal.yaml`
- [x] **R.10.2** Remove committed debug files from repo
- [x] **R.10.3** Create `tasks/lessons.md` (required by `agents.md` workflow)
- [x] **R.10.4** Fix duplicate `if __name__ == "__main__"` block in `kinetics.py`

---

### [ ] R.11: Uncertainty Quantification `[🔵 P10 | SOTA]`

> **Why:** Point estimates with no error bars are risky for guiding expensive wet-lab experiments (£1,500/sample). Scientists need to know "is this prediction reliable enough to skip the lab test?"

- [x] **R.11.1** Add `uncertainty_kcal` field to barrier lookups (heuristic: ±5, xTB: ±3, MLP: ±1, DFT: ±0.5)
- [x] **R.11.2** Propagate uncertainty through the Boltzmann scoring to produce OAV ranges
- [x] **R.11.3** Display confidence intervals in the sensory radar output
- [x] **R.11.4** Add confidence-weighted ranking: penalize formulations with high-uncertainty critical pathways

---

### [ ] R.12: Generalized Deamination Logic `[🔴 P1 | Correctness]`

> **Why (The Nitrogen Issue):** In Phase R.2, we implemented an **"N-Skip" guardrail** in the MFT pathway. 
> 1. **Root Cause:** The MFT template assumes a neutral sugar stoichiometry ($C_5H_8O_4 + 2H_2S \rightarrow MFT + 3H_2O + S$). 
> 2. **Failure:** If applied to an amino-deoxyosone (e.g., from a Ribose-Glycine Amadori product), the reaction becomes unbalanced because the Nitrogen has nowhere to go in the current products list, triggering an `AssertionError`.
> 3. **Gap:** This prevents the framework from predicting critical sulfur volatiles starting from the *actual* protein-bound or amino-acid intermediates unless they are manually deaminated first.

- [ ] **R.12.1** Implement `_deamination_template`: `Amino-Deoxyosone → Neutral Deoxyosone + Free Amine`.
- [ ] **R.12.2** Ensure deamination is triggered before cyclic sulfur volatile templates in `SmirksEngine.enumerate()`.
- [ ] **R.12.3** Remove the `if "N" in s.smiles: continue` guardrail from `_mft_pathway`.
- [ ] **R.12.4** Verify: run `ribose + cysteine` and confirm the `Amadori → Deamination → MFT` chain completes with valid mass balance.
- [ ] **R.12.5** Update `tests/unit/test_smirks_engine.py` with balanced deamination cases.

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

## ✅ Archived: Completed Work (Phases 0–Q)

<details>
<summary><b>Phase A: Literature Arrhenius Calibration ✅</b></summary>

- [x] A.1–A.6: Extracted (A, Ea) values, populated `arrhenius_params.yml`, updated `barrier_constants.py` with `get_arrhenius_params()`, updated Cantera exporter, verified integration, added tests.
- **Note:** Partial completion—A-factors still NaN for 6/11 families. See R.3 for follow-up.
</details>

<details>
<summary><b>Phase B: Public ML Potential Integration (MACE-OFF24) ✅</b></summary>

- [x] B.1–B.5: Installed MACE-OFF24, created `mlp_barrier.py`, benchmarked against xTB, integrated into `results_db.py`, added unit tests.
</details>

<details>
<summary><b>Phase C: Full Sensory Prediction Model ✅</b></summary>

- [x] C.1–C.5: Unified `SensoryDatabase`, Stevens' power-law OAV, radar chart output, integration with `InverseDesigner`, unit tests.
- **Note:** Synergy model (C.2) partially implemented—additive only. See R.6 for follow-up.
</details>

<details>
<summary><b>Phase D: Headspace & Volatility Partitioning ✅</b></summary>

- [x] D.1–D.5: `henry_constants.yml`, `HeadspaceModel`, Cantera→Headspace→Sensory chain, matrix composition in `ReactionConditions`, tests.
- **Note:** Henry database sparse (10 compounds). See R.5 for follow-up.
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
- **Note:** Thread-safety issue in state mutation. See R.8 for follow-up.
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
