# Maillard Reaction Computational Framework — Working Plan
## Status: SOTA Alignment & Production Scaling

---

## ⚡ ACTIVE EXECUTION & SOTA ROADMAP

The core Tier 0/1/2 pipeline is operational and the zero-DFT laptop pipeline is functional. All bug fixes and core infrastructure tasks are now archived. 

The immediate next priority is generating the required 500+ data points to unblock Machine Learning models, and optimizing MACE for sulfur chemistry to eliminate drift.

> [!IMPORTANT]
> **Environment Note:** Use `.venv` for full test coverage. The `conda_env` is missing `mace-torch` and will skip MLP/MACE related tests (`tests/qm/test_mlp_*`).

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **✅** | **Phase 12: Advanced Kinetic Features** | � | pH, Temperature, and Thermo gating completed & verified. |
| **🟠 2** | **Phase 13: Tier 2 Batch Execution (Dataset Gen)** | 🚧 | Generates 500+ DFT points required for fine-tuning ML potentials. |
| **� 3** | **Phase 14: MACE Fine-Tuning (Sulfur)** | 📋 | Fixes 1.11 Å drift on sulfur species (FFT, MFT) for accurate laptop runs. |
| **� 4** | **Phase 15: Δ-ML Network Scaling** | 📋 | Predicts DFT-quality energy correctly at MACE speed without full DFT. |
| **� 5** | **Phase 16: FFT Pathway Bottleneck** | 🚧 | Refine Ribose+Cys+Leu system investigation. |
| **🔵 6** | **Phase 17: Web Dashboard** | ⏳ | Deferred. GUI for food scientists. |
| **🔵 7** | **Phase 18: Experimental Validation Prep** | ⏳ | Deferred until pipeline is fully validated against ≥ 3 benchmarks. |

---

### [DONE] Phase 12: Advanced Kinetic Features `[✅ COMPLETE]`

> **Why:** Standard TST barriers are static. Real Maillard chemistry is heavily pH-dependent (e.g., Schiff Base is optimal at pH 5.5) and influenced by solvent. These features improve accuracy without requiring expensive DFT calculations.
- [x] **12.1 Barrier Auto-Scaling**: Implement pH and temperature-dependent activation scaling (Kirkwood-Onsager approximation) in `src/kinetics.py`.
- [x] **12.2 Multi-Reactant Integration**: Expand `CanteraExporter.add_reaction` to handle higher-order > 2 reactant/product rules.
- [x] **12.3 Dynamic Thermo-Gating**: Use `JobackEstimator` to pre-prune highly endergonic (ΔG > 30 kcal/mol) reactions before ODE integration.
- [x] **12.4 Solvent-Dependent Scaling**: Add empirical barrier scaling based on solvent dielectric constants.
- [x] **12.5 Verification**: Added `tests/unit/test_advanced_kinetics.py` to ensure long-term correctness of pH and solvent scaling.

---

### [ACTIVE] Phase 13: Tier 2 Batch Execution (Dataset Gen) `[🔴 CRITICAL | Diff: 7/10]`

> **Why:** The Machine Learning models (MACE fine-tuning and Δ-ML) are blocked on requiring 500–1000 high-fidelity DFT data points. We must batch-execute the `DFTRefiner` across our generated pathways to build this dataset.
- [ ] **13.1** Ensure all input geometries for the top 500 reactions exist in `data/geometries/xtb_inputs/`.
- [ ] **13.2** Run `src/dft_refiner.py` in batch mode (via multiprocessing or HPC) to converge these 500+ reactions using the PySCF/Sella bridge.
- [ ] **13.3** Save outputs consistently to `results/dft_tier2/*.json` to be consumed by training scripts.

---

### [ACTIVE] Phase 14: MACE Fine-Tuning for Sulfur Chemistry `[🔴 CRITICAL | Diff: 8/10]`

> **Why:** Pre-trained MACE has confirmed 1.11 Å RMSD drift on Maillard sulfur species (FFT, MFT, thiazoles). It cannot be used for production S-chemistry without fine-tuning on our explicit data.
> **Blocked on:** Phase 13 (Tier 2 Batch Execution).
- [ ] **14.1** Extend `scripts/generate_mace_training_data.py` to read `results/dft_tier2/*.json` and convert to `extxyz` format. Add sulfur species augmentation (FFT, MFT, thiazoles) to ensure domain coverage.
- [ ] **14.2** Add `fine_tune()` entrypoint to `src/mlp_optimizer.py` wrapping `mace_run_train` with Maillard-specific hyperparameters (cutoff 5.5 Å, 256 hidden channels).
- [ ] **14.3** Fine-tune MACE on Phase 13 data.
- [ ] **14.4** Re-run benchmark; verify RMSD drift falls below 0.05 Å on sulfur species and barrier MAE < 2 kcal/mol vs DFT.

---

### [ACTIVE] Phase 15: Δ-ML Network Scaling `[🟠 HIGH | Diff: 8/10]`

> **Why:** Train a correction model `E_DFT ≈ E_MACE + Δ_ML` on Phase 13 outputs to predict DFT-quality barriers at MACE speed across the entire network. E_MACE will be fast and structurally accurate (thanks to Phase 14), and Δ_ML fixes the remaining energy residual.
> **Blocked on:** Phase 13 and Phase 14.
- [ ] **15.1** Create `src/delta_ml.py` implementing a KRR or small NN regression on `(E_MACE, E_DFT)` pairs. Target: `E_DFT − E_MACE` residual.
- [ ] **15.2** Training script consuming `results/dft_tier2/*.json` as ground truth. Cross-validation with 80/20 train/test split.
- [ ] **15.3** Add `delta_ml` to the method priority list in `src/results_db.py` `find_barrier()`, positioned between `r2SCAN-3c` and `xtb`.
- [ ] **15.4** Write `tests/test_delta_ml.py` to verify MAE on held-out set is lower than raw MACE MAE.
- [ ] **15.5** Integrate into `DFTRefiner` as a new method tier (invoked after MACE guess, before full DFT SP).

---

### [ACTIVE] Phase 16: FFT Pathway Bottleneck Investigation `[🟡 MEDIUM | Diff: 5/10]`

> **Why:** The literature system (Ribose+Cys+Leu) does not optimally generate the flavor target FFT. Need to investigate the flux bottleneck in the current network constraint.
- [ ] **16.1** Analyze Cantera simulation logs for Ribose+Cys+Leu to identify which intermediate is stockpiling.
- [ ] **16.2** Adjust heuristic barrier or SMIRKS rule to alleviate bottleneck if it's a model artifact, or confirm if chemically accurate.
- [ ] **16.3** Run `test_regression.py` to verify rank improvement of FFT in the yield profile.

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

## 🏛️ ARCHIVE / HISTORICAL PROGRESS (Completed Phases)

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
