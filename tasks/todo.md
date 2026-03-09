# Maillard Reaction Computational Framework — Working Plan
## Status: SOTA Alignment & Production Scaling

---

## ⚡ ACTIVE EXECUTION & SOTA ROADMAP

The core Tier 0/1/2 pipeline is operational. The next objective is to **maximise usefulness on a laptop** by minimising DFT dependency, closing architectural gaps, and making the tool production-ready for food scientists.

> [!IMPORTANT]
> **Architecture Decision (2026-03-08):** DFT barriers are too expensive to compute on a laptop. The revised strategy uses the 17 literature-calibrated heuristic barriers in `barrier_constants.py` as the **primary production tier**. DFT results automatically override heuristics when available, but the tool works without them.

| Priority | Phase | Status | Impact |
|----------|-------|:------:|--------|
| **🔴 1** | **20 — Heuristic DB Population** | 📋 | Immediately enables full-network Cantera simulations (Gap 1, 5) |
| **🔴 2** | **21 — SmirksEngine → Cantera Bridge** | 📋 | Zero-DFT laptop pipeline; full mechanism from precursors (Gap 5) |
| **🔴 3** | **23 — Barrier Source Unification** | 📋 | InverseDesign & Cantera use different barrier sources (Gap 3) |
| **🟠 4** | **24 — NASA Polynomial Thermodynamics** | 📋 | constant-Cp makes reverse rates unphysical (Gap 4) |
| **🟠 5** | **18 — Regression Gate** | 📋 | Protect the heuristic baseline |
| **🟠 6** | **19 — Web Dashboard** | 📋 | Lower adoption barrier for food scientists |
| **� 7** | **3.3 — DFT Refinement Runs** | ⏳ | Optional cloud/HPC; auto-upgrades DB when available |
| **� 8** | **12 — Cantera Microkinetics** | ✅ | Predict concentration-vs-time (GC-MS comparable) |
| **🟢 9** | **15 — Temperature Ramp Modeling** | ✅ | Realistic extrusion profiles, not isothermal |
| **🟢 10** | **16 — Structured Results Database** | ✅ | Scale from 8 to 1000+ reactions |
| **🟢 11** | **17 — GC-MS Comparison Output** | ✅ | Direct overlay with experimental data |
| **🔧 8** | **9 — Explicit Solvation (CREST/QCG)** | 🔧 | Scaffolded; requires CREST binary (not bundled) |
| **🔧 9** | **10 — MLP Geometry Opt (MACE)** | 🔧 | Scaffolded; requires MACE weights (~500MB download) |
| **🔧 10** | **11 — Sella TS Search** | 🔧 | Scaffolded; requires `pip install sella` + ASE calculator |
| **�🟡 11** | **13 — Δ-ML Network Scaling** | 📋 | Blocked on 500+ DFT data points |
| **🟡 12** | **14 — React-TS Diffusion TS Guessing** | 📋 | Frontier: generative TS from 2D graphs |

### 🐛 ACTIVE BUGS — 14 Failing Tests (2026-03-09 Triage)

> **Summary:** 219 passing, 43 skipped, **14 failed** across 6 test files.
> Failures fall into **5 root-cause buckets**. Two are real bugs; three are test–code mismatches from Phase 20→21 architecture changes.

---

#### Bug A: Cantera Test Fixtures Use Unbalanced Dummy Reactions `[7 tests | 🔴 FIX TESTS]`

**Root Cause:** Phase 20 introduced strict mass-balance enforcement in `cantera_export.py` (line 77). Several test fixtures were written *before* this enforcement, using toy reactions like `C → CC` (methane → ethane) or `Ribose → FFT` that are not atom-balanced. The enforcement is correct — the *tests* need updating with balanced dummy reactions or a test-only bypass.

**Affected Tests:**
- `test_cantera_integration.py::TestIsothermalSimulation::test_simple_isothermal_simulation`
- `test_cantera_integration.py::TestIsothermalSimulation::test_isothermal_equilibration`
- `test_cantera_integration.py::TestTimeTemperatureProfile::test_time_temperature_profile`
- `test_cantera_integration.py::TestMultiplePathways::test_multiple_pathways_coexist`
- `test_cantera_integration.py::TestMultiplePathways::test_pathway_selectivity_chemistry`
- `test_cantera_integration.py::TestBarrierSensitivity::test_barrier_effect_on_kinetics`
- `test_temp_ramp.py::test_isothermal_vs_ramp`

**Fix:** Rewrite test fixtures to use atom-balanced dummy reactions that satisfy `add_reaction()`'s mass-balance check. For kinetics tests that only care about rate behaviour (not chemistry), use simple balanced isomerisation reactions like `A ⇌ B` where A and B have identical molecular formula (e.g. `CC=O ⇌ C=CO` — acetaldehyde ⇌ vinyl alcohol).

- [ ] **A.1** Fix `test_cantera_integration.py` — replace all unbalanced `[\"C\"] → [\"CC\"]` with balanced dummy reactions
- [ ] **A.2** Fix `test_temp_ramp.py` — replace `[\"C\"] → [\"CC\"]` with balanced isomerisation
- [ ] **A.3** Fix multi-step test assertions (`test_simple_isothermal_simulation`, `test_multiple_pathways_coexist`) — Ribose→FFT shortcuts are not atom-balanced; redesign with realistic multi-step intermediates or pre-build balanced mechanism

---

#### Bug B: Recommender `predict()` Fails to Match FFT Pathway `[2 tests | � FIX CODE]`

**Root Cause:** `_get_pathway_requirements()` in `recommend.py` computes the set of exogenous reactants needed for `C_S_Maillard_FFT`. Step 1 requires `[CYSTEINE, WATER]`, Step 3 requires `[FURFURAL, H2S, hydrogen]`. Furfural and H2S are produced in earlier steps, but **`water` and `hydrogen` (`[HH]`)** are treated as exogenous requirements. When the test calls `predict(["D-ribose", "L-cysteine"])`, the pathway never activates because `water` and `hydrogen` are missing from the pool.

**Affected Tests:**
- `test_recommend.py::test_recommender_canonical_systems`
- `test_recommend.py::test_recommender_penalties`

**Fix:** Treat ubiquitous small molecules (`water`, `hydrogen`, `ammonia`) as always-available in the reaction environment. Two options:
1. **Preferred:** Add implicit species set in `predict()` — `available_species |= {"water", "hydrogen", "ammonia"}` before the pathway loop.
2. **Alternative:** Remove `WATER` and `_s("hydrogen", "[HH]")` from reactant lists in `curated_pathways.py` and adjust mass balance accordingly (chemically less accurate).

- [ ] **B.1** Add implicit "solvent/ambient" species to `predict()` in `recommend.py`
- [ ] **B.2** Verify both `test_recommender_canonical_systems` and `test_recommender_penalties` pass

---

#### Bug C: Lysine Budget DHA — Scoping Bug `[1 test | 🟡 FIX CODE]`

**Root Cause:** In `recommend.py` `predict_from_steps()`, the DHA lysine budget calculation (lines 292–321) uses `max_r_dist` and `reachable` variables without re-initialising them for each step. These variables leak from the lipid trapping loop above, causing the DHA reachability check to use stale values. The test expects `budget_b > 0.0` but gets `0.0`.

**Affected Test:**
- `test_lysine_budget.py::test_lysine_budget_competition`

**Fix:** Re-initialise `max_r_dist = 0.0` and `reachable = True` inside the DHA loop for each step iteration (line ~298).

- [ ] **C.1** Fix variable scoping in `predict_from_steps()` lysine budget section

---

#### Bug D: MLP/MACE NameError `[3 tests | 🟠 FIX SKIPGUARD]`

**Root Cause:** `mlp_optimizer.py` imports `mace_mp` inside `try/except ImportError`, but when `mace-torch` is partially installed (importable module but broken internals) the `except` catches silently. The test's `@pytest.mark.skipif(MLPOptimizer is None, ...)` guard evaluates `True` because the class *is* importable, but when the constructor calls `mace_mp()` it raises `NameError: name 'mace_mp' is not defined` because the inner import failed silently.

**Affected Tests:**
- `test_mlp_optimizer.py::test_mlp_optimizer_smoke`
- `test_mlp_optimizer.py::test_mlp_optimizer_ts_fallback`
- `test_mlp_optimizer.py::test_dft_refiner_mace_backend`

**Fix:** Make `MLPOptimizer.__init__` more defensive: if `mace_mp` isn't available, raise `ImportError` explicitly. Or change the skip guard to check for `mace_mp` directly:
```python
_MACE_AVAILABLE = False
try:
    from mace.calculators import mace_mp
    _MACE_AVAILABLE = True
except ImportError:
    pass

@pytest.mark.skipif(not _MACE_AVAILABLE, reason="mace-torch not installed")
```

- [ ] **D.1** Fix test skip guard in `test_mlp_optimizer.py` to check `mace_mp` availability directly
- [ ] **D.2** Fix `mlp_optimizer.py` to set a module-level `_MACE_AVAILABLE` flag for reliable detection

---

#### Bug E: PySCF `Mole.to_ase()` Doesn't Exist `[1 test | � FIX CODE]`

**Root Cause:** `dft_refiner.py` line 236 calls `mol.to_ase()` and `mf.as_ase()`, but PySCF's `gto.Mole` has no `to_ase()` method and `mf.as_ase()` doesn't exist either. This is the Sella TS search integration path. The code was scaffolded but never tested end-to-end (noted in Phase 11.7: "PySCF↔ASE bridge is untested").

**Affected Test:**
- `test_explicit_solvation_integration.py::test_dft_refiner_explicit_solvation`

**Fix:** Two options:
1. **Implement the bridge** — convert PySCF Mole to ASE Atoms manually (reading coordinates + elements from `mol.atom_coords()` and `mol.atom_symbol()`).
2. **Skip the Sella path** when Sella is not available and fall through to geomeTRIC, which is the existing fallback. The test should still pass because it uses `is_ts=True` but doesn't require Sella. The real issue is that `self.ts_optimizer` is not `None` even though Sella may not work with PySCF directly.

- [ ] **E.1** Fix the Sella code path in `dft_refiner.py` to handle the PySCF↔ASE conversion properly, or make `ts_optimizer` initialization more defensive

---

### [DEPRECATED] Phase 20: Heuristic DB Population ~~`[🔴 CRITICAL]`~~ → Absorbed into Phase 21

> **ADR (Architectural Decision Record):** This phase was initially implemented and then deliberately reverted. The original plan was to populate `ResultsDB` with heuristic SMILES from `curated_pathways.py` using `literature_heuristic` barriers.
>
> **Why this was wrong:**
> 1. `ResultsDB` operates on *exact SMILES* strings. Populating it with 16 hardcoded molecule pairs misuses the DB as a generic lookup table — it will not match if a user substitutes even one precursor.
> 2. `barrier_constants.py` already provides the correct family-level constants. There is no need to copy them into the DB.
> 3. Pre-populating a DB to serve as a generic fallback is fundamentally inelegant and creates a non-generalisable system.
>
> **Correct approach (implemented in Phase 21):** The simulation script queries `ResultsDB` first for exact computed (DFT/xTB) barriers on specific molecules, and live-falls back to `barrier_constants.get_barrier(family)` otherwise. The DB grows organically from real computations only.

---

### [ACTIVE] Phase 21: SmirksEngine → Cantera Bridge `[🔴 CRITICAL | Diff: 5/10]` *(Closes Gap 1, 5)*

> **Why:** This is the *correct* zero-DFT laptop pipeline. `SmirksEngine` dynamically generates the exact, mass-balanced reaction network for any precursor combination. The bridge assigns barriers live: first checking `ResultsDB` for any existing DFT/xTB calculation for that exact SMILES pair, then falling back directly to `barrier_constants.get_barrier(reaction_family)`. No DB pre-population step is needed.

- [x] **21.1** Extend `run_cantera_kinetics.py` CLI with `--from-smirks` flag and adjust `run_simulation` entry point.
- [x] **21.2** Implement **SmirksEngine integration**: Initialize engine, enumerate network from `--precursors`, and deduplicate `ElementaryStep` objects.
- [x] **21.3** Implement **Dual-Lookup Barrier logic**:
    - Try `ResultsDB.find_barrier(reactants, products)` for exact calculated match.
    - Fallback to `barrier_constants.get_barrier(reaction_family)` for heuristic baseline.
    - Log which source was used for each reaction (transparency).
- [x] **21.4** Precursor Mapping: Ensure `--precursors` names (ribose, cysteine) map to the exactly balanced SMILES expected by SmirksEngine.
- [x] **21.5** End-to-end test: `python scripts/run_cantera_kinetics.py --from-smirks --precursors ribose:0.1,cysteine:0.1 --temp 150 --time 1800`.
- [x] **21.6** Database Override Verification: Manually add a fake low barrier to DB for a specific step and confirm Cantera simulation shifts its flux accordingly.
- [x] **21.7** Consistency Check: Verify Cantera product concentrations correlate with `inverse_design.py` pathway scores.

---


### [ACTIVE] Phase 23: Barrier Source Unification `[🟠 HIGH | Diff: 3/10]` *(Closes Gap 3)*

> **Why:** `inverse_design.py` uses `barrier_constants.get_barrier()` for scoring. `run_cantera_kinetics.py` queries `ResultsDB`. These are **different barrier values** for the same reaction families. A formulation ranked #1 in Inverse Design might rank #5 in Cantera mode. This erodes trust.

- [x] **23.1** Centralize lookup: Add `ResultsDB.get_best_barrier(reactants, products, family)` as the single source of truth.
- [x] **23.2** Refactor `run_cantera_kinetics.py` to use the centralized `ResultsDB.get_best_barrier()` method (removing duplicated logic).
- [x] **23.3** Refactor `inverse_design.py`: Initialize `ResultsDB` and query it for every step, ensuring fast-mode rankings "see" DFT results.
- [x] **23.4** Consistency Test: Run both Inverse Design and Cantera simulation for Ribose+Cysteine and verify identical barrier usage.
- [x] **23.5** Override Verification: Confirm that a high-level DB entry affects both ranking and kinetics flux.

---

### [ACTIVE] Phase 24: NASA Polynomial Thermodynamics `[🟠 HIGH | Diff: 6/10]` *(Closes Gap 4)*

> **Why:** `cantera_export.py` uses `constant-cp` with `h0=0, s0=0, cp0=75000` for ALL species. This makes equilibrium constants meaningless and reverse reaction rates unphysical. Any system where reversibility matters (Schiff base hydrolysis, Amadori equilibrium) will give wrong results.

- [x] **24.1** Implement `src/thermo.py`: Use Joback group additivity for $H_f$ and $C_p(T)$ estimation.
- [x] **24.2** SMILES-to-Groups: Add RDKit-based substructure matching to decompose Maillard intermediates into Joback groups.
- [x] **24.3** NASA-7 Polynomial Fitter: Fit $C_p(T)$ to 7 coefficients and integrate to obtain $H(T)$ and $S(T)$.
- [x] **24.4** Update `CanteraExporter`: Inject computed NASA polynomials into the YAML mechanism.
- [x] **24.5** Physical Validation: Verify that reverse rates in Cantera now reflect real thermodynamics (equilibrium constants).
- [x] **24.6** Non-Isothermal Test: Run a simulation with a temperature ramp and verify that thermo-dependent fluxes change realistically.

---

### [ACTIVE] Phase 18: Automated Regression Gate `[🟡 MEDIUM | Diff: 4/10]`

> **Why:** As we add MACE, CREST, Δ-ML, each change risks *regressing* accuracy on previously validated systems. An automated gate ensures we never break the Literature Validation.

- [x] **18.1** Define Ground Truth Data: Ensure `data/lit/` contains target JSONs for the 3 canonical systems (Ribose+Cysteine, Glucose+Glycine, Ribose+Cysteine+Leucine).
- [x] **18.2** Create `tests/test_regression.py`: A parametrized `pytest` suite marked with `@pytest.mark.slow`.
- [x] **18.3** End-to-End Test Logic: For each system, run `SmirksEngine` generation $\rightarrow$ Cantera Export $\rightarrow$ `KineticsEngine` simulation.
- [x] **18.4** Yield Assertions: Assert that the dominant experimental flavor compounds (e.g., FFT, Pyrazine) are correctly ranked in the top 3 most abundant volatiles produced by the simulation.
- [x] **18.5** CI Integration: Setup a `test-regression` npm/make script or document how to run the gate before merging major PRs.

---

### [ACTIVE] Phase 19: Web Dashboard for Food Scientists `[🟡 MEDIUM | Diff: 6/10]`

> **Why:** The CLI is powerful but alien to most food scientists. A web interface with visual outputs would dramatically lower the barrier to adoption for the alt-protein community.

- [ ] **19.1** Create `app/` directory with Streamlit app.
- [ ] **19.2** Input form: precursors, pH, temperature, target sensory profile.
- [ ] **19.3** Output: ranked formulation table, volatile concentration chart, sensory radar plot.
- [ ] **19.4** Connect to existing `SmirksEngine` + `InverseDesigner` backend.

---

### [ACTIVE] Phase 3: Tier 2 — DFT Refinement (r2SCAN-3c // wB97M-V)

- [x] **3.1 Workflow Implementation:** `DFTRefiner` in `src/dft_refiner.py`.
- [x] **3.2 Quasi-Harmonic Correction:** `QuasiHarmonicCorrector` in `src/thermo.py`.
- [/] **3.3 Compute barriers for key bifurcations:** `[Diff: 9/10]`
      *(Note: Inputs for all 8 reactions generated and mapped in `data/geometries/xtb_inputs/`).*
    - [x] 3.3a: Amadori rearrangement ✅ (Fast-mode verified: 1.09 kcal/mol)
    - [ ] 3.3b: 2,3-enolisation vs 1,2-enolisation bifurcation point (Inputs Ready)
    - [ ] 3.3c: Strecker decarboxylation (α-dicarbonyl + amino acid) (Inputs Ready)
    - [ ] 3.3d: Cysteine + ribose → FFT (via furfural + H₂S) (Inputs Ready)
    - [ ] 3.3e: Ribose retro-aldol → 1,4-dideoxyosone → MFT (Inputs Ready)
    - [ ] 3.3f: DHA β-elimination (Ser → dehydroalanine) (Inputs Ready)
    - [ ] 3.3g: Off-flavour trapping: hexanal + amino acid → Schiff base (Inputs Ready)
    - [ ] 3.3h: α-aminoketone dimerisation → pyrazine (Inputs Ready)
- [x] **3.4 Validate with IRC:** `[Diff: 6/10]` Native 'Displace + Optimize' engine in `src/dft_refiner.py`.
- [x] **3.5 Verification Study:** `[Diff: 4/10]` Scaffolded `verify_barrier` in `DFTRefiner`. Execution pending results.

---

### [ACTIVE] Phase 9: Explicit Solvation Automation (CREST/QCG) `[🔧 SCAFFOLDED | Diff: 7/10]`

> **Why:** Implicit solvation (ddCOSMO) systematically overestimates proton-transfer barriers by 10–25 kcal/mol. The Amadori rearrangement and enolisation steps are water-catalyzed; without explicit water molecules, our ±1 kcal/mol target is unreachable. SOTA §5 identifies this as the single largest source of error.
>
> **Status Note (2026-03-08):** Code is implemented with graceful degradation (`try/except`). Requires CREST binary to be separately installed and on `$PATH`. Not bundled with the project.

- [x] **9.1** Install Dependency: `conda_env/bin/crest` already present. Add `crest>=3.0` to `environment.yml`.
- [x] **9.2** Create core engine `src/solvation.py`: wrapper around `crest` binary (QCG mode).
    - [x] `generate_solvated_cluster(xyz_string, n_water=3, freeze_core=True)` method.
    - [x] **Elegance Requirement (per `agents.md`):** A naive approach will crash TS geometries into a local minimum. If `freeze_core=True`, generate a `.xcontrol` file that fixes the Cartesian coordinates of the solute atoms so only the water molecules conform.
    - [x] Run `crest input.xyz -qcg water -nsolv {n_water} -cinp .xcontrol`.
    - [x] Parse `crest_best.xyz` and return the cluster coordinates.
- [x] **9.3** Integrate seamlessly into `src/dft_refiner.py`:
    - [x] `use_explicit_solvent: bool` and `n_water: int` parameters already exist in `DFTRefiner`.
    - [x] TS-specific handling hook exists — calls `SolvationEngine.generate_solvated_cluster(freeze_core=True)` when enabled.
    - [x] Merged `Mole` object assembly with explicit water molecules already in place.
- [x] **9.4** Write `tests/test_solvation.py` (Verification Before Done):
    - [x] Assert that a 3-water cluster generates the correct atom count.
    - [x] **Crucial test:** Asserted solute atoms have not moved by more than 0.05 Å when `freeze_core=True` is provided (proving the TS won't be ruined) — against CREST's actual output.
- [x] **9.5** End-to-End Benchmark:
    - [x] Verified `DFTRefiner` with `use_explicit_solvent=True` correctly triggers the CREST/QCG workflow (or heuristic fallback) and passes integration tests.
- [ ] **9.6** ⚠️ **PRODUCTION GAP:** CREST binary is not pip-installable. Add installation instructions to README and verify on a clean machine.

---

### [ACTIVE] Phase 10: MLP-Accelerated Geometry Optimization (MACE) `[🔧 SCAFFOLDED | Diff: 8/10]`

> **Why:** SOTA §2 recommends replacing DFT-level optimization with MACE MLP, reducing optimization from hours to seconds. However, pre-trained MACE will fail on sulfur chemistry and sugar fragmentation without fine-tuning on Maillard-specific DFT data.
>
> **Status Note (2026-03-08):** Code is implemented with graceful degradation. Requires `mace-torch` and `torch` (~2GB). Pre-trained model has verified 1.11Å drift on sulfur species — cannot be used for production Maillard TS work without fine-tuning, which itself requires 500+ DFT data points not yet available.

- [x] **10.1** Add `mace-torch` and `torch` dependencies.
- [x] **10.2** Create `src/mlp_optimizer.py`: ASE-based wrapper around `mace-mp-0`.
- [x] **10.3** Refactor `DFTRefiner.optimize_geometry()` for Algorithmic Decoupling.
- [x] **10.4** Create Fine-Tuning Pipeline (`scripts/generate_mace_training_data.py`).
- [x] **10.5** Verification & Benchmarking (`tests/test_mlp_optimizer.py`).
- [x] **10.6** ⚠️ **VALIDATION GATE:** Pre-trained MACE fails on Maillard sulfur species (1.11 Å drift).
- [ ] **10.7** ⚠️ **PRODUCTION GAP:** Fine-tuning requires 500+ DFT data points. Blocked on Phase 3.3 batch execution.

---

### [ACTIVE] Phase 11: Sella Eigenvector-Following TS Search `[🔧 SCAFFOLDED | Diff: 6/10]`

> **Why:** SOTA §3 recommends Sella for direct saddle-point optimization. Instead of optimizing entire reaction paths (NEB), Sella walks directly toward the nearest saddle point, which is drastically faster when the initial xTB guess is good.
>
> **Status Note (2026-03-08):** Code is implemented with graceful degradation. Requires `pip install sella` and a working ASE calculator. The MACE+Sella path is blocked on MACE fine-tuning (Phase 10.7). The PySCF+Sella path requires the PySCF↔ASE bridge which is only scaffolded.

- [x] **11.1** Add `sella` to `environment.yml` and explicitly install it in the `.venv`.
- [x] **11.2** Create `src/ts_optimizer.py`: wrapper around `sella` Python package.
    - [x] `find_ts(atoms, calculator)` accepting ASE calculator (MACE or PySCF).
- [x] **11.3** Update `DFTRefiner.optimize_geometry()` to use `TSOptimizer` for TS searches when `is_ts=True`.
- [x] **11.4** Implement automatic fallback to `geomeTRIC` if Sella fails to converge.
- [x] **11.5** Write `tests/test_ts_optimizer.py`: verify saddle point found for H-transfer model system.
- [x] **11.6** Update `MLPOptimizer.optimize_ts()` to use Sella + MACE.
- [ ] **11.7** ⚠️ **PRODUCTION GAP:** PySCF↔ASE bridge for Sella is untested end-to-end.

---

### [ACTIVE] Phase 13: Δ-ML Network Scaling `[🟡 MEDIUM | Diff: 8/10]`

> **Why:** SOTA §7 proposes training a correction model `E_DFT ≈ E_MACE + Δ_ML` on 500–1000 DFT-calculated Maillard reactions to predict DFT-quality energies across the entire network without running expensive DFT on every step.

- [ ] **13.1** Create `src/delta_ml.py`: train simple regression (KRR or small NN) on `(MACE_energy, DFT_energy)` pairs.
- [ ] **13.2** Create training script consuming `results/dft_tier2/*.json` as ground truth.
- [ ] **13.3** Write `tests/test_delta_ml.py`: verify correction model reduces MAE vs. raw MACE energies.
- [ ] **13.4** ⚠️ **BLOCKED:** Requires 500+ diverse DFT data points (Phase 3.3 batch execution must complete first).

---

### [ACTIVE] Phase 14: React-TS Diffusion Model (Frontier) `[🟡 MEDIUM | Diff: 9/10]`

> **Why:** SOTA §3 identifies React-TS as the 2026 SOTA for generating 3D saddle points from 2D molecular graphs. Uses SE(3)-equivariant stochastic diffusion for sub-angstrom TS prediction. However, success rates are biased toward pharmaceutical datasets.

- [ ] **14.1** Create `src/diffusion_ts.py`: wrapper around React-TS inference.
    - `generate_ts_guess(reactant_smiles, product_smiles)` method.
- [ ] **14.2** Confidence scoring to trigger xTB fallback for out-of-distribution reactions.
- [ ] **14.3** ⚠️ **VALIDATION GATE:** Extensive validation against DFT-computed TSs required before replacing xTB.

---

## [DEFERRED] Phase 5: Experimental Validation Preparation

> **Context:** Physical wet-lab validation is deferred until the complete *in silico* pipeline is fully operational with SOTA methods.

- [ ] **5.1** Select 5–10 top-ranked novel formulations from computational output.
- [ ] **5.2** Generate a lab protocol summary.
- [ ] **5.3** Define GC-MS validation criteria (expected retention times from NIST).
- [ ] **5.4** Document comparison methodology for predicted vs. observed agreement.

---

## 🏛️ HISTORICAL PROGRESS (Completed Phases)

<details>
<summary><b>Phase 22: Documentation Sync ✅ (Expand to View)</b></summary>

### Phase 22: Documentation Sync `[🔴 CRITICAL | Diff: 2/10]` *(Closes Gap 2, 6)*
> **Why:** `docs/architecture.md` still references ORCA, Gaussian, and "Microkinetics out of scope." Phases 9-11 are marked ✅ in some places but require external binaries. A new user or collaborator would get a fundamentally wrong picture of the tool.
- [ ] **22.1** Rewrite `architecture.md` §3.3 to reference PySCF (not ORCA/Gaussian) and r2SCAN-3c // wB97M-V (not B3LYP-D3).
- [ ] **22.2** Update §4 (Out of Scope): remove Microkinetic ODE modelling (it's implemented). Add "Cloud/HPC DFT" as a deferred item.
- [ ] **22.3** Add a "Dependency Matrix" to README: for each feature, list what's required (pip-installable vs external binary).
- [ ] **22.4** Reconcile the README Roadmap section with this `todo.md` table.

- [x] **22.1** Updated `architecture.md` to reference PySCF and correct composite protocol.
- [x] **22.2** Updated `architecture.md` §4 to move Microkinetics to implemented status.
- [x] **22.3** Added Dependency Matrix to `README.md` (Gap 6).
- [x] **22.4** Reconciled README and roadmap tiers.
</details>

<details>
<summary><b>Phase 17: GC-MS Comparison Output ✅ (Expand to View)</b></summary>

- [x] **17.1 Automatic Mechanism Discovery:** Refactor `scripts/run_cantera_kinetics.py` to build mechanisms dynamically from `ResultsDB` rather than hardcoded SMILES.
- [x] **17.2 Literature Data Import:** Create `data/lit/ribose_glycine_2021.json` with experimental molarities for target volatiles (furfural, pyrazines).
- [x] **17.3 Implement `scripts/compare_sim_to_lit.py`:** Automated script that runs a simulation matching literature conditions and returns error metrics (MAE/R²).
- [x] **17.4 Comparison Visualization:** Generate bar charts comparing Predicted vs Experimental yields for the target aromatics.
- [x] **17.5 Threshold Calibration:** Adjust sensory mapping coefficients if simulation systematically over/under-predicts volatile yield.
</details>

<details>
<summary><b>Phase 16: Structured Results Database ✅ (Expand to View)</b></summary>

- [x] **16.1 Design SQLite Schema:** Reactions table (reactant/product SMILES, family), Barriers table (energy, method, basis, solvation, CPU time, timestamp).
- [x] **16.2 Implement `src/results_db.py`:** `ResultsDB` class focusing on fast lookup by reaction SMILES and provenance filtering.
- [x] **16.3 Integrate `DFTRefiner` logging:** Auto-write to DB after every successful converged calculation.
- [x] **16.4 Migration Script:** Convert existing `refinement_all.json` data into the new DB format.
- [x] **16.5 Query Interface:** Implement basic CLI/API for `get_best_barrier(reactants, products)` with automatic method prioritization (DFT > xTB).
</details>

<details>
<summary><b>Phase 15: Temperature Ramp Modeling ✅ (Expand to View)</b></summary>

- [x] **15.1 Extend `kinetics.py` support for $T(t)$:** Update `simulate_network_cantera` to accept `temperature_profile` list.
- [x] **15.2 Implement Reactor Callback:** Use Cantera's dynamic state updates to interpolate $T$ during integration.
- [x] **15.3 CLI Ramp Support:** Update `scripts/run_cantera_kinetics.py` to accept `--temp-ramp` CSV files.
- [x] **15.4 Verification:** Write `tests/test_temp_ramp.py` asserting divergence between isothermal and ramped results.
- [x] **15.5 Visualization:** Add temperature-vs-time overlaid on concentration plots in CLI output.
</details>

<details>
<summary><b>Phase 12: Cantera Microkinetic Integration ✅ (Expand to View)</b></summary>

- [x] **12.1 Extend `kinetics.py` with `simulate_network_cantera()`:** wraps Cantera Python API for batch reactor simulations.
- [x] **12.2 Create `src/cantera_export.py`:** exports Maillard network as Cantera-compatible `.yaml`.
- [x] **12.3 Setup Environment:** `cantera>=3.0` verified in environment.
- [x] **12.4 Verification Tests:** `tests/test_cantera_integration.py` and `tests/test_cantera_cli.py` (23/23 pass).
- [x] **12.5 Cantera CLI:** `scripts/run_cantera_cli.py` for direct kinetic simulations.
- [x] **12.6 Sensory Prediction:** Integrated with `--predict-sensory` in kinetics results.
</details>

<details>
<summary><b>Phase 8: Scientific Credibility & Formulation Depth ✅ (Expand to View)</b></summary>

- [x] **8.A — Expand Formulation Grid**
- [x] **8.B — Balance SmirksEngine Templates**
- [x] **8.C — Calibrate FAST-mode Heuristic Barriers**
    - [x] Literature Validation Gate passed.
- [x] **8.D — Concentrations + Scoring Function Redesign**
- [x] **8.E — Lipid-Maillard Synergy (Catalytic Role)**
</details>

<details>
<summary><b>Phase 7: Plant-Based Formulation Tools & PBMA Usability ✅ (Expand to View)</b></summary>

- [x] **7.1 Wire SmirksEngine to Recommender:** Created a unified pipeline script.
- [x] **7.2 Implement Core PBMA Precursors & Degradation Rules**
- [x] **7.3 Advanced CLI Formulation Interface**
- [x] **7.4 Sensory & Trapping Output Metrics**
- [x] **7.5 Inverse Design Mode**
- [x] **7.6 Lysine Budget & DHA Competition Metric**
- [x] **7.7 Water Activity CLI Flag**
</details>

<details>
<summary><b>Phase 6: Automated Pathway Generation ✅ (Expand to View)</b></summary>

- [x] **6.1 Deterministic Rule Enumeration (Hybrid SMIRKS + Templates)**
- [x] **6.2 Extended Template Coverage (Remaining 6 Families)**
</details>

<details>
<summary><b>Phase 4: Integration & Precursor Recommendation Prototype ✅ (Expand to View)</b></summary>

- [x] **4.0 Annotate curated pathways with target compounds**
- [x] **4.1 Rewrite `src/recommend.py`**
- [x] **4.2 Add multi-precursor comparison variants**
- [x] **4.3 Add penalty scoring for competing pathways**
- [x] **4.4 Add toxicity flags**
- [x] **4.5 Validate output on 3 canonical model systems**
- [x] **4.6 Write `tests/test_recommend.py`**
- [x] **4.7 Output format**
- [x] **4.8 Document known limitations**
</details>

<details>
<summary><b>Phase 2B: Physics & Implementation Improvements ✅ (Expand to View)</b></summary>

- [x] **2B.1 Replace fake barrier formula in `xtb_screener.py`:** Implemented xTB `--path` (NEB) mode for real approximate ΔE‡ values.
- [x] **2B.2 Add pH / Temperature / Water Activity (aᵥ) parametrisation:**
- [x] **2B.3 Add explicit water molecules in transition states for proton-transfer steps:**
- [x] **2B.4 Fix incorrect D-Ribose SMILES in `benchmark_xtb.py`:**
- [x] **2B.5 Validate Skala API against Microsoft documentation:**
- [x] **2B.6 Uncheck Phase 3 items that were only scaffolded, not executed** (3.3a–3.3h, 3.4, 3.5).
- [x] **2B.7 Multi-conformer sampling in `xtb_screener.py`:**
- [x] **2B.8 Add physics-level integration tests:**
- [x] **2B.9 Improve pathway ranking in `pathway_ranker.py`:**
</details>

<details>
<summary><b>Phase 2: Tier 1 — xTB Screening ✅ (Expand to View)</b></summary>
- [x] **2.1** Set up test framework.
- [x] **2.2** Write `src/pathway_extractor.py`.
- [x] **2.3** Write `tests/test_pathway_extractor.py`.
- [x] **2.4** Write `src/xtb_screener.py`.
- [x] **2.5** Write `tests/test_xtb_screener.py`.
- [x] **2.6** Parallelised via Python `multiprocessing`.
- [x] **2.7** Write `src/pathway_ranker.py`.
- [x] **2.8** Write `tests/test_pathway_ranker.py`.
- [x] **2.9** Verify Amadori barrier.
- [x] **2.10** Verify Ribose > Glucose.
- [x] **2.11** Documented limitations.
- [x] **2.12** Curated explicit pathways.
- [x] **2.13** Write `scripts/run_curated_screening.py`.
</details>

<details>
<summary><b>Phase 1: Knowledge Encoding ✅ (Expand to View)</b></summary>
- [x] **1A. Build the Target Compound Database**
- [x] **1B. Encode Domain-Specific Reaction Families for RMG**
- [x] **1C. RMG Rule Validation & Debugging**
</details>

<details>
<summary><b>Phase 0: Project Infrastructure ✅ (Expand to View)</b></summary>
- [x] **0.1** Initialised repo.
- [x] **0.2** Directories created.
- [x] **0.3** Conda env spec.
- [x] **0.4** Smoke-test.
</details>

---
*Last Updated: 2026-03-08*
