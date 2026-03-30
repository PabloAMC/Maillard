# Maillard Strategic Roadmap: Phase 2

## North Star
Evolve the Maillard framework from a static, heuristic-heavy decision surface into a dynamic, time-resolved kinetic engine, while aggressively paying down architectural debt that threatens long-term maintainability.

---

## Priority 0 — Test Velocity & Code Health (Immediate Wins)
**Goal:** Cut the 1h 39m test suite to ~30–40 min; eliminate known code quality hazards.
**Root cause audit (March 2026):** Full pipeline runs (`SmirksEngine.enumerate()` → `Recommender.predict_from_steps()`) are called repeatedly across tests with no cross-test caching. `_step_exists()` in `smirks_engine.py` is O(N) per call → O(N²) overall. Three `xfail` blind-spot tests run the full pipeline before unconditionally skipping. No session fixtures exist to amortize heavy benchmark evaluations.

### 0.1 Fix O(N²) `_step_exists` in `smirks_engine.py` — **[L1 | ✅ Done 2026-03-30]**
- Replace the `_step_exists(step, all_steps)` linear scan in Tier B with the `seen_step_keys: Set[str]` set that is already maintained for Tier A. Makes deduplication O(1) throughout enumerate().
- Estimated saving: 10–20% speedup on all pipeline-heavy tests.

- [x] Propagate `seen_step_keys` initialization to before the Tier B loop.
- [x] Replace all `if not _step_exists(step, all_steps):` calls with O(1) set lookup.
- [x] Remove now-dead `_step_exists` function (zero callers).
- [x] Add module-level `_ENUMERATE_CACHE` dict cache keyed on (frozenset SMILES, max_gen, pH, temp) — eliminates repeated enumeration for identical precursor pools across tests.

### 0.2 Fix wasted `xfail` runs in `test_blind_spots.py` — **[L1 | ✅ Done 2026-03-30]**
- Three tests call `evaluate_single()` (full pipeline, ~1–2 min each) **before** calling `pytest.xfail()`. Moving `xfail` to the first line of each test makes them skip instantly.
- Estimated saving: ~4–6 min per full run.

- [x] Move `pytest.xfail(reason)` to the first line of each blind-spot test body.

### 0.3 Add `@pytest.mark.slow` to full-pipeline scientific tests — **[L1 | ✅ Done 2026-03-30]**
- Tag `test_benchmarks.py` (11 parametrized pipeline runs), `test_time_propagation.py` (6 full runs), `test_accessibility_sensitivity.py` (integration test), `test_safety_and_flux.py`, and `test_bo_interventions.py` with `@pytest.mark.slow`.
- Enables `pytest -m "not slow"` for sub-5-minute dev iteration without deleting any tests.

- [x] Add `@pytest.mark.slow` to all tests that invoke `evaluate_single()` / `evaluate_benchmark()` / `evaluate_all()`.

### 0.4 Session-scoped fixtures for repeated benchmark evaluations — **[L2 | ✅ Done 2026-03-30]**
- The Mottram 1994 benchmark (`cys_ribose_150C_Mottram1994.json`) is evaluated independently in at least 4 test files. A `@pytest.fixture(scope="session")` runs it once and shares the result.
- Same pattern for Farmer 1999, and the pea/soy matrix benchmarks.
- Estimated saving: ~20–30 min on the full run (eliminates ~10–15 redundant pipeline executions).

- [x] Add `mottram_evaluation`, `farmer_evaluation`, `hofmann_meaty_result`, `mottram_meaty_result` session fixtures to `tests/conftest.py`.
- [x] Refactor `test_mottram_coverage.py`, `test_farmer_coverage.py`, `test_lipid_oxidation_guard.py` to consume the fixtures.

### 0.5 Pre-compile SMIRKS rule reactions in `smirks_engine.py` — **[L1 | ✅ Done 2026-03-30]**
- `_apply_smirks_rule` calls `AllChem.ReactionFromSmarts(smirks)` inside the hot loop (once per rule per SMIRKS generation). Pre-compiling at module load eliminates this repeated cost.

- [x] Create `_COMPILED_SMIRKS_RULES: List[Tuple[str, str, object, str]]` at module scope by compiling each rule at import time.
- [x] Pass compiled reactions into `_apply_smirks_rule` instead of the raw SMIRKS string.

### 0.6 Decompose `benchmark_reporting.py` (75KB) — **[L3 | Medium-term]**
- The Priority 2 refactor decomposed `benchmark_validation.py` but `benchmark_reporting.py` absorbed ~75KB of mixed responsibilities: statistics, markdown rendering, snapshot builders, target status builders, promotion artifacts, and closure audit — all in one file.
- Proposed split: `benchmark_summary.py` (stats + evaluation summary), `benchmark_markdown.py` (all `render_*_markdown()` functions), `benchmark_artifacts.py` (snapshot/promotion/closure builders).

- [ ] Extract `render_*_markdown()` functions → `benchmark_markdown.py`.
- [ ] Extract `snapshot_*`, `build_matrix_*`, `build_benchmark_index` → `benchmark_artifacts.py`.
- [ ] Update `benchmark_validation.py` facade to re-export from the new sub-modules.

---

## Priority 1 — Time-Resolved Microkinetics (The Major Scientific Upgrade)
**Goal:** Transition core prediction to a true time-resolved ODE engine.

- [ ] **1.1 Wire ODE solver into `pipeline.py`** — **[L4 | Claude 4.6]**
  - Implement `src/ode_kinetics.py` using `scipy.integrate.solve_ivp` to propagate species concentrations through time given the SMIRKS-enumerated reaction graph and computed barriers.
  - Implement a `mode="kinetic"` flag in the main pipeline that routes to the ODE solver instead of the budget-projection heuristic.
  - Bridge `smirks_engine.py` outputs into the ODE rate matrix dynamically.
- [ ] **1.2 Expose Temporal Dynamics in `recommend.py`** — **[L4 | Claude 4.6]**
  - Update `predict_from_steps` to accept optional pre-computed ODE time-series if provided.
  - Capture peak vs. plateau states for key markers from ODE output.
- [ ] **1.3 Dynamic Temperature & Moisture** — **[L3 | Claude 3.5]**
  - Allow `ReactionConditions` to accept time-temperature arrays for extrusion history.

---

## Priority 2 — Disarm the Monoliths (Refactoring)
**Goal:** Decompose 112K `benchmark_validation.py` and overloaded `recommend.py`.

- [/] **2.1 Decompose `benchmark_validation.py` (111K bytes)** — **[L4 | Claude 4.6]**
  - Move anchors to `benchmark_registry.py`.
  - Move execution to `benchmark_evaluator.py`, reporting to `benchmark_reporting.py`, and tests to `benchmark_assertions.py`.
- [/] **2.2 Decompose `recommend.py` (1455 lines)** — **[L3 | Claude 3.5]**
  - Extract target lookup/canonicalization to `target_registry.py`.
  - Isolate budget projection and headspace/matrix routing.

### P2 Stabilization Plan — 28 March 2026

- [x] Restore legacy `scripts/run_pipeline.py` entrypoint compatibility without regressing the new unified CLI.
- [x] Make `MaillardPipeline.evaluate_all()` accept both legacy dict formulations and typed `Formulation` objects.
- [x] Restore dict-like compatibility on `Formulation` for legacy tests and callers that still use `.get(...)`.
- [x] Preserve the historical pipeline → recommender concentration contract until speciation is owned consistently in one layer.
- [x] Re-run the targeted regression subset, then isolate and fix the remaining Hexanal ranking regression in the `recommend.py` / `projection.py` refactor.

Review — 28 March 2026

- Root cause of the last ranking failure was benchmark-path projection semantics, not direct drift in the extracted recommender math.
- Benchmark execution must pass explicit `process_metadata.state` into output projection.
- Exogenous injected targets must keep semantic compound names so calibration and matrix-family logic apply.
- Targeted regression subset passed after the fix: 38 tests green.

---

## Priority 3 — pH Speciation & Stealth Sinks (Scientific High-ROI)
**Goal:** Implement physiologically relevant reactivity without quantum overhead.

- [x] **3.1 Henderson-Hasselbalch Speciation** — **[L3 | Claude 3.5]**
  - Create `speciation.py` for pKa-based active fraction calculation (Lys/Cys).
- [x] **3.2 Melanoidin Sink Model** — **[L4 | Claude 4.6]**
  - Formulate continuous advanced-Maillard-coupled polymerization and trapping burden.
- [x] **3.3 Peptide / Intact-Protein Reactivity Prior** — **[L2 | Gemini 3.1]**
  - Apply heuristics for peptide-bound vs free amino reactivity.

---

## Priority 4 — Developer Experience & Usability
**Goal:** Make the codebase legible and library-ready.

- [x] **4.1 Typed Configuration Schema** — **[L2 | Gemini 3.1]**
  - Implement `@dataclass` or `Pydantic` models for `Formulation`.
- [x] **4.2 Unified CLI Framework** — **[L2 | Gemini 3.1]**
  - Migrate isolated scripts into a cohesive CLI app architecture.
- [x] **4.3 Public API Initialization** — **[L1 | Gemini 3 Flash]**
  - Create `src/__init__.py` exposing core objects.
- [x] **4.4 Code Quality Quick Wins** — **[L1 | Gemini 3 Flash]**
  - Full concrete type-hinting for `FormulationResult`.
  - Deduplicate sugar classification into `sugar_classifier.py`.
  - Archive root-level debugging scripts.

***

<details>
<summary><b>Archived Phase 1 Roadmap (Completed March 2026)</b></summary>
... [Archived content preserved]
</details>
