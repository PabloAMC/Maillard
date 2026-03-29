# Maillard Strategic Roadmap: Phase 2

## North Star
Evolve the Maillard framework from a static, heuristic-heavy decision surface into a dynamic, time-resolved kinetic engine, while aggressively paying down architectural debt that threatens long-term maintainability.

---

## Priority 1 — Time-Resolved Microkinetics (The Major Scientific Upgrade)
**Goal:** Transition core prediction to a true time-resolved ODE engine.

- [ ] **1.1 Wire Cantera deeply into `pipeline.py`** — **[L4 | Claude 4.6]**
  - Implement a `mode="kinetic"` flag in the main pipeline.
  - Bridge `smirks_engine.py` outputs into the Cantera mechanism generator dynamically.
- [ ] **1.2 Expose Temporal Dynamics in `recommend.py`** — **[L4 | Claude 4.6]**
  - Update `predict_from_steps` to parse Cantera time-series data.
  - Capture peak vs. plateau states for key markers.
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
