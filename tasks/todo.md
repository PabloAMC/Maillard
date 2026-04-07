# Maillard Strategic Roadmap

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact.

It is a problem of combining:

- a quantitatively credible free-precursor core
- matrix-aware observability and accessibility
- process-aware confidence boundaries
- scientist-facing reporting that states what is benchmarked, what is transferred, and what remains blocked by missing external evidence

The product question remains:

> Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

---

## P1 — Immediate Execution Slice (2026-04-05)

### S15. Literature Extraction Sprint — Top-20 Backlog Papers

**Rationale:** 164 of 172 scored citations (95%) remain in BACKLOG. Among those are ~30 papers scored 8/8 containing directly encodable constants (Ea values, partition coefficients, stoichiometries). The current runtime surface uses only 21 + 8 staged references. A focused extraction sprint on the highest-scored papers would approximately double the parametric surface with no wet-lab work required.

- [ ] Rank the 164 BACKLOG papers by score × family-coverage-gap to select the top-20 extraction targets.
- [ ] For each selected paper, extract deterministic constants (Ea, k, partition coefficients, stoichiometric ratios) into canonical runtime payload JSONs.
- [ ] Land extracted payloads into the appropriate operational registries (`process_state_calibrations.json`, `computational_priors.json`, `safety_reference_payloads.json`).
- [ ] Update the Deep Research backlog status from BACKLOG → RUNTIME_BOUND for each completed extraction.
- [ ] Validate that no existing benchmark regresses after each batch of new runtime payloads.

Papers with immediate high-value constants include:

- Glomb & Monnier 1995 — 3-DG retro-aldol fragmentation stoichiometry (MGO 41%, GO 28%, diacetyl 18%)
- Aliani & Farmer 2005 — ribose 3.8× MFT increase, G6P 3.2×
- Ordoudi et al. 2014 — HMF peak kinetics at pH 5.0, 125°C
- Hidalgo & Zamora 2004 — 4-HNE + Phe → 2-pentylpyrrole absolute concentrations
- Blank et al. 2001 — trans-4,5-epoxy-(E)-2-decenal from C20:4, ODT 0.07 ng/L

### S15.1. Land Staged Runtime Queue

- [x] Publish a machine-readable Deep Research runtime-first queue (8 candidates across process-state, safety, computational-prior lanes).
- [x] Land the first selected runtime-only citations from the curated queue into the operational registries without promoting any new benchmark payloads.
- [x] Advance the curated queue to the next runtime-first batch automatically once the previously staged citations are already landed.

### S15.2. Fix Cerny 2008 Benchmark Failure

**Rationale:** This is likely a data ingestion gap, not a modeling gap. The backlog contains `Cerny & Guntz-Dubini (2008)` scored 8/8 with detailed pH-resolved MFT data (thiamine-alone MFT at pH 4–8, mixed-system synergy factor 4.3×) that hasn't been encoded into the runtime. The 30.58× failure ratio is the only out-of-tolerance validation point dragging the entire surface.

- [x] Encode the Cerny 2008 key_values from the benchmark_intake_registry into the thiamine pathway priors.
- [x] Validate that the thiamine benchmark moves within the 1.5× acceptance band.
- [x] Confirm no regressions in the other 17 passing benchmarks.

### S15.3. Implement SLR-Identified Model Corrections

**Rationale:** The SLR (Section "Model corrections identified during the review") explicitly calls out 4 required corrections. Only 1 is implemented.

- [x] MFT/furfural ratio as quality constraint — implemented.
- [x] Make `volatile_retention[hexanal]` temperature-dependent (non-covalent, Ka 3.1×10²–3.1×10⁴ M⁻¹; currently scalar). Source: Ince et al. 2024.
- [x] Make `volatile_retention[MFT/FFT]` dynamic based on degree of browning (melanoidin trapping, 16× reduction). Source: Hofmann et al. 2001. (Currently deferred as 5.7.)
- [x] Add peptide-bound Cys reactivity distinction from free Cys. Source: Nishimura 2024 (17/18 consumed peptides contain Cys).

### S15.1b. Prepare The Next Runtime Queue Slice

- [x] Curate a third runtime-first batch behind the active batch so the queue can auto-advance again once batch 02 is landed.
- [x] Expose the prepared next batch in the generated `deep_research_runtime_queue` artifact instead of keeping it implicit in source only.

### S15.1c. Descriptive Test Naming And Runtime Registry Landing

- [x] Rename opaque scientific test files so filenames describe their actual benchmark or governance scope instead of using shorthand labels like `family03` or `p3`.
- [x] Land all six batch 02 runtime-only citations into the operational registries and update their Deep Research backlog state to `RUNTIME_BOUND`.
- [x] Land the first half of batch 03 with `Ordoudi et al. (2014 / PMC12484514)`, `Glomb & Monnier (1995)`, and `Aliani & Farmer (2005)` so batch 03 becomes the active queue.
- [x] Regenerate the Deep Research runtime queue artifact and validate the renamed scientific subset plus queue and registry landing tests in Docker.

Review 2026-04-05:
Batch 02 and batch 03 are now fully landed in runtime registries, the curated runtime-first queue is exhausted for the current three-batch set, and the scientific tests that previously used opaque filenames now use descriptive names.

---

### S16. README Rewrite for Food Scientists

**Rationale:** The current README (359 lines) is written for computational chemists. The target user — a food scientist in an alternative protein company — needs to see what goes in, what comes out, and how much to trust it, within 2 minutes.

- [ ] Rewrite README as a ~100-line landing page: "What this does → Install → Run → Example output → Where to learn more."
- [ ] Add a sample output snapshot (radar chart, ranked formulations table, confidence annotations) so scientists see value before installing.
- [ ] Replace internal vocabulary ("bench-neighborhood," "family-lane transparency," "validated envelope") with plain-language equivalents.
- [ ] Simplify the architecture diagram to a 3-box version ("Ingredients + Process → Maillard Engine → Predictions + Confidence") as the first visual. Move the detailed Mermaid diagram to `docs/architecture.md`.
- [ ] Consolidate the 16 validation artifact links into a single "click here to see what's validated today" reference.
- [ ] Move detailed content to dedicated documents: `docs/VALIDATION.md`, `docs/PHILOSOPHY.md`.
- [ ] Add explicit positioning statement against existing tools (NIST WebBook, RMG, manual literature review) so scientists understand the unique value.

---

### S13. Scientist-Facing Visual Output

**Rationale:** The pipeline produces rich numeric data but zero visual artifacts. Scientists need radar charts, kinetic traces, and safety dashboards to make decisions, not raw JSON. This is the single biggest product gap.

- [ ] Add matplotlib or plotly radar chart generation to `--report` output showing the multi-axis sensory profile (meaty, roasted, beany, sulfurous, sweet).
- [ ] Add a kinetic trace plot when `prediction_mode="kinetic"` to show temporal evolution of key volatiles.
- [ ] Add a safety-marker dashboard plot alongside the formulation report (acrylamide, furosine, CML/CEL).
- [ ] Add a parity plot export to the validation-figures command that a scientist can embed in presentations.
- [ ] Add a Pareto frontier visualization for multi-objective optimization (meaty-positive vs. safety vs. off-note suppression) instead of collapsing to a single blended ranking_score.

---

## P2 — High-Impact Medium-Effort

### S17. Extrusion Benchmark Validation

**Rationale:** Extrusion modeling is architecturally present (SME coupling, moisture-regime bifurcation, sequential isothermal zones, pre-extrusion damage baselines) but has zero benchmark validation. 0/2 extrusion matrices are closure-ready. For a tool aimed at alternative protein scientists, extrusion is literally the dominant production process.

#### S17.1. Extrusion Benchmark Experiment Design

- [ ] Specify the minimum viable extrusion benchmark: one protein type (PPI or SPI), two SME levels, one barrel temperature, measuring MFT + hexanal + furosine simultaneously.
- [ ] Generate a complete DoE protocol from `doe_generator.py` with real lab specifications (equipment model, SPME fiber type, exact internal standard concentrations).
- [ ] Publish the protocol as a shareable wet-lab request artifact in `results/validation/`.

#### S17.2. Extrusion Model Extensions

- [ ] Add volatile stripping correction at the die (flash-vaporization loss based on die temperature and compound vapor pressure).
- [ ] Add shear-volatile coupling beyond the simple linear SME→ΔT slope (cell-wall rupture → precursor release, protein aggregation → trapping landscape).
- [ ] Evaluate whether a simple RTD (residence time distribution) model is needed or if the sequential-zone model is sufficient for the target use case.

### S18. Selective xTB/DFT Unparking

**Rationale:** P3 refinement governance shows 0 approved jobs, but `why_not_closed.md` identifies 3 specific, narrow motif targets where xTB path search → r2SCAN-3c refinement is cost-effective and would meaningfully improve families 11, 12, and 14.

- [ ] Run xTB path search then r2SCAN-3c refinement for `hexanal_radical_quench` (Family 11: off-note suppression).
- [ ] Run r2SCAN-3c for `lysinoalanine_crosslink` (Family 12: AGE/ALE yield).
- [ ] Generate seed structure for `aa_ring_open_dicarbonyl` (Family 14: stealth browning).
- [ ] Evaluate asparagine-sugar transition state in explicit water cluster to computationally bound the matrix effect on acrylamide kinetics (narrows the wet-lab-only gap).

### S19. Web Interface (Minimal)

**Rationale:** Every interaction currently requires Docker + command-line. For food scientists, this is a major adoption barrier. A minimal web interface would increase adoption by an order of magnitude.

- [ ] Build a minimal Flask/FastAPI web interface with a formulation input form, "run prediction" button, and visual report output.
- [ ] Serve the radar chart, kinetic traces, and safety dashboard from S13 in the web response.
- [ ] Include a "download report" button for shareable PDF/HTML export.

### S14. Codebase Health & Maintainability

**Rationale:** `benchmark_validation.py` (117KB) and `recommend.py` (65KB) are monoliths that impede contribution and debugging. Test suite runtime (~1h40m) blocks iteration speed.

- [ ] Decompose `benchmark_validation.py` into modular components: registry, evaluation, reporting, assertion.
- [ ] Decompose `recommend.py` into modular components: concentration projection, observable mapping, scoring.
- [ ] Triage the test suite for performance: identify and optimize the slowest 10 tests, introduce pytest marks for fast/slow/full lanes.

---

## P3 — Strategic / Deferred

### S12. Scaling the Literature Pipeline & Uncertainty

#### S12.1: Formal Uncertainty Quantification (UQ)

- [ ] Replace narrative "trust heuristics" (e.g., Extrusion Exploratory Mode) with explicit mathematical confidence intervals (e.g., via parametric variance or Gaussian Processes) for out-of-domain predictions.
- [ ] Propagate UQ bounds into the predicted volatile headspace (ppb) figures so scientists know the exact variance of un-benchmarked estimates.

#### S12.2: Automated LLM-Assisted Payload Extraction

- [ ] Build an automated ingestion pipeline that parses eligible Deep Research summaries into canonical `benchmark_payload` JSONs to accelerate closing the ~150-paper backlog.
- [ ] Include a strict human-in-the-loop review interface to guarantee the 8-point SLR criteria are strictly maintained before merging into the main pipeline.

#### S12.3: Model-Guided Active Learning (DoE Feedback Loop)

- [ ] Formalize the "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
- [ ] Implement an API so that when the system identifies a critical gap (e.g., lack of MFT/FFT data in SPI extrudates), it auto-generates a precise wet-lab protocol optimised for maximal model calibration gain.

### Deferred Scientific Modeling Backlog

#### 5.7 Bidirectional Lipid-Maillard Crosstalk

- [ ] Add dicarbonyl-lipid oxidation promotion pathway in `lipid_oxidation.py`.
- [ ] Add melanoidin antioxidant capacity as a time-dependent LOPs suppressant.
- [ ] Validate against Report 11 crosstalk heuristics.

#### 5.8 Disulfide Bond Evolution / MFT Retention

- [ ] Model free-SH to disulfide kinetics as a function of SME and temperature.
- [ ] Link that state variable to MFT headspace recovery in the volatile retention model.

#### 5.10 Sunflower Chlorogenic Acid Off-Note

- [ ] Add temperature-gated 4-vinylguaiacol penalty for sunflower-containing formulations.
- [ ] Include chlorogenic acid to lysine covalent adduct formation as a lysine-accessibility sink.

#### 5.11 Transport / Diffusion Model for Volatile Release

- [ ] Design a 1D Fickian diffusion slab model for volatile release during cooling or serving.
- [ ] Integrate it with volatile retention factors as a compound-class-specific alternative to scalar correction.

### S9. Skipped Test Triage & QM Optionality

Status: supporting infrastructure, not current product bottleneck.

- [ ] Resume only after S14 if skipped-test cleanup blocks deterministic confidence in the active scientist workflow.

#### S9.1: Inventory and classify skipped tests

- [ ] Build a machine-readable skip registry from `pytest -rs` grouped by reason, file, and dependency class.
- [ ] Classify skips into: `not_implemented_module`, `missing_external_dataset`, `missing_optional_backend`, `long_running_campaign`.
- [ ] Add owner and unblock criteria for each skip cluster.

#### S9.2: Quasi-harmonic correction decision and implementation path

- [ ] Implement `src/quasi_harmonic_correction.py` with pure-Python deterministic numerics and no heavy backend coupling.
- [ ] Replace unconditional skips in `tests/benchmarks/test_quasi_harmonic_correction.py` with executable deterministic tests plus `skipif` only for optional integration points.

#### S9.3: Barrier and IRC benchmark skip policy

- [ ] Keep Phase 3.3/3.4 benchmark tests gated by explicit dataset/backend markers rather than unconditional `pytest.skip` inside the test body.
- [ ] Encode a run contract: default CI lane runs deterministic/unit and mock-backed checks; optional QM lane runs benchmark suites when datasets are mounted.

#### S9.4: DFT and ML-potential complement policy for skipped-test conversion

- [ ] Tie each formerly skipped QM test cluster to one of three execution lanes: `deterministic_helper_lane`, `optional_mlp_acceleration_lane`, `optional_dft_authority_lane`.
- [ ] Ensure report and validation artifacts expose lane provenance so users can distinguish deterministic numerics, MLP-assisted predictions, and DFT-confirmed evidence.

### P3–Refinement. Selective Mechanistic Refinement

Active only for matrix benchmarks that remain `mechanistic_priority` after observable closure review. Do not expand broad xTB or DFT activity beyond the specific targets in S18.

### P4. MLP Adoption

Offline accelerator lane only. External MLP evaluation must not substitute for missing matrix benchmarks or observable anchors.

### P6. Matrix-Family Expansion Beyond Pea and Soy

Keep matrix-family coverage explicit in artifacts. Do not broaden family-level scope faster than the evidence surface can support.

### Reproducibility & Provenance

- [ ] Add version-pinned reproducibility snapshots (`pip freeze` equivalent or Docker image tag) to every report artifact so predictions are exactly reproducible for peer review and regulatory contexts.
- [ ] Document the exact mapping between report provenance metadata and the reproducibility snapshot.

---

## Current Product Status

### Strong today

- Free-precursor screening is quantitatively credible inside the validated envelope (17/18 benchmarks pass, 10 strict-ready).
- 16 chemistry families tracked; 7 benchmark-linked, 4 with compound-level quantitative parity.
- 53 quantitative compound points with median ratio 1.019× in the core Maillard family.
- Family-aware ingestion, runtime, validation, and reporting are operational.
- Pea and soy matrix paths are executable and useful for directional prioritisation.
- Trust language, evidence posture, and family-lane transparency are already visible in reports.
- Intake-registry state model is normalised with explicit `triage_status`, `encoding_status`, and `runtime_artifacts_present` fields.
- Canonical literature backlog artifact exposes three exclusive queues: `ready_runtime`, `ready_benchmark`, `wet_lab_blocked`.
- Extrusion process modeling with SME coupling, moisture-regime bifurcation, and pre-extrusion damage baselines is operational.
- Safety markers (acrylamide, furosine, CML/CEL, LAL) are integrated into the pipeline.
- Cheap-first refinement screening pipeline (barrier offset sweeps with benchmark-visible decision gating) is operational.
- DoE generator templates exist for 5 gap types (blocking benchmark, missing anchor, missing kinetic, missing process-state, missing flavor anchor).
- SLR covers 990 lines with 16+ benchmark-eligible references and 5 structural gaps formally identified.

### Still blocking scientist value

- **Literature utilization**: 95% of scored citations (164/172) remain in BACKLOG; top-20 8/8 papers contain directly encodable constants.
- **No matrix benchmark** is yet external-decision-ready.
- **Cerny 2008 benchmark** (30.58× ratio) is the only out-of-tolerance validation point — likely fixable by encoding the already-scored paper's data.
- **No visual output** (charts, radar, kinetic traces) is generated from predictions — scientists must interpret raw JSON.
- **No web interface** — all interaction requires Docker + command-line, a major adoption barrier for food scientists.
- **Extrusion modeling is uncalibrated** — architecturally present but zero benchmark validation (0/2 closure-ready).
- **README is not accessible** to the target user (food scientist); uses internal vocabulary, no example output.
- Mixed pea and soy meaty-positive targets still rely on transferred or internal-candidate observable support.
- Six chemistry families still have zero benchmark-linked closure.
- SLR-identified model corrections (temperature-dependent volatile retention, melanoidin dynamic trapping, peptide-bound Cys) are not yet implemented.
- xTB/DFT pipeline is parked with 0 approved jobs despite 3 well-scoped motif targets identified.
- No version-pinned reproducibility for peer review.
- Test suite runtime (~1h40m) impedes development iteration.

---

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab or literature-backed matrix campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix usefulness improves.
- [ ] At least one matrix lane becomes meaningfully closer to external decision readiness.
- [ ] Reports and artifacts make promotion blockers explicit enough that synthetic closure is difficult.
- [ ] Expensive compute stays offline, sparse, and justified by benchmark-visible decisions.
- [ ] Visual output makes predictions actionable without requiring the scientist to interpret raw JSON.
- [ ] All benchmarks pass inside the 1.5× acceptance band.
- [ ] A food scientist can understand what the tool does and run a first prediction within 10 minutes of reading the README.
- [ ] The runtime parametric surface grows by ≥2× through literature extraction before any wet-lab work.
- [ ] Multi-objective tradeoffs (meaty vs. safe vs. off-note) are visible as Pareto frontiers, not collapsed scores.

---

## Completed Foundations

Sprints S0–S11 are complete as of 2026-04-05. All foundational architecture, family-aware ingestion and runtime, matrix observable closure, scientist experiment intake loop, family promotion contracts, intake-registry normalisation, Deep Research runtime queue publishing, extrusion process modeling (SME/moisture), safety marker integration (acrylamide, furosine, CML/CEL, LAL), cheap-first refinement screening, and DoE generator templates have been implemented and validated.
