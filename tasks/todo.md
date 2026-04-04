# Maillard Strategic Roadmap

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

## Current Execution Slice

### Active Slice: 2026-04-04 Literature Intake Priorities After Mixed-Matrix Cleanup

- [x] Confirm the visible scatter outlier is the Cerny 2008 reference-only anchor, not an experimental benchmark regression.
- [x] Standardize validation-overview benchmark labels and visually distinguish reference-only anchors from wet-lab measured comparators.
- [x] Reclassify `pmc_2026_hme_hexanal_baseline` onto its honest runtime surface: encode the extracted HME control points as Family 11 flavor-reference payloads, keep the executable-benchmark blockers explicit, and stop short of inventing final-blend pH or water activity.
- [x] Close the Family 12 intake loop for `acs_2022_pba_lysine_loss_benchmark` by wiring it to the already-encoded safety reference payload, so the learning loop reports it as encoded instead of template-required.
- [x] Run a strict literature triage for an external pea/soy mixed meaty-positive benchmark package; the SLR still resolves this as a wet-lab-only structural gap, so the blocker now stays explicit instead of being backfilled from internal lanes.
- [x] Reassess whether the next addition after those payload-state fixes should be runtime-only literature support or another benchmark lane, using the refreshed `literature_learning_loop` and objective-progress artifacts.

Strict triage verdict: no external pea/soy mixed meaty-positive benchmark package currently satisfies the benchmark contract. The next addition should therefore stay on runtime or reference-support lanes, or on benchmark lanes outside this mixed pea/soy closure gap, until a real wet-lab package exists.

### S12. Scaling the Literature Pipeline & Uncertainty

Based on a structural reflection of the matrix-literature pipeline, the strict separation of Runtime vs. Benchmark data (and the 8-point SLR) successfully prevents overfitting but fundamentally bottlenecks progress because perfect alternative-protein literature is sparse. To ensure robust high-quality predictions moving forward, we must mathematically formalize uncertainty when literature is missing, scale ingestion, and establish a closed feedback loop with wet-lab scientists.

#### S12.1: Formal Uncertainty Quantification (UQ)
- [ ] Replace narrative "trust heuristics" (e.g., Extrusion Exploratory Mode) with explicit mathematical confidence intervals (e.g., via parametric variance or Gaussian Processes) for out-of-domain predictions.
- [ ] Propagate UQ bounds into the predicted volatile headspace (ppb) figures so scientists know the exact variance of un-benchmarked estimates.

#### S12.2: Automated LLM-Assisted Payload Extraction
- [ ] Build an automated ingestion pipeline that parses eligible Deep Research summaries into canonical `benchmark_payload` JSONs to accelerate closing the ~150-paper backlog.
- [ ] Include a strict human-in-the-loop review interface to guarantee the 8-point SLR criteria are strictly maintained before merging into the main pipeline.

#### S12.3: Model-Guided Active Learning (DoE Feedback Loop)
- [ ] Formalize the "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
- [ ] Implement an API so that when the system identifies a critical gap (e.g., lack of MFT/FFT data in SPI extrudates), it auto-generates a precise wet-lab protocol optimized for maximal model calibration gain.

The product question remains:

Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

## Planning Review: 2026-03-24

This roadmap was re-audited against the current repository state before setting the next priority.

### What the current artifacts say

- `results/validation/family_deviation_audit.md` no longer shows severe residual tails in the families that already have quantitative points.
- `results/validation/family_validation_overview.md` shows 10 tracked families, but only 4 are benchmark-linked and only 4 have compound-level quantitative parity points.
- `results/validation/matrix_target_status.md` shows 5 matrix benchmarks covered, 3 with quantitative-support-ready compound closure, 0 external-decision-ready benchmarks, and 2 mechanistic-priority internal candidates.
- `results/validation/benchmark_summary.md` shows that free-precursor closure is strong, while matrix paths remain outside the strict release gate.

### Strategic conclusion

The primary bottleneck is no longer broad residual cleanup on already benchmark-backed families.

The primary bottleneck is matrix usefulness for alternative-protein scientists:

- closing observable support for decision-driving matrix compounds
- promoting at least one pea or soy mixed target benchmark from internal candidate toward externally defensible decision support
- giving scientists a machine-readable path to compare new measurements against the model and feed that evidence back into the trust surface

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact.

It is a problem of combining:

- a quantitatively credible free-precursor core
- matrix-aware observability and accessibility
- process-aware confidence boundaries
- scientist-facing reporting that states what is benchmarked, what is transferred, and what remains blocked by missing external evidence

## Current Product Status

### Strong today

- Free-precursor screening is quantitatively credible inside the validated envelope.
- Family-aware ingestion, runtime, validation, and reporting are operational.
- Pea and soy matrix paths are executable and useful for directional prioritization.
- Trust language, evidence posture, and family-lane transparency are already visible in reports.

### Still blocking scientist value

- No matrix benchmark is yet external-decision-ready.
- Mixed pea and soy meaty-positive targets still rely on transferred or internal-candidate observable support.
- Six chemistry families still have zero benchmark-linked closure.
- There is no first-class scientist data-ingestion loop that turns new measurements into reusable benchmark or calibration evidence.

## Primary Active Sprint

### Sprint S11. Matrix Observable Closure And Scientist Ingestion (Primary Active)

Goal: move the tool from matrix-directional usefulness toward matrix decision readiness for alternative-protein scientists.

This sprint replaces the earlier assumption that the highest-value work was broad residual-tail cleanup. Current evidence shows the larger product gap is matrix observability plus evidence ingestion.

### S11.A: Promote matrix target support from internal candidate to decision-usable

S11.A1: Freeze the promotion contract for matrix decision readiness
- [x] Define one explicit promotion rule from `internal_candidate` to `external_decision_ready` at the benchmark level.
- [x] Require the rule to distinguish external literature anchors, transferred calibration, internal synthetic references, and mechanistic-only support.
- [x] Publish the rule in machine-readable validation artifacts and in README trust language.

S11.A2: Close observable anchors for decision-driving mixed matrix compounds
- [x] Audit the current pea and soy mixed target panels compound by compound: sulfur positives, pyrazines, Strecker aldehydes, furans, and adverse lipid markers.
- [x] For each compound, classify the next closure action as one of: literature anchor available, class-level transfer acceptable, mechanistic blocker, or external-data blocker.
- [x] Prioritize compounds that actually drive scientist decisions rather than compounds that are merely easy to anchor.
- [x] Do not claim promotion closure for compounds that remain directional or internally constructed only.

S11.A3: Upgrade one benchmark lane, not many half-way
- [x] Choose one highest-impact matrix lane as the singular promotion target for this sprint.
- [x] Default target: pea or soy aqueous pre-extrusion mixed target ranking, because that lane combines desirable sulfur markers with adverse lipid markers in a scientist-relevant regime.
- [x] Improve that one lane until its support status is either genuinely promotable or explicitly blocked by missing external evidence.
- [x] If external evidence is insufficient, encode the blocker cleanly and stop instead of manufacturing synthetic closure.

S11.A4: Keep mechanistic refinement subordinate to the observable bottleneck
- [x] Use mechanistic refinement only where `results/validation/matrix_target_status.md` identifies a benchmark as `mechanistic_priority` after observable support is already classified.
- [x] Require any refinement candidate to name the target compounds and the expected benchmark-visible decision change before spending offline compute.
- [x] Do not reopen broad barrier retuning or generic QM cleanup unless the matrix observable audit proves it is the next limiting factor.

### S11.B: Build the scientist evidence-ingestion loop

S11.B1: Define a machine-readable experiment intake contract
- [x] Create one canonical schema for contributed matrix measurements: matrix family, process state, temperature, time, pH, formulation, measured volatiles, analytical context, and provenance.
- [x] Ensure the schema can represent partial compound panels without forcing fake completeness.
- [x] Distinguish external literature payloads from internal experiments and from synthetic diagnostic references.

S11.B2: Turn contributed measurements into trust-surface updates
- [x] Add an executable comparison artifact that answers: matched compounds, evidence posture, support delta, promotion blockers, and whether the new data changes benchmark readiness.
- [x] Make the output scientist-readable and machine-readable.
- [x] Ensure new evidence can land in benchmark payloads, calibration payloads, or explicit blocker registries rather than in narrative markdown only.

S11.B3: Expose the workflow to scientists
- [x] Add one documented scientist workflow for: run prediction, compare against a new measurement set, inspect support deltas, and decide whether the result strengthens calibration, promotion, or only hypothesis generation.
- [x] Keep this workflow lightweight enough to be used by collaborators without reading the entire repository internals.

### S11.C: Promote one non-quantitative family only if it changes matrix decisions

S11.C1: Narrow the promotion queue
- [x] Rank families 03, 04, 05, 06, 07, and 10 by direct impact on alternative-protein formulation decisions, not by narrative completeness.
- [x] Default first candidate: family 07 donor hierarchy, because it directly changes formulation choice and interacts with sulfur-positive yield.
- [x] Second candidate only if family 07 is blocked: family 03 thiamine degradation and sulfur support.

S11.C2: Require runtime landing, not literature summaries
- [x] For the chosen family, define the minimum machine-readable payload needed for first closure: benchmark payload, calibration payload, or bounded prior with explicit uncertainty.
- [x] Reject narrative-only literature additions.
- [x] Promote the family only if it changes the actual matrix decision surface or trust posture.

### S11 exit criteria

- [x] One explicit matrix benchmark promotion contract exists and is exposed in validation artifacts and trust language.
- [x] One pea or soy mixed target lane has a complete compound-level closure audit with next-action labels.
- [x] One machine-readable scientist experiment intake path exists end to end.
- [x] One new measurement set can be compared against the model with a generated support-delta artifact.
- [x] At least one non-quantitative family is either promoted with explicit bounds or documented as blocked for a concrete external-data reason.
- [x] The roadmap can explain clearly why the next bottleneck is still observable closure, or else it proves that mechanistic refinement has become the new limiter.

## Secondary Work, Explicitly Deferred

### S9. Skipped Test Triage And QM Optionality

Status: supporting infrastructure, not current product bottleneck.

- [ ] Resume only after S11 if skipped-test cleanup blocks deterministic confidence in the active scientist workflow.
- [ ] Keep quasi-harmonic correction and optional QM lane cleanup scoped as engineering hygiene, not as the main product step.

### P3. Selective Mechanistic Refinement

Status: keep active only on a narrow watchlist.

- [x] Continue only for matrix benchmarks that remain `mechanistic_priority` after observable closure review.
- [x] Do not expand broad xTB or DFT activity just because the infrastructure exists.

### P4. MLP Adoption

Status: offline accelerator lane, not main closure path.

- [x] Keep external molecular MLP evaluation bounded to approved offline roles.
- [x] Do not use MLP work as a substitute for missing matrix benchmarks or observable anchors.

### P6. Matrix-family expansion beyond pea and soy

Status: important, but second-order until one matrix lane becomes clearly decision-usable.

- [x] Keep matrix-family coverage explicit in artifacts.
- [x] Do not broaden family-level scope faster than the evidence surface can support.

## Execution Order

1. Audit the mixed pea and soy matrix target compounds against current observable support.
2. Freeze the promotion contract for `external_decision_ready`.
3. Build the scientist experiment intake schema and support-delta artifact.
4. Pick one singular matrix lane for promotion and stop splitting effort.
5. Promote one non-quantitative family only if it lands in runtime and changes the matrix decision surface.
6. Reassess whether the next blocker is still observable support or has shifted to mechanistic refinement.

## Completed Foundations

The foundational architecture remains complete and should be treated as the platform under this roadmap, not as unfinished backlog.

- [x] S0-S0d script quality and refactoring foundations
- [x] S1 chemistry-family strategy and scope policy
- [x] S2 literature generalization to family-aware payloads
- [x] S3 runtime refactor to family-aware registries and queries
- [x] S4 decision-panel expansion and family-lane metadata
- [x] S5 implementation of lipid, donor, pretreatment, and support lanes
- [x] S6 family-aware validation and reporting surfaces
- [x] S7 execution-sequence closure across ingestion, runtime, validation, and reporting
- [x] S8 definition-of-done closure for machine-readable family coverage
- [x] P1 accessibility states, process-state surrogates, and trust visibility foundations
- [x] P2 extrusion exploratory mode with explicit warnings
- [x] P5 trust surfaces across reports and comparison outputs

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab or literature-backed matrix campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix usefulness improves.
- [ ] At least one matrix lane becomes meaningfully closer to external decision readiness.
- [ ] Reports and artifacts make promotion blockers explicit enough that synthetic closure is difficult.
- [ ] Expensive compute stays offline, sparse, and justified by benchmark-visible decisions.

## Review

Review date: 2026-03-24

Why the roadmap changed:

- The previous active sprint overemphasized residual-tail reduction on benchmark-backed families.
- Current validation artifacts show that the larger product gap is no longer broad error tails inside the already quantitative surface.
- The highest-impact next step is to improve matrix observability, promotion logic, and evidence ingestion so alternative-protein scientists can trust and extend the tool in matrix-relevant regimes.

S10.A3: Promotion gate for family quantitative trust
- [x] Define family-level promotion rule: a family can be called near-quantitative only if it has minimum benchmark count plus bounded robust error.
- [x] Expose promotion state in validation markdown/json and README trust language.

Track B (second priority): Literature-to-runtime closure for families with zero quantitative parity

S10.B1: Literature closure queue by family
- [x] Rank families 03, 04, 05, 06, 07, 10 by expected decision impact and literature closability.
- [x] For each family, declare minimum viable payload needed for first benchmark-linked closure (benchmark payload, calibration payload, or bounded prior with explicit uncertainty).
- [x] Prevent narrative-only additions: every promoted literature item must land in machine-readable runtime payloads.

S10.B2: Alternative protein usability
- [x] Add explicit scientist workflow for alternative protein research (matrix choice, process state, interpretation gate, confidence posture).
- [x] Expand matrix-family coverage artifacts so users can see direct support vs transferred support by matrix family.

Track C (third priority): Optional QM/ML governance and skipped-test cleanup

- [ ] Keep S9 execution as supporting infrastructure lane.
- [ ] Do not treat S9 completion as sufficient for family predictive closure.
- [ ] Use QM/ML lanes only where Track A diagnostics prove benchmark-visible impact.

S10 measurable exit criteria:
- [x] Every benchmark-backed family has an outlier audit with top residual points and root-cause labels.
- [x] At least one high-residual family shows a verified reduction in robust error metrics after targeted fixes.
- [x] At least one currently non-quantitative family is promoted to benchmark-linked or calibration-grade with explicit uncertainty bounds.
- [x] README and validation artifacts present a scientist-operable alternative protein workflow and trust posture.

### Sprint S9. Skipped Test Triage and QM Optionality (Supporting Lane)

Goal: turn large skip clusters into an explicit policy: keep true research-stage tests gated, but promote deterministic physics helpers (like quasi-harmonic correction) into executable library tests.

S9.1: Inventory and classify skipped tests
- [ ] Build a machine-readable skip registry from `pytest -rs` grouped by reason, file, and dependency class.
- [ ] Classify skips into: `not_implemented_module`, `missing_external_dataset`, `missing_optional_backend`, `long_running_campaign`.
- [ ] Add owner and unblock criteria for each skip cluster.

S9.2: Quasi-harmonic correction decision and implementation path
- [ ] Treat quasi-harmonic correction as worthwhile if scoped to a deterministic helper layer (RRHO/Grimme/Truhlar entropy correction) rather than full QM workflow closure.
- [ ] Implement `src/quasi_harmonic_correction.py` with pure-Python deterministic numerics and no heavy backend coupling.
- [ ] Replace unconditional skips in `tests/benchmarks/test_quasi_harmonic_correction.py` with executable deterministic tests plus `skipif` only for optional integration points.
- [ ] Add a lightweight report hook that can mark whether quasi-harmonic correction was applied when barrier metadata is present.

S9.3: Barrier and IRC benchmark skip policy
- [ ] Keep Phase 3.3/3.4 benchmark tests gated by explicit dataset/backend markers rather than unconditional `pytest.skip` inside the test body.
- [ ] Encode a run contract: default CI lane runs deterministic/unit and mock-backed checks; optional QM lane runs benchmark suites when datasets are mounted.
- [ ] Add a validation note in reports/artifacts stating which QM benchmark lanes were executed versus skipped.

S9.4: DFT and ML-potential complement policy for skipped-test conversion
- [ ] Keep DFT as the authoritative physics lane for TS/IRC and high-sensitivity barrier validation; do not replace DFT acceptance criteria with MLP-only agreement.
- [ ] Keep ML potentials as bounded accelerators: candidate generation, geometry pre-relaxation, and offline triage before selective DFT confirmation.
- [ ] Tie each formerly skipped QM test cluster to one of three execution lanes: `deterministic_helper_lane`, `optional_mlp_acceleration_lane`, `optional_dft_authority_lane`.
- [ ] Ensure report and validation artifacts expose lane provenance so users can distinguish deterministic numerics, MLP-assisted predictions, and DFT-confirmed evidence.

S9 measurable exit criteria:
- [ ] `pytest -rs` shows no unconditional skip blocks for modules that are already implemented.
- [ ] Quasi-harmonic correction has executable deterministic coverage and clear activation metadata.
- [ ] Remaining skips are intentional, labeled, and linked to explicit unblock criteria.
- [ ] Every optional-QM or optional-ML skip maps to a declared lane policy and does not leak into default deterministic CI expectations.

---

## Completed Foundations (S0 through S8)

All foundational work has been implemented and integrated. Archive records below for traceability.

### Implementation Summary: S0–S8 Complete

All 8 major foundation sprints are complete as of 2026-03-23. The work comprises:

- **S0/S0b/S0c/S0d** (Script quality): All main orchestration scripts (`run_pipeline.py`, `run_cantera_kinetics.py`, `calibrate_barriers.py`) refactored, logging/config/exceptions centralized, dead code removed. ✓
- **S1** (Product positioning): Chemistry-family strategy explicit; lipid oxidation+crosstalk and donor hierarchy designated as next priority lanes; DFT and ML-potential role clarified. ✓
- **S2** (Literature generalization): Family ingestion plan frozen, all payload surfaces extended with family metadata, template vocabulary generalized to 7 canonical payload types. ✓
- **S3** (Runtime refactoring): Family-aware query layer implemented, legacy one-off getters replaced with `query_*` helpers, three runtime concepts (references/retention/priors) separated. ✓
- **S4** (Decision panel expansion): Observable panel split into 8 explicit family lanes, panel role/kind/regime metadata propagated through all projection and reporting surfaces. ✓
- **S5** (Family lanes implemented): Lipid oxidation/crosstalk, donor hierarchy, fermentation pretreatment, thiamine/nucleotide/glutathione support all wired into runtime and reporting. ✓
- **S6** (Family-aware validation): Benchmark metadata enriched with family tags, validation outputs grouped by lane, reporting surfaces expose family evidence ladder and active lanes. ✓
- **S7** (Execution record): All planned phases (A–E) completed; family ingestion/runtime/validation/reporting now unified under family-lane architecture. ✓
- **S8** (Definition of done): All criteria met; 10 numbered SLR families have explicit machine-readable lanes, no family relies on markdown-only narrative, reports explain active lanes with evidence strength, quantitative trunk intact, structural gaps remain explicit. ✓

*(Verbose closure artifacts for S0-S8 removed for brevity. See legacy commits for details on script quality refactors, early data model structures, and fundamental runtime refactorings).*

## The Three Modeling Regimes

1. Free precursors
   Goal: quantitative chemistry core
   Success means: strict-ready concentration-scale screening

2. Pea and soy matrices
   Goal: matrix target-ranking with explicit observable calibration
   Success means: useful near-quantitative prioritization for real plant proteins

3. Extrusion-heavy systems
   Goal: process-structured exploratory modeling
   Success means: useful experiment prioritization without pretending broad quantitative closure

These are not three independent backlogs. They are three trust regimes of the same tool.

---

## Program Record: Strategic Answer

### What problem should we solve?

- [ ] Make the tool reliably useful for deciding what to test next in plant-based meat flavor design, not merely for reproducing a narrow benchmark set.

### Does this reduce to modeling three classes correctly?

- [ ] Yes, but as trust regimes rather than isolated modules: free precursors, pea/soy matrices, and extrusion-heavy systems must form a coherent ladder of confidence.

### Are ML potentials the elegant main solution?

- [ ] No. The elegant main solution is benchmark-driven observable and accessibility modeling first, selective mechanistic refinement second, and external offline ML-potential acceleration only third.

## Program Record: Immediate Execution Order Under A Literature-Only Constraint

Phase L0: Freeze the scientist decision panel and evidence contract

- [x] Convert the current SLR intake into a canonical machine-readable decision panel covering meaty positives, Strecker aldehydes, pyrazines, furans/furanones, lipid off-notes, and safety markers.
- [x] Assign every panel compound one explicit evidence state: benchmarked, calibration-only, transferred prior, mechanistic surrogate, or missing.
- [x] Make reports, validation summaries, and optimizer outputs consume that same evidence contract so scientists can see what is truly supported.

Phase L1: Promote all literature-closable matrix surfaces before adding new mechanics

- [x] Turn the already eligible pea and soy literature records into executable benchmark or calibration artifacts wherever the SLR says the data are sufficient.
- [x] Encode Trikusuma, Pratap-Singh, Asen, Li, Squeo, Lincoln, and other already reviewed references into runtime-facing payloads instead of leaving them as narrative-only support.
- [x] Distinguish explicitly between what the literature closes today and the structural gaps that still require new primary data.

Phase L2: Close matrix accessibility and process priors without pretending they are benchmarks

- [x] Add peptide-bound and process-opened accessibility priors from literature where quantitative windows exist, even if they remain calibration-grade rather than benchmark-grade.
- [x] Add mycoprotein as a first-class matrix family if the literature supports bounded priors, even if early support is directional rather than quantitative.
- [x] Keep all such priors tagged with provenance, uncertainty posture, and process-state applicability.

Phase L3: Build a literature-only learning loop

- [x] Create a reproducible intake-to-runtime workflow so each newly accepted paper updates one or more of: benchmark payloads, calibration registries, evidence maps, or trust reports.
- [x] Add review artifacts that show which structural gaps are shrinking through literature ingestion and which remain impossible without new experiments.
- [x] Ensure the user can benefit from continuous literature updates even if no new wet-lab data arrive.

Phase L4: Open P3 only after L0-L3 produce a stable decision surface

- [x] Start mechanistic refinement only after the literature-backed matrix decision surface is explicit enough to say which remaining failures are truly mechanistic.
- [x] Require every expensive refinement candidate to show expected benchmark or ranking impact before execution.
- [x] Keep the daily workflow cheap and benchmark-facing while expensive physics remains offline and sparse.

### P0. Make Matrix Predictions Scientist-Useful

Goal: promote pea and soy from intake checks to matrix target-ranking systems.

Phase P0.1: Define the benchmark target panel

- [x] Freeze the compound panel for matrix target-ranking: sulfur positives, Strecker aldehydes, pyrazines, furans/furanones, and adverse lipid markers.
- [x] Map each panel compound to one of four evidence states: externally benchmarked, internally benchmarked, transferred prior, or still missing.
- [x] Publish the target panel and evidence map in machine-readable form so reports and validation use the same contract.

Phase P0.2: Build benchmark-grade pea and soy targets

- [x] Create at least one pea diagnostic candidate and one soy diagnostic candidate where desirable and adverse compounds are ranked together under a single process state.
- [ ] Ensure each benchmark has explicit measured values, ranking intent, adverse markers, and scale thresholds where justified.
- [x] Add focused tests proving that the benchmark executes, matches all target compounds, and preserves the intended ordering.
- [x] Treat mixed matrix candidates without external measurements as mechanistic triage inputs, not promotion evidence.

---

## Priority 5 — Strategic Gaps (Scientific Assessment — March 2026)
**Goal:** Close the gap between the tool's current validated envelope and the needs of working alternative protein scientists, using literature-closable computational methods.

### Tier 1 — High-impact, literature-closable (no wet-lab)

#### 5.1 Extrusion Process Model — **[L3 | Implemented with sequential isothermal zone model]**
The dominant commercial PBMA process (extrusion) is the weakest regime in the tool. Report 12 provides the calibration data.

- [x] **5.1a SME as independent process variable**: Added `sme_kj_per_kg` to `ReactionConditions`; effective extrusion temperature now applies an SME-driven correction with a 5–40°C bounded offset.
- [x] **5.1b Moisture regime classifier**: Added `moisture_regime: Literal['lme', 'hme']`-style runtime handling in conditions and flipped `predict_acrylamide()` moisture dependence between LME and HME.
- [x] **5.1c Pre-extrusion damage base load**: Added `pre_extrusion_damage_load(protein_type)` returning baseline furosine and LAL loads for concentrate/isolate regimes.
- [x] **5.1d Autoclave sterilization step**: Added a discrete sterilization damage increment (121–126°C, 15–30 min) layered separately from the extruder barrel.
- [x] **5.1e Spatial discretization**: Scoped to a sequential isothermal zone model for the barrel and surfaced through extrusion process metadata.

#### 5.2 Protein Source Registry — **[L3]**
Report 06 produced meaty potential multipliers for 14 protein sources. None are encoded at runtime.

- [x] Create `data/lit/protein_source_registry.json` with AA composition, meaty potential multiplier, off-note penalty, LOX activity flag, and methoxypyrazine ceiling per source.
- [x] Wire registry into `matrix_correction.py` so switching protein source auto-adjusts correction factors.
- [x] Add CLI `--protein-source` flag to `run_pipeline.py` that selects from the registry.
- [x] Include engineering heuristics from Report 06 (e.g., methoxypyrazine non-correctable flag for pea >50%, wheat gluten Cys advantage for MFT).

#### 5.3 Ingest Benchmark-Eligible Deep Research Data — **[L2]**
Over 150 benchmark-eligible datasets exist in the Gemini Deep Research corpus (~176 tracked) but aren't wired into the runtime.

- [x] **Deep Research Tracker Implementation**: Created `scripts/deep_research_tracker.py` to continuously parse `.md` files and track ingested vs backlog state.
- [x] **SPI-HVP + xylose**: Ingested as `data/benchmarks/spi_hvp_xylose_120C_PMC9905368.json` with runtime protein-source wiring and Docker-validated benchmark execution.
- [x] **Wheat gluten HVP + xylose**: Ingested as `data/benchmarks/wheat_gluten_hvp_xylose_120C_PMC9905368.json` and validated as the stronger MFT plant-source benchmark.
- [x] **Acrylamide fast kinetics**: Added an extrusion-context benchmark and routed safety benchmarking through effective extrusion temperature so the 22.36→62.62 µg/kg short-residence point is runtime-bound.
- [x] **CML/CEL commercial PBA ranges**: Added `predict_cml()` and `predict_cel()` proxy surfaces plus a Foods 2023 benchmark and safety reference payload bindings.
- [x] **Furosine formation-elimination crossover**: Added `predict_furosine()` proxy behavior, a crossover benchmark, and reference payload wiring so the post-150°C decline is explicit rather than hidden.

#### 5.3 Next Runtime Queue — **[L2 | Active]**
The tracker now needs a machine-readable family queue so we stop selecting the next ingestion tranche manually.

- [x] **5.3q1 Family-prioritized backlog artifact**: Extend `family_ingestion_plan` so it merges Deep Research backlog counts with family posture, implementation wave, and runtime module targets.
- [x] **5.3q2 Docker entrypoint**: Expose a named Docker command for generating `results/validation/family_ingestion_plan.{md,json}` from the current audit state.
- [x] **5.3q3 Next closure slice**: Execute Family 11 `lipid_maillard_crosstalk` as the next benchmark-facing slice, with Family 02 lipid adverse-marker vs carbonyl-competition cleanup as the prerequisite cut. Family 02 now exposes Trikusuma 2019 benchmark anchors and target values directly in the runtime lane, and Family 11 now exposes explicit hexanal/MFT competition priors plus xTB/HME anchor metadata in `build_flavor_axis_summary()`.
- [x] **5.3q3a Family 11 HME intake hardening**: Corrected the Foods 2026 HME intake entry to the exact paper (`10.3390/foods15050912`), encoded the confirmed pretreatment and extrusion method details from `data/lit/foods-15-00912.pdf`, and kept the executable-benchmark blockers explicit (`final-blend pH` and direct `water_activity/aw` still missing).
- [x] **5.3q3b Benchmark condition alias cleanup**: `src/benchmark_validation.py` now accepts both canonical `conditions.water_activity` and alias `conditions.aw`, so benchmark payloads and matrix-only paths do not silently diverge on moisture semantics. Docker validation passed on the touched subset.
- [ ] **5.3q3c Family 11 closure follow-up**: Find a companion source for the same SPI/wheat-gluten HME system that reports final-blend pH or direct `aw`; if none exists, keep `pmc_2026_hme_hexanal_baseline` intake-only and route it only as a runtime/reference anchor.
   - [x] Review candidate PDF `data/lit/1-s2.0-S0260877423001632-main.pdf` for exact-system pH/aw closure support. Result: useful process analogue for WG hydrolysate plasticization in SP-WG HME, but still no direct `aw` and no final extrudate pH for the exact Family 11 volatile anchor system.
   - [x] Review candidate PDF `data/lit/foods-12-00912.pdf` for exact-system pH/aw closure support. Result: useful HSPI/SP-WG HME structural-process analogue, but it does not report direct `aw` and does not close final-blend pH for the Family 11 control.
   - [x] If neither candidate closes Family 11, move immediately to the next literature-ingestion priority instead of stretching the HME claim.
- [x] **5.3q4 Follow-on support lane**: Family 03 is now promoted from availability-only to calibrated sulfur-support runtime logic. The runtime lane uses Hofmann 1996 as a conservative synergy floor, Cerny 2008 as the mixed-system pH-calibrated benchmark anchor, and De Leyn 2019 as the extrusion-survival attenuation anchor; focused Docker validation passed on `tests/scientific/test_family03_benchmark.py` and `tests/unit/test_literature_runtime.py` (35 passed).

Review 2026-04-03:
- The Foods 2026 HME paper is now encoded as a stronger Family 11 intake anchor, but still does not justify executable benchmark promotion.
- Local repo search did not surface a strong in-repo companion source already carrying the missing final-blend pH or water activity for this exact HME system.
- Family 03 support-lane promotion is already landed and now validated as complete at the runtime and benchmark level; the next unresolved literature-ingestion blocker remains Family 11 closure, not thiamine support.
- The two newly downloaded HME-adjacent candidates (`1-s2.0-S0260877423001632-main.pdf` and `foods-12-00912.pdf`) help on process analogues and moisture framing, but neither closes the exact Family 11 benchmark blockers (`final-blend pH`, direct `aw`).

#### 5.4 Kinetic Mode Documentation & User Exposure — **[L1]**
The ODE kinetics engine is a P1 deliverable but invisible to users.

- [ ] Add `prediction_mode="kinetic"` to `QUICKSTART.md` with example usage.
- [ ] Add a section to `food_scientist_walkthrough.md` showing kinetic mode for temporal extrusion dynamics.
- [ ] Update `SCIENTIFIC_RELIABILITY.md` to describe what kinetic mode improves over projection mode.
- [ ] Add `--prediction-mode kinetic` to the CLI help text prominently.

#### 5.5 Enhanced Onboarding for Food Scientists — **[L2]**
Current onboarding is too terse and assumes computational background.

- [ ] **Decision tree**: Create a visual flowchart "I have X protein, Y process → use these commands → interpret these numbers → trust them this much."
- [ ] **Annotated output walkthrough**: Add real (not hypothetical) example output to `food_scientist_walkthrough.md` with annotations explaining each field.
- [ ] **Glossary linkage**: Add inline glossary links from workflow guides to the README glossary.
- [ ] **Interactive notebook**: Create `docs/notebooks/formulation_comparison.ipynb` walking through a pea vs. soy comparison.

### Tier 2 — Medium-impact, improves scientific accuracy

#### 5.6 Explicit Carbonyl Donor Hierarchy — **[L3]**
Ribose, xylose, glucose are not interchangeable (ribose generates 5–10× more MFT). Currently collapsed to near-generic treatment.

- [x] Add sugar reactivity multipliers to `barrier_constants.py` keyed on donor identity.
- [x] Wire donor identity into the shared kinetic rate path (`src/conditions.py`, `src/kinetics.py`, `src/pipeline.py`, `src/benchmark_validation.py`) so reactant labels can penalize slow hexose routing without inflating already-calibrated pentose baselines.
- [x] Add benchmark test comparing ribose vs glucose MFT yield at equivalent conditions.

Review 2026-04-03:
- Carbonyl-donor identity is now passed from enumerated reactant labels into the rate-constant path, keeping ribose as the calibrated baseline and expressing hierarchy mainly as a glucose penalty plus a modest phosphorylated-donor uplift.
- Focused validation passed for the new donor-hierarchy slice: `tests/unit/test_barrier_constants.py::TestDonorReactivityMultipliers`, `tests/unit/test_conditions.py::test_rate_constant_respects_carbohydrate_donor_identity_context`, and `tests/scientific/test_free_aa_quantitative_regression.py::test_ribose_beats_glucose_for_equal_condition_mft_prediction`.
- A broader pre-existing regression still remains in `tests/scientific/test_free_aa_quantitative_regression.py::test_primary_free_amino_acid_benchmarks_stay_locally_calibrated` for `cys_ribose_150C_Mottram1994`; it was not reopened here because this slice was constrained to donor-discrimination wiring rather than full free-AA recalibration.

#### 5.7 Bidirectional Lipid-Maillard Crosstalk — **[L4]**
Currently one-directional (LOPs → Maillard). Missing: dicarbonyl → lipid oxidation catalysis, melanoidin → antioxidant protection.

- [ ] Add dicarbonyl–lipid oxidation promotion pathway in `lipid_oxidation.py`.
- [ ] Add melanoidin antioxidant capacity as a time-dependent LOPs suppressant.
- [ ] Validate against Report 11 crosstalk heuristics.

#### 5.8 Disulfide Bond Evolution / MFT Retention — **[L3]**
Free-SH consumption during extrusion constrains MFT flavor retention (disulfide interchange).

- [ ] Model free-SH → disulfide kinetics as function of SME and temperature.
- [ ] Link to MFT headspace recovery in volatile retention model.

### Tier 3 — Polish and completeness

#### 5.9 Visual Output in Main Prediction Workflow — **[L2]**
Scientists want charts, not JSON blobs.

- [ ] Add matplotlib/plotly radar chart generation to `--report` output.
- [ ] Add kinetic trace plot (concentration vs. time) when `prediction_mode="kinetic"`.
- [ ] Add safety marker dashboard plot alongside formulation report.

#### 5.10 Sunflower Chlorogenic Acid Off-Note — **[L2]**
Temperature-triggered 4-vinylguaiacol generation from chlorogenic acid degradation (>120°C). Documented in Report 06 but unmodeled.

- [ ] Add temperature-gated 4-vinylguaiacol penalty for sunflower-containing formulations.
- [ ] Include chlorogenic acid → lysine covalent adduct as a lysine accessibility sink.

#### 5.11 Transport / Diffusion Model for Volatile Release — **[L4 | Long-term]**
No mass transfer model exists for volatile release from protein matrices.

- [ ] Design a 1D Fickian diffusion slab model for volatile release during cooling/serving.
- [ ] Integrate with volatile retention factors as a compound-class-specific alternative to scalar correction.

***

- [ ] Add observable anchors for the compounds scientists actually use to decide, not only hexanal-like off-notes.
- [x] Separate compound classes that are quantitatively anchored from those that remain directional or transferred.
- [x] Make projection metadata expose which matrix compounds are directly anchored versus transferred.

Phase P0.4: Promote matrix validation from intake-only to target-ranking

- [x] Extend matrix benchmark contracts so they score multi-compound target-ranking behavior, not just executable intake behavior.
- [x] Add summary artifacts that show which matrix targets are quantitatively closed, directionally supported, or still open.
- [x] Define a measurable promotion rule for moving a matrix benchmark family from directional to near-quantitative support.

Phase P0.5: Make matrix-family coverage explicit beyond pea and soy

- [ ] Create a canonical matrix-family coverage registry that distinguishes chemistry-core support from matrix-family support for at least: free precursors, pea isolate, soy isolate, soy hydrolysate, mycoprotein, extrusion-heavy systems, and lipid-rich co-matrices such as coconut oil.
- [ ] Mark each family with one explicit runtime posture: quantitative_core, directional_matrix, qualitative_intake_only, indirect_generic_support, or open_gap.
- [ ] Encode the difference between explicit family support and indirect generic support so fat_fraction or generic lipid-trapping logic is not misread as coconut-oil validation.
- [ ] Link every family to its current evidence surface: executable benchmark, calibration prior, qualitative intake, mechanistic-only chemistry, or structural external-data gap.
- [ ] Add a scientist-facing artifact and reporting surface so the repo can answer "what matrix families are actually covered today?" without relying on chat summaries or institutional memory.
- [ ] Make coconut-oil and other lipid-rich co-matrix gaps explicit as product-scope gaps, not hidden inside generic lipid chemistry.

P0.5 measurable exit criteria:

- [ ] The repo can enumerate which matrix families are directly supported versus only indirectly approximated.
- [ ] Coconut oil or an equivalent lipid-rich co-matrix appears explicitly as a tracked open family-level gap.
- [ ] Reports and scientific references expose the coverage artifact so newcomers do not confuse pea/soy progress with broad PBMA matrix closure.

P0 measurable exit criteria:

- [ ] At least one pea and one soy matrix benchmark rank sulfur, Strecker, pyrazine, and adverse markers together.
- [x] All compounds in the panel have an explicit evidence label in reports.
- [x] Matrix summary artifacts distinguish quantitative support from directional support at the compound level.
- [x] Matrix summary artifacts explicitly separate external decision readiness from mechanistic-priority candidates.

### P1. Model Accessibility, Not Just Chemistry

Goal: capture the main reason plant matrices differ from free systems.

Phase P1.1: Encode accessibility states explicitly

- [x] Define a small set of accessibility states that the runtime understands: free, peptide-bound, protein-embedded, and process-opened.
- [x] Map current pea and soy assumptions onto these states instead of burying them in one-off correction factors.
- [x] Add tests showing that changing accessibility state changes only the intended observables.

Phase P1.2: Build reusable process-state surrogates

- [x] Add reusable surrogates for denaturation, release, retention, and browning severity that depend on temperature, time, pH, and matrix family.
- [x] Keep these surrogates cheap enough for daily runtime and explicit enough to surface in reporting.
- [ ] Validate that process-state changes improve matrix benchmark behavior without destabilizing free-precursor benchmarks.

Phase P1.3: Couple accessibility to prediction outputs

- [x] Propagate accessibility state and effective process state into the projection metadata used by reports.
- [x] Surface chemically reachable versus merely plausible compounds in scientist-facing artifacts.
- [x] Add domain-of-applicability warnings when accessibility assumptions dominate the prediction.

P1 measurable exit criteria:

- [x] Accessibility state is a first-class runtime concept rather than an implicit correction.
- [ ] Process-state surrogates are benchmarked for pea and soy on the compounds that drive decisions.
- [x] Reports show when accessibility, not chemistry, is the main source of uncertainty.

### P2. Open The Extrusion Regime Carefully

Goal: make extrusion-heavy predictions useful without overclaiming quantitative closure.

Phase P2.1: Define the extrusion neighborhood

- [x] Write an explicit extrusion-heavy benchmark neighborhood definition: what counts as in-domain, near-domain, and out-of-domain.
- [x] Identify the minimum observable panel needed for extrusion usefulness: at least key meaty positives, key off-notes, and one process-severity marker.
- [x] Encode default warnings so extrusion outputs cannot be misread as free-system quantitative predictions.

Phase P2.2: Add process-structured extrusion surrogates

- [x] Add cheap surrogates for moisture redistribution, accessibility loss/recovery, and volatile retention under high-severity processing.
- [x] Reuse P1 accessibility and process-state concepts instead of creating a separate extrusion-only stack.
- [x] Ensure extrusion surrogates degrade confidence cleanly rather than silently reusing free-precursor assumptions.

Phase P2.3: Create exploratory extrusion mode

- [x] Add an explicit exploratory recommendation mode for extrusion-heavy systems.
- [x] Report which outputs are still benchmark-backed from lower regimes and which are new extrusion extrapolations.
- [x] Add acceptance tests that confirm extrusion mode always carries the intended warnings and confidence posture.

P2 measurable exit criteria:

- [x] Extrusion-heavy runs use an explicit exploratory mode with visible warnings.
- [x] The runtime reuses shared accessibility/process abstractions rather than duplicating logic.
- [x] Scientists can distinguish reused lower-regime support from true extrusion-specific support in the default reports.

### P3. Refine Only Decisive Mechanistic Gaps

Goal: spend expensive compute only where it changes decisions.

Phase P3.0: Build the targeting layer first

- [x] Add motif-class tagging for elementary steps so refinement candidates are grouped by chemistry family rather than inspected manually one by one.
- [x] Implement global sensitivity analysis over formulation, process, and barrier families so the repo can rank which motif classes actually drive matrix ranking error or benchmark drift.
- [x] Produce a refinement watchlist that links each candidate motif class to observable impact, current evidence level, and expected value of refinement.
- [x] Use matrix target status next actions to distinguish external-data blockers from mechanistic-priority candidates before launching new offline refinement work.
- [x] Materialize an initial barrier-family sensitivity artifact over benchmark-visible systems so offline refinement can start from observed impact rather than occurrence counts alone.

Phase P3.1: Run cheap-first mechanistic triage

- [x] Run xTB or similarly cheap screening first on the top-ranked motif classes only.
- [x] Reject candidates that do not move the target observable panel, the benchmark ordering, or the calibration diagnostics.
- [x] Store all cheap-first outputs in reusable artifacts with provenance and uncertainty rather than ad hoc notebooks or one-off files.

Phase P3.2: Escalate only the decisive subset to selective DFT

- [x] Escalate only the small subset of candidates that remain benchmark-relevant after cheap screening.
- [x] Require an explicit rationale for solvent model, charge, spin, and IRC needs before launching DFT.
- [x] Compress every accepted DFT result into reusable barrier or delta artifacts consumable by the cheap daily workflow.

Phase P3.3: Prove decision impact before keeping the refinement

- [x] Re-run the affected benchmark and reporting surfaces after each accepted refinement.
- [x] Keep only refinements that improve ranking, scale, or confidence posture in a benchmark-visible way.
- [x] Revert or quarantine refinements that add cost without improving scientist-facing decisions.

P3 measurable exit criteria:

- [x] The repo can rank candidate motif classes for refinement using a reproducible sensitivity workflow.
- [x] Every accepted xTB or DFT artifact is written back into the cheap daily workflow with provenance.
- [ ] At least one benchmark-visible or panel-visible decision improves because of selective mechanistic refinement.

### P4. Keep MLPs In Their Proper Role

Goal: use ML potentials where they are elegant, not where they are fashionable.

Phase P4.0: Gate P4 against the real bottleneck

- [x] Treat P4 as an offline acceleration track, not as the main closure path for plant-matrix accuracy; benchmark-grade matrix observability and external data gaps remain separate blockers.
- [x] Only open a new P4 work item when P3 identifies a benchmark-visible mechanistic gap that is not better explained by missing matrix data, observability transfer, or process-state calibration.
- [x] Require every proposed MLP use case to name the target motif family, the benchmark-visible decision it could change, the expected speedup, the likely failure mode, and the non-ML comparator.

Phase P4.1: Build a chemistry-domain benchmark before any adoption

- [x] Replace endpoint-energy-only benchmarking with a reaction-centered benchmark set anchored to existing Tier 2 structures, TS searches, and barrier deltas for Maillard-relevant sulfur, carbonyl, proton-transfer, and fragmentation chemistry.
- [ ] Benchmark chemistry-first candidates such as MACE organic variants, OrbMol, and AIMNet2 before any materials-first or generic foundation model is considered.
- [x] Add explicit stop rules: if a candidate produces nonphysical energies, fails TS/geometry sanity checks, or cannot preserve barrier ordering on the benchmark set, quarantine it from workflow integration.

Phase P4.2: Restrict MLPs to narrow offline roles first

- [x] Restrict early MLP usage to conformer screening, geometry pre-optimization, TS initialization, or local surrogate ranking for already-approved motif families.
- [x] Keep MLP outputs out of the default prediction surface until they have been collapsed into validated cached artifacts with provenance and uncertainty.
- [x] Require every accepted MLP-assisted artifact to be reproducible with a higher-fidelity fallback path such as xTB plus selective DFT.

Phase P4.3: Adopt only benchmarked accelerators

- [x] Add a benchmarked adoption note for any selected MLP, including where it improves throughput, where it fails, and which chemistry families remain out of domain.
- [x] Expose in reports and provenance when a result depends on computational refinement or cached QM/ML surrogates rather than direct literature or benchmark support.
- [x] Keep a no-default-MLP posture until at least one candidate demonstrates stable value on the reaction-centered benchmark without degrading scientist-facing trust surfaces.

P4 measurable exit criteria:

- [x] No MLP is used as a substitute for missing matrix benchmarks.
- [x] Any adopted MLP serves a clearly bounded offline role with explicit provenance and a higher-fidelity fallback.
- [x] The chosen MLP shortlist is justified against chemistry-domain benchmark results rather than generic SOTA branding.
- [ ] At least one MLP-assisted offline role reduces compute cost or turnaround while preserving barrier ordering and geometry sanity on the accepted benchmark set.

### P6. Close Product-Scope Blind Spots

Goal: ensure the repo does not overfit its scientific product story to pea/soy isolates when important plant-based meat systems depend on other co-matrices and formulation families.

Phase P6.0: Audit matrix-family scope before expanding mechanics

- [ ] Audit the active runtime and validation surfaces to separate protein families, lipid-rich co-matrices, hydrolysate systems, and process-heavy systems.
- [ ] Distinguish chemistry families already modeled generically from matrix families that still lack explicit calibration or validation.
- [ ] Use that audit to prevent roadmap language from implying that pea/soy work equals broad plant-based meat coverage.

Phase P6.1: Prioritize co-matrices by decision impact

- [ ] Rank missing matrix/co-matrix families by how often they appear in real PBMA formulations and how strongly they alter observable tradeoffs.
- [ ] Treat coconut oil and similar fat phases as first-class scope questions because they change release, oxidation, and sulfur-versus-off-note tradeoffs even when they are not the main amino-acid source.
- [ ] Separate families that need family-specific lipid composition/retention models from those that only need generic directional warnings.

Phase P6.2: Expand one missing family explicitly before broadening further

- [ ] Start with one high-value missing family-level surface, likely a lipid-rich co-matrix lane, rather than diluting effort across many low-evidence proteins at once.
- [ ] Require the first expansion to encode both what the repo can already reuse from generic lipid chemistry and what still needs explicit family-specific evidence.
- [ ] Keep any new family outside the strict gate until it has its own benchmark or bounded calibration surface.

### P5. Make Trust Visible In Every Decision Surface

Goal: ensure a scientist can see why a prediction should or should not be trusted.

- [x] Surface benchmark neighborhood, calibration source, observable-layer assumptions, and extrapolation axes in all default reports.
- [x] Show which compounds are benchmark-backed, transferred, heuristic, or computationally refined.
- [x] Distinguish clearly between quantitative recommendation mode and directional hypothesis mode.

## Recent Verification

- [x] Chemistry-family scope now generates data/lit/chemistry_family_scope_registry.json plus results/validation/chemistry_family_scope.{md,json}, making explicit which families are first-class core, bounded lanes, high-priority partial lanes, or open gaps, and recommending lipid oxidation plus carbonylic crosstalk as the next product-scope expansion.
- [x] Docker validation passed on the chemistry-family scope suite: 11 tests green across the new artifact, reporting-surface exposure, and existing SLR payload exposure.
- [x] Explicit decision-panel evidence contract propagated through projection metadata, benchmark snapshots, confidence rows, and reporting surfaces.
- [x] Runtime reports now include a literature evidence summary derived from the existing benchmark intake registry and structural-gap payloads.
- [x] Docker validation passed on the affected suite: 36 tests green across benchmark targets, matrix delta reports, SLR payloads, usability reports, trust visibility, and projection contract coverage.
- [x] Literature-only learning loop now generates results/validation/literature_learning_loop.{md,json} and results/validation/literature_runtime_templates.json from the canonical intake registry.
- [x] Matrix priors now surface uncertainty posture and process-state applicability, and mycoprotein is encoded as a first-class bounded matrix family.
- [x] Docker validation passed on the expanded L1/L2/L3 suite: 66 tests green across the new learning-loop coverage, priors, reporting, matrix targets, and trust surfaces.
- [x] Matrix target support status now generates results/validation/matrix_target_status.{md,json} and separates quantitative closure, internal candidates, directional support, and open gaps at the compound level.
- [x] Matrix family coverage now generates data/lit/matrix_family_coverage_registry.json plus results/validation/matrix_family_coverage.{md,json}, explicitly separating direct matrix-family support from indirect generic support and making the coconut-oil co-matrix gap visible.
- [x] Refinement watchlist now generates results/validation/refinement_watchlist.{md,json} plus results/validation/offline_dft_jobs.json to keep expensive refinement benchmark-facing and offline.
- [x] Docker validation passed on the expanded L4/P0/P1 suite: 71 tests green across matrix status, refinement watchlist, benchmark evidence, target snapshots, accessibility, reporting, and projection contracts.
- [x] Reaction-centered benchmarking now materializes data/lit/reaction_benchmark_set.json plus results/validation/mlp_reaction_benchmark.{md,json}.
- [x] Geometry benchmarking now materializes data/lit/geometry_benchmark_set.json plus results/validation/mlp_geometry_benchmark.{md,json} and results/validation/mlp_geometry_assessment.{md,json}.
- [x] MLP governance now materializes data/lit/mlp_candidate_registry.json plus results/validation/mlp_assessment.{md,json} and results/validation/mlp_adoption_notes.{md,json}, with ResultsDB persistence for adoption decisions.
- [x] The external evidence lane now also materializes data/lit/mlp_external_benchmark_evidence.json plus results/validation/mlp_external_mlp_landscape.{md,json}, making external community evidence explicit as shortlist priority rather than as a substitute for local Maillard validation.
- [x] The current P4 outcome is explicit: mace_off_medium is quarantined for barrier-surrogate use, while mace_mp_small is accepted only for bounded offline geometry pre-optimization with fallback to r2SCAN-3c plus wB97M-V.
- [x] The current shortlist policy is explicit: AIMNet2, OrbMol, and MACE-OMOL are high-priority chemistry-first candidates because of external molecular evidence, but they remain deferred locally until geometry or TS benchmark evidence exists on Maillard-relevant systems.
- [x] Docker validation passed on the expanded P4 suite: 23 tests green across geometry benchmarks, reaction benchmarks, chemistry benchmark validation, reporting surface exposure, and MLP adoption persistence.
- [x] Docker validation passed on the external-evidence P4 suite: 17 tests green across the new external landscape artifact, chemistry benchmark validation, and reporting exposure.
- [x] Docker validation passed on the matrix-family coverage suite: 11 tests green across the new family-coverage registry, generated artifact, and reporting-surface exposure.
- [x] P3 global sensitivity now generates results/validation/p3_global_sensitivity.{md,json} and ranks barrier, process, and formulation axes on benchmark-visible systems.
- [x] Cheap-first refinement screening now generates results/validation/cheap_refinement_screening.{md,json} and explicitly quarantines candidates that are sensitivity-visible but not benchmark-improving.
- [x] Selective DFT planning now generates results/validation/selective_dft_plan.{md,json} plus results/validation/p3_offline_dft_jobs.json, including explicit no-escalation decisions when cheap-first does not justify QM.
- [x] Refinement impact now generates results/validation/refinement_impact.{md,json} and writes data/lit/refinement_surrogate_patches.json back into the cheap daily workflow.
- [x] Docker validation passed on the P3 suite: 24 tests green across refinement campaign, family sensitivity, matrix target status, refinement watchlist, DFT contract, and usability reporting.
- [x] Trust surfaces now expose decision mode, reachability, observable assumptions, calibration source, and extrapolation axes across CLI, single-run reports, comparison reports, and confidence rows.
- [x] Docker validation passed on the focused P5/P1 trust suite: 17 tests green across usability reports and trust visibility.

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix modeling expands.
- [ ] Pea and soy predictions become target-ranking useful for multi-compound flavor decisions.
- [ ] Extrusion-heavy systems become operationally useful as exploratory planning tools.
- [ ] Expensive computation remains offline, sparse, and reusable.
- [ ] Reports state confidence and extrapolation explicitly enough that overclaiming becomes difficult.
