# Maillard Strategic Roadmap

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

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

### S0–S0d: Script Quality Improvements (Completed)

**Goal**: Refactor main orchestration scripts and deep library code quality without altering scientific logic.

**Implementation closure**:
- The 2026-03-23 audit found no remaining unchecked S0/S0b/S0c/S0d implementation work in the active codebase.
- All three main scripts modularized with clear function boundaries, docstrings, type hints.
- Logger, config, and typed exceptions centralized in `src/logger.py`, `src/config.py`, `src/exceptions.py`.
- Removed `sys.exit(1)`, consolidated `_normalize_name` and `_canon` into `src/text_utils.py` and `src/chem_utils.py`.
- Removed dead interactive `__main__` stubs, deprecated legacy `predict()` flow.
- Scripts use logging instead of stdout; no hardcoded magic constants.
- Tests pass via `docker_maillard.sh pytest`.

### S1: Chemistry-Family Positioning (Completed)

**Goal**: Expand beyond amino acid–sugar without losing the core; explicitly name which families are core/partial/bounded/gap.

**Implementation closure**:
- `src/family_strategy_policy.py` derives machine-readable strategy, fixes lipid oxidation and carbonylic crosstalk as next expansion.
- Same artifact classifies families into first-class core, high-priority partial, bounded, and open gaps.
- P4 archived as infrastructure; DFT reserved for TS-sensitive gaps; MLPs kept as bounded offline accelerators.
- Amino acid–sugar core remains quantitative trunk.

### S2: Literature Generalization to Family-Aware Execution (Completed)

**Goal**: Every SLR family becomes first-class ingestible unit while preserving benchmark and reporting machinery.

**Implementation closure**:
- `data/lit/chemistry_family_scope_registry.json` and `data/lit/family_ingestion_plan.json` use canonical `family_id` keys.
- All payload surfaces extended with family metadata: `benchmark_intake_registry.json`, `flavor_reference_payloads.json`, `retention_reference_payloads.json`, `computational_priors.json`, `process_gap_registry.json`, `matrix_decision_panel.json`.
- `src/literature_learning_loop.py` generalized; template generation driven by declared payload role.
- Generalized payload vocabulary: `benchmark_payload`, `flavor_reference_payload`, `retention_payload`, `process_state_calibration`, `computational_prior`, `safety_payload`, `structural_gap_entry`.
- `results/validation/literature_learning_loop.{md,json}` and `results/validation/literature_runtime_templates.json` generated with family-aware summaries.

### S3: Runtime Refactoring to Family-Aware Registries (Completed)

**Goal**: Replace family-specific hardcoding with indexed family queries.

**Implementation closure**:
- `src/literature_runtime.py` exposes `query_family_runtime_priors(...)`, `query_flavor_reference_entries(...)`, `query_retention_reference_entries(...)`.
- Legacy getters (`get_pyrazine_control_priors`, `get_furanone_priors`, etc.) remain as thin compatibility wrappers.
- Three runtime concepts separated: reference payloads, retention payloads, family priors.
- `src/matrix_prior_registry.py` extended with `query_family_prior_entries(...)` and `summarize_family_prior_bundle(...)`.
- `build_flavor_axis_summary(...)` now exposes `family_prior_bundle` for scientist-facing diagnostics.

### S4: Decision Panel Expansion (Completed)

**Goal**: Observable panel now covers all family lanes, not just original targets.

**Implementation closure**:
- Decision panel split into 8 explicit family lanes with metadata: `panel_role`, `observable_kind`, `modeling_regimes`.
- Non-volatile state variables (nucleotide enrichment, fermentation precursor loading, thiamine, process severity) now first-class outputs.
- `build_flavor_axis_summary(...)` emits `family_state_markers` with influence modes (`upstream_state_only` vs `upstream_state_plus_marker_panel`).
- `src/presentation.py` surfaces panel metadata in projection tables and diagnostics.

### S5: New Family Lanes Implementation (Completed)

**Goal**: Sequence work to maximize product value quickly.

**Implementation closure** (lanes 1–5):
- **Lane 1 (Lipid oxidation/crosstalk)**: Explicit dual-lane with benchmark targets + retention and crosstalk priors. `maillard_closure_pressure` surfaced. Off-note markers and closure delta wired into projection.
- **Lane 2 (Donor hierarchy)**: Ribose, xylose, glucose, fructose, phosphorylated sugars, specialty donors no longer interchangeable. `src/precursor_resolver.py` updated; donor-family evidence in outputs.
- **Lane 3 (Fermentation pretreatment)**: Bounded pretreatment node upstream of cooking; modifies precursor pools, nucleotide support, thiamine, off-note burden, pH. Uses literature-backed fold-change payloads.
- **Lane 4 (Thiamine/nucleotides/glutathione)**: Support lanes once donor hierarchy and fermentation exist. Upstream availability modifiers + bounded thermal-routing priors.
- **Lane 5 (Caramelization/degradation)**: Severity and failure-mode lane; alternative proteins primarily in matrix-family coverage.
- `src/literature_runtime.py` exposes `build_family_upstream_contract(...)` for donor reweighting, pretreatment pH shifts, support activation.
- `src/pipeline.py` applies upstream contract before pathway enumeration.
- `src/presentation.py` surfaces effective runtime pH, donor routing, pretreatment interventions, upstream-added precursors.

### S6: Family-Lane Validation and Reporting (Completed)

**Goal**: User-facing trust surface shows what each family is doing.

**Implementation closure**:
- `src/benchmark_validation.py` enriched with chemistry-family, SLR-family, payload-role metadata; exposes `build_family_lane_validation_artifact(...)`.
- `src/presentation.py` shows chemistry-family and payload-role columns in benchmark summaries.
- `src/reporting.py` adds `family_evidence_ladder`, `family_runtime_support_summary`, `family_specific_open_gaps`, per-run `family_lane_sensitivity` to JSON and Markdown reports.
- `src/family_lane_sensitivity.py` provides toggle-impact artifact distinct from barrier-offset sensitivity.

### S7: Execution Sequence (Completed)

**Goal**: Record implemented phases A–E for traceability.

**Implementation closure** (all phases complete):
- **Phase A** (Data model/ingestion): Family-aware metadata flows through registries; `src/family_ingestion_plan.py` builds and renders plan; `src/literature_learning_loop.py` summarizes by family/payload-role.
- **Phase B** (Runtime refactor): Family-aware query helpers in `src/literature_runtime.py`; family-aware prior accessor in `src/matrix_prior_registry.py`; legacy wrappers stable.
- **Phase C** (First new family): Lipid oxidation/crosstalk wired; family-level outputs in reporting and recommendation explanations.
- **Phase D** (Upstream enabling): Donor hierarchy, fermentation pretreatment, thiamine/nucleotide/glutathione support all implemented.
- **Phase E** (Coverage/trust): Decision panel expanded, family-aware validation summaries, family open-gap reporting.

### S8: Definition of Done (Completed)

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

#### 5.1 Extrusion Process Model — **[L4 | Blocked on design decision]**
The dominant commercial PBMA process (extrusion) is the weakest regime in the tool. Report 12 provides the calibration data.

- [ ] **5.1a SME as independent process variable**: Add `sme_kj_per_kg` to `ReactionConditions`; compute `T_effective = T_jacket + f(SME)` correction from Report 12 SME–temperature offset data (5–40°C delta).
- [ ] **5.1b Moisture regime classifier**: Add `moisture_regime: Literal['lme', 'hme']` to conditions; flip sign of `∂acrylamide/∂moisture` in `predict_acrylamide()` for LME (<40%) vs HME (>50%).
- [ ] **5.1c Pre-extrusion damage base load**: Add `pre_extrusion_damage_load(protein_type)` returning baseline furosine and LAL from alkaline extraction + spray drying.
- [ ] **5.1d Autoclave sterilization step**: Model as a discrete additive damage increment (121–126°C, 15–30 min) separate from extruder.
- [ ] **5.1e Spatial discretization (pending scope decision)**: Either plug-flow reactor model or sequential isothermal zone model for barrel.

#### 5.2 Protein Source Registry — **[L3]**
Report 06 produced meaty potential multipliers for 14 protein sources. None are encoded at runtime.

- [ ] Create `data/lit/protein_source_registry.json` with AA composition, meaty potential multiplier, off-note penalty, LOX activity flag, and methoxypyrazine ceiling per source.
- [ ] Wire registry into `matrix_correction.py` so switching protein source auto-adjusts correction factors.
- [ ] Add CLI `--protein-source` flag to `run_pipeline.py` that selects from the registry.
- [ ] Include engineering heuristics from Report 06 (e.g., methoxypyrazine non-correctable flag for pea >50%, wheat gluten Cys advantage for MFT).

#### 5.3 Ingest Benchmark-Eligible Deep Research Data — **[L2]**
~20 benchmark-eligible datasets exist in the Gemini Deep Research corpus but aren't wired into the runtime.

- [ ] **SPI-HVP + xylose**: Ingest MFT OAV 450, FFT OAV 84 (120°C, 30 min, pH 6.0) as a matrix benchmark.
- [ ] **Wheat gluten HVP + xylose**: Ingest MFT OAV 850 as the highest-MFT plant-source reference.
- [ ] **Acrylamide fast kinetics**: Calibrate `predict_acrylamide()` against 22.36→62.62 µg/kg in 20–30s at 130°C.
- [ ] **CML/CEL commercial PBA ranges**: Validate `predict_cml()` and `predict_cel()` against 20+ product dataset (Foods 2023).
- [ ] **Furosine formation-elimination crossover**: Validate non-monotonic behavior (peak ~8.7 mg/100g at 140°C, fall above 150°C).

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

- [ ] Add sugar reactivity multipliers to `barrier_constants.py` keyed on donor identity.
- [ ] Wire donor identity into `ode_kinetics.py` rate constant computation.
- [ ] Add benchmark test comparing ribose vs glucose MFT yield at equivalent conditions.

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
- [x] P4 reaction-centered benchmarking now materializes data/lit/reaction_benchmark_set.json plus results/validation/p4_reaction_benchmark.{md,json}.
- [x] P4 geometry benchmarking now materializes data/lit/p4_geometry_benchmark_set.json plus results/validation/p4_geometry_benchmark.{md,json} and results/validation/p4_geometry_assessment.{md,json}.
- [x] P4 MLP governance now materializes data/lit/mlp_candidate_registry.json plus results/validation/p4_mlp_assessment.{md,json} and results/validation/p4_adoption_notes.{md,json}, with ResultsDB persistence for adoption decisions.
- [x] P4 now also materializes data/lit/mlp_external_benchmark_evidence.json plus results/validation/p4_external_mlp_landscape.{md,json}, making external community evidence explicit as shortlist priority rather than as a substitute for local Maillard validation.
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
