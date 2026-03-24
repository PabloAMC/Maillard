# Maillard Strategic Roadmap

## Mission

Build the most useful computational tool for scientists who want to imitate meat-like Maillard chemistry in plant-based systems.

The key product question is:

Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

## Product Thesis

This is not mainly a problem of making one chemistry engine more exact.

It is a problem of combining:

- a quantitatively credible free-precursor core
- matrix-aware observability and accessibility
- process-aware confidence boundaries
- scientist-facing reporting that states what is benchmarked and what is extrapolated

## Product Status Today

### What We Have

- **Family-aware runtime**: Amino acid–sugar core plus 9 additional chemistry families (lipid oxidation, fermentation pretreatment, donor hierarchy, thiamine support, nucleotides, caramelization, sulfur, off-notes, alternative matrices) all with machine-readable ingestion lanes.
- **Family evidence ladder**: Reports now show per-family evidence posture (core benchmarked, calibration-grade, directional prior, or structural-gap extrapolation).
- **Family-lane calibration**: Each family can have its own observable projection factors, retention models, and prior bundles; calibration is explicit per family, not monolithic.
- **Scientist-facing family transparency**: Reports disclose active family lanes, per-lane evidence strength, and open gaps so users know which parts of a recommendation are benchmark-backed vs. extrapolated.
- **Deterministic and optional QM support infrastructure**: Quasi-harmonic helpers, DFT authority lanes, and MLP acceleration lanes all have governance policies and dispatch contracts.

### What Guides This Work

The system combines three regimes of confidence (free precursors, pea/soy matrices, extrusion-heavy systems) with family-aware ingestion and runtime, so a scientist can:
- rank formulations before the wet lab,
- understand which predictions are benchmark-anchored and which are transferred or extrapolated,
- see how each chemistry family contributes to a recommendation,
- trust the confidence boundaries that are stated explicitly.

## Current Sprint

### Sprint S10. Family Predictive Closure and Deviation Reduction (Primary Active)

Goal: turn family-aware reporting into family-aware predictive closure by reducing large residuals, promoting high-impact families from directional to benchmark-linked, and exposing a scientist-operable workflow for alternative proteins.

Observed trigger for S10:

- `results/validation/family_validation_overview.{md,json,png}` shows that only a subset of families has quantitative parity and some family-level max-ratio tails remain very large.
- We must improve real predictive closure, not only metadata completeness.

Track A (highest priority): Quantitative deviation triage and closure on current benchmark-backed families

S10.A1: Build outlier diagnostics artifact from family validation points
- [x] Generate a machine-readable outlier audit listing worst compound-level ratio and log-error points by family, benchmark, and execution path.
- [x] Add robust metrics per family (`median`, `p90`, `max`, trimmed log-error) so one extreme point does not hide in aggregate means.
- [x] Publish a scientist-readable markdown artifact that explains which points dominate family-level error.

S10.A2: Convert diagnostics into benchmark-fix actions
- [x] For each family with large residual tails, classify root cause as one of: mapping mismatch, calibration mismatch, observable projection mismatch, benchmark data ambiguity.
- [x] Create one concrete fix ticket per root cause class and link expected impact metric (for example max-ratio drop in affected family).
- [x] Re-run validation artifacts after each fix and keep only changes that improve benchmark-visible error without regressions.

S10.A3: Promotion gate for family quantitative trust
- [ ] Define family-level promotion rule: a family can be called near-quantitative only if it has minimum benchmark count plus bounded robust error.
- [ ] Expose promotion state in validation markdown/json and README trust language.

Track B (second priority): Literature-to-runtime closure for families with zero quantitative parity

S10.B1: Literature closure queue by family
- [ ] Rank families 03, 04, 05, 06, 07, 10 by expected decision impact and literature closability.
- [ ] For each family, declare minimum viable payload needed for first benchmark-linked closure (benchmark payload, calibration payload, or bounded prior with explicit uncertainty).
- [ ] Prevent narrative-only additions: every promoted literature item must land in machine-readable runtime payloads.

S10.B2: Alternative protein usability
- [ ] Add explicit scientist workflow for alternative protein research (matrix choice, process state, interpretation gate, confidence posture).
- [ ] Expand matrix-family coverage artifacts so users can see direct support vs transferred support by matrix family.

Track C (third priority): Optional QM/ML governance and skipped-test cleanup

- [ ] Keep S9 execution as supporting infrastructure lane.
- [ ] Do not treat S9 completion as sufficient for family predictive closure.
- [ ] Use QM/ML lanes only where Track A diagnostics prove benchmark-visible impact.

S10 measurable exit criteria:
- [x] Every benchmark-backed family has an outlier audit with top residual points and root-cause labels.
- [x] At least one high-residual family shows a verified reduction in robust error metrics after targeted fixes.
- [ ] At least one currently non-quantitative family is promoted to benchmark-linked or calibration-grade with explicit uncertainty bounds.
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

Phase P0.3: Close matrix observability for decision-driving compounds

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
