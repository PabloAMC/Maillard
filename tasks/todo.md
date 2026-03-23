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

## Current Sprint

This sprint starts from a more honest state description:

- amino acid plus sugar Maillard chemistry is the current first-class core because it is the most transferable and benchmarkable trunk of the product, not because the rest of the chemistry is minor
- lipid oxidation and lipid-derived carbonyl crosstalk are already partially encoded, but not yet promoted to an equally explicit family-level validation and ingestion surface
- MLPs and DFT are support tooling for benchmark-visible mechanistic gaps, not the product roadmap themselves

Archive rule for this file:

- leave completed program foundations below as a compact record
- keep the active sprint clean by only listing unfinished or newly reopened work near the top

### Sprint S9. Skipped Test Triage and QM Optionality

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

### Sprint S0. Script Code Quality Improvements

Goal: Refactor main orchestration scripts (`scripts/`) to improve readability, modularity, and maintainability without altering scientific logic.

S0.1: Refactor `run_pipeline.py`
- [x] Extract argument parsing into a `parse_args()` function.
- [x] Separate the 250-line `main()` into `run_inverse_design()` and `run_forward_pipeline()`.
- [x] Consolidate shared validation/reporting logic into a helper function to DRY the code.
- [x] Add missing type hints to function signatures.

S0.2: Refactor `run_cantera_kinetics.py`
- [x] Break down the monolithic `run_simulation()` function into discrete steps: `build_mechanism()`, `run_integration()`, `plot_results()`.
- [x] Hoist inline hardcoded dictionaries (`NAME_MAP`, `lookup`) to module-level constants to prevent redeclaration inside loops.

S0.3: Refactor `calibrate_barriers.py`
- [x] Extract the target prediction evaluation logic from `run_sim()` into a standalone `compute_mae_for_benchmark()` helper.
- [x] Hoist `NAME_TO_SMILES` to a top-level constant.

S0 measurable exit criteria:
- [x] The three main scripts are modularized with clear function boundaries, docstrings, and type hints.
- [x] They execute identically to their previous monolithic versions.

Implemented closure:

- The 2026-03-23 audit found no remaining unchecked S0/S0b/S0c/S0d implementation work in the active codebase; the remaining quality debt in this file is markdown formatting, not unfinished script or library refactors.

### Sprint S0b. Script Robustness Improvements

Goal: Enhance reliability, validate inputs, and handle failures gracefully across orchestration scripts.

S0b.1: Foundation setup
- [x] Create `src/logger.py` for structured logging.
- [x] Create `src/config.py` for centralized defaults.
- [x] Create `src/exceptions.py` for specific `MaillardError` definitions.

S0b.2: Integration into CLI Scripts
- [x] Update `run_pipeline.py` with logger, config, and typed exceptions.
- [x] Update `run_cantera_kinetics.py` similarly.
- [x] Update `calibrate_barriers.py` similarly.

S0b.3: CLI Testing & Parallelisation
- [x] Add parallel execution (`ProcessPoolExecutor`) to `calibrate_barriers.py`.
- [x] Add integration smoke tests in `tests/unit/test_cli_scripts.py`.

S0b measurable exit criteria:
- [x] Scripts use logging instead of stdout prints (except for UI outputs).
- [x] Scripts no longer harbor internal hardcoded magic default constants.
- [x] Tests pass via `docker_maillard.sh pytest`.

### Sprint S0c: Deep Library Code Quality Refactor
Goal: Clean up fatal anti-patterns (sys.exit, mid-file imports, raw prints) migrating from CLI to the native `src/` libraries.

- [x] Remove `sys.exit(1)` from `src/recommend.py` and raise typed exceptions.
- [x] Migrate scattered `print()` logs to `src.logger` inside `src/recommend.py`, `src/kinetics.py`, `src/headspace.py`, `src/mlp_optimizer.py`, and `src/diffusion_ts.py`.
- [x] Remove duplicated definitions (e.g., `_weight`) and hoist mid-file PEP8 import violations in `src/recommend.py`.

### Sprint S0d: Code Audit Cleanups
Goal: Address remaining structural and quality issues surfaced during the deep audit.

- [x] Create `src/text_utils.py` and consolidate identical `_normalize_name` functions from `benchmark_validation.py` and `matrix_targets.py`.
- [x] Create `src/chem_utils.py` and consolidate `_canon` / `_canonical` from `recommend.py` and `smirks_engine.py`.
- [x] Extract repetitive step collection logic in `smirks_engine.py` into a helper function.
- [x] Replace remaining `print()` statements in `src/dft_refiner.py` with `logger` calls.
- [x] Remove dead interactive `__main__` stubs from library modules (`smirks_engine`, `headspace`, `kinetics`, `cantera_export`, `sensory`).
- [x] Deprecate/remove legacy static `predict()` flow from `src/recommend.py`.

### Sprint S1. Expand Beyond Amino Acid-Sugar Without Losing The Core

Goal: make the repo state of the art at the product level by explicitly promoting the next chemistry families that change real plant-based flavor decisions.

S1.1: Decide the next chemistry family lane

- [x] Treat lipid oxidation and carbonylic crosstalk as the default next expansion unless a stronger benchmark-visible family appears.
- [x] Keep amino acid plus sugar, Strecker, and sulfur chemistry as the benchmarked foundation rather than diffusing effort across many under-ingested families.
- [x] Name explicitly which additional families are open gaps versus bounded lanes versus first-class core.

Implemented closure:

- `src/family_strategy_policy.py` now derives a machine-readable strategy artifact from the chemistry-family scope and ingestion plan, fixing lipid oxidation and carbonylic crosstalk as the default next expansion lane.
- The same artifact makes the quantitative trunk explicit and classifies families into first-class core, high-priority partial lanes, bounded lanes, and open gaps.

S1.2: Generalize SLR ingestion beyond amino acid-sugar

- [x] Use the same ingestion contract already proven in the repo: intake registry for candidate papers, runtime payloads for closable evidence, and process-gap registries for non-closable scope.
- [x] Add a chemistry-family scope artifact that states, for each family, the preferred runtime payload type: benchmark payload, flavor reference payload, retention payload, computational prior, safety payload, or structural-gap registry.
- [x] Avoid creating a separate narrative-only markdown workflow for new families; new families should enter the same machine-readable ingestion path used by the current SLR loop.
- [x] For lipid oxidation and crosstalk specifically, structure ingestion as a dual lane: observable benchmark targets plus retention and competition payloads.

Implemented closure:

- The existing `chemistry_family_scope_registry`, `family_ingestion_plan`, and payload metadata remain the shared ingestion contract, and `family_strategy_policy` now states that contract explicitly as machine-readable policy.
- The new strategy artifact encodes a lipid dual-lane policy: observable benchmark-plus-retention payloads on one side and crosstalk/competition priors plus structural gaps on the other.

S1.3: Keep DFT and MLP policy clean for the next sprint

- [x] Archive the now-completed P4 setup work as infrastructure, not as the active scientific objective.
- [x] State clearly that selective DFT is reserved for benchmark-visible sulfur, carbonyl-transfer, and TS-sensitive gaps after cheap-first screening.
- [x] State clearly that MLPs remain bounded offline accelerators until local geometry or TS benchmarks accept them on Maillard-relevant systems.
- [x] Prevent roadmap language from implying that better MLP branding can substitute for missing family-level ingestion or benchmark closure.

Implemented closure:

- `family_strategy_policy` now captures P4 as infrastructure, keeps DFT as selective cheap-first escalation tooling, and keeps MLPs bounded as offline accelerators rather than a substitute for family-lane closure.

S1 measurable exit criteria:

- [x] The repo can explain why amino acid-sugar was first without implying that other families are unimportant.
- [x] The repo can identify lipid oxidation and crosstalk as either the active next family sprint or an explicitly rejected option with reasons.
- [x] The repo has a machine-readable chemistry-family scope artifact linked in reporting.
- [x] The top of this file reads as a clean active sprint rather than a mixed archive of finished and unfinished eras.

### Planning Notes. Multi-Family Runtime Extension

Goal: extend the existing amino-acid-sugar literature-to-runtime stack to the new 01-10 family corpus without creating a second parallel architecture.

Working conclusion from code review:

- the repo already has the right high-level ingestion contract:
   - intake registry for candidate references
   - executable benchmark payloads where closure exists
   - reference or retention payloads where runtime use is bounded
   - computational priors for directional or calibration-grade knowledge
   - process-gap registry for structurally non-closable scope
- the current implementation is still heavily shaped around the original family:
   - `src/literature_learning_loop.py` knows only a small set of template kinds and artifact classes
   - `src/literature_runtime.py` contains family-specific hardcoded getters and routing for pyrazines, furanones, thiamine, Strecker crosstalk, and a narrow set of retention rules
   - `src/matrix_prior_registry.py` is indexed by `protein_type` only and cannot yet represent family-specific priors cleanly
   - `src/benchmark_validation.py` and the benchmark metadata model are not yet explicitly family-lane aware
   - `data/lit/matrix_decision_panel.json` covers the original decision panel plus adverse markers, but not the broader family-by-family observable panel implied by 01-10
   - `src/family_sensitivity.py` screens kinetic reaction families, not literature families or ingestion lanes
- the repo therefore does not need a new architecture; it needs a generalized family-aware layer on top of the existing contracts.

Non-negotiable extension rules:

- do not parse the markdown SLR files directly at runtime
- do not create a narrative-only workflow for families 02-10
- keep the amino-acid-sugar core as the quantitative trunk
- add new families through machine-readable companion artifacts tied back to the numbered SLRs
- keep each new family explicit about whether it is:
   - benchmark-closable
   - calibration-grade
   - directional prior only
   - safety lane
   - structural gap only

### S2. Generalize The Literature Contract To Family-Aware Execution

Goal: make the repo capable of representing every SLR family as a first-class ingestible unit while preserving the existing benchmark and reporting machinery.

S2.1: Freeze the canonical family-to-runtime map

- [x] Create a machine-readable family ingestion plan artifact that maps each numbered SLR family 01-10 to:
   - chemistry family id
   - source SLR file
   - preferred payload types
   - target runtime modules
   - target compounds or state variables
   - benchmarkability status
   - next curation action
- [x] Use the current numbered corpus in `data/Gemini_Deep_Research/01.md` through `10.md` as the scientific source and the preserved synthetic summary in `data/Gemini_Deep_Research/cross_family_promotion_and_ingestion_priorities.md` as the ranking layer.
- [x] Decide one canonical identifier scheme for families so the same key appears in:
   - `data/lit/chemistry_family_scope_registry.json`
   - intake and payload artifacts
   - reporting summaries
   - future validation outputs
- [x] Keep matrix family and chemistry family as separate axes:
   - chemistry family answers what reaction knowledge is encoded
   - matrix family answers where that knowledge is supported

Implemented closure:

- `data/lit/chemistry_family_scope_registry.json` now uses the same canonical `family_id` keys as `data/lit/family_ingestion_plan.json`.
- `results/validation/family_identifier_contract.{json,md}` verifies scope-plan alignment, payload-family validity, and chemistry-vs-matrix axis separation.

S2.2: Add family metadata to all literature payload surfaces

- [x] Extend `data/lit/benchmark_intake_registry.json` entries to carry at least:
   - `chemistry_family`
   - `slr_family_source`
   - `payload_role`
   - `observable_panel_tags`
   - `process_state_scope`
- [x] Extend `data/lit/flavor_reference_payloads.json` entries with family metadata so sulfur, nucleotide, caramelization, and fermentation-derived references can coexist without ad hoc code paths.
- [x] Extend `data/lit/retention_reference_payloads.json` so retention and release rules can be grouped by chemistry family, not only by matrix and compound.
- [x] Extend `data/lit/computational_priors.json` so priors can be indexed by both `protein_type` and `chemistry_family` where appropriate.
- [x] Extend `data/lit/process_gap_registry.json` so structural blockers are attributable to a family lane and not just to a benchmark id.
- [x] Extend `data/lit/matrix_decision_panel.json` with family tags so the decision panel can say which compounds are part of which lane and whether they are core targets, adverse markers, or process severity markers.

S2.3: Generalize the ingestion compiler instead of multiplying hardcoded templates

- [x] Refactor `src/literature_learning_loop.py` so template generation is driven by declared payload role rather than a short hardcoded map in `READY_TEMPLATE_KIND`.
- [x] Add support for family-specific output targets beyond the current set of:
   - benchmark
   - process state calibration
   - directional prior
   - safety reference
- [x] Introduce a generalized template vocabulary such as:
   - `benchmark_payload`
   - `flavor_reference_payload`
   - `retention_payload`
   - `process_state_calibration`
   - `computational_prior`
   - `safety_payload`
   - `structural_gap_entry`
- [x] Generate a new review artifact showing the queue by chemistry family and payload type, not only by encoding status.

Implemented closure:

- `src/literature_learning_loop.py` now separates `source_payload_role` from `target_payload_types` and derives the primary `template_kind` from declared runtime artifact payload types, with readiness status used only as a fallback.
- `results/validation/literature_learning_loop.{md,json}` now include `payload_queue_review` and summary counts by generalized payload type.
- `results/validation/literature_runtime_templates.json` now uses the generalized payload vocabulary (`benchmark_payload`, `computational_prior`, `process_state_calibration`, `safety_payload`) and carries family-aware metadata needed for future `flavor_reference_payload` and `retention_payload` templates.

### S3. Refactor Runtime Access From Family-Specific Hardcoding To Family-Aware Registries

Goal: stop baking the original family assumptions directly into `src/literature_runtime.py` and related modules.

S3.1: Replace singleton getters with indexed family queries

- [x] Refactor `src/literature_runtime.py` so current functions like `get_pyrazine_control_priors`, `get_furanone_priors`, `get_thiamine_priors`, and `get_strecker_crosstalk_priors` become family-aware query helpers rather than one-off retrieval functions.
- [x] Add generic lookup helpers that can answer:
   - which priors apply to a compound, family, matrix family, and process state
   - which retention references apply to a given compound and family lane
   - which flavor references are scoring targets versus diagnostics versus constraints
- [x] Keep the legacy helpers as compatibility wrappers until the new query layer is stable.

Implemented closure:

- `src/literature_runtime.py` now exposes `query_family_runtime_priors(...)`, `query_flavor_reference_entries(...)`, and `query_retention_reference_entries(...)` as the query-first family-aware interface.
- Legacy getters (`get_pyrazine_control_priors`, `get_furanone_priors`, `get_thiamine_priors`, `get_strecker_crosstalk_priors`) remain as thin compatibility wrappers over the new prior query layer.

S3.2: Separate three runtime concepts that are currently partially conflated

- [x] Keep `reference payloads` for scoring or reporting targets
- [x] Keep `retention payloads` for observable release or attenuation behavior
- [x] Keep `family priors` for mechanistic directional modifiers and pathway plausibility
- [x] Ensure that no single payload file silently becomes all three.

Implemented closure:

- `build_flavor_axis_summary(...)` continues to consume compatibility getters, but those getters now route through the dedicated prior query path rather than section-specific hardcoding.
- flavor target lookup and policy summary now flow through `query_flavor_reference_entries(...)`.
- retention routing helpers now flow through `query_retention_reference_entries(...)`, preserving the current runtime behaviour while separating retention evidence from flavor references and priors.

S3.3: Make priors family-aware without losing the existing protein-type bundle interface

- [x] Preserve `src/matrix_prior_registry.py` for matrix-state summaries.
- [x] Add a companion family-aware registry or accessor so the runtime can answer questions such as:
   - what sulfur-family priors apply in soy isolate
   - what fermentation-family priors modify the precursor pool before cooking
   - what off-note-family priors trap dicarbonyls or block amino groups
- [x] Do not overload the current matrix prior bundle with every chemistry concept; add a separate layer if necessary.

Implemented closure:

- `src/matrix_prior_registry.py` now keeps the existing matrix-state bundle intact while adding `query_family_prior_entries(...)` and `summarize_family_prior_bundle(...)` for chemistry-family-aware prior access.
- `build_flavor_axis_summary(...)` now exposes `family_prior_bundle` so scientist-facing diagnostics can distinguish matrix-state support from chemistry-family priors.

### S4. Extend The Scientist Decision Panel Beyond The Original Family

Goal: make the tool useful for real formulation decisions across the full family map rather than only the original target panel.

S4.1: Freeze the expanded observable panel by lane

- [x] Split the decision panel into explicit lane groups:
   - sulfur positives
   - Strecker aldehydes
   - pyrazines
   - furanones and caramelization severity markers
   - lipid oxidation adverse markers
   - nucleotide and umami support markers
   - fermentation pretreatment state markers
   - safety markers
- [x] For each panel compound or state variable, encode:
   - family ownership
   - evidence state
   - whether it is scored, constrained, diagnostic, or report-only
   - which modeling regimes it applies to

Implemented closure:

- `data/lit/matrix_decision_panel.json` now carries explicit `panel_role`, `observable_kind`, and `modeling_regimes` metadata by lane, plus first-class pretreatment/state-marker entries.
- `src/matrix_targets.py`, `src/projection_metadata.py`, and `src/projection_utils.py` now preserve and propagate those panel fields into downstream reporting.

S4.2: Keep family-specific observables from being overclaimed

- [x] Do not force every family into the same standard as MFT or FFT.
- [x] Allow non-volatile but decision-relevant state variables to be first-class outputs where appropriate, especially for:
   - nucleotide enrichment
   - fermentation-derived precursor loading
   - thiamine availability
   - process severity markers such as HMF or furfural
- [x] Make the reports explicit when a family influences the recommendation through upstream state changes rather than through direct volatile prediction.

Implemented closure:

- `build_flavor_axis_summary(...)` now emits `family_state_markers` with explicit influence modes (`upstream_state_only` vs `upstream_state_plus_marker_panel`) for thiamine, nucleotide, pretreatment, and caramelization state support.
- `src/presentation.py` now surfaces panel role/kind/regime metadata in projection tables and prints state-marker diagnostics directly in the flavor-axis section.

### S5. Implement New Family Lanes In The Order That Maximizes Product Value

Goal: sequence the work so the next code changes improve recommendation quality quickly instead of scattering effort.

S5.1: Priority lane 1. Lipid oxidation and carbonylic crosstalk

- [x] Promote `lipid_oxidation_and_carbonylic_crosstalk` to the first new family lane because it has the highest immediate impact on PBMA decisions.
- [x] Split this lane into two runtime sub-lanes:
   - adverse marker generation and retention
   - carbonyl competition and crosstalk priors
- [x] Extend `src/lipid_oxidation.py`, `src/literature_runtime.py`, `src/projection.py`, and reporting so lipid-derived aldehydes are not only adverse outputs but also modifiers of Maillard closure.
- [x] Define benchmark-ready targets from the existing off-note panel and calibration-ready crosstalk priors from the literature.

Implemented closure:

- `src/lipid_oxidation.py` now resolves generic oil inputs into proxy lipid loads and exposes named benchmark-ready off-note markers instead of only raw SMILES outputs.
- `src/literature_runtime.py` now treats family 02 as an explicit dual lane: adverse marker generation plus retention, and carbonyl competition plus crosstalk priors, with `maillard_closure_pressure` surfaced alongside benchmark-ready targets and Lincoln 2025 prior ids.
- The pipeline projection path now lets family 02 modify target closure explicitly through `maillard_closure_delta` instead of only inflating off-note risk.
- `src/presentation.py` now reports lipid benchmark targets, crosstalk priors, and closure pressure directly in the flavor-axis diagnostics.

S5.2: Priority lane 2. Carbonyl donor hierarchy

- [x] Extend precursor handling so ribose, xylose, glucose, fructose, phosphorylated sugars, and specialty donors are no longer treated as interchangeable sugar inputs.
- [x] Update `src/precursor_resolver.py`, `src/recommend.py`, and runtime priors to encode donor identity as a first-class variable.
- [x] Add donor-family evidence to recommendation outputs so the tool can explain why a formulation is donor-limited rather than just sugar-limited.

S5.3: Priority lane 3. Fermentation pretreatment

- [x] Add a new `fermentation_pretreatment_node` concept that sits upstream of cooking chemistry and can modify:
   - precursor pools
   - nucleotide support
   - thiamine availability
   - off-note burden
   - pH routing
- [x] Implement this as a bounded pretreatment layer first, not as a full fermentation simulator.
- [x] Start by using literature-backed fold-change or enrichment payloads rather than detailed microbial kinetics.

S5.4: Priority lane 4. Thiamine, nucleotides, and glutathione or peptide support

- [x] Treat thiamine, nucleotide, and glutathione families as explicit additive or support lanes once donor hierarchy and fermentation pretreatment exist.
- [x] Represent them as upstream availability modifiers plus bounded thermal-routing priors, not as fully independent kinetic solvers on day one.
- [x] Reuse existing sulfur and Strecker scoring surfaces wherever possible.

S5.5: Priority lane 5. Carbohydrate pyrolysis and alternative matrix scope

- [x] Keep caramelization and carbohydrate thermal degradation as a severity and failure-mode lane before promoting it to a major optimization axis.
- [x] Keep alternative proteins primarily in `matrix_family_coverage` and bounded priors until enough direct quantitative closure exists to justify family-specific benchmarks.

Implemented closure:

- `src/literature_runtime.py` now exposes a reusable `build_family_upstream_contract(...)` layer that reweights donor pools, applies bounded pretreatment pH shifts, carries support-lane activation, and injects bounded thiamine support when formulation metadata says it is available but not explicitly listed.
- `src/pipeline.py` now applies that upstream contract before pathway enumeration so donor hierarchy and fermentation pretreatment modify effective precursor loading and runtime pH instead of remaining only diagnostic lanes.
- `src/presentation.py` now surfaces the effective runtime pH, donor-class routing, donor pool factors, pretreatment interventions, and upstream-added precursors in the scientist-facing flavor-axis report.
- Focused runtime, pipeline, usability, pre-processor, and integration subsets passed in Docker after the change.

### S6. Make Validation And Reporting Family-Lane Aware

Goal: the user-facing trust surface must show what each family is doing, not just what the free core predicts.

S6.1: Extend benchmark metadata and validation outputs

- [x] Add chemistry-family and payload-role metadata to benchmark and calibration artifacts.
- [x] Update `src/benchmark_validation.py` so benchmark summaries can be grouped by chemistry family and by lane.
- [x] Add family-lane summaries to the validation outputs analogous to the current chemistry-family scope and matrix-family coverage artifacts.

S6.2: Extend reporting surfaces

- [x] Update `src/reporting.py` to expose:
   - family ingestion plan artifact
   - family runtime support summary
   - family-specific open gaps
   - family-specific evidence ladders in recommendation outputs
- [x] Make it obvious in scientist-facing reports whether a recommendation was driven by:
   - core benchmarked chemistry
   - calibration-grade family payloads
   - directional priors
   - structural-gap extrapolation

S6.3: Keep family sensitivity separate from barrier sensitivity

- [x] Leave `src/family_sensitivity.py` focused on kinetic reaction-family perturbation.
- [x] Add a new family-lane sensitivity artifact if needed that asks a different question:
   - how much does enabling or disabling a literature family lane change the recommendation or ranking outcome
- [x] Do not overload the existing barrier-offset sensitivity tool with literature-family semantics.

Implemented closure:

- `src/benchmark_validation.py` now enriches `BenchmarkSummary` with chemistry-family, SLR-family, and payload-role metadata and exposes `build_family_lane_validation_artifact(...)` plus a Markdown renderer and generator script.
- `src/presentation.py` now shows chemistry-family and payload-role columns in the benchmark summary surface so family-aware validation is visible instead of buried in JSON.
- `src/reporting.py` now adds `family_evidence_ladder`, `family_runtime_support_summary`, `family_specific_open_gaps`, and per-run `family_lane_sensitivity` to the JSON and Markdown reports.
- `src/family_lane_sensitivity.py` now provides a runtime toggle-impact artifact distinct from `src/family_sensitivity.py`, keeping literature family-lane semantics separate from barrier-offset perturbations.

### S7. Implemented Execution Sequence Record

This sequence is now historical record. The originally proposed S7 work was absorbed into S2, S3, S5, and S6 and is no longer an open sprint.

Phase A. Data model and ingestion generalization

- [x] Add chemistry-family metadata to all literature payload files.
- [x] Build the family ingestion plan artifact and generator.
- [x] Extend the literature learning loop to summarize by family and payload role.

Implemented closure:

- Family-aware metadata now flows across the literature payload registries consumed by `src/literature_family_registry.py` and `src/literature_learning_loop.py`.
- `src/family_ingestion_plan.py` builds and renders the family ingestion plan artifact, and `scripts/generate_family_ingestion_plan.py` now emits the validation markdown/json pair.
- `src/family_ingestion_plan.py` builds and renders the family ingestion plan artifact, and `scripts/generators/generate_family_ingestion_plan.py` now emits the validation markdown/json pair.
- `src/literature_learning_loop.py` summarizes ready references by chemistry family, payload role, and generalized payload type.

Phase B. Runtime refactor

- [x] Add family-aware query helpers in `src/literature_runtime.py`.
- [x] Introduce a family-aware prior accessor alongside `src/matrix_prior_registry.py`.
- [x] Keep old hardcoded accessors as compatibility wrappers until tests are migrated.

Implemented closure:

- `src/literature_runtime.py` now exposes query-first helpers for priors, flavor references, and retention references.
- `src/matrix_prior_registry.py` now supports family-aware prior queries and bundle summaries in parallel with the original protein-centric accessors.
- Legacy one-off getters remain as thin compatibility wrappers over the family-aware layer.

Phase C. First new family implementation

- [x] Implement lipid oxidation and crosstalk as the first explicit multi-payload family lane.
- [x] Wire family-level outputs into reporting and recommendation explanations.

Implemented closure:

- `src/literature_runtime.py` and `src/pipeline.py` now treat lipid oxidation and carbonylic crosstalk as an explicit runtime lane with causal upstream and scoring effects.
- Scientist-facing explanation surfaces now expose active family lanes, lane-specific summaries, and family-level score deltas.

Phase D. Upstream enabling families

- [x] Implement donor hierarchy.
- [x] Implement fermentation pretreatment.
- [x] Add thiamine, nucleotide, and glutathione support lanes.

Phase E. Coverage and trust surfaces

- [x] Expand the decision panel.
- [x] Add family-aware validation summaries.
- [x] Add family-lane open-gap reporting.

Implemented closure:

- The decision-panel layer now carries family-aware metadata through `data/lit/matrix_decision_panel.json`, `src/matrix_targets.py`, and projection metadata.
- `src/benchmark_validation.py` and `src/presentation.py` expose family-aware validation summaries, including chemistry families, payload roles, and family-lane validation artifacts.
- `src/reporting.py` now emits family evidence ladders, runtime support summaries, and family-specific open-gap reporting for each run and comparison report.

### S8. Definition Of Done For This Extension Program

- [x] Every numbered SLR family 01-10 maps to one explicit machine-readable ingestion lane.
- [x] No family is represented only by markdown if it is supposed to affect runtime behavior.
- [x] Reports can explain which family lanes were active in a recommendation and with what evidence strength.
- [x] The amino-acid-sugar core remains the quantitative trunk rather than being diluted into a generic all-families model.
- [x] Lipid oxidation and crosstalk, donor hierarchy, and fermentation pretreatment are all visible as explicit runtime concepts.
- [x] Structural gaps remain explicit where no benchmark-grade closure exists.

Implemented closure:

- `data/lit/family_ingestion_plan.json`, `src/literature_family_registry.py`, and family-aware payload metadata provide explicit machine-readable lanes for SLR families 01-10.
- Runtime behavior no longer depends on markdown-only SLR narratives; family behavior is wired through structured registries, priors, and runtime query helpers.
- `src/reporting.py` now exposes active family lanes, evidence posture, open gaps, and lane sensitivity so recommendations disclose evidence strength by family.
- The quantitative trunk remains anchored on amino-acid plus sugar chemistry while additional families are layered as bounded lanes with explicit uncertainty posture.
- Structural gaps remain first-class via `data/lit/process_gap_registry.json`, family-specific open-gap reporting, and non-promotional lane policy.

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

## Completed Foundations (Collapsed)

- [x] Docker plus conda is the authoritative validation path.
- [x] Free-precursor strict-ready benchmark surface exists and is versioned.
- [x] Matrix-only pea/soy intake-headspace path is executable and benchmarked directionally.
- [x] Validation overview now shows all quantitative benchmarks, including matrix-only and matrix-augmented lanes.
- [x] Reporting exposes provenance, evidence ladders, confidence warnings, MFT/furfural ratio, meaty quality penalty, and sulfur trapping.
- [x] Matrix observable calibration is applied exactly once in the matrix-only path.
- [x] Validation contract now includes both max ratio and mean absolute log-scale error.
- [x] Per-benchmark scale thresholds are supported for the authoritative free-precursor set.
- [x] xTB, selective DFT, and offline ML-potential integration hooks exist as refinement tooling.

## Strategic Answer

### What problem should we solve?

- [ ] Make the tool reliably useful for deciding what to test next in plant-based meat flavor design, not merely for reproducing a narrow benchmark set.

### Does this reduce to modeling three classes correctly?

- [ ] Yes, but as trust regimes rather than isolated modules: free precursors, pea/soy matrices, and extrusion-heavy systems must form a coherent ladder of confidence.

### Are ML potentials the elegant main solution?

- [ ] No. The elegant main solution is benchmark-driven observable and accessibility modeling first, selective mechanistic refinement second, and external offline ML-potential acceleration only third.

## Program Record

The sections below are the program record and archived foundations for the current branch. Keep them for traceability, but prefer the Current Sprint section above for active execution.

### Current Operating Constraint

Assume no internal wet-lab loop is available unless explicitly stated otherwise.

That changes the execution order:

- use the literature and machine-readable intake registry to close every benchmark or calibration surface that is actually closable without new experiments
- represent structural literature gaps explicitly instead of pretending they can be calibrated away
- do not count internally constructed mixed matrix benchmarks as external closure; if a mixed matrix surface lacks external measurements, route it into mechanistic triage rather than P0 promotion
- use offline xTB, selective DFT, and external ML potentials only for narrow mechanistic gaps that remain decision-relevant after the literature closure work

### Immediate Execution Order Under A Literature-Only Constraint

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
