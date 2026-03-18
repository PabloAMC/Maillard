# Maillard Strategic Backlog

## Product Goal

Convert Maillard from a scientifically promising validation harness into a decision-grade formulation tool for alternative-protein scientists.

For this repository to be genuinely state of the art, it must simultaneously provide:

- trustworthy quantitative envelopes, not only plausible rankings;
- scientist-facing outputs that explain what to do next in formulation terms;
- reproducible benchmark, uncertainty, and domain-of-validity surfaces;
- a workflow that is faster and more informative than spreadsheet-plus-wet-lab iteration.

## Truthful Current Status — 2026-03-17

### Documentation and onboarding review — 2026-03-18

- [x] Rework the entry documentation so a new scientist can understand what the tool does, how to run it, and what not to claim from it.
- [x] Add a dedicated trust-and-limitations guide and a fast quickstart path.
- [x] Add a command reference and a project-structure reference.
- [x] Add a reproducible validation-figure artifact so reliability and limitations can be shown graphically.
- [x] Group development and research notes behind explicit documentation indexes instead of leaving them as an unstructured reading path.

### Validation communication cleanup — 2026-03-18

- [x] Remove the fold-error-by-compound panel from benchmark comparison figures because it is visually confusing for first-pass scientific review.
- [x] Align the validation deep dive with the current proof model: strict quantitative free-precursor evidence vs. non-authoritative matrix readiness.
- [x] Rewrite the root README so a first-time scientist can understand trust levels, main workflows, and how to get useful results without overclaiming.
- [x] Retire legacy validation plot artifacts that no longer represent the current benchmark narrative.

### External matrix benchmark candidate curation — 2026-03-18

- [x] Draft literature-backed pea and soy meaty-positive benchmark candidate reports so the missing data package is concrete rather than abstract.
- [x] Promote those candidate reports into visible documentation instead of leaving them buried in research notes.

### SLR benchmark triage — 2026-03-18

- [x] Incorporate the SLR verdict that the missing pea/SPI meaty-positive matrix benchmarks are structural literature gaps, not a search gap.
- [x] Register the benchmark-intake status of Nishimura 2024, Squeo 2023, Asen 2022, and Malia 2025 in a machine-readable registry instead of leaving them only in prose.
- [x] Re-anchor pea process-state calibration provenance to Asen 2022 and Malia 2025 in the matrix correction layer.
- [x] Create explicit machine-readable payloads for the pea process-state calibration anchors and the Squeo 2023 safety reference anchor.
- [x] Create a canonical human-readable scientific reference document that links pathways, validated articles, numeric anchors, and comments in one place.
- [x] Recover full text for Nishimura & Abe (2024) and verify whether absolute MFT/FFT concentrations plus internal-standard details are present. Result: full text confirms soy-hydrolysate qualitative sulfur chemistry at 95 C / 90 min, but only as relative peak-area output without absolute ppb or internal-standard benchmark quantitation. (Dificultad: Media — Modelo recomendado: GPT-5.4)
- [x] Convert the candidate reports into a benchmark-intake checklist with only fully curated references before creating benchmark JSON entries. (Dificultad: Alta — Modelo recomendado: GPT-5.4)

### What is already strong

- [x] Free-amino-acid PRIMARY benchmarks are in a credible Docker-validated envelope.
- [x] FAST proxy vs observable projection semantics are explicit and benchmark artifacts expose the projection provenance.
- [x] Pea/soy matrix handling is materially better than before: accessibility, denaturation inference, headspace fallback, and compound-class-aware retention now exist.
- [x] The repo now has validated-envelope and formulation explainability artifacts.
- [x] The Docker-first validation workflow is coherent and reproducible.

### What is still not enough for a best-in-class scientific tool

- [ ] External meaty-positive matrix benchmark data is still missing for pea and soy, and the SLR now confirms this is a structural literature gap rather than an unresolved search problem, so matrix trust cannot advance beyond directional support. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Matrix benchmarks are executable intake checks, not yet benchmark-facing target-ranking validations. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Process-state calibration for accessibility, retention, and release is not yet benchmark-closed across real plant matrices; the SLR only upgrades pea denaturation/SH anchors, not the full PPI/SPI benchmark conditions. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Benchmark breadth is narrow relative to the scientific space alternative-protein teams actually care about. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] There is still a gap between “validated model component” and “formulation operating system for scientists”, especially for experiment ingestion and campaign-scale comparison. (Dificultad: Alta — Modelo recomendado: Gemini 3 Flash)

## Strategic Assessment

The previous summary was directionally correct but not sufficient as a roadmap to a genuinely state-of-the-art repository.

It correctly identified immediate usability gaps:

- inline CLI explainability;
- domain-of-validity warnings;
- compound-aware matrix observability;
- better user-facing diagnostics.

But a best possible tool for scientists also needs broader priorities that are at least as important:

- benchmark expansion beyond the current narrow envelope;
- uncertainty and confidence reporting;
- stronger product ergonomics for comparing candidate formulations;
- scientist-centered outputs that connect predictions to concrete intervention decisions;
- data-ingestion and calibration workflows that allow the model to improve from new experiments.

## Priority Order

## P0 — Make The Existing Science Usable In Daily Work

These are the highest-leverage items because they convert existing validated machinery into something scientists can actually trust and use without reading internal docs.

- [x] Integrate deep explainability into CLI output (bottlenecks, penalties, intervention hints) (P0)
- [x] Print a concise default decision summary for every run: dominant desirable compounds, dominant penalties, matrix state, validated-envelope status, and key caveats. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [x] Surface domain-of-validity warnings inline whenever a run leaves the trusted envelope (`matrix_only`, unsupported matrices, peptide-bound systems, aggressive process assumptions, sparse benchmark analogies). (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)
- [x] Add a single scientist-facing report mode that outputs both machine-readable JSON and human-readable Markdown from the main CLI, not only from side scripts. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [x] Add recommendation diagnostics that explicitly answer: which precursors helped, which penalties dominated, which compounds were suppressed by matrix retention/headspace, and what intervention is most likely to move the outcome. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [x] Create side-by-side comparison outputs for multiple candidate formulations so scientists can compare tradeoffs instead of inspecting one run at a time. (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)

### P0 acceptance criteria

- [x] A scientist can run one command and immediately see what the model predicts, why, how trustworthy it is, and what to try next. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [x] The default CLI no longer requires reading `results/validation/*.md` to understand whether a result is meaningful. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)

## P1 — Turn Matrix Handling Into A Real Validation Surface

This is the biggest scientific-product gap. Without it, the repo is credible for a narrow envelope but not yet best-in-class for alternative proteins.

- [ ] Promote plant-matrix benchmarks from intake checks toward target-ranking benchmarks with benchmark-facing observable outputs. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Expand matrix observability from compound classes toward compound-specific calibration where literature or internal data exists. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Add matrix benchmark families beyond the current pea/soy intake cases: broader protein systems, process states, and off-flavour/meaty tradeoff regimes relevant to alternative proteins. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Validate matrix predictions against benchmarks that include both beneficial sulfur/meaty targets and adverse lipid-oxidation targets in the same system; the SLR confirms this tradeoff benchmark does not currently exist in the literature for PPI/SPI isolate systems. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Add process-state realism for extrusion/heating history where that state materially changes accessibility or release. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Keep strict-gate promotion blocked until matrix target ranking is reproducible in Docker and visible in reports. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)

### P1 acceptance criteria

- [ ] Matrix systems can be benchmarked on ranked observable targets, not just “executable intake”. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] A pea/soy recommendation is scientifically interpretable as a formulation signal, not only as a conservative heuristic. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)

## P2 — Add Confidence, Uncertainty, And Decision Boundaries

Scientists need to know not only the prediction, but how much to trust it.

- [x] Add per-result confidence metadata tied to benchmark neighborhood, matrix support level, and projection assumptions. (Dificultad: Media — Modelo recomendado: Gemini 3 Pro)
- [x] Report uncertainty bands or at least confidence tiers for major predicted compounds and aggregate sensory scores. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
- [x] Add calibration diagnostics that show when a recommendation is extrapolating beyond supported chemistry or process conditions. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [x] Distinguish clearly between benchmark-supported quantitative predictions, directional heuristics, and speculative outputs. (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)
- [x] Add sensitivity summaries showing which inputs most strongly change the ranking or safety outcome. (Dificultad: Media — Modelo recomendado: Gemini 3 Pro)

### P2 acceptance criteria

- [x] Every recommendation carries an explicit confidence story. (Dificultad: Media — Modelo recomendado: Gemini 3 Pro)
- [x] Scientists can tell the difference between “ship this to the wet lab first” and “interesting but speculative”. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)

### P1 detailed execution plan

#### Phase 1 — Benchmark contract expansion

1. [x] Define a matrix target-ranking contract that every candidate benchmark must satisfy: observable targets, expected ordering, adverse markers, process metadata, and citation provenance. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
2. [x] Add contract fields to the benchmark artifacts and validation loaders so matrix cases expose ranked desirable and adverse targets, not just intake executability. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
3. [x] Create a first promotion set for pea and soy with at least one meaty-positive and one off-flavour-negative benchmark per family. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)

#### Phase 2 — Observable calibration surface

1. [x] Replace class-level matrix observability where evidence exists with compound-level calibration entries keyed by compound, matrix family, and process state. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
2. [x] Extend projection metadata so every ranked compound exposes its calibration source, fallback mode, and evidence strength in CLI/report JSON. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
3. [ ] Add explicit adverse-target calibration for lipid-coupled notes so matrix validation covers the meaty/off-flavour tradeoff in the same experiment. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
4. [ ] Encode Asen 2022 and Malia 2025 as explicit process-state calibration payloads instead of leaving them only as provenance strings in the pea thermal envelope. (Dificultad: Alta — Modelo recomendado: GPT-5.4)

#### Phase 3 — Process realism

1. [x] Introduce process-state descriptors for extrusion or pre-heating history into matrix benchmark inputs and explainability payloads. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
2. [ ] Calibrate accessibility and release heuristics against those states before promoting any matrix family into a stricter gate. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)

#### Phase 4 — Validation and promotion

1. [x] Add Docker-visible ranking assertions for matrix benchmarks: top-k desirable hits, adverse marker ordering, and tolerance thresholds. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
2. [x] Add report outputs that show benchmark deltas for every promoted matrix case, including successes, misses, and fallback paths used. (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)
3. [x] Keep strict-gate promotion disabled until the first pea/soy candidate passes reproducibly in Docker on ranking metrics, not only execution checks. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)

#### P1 deliverables

- [x] New/updated matrix benchmark YAML cases with ranked targets and adverse markers. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
- [x] Projection/calibration registry for compound-level matrix observability. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [x] Validation code and reports that surface ranking deltas in Docker artifacts. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)

#### P1 verification plan

- [x] Add focused pytest coverage for contract parsing, calibration fallback behavior, and ranking assertions. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Run Docker benchmark lanes for the first promoted matrix family before merge. (Dificultad: Media — Modelo recomendado: GPT-5 mini)
- [x] Compare default-branch vs feature-branch benchmark deltas before any strict-gate claim. (Dificultad: Media — Modelo recomendado: GPT-5 mini)

## P3 — Broaden The Scientific Envelope That Matters To Alt-Protein R&D

This is how the repo stops being a strong niche engine and becomes genuinely state of the art.

- [ ] Expand benchmark coverage for precursor families, protein systems, process regimes, and target compound families that matter in alternative proteins. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Add peptide-bound and intact-protein reactivity where it materially changes accessible chemistry. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Add broader carbohydrate realism beyond the current free-sugar emphasis when relevant to commercial formulations. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Strengthen safety validation with literature-backed dynamic acrylamide and related risk markers across realistic process trajectories; Squeo 2023 now provides an industrial endpoint anchor but not dynamic T/t/[Asn] kinetics. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Improve temporal/process-path validation so the model can compare alternative heating profiles, not only endpoint conditions. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Add better coverage for lipid/Maillard coupling in realistic plant-fat systems. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)

### P3 acceptance criteria

- [ ] The validated envelope covers a materially broader slice of real alternative-protein formulation space. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Scientists can use the tool for more than free-precursor exploration plus narrow pea/soy intake modeling. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)

## P4 — Build The Scientific Operating System Layer

This is the difference between a strong model repo and a tool scientists keep using.

- [ ] Add experiment-ingestion workflows so internal wet-lab data can be added as calibration/validation cases without bespoke code edits. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
- [ ] Add standardized project/report templates for comparing candidate recipes, interventions, and benchmark deltas over time. (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)
- [ ] Add versioned result provenance so scientists can trace a recommendation to the exact benchmark contract, calibration state, and code revision. (Dificultad: Media — Modelo recomendado: Gemini 3 Pro)
- [ ] Add batch campaign workflows for screening formulation sets and exporting ranked decisions with rationale. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)
- [ ] Add notebook- and report-friendly APIs for external scientific teams. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [ ] Add explicit support for comparing model output against newly collected GC-MS or sensory panel data. (Dificultad: Alta — Modelo recomendado: Gemini 3 Pro)

### P4 acceptance criteria

- [ ] A team can use the repo as a repeatable formulation platform, not just as a codebase with scripts. (Dificultad: Media — Modelo recomendado: Gemini 3 Flash)

## Guardrails

- [ ] Do not chase new chemistry breadth before the user-facing decision surface is honest and usable. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [ ] Do not promote matrix systems into the strict gate before observable target ranking is reproducible. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [ ] Do not claim quantitative confidence where only directional support exists. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)
- [ ] Keep Docker as the authoritative validation environment. (Dificultad: Baja — Modelo recomendado: GPT-5 mini)

## Immediate Execution Sequence

1. [x] Recover full text for Nishimura & Abe (2024) and decide whether it can be promoted from soy-hydrolysate conditional calibration to an encoded literature anchor. Result: it can be promoted only to a qualitative intake anchor, not to an absolute benchmark or calibrated ppb reference. (Dificultad: Media — Modelo recomendado: GPT-5.4)
2. [x] Encode Squeo 2023 as a formal safety-reference payload and expose its upper-bound acrylamide range in a code-consumable safety reference surface without pretending it is a kinetic benchmark. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
3. [x] Encode Asen 2022 and Malia 2025 as explicit pea process-state calibration payloads for `denaturation_state` and `cysteine_accessibility`. (Dificultad: Alta — Modelo recomendado: GPT-5.4)
4. [ ] Design the first primary cross-benchmark experiment for PPI/SPI + ribose + cysteine because the SLR confirms the key meaty/off-flavour/time-series gaps cannot be closed with more literature search; use the intake checklist plus the minimum primary experiment spec as the starting contract. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)
5. [ ] Build the ingestion path for new internal calibration and validation experiments. (Dificultad: Muy alta — Modelo recomendado: GPT-5.4)

## Review

- This backlog replaces the older phase-heavy plan with a scientist-product roadmap.
- The current repo is already scientifically serious and more usable than a typical research prototype.
- It is not yet state of the art for alternative-protein scientists because the decisive missing layers are user-facing decision support, confidence reporting, matrix-ranking validation, and experiment-to-model iteration.
- P2 is now functionally complete in the main scientist-facing surfaces: CLI/report outputs include prediction mode, calibration diagnostics, compound and aggregate confidence tiers, and local sensitivity summaries.
- P1 is now the dominant scientific-product blocker and should start with benchmark contract design before any broader chemistry expansion.
- P1 progress on 2026-03-17: pea/soy matrix benchmarks now carry an explicit ranking contract with process metadata and citation provenance; validation/report artifacts surface ranking-contract status and compound-level calibration provenance while the strict gate remains blocked.
- P1 evidence status on 2026-03-17: the repo now distinguishes external matrix evidence from internal reproducibility harnesses explicitly; external pea/soy meaty-positive promotion remains blocked because no wet-lab quantitative benchmark is currently available in-repo.
- P1 completion status on 2026-03-17: the repo now emits matrix assertions, evidence audits, family-readiness summaries, and branch-vs-main deltas in Docker. Pea and soy each have an off-flavour anchor and a meaty candidate, but external meaty-positive assessment remains blocked pending wet-lab data.
- SLR update on 2026-03-18: the literature search closes the planning uncertainty. Nishimura 2024 is the only soy-protein matrix chemistry candidate for MFT/FFT, Squeo 2023 is the only benchmark-eligible acrylamide endpoint anchor, Asen 2022 and Malia 2025 are the best open pea process-state calibration anchors, and the core PPI/SPI meaty/off-flavour/time-series benchmark gaps are structural and require primary data.
- Nishimura full-text update on 2026-03-18: the paper is now resolved as a qualitative soy-hydrolysate intake anchor. It confirms the soy matrix chemistry workflow and MRP conditions, but its volatile outputs are relative peak-area signals rather than absolute benchmark-ready concentrations.

## Model guidance (resumen)

- `GPT-5 mini`: cambios tácticos de CLI, scripts, y reportes; pruebas y ajustes triviales.
- `Gemini 3 Flash`: integraciones de ingeniería productiva y ayudantía para tareas de complejidad media.
- `Gemini 3 Pro`: diseño de pipelines, validación de proyección/headspace, y pruebas de integración.
- `GPT-5.4` / `Claude Opus 4.5`: calibración científica profunda, diseño de benchmarks y pipelines de ingestión de datos experimentales.
