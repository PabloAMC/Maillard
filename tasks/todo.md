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

## Active Roadmap

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

Phase P4.0: Adopt explicit model-selection rules

- [ ] Do not use MLPs as the primary closure path for plant-matrix accuracy; literature-backed observable closure stays first.
- [ ] Restrict default consideration to chemistry-domain external models before any materials-first model is allowed into the molecular workflow.
- [ ] Require each proposed MLP use case to declare the target motif class, the expected speedup, the expected failure mode, and the validation comparator.

Phase P4.1: Use MLPs only in narrow offline roles

- [ ] Restrict ML-potential usage to conformer screening, geometry pre-optimization, local barrier surrogates, or reusable offline acceleration artifacts.
- [ ] Keep MLP outputs out of the default prediction surface until they have been reduced to validated cached artifacts.
- [ ] Expose in reports when a result depends on computational refinement rather than direct literature or benchmark support.

Phase P4.2: Start from the strongest chemistry-domain shortlist

- [ ] Evaluate chemistry-first external models such as MACE organic variants, OrbMol, and AIMNet2 before considering materials-first defaults.
- [ ] Avoid using materials-only foundation models as the default molecular refinement path for Maillard chemistry.
- [ ] Add a benchmarked adoption note for any selected MLP, including where it helps and where it is known to fail.

P4 measurable exit criteria:

- [ ] No MLP is used as a substitute for missing matrix benchmarks.
- [ ] Any adopted MLP serves a clearly bounded offline role with explicit provenance.
- [ ] The chosen MLP shortlist is justified against chemistry-domain needs rather than generic SOTA branding.

### P5. Make Trust Visible In Every Decision Surface

Goal: ensure a scientist can see why a prediction should or should not be trusted.

- [x] Surface benchmark neighborhood, calibration source, observable-layer assumptions, and extrapolation axes in all default reports.
- [x] Show which compounds are benchmark-backed, transferred, heuristic, or computationally refined.
- [x] Distinguish clearly between quantitative recommendation mode and directional hypothesis mode.

## Recent Verification

- [x] Explicit decision-panel evidence contract propagated through projection metadata, benchmark snapshots, confidence rows, and reporting surfaces.
- [x] Runtime reports now include a literature evidence summary derived from the existing benchmark intake registry and structural-gap payloads.
- [x] Docker validation passed on the affected suite: 36 tests green across benchmark targets, matrix delta reports, SLR payloads, usability reports, trust visibility, and projection contract coverage.
- [x] Literature-only learning loop now generates results/validation/literature_learning_loop.{md,json} and results/validation/literature_runtime_templates.json from the canonical intake registry.
- [x] Matrix priors now surface uncertainty posture and process-state applicability, and mycoprotein is encoded as a first-class bounded matrix family.
- [x] Docker validation passed on the expanded L1/L2/L3 suite: 66 tests green across the new learning-loop coverage, priors, reporting, matrix targets, and trust surfaces.
- [x] Matrix target support status now generates results/validation/matrix_target_status.{md,json} and separates quantitative closure, internal candidates, directional support, and open gaps at the compound level.
- [x] Refinement watchlist now generates results/validation/refinement_watchlist.{md,json} plus results/validation/offline_dft_jobs.json to keep expensive refinement benchmark-facing and offline.
- [x] Docker validation passed on the expanded L4/P0/P1 suite: 71 tests green across matrix status, refinement watchlist, benchmark evidence, target snapshots, accessibility, reporting, and projection contracts.
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
