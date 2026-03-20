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

### P0. Make Matrix Predictions Scientist-Useful

Goal: promote pea and soy from intake checks to matrix target-ranking systems.

Phase P0.1: Define the benchmark target panel

- [ ] Freeze the compound panel for matrix target-ranking: sulfur positives, Strecker aldehydes, pyrazines, furans/furanones, and adverse lipid markers.
- [ ] Map each panel compound to one of four evidence states: externally benchmarked, internally benchmarked, transferred prior, or still missing.
- [ ] Publish the target panel and evidence map in machine-readable form so reports and validation use the same contract.

Phase P0.2: Build benchmark-grade pea and soy targets

- [ ] Create at least one pea benchmark and one soy benchmark where desirable and adverse compounds are ranked together under a single process state.
- [ ] Ensure each benchmark has explicit measured values, ranking intent, adverse markers, and scale thresholds where justified.
- [ ] Add focused tests proving that the benchmark executes, matches all target compounds, and preserves the intended ordering.

Phase P0.3: Close matrix observability for decision-driving compounds

- [ ] Add observable anchors for the compounds scientists actually use to decide, not only hexanal-like off-notes.
- [ ] Separate compound classes that are quantitatively anchored from those that remain directional or transferred.
- [ ] Make projection metadata expose which matrix compounds are directly anchored versus transferred.

Phase P0.4: Promote matrix validation from intake-only to target-ranking

- [ ] Extend matrix benchmark contracts so they score multi-compound target-ranking behavior, not just executable intake behavior.
- [ ] Add summary artifacts that show which matrix targets are quantitatively closed, directionally supported, or still open.
- [ ] Define a measurable promotion rule for moving a matrix benchmark family from directional to near-quantitative support.

P0 measurable exit criteria:

- [ ] At least one pea and one soy matrix benchmark rank sulfur, Strecker, pyrazine, and adverse markers together.
- [ ] All compounds in the panel have an explicit evidence label in reports.
- [ ] Matrix summary artifacts distinguish quantitative support from directional support at the compound level.

### P1. Model Accessibility, Not Just Chemistry

Goal: capture the main reason plant matrices differ from free systems.

Phase P1.1: Encode accessibility states explicitly

- [ ] Define a small set of accessibility states that the runtime understands: free, peptide-bound, protein-embedded, and process-opened.
- [ ] Map current pea and soy assumptions onto these states instead of burying them in one-off correction factors.
- [ ] Add tests showing that changing accessibility state changes only the intended observables.

Phase P1.2: Build reusable process-state surrogates

- [ ] Add reusable surrogates for denaturation, release, retention, and browning severity that depend on temperature, time, pH, and matrix family.
- [ ] Keep these surrogates cheap enough for daily runtime and explicit enough to surface in reporting.
- [ ] Validate that process-state changes improve matrix benchmark behavior without destabilizing free-precursor benchmarks.

Phase P1.3: Couple accessibility to prediction outputs

- [ ] Propagate accessibility state and effective process state into the projection metadata used by reports.
- [ ] Surface chemically reachable versus merely plausible compounds in scientist-facing artifacts.
- [ ] Add domain-of-applicability warnings when accessibility assumptions dominate the prediction.

P1 measurable exit criteria:

- [ ] Accessibility state is a first-class runtime concept rather than an implicit correction.
- [ ] Process-state surrogates are benchmarked for pea and soy on the compounds that drive decisions.
- [ ] Reports show when accessibility, not chemistry, is the main source of uncertainty.

### P2. Open The Extrusion Regime Carefully

Goal: make extrusion-heavy predictions useful without overclaiming quantitative closure.

Phase P2.1: Define the extrusion neighborhood

- [ ] Write an explicit extrusion-heavy benchmark neighborhood definition: what counts as in-domain, near-domain, and out-of-domain.
- [ ] Identify the minimum observable panel needed for extrusion usefulness: at least key meaty positives, key off-notes, and one process-severity marker.
- [ ] Encode default warnings so extrusion outputs cannot be misread as free-system quantitative predictions.

Phase P2.2: Add process-structured extrusion surrogates

- [ ] Add cheap surrogates for moisture redistribution, accessibility loss/recovery, and volatile retention under high-severity processing.
- [ ] Reuse P1 accessibility and process-state concepts instead of creating a separate extrusion-only stack.
- [ ] Ensure extrusion surrogates degrade confidence cleanly rather than silently reusing free-precursor assumptions.

Phase P2.3: Create exploratory extrusion mode

- [ ] Add an explicit exploratory recommendation mode for extrusion-heavy systems.
- [ ] Report which outputs are still benchmark-backed from lower regimes and which are new extrusion extrapolations.
- [ ] Add acceptance tests that confirm extrusion mode always carries the intended warnings and confidence posture.

P2 measurable exit criteria:

- [ ] Extrusion-heavy runs use an explicit exploratory mode with visible warnings.
- [ ] The runtime reuses shared accessibility/process abstractions rather than duplicating logic.
- [ ] Scientists can distinguish reused lower-regime support from true extrusion-specific support in the default reports.

### P3. Refine Only Decisive Mechanistic Gaps

Goal: spend expensive compute only where it changes decisions.

- [ ] Use sensitivity analysis to identify which motif classes actually limit matrix prediction quality.
- [ ] Run xTB first and selective DFT second on only those decisive motif classes.
- [ ] Write every expensive correction back as a reusable artifact consumed by the cheap daily workflow.

### P4. Keep MLPs In Their Proper Role

Goal: use ML potentials where they are elegant, not where they are fashionable.

- [ ] Do not use MLPs as the primary closure path for plant-matrix accuracy.
- [ ] Use external state-of-the-art ML potentials only in offline roles once the target motif class is clearly defined.
- [ ] Restrict ML-potential usage to conformer screening, local delta models, or reusable offline acceleration artifacts.

### P5. Make Trust Visible In Every Decision Surface

Goal: ensure a scientist can see why a prediction should or should not be trusted.

- [ ] Surface benchmark neighborhood, calibration source, observable-layer assumptions, and extrapolation axes in all default reports.
- [ ] Show which compounds are benchmark-backed, transferred, heuristic, or computationally refined.
- [ ] Distinguish clearly between quantitative recommendation mode and directional hypothesis mode.

## Success Criteria

- [ ] A scientist can use the tool to narrow a wet-lab campaign, not just inspect simulated chemistry.
- [ ] Free-precursor predictions remain quantitatively stable while matrix modeling expands.
- [ ] Pea and soy predictions become target-ranking useful for multi-compound flavor decisions.
- [ ] Extrusion-heavy systems become operationally useful as exploratory planning tools.
- [ ] Expensive computation remains offline, sparse, and reusable.
- [ ] Reports state confidence and extrapolation explicitly enough that overclaiming becomes difficult.
