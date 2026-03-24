# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Recommended](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Maillard is a computational screening framework for scientists who want to imitate meat-like Maillard chemistry in plant-based systems before running wet-lab campaigns.

The point is not to predict every food matrix equally well. The point is to help scientists answer the most useful pre-lab question:

Which formulation and process changes are most worth testing next if the goal is meat-like aroma under plant-matrix constraints?

Important note: this project requires conda or Docker. Pip-only installation is not supported.

## Quick Start

If you want to run the tool immediately, start with [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md).

If you want definitions for terms such as FAST mode, validated envelope, or benchmark neighborhood, see [docs/guides/GLOSSARY.md](docs/guides/GLOSSARY.md).

If you want the current matrix-family scope and the experiment-ingestion workflow, see [results/validation/matrix_family_coverage.md](results/validation/matrix_family_coverage.md) and [docs/guides/MATRIX_EXPERIMENT_INGESTION.md](docs/guides/MATRIX_EXPERIMENT_INGESTION.md).

## The Real Problem

The most useful version of this tool is not a universal chemistry oracle. It is a scientist-facing decision system that can:

- rank candidate formulations before the wet lab
- explain why a prediction should be trusted or discounted
- separate benchmark-backed claims from directional extrapolation
- show which compounds are driven by free chemistry, matrix observability, or transferred priors

That means the central problem is broader than matching one benchmark table. We need a model that is useful across three regimes of scientific confidence.

## The Three Regimes

| Regime                  | What it represents                                                         | Current trust | What the tool should do well                                                  |
| ----------------------- | -------------------------------------------------------------------------- | ------------- | ----------------------------------------------------------------------------- |
| Free precursors         | Buffer-like systems where sugars and amino acids are directly accessible   | High          | Quantitative ranking and concentration-scale screening                        |
| Pea / soy matrices      | Real protein matrices where accessibility, retention, and headspace matter | Moderate      | Directional and near-quantitative prioritization where explicit anchors exist |
| Extrusion-heavy systems | Highly processed systems with strong structural and transport effects      | Low           | Exploratory hypothesis generation and experiment prioritization               |

These regimes are not three unrelated products. They are three layers of increasing structural complexity.

- Free precursors test whether the chemistry core is right.
- Pea and soy matrices test whether observable reality is right.
- Extrusion-heavy systems test whether process-structured accessibility and transport are right.

## Family Validation Surface

The important question for a scientist is not whether the repository names many chemistry families. The important question is whether each family has an explicit experimental prediction surface.

This repo therefore separates three cases:

- families with **compound-level quantitative parity** against experiments
- families with **benchmark-linked calibration support** but not yet strict quantitative closure
- families that remain **directional or gap-limited**, and are reported as such rather than being overclaimed

The figures below are generated with `./scripts/docker_maillard.sh validation-figures` and are provided as three standalone PNGs so they can be embedded, cropped, or cited independently.

![Compound parity](docs/assets/family_parity.png)

![Per-benchmark accuracy](docs/assets/family_benchmark_accuracy.png)

Captions:

- **Compound parity:** Predicted vs measured concentrations (log–log). Colour = chemistry family; marker shape = execution lane. Green/yellow bands denote 1.5× and 2× tolerances.
- **Per-benchmark accuracy:** Worst-case predicted/measured ratio per benchmark (human-readable study labels). Vertical lines mark strict-gate (1.5×) and matrix tolerance (2×).
- **Family coverage:** Counts of matched quantitative compound points across all 10 families; families with no benchmark-backed points are annotated as explicit gaps.

If your markdown viewer does not render images inline, open the files directly in `docs/assets`.

For machine-readable artifacts, see [results/validation/family_validation_overview.md](results/validation/family_validation_overview.md). For the detailed per-benchmark drill-down see [docs/assets/validation_overview.png](docs/assets/validation_overview.png).

## How To Use This For Alternative Protein Research

For a scientist evaluating alternative proteins, the operational workflow is:

1. Select matrix family and process state explicitly (for example pea isolate, soy isolate, mycoprotein, extrusion-heavy process).
2. Run a forward prediction and generate report artifacts.
3. Check the matrix-family coverage artifact so you know whether your family is explicit support, indirect support, or an open gap.
4. Read family evidence ladder and family lane sensitivity before trusting absolute concentrations.
5. Use benchmark-backed families for quantitative decisions and directional families for experiment prioritization.
6. Promote a family lane only after adding benchmark or calibration evidence, not by tuning barriers alone.
7. When you obtain a new measurement set, compare it against the model through the experiment-intake workflow before treating it as promotion evidence.

Practical commands:

```bash
./scripts/docker_maillard.sh validation-figures
python scripts/run_pipeline.py --protein-type pea_iso --target meaty --report --output-dir results/first_run
```

Primary artifacts for this workflow:

- [results/validation/family_validation_overview.md](results/validation/family_validation_overview.md)
- [results/validation/family_lane_validation.md](results/validation/family_lane_validation.md)
- [results/validation/family_deviation_audit.md](results/validation/family_deviation_audit.md)
- [results/validation/benchmark_summary.md](results/validation/benchmark_summary.md)
- [results/validation/validated_envelope.md](results/validation/validated_envelope.md)
- [results/validation/matrix_family_coverage.md](results/validation/matrix_family_coverage.md)
- [results/validation/matrix_target_status.md](results/validation/matrix_target_status.md)
- [results/validation/matrix_promotion_contract.md](results/validation/matrix_promotion_contract.md)
- [results/validation/matrix_observable_closure_audit.md](results/validation/matrix_observable_closure_audit.md)
- [results/validation/matrix_experiment_intake_schema.md](results/validation/matrix_experiment_intake_schema.md)
- [results/validation/p3_refinement_governance.md](results/validation/p3_refinement_governance.md)
- [results/validation/p4_mlp_assessment.md](results/validation/p4_mlp_assessment.md)
- [results/validation/p4_adoption_notes.md](results/validation/p4_adoption_notes.md)
- [results/validation/literature_learning_loop.md](results/validation/literature_learning_loop.md)
- [results/validation/family_promotion_state.md](results/validation/family_promotion_state.md)

## Predictive Accuracy: What Is Quantitative Today

Current artifact-backed status (from results/validation):

- 10 chemistry families tracked in runtime scope
- 4 families currently benchmark-linked
- 4 families currently have compound-level quantitative parity points
- 32 quantitative compound points in the current validation surface

Interpretation:

- Quantitative prediction claims are strongest in benchmark-backed families and free-precursor strict-ready lanes.
- Matrix and support families can still be decision-useful, but some remain calibration-grade or directional.
- This is a deliberate trust posture: explicit boundaries are preferred over overclaiming broad quantitative closure.

## Current Validation Snapshot

Current in-repo validation summary:

- 9 supported benchmarks are tracked in Docker-validated artifacts
- 4 benchmarks are strict-ready today, all inside the free-precursor envelope
- 9 matched compounds define the current authoritative quantitative proof surface
- the median matched-compound ratio in that proof surface is 1.118x

The authoritative benchmark-level parity plot (per compound, per benchmark, coloured by study) is shown below.

![Validation Overview](docs/assets/validation_overview.png)

How to interpret this trust surface:

- Free-precursor benchmarks are the quantitative proof surface.
- Pea and soy matrix paths are useful for prioritization, not yet for broad release-grade claims.
- Soy and pea are not the only matrix families tracked by the repo, but they are the only plant-protein matrix lanes with executable benchmark-plus-calibration support today.
- Family 07 carbonyl donor hierarchy is now promoted to benchmark-linked support, meaning sugar identity is no longer only a heuristic lane: existing benchmark-linked compounds constrain it with explicit uncertainty, but it is not yet near-quantitative as a standalone family.
- P3 mechanistic refinement is now explicitly gate-kept by benchmark-visible compounds and cheap-first screening; if [results/validation/p3_refinement_governance.md](results/validation/p3_refinement_governance.md) shows zero approved jobs, offline QM stays parked.
- P4 MLP work remains an offline accelerator lane only. The current policy in [results/validation/p4_mlp_assessment.md](results/validation/p4_mlp_assessment.md) and [results/validation/p4_adoption_notes.md](results/validation/p4_adoption_notes.md) is no default MLP adoption until the reaction benchmark passes.
- Mycoprotein is currently bounded prior support, soy hydrolysate remains qualitative intake support, and other plant proteins remain explicit scope gaps until elevated into runtime-facing evidence.
- P6 matrix expansion is intentionally bounded: [results/validation/matrix_family_coverage.md](results/validation/matrix_family_coverage.md) now separates bounded expansion candidates from evidence-blocked matrix families so scope cannot drift faster than the evidence surface.
- Extrusion-heavy systems remain exploratory until benchmarked directly.
- For the full 10-family strategic view including coverage gaps, see the **Family Validation Surface** section above.

## What Accuracy Depends On

The tool does not have a single accuracy mode. Its precision depends on which layer is doing most of the work.

| Layer                                    | Main inputs                                                                                        | What it mainly controls                                         |
| ---------------------------------------- | -------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| Core FAST kinetics                       | Literature barriers, cached barrier records, calibrated Arrhenius heuristics                       | Pathway ordering and relative competitiveness                   |
| Observable projection                    | Proxy-to-observable mapping, volatility, retention, process-sensitive projection factors           | Absolute ppb-scale fit                                          |
| Matrix calibration and headspace release | Protein type, process state, compound-specific observable factors, transferred anchors             | Whether matrix predictions stay physically plausible and useful |
| Offline refinement                       | xTB, selective DFT, and optional external ML-potential refinement written back as cached artifacts | Narrow mechanistic correction of decisive motif classes         |

Two rules follow from this:

- We should not retune energy barriers just to make one benchmark concentration match.
- When ranking looks right but ppb scale is wrong, the preferred fix is usually in the observable or matrix-calibration layer.

The current validation contract therefore uses two scale checks:

- max ratio, to catch worst-case outliers
- mean absolute log-scale error, to catch broad multiplicative drift

## Scientific Architecture

The scientific stack is deliberately layered. Different tools serve different roles.

| Layer                          | Main tool or model                                         | Role in the system                                                                                                      | Main dependency                                                             |
| ------------------------------ | ---------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Reaction enumeration           | SMIRKS rules and curated pathway families                  | Generate plausible Maillard and lipid-derived chemistry                                                                 | Coverage and correctness of the encoded reaction families                   |
| Fast prediction core           | FAST observable path with empirical or cached barriers     | Daily screening, ranking, and benchmark-facing concentrations                                                           | Calibration quality of barrier tables and projection layers                 |
| Diagnostic thermochemistry     | Cantera and thermodynamic gating                           | Diagnose whether a pathway family is physically plausible and whether gating would materially change benchmark behavior | Thermodynamic data and gating policy                                        |
| Observable projection          | Headspace, retention, matrix, and process-state surrogates | Convert pathway signal into what a scientist would actually measure or smell                                            | Benchmark-anchored observable calibration                                   |
| Offline mechanistic refinement | xTB first, selective DFT second                            | Refine decisive motif classes when cheap surrogates stop improving decisions                                            | Availability of narrow, benchmark-relevant refinement targets               |
| Optional ML acceleration       | External state-of-the-art ML potentials, only offline      | Accelerate conformer or local motif refinement after the refinement task is well-defined                                | Quality and relevance of the external model for the motif class of interest |

What this means in practice:

- FAST is the day-to-day engine because it is cheap enough to screen formulations and calibrate against benchmarks.
- Cantera is a diagnostic lane, not the main benchmark-facing prediction surface.
- xTB and selective DFT are for targeted refinement, not routine production inference.
- If ML potentials are used, the right default is to consume a state-of-the-art external model in an offline role, not to train a repository-specific model by default.

## Should ML Potentials Be The Main Missing Piece?

No.

ML potentials are useful, but they are not the elegant primary solution to the current product gap.

The main missing problem is not raw quantum precision. The main missing problem is that plant-matrix and extrusion-heavy systems are under-benchmarked in accessibility, retention, headspace release, and process-state dependence.

That means the clean modeling order is:

1. benchmark and calibrate observable-layer behavior where literature or internal data exist
2. add process-aware matrix accessibility and transport surrogates
3. use xTB and selective DFT to refine narrow motif classes only after cheap closure stops improving decisions
4. use external state-of-the-art ML potentials later, offline, only if repeated motif classes justify them and the target motif class is already well-scoped

ML potentials should therefore be a tertiary offline accelerator, not the main path to making the tool useful, and the default assumption should be to reuse the best available external model rather than train a new one ourselves.

## What The Best Version Of This Tool Looks Like

The best version of Maillard is the tool that helps a scientist decide what to make next, not the tool that claims uniform absolute accuracy everywhere.

To get there, the repo needs to do five things well:

- keep the free-precursor chemistry quantitatively solid
- promote pea and soy from intake checks to matrix target-ranking benchmarks
- model extrusion-heavy systems as process-structured accessibility problems rather than pretending they are just bigger free systems
- expose confidence, extrapolation axes, and calibration sources clearly in reports
- keep expensive computation offline and reusable

## Recommended Workflow

1. Check the validation surface first.
2. Run a forward prediction for your candidate formulation.
3. Save a report bundle with provenance.
4. Use the optimizer to search nearby operating points.
5. Treat matrix-heavy and extrusion-heavy outputs according to their trust regime.

## Key Commands

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-promotion-contract
./scripts/docker_maillard.sh matrix-closure-audit
./scripts/docker_maillard.sh experiment-intake-schema
./scripts/docker_maillard.sh matrix-family-coverage
./scripts/docker_maillard.sh p3-refinement
./scripts/docker_maillard.sh p4-mlp-assessment
./scripts/docker_maillard.sh literature-learning-loop
./scripts/docker_maillard.sh family-promotion-state
./scripts/docker_maillard.sh compare-experiment data/protocols/example_matrix_experiment_intake.yaml
```

Generate one prediction:

```bash
python scripts/run_pipeline.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --time-minutes 45 \
  --protein-type pea_iso \
  --target meaty \
  --minimize beany \
  --report \
  --output-dir results/first_run
```

Optimize a formulation:

```bash
python scripts/optimize_formulation.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --target-tag meaty \
  --minimize-tag beany \
  --protein-type pea_iso \
  --n-iterations 50
```

## Scientific References

The canonical human-readable reference list is [docs/reference/SCIENTIFIC_REFERENCE.md](docs/reference/SCIENTIFIC_REFERENCE.md).

For the broader literature-screening and ingestion view, see [docs/slr_benchmark_evaluation.md](docs/slr_benchmark_evaluation.md).

## Current Boundary

Use the generated validation artifacts when you need the exact benchmark-by-benchmark contract:

- [results/validation/benchmark_summary.md](results/validation/benchmark_summary.md)
- [results/validation/validation_overview.md](results/validation/validation_overview.md)
- [results/validation/family_validation_overview.md](results/validation/family_validation_overview.md)
- [results/validation/family_deviation_audit.md](results/validation/family_deviation_audit.md)
- [results/validation/validated_envelope.md](results/validation/validated_envelope.md)
- [results/validation/matrix_family_coverage.md](results/validation/matrix_family_coverage.md)
- [results/validation/matrix_target_status.md](results/validation/matrix_target_status.md)
- [results/validation/matrix_promotion_contract.md](results/validation/matrix_promotion_contract.md)
- [results/validation/matrix_observable_closure_audit.md](results/validation/matrix_observable_closure_audit.md)
- [results/validation/matrix_experiment_intake_schema.md](results/validation/matrix_experiment_intake_schema.md)
- [results/validation/literature_learning_loop.md](results/validation/literature_learning_loop.md)
- [results/validation/family_promotion_state.md](results/validation/family_promotion_state.md)
