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

## The Real Problem

The most useful version of this tool is not a universal chemistry oracle. It is a scientist-facing decision system that can:

- rank candidate formulations before the wet lab
- explain why a prediction should be trusted or discounted
- separate benchmark-backed claims from directional extrapolation
- show which compounds are driven by free chemistry, matrix observability, or transferred priors

That means the central problem is broader than matching one benchmark table. We need a model that is useful across three regimes of scientific confidence.

## The Three Regimes

| Regime | What it represents | Current trust | What the tool should do well |
| --- | --- | --- | --- |
| Free precursors | Buffer-like systems where sugars and amino acids are directly accessible | High | Quantitative ranking and concentration-scale screening |
| Pea / soy matrices | Real protein matrices where accessibility, retention, and headspace matter | Moderate | Directional and near-quantitative prioritization where explicit anchors exist |
| Extrusion-heavy systems | Highly processed systems with strong structural and transport effects | Low | Exploratory hypothesis generation and experiment prioritization |

These regimes are not three unrelated products. They are three layers of increasing structural complexity.

- Free precursors test whether the chemistry core is right.
- Pea and soy matrices test whether observable reality is right.
- Extrusion-heavy systems test whether process-structured accessibility and transport are right.

## What Works Today

Current in-repo validation summary:

- 9 supported benchmarks are tracked in Docker-validated artifacts
- 4 benchmarks are strict-ready today, all inside the free-precursor envelope
- 9 matched compounds define the current authoritative quantitative proof surface
- the median matched-compound ratio in that proof surface is 1.118x

The current validation overview is shown below. The authoritative artifacts are generated locally with `./scripts/docker_maillard.sh validation-figures`.

![Validation Overview](docs/assets/validation_overview.png)

How to interpret this trust surface:

- Free-precursor benchmarks are the quantitative proof surface.
- Pea and soy matrix paths are useful for prioritization, not yet for broad release-grade claims.
- Extrusion-heavy systems remain exploratory until benchmarked directly.

## What Accuracy Depends On

The tool does not have a single accuracy mode. Its precision depends on which layer is doing most of the work.

| Layer | Main inputs | What it mainly controls |
| --- | --- | --- |
| Core FAST kinetics | Literature barriers, cached barrier records, calibrated Arrhenius heuristics | Pathway ordering and relative competitiveness |
| Observable projection | Proxy-to-observable mapping, volatility, retention, process-sensitive projection factors | Absolute ppb-scale fit |
| Matrix calibration and headspace release | Protein type, process state, compound-specific observable factors, transferred anchors | Whether matrix predictions stay physically plausible and useful |
| Offline refinement | xTB, selective DFT, and optional external ML-potential refinement written back as cached artifacts | Narrow mechanistic correction of decisive motif classes |

Two rules follow from this:

- We should not retune energy barriers just to make one benchmark concentration match.
- When ranking looks right but ppb scale is wrong, the preferred fix is usually in the observable or matrix-calibration layer.

The current validation contract therefore uses two scale checks:

- max ratio, to catch worst-case outliers
- mean absolute log-scale error, to catch broad multiplicative drift

## Scientific Architecture

The scientific stack is deliberately layered. Different tools serve different roles.

| Layer | Main tool or model | Role in the system | Main dependency |
| --- | --- | --- | --- |
| Reaction enumeration | SMIRKS rules and curated pathway families | Generate plausible Maillard and lipid-derived chemistry | Coverage and correctness of the encoded reaction families |
| Fast prediction core | FAST observable path with empirical or cached barriers | Daily screening, ranking, and benchmark-facing concentrations | Calibration quality of barrier tables and projection layers |
| Diagnostic thermochemistry | Cantera and thermodynamic gating | Diagnose whether a pathway family is physically plausible and whether gating would materially change benchmark behavior | Thermodynamic data and gating policy |
| Observable projection | Headspace, retention, matrix, and process-state surrogates | Convert pathway signal into what a scientist would actually measure or smell | Benchmark-anchored observable calibration |
| Offline mechanistic refinement | xTB first, selective DFT second | Refine decisive motif classes when cheap surrogates stop improving decisions | Availability of narrow, benchmark-relevant refinement targets |
| Optional ML acceleration | External state-of-the-art ML potentials, only offline | Accelerate conformer or local motif refinement after the refinement task is well-defined | Quality and relevance of the external model for the motif class of interest |

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
- [results/validation/validated_envelope.md](results/validation/validated_envelope.md)
