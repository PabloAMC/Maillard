# Maillard Framework Philosophy

## The Goal and the Problem

The framework operates as a rigorous forward-prediction and optimization engine, rather than a static lookup database. It takes initial formulation conditions (precursors, protein matrix type, physical state) and deterministic process variables (pH, temperature, time), and computationally resolves the competitive kinetic pathways of the Maillard reaction. 

The most useful version of this tool is a scientist-facing decision system that can:

- rank candidate formulations prior to wet-lab campaigns
- optimize targeted aromatic profiles (e.g., maximizing meaty thiols while minimizing beany off-notes) by computationally navigating the complex reaction network
- separate benchmark-backed certainty from directional extrapolation
- isolate which volatile markers are driven by free chemistry versus matrix retention

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

Why all families are not quantitative today:
- A quantitative family needs at least one executable benchmark with direct measured compounds or markers that can be matched against predictions.
- Several integrated families are upstream modifiers, guardrails, or support lanes rather than endpoint volatile panels, so they change the prediction context without yet having their own matched numeric compound surface.
- Those families are still tracked explicitly in validation reports instead of being over-claimed as quantitative.

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

The validation contract therefore uses two scale checks:
- max ratio, to catch worst-case outliers
- mean absolute log-scale error, to catch broad multiplicative drift

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

## Literature Integration: Runtime vs. Benchmark

Literature serves as the strict deterministic control system governing predictive accuracy, rather than mere citations. Every paper is rigorously evaluated through an 8-point Systematic Literature Review (SLR) checklist. Literature data maps into two mutually exclusive functions:

### 1. Runtime Integration (Parametrization)
Runtime data drives the underlying scientific engine. These papers supply the chemical constants, activation configurations, and thermodynamic limits used to compute the paths.
- **Data type:** Rate constants, activation energies ($E_a$), matrix retention multipliers.
- **Example:** A structural study providing the dynamic, reversible binding affinity ($K_a$) of hexanal to a soy glycinin isolate. This is parameterized into the retention module to accurately constrain expected headspace release levels.

### 2. Benchmark Integration (Validation)
Benchmark data establishes the isolated "ground truth" used exclusively for quality control and structural validation. To strictly prevent overfitting ("circular validation"), benchmark targets are fully independent of runtime calibration parameters.
- **Data type:** Absolute end-state product concentrations (e.g., ppb) measured under definitively constrained starting conditions.
- **Example:** A study producing 45 ppb of 2-methyl-3-furanthiol isolated from 10 mM ribose and 10 mM cysteine in a pea protein isolate under 120°C for 30 minutes.
