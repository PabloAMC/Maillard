# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Recommended](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Maillard is a scientist-facing computational framework for designing and screening meat-like flavour chemistry in plant-protein systems before committing to wet-lab campaigns.

The repository is built around one practical question: if a scientist changes precursors, matrix, or process severity, which flavour outcomes are worth testing next, what tradeoffs are likely, and how much of that prediction is actually supported by evidence?

## 1. The Problem

The Maillard reaction is one of the main reasons heated foods smell roasted, savory, sulfurous, nutty, or meaty. In a simple precursor system, the chemistry is already complex. In plant-based systems it becomes harder for four reasons:

- the reactive sugars and amino acids are not equally accessible inside the matrix;
- proteins, lipids, and water activity alter which compounds actually form and escape into headspace;
- severe processing steps such as extrusion create additional structure, transport, and damage effects;
- the strongest literature evidence is unevenly distributed across flavour families, matrices, and process states.

That means the real scientific problem is not just to predict Maillard chemistry. It is to make useful formulation decisions under uneven evidence.

In practice, a scientist usually wants to answer questions like these:

- Which formulation is the best next experiment for a meaty target?
- Which compounds are likely to improve at the same time as off-notes?
- Is a promising signal benchmark-backed, only directional, or still exploratory?
- What exact external experiment would move a blocked claim from plausible to decision-ready?

This repository exists to make those decisions explicit.

## 2. The Solution

The solution in this repository is not a single model. It is a layered decision system that separates what can be predicted mechanistically from what must be treated as matrix-limited or evidence-limited.

### 2.1 Core Scientific Strategy

The framework combines four layers.

| Layer | What it does | Why it matters |
| :--- | :--- | :--- |
| Chemistry layer | Encodes Maillard, Strecker, sulfur, carbonyl, and related reaction logic | Gives a mechanistic base for ranking likely flavour outcomes |
| Matrix/process layer | Adjusts interpretation using matrix accessibility, retention, denaturation, headspace, and process-state context | Prevents free-precursor chemistry from being over-transferred into plant matrices or extrusion |
| Evidence layer | Tracks whether each claim is benchmark-backed, transferred, directional, or blocked | Prevents quantitative overclaiming |
| Closure layer | States the smallest external package needed to advance a blocked decision | Turns “we need more data” into an exact experiment bundle |

### 2.2 What the Framework Is Trying To Achieve

The repository is trying to solve the plant-based flavour problem in a specific order:

1. Model the chemistry well enough to rank candidate flavour outcomes.
2. Correct that ranking for matrix and process effects that change real observables.
3. Show how much of the result is supported by benchmarked evidence.
4. Name the exact experiment needed when the evidence is not yet strong enough.

That ordering matters. The framework is useful because it does not confuse mechanism, calibration, and external closure.

### 2.3 The Three Trust Regimes

The same predicted compound means different things in different regimes.

| Regime | What dominates | What the tool is good for | What you should not claim |
| :--- | :--- | :--- | :--- |
| Free precursors | Core Maillard and Strecker chemistry | Quantitative screening and ranking | Direct transfer to real matrices without matrix evidence |
| Pea and soy matrices | Accessibility, retention, denaturation, headspace | Directional prioritization and bounded calibration | Promotion-grade quantitative closure without external mixed benchmarks |
| Extrusion-heavy systems | Thermal severity, residence time, structure, damage markers | Mechanistic staging and blocker discovery | Closed direct-damage claims without external reactive-lysine, furosine, and lysinoalanine data |

There is no single global accuracy claim. There are different decision postures for different scientific regimes.

## 3. The Tool

This repository turns that scientific strategy into a set of concrete workflows.

### 3.1 What the Tool Does

Maillard is:

- a formulation-screening engine for flavour targets and tradeoffs;
- a matrix-aware interpretation layer;
- an evidence and confidence surface for scientific claims;
- an experiment-design surface that states what evidence is still missing.

Maillard is not:

- a substitute for GC-MS, LC-MS, sensory work, or process trials;
- a universal quantitative predictor across every matrix and process condition;
- a justification for promotion-grade claims when the required external package is still open.

### 3.2 Who This Tool Is For

This repository is primarily for:

- scientists designing or prioritizing plant-based flavour experiments;
- computational researchers curating Maillard pathways, benchmarks, and evidence surfaces;
- teams that need a transparent bridge between mechanistic chemistry and experiment planning.

It is less useful if the goal is only to run a black-box aroma predictor without caring how strongly the result is supported.

### 3.3 What a Typical Run Produces

A standard workflow can produce:

- predicted volatile concentrations in ppb;
- flavour-target and off-note tradeoff summaries;
- matrix and process-aware interpretation panels;
- a claim posture such as benchmark-backed, directional, or exploratory;
- explicit artifacts that say what evidence would change the current decision.

Typical report shape:

```text
Formulation: Cysteine Enrichment (Pea Isolate, 120°C, 30 min)
─────────────────────────────────────────────────────────────
Compound                     Predicted    Trust         Note
2-Methyl-3-furanthiol        18.3 ppb     Benchmark     Core meaty marker
2-Furfurylthiol               9.1 ppb     Benchmark     Roasted / savory
2,5-Dimethylpyrazine          4.4 ppb     Directional   Nutty / roasted
Hexanal                      12.7 ppb     Benchmark     Off-note risk
─────────────────────────────────────────────────────────────
Matrix: Pea isolate
Decision posture: Directional ranking supported
```

### 3.4 How To Use The Tool To Implement The Solution

Use the repository by starting from the scientific task.

| If you need to do this part of the solution... | Use these workflows | What they implement |
| :--- | :--- | :--- |
| Rank candidate formulations | Forward prediction, formulation comparison, optimization | Chemistry layer plus practical ranking |
| Interpret matrix or process effects | Standard reports plus matrix and extrusion artifacts | Matrix/process layer |
| Check whether a claim is strong enough | Objective progress, family coverage, external evidence artifacts | Evidence layer |
| Decide what new experiment would help most | Primary or dual external package artifacts, extrusion package artifacts | Closure layer |
| Compare a new dataset to current support | Experiment intake and benchmark materialization workflows | Evidence update loop |

### 3.5 Quick Start

Docker is the recommended path because the repository is routinely validated there.

#### 1. Bootstrap the environment

```bash
./scripts/docker_maillard.sh bootstrap
```

#### 2. Run one prediction

```bash
./scripts/docker_maillard.sh run python scripts/run_pipeline.py \
  --protein-type pea_iso \
  --target meaty \
  --report \
  --output-dir results/first_run
```

#### 3. Compare or optimize candidates

```bash
./scripts/docker_maillard.sh run python scripts/compare_formulations.py \
  --baseline "Soy/Pea Base (Untreated)" \
  --candidate "Cysteine Enrichment (Basic)" \
  --protein-type pea_iso

./scripts/docker_maillard.sh run python scripts/optimize_formulation.py \
  --target-tag meaty \
  --minimize-tag beany \
  --protein-type pea_iso
```

#### 4. Inspect the evidence before making a strong claim

```bash
./scripts/docker_maillard.sh objective-progress
./scripts/docker_maillard.sh dft-coverage-map
./scripts/docker_maillard.sh pea-soy-external-evidence
./scripts/docker_maillard.sh pea-soy-mixed-external-package
```

#### 5. If extrusion matters, inspect the explicit blocker package

```bash
./scripts/docker_maillard.sh extrusion-external-closure
./scripts/docker_maillard.sh dha-lysinoalanine-external-package
```

### 3.6 How To Read The Output

The key question is not whether the number looks chemically plausible. The key question is what kind of claim that number can support.

| Output feature | What it means | How to use it |
| :--- | :--- | :--- |
| Benchmark-backed | The signal is tied to explicit literature-linked or equivalent benchmark support | Use for stronger decisions inside the validated regime |
| Directional | Relative ranking is more trustworthy than the absolute ppb | Use to prioritize experiments, not to claim closure |
| Exploratory | Mechanistically useful but not externally closed | Use as a hypothesis |
| Off-note or damage flag | Adverse markers or tradeoffs are rising | Use to reject or redesign the candidate |
| Matrix/process panel | Context changes interpretation of the same chemistry | Read this before trusting absolute concentration output |

Practical rule:

- in free chemistry, trust the ranking and often the scale;
- in pea and soy matrices, trust prioritization before absolute closure;
- in extrusion, trust the blocker analysis before the concentration estimate.

## 4. Why We Trust The Tool

Trust in this repository does not come from a single benchmark score. It comes from explicit evidence handling.

### 4.1 What the Trust Claim Is

The trust claim is narrow and concrete:

- the framework is strongest for free-precursor chemistry;
- it is useful for directional prioritization and bounded calibration in pea and soy matrices;
- it is useful for blocker discovery and experiment design in extrusion-heavy systems;
- it explicitly distinguishes what is benchmark-closed from what is still waiting on external data.

### 4.2 Where the Evidence Lives

The main evidence surfaces are these artifacts:

- [results/validation/objective_progress.md](results/validation/objective_progress.md): high-level objective closure without mixing internal calibration with external promotion.
- [results/validation/dft_coverage_map.md](results/validation/dft_coverage_map.md): chemistry-family coverage and where literature, xTB, DFT, and MLP do or do not matter.
- [results/validation/pea_soy_external_evidence.md](results/validation/pea_soy_external_evidence.md): what still blocks promotion-grade mixed claims in pea and soy.
- [results/validation/pea_soy_mixed_external_package.md](results/validation/pea_soy_mixed_external_package.md): the explicit dual external package needed to advance mixed pea/soy closure.
- [results/validation/primary_matrix_external_package.md](results/validation/primary_matrix_external_package.md): the prioritized single-lane package view for the next matrix move.
- [results/validation/extrusion_external_closure.md](results/validation/extrusion_external_closure.md): the explicit extrusion blocker surface.
- [results/validation/dha_lysinoalanine_external_package.md](results/validation/dha_lysinoalanine_external_package.md): the direct-damage external package needed to move DHA/LAL closure.

These artifacts are evidence maps. They are intentionally more useful than a single summary figure because they tell you exactly what is supported and exactly what is still missing.

![Benchmark-level quantitative validation](docs/assets/validation_overview.png)

This figure is the compact quantitative trust surface for the current benchmark set. The left panel shows predicted versus measured concentrations across benchmark compounds, and the right panel shows worst-case ratios per benchmark against explicit tolerance gates. Read it as evidence that the strongest regime in the repository is not just mechanistically plausible, but quantitatively anchored on the benchmark surface that the repo exposes.

### 4.3 What the Current Evidence Supports

At a high level, the repository currently supports these claims:

- free-precursor chemistry is the strongest and most quantitative part of the current system;
- matrix-aware work in pea and soy is useful for prioritization and bounded calibration;
- the repo now states the exact external packages needed for mixed pea/soy closure and for extrusion direct-damage closure;
- selective xTB and DFT are governance tools for decisive missing steps, not blanket replacements for literature-calibrated kinetics;
- MLP support is limited to offline acceleration and geometry staging, not barrier prediction.

### 4.4 What the Evidence Does Not Yet Support

The repository does not currently support these stronger claims:

- that pea or soy mixed meaty-positive lanes are already promotion-closed externally;
- that extrusion direct-damage claims are closed without reactive lysine, furosine, and lysinoalanine measurements;
- that all chemistry families are equally benchmarked;
- that DFT or MLP can replace literature-grounded calibration everywhere.

That boundary is part of the product, not an inconvenience. The repository is useful because it makes those limits explicit.

## 5. What To Do When The Tool Is Not Enough

If a result is scientifically interesting but still weak, the next step depends on why it is weak.

| If the limitation is... | The right next step is... |
| :--- | :--- |
| Missing benchmarked chemistry family | Extend literature and benchmark coverage before deeper refinement |
| Directional-only matrix support | Run the smallest external package that closes the matrix question |
| Extrusion direct-damage blocker | Measure reactive lysine, furosine, and lysinoalanine under controlled extrusion conditions |
| Uncertain mechanistic plausibility | Use xTB or selective DFT only for the decisive step, not the whole reaction graph |
| New experimental dataset | Intake and compare it before retuning the chemistry core |

The repository is designed to make those next steps concrete instead of vague.

![Objective progress and open blockers](docs/assets/objective_progress.png)

This figure is the fastest way to see when the tool is not enough yet. The left panel shows which tracked objectives are closed versus still blocked, and the right panel shows that one internal calibration lane is already in band while the external mixed-package and extrusion direct-damage lanes remain open. Read it as a decision map: if a lane is still blocked here, the next step is not more interpretation of the same prediction, but the specific external package named in the corresponding validation artifacts.

## 6. Repository Guide

If you are new to the codebase, these are the most useful starting points.

- Forward prediction and optimization scripts live in [scripts](scripts).
- Core implementation lives in [src](src).
- Validation artifacts are written to [results/validation](results/validation).
- Scientific and architectural documentation lives in [docs](docs).
- Tests live in [tests](tests).

Useful references:

- Scientific references: [docs/reference/SCIENTIFIC_REFERENCE.md](docs/reference/SCIENTIFIC_REFERENCE.md)
- Development documentation: [docs/development/README.md](docs/development/README.md)

## 7. Minimal Command Reference

### Prediction and optimization

```bash
./scripts/docker_maillard.sh run python scripts/run_pipeline.py --protein-type pea_iso --target meaty --report
./scripts/docker_maillard.sh run python scripts/compare_formulations.py --baseline BASE --candidate CANDIDATE --protein-type pea_iso
./scripts/docker_maillard.sh run python scripts/optimize_formulation.py --target-tag meaty --minimize-tag beany --protein-type pea_iso
```

### Trust and evidence

```bash
./scripts/docker_maillard.sh objective-progress
./scripts/docker_maillard.sh dft-coverage-map
./scripts/docker_maillard.sh pea-soy-external-evidence
./scripts/docker_maillard.sh pea-soy-mixed-external-package
./scripts/docker_maillard.sh primary-matrix-external-package
./scripts/docker_maillard.sh extrusion-external-closure
./scripts/docker_maillard.sh dha-lysinoalanine-external-package
```

### Experiment intake

```bash
./scripts/docker_maillard.sh compare-experiment path/to/intake.json
./scripts/docker_maillard.sh materialize-experiment-benchmark path/to/intake.json data/benchmarks/new_matrix_benchmark.json
```

### Reporting bundles

```bash
./scripts/docker_maillard.sh reporting-fast
./scripts/docker_maillard.sh scientific
```

## 8. Limits

- The strongest predictive regime today is still free-precursor chemistry.
- Pea and soy support is useful, but mixed promotion-grade closure still depends on external packages that are now specified rather than assumed.
- Extrusion is still a structured blocker-discovery regime, not a fully closed quantitative one.
- Selective DFT and xTB are refinement tools for decisive gaps, not full replacement engines.
- The repo is designed to support scientific decisions, not to hide uncertainty.

## 9. Glossary

| Term | Meaning in this repository |
| :--- | :--- |
| Maillard reaction | Reaction network initiated when reducing sugars and amino groups condense and rearrange under heat |
| Strecker degradation | Reaction between dicarbonyl species and amino acids that forms aroma aldehydes and heterocycle precursors |
| FAST screening | The lightweight runtime screening layer used for day-to-day formulation ranking |
| SMIRKS | Reaction transformation notation used to encode candidate chemistry |
| Headspace partitioning | Fraction of a volatile that actually escapes into the gas phase above the matrix |
| Benchmark-backed | Directly tied to a literature-linked measurement in a comparable context |
| Directional | Useful for relative prioritization even when absolute concentration is not closed |
| Exploratory | Mechanistically plausible but not externally closed enough for a stronger claim |
| xTB | Cheap semi-empirical quantum method used as a plausibility gate for selected steps |
| DFT | Higher-cost quantum chemistry used selectively for decisive missing barriers |
| MLP | Machine learning potential used here as an offline accelerator, not a truth engine |

## License

Apache License 2.0
