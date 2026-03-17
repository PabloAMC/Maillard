# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Recommended](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Maillard is a computational formulation and validation framework for Maillard chemistry in alternative-protein systems. It is designed to help scientists understand what to try next in the wet lab, what the current model can support, and where the scientific limits still are.

## Start Here

If this is your first visit to the repository, read these in order:

1. [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md)
2. [docs/guides/SCIENTIFIC_RELIABILITY.md](docs/guides/SCIENTIFIC_RELIABILITY.md)
3. [docs/reference/COMMAND_REFERENCE.md](docs/reference/COMMAND_REFERENCE.md)
4. [docs/README.md](docs/README.md)

## What The Tool Can Do Today

- Explore Maillard reaction networks and rank likely flavor-active endpoints.
- Predict observable volatile outputs for free-precursor systems in a Docker-validated benchmark envelope.
- Score formulations against desirable, adverse, and safety targets.
- Optimize formulations over pH, temperature, and precursor choices.
- Run narrow matrix-aware pea and soy workflows with process-state and calibration metadata.
- Generate reproducible validation artifacts, including benchmark summaries, target snapshots, matrix assertions, readiness tables, branch deltas, and coverage gaps.

## What The Tool Cannot Honestly Claim Yet

- It is not yet a general-purpose predictor for intact protein foods.
- It does not yet have external meaty-positive matrix benchmarks for pea or soy.
- Matrix strict-gate promotion is still blocked pending wet-lab evidence.
- It should not be used as a sole source for absolute concentration claims in unbenchmarked systems.
- It does not replace wet-lab confirmation for extrusion-heavy, peptide-bound, or highly formulation-specific processes.

## Quick Trust Summary

- High trust: free-precursor ranking and comparative formulation work inside the current benchmark envelope.
- Moderate trust: matrix-aware pea and soy directional comparisons where the repo surfaces explicit caveats.
- Low trust: new protein families, severe processing states, or absolute claims without a nearby benchmark.

The detailed trust model, limits, and current benchmark table are in [docs/guides/SCIENTIFIC_RELIABILITY.md](docs/guides/SCIENTIFIC_RELIABILITY.md).

## Five-Minute Quick Start

On macOS, use Docker:

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validation-figures
```

This gives you two immediate outputs:

- a textual benchmark snapshot in [results/validation/benchmark_summary.md](results/validation/benchmark_summary.md)
- a graphical reliability and limitation snapshot in [results/validation/validation_overview.md](results/validation/validation_overview.md)

For platform-specific setup details, see [Installation.md](Installation.md).

## Typical Workflows

### 1. Check Whether The Current Model Is Trustworthy For Your Question

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
```

### 2. Inspect A Literature Benchmark

```bash
./scripts/docker_maillard.sh targets data/benchmarks/cys_glucose_150C_Farmer1999.json
```

### 3. Inspect Matrix Readiness And Current Matrix Limits

```bash
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh coverage-gaps
```

### 4. Compare Your Branch Against Main Before Claiming Progress

```bash
./scripts/docker_maillard.sh matrix-branch-deltas main
```

### 5. Run A Forward Prediction

```bash
python scripts/run_pipeline.py \
  --sugars ribose:0.5 \
  --amino-acids cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --protein-type pea_iso
```

### 6. Run A Formulation Optimizer

```bash
python scripts/optimize_formulation.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --target-tag meaty \
  --minimize-tag beany \
  --n-iterations 50
```

## What The Outputs Mean

- Sensory radar: a relative sensory profile, not a direct instrument reading.
- Predicted ppb: benchmark-facing observable concentration estimate where supported.
- Confidence metadata: whether the result is benchmark-backed, directional, or speculative.
- Projection metadata: which matrix and headspace corrections were applied.
- Validated-envelope warnings: when you leave the current scientific support zone.

## Architecture In One Table

| Layer | Role | Typical use |
| --- | --- | --- |
| Tier 0 | Reaction network generation | Enumerate plausible pathways |
| FAST | Observable concentration and flavor ranking | Default formulation work |
| Cantera | Diagnostic kinetics lane | Mechanism inspection and slower validation |
| MACE / xTB / DFT | Barrier refinement | Research-grade deep dives |

For a deeper architectural walkthrough, see [docs/architecture.md](docs/architecture.md) and [docs/reference/PROJECT_STRUCTURE.md](docs/reference/PROJECT_STRUCTURE.md).

## Validation And Reproducibility

The repository is Docker-first for scientific verification. The main reproducible lanes are:

- `./scripts/docker_maillard.sh core`
- `./scripts/docker_maillard.sh scientific-fast`
- `./scripts/docker_maillard.sh kinetics-validation`
- `./scripts/docker_maillard.sh scientific`
- `./scripts/docker_maillard.sh validation-figures`

If you need the strict free-precursor gate:

```bash
MAILLARD_STRICT_BENCHMARKS=1 ./scripts/docker_maillard.sh pytest tests/scientific/test_benchmarks.py -q
```

## Documentation Map

- [docs/README.md](docs/README.md): full document index by audience and purpose.
- [docs/guides/QUICKSTART.md](docs/guides/QUICKSTART.md): fastest path to a successful first run.
- [docs/guides/SCIENTIFIC_RELIABILITY.md](docs/guides/SCIENTIFIC_RELIABILITY.md): validated envelope, limits, and trust guidance.
- [docs/reference/COMMAND_REFERENCE.md](docs/reference/COMMAND_REFERENCE.md): Docker and validation command reference.
- [docs/reference/PROJECT_STRUCTURE.md](docs/reference/PROJECT_STRUCTURE.md): where code, data, tests, and artifacts live.
- [docs/VALIDATION_GUIDE.md](docs/VALIDATION_GUIDE.md): detailed benchmark contract and validation methodology.
- [docs/development/README.md](docs/development/README.md): backlog, lessons learned, and decision logs.
- [docs/research/README.md](docs/research/README.md): research notes, literature syntheses, and external reports.
