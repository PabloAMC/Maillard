# Quick Start

This guide is for a scientist who wants to understand the tool quickly and get one trustworthy result out of it.

## Goal

In under 10 minutes, you should be able to:

- start the validated environment
- inspect the current benchmark envelope
- generate a graphical trust snapshot
- run one formulation prediction

## macOS Recommended Workflow

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validation-figures
```

Then open:

- [../../results/validation/benchmark_summary.md](../../results/validation/benchmark_summary.md)
- [../../results/validation/validation_overview.md](../../results/validation/validation_overview.md)

## First Prediction

```bash
python scripts/run_pipeline.py \
  --sugars ribose:0.5 \
  --amino-acids cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --protein-type pea_iso
```

What you should look for:

- `predicted_ppb`: observable concentration estimate
- `confidence_metadata`: how strongly the result is benchmark-backed
- `projection_metadata`: matrix/headspace corrections and calibration provenance
- any validated-envelope warnings that say the run is speculative

## If You Only Want To Check Scientific Trust

Run:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
```

Interpretation:

- summary: exact benchmark-by-benchmark status
- validated-envelope: plain-language boundary of what is supported
- validation-overview figure: one-page view of reliability and current limitations

## If You Care About Matrix Systems Specifically

Run:

```bash
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh coverage-gaps
```

These commands tell you three different things:

- what matrix data exists
- whether the current matrix benchmarks pass their own ranking assertions
- what still blocks external scientific assessment

## Before You Trust A Result

Read [SCIENTIFIC_RELIABILITY.md](SCIENTIFIC_RELIABILITY.md).

That document explains what the tool can support today, what remains directional only, and where wet-lab confirmation is still mandatory.
