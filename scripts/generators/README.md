# Generators

This folder contains report and validation artifact generators (`generate_*.py`).

## Purpose

Each generator is a thin CLI wrapper that:

1. Builds a payload from `src/*` artifact builders.
2. Renders markdown and JSON outputs.
3. Writes outputs under `results/validation` (or a custom output dir when supported).

## Usage

Run a generator from the repository root, for example:

```bash
python scripts/generators/generate_family_ingestion_plan.py
python scripts/generators/generate_family_lane_validation.py --output-dir results/validation
python scripts/generators/generate_family_deviation_audit.py --output-dir results/validation
python scripts/generators/generate_family_validation_figures.py --output-dir results/validation --docs-asset-dir docs/assets
```

## Conventions

- Keep scientific logic in `src/`; keep these scripts orchestration-only.
- Prefer one script per artifact family for reproducibility and CI traceability.
- New scripts should follow existing patterns: argparse, deterministic payload build, markdown/json emit.

## The S3 sensitivity of the pH balance

`build_ph_balance_s3_sensitivity.py` (2026-11-07) fixes which shipped outputs depend on the pH balance that M1b part (i)
falsified (prereg secs. 20 - 21; CF-084, CF-134): every sulfur-lane output path at B0 and at the carried rivals'
representative points (B1 at phi_cys = 1, B13 at its best draw), read against omega 0.419 dex and the prereg's pH terms.
`--time-one` times one recipe and projects the run; `--workers 4` writes `results/validation/ph_balance_s3_sensitivity.{json,md}`,
the audit of the rule `src/kinetic_core/ph_balance_refusal.py` applies live on every balance-solved pot.
