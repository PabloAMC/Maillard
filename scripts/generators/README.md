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
```

## Conventions

- Keep scientific logic in `src/`; keep these scripts orchestration-only.
- Prefer one script per artifact family for reproducibility and CI traceability.
- New scripts should follow existing patterns: argparse, deterministic payload build, markdown/json emit.
