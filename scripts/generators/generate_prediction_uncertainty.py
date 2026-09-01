#!/usr/bin/env python3
"""S20.1 CLI — propagate barrier-family priors through the benchmark panel."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.uncertainty_propagation import propagate_benchmarks, write_artifact


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--basename", default="prediction_uncertainty")
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    payload = propagate_benchmarks(
        n_samples=args.n_samples,
        seed=args.seed,
        target_tag=args.target_tag,
    )
    paths = write_artifact(payload, output_dir=args.output_dir, basename=args.basename)
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
    summary = payload.get("summary", {})
    coverage = summary.get("ci_coverage_rate")
    if coverage is not None:
        print(
            f"90% CI coverage: {coverage * 100:.1f}% "
            f"({summary.get('ci_coverage_hits')}/{summary.get('matched_compound_count')})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
