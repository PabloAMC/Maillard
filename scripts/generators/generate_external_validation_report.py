#!/usr/bin/env python3
"""S20.5 Lane A.3 — score the isolated hold-out bundles against the frozen calibration."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.external_validation import build_external_validation_report, write_external_validation_artifact


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--basename", default="external_validation_report")
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    payload = build_external_validation_report(
        n_samples=args.n_samples,
        seed=args.seed,
        target_tag=args.target_tag,
    )
    paths = write_external_validation_artifact(
        payload,
        output_dir=Path(args.output_dir),
        basename=args.basename,
    )
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
    summary = payload.get("summary", {})
    coverage = summary.get("ci_coverage_rate")
    coverage_str = f"{coverage * 100:.1f}%" if coverage is not None else "n/a"
    median_accuracy = summary.get("median_accuracy_fold")
    median_str = f"{float(median_accuracy):.2f}x" if median_accuracy is not None else "n/a"
    print(
        f"External hold-out coverage: {summary.get('ci_coverage_hits', 0)}/{summary.get('matched_compound_count', 0)} "
        f"({coverage_str}) | median accuracy {median_str}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())