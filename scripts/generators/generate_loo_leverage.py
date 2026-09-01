#!/usr/bin/env python3
"""S20.2 CLI — leave-one-benchmark-out leverage diagnostics."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.cross_validation import (
    PREDICTION_UNCERTAINTY_PATH,
    compute_leverage,
    load_prediction_payload,
    write_artifact,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-path", default=str(PREDICTION_UNCERTAINTY_PATH))
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--basename", default="loo_leverage")
    args = parser.parse_args()

    payload = load_prediction_payload(args.prediction_path)
    report = compute_leverage(payload)
    paths = write_artifact(report, output_dir=args.output_dir, basename=args.basename)
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
    panel = report.get("panel", {})
    print(
        f"Panel coverage: {panel.get('inside_ci_count')}/{panel.get('matched_compound_count')}"
        f" ({(panel.get('coverage') or 0.0) * 100:.1f}%)"
    )
    print("Top drag (largest negative leverage ⇒ best LOO target):")
    for entry in sorted(
        report.get("benchmarks", []),
        key=lambda r: (r.get("leverage") if r.get("leverage") is not None else 0.0),
    )[:5]:
        lev = entry.get("leverage")
        lev_str = f"{lev:+.3f}" if lev is not None else "n/a"
        print(
            f"  {lev_str}  {entry.get('benchmark_id')}"
            f"  (mean|Δlog₁₀|={entry.get('mean_abs_log10_residual')})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
