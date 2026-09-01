#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.matrix_recalibration import run_matrix_recalibration, write_refit_artifacts


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the matrix observable recalibration guard and emit an accept/reject artifact.")
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--write-changes", action="store_true", help="Persist accepted observable multipliers to data/lit/matrix_calibration_offsets.json.")
    parser.add_argument("--uncertainty-samples", type=int, default=60)
    parser.add_argument("--uncertainty-seed", type=int, default=0)
    args = parser.parse_args()

    report = run_matrix_recalibration(
        write_changes=bool(args.write_changes),
        uncertainty_samples=int(args.uncertainty_samples),
        uncertainty_seed=int(args.uncertainty_seed),
    )
    paths = write_refit_artifacts(report, output_dir=Path(args.output_dir))
    print(f"Decision: {report['decision']}")
    print(f"Wrote {paths['md']}")
    print(f"Wrote {paths['json']}")
    if report.get("persistence_path"):
        print(f"Persisted overrides: {report['persistence_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())