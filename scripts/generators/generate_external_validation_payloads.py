#!/usr/bin/env python3
"""S20.5 Lane A.2 — materialize external-validation-only hold-out payloads."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.external_validation import build_holdout_bundles, write_holdout_bundles


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--protocol-dir",
        default="data/protocols/external_validation",
        help="Destination directory for synthesized intake YAML files.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="data/benchmarks/external_validation",
        help="Destination directory for isolated hold-out benchmark JSON files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Rewrite bundles that already exist. By default existing bundles are frozen "
            "evidence and are skipped; with --overwrite, keys the committed file has and the "
            "regenerated payload lacks are still carried forward."
        ),
    )
    args = parser.parse_args()

    bundles = build_holdout_bundles()
    written = write_holdout_bundles(
        bundles,
        protocol_dir=Path(args.protocol_dir),
        benchmark_dir=Path(args.benchmark_dir),
        overwrite=args.overwrite,
    )
    payload = {
        "bundle_count": len(bundles),
        "matched_compound_count": sum(bundle.matched_compound_count() for bundle in bundles),
        "protocols": [str(path) for path in written["protocols"]],
        "benchmarks": [str(path) for path in written["benchmarks"]],
        "skipped_existing": [str(path) for path in written["skipped"]],
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())