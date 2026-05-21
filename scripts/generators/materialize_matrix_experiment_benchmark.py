#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_experiment_intake import materialize_matrix_experiment_benchmark


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_payload")
    parser.add_argument("output_benchmark")
    args = parser.parse_args()

    payload_path = Path(args.input_payload)
    output_path = Path(args.output_benchmark)
    benchmark = materialize_matrix_experiment_benchmark(payload_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(benchmark, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(benchmark, indent=2))
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())