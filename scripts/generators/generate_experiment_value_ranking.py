#!/usr/bin/env python3
"""S22.1 CLI — rank experiments by value-of-information."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiment_value import (
    PREDICTION_UNCERTAINTY_PATH,
    build_ranking_payload,
    load_prediction_payload,
    rank_experiments,
    write_artifact,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-path", default=str(PREDICTION_UNCERTAINTY_PATH))
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--basename", default="experiment_value_ranking")
    parser.add_argument("--top-n", type=int, default=None)
    args = parser.parse_args()

    source = Path(args.prediction_path)
    prediction_payload = load_prediction_payload(source)
    candidates = rank_experiments(prediction_payload, top_n=args.top_n)
    payload = build_ranking_payload(candidates, source_path=source)
    paths = write_artifact(payload, output_dir=args.output_dir, basename=args.basename)
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
    miss = payload["miss_count"]
    total = payload["candidate_count"]
    print(f"Out-of-CI rows: {miss}/{total}")
    if candidates:
        print("Top 5:")
        for c in candidates[:5]:
            print(
                f"  #{c.rank} VoI={c.voi_score:.2f} {c.benchmark_id} :: {c.compound} "
                f"(template={c.suggested_doe_template})"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
