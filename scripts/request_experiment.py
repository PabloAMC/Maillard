#!/usr/bin/env python3
"""S22.2 CLI — request concrete experiments to close model gaps.

Usage:
  ./scripts/docker_maillard.sh run "python scripts/request_experiment.py \\
      --top 5 --protein-type soy --goal 'meaty aroma' --budget '3 lab days'"

For each top-N candidate, writes a pre-filled intake YAML under
`data/protocols/requested_*.yaml` and a human protocol under
`results/validation/experiment_requests/*.md`, plus an index file.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.experiment_request import build_requests, write_index
from src.experiment_value import (
    PREDICTION_UNCERTAINTY_PATH,
    load_prediction_payload,
    rank_experiments,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-path", default=str(PREDICTION_UNCERTAINTY_PATH))
    parser.add_argument("--top", type=int, default=5, help="Number of requests to emit.")
    parser.add_argument(
        "--protein-type",
        default=None,
        help="Substring filter against benchmark.protein_type (e.g. 'soy', 'pea').",
    )
    parser.add_argument(
        "--matrix",
        default=None,
        help=(
            "Comma-separated matrix families inferred from benchmark IDs "
            "(e.g. 'soy_iso,wheat_gluten'). Useful when benchmarks lack an "
            "explicit protein_type field."
        ),
    )
    parser.add_argument(
        "--goal",
        default=None,
        help="Free-text scientist intent recorded with the request (e.g. 'meaty aroma').",
    )
    parser.add_argument(
        "--budget",
        default=None,
        help="Free-text budget label (e.g. '3 lab days', '1 week SIDA slot').",
    )
    args = parser.parse_args()

    payload = load_prediction_payload(args.prediction_path)
    candidates = rank_experiments(payload)
    matrix_filter = (
        [m.strip() for m in args.matrix.split(",") if m.strip()]
        if args.matrix
        else None
    )
    requests = build_requests(
        candidates,
        top_n=args.top,
        protein_type=args.protein_type,
        matrix_filter=matrix_filter,
        goal=args.goal,
        budget_label=args.budget,
    )
    index_path = write_index(requests)

    print(f"Generated {len(requests)} request(s). Index: {index_path}")
    for r in requests:
        print(
            f"  #{r.candidate_rank} VoI={r.voi_score:.2f} {r.benchmark_id} :: {r.compound}"
        )
        print(f"      intake : {r.intake_yaml_path}")
        print(f"      protocol: {r.protocol_md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
