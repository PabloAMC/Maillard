#!/usr/bin/env python
"""
THE CORE MONTE-CARLO ENVELOPE (retirement step B2).

Writes ``results/validation/core_prediction_uncertainty.{json,md}``: per
(benchmark, compound) P5/P50/P95 predicted values from re-integrating the
KINETIC CORE under draws of its own fitted coordinates and declared bands,
over the union panel (trust loop + maillard_path hold-outs + external matrix
bundles). See ``src/kinetic_core/uncertainty.py`` for what is sampled and
from where -- and for what is NOT (the sulfur lane's fit report carries no
uncertainty, and none is invented for it).

This is a NEW artifact. The legacy ``prediction_uncertainty.json`` is left
untouched; the re-pin of README / model card / the headline guards onto this
file is a later step (B4).

Usage:
    python scripts/generators/generate_core_prediction_uncertainty.py \
        --n-samples 200 --seed 0 --output results/validation/core_prediction_uncertainty.json
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

DEFAULT_OUTPUT = data_paths.VALIDATION_DIR / "core_prediction_uncertainty.json"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--output", default=str(DEFAULT_OUTPUT),
        help="JSON path; the markdown table is written beside it with .md",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="worker processes over draws (result is independent of this)",
    )
    args = parser.parse_args(argv)

    from src.kinetic_core.uncertainty import propagate_panel, write_artifact

    started = time.perf_counter()
    payload = propagate_panel(
        n_samples=args.n_samples, seed=args.seed, workers=args.workers
    )
    wall = time.perf_counter() - started
    payload["summary"]["wall_seconds"] = round(wall, 1)
    payload["summary"]["workers"] = int(args.workers)
    json_path, md_path = write_artifact(payload, Path(args.output))

    s = payload["summary"]
    lit = s["honest_literature_coverage"]
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(
        f"n={s['n_samples']} seed={s['seed']} workers={args.workers} wall={wall:.1f}s | "
        f"benchmarks {s['benchmark_count']}/{s['panel_benchmark_count']}, rows "
        f"{s['matched_compound_count']} (refused {s['refused_compound_count']}) | "
        f"honest literature coverage {lit['hits']}/{lit['total']} "
        f"median width {lit['median_ci_width_log10']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
