#!/usr/bin/env python
"""
THE CORE PANEL SCORECARD (retirement step B3).

Writes ``results/validation/core_panel_scores.{json,md}``: the kinetic core
scored deterministically on the union panel (trust loop + maillard_path
hold-outs + external matrix bundles) in the legacy harness's vocabulary --
per-benchmark contract status, strict-ready, evidence roles -- plus the 3x
band the hold-out was pre-registered against. See ``src/kinetic_core/scoring.py``.

This is a NEW artifact. The legacy ``benchmark_summary.json`` and the headline
guards are untouched; putting the two side by side and re-pinning is step B4.

Usage:
    python scripts/generators/generate_core_panel_scores.py \
        --output results/validation/core_panel_scores.json
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

DEFAULT_OUTPUT = data_paths.VALIDATION_DIR / "core_panel_scores.json"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", default=str(DEFAULT_OUTPUT),
        help="JSON path; the markdown table is written beside it with .md",
    )
    parser.add_argument("--pass-band", type=float, default=None,
                        help="fold band counted as 'within band' (default: the exam's 3x)")
    args = parser.parse_args(argv)

    from src.kinetic_core.scoring import PASS_BAND_LEVEL, score_panel, write_artifact

    started = time.perf_counter()
    payload = score_panel(pass_band=args.pass_band or PASS_BAND_LEVEL)
    payload["summary"]["wall_seconds"] = round(time.perf_counter() - started, 1)
    json_path, md_path = write_artifact(payload, Path(args.output))
    s = payload["summary"]
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(
        f"panel {s['panel_benchmark_count']} | scored {s['scored_benchmark_count']} | "
        f"rows {s['matched_compound_count']} (refused {s['refused_compound_count']}) | "
        f"within {payload['pass_band_level']:.0f}x {s['within_band_count']}/{s['matched_compound_count']} | "
        f"honest {s['honest_literature']['within_band']}/{s['honest_literature']['rows']} | "
        f"out-of-sample {s['out_of_sample']['within_band']}/{s['out_of_sample']['rows']} | "
        f"predictive passes {len(s['predictive_passes'])} | strict-ready {len(s['strict_ready'])} | "
        f"wall={s['wall_seconds']}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
