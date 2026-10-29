"""
C5 kind-1 check (independent review, 2026-10-29): the whole scored panel under the inline pH coupling against B60's
converged node fixed point, row by row. Ships the coupling only if every scored row agrees within the numerics
tolerance B60 itself used for its node-grid verification (0.03 dex) -- a number B60 pre-registered, not chosen here.
Writes results/validation/inline_coupling_panel_diff.json.

    python scripts/generators/generate_inline_coupling_panel_diff.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import sulfur  # noqa: E402
from src.kinetic_core.scoring import score_panel  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "inline_coupling_panel_diff.json"
TOL_DEX = 0.03  # B60 verification's V1 / quick-setting tolerance (kinetic_core_b60_verification.py)


def _rows(payload):
    return {(b["benchmark_id"], r["compound"]): r["predicted"] for b in payload["benchmarks"] for r in b["compounds"]
            if r.get("predicted")}


def main() -> int:
    sulfur.PH_COUPLING = "fixed_point_nodes"
    fixed = _rows(score_panel())
    sulfur.PH_COUPLING = "inline"
    inline = _rows(score_panel())
    diffs = {f"{k[0]} | {k[1]}": abs(math.log10(inline[k] / fixed[k])) for k in fixed if k in inline}
    worst = max(diffs, key=diffs.get)
    report = {"artifact": "inline_coupling_panel_diff", "tolerance_dex": TOL_DEX, "rows": len(diffs),
              "max_dex": diffs[worst], "worst_row": worst, "rows_over_tolerance": sorted(k for k, v in diffs.items() if v > TOL_DEX),
              "missing_in_one": sorted(set(map(str, fixed)) ^ set(map(str, inline))), "pass": diffs[worst] <= TOL_DEX}
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in ("rows", "max_dex", "worst_row", "pass")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
