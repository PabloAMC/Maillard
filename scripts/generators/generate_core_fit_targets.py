#!/usr/bin/env python
"""
THE KINETIC CORE'S FIT-ROW DECLARATIONS, in the fit-target index's vocabulary (2026-09-03).

Reads the sulfur fit's row table (``generate_kinetic_core_b2_3_fit.FIT_ROWS`` plus the
four rows B8 installs) -- each level row now carries ``benchmark_id`` and
``benchmark_compound`` next to its source anchor -- and the B8 report's free set, and
writes ``results/validation/kinetic_core_b8_fit_targets.json``:

    fit_target_ids   the panel bundles whose level rows the fit read
    fit_leverage     {free_parameters, fitted_rows}  (23 / 62 -> global_low_leverage)
    rows             [{id, benchmark_id, compound, kind, anchor}, ...]

``src/fit_target_index.py`` and ``scripts/ci/fit_target_gate.py`` glob ``*_fit_targets.json``,
so the index, the hold-out guard (check 4: no fit record names a hold-out bundle) and
``src/kinetic_core/fit_targets.py`` all read THIS file. The hand-typed table that stood
in for it from B3 to this step is gone. Nothing here fits anything.

Usage:
    python scripts/generators/generate_core_fit_targets.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import data_paths  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b8_fit_targets.json"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(OUT))
    args = parser.parse_args(argv)

    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402  (installs nothing)
    import generate_kinetic_core_b8_fit as B8  # noqa: E402  (installs B8's four rows into B23)

    report = json.loads(Path(B8.OUT_FIT_REPORT).read_text(encoding="utf-8"))
    n_free = int(report["free_set"]["n_free"])
    n_rows = len(B23.ACTIVE_FIT_ROWS)
    rows: List[Dict[str, Any]] = []
    for row in B23.ACTIVE_FIT_ROWS:
        if row.get("benchmark_id"):
            rows.append(
                {
                    "id": row["id"],
                    "benchmark_id": row["benchmark_id"],
                    "compound": row["benchmark_compound"],
                    "kind": row["kind"],
                    "anchor": row.get("anchor", ""),
                }
            )
    payload = {
        "artifact": "kinetic_core_fit_targets",
        "lane": "sulfur",
        "fit_report": data_paths.rel(Path(B8.OUT_FIT_REPORT)),
        "generated_by": "scripts/generators/generate_core_fit_targets.py",
        "declaration": (
            "Level rows of the sulfur objective that are ALSO scored panel bundles. Ratio, share "
            "and conversion rows constrain nothing a bundle scores and are not listed. Leverage "
            f"is {n_free} free parameters over {n_rows} rows, below the per-row-recovery "
            "threshold: these bundles stay in the coverage counts, annotated in_core_fit."
        ),
        "fit_target_ids": sorted({r["benchmark_id"] for r in rows}),
        "fit_leverage": {"free_parameters": n_free, "fitted_rows": n_rows},
        "rows": rows,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}: {len(rows)} level rows over {len(payload['fit_target_ids'])} bundles; leverage {n_free}/{n_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
