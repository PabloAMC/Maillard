#!/usr/bin/env python
"""
THE KINETIC CORE'S FIT-ROW DECLARATIONS, in the fit-target index's vocabulary (2026-09-03).

Reads the sulfur fit's row table (``generate_kinetic_core_b2_3_fit.FIT_ROWS`` plus the
four rows B8 installs) -- each level row now carries ``benchmark_id`` and
``benchmark_compound`` next to its source anchor -- and the B8 report's free set, and
writes ``results/validation/kinetic_core_<wave>_fit_targets.json`` (default wave b9):

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

OUT = None  # kinetic_core_<wave>_fit_targets.json


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=None)
    parser.add_argument("--wave", default="b9", help="the sulfur wave whose rows are declared (b8 or b9)")
    args = parser.parse_args(argv)

    import importlib

    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402  (installs nothing)
    full_rows = list(B23.FIT_ROWS)  # before a wave import removes level rows (B9) and with them the pot -> bundle links
    B8 = importlib.import_module(f"generate_kinetic_core_{args.wave}_fit")  # installs / removes rows in B23
    if not Path(B8.OUT_FIT_REPORT).exists() and args.wave != "b8":
        B8 = importlib.import_module("generate_kinetic_core_b8_fit")
    wave = B8.WAVE.lower() if hasattr(B8, "WAVE") else "b8"
    out = Path(args.output) if args.output else data_paths.VALIDATION_DIR / f"kinetic_core_{wave}_fit_targets.json"

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
        if row.get("benchmark_id_b"):
            # B10: a within-study FOLD reads two bundles' values (numerator and
            # denominator); both leave the strict out-of-sample count.
            rows.append(
                {
                    "id": row["id"] + "__denominator",
                    "benchmark_id": row["benchmark_id_b"],
                    "compound": row["benchmark_compound"],
                    "kind": row["kind"],
                    "anchor": row.get("anchor", ""),
                }
            )
    # 2026-10-28 (independent review; north-star decisions): a fit row that reads NO bundle compound but reads
    # the SAME POT as a bundle -- Hofmann's in-situ norfuraneol and furfural, MFT's and FFT's precursors in the
    # ribose pH 5 pot -- still informs that bundle's predictions. Link it at the pot level, through the pot
    # (`system`) -> bundle map that the full row table declares, removed level rows included.
    pot_bundle: Dict[str, set] = {}
    for row in full_rows + list(B23.ACTIVE_FIT_ROWS):
        if row.get("benchmark_id") and row.get("system"):
            pot_bundle.setdefault(str(row["system"]), set()).add(str(row["benchmark_id"]))
    from src.kinetic_core.panel import panel_bundles

    panel_ids = {Path(p).stem for p, _ in panel_bundles()[0]}

    def _resolve(bundle: str) -> str:
        """A row table written before a bundle moved to the hold-out panel names it by its old id."""
        return bundle if bundle in panel_ids or f"mp_holdout_{bundle}" not in panel_ids else f"mp_holdout_{bundle}"

    for row in B23.ACTIVE_FIT_ROWS:
        if row.get("benchmark_id") or not row.get("system"):
            continue
        for bundle in sorted({_resolve(b) for b in pot_bundle.get(str(row["system"]), ())}):
            rows.append(
                {
                    "id": row["id"],
                    "benchmark_id": bundle,
                    "compound": f"[same pot] {row.get('species', '')}",
                    "kind": row["kind"],
                    "anchor": row.get("anchor", ""),
                    "pot_level": True,
                }
            )
    payload = {
        "artifact": "kinetic_core_fit_targets",
        "lane": "sulfur",
        "fit_report": data_paths.rel(Path(B8.OUT_FIT_REPORT)),
        "generated_by": "scripts/generators/generate_core_fit_targets.py",
        "declaration": (
            "Level rows of the sulfur objective that are ALSO scored panel bundles; (B10) the "
            "within-study FOLD rows that read two bundles' values; and (2026-10-28) POT-LEVEL rows, "
            "which read an intermediate of the same pot a bundle scores (pot_level: true), so every row "
            "of that bundle leaves the out-of-sample count. Rows that read no bundle's pot are not listed. Leverage "
            f"is {n_free} free parameters over {n_rows} rows, below the per-row-recovery "
            "threshold: these bundles stay in the coverage counts, annotated in_core_fit."
        ),
        "fit_target_ids": sorted({r["benchmark_id"] for r in rows}),
        "fit_leverage": {"free_parameters": n_free, "fitted_rows": n_rows},
        "rows": rows,
    }
    payload["wave"] = wave
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}: {len(rows)} level rows over {len(payload['fit_target_ids'])} bundles; leverage {n_free}/{n_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
