#!/usr/bin/env python3
"""
scripts/generators/generate_trunk_calibration_holdout_prepost.py

Scores Wave S3's trunk rate calibration against the PRE-REGISTERED hold-out
baseline, and writes `results/validation/maillard_path_holdout_S3_prepost.{json,md}`.

THE BASELINE IS A FILE, NOT A RE-RUN.
-------------------------------------
`results/validation/maillard_path_holdout_frozen_predictions.json` was
generated at git HEAD `12f43dd`, before any rate-calibration wave saw these
twelve bundles. This script READS that file and never regenerates it. Comparing
a wave against a baseline the same wave produced is the circularity the audit
campaign exists to stop, and it is why the artifact is committed and
append-only.

THREE COLUMNS, AND THE MIDDLE ONE IS THE HONEST ANSWER
------------------------------------------------------
  * FROZEN    -- the committed pre-registration at 12f43dd.
  * AS SHIPPED -- the hold-out re-run at the current tree. Wave S3 changes no
    barrier, no projection constant and no observability factor, so this column
    is EXPECTED to reproduce the frozen one exactly. Any drift is a finding
    about something else and is reported as such rather than absorbed.
  * COUNTERFACTUAL -- the hold-out re-run with the derived screening-lane
    barriers from `results/validation/trunk_rate_calibration_refit.json`
    applied in memory. This is NOT shipped and is NOT written to
    `src/barrier_constants.py`; it exists so the owner can see the size of the
    decision they are being asked to make, measured rather than argued.

Usage:
    python scripts/generators/generate_trunk_calibration_holdout_prepost.py
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths

VALIDATION = data_paths.VALIDATION_DIR
FROZEN = VALIDATION / "maillard_path_holdout_frozen_predictions.json"
CALIBRATION = VALIDATION / "trunk_rate_calibration_refit.json"
FROZEN_COMMIT = "12f43dd"


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def flatten(payload: Dict[str, Any]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """Key every scored target by (benchmark_id, compound)."""
    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for bench in payload.get("benchmarks", []):
        for compound in bench.get("compounds", []):
            out[(bench["benchmark_id"], compound["compound"])] = compound
    return out


def flatten_series(payload: Dict[str, Any]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for series in payload.get("series", []):
        for entry in series.get("entries", []):
            out[(series["series_id"], entry["compound"])] = entry
    return out


def run_holdout(barrier_overrides: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """
    Re-run the twelve bundles, optionally with FAST_BARRIERS patched in memory.

    The patch is applied to the live dict and reverted in a `finally`, so a
    failure cannot leave the process holding modified barriers.
    """
    from src import barrier_constants
    from scripts.generators.generate_maillard_path_holdout_frozen_predictions import (
        build_frozen_predictions,
    )

    saved: Dict[str, Any] = {}
    try:
        if barrier_overrides:
            for family, value in barrier_overrides.items():
                if family not in barrier_constants.FAST_BARRIERS:
                    raise KeyError(f"unknown FAST_BARRIERS family {family!r}")
                saved[family] = barrier_constants.FAST_BARRIERS[family]
                barrier_constants.FAST_BARRIERS[family] = (
                    float(value), "WAVE S3 COUNTERFACTUAL -- in memory only, never shipped"
                )
        return build_frozen_predictions(target_tag="meaty")
    finally:
        for family, original in saved.items():
            barrier_constants.FAST_BARRIERS[family] = original


def derived_overrides() -> Tuple[Dict[str, float], List[Dict[str, Any]]]:
    calibration = load_json(CALIBRATION)
    propagation = calibration["screening_lane_propagation"]
    overrides: Dict[str, float] = {}
    rows: List[Dict[str, Any]] = []
    for row in propagation["rows"]:
        if row.get("derived_barrier_kcal_mol") is None:
            rows.append(row)
            continue
        overrides[row["family"]] = float(row["derived_barrier_kcal_mol"])
        rows.append(row)
    return overrides, rows


def summarise(payload: Dict[str, Any]) -> Dict[str, Any]:
    return dict(payload.get("summary", {}))


def compare(frozen: Dict[str, Any], other: Dict[str, Any], label: str) -> Dict[str, Any]:
    a, b = flatten(frozen), flatten(other)
    rows = []
    moved = 0
    improved = 0
    worsened = 0
    for key in sorted(set(a) | set(b)):
        left, right = a.get(key), b.get(key)
        lf = (left or {}).get("fold_error")
        rf = (right or {}).get("fold_error")
        lp = (left or {}).get("predicted_native_unit")
        rp = (right or {}).get("predicted_native_unit")
        changed = not (
            lp is not None and rp is not None
            and (lp == rp or (rp != 0 and abs(lp - rp) <= 1e-9 * max(abs(lp), abs(rp))))
        )
        if lp is None and rp is None:
            changed = False
        if changed:
            moved += 1
            if lf is not None and rf is not None:
                if rf < lf:
                    improved += 1
                elif rf > lf:
                    worsened += 1
        rows.append({
            "benchmark_id": key[0],
            "compound": key[1],
            "target_value": (left or right or {}).get("target_value"),
            "target_unit": (left or right or {}).get("target_unit"),
            "frozen_predicted": lp,
            f"{label}_predicted": rp,
            "frozen_fold_error": lf,
            f"{label}_fold_error": rf,
            "changed": changed,
            "scoring": (right or left or {}).get("scoring"),
        })

    sa, sb = flatten_series(frozen), flatten_series(other)
    series_rows = []
    for key in sorted(set(sa) | set(sb)):
        left, right = sa.get(key, {}), sb.get(key, {})
        series_rows.append({
            "series_id": key[0],
            "compound": key[1],
            "measured_ratio": left.get("measured_ratio_high_over_low"),
            "frozen_predicted_ratio": left.get("predicted_ratio_high_over_low"),
            f"{label}_predicted_ratio": right.get("predicted_ratio_high_over_low"),
            "frozen_direction_correct": left.get("direction_correct"),
            f"{label}_direction_correct": right.get("direction_correct"),
        })

    return {
        "targets_changed": moved,
        "targets_total": len(rows),
        "improved": improved,
        "worsened": worsened,
        "frozen_summary": summarise(frozen),
        f"{label}_summary": summarise(other),
        "rows": rows,
        "series": series_rows,
    }


def build() -> Dict[str, Any]:
    frozen = load_json(FROZEN)
    if not frozen.get("git", {}).get("short", "").startswith(FROZEN_COMMIT):
        raise SystemExit(
            f"the frozen artifact does not name commit {FROZEN_COMMIT}; refusing to score "
            f"against a baseline that is not the pre-registration"
        )

    as_shipped = run_holdout(None)
    overrides, propagation_rows = derived_overrides()
    counterfactual = run_holdout(overrides)

    return {
        "artifact": "maillard_path_holdout_S3_prepost",
        "wave": "S3",
        "generated_on": date.today().isoformat(),
        "baseline": {
            "file": data_paths.rel(data_paths.VALIDATION_DIR / "maillard_path_holdout_frozen_predictions.json"),
            "git_commit": frozen["git"]["short"],
            "note": "READ, never regenerated. The pre-registration is the file, not a re-run.",
        },
        "as_shipped": compare(frozen, as_shipped, "shipped"),
        "counterfactual_derived_barriers": {
            "applied_in_memory_only": True,
            "overrides_kcal_mol": overrides,
            "propagation_rows": propagation_rows,
            "comparison": compare(frozen, counterfactual, "counterfactual"),
        },
    }


def render(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    add = lines.append
    add("# Wave S3 — hold-out pre/post against the pre-registered baseline")
    add("")
    add(f"Generated {payload['generated_on']}. Baseline: "
        f"`{payload['baseline']['file']}` at commit `{payload['baseline']['git_commit']}`. "
        f"{payload['baseline']['note']}")
    add("")

    for section, title, label in (
        ("as_shipped", "AS SHIPPED — what this wave actually changes", "shipped"),
        (None, "COUNTERFACTUAL — what applying the derived barriers WOULD do "
               "(measured; not shipped)", "counterfactual"),
    ):
        if section is None:
            block = payload["counterfactual_derived_barriers"]["comparison"]
            add(f"## {title}")
            add("")
            add("Barrier overrides applied in memory only:")
            add("")
            for family, value in payload["counterfactual_derived_barriers"][
                    "overrides_kcal_mol"].items():
                add(f"* `{family}` → {value:.2f} kcal/mol")
            add("")
        else:
            block = payload[section]
            add(f"## {title}")
            add("")

        frozen_summary = block["frozen_summary"]
        other_summary = block[f"{label}_summary"]
        add(f"**{block['targets_changed']}/{block['targets_total']} targets moved** "
            f"({block['improved']} improved, {block['worsened']} worsened).")
        add("")
        add("| measure | frozen | " + label + " |")
        add("|---|---:|---:|")
        for key in ("median_fold_error", "median_abs_log10_error", "within_10x",
                    "worst_fold_error", "best_fold_error", "structural_zero_count",
                    "ordinal_pairs_correct", "series_directions_correct"):
            left = frozen_summary.get(key)
            right = other_summary.get(key)
            fmt = (lambda v: f"{v:.4f}" if isinstance(v, float) else str(v))
            add(f"| `{key}` | {fmt(left)} | {fmt(right)} |")
        add("")
        add("### Series directions (the three Wave U structural errors live here)")
        add("")
        add("| series | compound | measured ratio | frozen predicted | "
            + label + " predicted | frozen ok? | " + label + " ok? |")
        add("|---|---|---:|---:|---:|---|---|")
        for row in block["series"]:
            def num(v: Any) -> str:
                return f"{v:.4g}" if isinstance(v, (int, float)) else "—"
            add(f"| {row['series_id']} | {row['compound']} | "
                f"{num(row['measured_ratio'])} | "
                f"{num(row['frozen_predicted_ratio'])} | "
                f"{num(row[f'{label}_predicted_ratio'])} | "
                f"{row['frozen_direction_correct']} | "
                f"{row[f'{label}_direction_correct']} |")
        add("")
        changed = [r for r in block["rows"] if r["changed"]]
        if changed:
            add("### Targets that moved")
            add("")
            add("| benchmark | compound | target | frozen fold | " + label + " fold |")
            add("|---|---|---:|---:|---:|")
            for row in changed:
                def num(v: Any) -> str:
                    return f"{v:.4g}" if isinstance(v, (int, float)) else "—"
                add(f"| {row['benchmark_id']} | {row['compound']} | "
                    f"{num(row['target_value'])} {row['target_unit']} | "
                    f"{num(row['frozen_fold_error'])} | "
                    f"{num(row[f'{label}_fold_error'])} |")
        else:
            add("**No target moved. Every prediction is bit-identical to the "
                "pre-registration.**")
        add("")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Score Wave S3's trunk rate calibration against the frozen hold-out "
            "baseline in three columns (frozen / as shipped / counterfactual); "
            "writes "
            "results/validation/maillard_path_holdout_S3_prepost.{json,md}."
        )
    )
    parser.parse_args(argv)

    payload = build()
    (VALIDATION / "maillard_path_holdout_S3_prepost.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (VALIDATION / "maillard_path_holdout_S3_prepost.md").write_text(
        render(payload), encoding="utf-8")
    shipped = payload["as_shipped"]
    counter = payload["counterfactual_derived_barriers"]["comparison"]
    print(f"as shipped: {shipped['targets_changed']}/{shipped['targets_total']} targets moved")
    print(f"counterfactual: {counter['targets_changed']}/{counter['targets_total']} "
          f"targets moved ({counter['improved']} better, {counter['worsened']} worse)")
    print(f"frozen median  {shipped['frozen_summary'].get('median_fold_error')}")
    print(f"shipped median {shipped['shipped_summary'].get('median_fold_error')}")
    print(f"counter median {counter['counterfactual_summary'].get('median_fold_error')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
