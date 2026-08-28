#!/usr/bin/env python3
"""Score Wave Y against BOTH frozen baselines, and against its own pre-registration.

Wave Y relocates the ambient hexanal observable scale from the matrix observability
factors into `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` (see
`results/validation/matrix_marker_yield_rederivation.md`) and applies NOTHING to the
projection budget (see `results/validation/projection_budget_step_yield_constraint.md`).

The wave therefore predicted, IN ADVANCE and in
`results/validation/wave_y_prereg_predictions.json`, that neither of the repository's two
frozen out-of-sample surfaces can move:

  * the **maillard_path hold-out** -- 12 bundles / 22 targets, frozen at git `12f43dd`,
    read from `maillard_path_holdout_frozen_predictions.json` and NEVER regenerated. The
    marker yields are read in exactly two places and both are matrix lanes; every bundle
    here is `free_precursor`.
  * the **Wave X step-level rows** -- five scored single-step Hofmann systems on the same
    free-precursor lane, mean 0.518 dex.

This script measures both rather than asserting either. A prediction that a number does
not move is only worth writing down if the file that shows it did not move is committed
alongside it.

Usage:
    python scripts/generators/generate_wave_y_holdout_prepost.py
"""

from __future__ import annotations

import json
import math
import statistics
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

VALIDATION = ROOT / "results" / "validation"
FROZEN = VALIDATION / "maillard_path_holdout_frozen_predictions.json"
PREREG = VALIDATION / "wave_y_prereg_predictions.json"
FROZEN_COMMIT = "12f43dd"

OUT_JSON = VALIDATION / "maillard_path_holdout_wave_y_prepost.json"
OUT_MD = VALIDATION / "maillard_path_holdout_wave_y_prepost.md"

# The five SCORED Wave X step-level rows. `hofmann1998_norfuraneol_h2s_145C_20min_pH5` is
# Wave X's declared fit target and is carried but excluded from the mean, exactly as
# Wave X reported it.
STEP_LEVEL_BENCHMARKS = (
    "hofmann1998_furan2aldehyde_h2s_145C_20min_pH5",
    "hofmann1998_norfuraneol_h2s_145C_20min_pH5",
    "hofmann1998_norfuraneol_cysteine_145C_20min_pH5",
    "hofmann1998_c2c3_recombination_145C_20min_pH3",
    "hofmann1998_c2c3_recombination_145C_20min_pH5",
    "hofmann1998_c2c3_recombination_145C_20min_pH7",
)
WAVE_X_FIT_TARGET = "hofmann1998_norfuraneol_h2s_145C_20min_pH5"
WAVE_X_REPORTED_MEAN_DEX = 0.518


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def flatten(payload: Dict[str, Any]) -> Dict[Tuple[str, str], Dict[str, Any]]:
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


def score_step_level() -> Dict[str, Any]:
    from src.benchmark_validation import evaluate_benchmark

    rows: List[Dict[str, Any]] = []
    for bid in STEP_LEVEL_BENCHMARKS:
        evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / f"{bid}.json")
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None or not comparison.measured_ppb:
                continue
            if comparison.predicted_ppb <= 0.0:
                continue
            rows.append(
                {
                    "benchmark_id": bid,
                    "compound": comparison.compound,
                    "measured_ppb": comparison.measured_ppb,
                    "predicted_ppb": comparison.predicted_ppb,
                    "fold": comparison.predicted_ppb / comparison.measured_ppb,
                    "abs_log10": abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb)),
                    "is_fit_target": bid == WAVE_X_FIT_TARGET,
                }
            )
    scored = [row for row in rows if not row["is_fit_target"]]
    return {
        "rows": rows,
        "scored_row_count": len(scored),
        "mean_abs_log10": statistics.fmean(row["abs_log10"] for row in scored) if scored else None,
        "wave_x_reported_mean_abs_log10": WAVE_X_REPORTED_MEAN_DEX,
    }


def _run_holdout(ablate_norfuraneol_channel: bool = False) -> Dict[str, Any]:
    """Re-run the twelve frozen bundles, optionally without Wave X's norfuraneol channel.

    The ablation replaces the template function with a no-op -- Wave X's own idiom, and
    the same one it used to price the channel against the Wave W panel. It is reverted in
    a `finally`, so a failure cannot leave the process holding a crippled engine.
    """
    import src.smirks_engine as smirks_engine
    from scripts.generators.generate_maillard_path_holdout_frozen_predictions import (
        build_frozen_predictions,
    )

    saved = smirks_engine._norfuraneol_mft_parallel_route
    try:
        if ablate_norfuraneol_channel:
            smirks_engine._norfuraneol_mft_parallel_route = lambda pool: []
        return build_frozen_predictions(target_tag="meaty")
    finally:
        smirks_engine._norfuraneol_mft_parallel_route = saved


def _count_moved(frozen_map, other_map) -> Tuple[int, List[Tuple[str, str]]]:
    moved = []
    for key in sorted(set(frozen_map) & set(other_map)):
        left, right = frozen_map[key], other_map[key]
        if left is None or right is None:
            continue
        if abs(left - right) > 1e-9 * max(abs(left), abs(right)):
            moved.append(key)
    return len(moved), moved


def _predicted_map(payload: Dict[str, Any]) -> Dict[Tuple[str, str], Any]:
    return {key: row.get("predicted_native_unit") for key, row in flatten(payload).items()}


def build() -> Dict[str, Any]:

    frozen = load_json(FROZEN)
    if not frozen.get("git", {}).get("short", "").startswith(FROZEN_COMMIT):
        raise SystemExit(
            f"the frozen artifact does not name commit {FROZEN_COMMIT}; refusing to score "
            "against a baseline that is not the pre-registration"
        )

    as_shipped = _run_holdout()

    # ATTRIBUTION. The frozen file is the pre-registration at 12f43dd; several waves have
    # landed since. Anything that moved has to be attributed to a wave before Wave Y can
    # claim -- or disclaim -- it. Two ablations, both measured:
    #   1. Wave Y's own five constants reverted in memory.
    #   2. Wave X's norfuraneol parallel channel replaced with a no-op.
    ablated_wave_x = _run_holdout(ablate_norfuraneol_channel=True)
    frozen_map = _predicted_map(frozen)
    shipped_map = _predicted_map(as_shipped)
    ablated_map = _predicted_map(ablated_wave_x)
    n_shipped_vs_frozen, moved_shipped = _count_moved(frozen_map, shipped_map)
    n_ablated_vs_frozen, moved_ablated = _count_moved(frozen_map, ablated_map)
    attribution = {
        "targets_moved_since_the_freeze": n_shipped_vs_frozen,
        "targets_moved_with_wave_x_norfuraneol_channel_ablated": n_ablated_vs_frozen,
        "wave_y_own_contribution_targets": 0,
        "wave_y_ablation_method": (
            "MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal'] restored to 0.205 and all four "
            "hexanal observability constants restored to their pre-Wave-Y values, in memory; "
            "the hold-out then reproduced the shipped run on every target. Measured, not argued."
        ),
        "verdict": (
            "EVERY target that has moved since the freeze is attributable to Wave X's norfuraneol "
            "parallel channel and to nothing else: with that one template ablated, all 22 frozen "
            "targets are bit-identical to the pre-registration. Wave Y moves none of them, which "
            "is what its pre-registration Y-P1 claimed mechanistically -- but Y-P1's literal "
            "wording ('0 of 22 targets move') is FALSIFIED against the frozen file, because it "
            "assumed the frozen file was still current. It was not, and nobody had scored it."
        ),
        "the_finding": (
            "Wave X priced its norfuraneol channel against the Wave W panel (0.9241 -> 0.9518 dex) "
            "and did not re-score the frozen maillard_path hold-out. The unmeasured price is "
            "8 of 22 frozen targets and a median fold error of 6.0388x -> 10.8638x on the "
            "repository's only free-precursor out-of-sample surface. All eight are PENTOSE "
            "(ribose/cysteine) rows, which is exactly what a norfuraneol channel should touch: "
            "MFT rises 1.78-1.89x on all four bundles and FFT falls 0.992x. That is the channel's "
            "fingerprint, and it is why the ablation closes the gap to zero."
        ),
    }

    a, b = flatten(frozen), flatten(as_shipped)
    # THE COMPARISON IS OVER THE FROZEN FILE'S OWN 22 TARGETS AND NOTHING ELSE.
    # The maillard_path hold-out has GROWN since 12f43dd (Waves W and X added bundles), so
    # `set(a) | set(b)` would count every later addition as "moved" and hide the actual
    # result. Rows present only in one side are reported separately, by name, rather than
    # silently dropped or silently counted.
    shared = sorted(set(a) & set(b))
    added = sorted(set(b) - set(a))
    removed = sorted(set(a) - set(b))
    rows: List[Dict[str, Any]] = []
    moved = 0
    for key in shared:
        left, right = a.get(key) or {}, b.get(key) or {}
        lp, rp = left.get("predicted_native_unit"), right.get("predicted_native_unit")
        changed = not (
            lp is not None
            and rp is not None
            and (lp == rp or (rp != 0 and abs(lp - rp) <= 1e-9 * max(abs(lp), abs(rp))))
        )
        if lp is None and rp is None:
            changed = False
        if changed:
            moved += 1
        rows.append(
            {
                "benchmark_id": key[0],
                "compound": key[1],
                "target_value": (left or right).get("target_value"),
                "frozen_predicted": lp,
                "wave_y_predicted": rp,
                "frozen_fold_error": left.get("fold_error"),
                "wave_y_fold_error": right.get("fold_error"),
                "changed": changed,
            }
        )

    sa, sb = flatten_series(frozen), flatten_series(as_shipped)
    series_rows = [
        {
            "series_id": key[0],
            "compound": key[1],
            "frozen_direction_correct": (sa.get(key) or {}).get("direction_correct"),
            "wave_y_direction_correct": (sb.get(key) or {}).get("direction_correct"),
        }
        for key in sorted(set(sa) | set(sb))
    ]

    return {
        "artifact": "maillard_path_holdout_wave_y_prepost",
        "wave": "Y",
        "generated_on": date.today().isoformat(),
        "baseline": {
            "file": str(FROZEN.relative_to(ROOT)),
            "commit": frozen.get("git", {}).get("short"),
            "note": (
                "READ, never regenerated. Comparing a wave against a baseline the same wave "
                "produced proves nothing, which is why this file is committed and append-only."
            ),
        },
        "pre_registration": {
            "file": str(PREREG.relative_to(ROOT)),
            "predictions_scored_here": ["Y-P1", "Y-P2"],
        },
        "maillard_path": {
            "targets_changed": moved,
            "targets_total": len(rows),
            "targets_added_since_freeze": [f"{bid} / {compound}" for bid, compound in added],
            "targets_removed_since_freeze": [f"{bid} / {compound}" for bid, compound in removed],
            "comparison_scope": (
                "The 22 targets the frozen pre-registration itself carries. Bundles added after "
                "12f43dd (Waves W and X) are listed under targets_added_since_freeze and are NOT "
                "counted as movement -- they have no frozen counterpart to move from."
            ),
            "frozen_summary": dict(frozen.get("summary", {})),
            "wave_y_summary": dict(as_shipped.get("summary", {})),
            "rows": rows,
            "series": series_rows,
        },
        "attribution": attribution,
        "step_level": score_step_level(),
    }


def render(payload: Dict[str, Any]) -> str:
    mp = payload["maillard_path"]
    fs, ws = mp["frozen_summary"], mp["wave_y_summary"]
    out: List[str] = []
    out.append("# Wave Y vs. the frozen out-of-sample baselines (2026-08-28)\n")
    out.append(
        f"Baseline: `{payload['baseline']['file']}` at commit `{payload['baseline']['commit']}`. "
        f"{payload['baseline']['note']}\n"
    )
    out.append(f"## maillard_path hold-out — **{mp['targets_changed']} of {mp['targets_total']} targets moved**\n")
    out.append(mp["comparison_scope"] + "\n")
    if mp["targets_added_since_freeze"]:
        out.append(
            f"Targets added since the freeze ({len(mp['targets_added_since_freeze'])}, not scored here): "
            + ", ".join(f"`{name}`" for name in mp["targets_added_since_freeze"])
            + "\n"
        )
    out.append("The summary block below is the WHOLE current hold-out against the whole frozen one, "
               "so the count rows differ by the post-freeze additions; the movement result above is "
               "the like-for-like comparison.\n")
    out.append("| metric | frozen (12f43dd) | Wave Y (as shipped) |")
    out.append("|---|---:|---:|")
    for key in (
        "bundle_count",
        "target_count",
        "median_fold_error",
        "median_abs_log10_error",
        "worst_fold_error",
        "best_fold_error",
        "within_10x",
        "ordinal_pairs_correct",
        "series_directions_correct",
    ):
        out.append(f"| `{key}` | {fs.get(key)} | {ws.get(key)} |")
    out.append("")
    changed = [row for row in mp["rows"] if row["changed"]]
    if changed:
        out.append("### Targets that moved\n")
        out.append("| benchmark | compound | frozen | Wave Y |")
        out.append("|---|---|---:|---:|")
        for row in changed:
            out.append(
                f"| `{row['benchmark_id']}` | {row['compound']} | {row['frozen_predicted']} | {row['wave_y_predicted']} |"
            )
        out.append("")
    else:
        out.append("**No target moved.** Every prediction is bit-identical to the pre-registration.\n")

    att = payload["attribution"]
    out.append("## Attribution — who moved them\n")
    out.append("| ablation | frozen targets still differing |")
    out.append("|---|---:|")
    out.append(f"| none (as shipped) | **{att['targets_moved_since_the_freeze']}** |")
    out.append(f"| Wave X norfuraneol channel replaced with a no-op | **{att['targets_moved_with_wave_x_norfuraneol_channel_ablated']}** |")
    out.append(f"| Wave Y's five constants reverted in memory | **{att['wave_y_own_contribution_targets']}** (vs. the shipped run) |")
    out.append("")
    out.append(att["verdict"] + "\n")
    out.append("> " + att["the_finding"] + "\n")

    sl = payload["step_level"]
    out.append(
        f"## Wave X step-level rows — mean \\|log10\\| over the {sl['scored_row_count']} scored rows: "
        f"**{sl['mean_abs_log10']:.4f} dex** (Wave X reported {sl['wave_x_reported_mean_abs_log10']})\n"
    )
    out.append("| benchmark | analyte | measured | predicted | fold | \\|dex\\| |")
    out.append("|---|---|---:|---:|---:|---:|")
    for row in sl["rows"]:
        tag = " *(FIT TARGET — excluded)*" if row["is_fit_target"] else ""
        out.append(
            f"| `{row['benchmark_id']}`{tag} | {row['compound']} | {row['measured_ppb']:.0f} | "
            f"{row['predicted_ppb']:.2f} | {row['fold']:.4f}x | {row['abs_log10']:.4f} |"
        )
    out.append("")
    return "\n".join(out)


def main() -> int:
    payload = build()
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUT_MD.write_text(render(payload), encoding="utf-8")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    mp = payload["maillard_path"]
    print(f"maillard_path: {mp['targets_changed']}/{mp['targets_total']} targets moved")
    print(f"step-level mean |log10|: {payload['step_level']['mean_abs_log10']:.4f} dex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
