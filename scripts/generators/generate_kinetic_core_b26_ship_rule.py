#!/usr/bin/env python
"""
Wave B26 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-09).

`results/validation/kinetic_core_b26_prereg.md` section 4, computed from the LIVE matrix registry
against wave B4's FROZEN blind predictions, written to
`results/validation/kinetic_core_b26_ship_rule.{json,md}`.

  T1  arithmetic: the pooled n_alkanal constant is the geometric mean of the two rows, the branched
      surrogate follows the 2.81x/CH2 slope, the reference loading is unchanged, no other class moves
  T2  the flagship hold-out, DECISIVE: on the three Hong 2020 rows where the binding term is active
      every fold error falls, no sign inverts, and the ten-row verdict does not get worse
  T3  the evidence ceiling: the reversible term's share of each row's log-shift, before and after
  T4  the loading and the within-paper method spread, as a band on the hexanal prediction
  T5  nothing else moves: sealed keys still sealed, the unsaturation penalty exactly sqrt(2.81*4.95),
      the kinetic panel scorecard bit-for-bit unchanged
Ship rule: SHIP if T1, T2 and T5 hold; T3 and T4 are reported.

WHY THIS DOES NOT RE-RUN WAVE B4. `kinetic_core_b4_frozen_predictions.json` is the record of a BLIND
prediction, written before that wave read the paired thresholds. Re-running B4's generator would
overwrite it with a prediction made by somebody who has now seen the answer, and would silently
replace a pre-registration with a post-hoc fit. That generator now refuses to do so without an
explicit --refreeze; this ship rule reads the frozen file and writes only its own artifact.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
FROZEN = V / "kinetic_core_b4_frozen_predictions.json"
HOLDOUT_VALUES = V / "holdout_frozen/hong2020_paired_thresholds.json"
OUT = V / "kinetic_core_b26_ship_rule.json"
#: The prereg's decisive rows: the only Hong compounds whose class carries a binding constant.
ACTIVE = ("hexanal", "3_methylbutanal", "2_methylbutanal")
#: The two rows the wave pools, and the two it declines to pool.
SHIPPED_KEYS = ("kg_hexanal_pea", "kg_z_2_penten_1_ol_pea")
QUARANTINED_KEY = "kg_t_2_octenal_pea"


def _read(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _name_key(name: str) -> str:
    from generate_kinetic_core_b4_holdout import _name_key as k  # noqa: E402

    return k(name)


def _measured() -> Dict[str, float]:
    payload = _read(HOLDOUT_VALUES)
    out: Dict[str, float] = {}
    for entry in payload["compounds"]:
        ratio = entry.get("ratio_soy_over_water")
        if ratio is None and entry.get("threshold_water"):
            ratio = entry["threshold_soy"] / entry["threshold_water"]
        out[_name_key(entry["name"])] = float(ratio)
    return out


def t1(frozen: Dict[str, Any]) -> Dict[str, Any]:
    from src.kinetic_core.matrix_oav import fit_class_binding_constants
    from src.kinetic_core.parameters_matrix import CHAIN_LENGTH_SLOPE_PER_CH2, REVERSIBLE_BINDING

    before = frozen["class_binding_constants"]
    after = fit_class_binding_constants()
    rows = {p.key: p for p in REVERSIBLE_BINDING}
    pooled = math.exp(
        (math.log(rows["kg_hexanal_dairy"].value) + math.log(rows["kg_hexanal_pea"].value)) / 2.0)
    moved, unmoved = {}, []
    for name, entry in after.items():
        old = before.get(name)
        if old is None:
            moved[name] = {"before": None, "after": entry["k_g_l_per_g"], "new_class": True}
            continue
        delta = abs(entry["k_g_l_per_g"] - old["k_g_l_per_g"])
        if delta > 1e-12:
            moved[name] = {"before": old["k_g_l_per_g"], "after": entry["k_g_l_per_g"],
                           "x": entry["k_g_l_per_g"] / old["k_g_l_per_g"] if old["k_g_l_per_g"] else None}
        else:
            unmoved.append(name)
    n_alk = after["n_alkanal"]
    branched = after["branched_alkanal"]
    ok = (abs(n_alk["k_g_l_per_g"] - pooled) < 1e-12
          and abs(branched["k_g_l_per_g"] - pooled / CHAIN_LENGTH_SLOPE_PER_CH2) < 1e-12
          and n_alk["reference_loading_g_per_l"] == before["n_alkanal"]["reference_loading_g_per_l"]
          and set(moved) == {"n_alkanal", "branched_alkanal", "alkenol"})
    return {"geometric_mean_expected": pooled, "n_alkanal_after": n_alk["k_g_l_per_g"],
            "n_alkanal_members": n_alk["members"], "branched_after": branched["k_g_l_per_g"],
            "reference_loading_before_after": [before["n_alkanal"]["reference_loading_g_per_l"],
                                               n_alk["reference_loading_g_per_l"]],
            "classes_moved": moved, "classes_unmoved": sorted(unmoved), "pass": bool(ok)}


def t2_t3(frozen: Dict[str, Any]) -> Dict[str, Any]:
    from src.kinetic_core.matrix_oav import decompose_residual, predict_matrix_shift

    measured = _measured()
    frozen_rows = {r["compound"]: r for r in frozen["predictions"]}
    rows: List[Dict[str, Any]] = []
    worse: List[str] = []
    for key, row in frozen_rows.items():
        obs = measured.get(key)
        if obs is None:
            continue
        old_pred = float(row["predicted_matrix_over_water_ratio"])
        new = predict_matrix_shift(key, row["matrix"])
        new_pred = new.predicted_ratio
        old_fold = max(obs, old_pred) / min(obs, old_pred)
        new_fold = max(obs, new_pred) / min(obs, new_pred)
        decomposition = decompose_residual(new, obs)
        share_after = (decomposition.per_term_decades["reversible_binding"]
                       / decomposition.measured_decades if decomposition.measured_decades else 0.0)
        share_before = (math.log10(old_pred) / math.log10(obs)
                        if old_pred > 0 and obs > 0 and obs != 1 else 0.0)
        sign_ok = ((obs > 1.0) == (new_pred > 1.0)) or new_pred == 1.0 and old_pred == 1.0
        active = key in ACTIVE
        if active and (new_fold > old_fold + 1e-12 or not sign_ok):
            worse.append(key)
        rows.append({"compound": key, "active_term": active, "measured": obs,
                     "predicted_before": old_pred, "predicted_after": new_pred,
                     "fold_before": old_fold, "fold_after": new_fold,
                     "sign_ok": bool(sign_ok), "state": new.state,
                     "explained_share_before": share_before, "explained_share_after": share_after,
                     "ceiling_flag": [f for f in decomposition.flags if "CEILING" in f]})
    within_before = sum(1 for r in rows if r["fold_before"] <= 5.0)
    within_after = sum(1 for r in rows if r["fold_after"] <= 5.0)
    signs_before = sum(1 for r in rows if (r["measured"] > 1.0) == (r["predicted_before"] > 1.0)
                       and r["predicted_before"] != 1.0)
    signs_after = sum(1 for r in rows if (r["measured"] > 1.0) == (r["predicted_after"] > 1.0)
                      and r["predicted_after"] != 1.0)
    t2 = {"rows": rows, "rows_made_worse": worse,
          "within_5x_before_after": [within_before, within_after],
          "signs_correct_before_after": [signs_before, signs_after],
          "pass": bool(not worse and within_after >= within_before and signs_after >= signs_before)}
    active_rows = [r for r in rows if r["active_term"]]
    t3 = {"ceiling": 0.25,
          "explained_share_by_row": {r["compound"]: [r["explained_share_before"],
                                                     r["explained_share_after"]] for r in active_rows},
          "rows_now_over_the_ceiling": [r["compound"] for r in active_rows
                                        if r["explained_share_after"] > 0.25],
          "reading_declared_in_the_prereg": (
              "The cap was computed from ONE compound in beef and one dairy protein. A plant "
              "isolate that binds an alkanal 22x harder than cow's milk does is a reason to "
              "doubt that the cap transfers, not a reason to shrink a measured constant. The "
              "layer's flag stays and starts firing; that is the flag doing its job."),
          "reported_only": True}
    return {"T2": t2, "T3": t3}


def t4() -> Dict[str, Any]:
    """The loading is inherited, not restated; and the paper's two headspace routes disagree 2.5x."""
    from src.kinetic_core.matrix_oav import predict_matrix_shift
    from src.kinetic_core.parameters_matrix import MATRIX_LOADING, REVERSIBLE_BINDING

    rows = {p.key: p for p in REVERSIBLE_BINDING}
    k_pea = rows["kg_hexanal_pea"].value
    k_dairy = rows["kg_hexanal_dairy"].value
    hong_loading = MATRIX_LOADING["soy_paste_hong"].protein_g_per_l
    out = {}
    for label, factor in (("loading_half_5_g_per_L", 2.0), ("as_printed_10_g_per_L", 1.0),
                          ("loading_double_20_g_per_L", 0.5),
                          ("method_spread_low_2.5x", 1 / 2.5), ("method_spread_high_2.5x", 2.5)):
        k = k_pea * factor
        pooled = math.exp((math.log(k_dairy) + math.log(k)) / 2.0)
        out[label] = {"kg_pea_l_per_g": k, "pooled_n_alkanal": pooled,
                      "hong_hexanal_prediction_x": 1.0 + pooled * hong_loading}
    out["live"] = {"hong_hexanal_prediction_x": predict_matrix_shift("hexanal", "soy_paste_hong").predicted_ratio}
    out["reported_only"] = True
    return out


def t5(frozen: Dict[str, Any]) -> Dict[str, Any]:
    from src.kinetic_core.matrix_oav import fit_unsaturation_penalty
    from src.kinetic_core.parameters_matrix import HOLDOUT_SEALED_BINDING, REVERSIBLE_BINDING

    carried = {p.key for p in REVERSIBLE_BINDING}
    penalty = fit_unsaturation_penalty()
    expected = math.sqrt(2.81 * 4.95)
    try:
        tracked = json.loads(subprocess.check_output(
            ["git", "show", "HEAD:" + data_paths.rel(V / "core_panel_scores.json")], cwd=ROOT, text=True))
    except Exception as exc:  # pragma: no cover
        return {"status": f"tracked scorecard unavailable: {exc}", "pass": False}
    from src.kinetic_core import scoring

    live = scoring.score_panel()

    def strip(payload):
        return json.dumps({b["benchmark_id"]: b.get("compounds") for b in payload["benchmarks"]},
                          sort_keys=True, default=str)

    panel_identical = strip(tracked) == strip(live)
    return {"sealed_keys_still_valueless": bool(carried.isdisjoint(HOLDOUT_SEALED_BINDING)),
            "n_sealed": len(HOLDOUT_SEALED_BINDING),
            "unsaturation_penalty_x": penalty["penalty_x"], "unsaturation_n_fit_rows": penalty["n_fit_rows"],
            "unsaturation_excluded": penalty["excluded"],
            "kinetic_panel_identical": bool(panel_identical),
            "pass": bool(carried.isdisjoint(HOLDOUT_SEALED_BINDING)
                         and abs(penalty["penalty_x"] - expected) < 1e-9
                         and penalty["n_fit_rows"] == 2 and panel_identical)}


def render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4, T5 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"]
    L = [f"# Wave B26 ship rule: {p['verdict']}", "",
         f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         f"Rows installed: {', '.join(p['shipped_keys'])} (pooled) and {p['quarantined_key']} "
         f"(quarantined as a binding constant and excluded from the unsaturation fit).", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 arithmetic | n_alkanal {T1['n_alkanal_after']:.5g} L/g from "
         f"{', '.join(T1['n_alkanal_members'])}; branched {T1['branched_after']:.5g}; reference "
         f"loading {T1['reference_loading_before_after'][0]:g} -> "
         f"{T1['reference_loading_before_after'][1]:g} g/L; classes moved "
         f"{sorted(T1['classes_moved'])}; unmoved {T1['classes_unmoved']} | {T1['pass']} |",
         f"| T2 flagship hold-out | rows made worse: {T2['rows_made_worse'] or 'none'}; within 5x "
         f"{T2['within_5x_before_after'][0]} -> {T2['within_5x_before_after'][1]} of "
         f"{len(T2['rows'])}; signs correct {T2['signs_correct_before_after'][0]} -> "
         f"{T2['signs_correct_before_after'][1]} | {T2['pass']} |",
         f"| T3 evidence ceiling | rows now over the {T3['ceiling']:.0%} cap: "
         f"{T3['rows_now_over_the_ceiling'] or 'none'} | reported |",
         f"| T4 the loading | hexanal in the Hong paste from "
         f"{T4['loading_double_20_g_per_L']['hong_hexanal_prediction_x']:.2f}x to "
         f"{T4['loading_half_5_g_per_L']['hong_hexanal_prediction_x']:.2f}x across a factor of four "
         f"in the assumed loading; live {T4['live']['hong_hexanal_prediction_x']:.2f}x | reported |",
         f"| T5 nothing else moves | {T5['n_sealed']} sealed keys still valueless; penalty "
         f"{T5['unsaturation_penalty_x']:.4f}x on {T5['unsaturation_n_fit_rows']} rows; kinetic "
         f"panel identical {T5['kinetic_panel_identical']} | {T5['pass']} |", "",
         "## The three rows that carry a binding term", "",
         "| compound | measured | predicted before | predicted after | fold before | fold after | "
         "explained share before | after |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in T2["rows"]:
        if not row["active_term"]:
            continue
        L.append(f"| {row['compound']} | {row['measured']:.4g}x | {row['predicted_before']:.4g}x | "
                 f"{row['predicted_after']:.4g}x | {row['fold_before']:.4g}x | {row['fold_after']:.4g}x | "
                 f"{row['explained_share_before']:.1%} | {row['explained_share_after']:.1%} |")
    L += ["", "## The seven rows that do not", "",
          "Unchanged, and unchanged for the same reason as before: the corpus supplies no binding "
          "constant for their classes, so the layer emits exactly 1.0 and reports the whole shift "
          "as unexplained residual. This wave adds a protein, not a class.", "",
          "## What the ceiling says now", "", f"> {T3['reading_declared_in_the_prereg']}", ""]
    return "\n".join(L)


def main() -> int:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    if not FROZEN.exists():
        raise SystemExit(f"{FROZEN} missing")
    frozen = _read(FROZEN)
    T1 = t1(frozen)
    both = t2_t3(frozen)
    T2, T3 = both["T2"], both["T3"]
    T4, T5 = t4(), t5(frozen)
    ships = bool(T1["pass"] and T2["pass"] and T5["pass"])
    payload = {
        "artifact": "kinetic_core_b26_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b26_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b26_ship_rule.py",
            inputs=[FROZEN, HOLDOUT_VALUES]),
        "prereg": data_paths.rel(V / "kinetic_core_b26_prereg.md"),
        "rule": ("SHIP if T1 (the arithmetic), T2 (every active hold-out row improves and no sign "
                 "inverts) and T5 (nothing else moves) hold; T3 and T4 reported"),
        "shipped_keys": list(SHIPPED_KEYS), "quarantined_key": QUARANTINED_KEY,
        "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], "| T2", T2["pass"], T2["rows_made_worse"],
          "| T3 over-ceiling", T3["rows_now_over_the_ceiling"], "| T5", T5["pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
