#!/usr/bin/env python
"""
Wave B18 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-08).

`results/validation/kinetic_core_b18_prereg.md` sec. 4, computed from the frozen B18 fit report (and the
information-only `nosink` variant) against the engine as it stands, and written to
`results/validation/kinetic_core_b18_ship_rule.{json,md}`.

  T1  each of Zhou 2024's six rates within 0.3 dex; each barrier inside its printed-to-refit band
  T2  every scored panel row's predicted values move by less than 0.05 dex against the tracked scorecard
  T3  Leahy 1989's 95 C / 2 h distribution: pyrazine : 2,5-dimethylpyrazine within 0.5 dex (glycine for lysine);
      methylpyrazine's share reported (within 0.5 dex -> the mixed route ships)
  T4  Leahy's 2 h total within tenfold (reported)
  T5  Laplace: sigma below one decade on the two log10 constants (the barriers' sigma and bound state reported)
  T6  k(pH 5) / k(pH 9) for pyrazine at 95 C between 1/60 and 1/20 (reported)
  extra: Yu 2018's and Leahy's whole-cascade barriers against the model's apparent ones; the glyoxal-sink
         conditionality sized from the nosink variant
Ship rule: SHIP if T1, T2 and T5 hold.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import operative_parameters, trunk_conditions  # noqa: E402
from src.kinetic_core.engine import b1_fitted  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_pyrazine import with_fitted_pyrazine  # noqa: E402

V = data_paths.VALIDATION_DIR
REPORT = V / "kinetic_core_b18_fit_report.json"
NOSINK = V / "kinetic_core_b18_nosink_fit_report.json"
SCORECARD = V / "core_panel_scores.json"
OUT = V / "kinetic_core_b18_ship_rule.json"
CELSIUS = 273.15
MW = {"PZ": 80.09, "MPZ": 94.12, "DMP": 108.14}
#: Leahy 1989 ch. 7 Table III, lysine + glucose, 95 C, 2 h, pH 9: percent of total pyrazines and the total
LEAHY_TABLE_III = {"PZ": 55.8, "MPZ": 41.5, "DMP": 2.4, "total_ppm": 13.1}
#: Yu 2018 Table 1 (thermal arm, glucose 100 + glycine 100 mM, pH 10, 70-90 C): Ea kJ/mol +/- SE
YU_EA = {"DMP": (99.8, 6.7), "TMP": (117.8, 18.0), "TETRA": (104.1, 5.1)}
#: Leahy 1989 ch. 7 Table II, lysine + glucose, pH 9, 75-95 C: whole-cascade Ea kJ/mol
LEAHY_EA = {"PZ": 149.8, "MPZ": 153.1, "DMP": 177.0}


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def _params(fr: Dict[str, float]) -> Dict[str, Any]:
    p = dict(operative_parameters(b1_fitted()))
    p.update(with_fitted_pyrazine(fr["log10_k_go_ak_100C"], fr["ea_go_ak_kj_mol"], fr["log10_k_mgo_ak_100C"], fr["ea_mgo_ak_kj_mol"]))
    return p


def _run(fr, initial, t_c, minutes, ph, points=13):
    slopes = (fr["ph_slope_above_7_decades_per_unit"], fr["ph_slope_below_7_decades_per_unit"])
    params, _ = trunk_conditions.apply(_params(fr), SimpleNamespace(ph=ph, water_activity=None), pyrazine_slopes=slopes)
    return integrate(params, t_c + CELSIUS, initial, np.linspace(0.0, minutes, points), rtol=1e-8, atol=1e-14)


def _mean_rate(run, species, minutes):
    c = run.series(species)
    return float(c[-1] - c[0]) / minutes * 1000.0


def t1(report):
    res = report["residual_by_row_dex"]
    zhou = {k: v for k, v in res.items() if k.startswith("zhou_")}
    fr = report["frozen_parameters"]["pyrazine"]
    b = report["bounds"]
    inside = {k: bool(b[k][0] - 1e-9 <= fr[k] <= b[k][1] + 1e-9) for k in ("ea_go_ak_kj_mol", "ea_mgo_ak_kj_mol")}
    worst = max(zhou.items(), key=lambda kv: abs(kv[1]))
    return {"worst_row": worst[0], "worst_dex": worst[1], "rows_dex": zhou, "barriers_inside_band": inside,
            "pass": bool(abs(worst[1]) < 0.3 and all(inside.values()))}


def _predicted_leaves(o, path=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(_predicted_leaves(v, f"{path}/{k}"))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(_predicted_leaves(v, f"{path}[{i}]"))
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        if "pred" in path.lower():
            out[path] = float(o)
    return out


def t2():
    """The tracked scorecard (git HEAD) against the live one: every predicted number's log10 change."""
    from src.kinetic_core import scoring

    try:
        tracked = json.loads(subprocess.check_output(["git", "show", "HEAD:" + data_paths.rel(SCORECARD)], cwd=ROOT, text=True))
    except Exception as exc:  # pragma: no cover
        return {"status": f"tracked scorecard unavailable: {exc}", "pass": False}
    live = scoring.score_panel()
    old = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in tracked["benchmarks"]}
    new = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in live["benchmarks"]}
    deltas = []
    for bid, leaves in old.items():
        for path, v in leaves.items():
            w = new.get(bid, {}).get(path)
            if w is None or v <= 0 or w <= 0:
                continue
            deltas.append((abs(math.log10(w / v)), bid, path))
    deltas.sort(reverse=True)
    worst = deltas[0] if deltas else (0.0, None, None)
    return {"n_compared": len(deltas), "worst_dex": worst[0], "worst_benchmark": worst[1], "worst_path": worst[2],
            "n_rows_moved_over_0.01_dex": sum(1 for d in deltas if d[0] > 0.01), "pass": bool(worst[0] < 0.05)}


def leahy_holdouts(fr):
    run = _run(fr, {"Glc": 100.0, "Gly": 100.0}, 95.0, 120.0, 9.0)
    ug = {s: float(run.series(s)[-1]) * MW[s] * 1000.0 for s in ("PZ", "MPZ", "DMP")}
    total_ug = sum(ug.values())
    share = {s: 100.0 * ug[s] / total_ug if total_ug > 0 else float("nan") for s in ug}
    model_pz_over_dmp = ug["PZ"] / ug["DMP"] if ug["DMP"] > 0 else float("inf")
    leahy_pz_over_dmp = LEAHY_TABLE_III["PZ"] / LEAHY_TABLE_III["DMP"]
    model_mpz_over_pz = ug["MPZ"] / ug["PZ"] if ug["PZ"] > 0 else float("inf")
    leahy_mpz_over_pz = LEAHY_TABLE_III["MPZ"] / LEAHY_TABLE_III["PZ"]
    d_pz_dmp = math.log10(model_pz_over_dmp / leahy_pz_over_dmp) if math.isfinite(model_pz_over_dmp) and model_pz_over_dmp > 0 else float("inf")
    d_mpz = math.log10(model_mpz_over_pz / leahy_mpz_over_pz) if math.isfinite(model_mpz_over_pz) and model_mpz_over_pz > 0 else float("inf")
    total_dex = math.log10(total_ug / (LEAHY_TABLE_III["total_ppm"] * 1000.0)) if total_ug > 0 else float("-inf")
    return {
        "model_ug_per_l": ug, "model_share_pct": share, "leahy_share_pct": {k: v for k, v in LEAHY_TABLE_III.items() if k != "total_ppm"},
        "T3_pz_over_dmp": {"model": model_pz_over_dmp, "leahy": leahy_pz_over_dmp, "dex": d_pz_dmp, "pass": bool(abs(d_pz_dmp) < 0.5)},
        "T3_mpz_share": {"model_mpz_over_pz": model_mpz_over_pz, "leahy": leahy_mpz_over_pz, "dex": d_mpz, "mixed_route_ships": bool(abs(d_mpz) < 0.5)},
        "T4_total": {"model_ug_per_l": total_ug, "leahy_ug_per_l": LEAHY_TABLE_III["total_ppm"] * 1000.0, "dex": total_dex, "pass": bool(abs(total_dex) < 1.0)},
        "note": "glycine stands in for lysine (declared); Leahy's Table III is percent of total pyrazines after 2 h at 95 C, pH 9 borate",
    }


def t6(fr):
    r = {}
    for ph in (9.0, 5.0):
        run = _run(fr, {"Glc": 100.0, "Gly": 100.0}, 95.0, 120.0, ph)
        r[ph] = _mean_rate(run, "PZ", 120.0)
    ratio = r[5.0] / r[9.0] if r[9.0] > 0 else float("nan")
    return {"k_ph5_over_k_ph9": ratio, "band": [1 / 60.0, 1 / 20.0], "pass": bool(1 / 60.0 <= ratio <= 1 / 20.0)}


def _apparent_ea(fr, initial, temps, minutes, ph, species):
    xs, ys = [], []
    for t_c in temps:
        rate = _mean_rate(_run(fr, initial, t_c, minutes, ph), species, minutes)
        if rate > 0:
            xs.append(1.0 / (t_c + CELSIUS)); ys.append(math.log(rate))
    if len(xs) < 2:
        return None
    slope = np.polyfit(xs, ys, 1)[0]
    return float(-slope * 8.314462618e-3)


def barriers(fr):
    out = {"yu2018": {}, "leahy1989": {}}
    ea = _apparent_ea(fr, {"Glc": 100.0, "Gly": 100.0}, (70.0, 80.0, 90.0), 90.0, 10.0, "DMP")
    out["yu2018"]["DMP"] = {"model_apparent_kj_mol": ea, "yu": YU_EA["DMP"], "note": "glucose 100 + glycine 100 mM, pH 10, 70-90 C, 90 min window; Yu's k time unit is unprinted so only the barrier is compared"}
    for species in ("PZ", "MPZ", "DMP"):
        ea = _apparent_ea(fr, {"Glc": 100.0, "Gly": 100.0}, (75.0, 85.0, 95.0), 120.0, 9.0, species)
        out["leahy1989"][species] = {"model_apparent_kj_mol": ea, "leahy_whole_cascade_kj_mol": LEAHY_EA[species]}
    return out


def sink_conditionality(report, nosink):
    if nosink is None:
        return {"status": "nosink variant not run"}
    a, b = report["frozen_parameters"]["pyrazine"], nosink["frozen_parameters"]["pyrazine"]
    d = report["diagnostics"]
    return {
        "glyoxal_remaining_fraction_at_120min_in_zhou_pot": {k: v for k, v in d.items() if "GO_remaining" in k},
        "methylglyoxal_remaining_fraction_at_120min_in_zhou_pot": {k: v for k, v in d.items() if "MGO_remaining" in k},
        "pyrazine_rate_0_60_over_0_120_min": d.get("zhou_go_120_PZ_rate_0_60_over_0_120"),
        "log10_k_go_ak_shift_when_the_dry_glass_glyoxal_sink_is_zeroed": b["log10_k_go_ak_100C"] - a["log10_k_go_ak_100C"],
        "log10_k_mgo_ak_shift": b["log10_k_mgo_ak_100C"] - a["log10_k_mgo_ak_100C"],
        "nosink_cost": nosink["objective"]["final_cost"], "b18_cost": report["objective"]["final_cost"],
        "reading": ("The B13 glyoxal sink (Kocadagli's dry glass at 180 C, barrier fixed to zero) removes 98 % of a fed 20 mM glyoxal "
                    "in two hours at 100 C, so the pyrazine growth in the modelled Zhou pot is not linear (Zhou's is, figure-only) and "
                    "the fitted glyoxal Strecker constant is higher by the shift above than it would be over a constant pool. The "
                    "methylglyoxal pot loses its dicarbonyl too, through the trunk's own B1 melanoidin sink and the B7 furanone step; "
                    "those are fitted or measured constants and are not zeroed here. Recorded as a conditionality on every pyrazine "
                    "answer; the fix is a wave on the dicarbonyl sinks in water, not a move of this one."),
    }


def render(payload: Dict[str, Any]) -> str:
    t = payload
    L = [f"# Wave B18 ship rule: {t['verdict']}", "",
         f"*Rule: {t['rule']}. Pre-registration `{t['prereg']}`.*", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 Zhou rates | worst {t['T1']['worst_row']} {t['T1']['worst_dex']:+.3f} dex; barriers inside band {t['T1']['barriers_inside_band']} | {t['T1']['pass']} |",
         f"| T2 panel unchanged | {t['T2'].get('n_compared')} predicted numbers compared; worst {t['T2'].get('worst_dex', float('nan')):.4f} dex ({t['T2'].get('worst_benchmark')}) | {t['T2']['pass']} |",
         f"| T3 Leahy distribution | pyrazine : 2,5-DMP model {t['leahy']['T3_pz_over_dmp']['model']:.3g} vs {t['leahy']['T3_pz_over_dmp']['leahy']:.3g} ({t['leahy']['T3_pz_over_dmp']['dex']:+.2f} dex); methylpyrazine/pyrazine model {t['leahy']['T3_mpz_share']['model_mpz_over_pz']:.3g} vs {t['leahy']['T3_mpz_share']['leahy']:.3g} ({t['leahy']['T3_mpz_share']['dex']:+.2f} dex) | {t['leahy']['T3_pz_over_dmp']['pass']} (mixed route ships: {t['leahy']['T3_mpz_share']['mixed_route_ships']}) |",
         f"| T4 Leahy total | model {t['leahy']['T4_total']['model_ug_per_l']:.3g} vs {t['leahy']['T4_total']['leahy_ug_per_l']:.3g} ug/L ({t['leahy']['T4_total']['dex']:+.2f} dex) | {t['leahy']['T4_total']['pass']} |",
         f"| T5 identification | log10 k sigma {t['T5']['log10_sigma']}; barrier sigma {t['T5']['ea_sigma']} kJ/mol, on bound {t['T5']['ea_on_bound']} | {t['T5']['pass']} |",
         f"| T6 pH direction | k(5)/k(9) = {t['T6']['k_ph5_over_k_ph9']:.4f} (band 1/60 to 1/20) | {t['T6']['pass']} |",
         "", "## Barriers against the hold-out laboratories", ""]
    b = t["barriers"]
    L.append(f"- Yu 2018, 2,5-dimethylpyrazine, glucose + glycine pH 10, 70-90 C: model apparent {b['yu2018']['DMP']['model_apparent_kj_mol']:.1f} kJ/mol vs {b['yu2018']['DMP']['yu'][0]} +/- {b['yu2018']['DMP']['yu'][1]}")
    for s, v in b["leahy1989"].items():
        L.append(f"- Leahy 1989, {s}, pH 9, 75-95 C: model apparent {v['model_apparent_kj_mol']:.1f} vs whole-cascade {v['leahy_whole_cascade_kj_mol']}")
    sc = t["sink_conditionality"]
    L += ["", "## The glyoxal-sink conditionality", "", f"- {sc.get('reading', sc.get('status'))}", ""]
    for k, v in sc.items():
        if k != "reading":
            L.append(f"- {k}: {v}")
    L += ["", "## Leahy 95 C / 2 h, pH 9 (glycine for lysine)", "",
          f"- model ug/L: {t['leahy']['model_ug_per_l']}", f"- model shares %: {t['leahy']['model_share_pct']}", f"- Leahy shares %: {t['leahy']['leahy_share_pct']}"]
    return "\n".join(L) + "\n"


def main() -> int:
    report = _read(REPORT)
    nosink = _read(NOSINK) if NOSINK.exists() else None
    fr = report["frozen_parameters"]["pyrazine"]
    lap = report["laplace"]
    names = lap["coordinates"]
    sig = dict(zip(names, lap["sigma"]))
    onb = dict(zip(names, lap["on_bound"]))
    T5 = {"log10_sigma": {k: sig[k] for k in ("log10_k_go_ak_100C", "log10_k_mgo_ak_100C")},
          "ea_sigma": {k: sig[k] for k in ("ea_go_ak_kj_mol", "ea_mgo_ak_kj_mol")},
          "ea_on_bound": {k: onb[k] for k in ("ea_go_ak_kj_mol", "ea_mgo_ak_kj_mol")},
          "slope_sigma": {k: sig[k] for k in names if k.startswith("ph_slope")},
          "pass": bool(all(sig[k] < 1.0 for k in ("log10_k_go_ak_100C", "log10_k_mgo_ak_100C")))}
    T1 = t1(report)
    T2 = t2()
    leahy = leahy_holdouts(fr)
    T6 = t6(fr)
    ships = bool(T1["pass"] and T2["pass"] and T5["pass"])
    payload = {
        "artifact": "kinetic_core_b18_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b18_ship_rule",
                                                  generated_by="scripts/generators/generate_kinetic_core_b18_ship_rule.py",
                                                  inputs=[p for p in (REPORT, NOSINK) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b18_prereg.md"),
        "rule": "SHIP if T1 (six Zhou rates within 0.3 dex, barriers inside band), T2 (no scored panel value moves 0.05 dex) and T5 (log10 k sigma below one decade) hold",
        "T1": T1, "T2": T2, "leahy": leahy, "T5": T5, "T6": T6, "barriers": barriers(fr),
        "sink_conditionality": sink_conditionality(report, nosink),
        "frozen_parameters": {"pyrazine": fr},
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], T1["worst_dex"], "| T2", T2["pass"], T2.get("worst_dex"), "| T3", leahy["T3_pz_over_dmp"]["dex"], "mixed", leahy["T3_mpz_share"]["dex"],
          "| T4", leahy["T4_total"]["dex"], "| T5", T5["pass"], "| T6", T6["k_ph5_over_k_ph9"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
