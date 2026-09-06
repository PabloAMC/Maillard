#!/usr/bin/env python
"""
Wave B10 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-06).

`results/validation/kinetic_core_b10_prereg.md` sec. 5 names six tests and a rule.
This script computes T1-T5 from the frozen B10, B10-noyil and B9 artifacts WITHOUT
switching the engine (B10's parameters are handed to `predict` as an override), and
writes `results/validation/kinetic_core_b10_ship_rule.{json,md}`. T6 (the directional
panel) is read from the regenerated scorecard after the decision, not here.

  T1  both route barriers identified by the Laplace (finite sigma, not rank-deficient)
  T2  no B9 objective row's |residual| grows by more than 0.3 dex
  T3  the eight Yiltirak LEVEL rows (fit-adjacent): median fold error below 10x
  T4  hold-outs: Kang 140 C rung direction (MFT and FFT rise 120 -> 140), Hofmann 2002
      brew FFT loss at 80 C not worse than under B9
  T5  the leave-Yiltirak-out refit, reported alongside
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import artifact_io, data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, panel  # noqa: E402
from src.kinetic_core.engine import SULFUR, CoreDraw, core_parameters  # noqa: E402
from src.kinetic_core.ph_state import PhDrift  # noqa: E402
from src.kinetic_core.sulfur import integrate_sulfur  # noqa: E402

V = data_paths.VALIDATION_DIR
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
B10_REPORT = V / "kinetic_core_b10_fit_report.json"
B10_NOYIL_REPORT = V / "kinetic_core_b10_noyil_fit_report.json"
B10_LAPLACE = V / "kinetic_core_b10_laplace_covariance.json"
OUT = V / "kinetic_core_b10_ship_rule.json"
CELSIUS = 273.15
YILTIRAK_DIR = data_paths.MAILLARD_PATH_HOLDOUT_DIR
YILTIRAK_STEMS = (
    "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026",
)
KANG_OBSERVED_FOLD = {"MFT": 5.907 / 1.388, "FFT": 11.439 / 4.107}   # Kang 2026 SI Table S4, 140 over 120 C
BREW_K_OBSERVED = 0.023   # Hofmann 2002 Fig. 1, /min at 80 C


def _read(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _operative_and_drift(report: Dict[str, Any]):
    fr = report["frozen_parameters"]
    block = {
        "log10_k_ref_at_145C": dict(fr["log10_k_ref_at_145C"]),
        "lumped_formation_Ea_kJ_mol": float(fr["lumped_formation_Ea_kJ_mol"]),
        "decay_Ea_kJ_mol": dict(fr["decay_Ea_kJ_mol"]),
    }
    if fr.get("formation_Ea_by_route_kJ_mol"):
        block["formation_Ea_by_route_kJ_mol"] = dict(fr["formation_Ea_by_route_kJ_mol"])
    drift = PhDrift(acid_yield=float(fr["ph_drift"]["acid_yield_per_sink_event"]),
                    arp_amine_pka=float(fr["ph_drift"]["arp_secondary_ammonium_pKa"]))
    return core_parameters(SULFUR, frozen=block), drift


def _yiltirak_folds(operative, drift) -> List[Dict[str, Any]]:
    rows = []
    for stem in YILTIRAK_STEMS:
        bench = panel.load_bundle(YILTIRAK_DIR / f"{stem}.json")
        spec = panel.core_spec(bench, use_buffer=True)
        _, limiting = panel.limiting_precursor_molar(bench)
        for compound, target in panel.bundle_targets(bench).items():
            measured = panel.measured_value(bench, compound, target)
            run = engine.predict(spec, [compound], parameters=operative, draw=CoreDraw(ph_drift=drift))
            predicted = panel.core_native_value(run, compound, "ppb", limiting)
            fold = max(predicted / measured, measured / predicted) if predicted and measured else float("inf")
            rows.append({"bundle": stem, "compound": compound, "measured": measured,
                         "predicted": predicted, "fold": fold})
    return rows


def _kang_rung(operative, drift, B23) -> Dict[str, Any]:
    spec = B23.SYSTEMS["kang_ttca_120"]
    out = {}
    for t_c in (120.0, 140.0):
        run = integrate_sulfur(operative, t_c + CELSIUS, spec["initial"], np.array([0.0, float(spec["minutes"])]),
                               ph=float(spec["ph"]), buffer_spec=spec.get("buffer"), ph_drift=drift,
                               rtol=1e-8, atol=1e-14)
        out[t_c] = {s: run.final(s) for s in ("MFT", "FFT")}
    return {
        s: {"predicted_fold_140_over_120": (out[140.0][s] / out[120.0][s]) if out[120.0][s] > 0 else float("nan"),
            "observed_fold": KANG_OBSERVED_FOLD[s]}
        for s in ("MFT", "FFT")
    }


def _brew(operative) -> float:
    run = integrate_sulfur(operative, 80.0 + CELSIUS, {"FFT": 0.5, "MELE": 9.0 * 50.0, "OX": 1.0},
                           np.array([0.0, 60.0]), ph=5.2, rtol=1e-8, atol=1e-14)
    remaining = run.final("FFT") / 0.5
    return -math.log(max(remaining, 1e-12)) / 60.0


def main() -> int:
    for p in (B9_REPORT, B10_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    import generate_kinetic_core_b10_fit as B10  # noqa: E402  (configure(True) at import)
    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
    B10.configure(True)
    sigma = {r["id"]: float(r["sigma_log"]) for r in B23.ACTIVE_FIT_ROWS if r["kind"] != "ph_endpoint"}

    b9, b10 = _read(B9_REPORT), _read(B10_REPORT)
    noyil = _read(B10_NOYIL_REPORT) if B10_NOYIL_REPORT.exists() else None
    lap = _read(B10_LAPLACE) if B10_LAPLACE.exists() else None
    routes = b10["frozen_parameters"]["formation_Ea_by_route_kJ_mol"]

    # T1 ---------------------------------------------------------------
    t1 = {"status": "no laplace artifact"}
    if lap:
        ident = {c["key"] or c["block"]: bool(ok) for c, ok in zip(lap["coordinates"], lap["identified"])}
        sig = {c["key"] or c["block"]: s for c, s in zip(lap["coordinates"], lap["sigma"])}
        t1 = {"sugar_trunk": {"identified": ident.get("sugar_trunk"), "sigma_kj_mol": sig.get("sugar_trunk")},
              "thiol_assembly": {"identified": ident.get("thiol_assembly"), "sigma_kj_mol": sig.get("thiol_assembly")}}
        t1["pass"] = bool(t1["sugar_trunk"]["identified"] and t1["thiol_assembly"]["identified"])

    # T2 ---------------------------------------------------------------
    r9 = b9["members"][b9["best_start"]]["residual_by_row"]
    r10 = b10["members"][b10["best_start"]]["residual_by_row"]
    growth = []
    for rid, v9 in r9.items():
        if rid in r10 and rid in sigma:
            dex = (abs(float(r10[rid])) - abs(float(v9))) * sigma[rid]
            growth.append({"row": rid, "b9_residual": float(v9), "b10_residual": float(r10[rid]), "growth_dex": dex})
    worst = max(growth, key=lambda g: g["growth_dex"])
    t2 = {"rows_compared": len(growth), "worst": worst, "n_over_0_3_dex": sum(1 for g in growth if g["growth_dex"] > 0.3),
          "n_over_0_5_dex": sum(1 for g in growth if g["growth_dex"] > 0.5), "pass": worst["growth_dex"] <= 0.3}

    # T3 ---------------------------------------------------------------
    op9, d9 = _operative_and_drift(b9)
    op10, d10 = _operative_and_drift(b10)
    y9, y10 = _yiltirak_folds(op9, d9), _yiltirak_folds(op10, d10)
    med9, med10 = statistics.median(r["fold"] for r in y9), statistics.median(r["fold"] for r in y10)
    t3 = {"b9_median_fold": med9, "b10_median_fold": med10, "rows_b10": y10, "pass": med10 < 10.0,
          "improved": med10 < med9, "within_3x_b9": sum(r["fold"] <= 3 for r in y9), "within_3x_b10": sum(r["fold"] <= 3 for r in y10)}

    # T4 ---------------------------------------------------------------
    kang9, kang10 = _kang_rung(op9, d9, B23), _kang_rung(op10, d10, B23)
    brew9, brew10 = _brew(op9), _brew(op10)
    brew_fold9 = max(brew9 / BREW_K_OBSERVED, BREW_K_OBSERVED / brew9)
    brew_fold10 = max(brew10 / BREW_K_OBSERVED, BREW_K_OBSERVED / brew10)
    t4 = {"kang_140_over_120": {"b9": kang9, "b10": kang10,
                                "direction_pass": all(kang10[s]["predicted_fold_140_over_120"] > 1.0 for s in ("MFT", "FFT"))},
          "hofmann2002_brew_80C": {"observed_k_per_min": BREW_K_OBSERVED, "b9_k_per_min": brew9, "b10_k_per_min": brew10,
                                   "b9_fold": brew_fold9, "b10_fold": brew_fold10, "not_worse": brew_fold10 <= brew_fold9 * 1.05}}
    t4["pass"] = t4["kang_140_over_120"]["direction_pass"]

    # T5 ---------------------------------------------------------------
    t5 = {"status": "no leave-Yiltirak-out report"}
    if noyil:
        nr = noyil["frozen_parameters"]["formation_Ea_by_route_kJ_mol"]
        t5 = {"with_folds": routes, "without_folds": nr,
              "difference_kj_mol": {k: routes[k] - nr[k] for k in routes},
              "cost_with_folds_on_shared_rows": b10["objective"]["sum_r2_level_shared_with_b2_4"],
              "cost_without_folds_on_shared_rows": noyil["objective"]["sum_r2_level_shared_with_b2_4"]}
        if lap:
            t5["difference_in_laplace_sigmas"] = {
                k: (abs(routes[k] - nr[k]) / t1[k]["sigma_kj_mol"]) if t1.get(k, {}).get("sigma_kj_mol") else None
                for k in routes}

    # the rule --------------------------------------------------------
    ships = bool(t1.get("pass") and t2["pass"] and t4["pass"] and t3["improved"])
    verdict = "SHIP" if ships else ("RE-MERGE: barriers not identified" if lap and not t1.get("pass") else "DO NOT SHIP")
    payload = {
        "artifact": "kinetic_core_b10_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b10_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b10_ship_rule.py",
            inputs=[p for p in (B9_REPORT, B10_REPORT, B10_NOYIL_REPORT, B10_LAPLACE) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b10_prereg.md"),
        "route_barriers_kj_mol": routes,
        "b9_lumped_barrier_kj_mol": b9["frozen_parameters"]["lumped_formation_Ea_kJ_mol"],
        "active_bounds": b10.get("active_bounds"),
        "T1_laplace_identification": t1, "T2_in_sample_discipline": t2, "T3_yiltirak_levels_fit_adjacent": t3,
        "T4_holdouts": t4, "T5_leave_yiltirak_out": t5,
        "verdict": verdict,
        "rule": "SHIP if T1 and T2 and T4 (Kang direction) hold and T3 improved at all (prereg sec. 5)",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: routes {routes}; T1 {t1.get('pass')} T2 {t2['pass']} (worst {worst['row']} {worst['growth_dex']:+.2f} dex) "
          f"T3 {med9:.1f}x -> {med10:.1f}x T4 kang {t4['pass']} brew {brew_fold9:.1f}x -> {brew_fold10:.1f}x")
    return 0


def render(p: Dict[str, Any]) -> str:
    t1, t2, t3, t4, t5 = (p["T1_laplace_identification"], p["T2_in_sample_discipline"], p["T3_yiltirak_levels_fit_adjacent"],
                          p["T4_holdouts"], p["T5_leave_yiltirak_out"])
    out = [f"# Wave B10 ship rule -- **{p['verdict']}**", "",
           f"Prereg: `{p['prereg']}` sec. 5. Rule: {p['rule']}.", "",
           f"* route barriers: sugar trunk **{p['route_barriers_kj_mol']['sugar_trunk']:.1f}**, thiol assembly "
           f"**{p['route_barriers_kj_mol']['thiol_assembly']:.1f}** kJ/mol (B9's single lumped barrier: {p['b9_lumped_barrier_kj_mol']:.1f})",
           f"* active bounds at the optimum: {[a['key'] for a in (p.get('active_bounds') or [])] or 'none'}", "",
           "## T1 -- Laplace identification of the two barriers", ""]
    if "pass" in t1:
        for k in ("sugar_trunk", "thiol_assembly"):
            out.append(f"* {k}: identified = {t1[k]['identified']}, sigma = {t1[k]['sigma_kj_mol']}")
        out.append(f"* **{'PASS' if t1['pass'] else 'FAIL'}**")
    else:
        out.append(f"* {t1['status']}")
    out += ["", "## T2 -- no B9 row moved more than 0.3 dex", "",
            f"* {t2['rows_compared']} shared rows; worst growth {t2['worst']['growth_dex']:+.2f} dex on `{t2['worst']['row']}`; "
            f"{t2['n_over_0_3_dex']} rows over 0.3 dex, {t2['n_over_0_5_dex']} over 0.5 dex -> **{'PASS' if t2['pass'] else 'FAIL'}**",
            "", "## T3 -- the eight Yiltirak levels (FIT-ADJACENT, not out of sample)", "",
            f"* median fold error {t3['b9_median_fold']:.1f}x (B9) -> **{t3['b10_median_fold']:.1f}x** (B10); within 3x "
            f"{t3['within_3x_b9']} -> {t3['within_3x_b10']} of 8 -> below 10x: **{'PASS' if t3['pass'] else 'FAIL'}**; improved: {t3['improved']}",
            "", "| bundle | compound | measured | B10 predicted | fold |", "|---|---|---:|---:|---:|"]
    for r in t3["rows_b10"]:
        out.append(f"| {r['bundle']} | {r['compound']} | {r['measured']:.3g} | {r['predicted']:.3g} | {r['fold']:.1f} |")
    k = t4["kang_140_over_120"]; b = t4["hofmann2002_brew_80C"]
    out += ["", "## T4 -- hold-outs", "",
            f"* Kang 140 C over 120 C: MFT predicted x{k['b10']['MFT']['predicted_fold_140_over_120']:.2f} (B9 x{k['b9']['MFT']['predicted_fold_140_over_120']:.2f}; "
            f"observed x{KANG_OBSERVED_FOLD['MFT']:.2f}), FFT x{k['b10']['FFT']['predicted_fold_140_over_120']:.2f} (B9 x{k['b9']['FFT']['predicted_fold_140_over_120']:.2f}; "
            f"observed x{KANG_OBSERVED_FOLD['FFT']:.2f}) -> direction **{'PASS' if k['direction_pass'] else 'FAIL'}**",
            f"* Hofmann 2002 brew, FFT loss at 80 C: observed {b['observed_k_per_min']:.3f} /min; B9 {b['b9_k_per_min']:.4f} ({b['b9_fold']:.1f}x), "
            f"B10 {b['b10_k_per_min']:.4f} ({b['b10_fold']:.1f}x) -> not worse: {b['not_worse']}",
            "", "## T5 -- leave-Yiltirak-out", ""]
    if "with_folds" in t5:
        out.append(f"* with the folds: {t5['with_folds']}; without: {t5['without_folds']}; difference: "
                   f"{ {k: round(v, 1) for k, v in t5['difference_kj_mol'].items()} } kJ/mol"
                   + (f"; in Laplace sigmas: {t5.get('difference_in_laplace_sigmas')}" if t5.get('difference_in_laplace_sigmas') else ""))
        out.append(f"* shared-row cost: with folds {t5['cost_with_folds_on_shared_rows']:.2f}, without {t5['cost_without_folds_on_shared_rows']:.2f}")
    else:
        out.append(f"* {t5['status']}")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
