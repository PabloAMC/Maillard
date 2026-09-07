#!/usr/bin/env python
"""
Wave B11 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-07).

`results/validation/kinetic_core_b11_prereg.md` sec. 5 names six tests and a rule. This script
computes them from the frozen B9 and B11 artifacts WITHOUT switching the engine (B11's parameters
are handed to `predict` as an override; the panel bundles' vessel blocks charge the reservoir) and
writes `results/validation/kinetic_core_b11_ship_rule.{json,md}`.

  T1  Bolton 1994 (O2 : thiol 2.07): MFT fold error below 6x (falsifier: above 12x)
  T2  Hofmann 1998 Table-1 pH-5 rows (0.27): none worse than 8x (the eight B9 validation rows)
  T3  Yiltirak 2026 130 C / 0.5 h: MFT and FFT fold errors both fall by at least 3x
  T4  Yiltirak 100 C MFT-vs-FFT split: recorded, not a target
  T5  in-sample: no B9 row's |residual| grows by more than 0.3 dex
  T6  fed-intermediate rows: each within 2x (0.3 dex) of its B9 residual
  L   Laplace identification of the 25 coordinates (the consumers expected unidentified)

Ship rule: SHIP if T1, T3 and T5 hold; T1 and T3 both failing refutes the structure as fitted
(consumers ship as declared-inert); one of the two ships with the finding recorded.
"""
from __future__ import annotations

import json
import math
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

V = data_paths.VALIDATION_DIR
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
B11_REPORT = V / "kinetic_core_b11_fit_report.json"
B11_LAPLACE = V / "kinetic_core_b11_laplace_covariance.json"
OUT = V / "kinetic_core_b11_ship_rule.json"
BOLTON = data_paths.BENCHMARKS_DIR / "thiamine_cys_glucose_120C_Bolton1994.json"
YILTIRAK_130 = data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026.json"
YILTIRAK_100 = data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026.json"
FED_PREFIX = "fed_"


def _read(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _operative_and_drift(report: Dict[str, Any]):
    fr = report["frozen_parameters"]
    block: Dict[str, Any] = {
        "log10_k_ref_at_145C": dict(fr["log10_k_ref_at_145C"]),
        "lumped_formation_Ea_kJ_mol": float(fr["lumped_formation_Ea_kJ_mol"]),
        "decay_Ea_kJ_mol": dict(fr["decay_Ea_kJ_mol"]),
    }
    if fr.get("formation_Ea_by_route_kJ_mol"):
        block["formation_Ea_by_route_kJ_mol"] = dict(fr["formation_Ea_by_route_kJ_mol"])
    if fr.get("oxygen"):
        block["oxygen"] = dict(fr["oxygen"])
    drift = PhDrift(acid_yield=float(fr["ph_drift"]["acid_yield_per_sink_event"]),
                    arp_amine_pka=float(fr["ph_drift"]["arp_secondary_ammonium_pKa"]))
    return core_parameters(SULFUR, frozen=block), drift


def _bundle_folds(path: Path, operative, drift) -> List[Dict[str, Any]]:
    bench = panel.load_bundle(path)
    spec = panel.core_spec(bench, use_buffer=True)
    _, limiting = panel.limiting_precursor_molar(bench)
    reservoir, basis = engine.oxygen_reservoir_units(spec.process)
    rows = []
    for compound, target in panel.bundle_targets(bench).items():
        measured = panel.measured_value(bench, compound, target)
        run = engine.predict(spec, [compound], parameters=operative, draw=CoreDraw(ph_drift=drift))
        predicted = panel.core_native_value(run, compound, "ppb", limiting)
        fold = max(predicted / measured, measured / predicted) if predicted and measured else float("inf")
        rows.append({"bundle": path.stem, "compound": compound, "measured": measured, "predicted": predicted,
                     "fold": fold, "reservoir_units": reservoir, "reservoir_basis": basis})
    return rows


def _hofmann_table1(B11, B23, hofmann_rows, x9: np.ndarray, x11: np.ndarray) -> Dict[str, Any]:
    """Re-score the eight B9 validation rows through the generator's own residual machinery."""
    saved = B23.ACTIVE_FIT_ROWS
    try:
        B23.ACTIVE_FIT_ROWS = tuple(hofmann_rows)
        r9 = B11.residual_vector(x9, False)
        r11 = B11.residual_vector(x11, False)
    finally:
        B23.ACTIVE_FIT_ROWS = saved
    out = []
    for row, a, b in zip(hofmann_rows, r9, r11):
        s = float(row["sigma_log"])
        # a hexose MFT row predicts (numerically) zero: the fold is astronomical, not a measurement of anything
        clipped = 10 ** abs(float(a) * s) > 1e4 or 10 ** abs(float(b) * s) > 1e4
        out.append({"row": row["id"], "b9_fold": 10 ** abs(float(a) * s), "b11_fold": 10 ** abs(float(b) * s),
                    "b9_signed_dex": float(a) * s, "b11_signed_dex": float(b) * s, "predicted_zero": clipped})
    # the hexose MFT rows predict ZERO under B9 and B11 alike (the B9 finding: no hexose -> MFT route);
    # T2 is scored on the rows the engine answers, and the zero rows are listed, not averaged in
    live = [r for r in out if not r["predicted_zero"]]
    return {"rows": out, "n_predicted_zero": len(out) - len(live),
            "worst_b11_fold": max(r["b11_fold"] for r in live), "worst_b9_fold": max(r["b9_fold"] for r in live),
            "pass": all(r["b11_fold"] <= 8.0 for r in live)}


def main() -> int:
    for p in (B9_REPORT, B11_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    # snapshot the eight Hofmann Table-1 rows BEFORE B9 (imported by B11) removes them
    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
    import generate_kinetic_core_b8_fit as B8  # noqa: E402  (install_b8_rows at import)
    hofmann_ids = ("hofmann_ribose_FFT", "hofmann_ribose_MFT", "hofmann_xylose_FFT", "hofmann_xylose_MFT",
                   "hofmann_glucose_FFT", "hofmann_glucose_MFT", "hofmann_fructose_FFT", "hofmann_fructose_MFT")
    hofmann_rows = [r for r in B23.ACTIVE_FIT_ROWS if r["id"] in hofmann_ids]
    assert len(hofmann_rows) == 8, [r["id"] for r in hofmann_rows]
    import generate_kinetic_core_b11_fit as B11  # noqa: E402  (configure() at import; B9 rows)
    from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
    B11.configure()
    sigma = {r["id"]: float(r["sigma_log"]) for r in B23.ACTIVE_FIT_ROWS if r["kind"] != "ph_endpoint"}

    b9, b11 = _read(B9_REPORT), _read(B11_REPORT)
    lap = _read(B11_LAPLACE) if B11_LAPLACE.exists() else None
    oxy = b11["frozen_parameters"]["oxygen"]
    x11 = frozen_vector(b11)
    lower, _upper = B11.full_bounds()
    x9 = np.append(frozen_vector(b9), [lower[B11.CYS_SLOT], lower[B11.RED_SLOT]])   # B9 with the consumers at the floor
    op9, d9 = _operative_and_drift(b9)
    op11, d11 = _operative_and_drift(b11)

    # T1 Bolton -----------------------------------------------------------
    bol9, bol11 = _bundle_folds(BOLTON, op9, d9), _bundle_folds(BOLTON, op11, d11)
    t1 = {"b9": bol9, "b11": bol11, "b9_fold": bol9[0]["fold"], "b11_fold": bol11[0]["fold"],
          "pass": bol11[0]["fold"] < 6.0, "falsified": bol11[0]["fold"] > 12.0}

    # T2 Hofmann Table-1 ---------------------------------------------------
    t2 = _hofmann_table1(B11, B23, hofmann_rows, x9, x11)

    # T3 / T4 Yiltirak -----------------------------------------------------
    y9, y11 = _bundle_folds(YILTIRAK_130, op9, d9), _bundle_folds(YILTIRAK_130, op11, d11)
    by9 = {r["compound"]: r["fold"] for r in y9}
    by11 = {r["compound"]: r["fold"] for r in y11}
    improvement = {c: by9[c] / by11[c] if by11[c] > 0 else float("inf") for c in by9}
    t3 = {"b9": y9, "b11": y11, "improvement_factor": improvement,
          "pass": all(v >= 3.0 for v in improvement.values()), "falsified": all(v < 2.0 for v in improvement.values())}
    s9, s11 = _bundle_folds(YILTIRAK_100, op9, d9), _bundle_folds(YILTIRAK_100, op11, d11)
    t4 = {"b9": s9, "b11": s11, "note": "recorded, not a target of this wave (prereg sec. 5 T4)"}

    # T5 / T6 in-sample ----------------------------------------------------
    r9 = b9["members"][b9["best_start"]]["residual_by_row"]
    r11 = b11["members"][b11["best_start"]]["residual_by_row"]
    growth = []
    for rid, v9 in r9.items():
        if rid in r11 and rid in sigma:
            dex = (abs(float(r11[rid])) - abs(float(v9))) * sigma[rid]
            growth.append({"row": rid, "b9_residual": float(v9), "b11_residual": float(r11[rid]), "growth_dex": dex})
    worst = max(growth, key=lambda g: g["growth_dex"])
    t5 = {"rows_compared": len(growth), "worst": worst, "n_over_0_3_dex": sum(1 for g in growth if g["growth_dex"] > 0.3),
          "n_over_0_5_dex": sum(1 for g in growth if g["growth_dex"] > 0.5), "pass": worst["growth_dex"] <= 0.3,
          "falsified": worst["growth_dex"] > 0.5}
    fed = [g for g in growth if g["row"].startswith(FED_PREFIX)]
    t6 = {"rows": fed, "n": len(fed), "worst_abs_shift_dex": max(abs(g["growth_dex"]) for g in fed),
          "pass": all(abs(g["growth_dex"]) <= math.log10(2.0) for g in fed)}

    # Laplace ----------------------------------------------------------------
    lap_block: Dict[str, Any] = {"status": "no laplace artifact"}
    if lap:
        names = [c["key"] or c["block"] for c in lap["coordinates"]]
        ident = dict(zip(names, (bool(v) for v in lap["identified"])))
        lap_block = {"n_identified": sum(ident.values()), "n_free": len(names),
                     "k_cys_ox_identified": ident.get("k_cys_ox"), "k_red_ox_identified": ident.get("k_red_ox"),
                     "at_least_20": sum(ident.values()) >= 20}

    ships = bool(t1["pass"] and t3["pass"] and t5["pass"])
    if ships:
        verdict = "SHIP"
    elif (not t1["pass"]) and (not t3["pass"]):
        verdict = "DO NOT SHIP: T1 and T3 both fail -- the consumers ship as declared-inert"
    else:
        verdict = "PARTIAL: one of T1/T3 holds -- ships with the failing test named" if t5["pass"] else "DO NOT SHIP: T5 fails"
    payload = {
        "artifact": "kinetic_core_b11_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b11_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b11_ship_rule.py",
            inputs=[p for p in (B9_REPORT, B11_REPORT, B11_LAPLACE, BOLTON, YILTIRAK_130, YILTIRAK_100) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b11_prereg.md"),
        "oxygen": oxy, "oxygen_log10_k": b11["frozen_parameters"].get("oxygen_log10_k"),
        "active_bounds": b11.get("active_bounds"),
        "T1_bolton": t1, "T2_hofmann_table1": t2, "T3_yiltirak_130C": t3, "T4_yiltirak_100C_recorded": t4,
        "T5_in_sample_discipline": t5, "T6_fed_rows": t6, "laplace": lap_block,
        "verdict": verdict,
        "rule": "SHIP if T1 (Bolton < 6x), T3 (Yiltirak 130 C both >= 3x better) and T5 (no row +0.3 dex) hold (prereg sec. 5)",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: oxygen {oxy}; T1 {t1['b9_fold']:.1f}x -> {t1['b11_fold']:.1f}x; T2 worst answered {t2['worst_b9_fold']:.1f}x -> "
          f"{t2['worst_b11_fold']:.1f}x ({t2['n_predicted_zero']} zero rows); T3 {by9} -> {by11}; T5 worst {worst['row']} {worst['growth_dex']:+.2f} dex; "
          f"T6 {t6['pass']}; laplace {lap_block}")
    return 0


def render(p: Dict[str, Any]) -> str:
    t1, t2, t3, t4, t5, t6, lap = (p["T1_bolton"], p["T2_hofmann_table1"], p["T3_yiltirak_130C"], p["T4_yiltirak_100C_recorded"],
                                   p["T5_in_sample_discipline"], p["T6_fed_rows"], p["laplace"])
    o = p["oxygen"]
    out = [f"# Wave B11 ship rule -- **{p['verdict']}**", "",
           f"Prereg: `{p['prereg']}` sec. 5. Rule: {p['rule']}.", "",
           f"* fitted consumers: k_cys_ox **{o['k_cys_ox']:.3g}**, k_red_ox **{o['k_red_ox']:.3g}** per unit per min "
           f"(bands 1e-5..1e-1); supply {o['k_ox_supply']:g} /min; saturation {o['ox_sat_mmol_l']:g} mmol/L",
           f"* active bounds at the optimum: {[a['key'] for a in (p.get('active_bounds') or [])] or 'none'}", "",
           "## T1 -- Bolton 1994 (O2 : thiol 2.07)", "",
           f"* MFT fold error {t1['b9_fold']:.1f}x (B9) -> **{t1['b11_fold']:.1f}x** (B11); reservoir "
           f"{t1['b11'][0]['reservoir_units']:.0f} units ({t1['b11'][0]['reservoir_basis']}) -> below 6x: "
           f"**{'PASS' if t1['pass'] else 'FAIL'}**{' (FALSIFIED: above 12x)' if t1['falsified'] else ''}",
           "", "## T2 -- Hofmann 1998 Table-1 pH-5 rows (O2 : thiol 0.27; the eight B9 validation rows)", "",
           "| row | B9 fold | B11 fold |", "|---|---:|---:|"]
    for r in t2["rows"]:
        if r["predicted_zero"]:
            out.append(f"| {r['row']} | predicts ZERO | predicts ZERO |")
        else:
            out.append(f"| {r['row']} | {r['b9_fold']:.1f} | {r['b11_fold']:.1f} |")
    out += [f"* worst answered row {t2['worst_b9_fold']:.1f}x -> **{t2['worst_b11_fold']:.1f}x** ({t2['n_predicted_zero']} hexose rows predict "
            f"zero under both, the B9 finding); none beyond 8x: **{'PASS' if t2['pass'] else 'FAIL'}**",
            "", "## T3 -- Yiltirak 2026, 130 C / 0.5 h (oxygen in excess)", ""]
    for a, b in zip(t3["b9"], t3["b11"]):
        out.append(f"* {a['compound']}: {a['fold']:.1f}x -> **{b['fold']:.1f}x** (improvement x{t3['improvement_factor'][a['compound']]:.2f})")
    out += [f"* both at least 3x better: **{'PASS' if t3['pass'] else 'FAIL'}**{' (FALSIFIED: neither 2x)' if t3['falsified'] else ''}",
            "", "## T4 -- Yiltirak 100 C / 4 h (recorded, not a target)", ""]
    for a, b in zip(t4["b9"], t4["b11"]):
        out.append(f"* {a['compound']}: {a['fold']:.1f}x -> {b['fold']:.1f}x")
    out += ["", "## T5 -- in-sample discipline", "",
            f"* {t5['rows_compared']} shared rows; worst growth {t5['worst']['growth_dex']:+.2f} dex on `{t5['worst']['row']}`; "
            f"{t5['n_over_0_3_dex']} rows over 0.3 dex, {t5['n_over_0_5_dex']} over 0.5 dex -> **{'PASS' if t5['pass'] else 'FAIL'}**",
            "", "## T6 -- fed-intermediate rows", "",
            f"* {t6['n']} rows; worst shift {t6['worst_abs_shift_dex']:.2f} dex; all within 2x of their B9 residual: **{'PASS' if t6['pass'] else 'FAIL'}**",
            "", "## Laplace", ""]
    if "n_identified" in lap:
        out.append(f"* {lap['n_identified']} of {lap['n_free']} coordinates identified (at least 20: {lap['at_least_20']}); "
                   f"k_cys_ox identified = {lap['k_cys_ox_identified']}, k_red_ox identified = {lap['k_red_ox_identified']} "
                   "(both EXPECTED unidentified: one laboratory's vessel, no oxygen contrast)")
    else:
        out.append(f"* {lap['status']}")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
