#!/usr/bin/env python
"""
Wave B16 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-07).

`results/validation/kinetic_core_b16_prereg.md` sec. 4. Computed from the frozen B9, B16 and B16-lift
artifacts WITHOUT switching the engine (parameters handed to `predict` as an override); writes
`results/validation/kinetic_core_b16_ship_rule.{json,md}`.

  T1  shape: SCH-T-01 agrees under B16 (MFT monotone 30 -> 720 min at 100 C) and the seven fold rows within 0.3 dex
  T2  in-sample: no B9 row's |residual| grows by more than 0.3 dex
  T3  Yiltirak 100 C / 4 h and 110 C / 2 h levels not worse than B9
  T4  the three TTCA rows within 0.1 dex
  T5  Laplace: Ea_decay_thiol_sink identified or off its bound
  T6  recorded: RIB-T-01 / RIB-T-02 (Liu 2023, 168 C) evaluable and their status
Ship rule: SHIP if T1, T2 and T4 hold (variant b16 only; the lift variant is information).
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List

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
B16_REPORT = V / "kinetic_core_b16_fit_report.json"
B16_LIFT_REPORT = V / "kinetic_core_b16_lift_fit_report.json"
B16_LAPLACE = V / "kinetic_core_b16_laplace_covariance.json"
OUT = V / "kinetic_core_b16_ship_rule.json"
YIL = data_paths.MAILLARD_PATH_HOLDOUT_DIR
YIL_STEMS = ("mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026", "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026")
SCH_MINUTES = (30.0, 60.0, 360.0, 720.0)
RIB = {"D-Ribose": 26.1, "L-Cysteine": 326.6}


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def _operative_and_drift(report):
    fr = report["frozen_parameters"]
    block = {"log10_k_ref_at_145C": dict(fr["log10_k_ref_at_145C"]),
             "lumped_formation_Ea_kJ_mol": float(fr["lumped_formation_Ea_kJ_mol"]),
             "decay_Ea_kJ_mol": dict(fr["decay_Ea_kJ_mol"])}
    if fr.get("formation_Ea_by_route_kJ_mol"):
        block["formation_Ea_by_route_kJ_mol"] = dict(fr["formation_Ea_by_route_kJ_mol"])
    if fr.get("oxygen"):
        block["oxygen"] = dict(fr["oxygen"])
    drift = PhDrift(acid_yield=float(fr["ph_drift"]["acid_yield_per_sink_event"]),
                    arp_amine_pka=float(fr["ph_drift"]["arp_secondary_ammonium_pKa"]))
    return core_parameters(SULFUR, frozen=block), drift


def _series(precursors, t_c, minutes, ph, buffer, operative, drift, targets):
    from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram

    out = {}
    for m in minutes:
        spec = FormulationSpec("pot", dict(precursors), ProcessSpec(ThermalProgram.isothermal(float(t_c), float(m)), ph=ph, buffer=buffer))
        run = engine.predict(spec, list(targets), parameters=operative, draw=CoreDraw(ph_drift=drift))
        out[float(m)] = {t: (run.concentrations_ug_per_l.get(t) if run.answered else None) for t in targets}
    return out


def _monotone_up(values: List[float]) -> bool:
    return all(b > a for a, b in zip(values, values[1:]))


def _yiltirak(operative, drift):
    rows = []
    for stem in YIL_STEMS:
        bench = panel.load_bundle(YIL / f"{stem}.json")
        spec = panel.core_spec(bench, use_buffer=True)
        _, limiting = panel.limiting_precursor_molar(bench)
        for compound, target in panel.bundle_targets(bench).items():
            measured = panel.measured_value(bench, compound, target)
            run = engine.predict(spec, [compound], parameters=operative, draw=CoreDraw(ph_drift=drift))
            predicted = panel.core_native_value(run, compound, "ppb", limiting)
            fold = max(predicted / measured, measured / predicted) if predicted and measured else float("inf")
            rows.append({"bundle": stem, "compound": compound, "measured": measured, "predicted": predicted, "fold": fold})
    return rows


def evaluate(report_path: Path, tag: str, b9, sigma, hof_buffer, B16) -> Dict[str, Any]:
    rep = _read(report_path)
    op, drift = _operative_and_drift(rep)
    op9, d9 = _operative_and_drift(b9)
    # T1 ---------------------------------------------------------------
    sch = _series({"D-Ribose": 100.0, "L-Cysteine": 33.0}, 100.0, SCH_MINUTES, 5.0, hof_buffer, op, drift, ("MFT", "FFT"))
    mft = [sch[m]["MFT"] for m in SCH_MINUTES]
    fft = [sch[m]["FFT"] for m in SCH_MINUTES]
    resid = rep["members"][rep["best_start"]]["residual_by_row"]
    folds = {r["id"]: float(resid[r["id"]]) * float(r["sigma_log"]) for r in B16.B16_FIT_ROWS if r["id"] in resid}
    fold_rows = {k: v for k, v in folds.items() if k.startswith("schieberle")}
    ttca_rows = {k: v for k, v in folds.items() if k.startswith("zhai2021")}
    t1 = {"mft_ug_per_l": mft, "fft_ug_per_l": fft, "mft_monotone": _monotone_up(mft), "fft_monotone": _monotone_up(fft),
          "fold_rows_dex": fold_rows, "worst_fold_dex": max(abs(v) for v in fold_rows.values()),
          "pass": _monotone_up(mft) and all(abs(v) <= 0.3 for v in fold_rows.values())}
    # T2 ---------------------------------------------------------------
    r9 = b9["members"][b9["best_start"]]["residual_by_row"]
    growth = []
    for rid, v9 in r9.items():
        if rid in resid and rid in sigma:
            growth.append({"row": rid, "b9_residual": float(v9), "b16_residual": float(resid[rid]),
                           "growth_dex": (abs(float(resid[rid])) - abs(float(v9))) * sigma[rid]})
    worst = max(growth, key=lambda g: g["growth_dex"])
    t2 = {"rows_compared": len(growth), "worst": worst, "n_over_0_3_dex": sum(1 for g in growth if g["growth_dex"] > 0.3),
          "n_over_0_5_dex": sum(1 for g in growth if g["growth_dex"] > 0.5), "pass": worst["growth_dex"] <= 0.3}
    # T3 ---------------------------------------------------------------
    y9, y16 = _yiltirak(op9, d9), _yiltirak(op, drift)
    t3 = {"b9": y9, "b16": y16, "b9_median_fold": statistics.median(r["fold"] for r in y9),
          "b16_median_fold": statistics.median(r["fold"] for r in y16)}
    t3["pass"] = t3["b16_median_fold"] <= t3["b9_median_fold"] * 1.05
    # T4 ---------------------------------------------------------------
    t4 = {"rows_dex": ttca_rows, "worst_dex": max(abs(v) for v in ttca_rows.values()), "pass": all(abs(v) <= 0.1 for v in ttca_rows.values())}
    # T6 ---------------------------------------------------------------
    from src.kinetic_core.ph_state import BufferSpec

    rib_buf = BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="Liu 2023 sec. 2.2")
    rib = _series(RIB, 168.0, (20.0, 40.0, 60.0), 5.88, rib_buf, op, drift, ("MFT", "FFT"))
    rmft = [rib[m]["MFT"] for m in (20.0, 40.0, 60.0)]
    rfft = [rib[m]["FFT"] for m in (20.0, 60.0)]
    evaluable = all(v is not None and v > 1e-6 for v in rmft + rfft)
    t6 = {"mft_ug_per_l": rmft, "fft_20_60": rfft, "evaluable": evaluable,
          "rib_t_01_decreasing": evaluable and rmft[0] > rmft[1] > rmft[2],
          "rib_t_02_flat": evaluable and abs(math.log10(rfft[1] / rfft[0])) <= math.log10(1.05)}
    fr = rep["frozen_parameters"]
    return {"variant": tag, "thiol_sink_ea": fr["decay_Ea_kJ_mol"].get("thiol_sink"), "ceiling": rep.get("thiol_sink_ceiling_kj_mol"),
            "cost": rep["objective"]["final_cost"], "active_bounds": [a["key"] for a in rep.get("active_bounds") or []],
            "T1_shape": t1, "T2_in_sample": t2, "T3_yiltirak_levels": t3, "T4_ttca_rows": t4, "T6_liu2023_168C": t6}


def main() -> int:
    for p in (B9_REPORT, B16_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    import generate_kinetic_core_b16_fit as B16  # noqa: E402  (configure(False) at import)
    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
    B16.configure(False)
    sigma = {r["id"]: float(r["sigma_log"]) for r in B23.ACTIVE_FIT_ROWS if r["kind"] != "ph_endpoint"}
    b9 = _read(B9_REPORT)
    hof_buffer = B23.BUFFER_HOFMANN
    main_v = evaluate(B16_REPORT, "b16", b9, sigma, hof_buffer, B16)
    lift_v = evaluate(B16_LIFT_REPORT, "b16_lift", b9, sigma, hof_buffer, B16) if B16_LIFT_REPORT.exists() else None
    lap = _read(B16_LAPLACE) if B16_LAPLACE.exists() else None
    t5: Dict[str, Any] = {"status": "no laplace artifact"}
    if lap:
        names = [c["key"] or c["block"] for c in lap["coordinates"]]
        ident = dict(zip(names, (bool(v) for v in lap["identified"])))
        sig = dict(zip(names, lap["sigma"]))
        off_bound = "Ea_decay_thiol_sink" not in main_v["active_bounds"]
        t5 = {"thiol_sink_identified": ident.get("thiol_sink"), "thiol_sink_sigma": sig.get("thiol_sink"),
              "off_bound": off_bound, "n_identified": sum(ident.values()), "n_free": len(names),
              "pass": bool(ident.get("thiol_sink")) or off_bound}
    ships = bool(main_v["T1_shape"]["pass"] and main_v["T2_in_sample"]["pass"] and main_v["T4_ttca_rows"]["pass"])
    if ships:
        verdict = "SHIP"
    elif lift_v and lift_v["T1_shape"]["pass"] and not main_v["T1_shape"]["pass"]:
        verdict = "DO NOT SHIP: the 100 C series is reproduced only with the thiol-sink ceiling lifted -- the owner's decision"
    else:
        verdict = "DO NOT SHIP"
    payload = {
        "artifact": "kinetic_core_b16_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b16_ship_rule",
                                                  generated_by="scripts/generators/generate_kinetic_core_b16_ship_rule.py",
                                                  inputs=[p for p in (B9_REPORT, B16_REPORT, B16_LIFT_REPORT, B16_LAPLACE) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b16_prereg.md"),
        "b16": main_v, "b16_lift": lift_v, "T5_laplace": t5, "verdict": verdict,
        "rule": "SHIP if T1 (shape + folds within 0.3 dex), T2 (no B9 row +0.3 dex) and T4 (TTCA rows within 0.1 dex) hold on variant b16",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    m = main_v
    print(f"{verdict}: b16 thiol-sink Ea {m['thiol_sink_ea']:.1f}; T1 {m['T1_shape']['pass']} (MFT {[round(v, 1) for v in m['T1_shape']['mft_ug_per_l']]}, "
          f"worst fold {m['T1_shape']['worst_fold_dex']:.2f} dex); T2 {m['T2_in_sample']['pass']} (worst {m['T2_in_sample']['worst']['row']} "
          f"{m['T2_in_sample']['worst']['growth_dex']:+.2f}); T3 {m['T3_yiltirak_levels']['b9_median_fold']:.0f}x -> {m['T3_yiltirak_levels']['b16_median_fold']:.0f}x; "
          f"T4 {m['T4_ttca_rows']['pass']} (worst {m['T4_ttca_rows']['worst_dex']:.2f}); T6 {m['T6_liu2023_168C']}")
    if lift_v:
        l = lift_v
        print(f"  lift: Ea {l['thiol_sink_ea']:.1f}; T1 {l['T1_shape']['pass']} (MFT {[round(v, 1) for v in l['T1_shape']['mft_ug_per_l']]}); "
              f"T2 {l['T2_in_sample']['pass']} ({l['T2_in_sample']['worst']['growth_dex']:+.2f}); T4 {l['T4_ttca_rows']['pass']}; T6 {l['T6_liu2023_168C']}")
    return 0


def _variant_md(v: Dict[str, Any]) -> List[str]:
    t1, t2, t3, t4, t6 = v["T1_shape"], v["T2_in_sample"], v["T3_yiltirak_levels"], v["T4_ttca_rows"], v["T6_liu2023_168C"]
    out = [f"## Variant `{v['variant']}` -- thiol-sink Ea {v['thiol_sink_ea']:.1f} kJ/mol (ceiling {v['ceiling']:.0f}); cost {v['cost']:.2f}; "
           f"active bounds {v['active_bounds'] or 'none'}", "",
           f"* **T1 shape:** MFT at 100 C over 30/60/360/720 min = {[round(x, 1) for x in t1['mft_ug_per_l']]} ug/L "
           f"(monotone: {t1['mft_monotone']}); FFT {[round(x, 1) for x in t1['fft_ug_per_l']]} (monotone: {t1['fft_monotone']}); "
           f"fold rows worst {t1['worst_fold_dex']:.2f} dex -> **{'PASS' if t1['pass'] else 'FAIL'}**",
           "  * " + ", ".join(f"{k.replace('schieberle_', '')} {v:+.2f}" for k, v in t1["fold_rows_dex"].items()),
           f"* **T2 in-sample:** worst growth {t2['worst']['growth_dex']:+.2f} dex on `{t2['worst']['row']}`; {t2['n_over_0_3_dex']} rows over 0.3, "
           f"{t2['n_over_0_5_dex']} over 0.5 -> **{'PASS' if t2['pass'] else 'FAIL'}**",
           f"* **T3 Yiltirak 100/110 C levels:** median fold {t3['b9_median_fold']:.1f}x (B9) -> {t3['b16_median_fold']:.1f}x -> {'PASS' if t3['pass'] else 'FAIL'}",
           f"* **T4 TTCA rows:** " + ", ".join(f"{k.split('_')[3]} {v:+.2f}" for k, v in t4["rows_dex"].items()) + f" dex -> **{'PASS' if t4['pass'] else 'FAIL'}**",
           f"* **T6 Liu 2023 at 168 C (recorded):** MFT {[('%.3g' % x) if x is not None else None for x in t6['mft_ug_per_l']]} ug/L at 20/40/60 min; "
           f"evaluable {t6['evaluable']}; RIB-T-01 decreasing {t6['rib_t_01_decreasing']}; RIB-T-02 flat {t6['rib_t_02_flat']}", ""]
    return out


def render(p: Dict[str, Any]) -> str:
    out = [f"# Wave B16 ship rule -- **{p['verdict']}**", "", f"Prereg: `{p['prereg']}` sec. 4. Rule: {p['rule']}.", ""]
    out += _variant_md(p["b16"])
    if p.get("b16_lift"):
        out += _variant_md(p["b16_lift"])
    t5 = p["T5_laplace"]
    out += ["## T5 -- Laplace (variant b16)", ""]
    if "pass" in t5:
        out.append(f"* thiol-sink barrier identified = {t5['thiol_sink_identified']}, sigma = {t5['thiol_sink_sigma']}, off its bound = {t5['off_bound']}; "
                   f"{t5['n_identified']} of {t5['n_free']} identified -> **{'PASS' if t5['pass'] else 'FAIL'}**")
    else:
        out.append(f"* {t5['status']}")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
