#!/usr/bin/env python
"""
Wave B17 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-08).

`results/validation/kinetic_core_b17_prereg.md` sec. 4, computed from the frozen B9 and B17 artifacts
WITHOUT switching the engine (the candidate vector is handed to `core_parameters(SULFUR, frozen=...)`),
written to `results/validation/kinetic_core_b17_ship_rule.{json,md}`.

  T1  the reference pot at 100 C: MFT and FFT at 360 and 720 min within 0.5 dex of Schieberle 2000 Table IV
      as within-study ratios to the 30 min point, and both still rising between 6 and 12 h
  T2  the 145 C fed pots do not break: every B9 fit row within 0.3 dex of its B9 residual
  T3  the dimer shares: Zhou 2023 (ARP 20 + Cys 20 mM, water, pH 6/7/8, 120 C / 60 min) and Zhang 2024
      (thiamine + xylose + cysteine, pH 4.9 phosphate, 115 C / 60 min) thiol-equivalents-in-dimer over free
      thiol within 0.3 dex
  T4  Yiltirak 2026's four pots: median fold error below 20 (B9: 115)
  T5  Wang 2022's 140 C pot: MFT and FFT decline from their peak to 180 min by less than one decade
  T6  identification: Laplace sigma on log10 k_dimer_release below one decade, off its bound; the cost slice
      along the coordinate is not bound-limited
Ship rule: SHIP if T1, T2 and T6 hold; T3 to T5 are reported.
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
from src.kinetic_core.engine import SULFUR, CoreDraw, FormulationSpec, ProcessSpec, ThermalProgram, core_parameters  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec, PhDrift  # noqa: E402

V = data_paths.VALIDATION_DIR
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
B17_REPORT = V / "kinetic_core_b17_fit_report.json"
B17_LAPLACE = V / "kinetic_core_b17_laplace_covariance.json"
OUT = V / "kinetic_core_b17_ship_rule.json"
YIL = data_paths.MAILLARD_PATH_HOLDOUT_DIR
YIL_STEMS = tuple(sorted(p.stem for p in YIL.glob("mp_holdout_ribose_cysteine_buffer_*_Yiltirak2026.json")))
SCH_MINUTES = (30.0, 60.0, 360.0, 720.0)
TABLE_IV = {"MFT": {30.0: 4.5, 60.0: 13.8, 360.0: 156.0, 720.0: 179.0}, "FFT": {30.0: 2.0, 60.0: 3.1, 360.0: 110.0, 720.0: 132.0}}
#: Zhou 2023 sec. 2.1: thiol-equivalents in the dimer over the free monomer, percent, pH 6 / 7 / 8
ZHOU_SHARE = {"MFT": {6.0: 8.6, 7.0: 6.5, 8.0: 9.6}, "FFT": {6.0: 4.4, 7.0: 2.4, 8.0: 0.0}}
#: Zhang 2024 Table (Cys arm): MFT-MFT / MFT mass ratio 0.086 -> molar dimer/monomer x 114.17/226.32, x2 thiol equivalents
ZHANG_SHARE_PCT = 0.086 * 114.17 / 226.32 * 2.0 * 100.0
MW = {"MFT": 114.17, "FFT": 114.17, "MFT dimer": 226.32}


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def _operative_and_drift(report, with_release: bool):
    fr = report["frozen_parameters"]
    block = {"log10_k_ref_at_145C": dict(fr["log10_k_ref_at_145C"]),
             "lumped_formation_Ea_kJ_mol": float(fr["lumped_formation_Ea_kJ_mol"]),
             "decay_Ea_kJ_mol": dict(fr["decay_Ea_kJ_mol"])}
    if fr.get("formation_Ea_by_route_kJ_mol"):
        block["formation_Ea_by_route_kJ_mol"] = dict(fr["formation_Ea_by_route_kJ_mol"])
    if fr.get("oxygen"):
        block["oxygen"] = dict(fr["oxygen"])
    if with_release and fr.get("dimer_release_log10_k"):
        block["dimer_release_log10_k"] = dict(fr["dimer_release_log10_k"])
    drift = PhDrift(acid_yield=float(fr["ph_drift"]["acid_yield_per_sink_event"]),
                    arp_amine_pka=float(fr["ph_drift"]["arp_secondary_ammonium_pKa"]))
    return core_parameters(SULFUR, frozen=block), drift


def _predict(precursors, t_c, minutes, ph, buffer, operative, drift, targets):
    spec = FormulationSpec("pot", dict(precursors), ProcessSpec(ThermalProgram.isothermal(float(t_c), float(minutes)), ph=ph, buffer=buffer))
    run = engine.predict(spec, list(targets), parameters=operative, draw=CoreDraw(ph_drift=drift))
    return {t: (run.concentrations_ug_per_l.get(t) if run.answered else None) for t in targets}, run


def _series(precursors, t_c, minutes, ph, buffer, operative, drift, targets):
    return {float(m): _predict(precursors, t_c, m, ph, buffer, operative, drift, targets)[0] for m in minutes}


def t1_reference_pot(op, drift, hof_buffer):
    sch = _series({"D-Ribose": 100.0, "L-Cysteine": 33.0}, 100.0, SCH_MINUTES, 5.0, hof_buffer, op, drift, ("MFT", "FFT"))
    out = {"mft_ug_per_l": [sch[m]["MFT"] for m in SCH_MINUTES], "fft_ug_per_l": [sch[m]["FFT"] for m in SCH_MINUTES], "ratios_dex": {}}
    ok = True
    for s in ("MFT", "FFT"):
        for m in (360.0, 720.0):
            model = sch[m][s] / sch[30.0][s] if sch[30.0][s] else float("inf")
            table = TABLE_IV[s][m] / TABLE_IV[s][30.0]
            d = math.log10(model / table) if model > 0 and math.isfinite(model) else float("inf")
            out["ratios_dex"][f"{s}_{int(m)}_over_30"] = d
            ok = ok and abs(d) <= 0.5
    rising = {s: bool(sch[720.0][s] > sch[360.0][s]) for s in ("MFT", "FFT")}
    out["still_rising_6_to_12_h"] = rising
    out["pass"] = bool(ok and all(rising.values()))
    return out


def t2_b9_rows(rep, b9, sigma):
    resid = rep["members"][rep["best_start"]]["residual_by_row"]
    r9 = b9["members"][b9["best_start"]]["residual_by_row"]
    growth = []
    for rid, v9 in r9.items():
        if rid in resid and rid in sigma:
            growth.append({"row": rid, "b9_residual": float(v9), "b17_residual": float(resid[rid]),
                           "growth_dex": (abs(float(resid[rid])) - abs(float(v9))) * sigma[rid]})
    worst = max(growth, key=lambda g: g["growth_dex"])
    return {"rows_compared": len(growth), "worst": worst, "n_over_0_3_dex": sum(1 for g in growth if g["growth_dex"] > 0.3),
            "pass": worst["growth_dex"] <= 0.3}


def _share_pct(conc, thiol, dimer):
    """Thiol equivalents held in the dimer over the free thiol, percent (2 mol thiol per mol dimer)."""
    t, d = conc.get(thiol), conc.get(dimer)
    if not t or d is None:
        return None
    return 100.0 * 2.0 * (d / MW[dimer]) / (t / MW[thiol])


def t3_dimer_shares(op, drift):
    out = {"zhou2023": {}, "zhang2024": {}}
    worst = 0.0
    water = BufferSpec(kind="none", declared=True, source="Zhou 2023 sec. 2: deionized water, no buffer, initial pH by NaOH")
    for ph in (6.0, 7.0, 8.0):
        try:
            conc, run = _predict({"ARP": 20.0, "L-Cysteine": 20.0}, 120.0, 60.0, ph, water, op, drift, ("MFT", "MFT dimer"))
        except Exception as exc:  # pragma: no cover
            out["zhou2023"][str(ph)] = {"error": str(exc)}
            continue
        mft_share = _share_pct(conc, "MFT", "MFT dimer")
        row = {"model_mft_share_pct": mft_share, "zhou_mft_share_pct": ZHOU_SHARE["MFT"][ph], "answered": run.answered}
        if mft_share and mft_share > 0:
            row["dex"] = math.log10(mft_share / ZHOU_SHARE["MFT"][ph])
            worst = max(worst, abs(row["dex"]))
        out["zhou2023"][str(ph)] = row
    phosphate = BufferSpec(kind="phosphate", phosphate_mol_l=0.1, declared=False, source="Zhang 2024: 'pH 4.9 phosphate buffered solution', molarity unstated; 0.1 M assumed")
    conc, run = _predict({"thiamine": 44.5, "D-Xylose": 99.9, "L-Cysteine": 123.8}, 115.0, 60.0, 4.9, phosphate, op, drift, ("MFT", "MFT dimer"))
    share = _share_pct(conc, "MFT", "MFT dimer")
    out["zhang2024"] = {"model_mft_share_pct": share, "zhang_mft_share_pct": ZHANG_SHARE_PCT, "answered": run.answered,
                        "note": "the cysteine arm; methionine (the paper's methanethiol source) is not a precursor the core charges"}
    if share and share > 0:
        out["zhang2024"]["dex"] = math.log10(share / ZHANG_SHARE_PCT)
        worst = max(worst, abs(out["zhang2024"]["dex"]))
    out["worst_dex"] = worst
    out["pass"] = bool(worst <= 0.3 and all("dex" in r for r in out["zhou2023"].values()) and "dex" in out["zhang2024"])
    return out


def t4_yiltirak(op, drift):
    rows = []
    for stem in YIL_STEMS:
        bench = panel.load_bundle(YIL / f"{stem}.json")
        spec = panel.core_spec(bench, use_buffer=True)
        _, limiting = panel.limiting_precursor_molar(bench)
        for compound, target in panel.bundle_targets(bench).items():
            measured = panel.measured_value(bench, compound, target)
            run = engine.predict(spec, [compound], parameters=op, draw=CoreDraw(ph_drift=drift))
            predicted = panel.core_native_value(run, compound, "ppb", limiting)
            fold = max(predicted / measured, measured / predicted) if predicted and measured else float("inf")
            rows.append({"bundle": stem, "compound": compound, "measured": measured, "predicted": predicted, "fold": fold})
    med = statistics.median(r["fold"] for r in rows)
    return {"rows": rows, "median_fold": med, "pass": bool(med < 20.0)}


def t5_wang_shape(op, drift):
    buf = BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Wang 2022 sec. 2: sodium phosphate buffer pH 5.5, 0.2 mol/L")
    minutes = (30.0, 60.0, 90.0, 120.0, 150.0, 180.0)
    ser = _series({"L-Cysteine": 400.0, "D-Xylose": 400.0}, 140.0, minutes, 5.5, buf, op, drift, ("MFT", "FFT"))
    out = {"minutes": list(minutes)}
    ok = True
    for s in ("MFT", "FFT"):
        vals = [ser[m][s] for m in minutes]
        peak = max(vals)
        decline = math.log10(peak / vals[-1]) if vals[-1] and peak else float("inf")
        out[s] = {"ug_per_l": vals, "peak_min": minutes[int(np.argmax(vals))], "decline_from_peak_dex": decline}
        ok = ok and decline < 1.0
    out["pass"] = bool(ok)
    return out


def t6_identification(rep, lap, B17):
    fr = rep["frozen_parameters"]
    x = float(fr["dimer_release_log10_k"]["k_dimer_release"])
    lo, hi = rep["release_band_log10k"]
    on_bound = bool(x - lo <= 1e-3 * (hi - lo) or hi - x <= 1e-3 * (hi - lo))
    out: Dict[str, Any] = {"log10_k_dimer_release": x, "band": [lo, hi], "on_bound": on_bound}
    if lap:
        names = [c["key"] or c["block"] for c in lap["coordinates"]]
        sig = dict(zip(names, lap["sigma"]))
        ident = dict(zip(names, (bool(v) for v in lap["identified"])))
        out["sigma_dex"] = sig.get("k_dimer_release")
        out["identified"] = ident.get("k_dimer_release")
        out["n_identified"] = sum(ident.values())
        out["n_free"] = len(names)
    # the cost slice along the release coordinate (quick mode), the profile the prereg asks for
    from generate_kinetic_core_b8_laplace import frozen_vector

    x_opt = frozen_vector(rep)
    slice_pts = {}
    for d in (-1.0, -0.5, 0.0, 0.5, 1.0):
        xx = x_opt.copy()
        xx[B17.RELEASE_SLOT] = min(max(x + d, lo), hi)
        r = B17.residual_vector(xx, True)
        slice_pts[f"{d:+.1f}"] = 0.5 * float(np.dot(r, r))
    c0 = slice_pts["+0.0"]
    out["cost_slice"] = slice_pts
    out["slice_verdict"] = ("bound_limited" if on_bound else
                            ("flat" if max(slice_pts.values()) - c0 < 0.5 else
                             ("quadratic" if slice_pts["-1.0"] > c0 and slice_pts["+1.0"] > c0 else "shifted")))
    out["pass"] = bool((out.get("sigma_dex") is not None and out["sigma_dex"] < 1.0) and not on_bound and out["slice_verdict"] not in ("bound_limited", "flat"))
    return out


def render(p: Dict[str, Any]) -> str:
    t1, t2, t3, t4, t5, t6 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"], p["T6"]
    L = [f"# Wave B17 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         f"log10 k_dimer_release at 145 C: {t6['log10_k_dimer_release']:.3f} (band {t6['band']}); cost {p['cost']:.2f}; active bounds {p['active_bounds'] or 'none'}", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 reference pot 100 C | MFT {[round(v, 1) for v in t1['mft_ug_per_l']]} ug/L at 30/60/360/720 min; ratios (dex) {t1['ratios_dex']}; rising 6->12 h {t1['still_rising_6_to_12_h']} | {t1['pass']} |",
         f"| T2 B9 rows | worst {t2['worst']['row']} {t2['worst']['growth_dex']:+.2f} dex; {t2['n_over_0_3_dex']} rows over 0.3 | {t2['pass']} |",
         f"| T3 dimer shares | Zhou 2023: {{ {', '.join(f'pH {k}: model {v.get('model_mft_share_pct', float('nan')):.2f} % vs {v.get('zhou_mft_share_pct')} %' for k, v in t3['zhou2023'].items())} }}; Zhang 2024 Cys arm: model {t3['zhang2024'].get('model_mft_share_pct', float('nan')):.2f} % vs {t3['zhang2024']['zhang_mft_share_pct']:.1f} %; worst {t3['worst_dex']:.2f} dex | {t3['pass']} |",
         f"| T4 Yiltirak | median fold {t4['median_fold']:.1f} (B9 115; B9 under this rule: {p['b9_reference']['T4_median_fold']:.1f}) | {t4['pass']} |",
         f"| T5 Wang 2022 140 C shape | MFT peak {t5['MFT']['peak_min']:.0f} min, decline {t5['MFT']['decline_from_peak_dex']:.2f} dex; FFT peak {t5['FFT']['peak_min']:.0f} min, decline {t5['FFT']['decline_from_peak_dex']:.2f} dex | {t5['pass']} |",
         f"| T6 identification | sigma {t6.get('sigma_dex')} dex, identified {t6.get('identified')}, on bound {t6['on_bound']}, slice {t6['slice_verdict']} {t6['cost_slice']} | {t6['pass']} |",
         "", "## The reference pot under B9, for comparison", "",
         f"- MFT {[round(v, 1) for v in p['b9_reference']['T1']['mft_ug_per_l']]} ug/L; ratios {p['b9_reference']['T1']['ratios_dex']}; rising {p['b9_reference']['T1']['still_rising_6_to_12_h']}",
         f"- Wang 140 C under B9: MFT decline {p['b9_reference']['T5']['MFT']['decline_from_peak_dex']:.2f} dex, FFT {p['b9_reference']['T5']['FFT']['decline_from_peak_dex']:.2f} dex",
         f"- dimer shares under B9: Zhou {{ {', '.join(f'pH {k}: {v.get('model_mft_share_pct', float('nan')):.2f} %' for k, v in p['b9_reference']['T3']['zhou2023'].items())} }}, Zhang {p['b9_reference']['T3']['zhang2024'].get('model_mft_share_pct', float('nan')):.2f} %"]
    return "\n".join(L) + "\n"


def main() -> int:
    for p in (B9_REPORT, B17_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    import generate_kinetic_core_b17_fit as B17  # noqa: E402  (configure() at import)
    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
    sigma = {r["id"]: float(r["sigma_log"]) for r in B23.ACTIVE_FIT_ROWS if r["kind"] != "ph_endpoint"}
    b9, rep = _read(B9_REPORT), _read(B17_REPORT)
    lap = _read(B17_LAPLACE) if B17_LAPLACE.exists() else None
    hof_buffer = B23.BUFFER_HOFMANN
    op, drift = _operative_and_drift(rep, True)
    op9, d9 = _operative_and_drift(b9, False)
    T1 = t1_reference_pot(op, drift, hof_buffer)
    T2 = t2_b9_rows(rep, b9, sigma)
    T3 = t3_dimer_shares(op, drift)
    T4 = t4_yiltirak(op, drift)
    T5 = t5_wang_shape(op, drift)
    T6 = t6_identification(rep, lap, B17)
    ref = {"T1": t1_reference_pot(op9, d9, hof_buffer), "T3": t3_dimer_shares(op9, d9), "T4_median_fold": t4_yiltirak(op9, d9)["median_fold"], "T5": t5_wang_shape(op9, d9)}
    ships = bool(T1["pass"] and T2["pass"] and T6["pass"])
    payload = {
        "artifact": "kinetic_core_b17_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b17_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b17_ship_rule.py",
                                                  inputs=[p for p in (B9_REPORT, B17_REPORT, B17_LAPLACE) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b17_prereg.md"),
        "rule": "SHIP if T1 (reference pot ratios within 0.5 dex and still rising 6 -> 12 h), T2 (no B9 row +0.3 dex) and T6 (release identified, off its bound, slice not flat) hold; T3-T5 reported",
        "cost": rep["objective"]["final_cost"], "active_bounds": [a["key"] for a in rep.get("active_bounds") or []],
        "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5, "T6": T6, "b9_reference": ref,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], T1["ratios_dex"], T1["still_rising_6_to_12_h"], "| T2", T2["pass"], T2["worst"]["growth_dex"],
          "| T3", T3["pass"], T3["worst_dex"], "| T4", T4["median_fold"], "| T5", T5["pass"], "| T6", T6["pass"], T6.get("sigma_dex"), T6["slice_verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
