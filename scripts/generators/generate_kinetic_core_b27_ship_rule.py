#!/usr/bin/env python
"""
Wave B27 -- THE PRE-REGISTERED SHIP RULE, evaluated BEFORE ANY CONSTANT WAS FITTED (2026-09-11).

`results/validation/kinetic_core_b27_prereg.md` sec. 5 says SHIP if T1, T2, T3 and T6 hold. This generator asks,
at the shipped B9 vector with the one new coordinate phi swept over its whole band, whether T3 CAN hold -- and it
cannot, for a reason no fit can change. So the fit was not run (the B30 discipline), and this artifact is the
record of why, written to `results/validation/kinetic_core_b27_ship_rule.{json,md}`.

  G1  T3 REACHABILITY. Zhou 2023's and Zhang 2024's pots carry the ambient oxidant, 1.0. The dimer step is first
      order in it. How much of that pool do the consumers actually use, and how much could the whole
      mercaptoketone flux add at phi = 1? A decade on the dimer share needs the pool to rise about tenfold.
  G2  THE WHITFIELD POT. The one pot whose oxidant is genuinely zero: the MFT disulfide share against phi, and
      the phi at which the 35 % floor is first reached. Prereg sec. 10: phi pinned at its ceiling of 1 is
      evidence AGAINST the structure.
  G3  FIT-FREE T2. The residual on every row at phi = 1 against phi -> 0, in dex; the Kumazawa four separately
      (prereg sec. 9: they must move by under 0.05 dex).
  T1, T5  the reference pot and the Wang shape, at phi = 1 (reported).
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
from src.kinetic_core import engine  # noqa: E402
from src.kinetic_core.engine import SULFUR, CoreDraw, FormulationSpec, ProcessSpec, ThermalProgram, core_parameters  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec, PhDrift  # noqa: E402
from src.kinetic_core.sulfur import integrate_sulfur, sulfur_flux_budget  # noqa: E402

V = data_paths.VALIDATION_DIR
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
PREREG = V / "kinetic_core_b27_prereg.md"
OUT = V / "kinetic_core_b27_ship_rule.json"
SCH_MINUTES = (30.0, 60.0, 360.0, 720.0)
TABLE_IV = {"MFT": {30.0: 4.5, 60.0: 13.8, 360.0: 156.0, 720.0: 179.0}, "FFT": {30.0: 2.0, 60.0: 3.1, 360.0: 110.0, 720.0: 132.0}}
AMBIENT_POTS = ("zhou_arp_cys_pH7", "zhang_fig1_cys", "zhang_fig1_gcys")
OX_CONSUMERS = ("ch_dimer_mft", "ch_dimer_fft", "ch_cys_ox", "ch_red_ox_nf", "ch_red_ox_dpo")
MERCAPTOKETONE = ("r_nf_mp3p", "r_mgo_mp", "ch_redox_mp3p")
PHI_GRID = (1e-4, 1e-2, 0.1, 0.2, 0.316, 0.5, 0.708, 1.0)
KUMAZAWA = ("kumazawa_FFT_retention_pH5_0", "kumazawa_FFT_retention_pH5_4", "kumazawa_FFT_retention_pH6_0", "kumazawa_FFT_retention_pH6_4")
CELSIUS = 273.15


def _operative_and_drift(report, log10_phi):
    fr = report["frozen_parameters"]
    block = {"log10_k_ref_at_145C": dict(fr["log10_k_ref_at_145C"]),
             "lumped_formation_Ea_kJ_mol": float(fr["lumped_formation_Ea_kJ_mol"]),
             "decay_Ea_kJ_mol": dict(fr["decay_Ea_kJ_mol"])}
    if fr.get("formation_Ea_by_route_kJ_mol"):
        block["formation_Ea_by_route_kJ_mol"] = dict(fr["formation_Ea_by_route_kJ_mol"])
    if fr.get("oxygen"):
        block["oxygen"] = dict(fr["oxygen"])
    if log10_phi is not None:
        block["dicarbonyl_redox"] = {"log10_ox_yield_per_mercaptoketone": float(log10_phi)}
    drift = PhDrift(acid_yield=float(fr["ph_drift"]["acid_yield_per_sink_event"]),
                    arp_amine_pka=float(fr["ph_drift"]["arp_secondary_ammonium_pKa"]))
    return core_parameters(SULFUR, frozen=block), drift


def _predict(precursors, t_c, minutes, ph, buffer, operative, drift, targets):
    spec = FormulationSpec("pot", dict(precursors), ProcessSpec(ThermalProgram.isothermal(float(t_c), float(minutes)), ph=ph, buffer=buffer))
    run = engine.predict(spec, list(targets), parameters=operative, draw=CoreDraw(ph_drift=drift))
    return {t: (run.concentrations_ug_per_l.get(t) if run.answered else None) for t in targets}


def t1_reference_pot(op, drift, hof_buffer):
    sch = {m: _predict({"D-Ribose": 100.0, "L-Cysteine": 33.0}, 100.0, m, 5.0, hof_buffer, op, drift, ("MFT", "FFT")) for m in SCH_MINUTES}
    out = {"mft_ug_per_l": [sch[m]["MFT"] for m in SCH_MINUTES], "fft_ug_per_l": [sch[m]["FFT"] for m in SCH_MINUTES], "ratios_dex": {}}
    ok = True
    for s in ("MFT", "FFT"):
        for m in (360.0, 720.0):
            model = sch[m][s] / sch[30.0][s] if sch[30.0][s] else float("inf")
            d = math.log10(model / (TABLE_IV[s][m] / TABLE_IV[s][30.0])) if model > 0 and math.isfinite(model) else float("inf")
            out["ratios_dex"][f"{s}_{int(m)}_over_30"] = d
            ok = ok and abs(d) <= 0.5
    rising = {s: bool(sch[720.0][s] > sch[360.0][s]) for s in ("MFT", "FFT")}
    out["still_rising_6_to_12_h"] = rising
    out["pass"] = bool(ok and all(rising.values()))
    return out


def t5_wang_shape(op, drift):
    buf = BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Wang 2022 sec. 2: sodium phosphate buffer pH 5.5, 0.2 mol/L")
    minutes = (30.0, 60.0, 90.0, 120.0, 150.0, 180.0)
    ser = {m: _predict({"L-Cysteine": 400.0, "D-Xylose": 400.0}, 140.0, m, 5.5, buf, op, drift, ("MFT", "FFT")) for m in minutes}
    out: Dict[str, Any] = {"minutes": list(minutes)}
    ok = True
    for s in ("MFT", "FFT"):
        vals = [ser[m][s] for m in minutes]
        peak = max(vals)
        decline = math.log10(peak / vals[-1]) if vals[-1] and peak else float("inf")
        out[s] = {"ug_per_l": vals, "peak_min": minutes[int(np.argmax(vals))], "decline_from_peak_dex": decline}
        ok = ok and decline < 1.0
    out["pass"] = bool(ok)
    return out


def g1_t3_reachability(B27, B23, x) -> Dict[str, Any]:
    p = B27.build_parameters(x)
    _f, _e, _d, drift = B23.unpack(x)
    pots = {}
    for name in AMBIENT_POTS:
        s = B23.SYSTEMS[name]
        fl = sulfur_flux_budget(p, s["t_c"] + CELSIUS, s["initial"], s["minutes"], ph=s["ph"], buffer_spec=s.get("buffer"), ph_drift=drift)
        ox0 = float(s["initial"].get("OX", 0.0))
        used = sum(float(fl.get(k, 0.0)) for k in OX_CONSUMERS)
        supply = sum(float(fl.get(k, 0.0)) for k in MERCAPTOKETONE)
        pots[name] = {"ox_charged": ox0, "ox_consumed": used, "consumed_fraction": used / ox0 if ox0 else None,
                      "mercaptoketone_flux_mmol_l": supply, "supply_at_phi_1_over_pool": supply / ox0 if ox0 else None,
                      "dimer_flux": float(fl.get("ch_dimer_mft", 0.0)) + float(fl.get("ch_dimer_fft", 0.0))}
    # A decade on a first-order-in-OX share needs the pool to rise about tenfold: supply / pool >= 9.
    best = max(v["supply_at_phi_1_over_pool"] or 0.0 for v in pots.values())
    return {"pots": pots, "needed_supply_over_pool_for_one_decade": 9.0, "best_available_at_phi_1": best,
            "reachable": bool(best >= 9.0),
            "reading": ("In the ambient pots the oxidant is NOT a budget: the consumers use well under one per cent of it, so it "
                        "acts as a constant multiplier on the dimer rate. The whole mercaptoketone flux at phi = 1 would raise "
                        "that multiplier by under one per cent. What those rows need is the dimer RATE CONSTANT, which sits at "
                        "its band ceiling and is opposed by Kumazawa's retention rows. No oxidant SOURCE reaches T3.")}


def g2_whitfield_share(B27, B23) -> Dict[str, Any]:
    x0 = B27.incumbent_vector()
    s = B23.SYSTEMS["whitfield_nf_cys"]
    curve = []
    first_reached = None
    for phi in PHI_GRID:
        xx = x0.copy(); xx[B27.K_SLOT] = math.log10(phi)
        p = B27.build_parameters(xx); _f, _e, _d, drift = B23.unpack(xx)
        run = integrate_sulfur(p, s["t_c"] + CELSIUS, s["initial"], np.array([0.0, s["minutes"]]), ph=s["ph"],
                               buffer_spec=s["buffer"], ph_drift=drift, rtol=1e-6, atol=1e-14)
        mft, d = run.final("MFT"), run.final("MFTD")
        share = 2.0 * d / (mft + 2.0 * d) if (mft + 2.0 * d) > 0 else 0.0
        curve.append({"phi": phi, "mft_disulfide_share": share, "total_mft_molpct": 100.0 * (mft + 2.0 * d) / float(B27.NF_MMOL_L)})
        if first_reached is None and share >= 0.35:
            first_reached = phi
    return {"charge": dict(s["initial"]), "curve": curve, "floor": 0.35, "phi_first_reaching_floor": first_reached,
            "at_ceiling": bool(first_reached is not None and first_reached >= 1.0),
            "reading": ("The pot whose oxidant is genuinely zero. The share climbs with phi and reaches the 35 % floor only at "
                        "phi = 1.0 exactly -- every mercaptoketone-forming event oxidising a thiol. The pre-registration (sec. 10) "
                        "declared that a phi pinned at its ceiling means the objective is asking the mercaptoketone flux for more "
                        "than it can supply, and is evidence against the structure rather than a fitted value.")}


def g3_fit_free_t2(B27, B23) -> Dict[str, Any]:
    x0 = B27.incumbent_vector()
    ids = [r["id"] for r in B23.ACTIVE_FIT_ROWS]
    sigma = {r["id"]: float(r["sigma_log"]) for r in B23.ACTIVE_FIT_ROWS if r["kind"] != "ph_endpoint"}
    lo, hi = x0.copy(), x0.copy()
    lo[B27.K_SLOT], hi[B27.K_SLOT] = -4.0, 0.0
    r_lo, r_hi = B27.residual_vector(lo, True), B27.residual_vector(hi, True)
    growth = []
    for i, rid in enumerate(ids):
        if rid not in sigma:
            continue
        growth.append({"row": rid, "residual_phi_0": float(r_lo[i]), "residual_phi_1": float(r_hi[i]),
                       "growth_dex": (abs(float(r_hi[i])) - abs(float(r_lo[i]))) * sigma[rid]})
    worst = max(growth, key=lambda g: g["growth_dex"])
    kum = [g for g in growth if g["row"] in KUMAZAWA]
    return {"cost_phi_0": 0.5 * float(np.dot(r_lo, r_lo)), "cost_phi_1": 0.5 * float(np.dot(r_hi, r_hi)),
            "rows_compared": len(growth), "worst": worst, "n_over_0_3_dex": sum(1 for g in growth if g["growth_dex"] > 0.3),
            "kumazawa_max_abs_growth_dex": max(abs(g["growth_dex"]) for g in kum), "kumazawa": kum,
            "improved": sorted((g for g in growth if g["growth_dex"] < -0.1), key=lambda g: g["growth_dex"]),
            "pass": worst["growth_dex"] <= 0.3}


def render(p: Dict[str, Any]) -> str:
    g1, g2, g3 = p["G1"], p["G2"], p["G3"]
    L = [f"# Wave B27 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         "**The fit was not run.** The decisive test cannot be reached by this structure or by any oxidant source, and the one pot "
         "the structure does fix needs its single coordinate at the physical ceiling the pre-registration declared disqualifying.", "",
         "| gate | result | pass |", "|---|---|---|",
         f"| G1 T3 reachability | " + "; ".join(f"{k}: consumers use {100 * v['consumed_fraction']:.3f} % of the pool, whole mercaptoketone flux at phi = 1 adds {100 * v['supply_at_phi_1_over_pool']:.2f} %" for k, v in g1["pots"].items())
         + f"; a decade needs about {g1['needed_supply_over_pool_for_one_decade']:.0f}x the pool | {g1['reachable']} |",
         f"| G2 Whitfield share | " + ", ".join(f"phi {c['phi']:g}: {100 * c['mft_disulfide_share']:.1f} %" for c in g2["curve"]) + f"; the 35 % floor is first reached at phi = {g2['phi_first_reaching_floor']} | {not g2['at_ceiling']} |",
         f"| G3 fit-free T2 | cost {g3['cost_phi_0']:.1f} (phi -> 0) -> {g3['cost_phi_1']:.1f} (phi = 1); worst growth {g3['worst']['row']} {g3['worst']['growth_dex']:+.2f} dex; {g3['n_over_0_3_dex']} rows over 0.3; Kumazawa max |growth| {g3['kumazawa_max_abs_growth_dex']:.4f} dex | {g3['pass']} |",
         f"| T1 reference pot at phi = 1 | ratios (dex) {p['T1']['ratios_dex']}; rising 6->12 h {p['T1']['still_rising_6_to_12_h']} | {p['T1']['pass']} |",
         f"| T5 Wang shape at phi = 1 | MFT decline {p['T5']['MFT']['decline_from_peak_dex']:.2f} dex, FFT {p['T5']['FFT']['decline_from_peak_dex']:.2f} dex | {p['T5']['pass']} |",
         f"| the same two under B9 itself (phi absent) | T1 ratios {p['b9_reference']['T1']['ratios_dex']}, rising {p['b9_reference']['T1']['still_rising_6_to_12_h']}, pass {p['b9_reference']['T1']['pass']}; T5 MFT decline {p['b9_reference']['T5']['MFT']['decline_from_peak_dex']:.2f}, FFT {p['b9_reference']['T5']['FFT']['decline_from_peak_dex']:.2f}, pass {p['b9_reference']['T5']['pass']} | reference |",
         "", "## G1, read plainly", "", g1["reading"], "", "## G2, read plainly", "", g2["reading"], "",
         "## What phi = 1 does to the rest of the objective", "",
         "Rows that improve by more than 0.1 dex: " + (", ".join(f"{g['row']} ({g['growth_dex']:+.2f})" for g in g3["improved"]) or "none") + ".",
         f"The four Kumazawa rows move by at most {g3['kumazawa_max_abs_growth_dex']:.4f} dex, as section 9 predicted (they carry no norfuraneol).", "",
         "## What is kept", "",
         "The step `ch_redox_mp3p`, the parameter `k_redox_mp3p` and the engine hook stay in the code, INERT at phi = 0, the way B17's and "
         "B25's refused structures do. The three Whitfield charge corrections are installed only by this wave's generator and are restored "
         "on exit, so the shipped objective is unchanged -- and that is a named debt: the next sulfur refit must carry the printed charges "
         "(norfuraneol 50, cysteine 50, H2S 97 mmol/L; 0.5 M phosphate pH 4.5; total MFT 0.230 mol %)."]
    return "\n".join(L) + "\n"


def main() -> int:
    if not B9_REPORT.exists():
        raise SystemExit(f"{B9_REPORT} missing")
    import generate_kinetic_core_b27_fit as B27  # noqa: E402  (configure() at import: the corrections are installed)
    import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
    try:
        x9 = B27.incumbent_vector()
        G1 = g1_t3_reachability(B27, B23, x9)
        G2 = g2_whitfield_share(B27, B23)
        G3 = g3_fit_free_t2(B27, B23)
        hof = B23.BUFFER_HOFMANN
    finally:
        B27.restore()
    rep = json.loads(B9_REPORT.read_text(encoding="utf-8"))
    op1, d1 = _operative_and_drift(rep, 0.0)
    T1, T5 = t1_reference_pot(op1, d1, hof), t5_wang_shape(op1, d1)
    op9, d9 = _operative_and_drift(rep, None)
    ref = {"T1": t1_reference_pot(op9, d9, hof), "T5": t5_wang_shape(op9, d9)}
    payload = {
        "artifact": "kinetic_core_b27_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b27_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b27_ship_rule.py",
                                                  inputs=[B9_REPORT, PREREG]),
        "prereg": data_paths.rel(PREREG),
        "rule": "SHIP if T1, T2, T3 (Zhou 2023's three dimer shares within 0.3 dex) and T6 hold; evaluated here as a GATE before the fit, at the shipped B9 vector with phi swept over its band",
        "fitted": False,
        "G1": G1, "G2": G2, "G3": G3, "T1": T1, "T5": T5, "b9_reference": ref,
        "verdict": "NOT FITTED -- DO NOT SHIP",
        "why": ["T3 is unreachable: the ambient pots are not oxidant-limited (G1)",
                "the one pot the structure fixes needs phi at its physical ceiling, which the pre-registration declared disqualifying (G2)"],
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| G1 reachable", G1["reachable"], f"best {G1['best_available_at_phi_1']:.4f} of the 9x needed",
          "| G2 floor at phi", G2["phi_first_reaching_floor"], "| G3", G3["pass"], f"worst {G3['worst']['growth_dex']:+.2f}",
          f"kumazawa {G3['kumazawa_max_abs_growth_dex']:.4f}", "| T1", T1["pass"], "(B9:", ref["T1"]["pass"], ") | T5", T5["pass"], "(B9:", ref["T5"]["pass"], ")")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
