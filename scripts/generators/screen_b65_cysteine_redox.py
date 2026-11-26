#!/usr/bin/env python
"""
screen_b65_cysteine_redox.py -- B65's plug-in mechanism-signature screen (lesson 2026-11-05), 2026-11-25.

Read-only diagnostic: the SHIPPED sulfur lane (``engine.core_parameters``), a handful of runs at shipped values and at
band corners, and three post-hoc / per-pot plug-ins where the engine lacks the structure. Nothing is fitted or sampled,
no engine file, constant or data file is touched. Every pot runs with the pH CLAMPED at its stated value (no drift,
no buffer model), isothermal.

    ./scripts/docker_maillard.sh run "python -m scripts.generators.screen_b65_cysteine_redox"

Writes results/validation/kinetic_core_b65_plugin_screen.json; the record is kinetic_core_b65_plugin_screen.md.
"""
from __future__ import annotations

import dataclasses
import json
import math
import os
from multiprocessing import Pool
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import brentq

from src import data_paths
from src.kinetic_core import engine
from src.kinetic_core.parameters_sulfur import OX_RESERVOIR_DEFAULT_UNITS, oxygen_parameters
from src.kinetic_core.sulfur import integrate_sulfur, sulfur_flux_budget

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b65_plugin_screen.json"
CELSIUS = 273.15
MW = {"MFT": 114.17, "MFTD": 226.32, "THI_HCl": 337.27, "XYL": 150.13, "CYS": 121.16}
WORKERS = min(8, int(os.environ.get("MAILLARD_WORKERS", "8")))

# ----------------------------------------------------------------------------------------------------------------
# Pots (mmol/L). The engine's Zhang Fig. 1 convention is generate_kinetic_core_b2_3_fit.py:404-417 (ambient OX 1.0;
# cystine = Cys 124.8 + OX 62.4).
# ----------------------------------------------------------------------------------------------------------------
OX_AMB = 1.0
ZHANG_T, ZHANG_MIN, ZHANG_PH = 115.0, 60.0, 4.9
GROSCH_T, GROSCH_MIN, GROSCH_PH = 100.0, 120.0, 5.7
FIG1 = {
    "z1_control": {"THI": 44.5, "PENT": 99.9, "OX": OX_AMB},
    "z1_cys": {"THI": 44.5, "PENT": 99.9, "Cys": 123.8, "OX": OX_AMB},
    "z1_gsh": {"THI": 44.5, "PENT": 99.9, "GSH": 48.8, "OX": OX_AMB},
    "z1_cystine": {"THI": 44.5, "PENT": 99.9, "Cys": 124.8, "OX": 62.4},
}
#: Not a Zhang arm: GSH's thiol charged as Cys (48.8 mM), a proxy for "GSH would scavenge oxidant as Cys does",
#: because the engine has no GSH + OX step. Reported beside S2 only.
GSH_PROXY = {"z1_gsh_as_cys": {"THI": 44.5, "PENT": 99.9, "Cys": 48.8, "OX": OX_AMB}}
FIG2_DOSES_MG_ML = (0.5, 2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 62.5, 75.0)
FIG2_THI_MM = (10.0, 30.0, 100.0, 133.0)   # absolute thiamine not printed; 133 mM = Met 15 mg/mL at 1:3:3
GROSCH_30 = (0.0, 0.17, 0.34, 0.51, 0.65, 0.68, 1.36, 2.72, 5.44)   # Cys mM, thiamine 30 umol / 0.5 L = 60 uM
GROSCH_6 = (0.34, 0.68, 1.36, 2.72)                                   # thiamine 6 umol / 0.5 L = 12 uM

MEASURED = {
    "S1_ratio": {"z1_cys": 1.34 / 0.07, "z1_gsh": 1.05 / 0.07, "z1_cystine": 2.01 / 0.07},
    "S2_share": {"z1_control": 0.04 / 0.11, "z1_cys": 0.115 / 1.455, "z1_gsh": 0.07 / 1.12, "z1_cystine": 1.09 / 3.10},
    "S3_mft_cys": dict(zip(FIG2_DOSES_MG_ML, (0.4, 0.75, 1.35, 1.6, 2.4, 2.5, 4.0, 5.7, 11.0, 11.0))),
    "S3_dimer_cys": dict(zip(FIG2_DOSES_MG_ML, (0.02, 0.03, 0.04, 0.04, 0.05, 0.10, 0.55, 1.00, 0.68, 0.55))),
    "S4_mft_ug": dict(zip(GROSCH_30, (5.2, 14, 18, 18, 19, 14, 10, 4.8, 2.0))),
    "S5_mft_ug": dict(zip(GROSCH_6, (2.1, 2.8, 1.6, 0.5))),
}


def pots() -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for name, init in {**FIG1, **GSH_PROXY}.items():
        out[name] = dict(init=init, t_c=ZHANG_T, minutes=ZHANG_MIN, ph=ZHANG_PH)
    for thi in FIG2_THI_MM:
        xyl = thi * MW["THI_HCl"] / MW["XYL"]   # xylose = thiamine by weight (1:3:3 Met:VB1:Xyl; Met not in the lane)
        for d in FIG2_DOSES_MG_ML:
            out[f"z2_thi{thi:g}_cys{d:g}"] = dict(
                init={"THI": thi, "PENT": xyl, "Cys": d / MW["CYS"] * 1000.0, "OX": OX_AMB},
                t_c=ZHANG_T, minutes=ZHANG_MIN, ph=ZHANG_PH)
    for c in GROSCH_30:
        init = {"THI": 0.060, "PENT": 1.34, "OX": OX_AMB}
        if c > 0:
            init["Cys"] = c
        out[f"g30_cys{c:g}"] = dict(init=init, t_c=GROSCH_T, minutes=GROSCH_MIN, ph=GROSCH_PH)
    for c in GROSCH_6:
        out[f"g6_cys{c:g}"] = dict(init={"THI": 0.012, "PENT": 1.34, "Cys": c, "OX": OX_AMB},
                                   t_c=GROSCH_T, minutes=GROSCH_MIN, ph=GROSCH_PH)
    return out


def thiol(init: Dict[str, float]) -> float:
    """Charged total free thiol (mM): Cys + GSH. A plug-in: constant over the run (cysteine loss ignored)."""
    return float(init.get("Cys", 0.0)) + float(init.get("GSH", 0.0))


# ----------------------------------------------------------------------------------------------------------------
# Variants: each maps (shipped params, pot) -> (params, initial state)
# ----------------------------------------------------------------------------------------------------------------
SUGAR_KEYS = ("k_nf_mft", "k_ddp_mft", "k_ddp_mft_hs", "k_ha_mp_mft")


def _scale(p: Dict[str, Any], key: str, f: float) -> None:
    p[key] = dataclasses.replace(p[key], k_ref=float(p[key].k_ref) * float(f))


def _oxygen(p: Dict[str, Any], k_cys_ox: float, supply_on: bool) -> None:
    ox = oxygen_parameters(k_cys_ox=k_cys_ox, k_ox_supply=10.0 if supply_on else 0.0)
    p["k_cys_ox"], p["k_ox_supply"] = ox["k_cys_ox"], ox["k_ox_supply"]


def g_cf089(c: float, a: float, kp: float, kd: float) -> float:
    return (1.0 + a * c / (c + kp)) / (1.0 + (c / kd) ** 2)


def variants() -> Dict[str, Tuple[str, Callable]]:
    v: Dict[str, Tuple[str, Callable]] = {}

    def shipped(p, pot):
        return p, dict(pot["init"])
    v["V1_shipped"] = ("1 engine as shipped", shipped)

    for d in (-1.9, -0.55, 0.55, 1.9):
        def f(p, pot, d=d):
            for k in SUGAR_KEYS:
                _scale(p, k, 10.0 ** d)
            return p, dict(pot["init"])
        v[f"V2_sugar_{d:+.2f}dex"] = ("2 sugar route x10^d", f)

    for lk in (-5.0, -3.0, -1.0):
        for res in (False, True):
            def f(p, pot, lk=lk, res=res):
                _oxygen(p, 10.0 ** lk, supply_on=res)
                init = dict(pot["init"])
                if res:
                    init["OXR"] = OX_RESERVOIR_DEFAULT_UNITS
                return p, init
            v[f"V3_kcysox_{lk:g}_{'B11res' if res else 'closed'}"] = ("3 k_cys_ox corner", f)

    for ox0 in (0.01, 0.1, 1.0, 10.0):
        for s in (0.0, 0.1, 1.0):
            for lk in (None, -1.0):
                def f(p, pot, ox0=ox0, s=s, lk=lk):
                    _oxygen(p, 0.0 if lk is None else 10.0 ** lk, supply_on=s > 0)
                    init = dict(pot["init"])
                    init["OX"] = init.get("OX", OX_AMB) - OX_AMB + ox0   # cystine's 62.4 is kept on top
                    if init["OX"] <= 0:
                        init["OX"] = ox0
                    if s > 0:
                        init["OXR"] = s * pot["minutes"]   # supply budget = rate x hold, held at OX0 while it lasts
                    return p, init
                tag = "" if lk is None else f"_kcysox{lk:g}"
                v[f"V4_ox{ox0:g}_sup{s:g}{tag}"] = ("4 oxidant charge / supply", f)

    corners = [(a, kp, kd) for a in (1.0, 10.0) for kp in (0.03, 3.0) for kd in (0.3, 30.0)]
    corners += [(10.0, 3.0, 1000.0), (10.0, 0.03, 1000.0)]   # the K_d a B1T would need at 115 C (premortem 2b)
    for a, kp, kd in corners:
        def f(p, pot, a=a, kp=kp, kd=kd):
            _scale(p, "k_hmp_mft", g_cf089(thiol(pot["init"]), a, kp, kd))
            return p, dict(pot["init"])
        v[f"V5b_g_a{a:g}_Kp{kp:g}_Kd{kd:g}"] = ("5b CF-089 g(c) on k_hmp_mft", f)

    for kc in (1.0, 10.0, 100.0):
        def f(p, pot, kc=kc):
            fac = 1.0 / (1.0 + thiol(pot["init"]) / kc)
            _scale(p, "k_mft_decay", fac)
            _scale(p, "k_thiolate_loss", fac)
            return p, dict(pot["init"])
        v[f"V5c_slowloss_Kc{kc:g}"] = ("5c cysteine-slowed MFT loss", f)

    # 6: combinations of the corners that each moved one signature, and the sugar route alone
    def combo(p, pot, sugar=0.0, kc=None, ox=None, thiamine_off=False):
        if sugar:
            for k in SUGAR_KEYS:
                _scale(p, k, 10.0 ** sugar)
        if thiamine_off:
            _scale(p, "k_hmp_mft", 0.0)
        if kc is not None:
            fac = 1.0 / (1.0 + thiol(pot["init"]) / kc)
            _scale(p, "k_mft_decay", fac)
            _scale(p, "k_thiolate_loss", fac)
        init = dict(pot["init"])
        if ox is not None:
            ox0, s, lk = ox
            _oxygen(p, 10.0 ** lk, supply_on=s > 0)
            init["OX"] = init.get("OX", OX_AMB) - OX_AMB + ox0
            if s > 0:
                init["OXR"] = s * pot["minutes"]
        return p, init
    six = {
        "V6_sugar_only_shipped": dict(thiamine_off=True),
        "V6_sugar_only_+1.90dex": dict(thiamine_off=True, sugar=1.9),
        "V6_ox10_sup1_kcysox-3": dict(ox=(10.0, 1.0, -3.0)),
        "V6_ox10_sup1_kcysox-2": dict(ox=(10.0, 1.0, -2.0)),
        "V6_ox10_sup0_kcysox-2": dict(ox=(10.0, 0.0, -2.0)),
        "V6_ox10_sup0_kcysox-2.5": dict(ox=(10.0, 0.0, -2.5)),
        "V6_ox10_sup0.5_kcysox-2.5": dict(ox=(10.0, 0.5, -2.5)),
        "V6_ox30_sup0_kcysox-2.5": dict(ox=(30.0, 0.0, -2.5)),
        "V6_sugar+1.9_Kc10": dict(sugar=1.9, kc=10.0),
        "V6_sugar+1.9_Kc1": dict(sugar=1.9, kc=1.0),
        "V6_sugar+0.55_Kc1": dict(sugar=0.55, kc=1.0),
        "V6_sugar+1.9_Kc10_ox10_sup1_kcysox-1": dict(sugar=1.9, kc=10.0, ox=(10.0, 1.0, -1.0)),
        "V6_sugar+1.9_Kc1_ox10_sup1_kcysox-1": dict(sugar=1.9, kc=1.0, ox=(10.0, 1.0, -1.0)),
    }
    # 7: a small grid around the scavenging edge, to test whether S2's four arms (GSH as the Cys proxy) can pass jointly
    for ox0 in (5.0, 10.0, 20.0):
        for sup in (0.0, 0.1, 0.25):
            for lk in (-2.75, -2.5, -2.25):
                six[f"V7_ox{ox0:g}_sup{sup:g}_kcysox{lk:g}"] = dict(ox=(ox0, sup, lk))
    for name, kw in six.items():
        v[name] = ("6 combination", lambda p, pot, kw=kw: combo(p, pot, **kw))
    return v


# ----------------------------------------------------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------------------------------------------------
_SHIPPED: Optional[Dict[str, Any]] = None
_POTS = pots()
_VARIANTS = variants()


def _run_one(job: Tuple[str, str]) -> Tuple[str, str, Dict[str, float]]:
    global _SHIPPED
    if _SHIPPED is None:
        _SHIPPED = engine.core_parameters(engine.SULFUR)
    vname, pname = job
    pot = _POTS[pname]
    params, init = _VARIANTS[vname][1](dict(_SHIPPED), pot)
    run = integrate_sulfur(params, pot["t_c"] + CELSIUS, init, np.array([0.0, pot["minutes"]]), ph=pot["ph"],
                           rtol=1e-8, atol=1e-16)
    res = {s: float(run.final(s)) for s in ("MFT", "MFTD", "Cys", "GSH", "OX", "OXR", "HMP", "H2S")}
    if vname == "V1_shipped" and pname in FIG1:
        b = sulfur_flux_budget(params, pot["t_c"] + CELSIUS, init, pot["minutes"], ph=pot["ph"])
        res["flux_thiamine_mft"] = float(b.get("r_hmp_mft", 0.0))
        res["flux_sugar_mft"] = float(sum(b.get(k, 0.0) for k in ("r_nf_mft", "r_ddp_mft", "r_ddp_mft_hs",
                                                                  "r_ha_mp_mft")))
    return vname, pname, res


def share(r: Dict[str, float]) -> float:
    m, d = r["MFT"] * MW["MFT"], r["MFTD"] * MW["MFTD"]
    return d / (m + d) if (m + d) > 0 else float("nan")


def exchange_share(r: Dict[str, float], partner_ssr: float, partner_sh: float, k_eq: float) -> float:
    """MFTD + 2 RSH <-> 2 MFT + RSSR at equilibrium, K = [MFT]^2[RSSR] / ([MFTD][RSH]^2), from the end state."""
    m, d = r["MFT"], r["MFTD"]
    if partner_sh <= 0 or d <= 0:
        return share(r)

    def f(x):
        return (2 * math.log(m + 2 * x) + math.log(partner_ssr + x)
                - math.log(k_eq) - math.log(d - x) - 2 * math.log(partner_sh - 2 * x))
    lo = -min(m / 2, partner_ssr) * (1 - 1e-12) if partner_ssr > 0 else 1e-300
    if partner_ssr <= 0:
        lo = 1e-30 * d
    hi = d * (1 - 1e-12)
    try:
        x = brentq(f, lo, hi, xtol=1e-30, rtol=1e-14, maxiter=500)
    except ValueError:
        x = hi if f(hi) < 0 else lo
    return share({"MFT": m + 2 * x, "MFTD": d - x})


def judge(res: Dict[str, Dict[str, Dict[str, float]]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for vname, rows in res.items():
        j: Dict[str, Any] = {}
        ctrl = max(rows["z1_control"]["MFT"], 1e-30)   # the sugar-only variants have no control MFT
        ratios = {a: rows[a]["MFT"] / ctrl for a in ("z1_cys", "z1_gsh", "z1_cystine")}
        j["S1_ratio"] = ratios
        j["S1_pass_0p3dex"] = all(abs(math.log10(ratios[a] / MEASURED["S1_ratio"][a])) <= 0.3 for a in ratios)
        j["S1_pass_floor5"] = all(ratios[a] >= 5.0 for a in ratios)
        shares = {a: share(rows[a]) for a in FIG1}
        j["S2_share"] = shares
        j["S2_gsh_proxy_share"] = share(rows["z1_gsh_as_cys"])
        j["S2_pass"] = (shares["z1_cys"] < 0.15 and shares["z1_gsh"] < 0.15 and shares["z1_cystine"] > 0.2
                        and shares["z1_control"] > 0.2)
        s3 = {}
        for thi in FIG2_THI_MM:
            mft = [max(rows[f"z2_thi{thi:g}_cys{d:g}"]["MFT"], 1e-30) for d in FIG2_DOSES_MG_ML]
            dim = [share(rows[f"z2_thi{thi:g}_cys{d:g}"]) for d in FIG2_DOSES_MG_ML]
            upto = FIG2_DOSES_MG_ML.index(62.5)
            steps = [math.log10(mft[i + 1] / mft[i]) for i in range(upto)]
            rise = mft[upto] / mft[0]
            s3[f"{thi:g}"] = {"rise_62p5_over_0p5": rise, "worst_step_dex": min(steps),
                              "tail_75_over_62p5": mft[-1] / mft[upto],
                              "pass": rise >= 9.0 and min(steps) >= -0.1,
                              "dimer_share_0p5_to_62p5": (dim[0], dim[upto]), "mft_mM": mft}
        j["S3"] = s3
        j["S3_pass_any_thi"] = any(x["pass"] for x in s3.values())
        j["S3_pass_all_thi"] = all(x["pass"] for x in s3.values())
        g30 = [rows[f"g30_cys{c:g}"]["MFT"] for c in GROSCH_30]
        g6 = [rows[f"g6_cys{c:g}"]["MFT"] for c in GROSCH_6]
        opt30 = GROSCH_30[int(np.argmax(g30))]
        opt6 = GROSCH_6[int(np.argmax(g6))]
        g0 = max(g30[0], 1e-30)
        j["S4"] = {"optimum_mM": opt30, "max_over_cysfree": max(g30) / g0, "fall_5p44_over_0": g30[-1] / g0,
                   "pass_optimum": 0.17 <= opt30 <= 1.36, "pass_fall": g30[-1] < g30[0], "mft_mM": g30}
        j["S5"] = {"optimum_mM": opt6, "pass": opt6 in (0.34, 0.68, 1.36) and 0 < int(np.argmax(g6)) < 3,
                   "mft_mM": g6}
        out[vname] = j
    return out


def post_hoc_exchange(shipped_rows: Dict[str, Dict[str, float]], tag: str) -> Dict[str, Any]:
    """Variant 5a on one variant's Fig. 1 end states. Partner disulfide: none in the Cys / GSH arms (engine has no
    cystine; k_cys_ox makes fragments); in the cystine arm the unconsumed oxidant (= the cystine charge, engine
    convention) or a solubility-limited 2 mM."""
    out: Dict[str, Any] = {"variant": tag}
    for k_eq in (0.1, 1.0, 10.0):
        row = {"z1_control": share(shipped_rows["z1_control"])}
        r = shipped_rows["z1_cys"]
        row["z1_cys"] = exchange_share(r, 0.0, r["Cys"], k_eq)
        r = shipped_rows["z1_gsh"]
        row["z1_gsh"] = exchange_share(r, 0.0, r["GSH"], k_eq)
        r = shipped_rows["z1_cystine"]
        row["z1_cystine_ssr_is_OX"] = exchange_share(r, r["OX"], r["Cys"], k_eq)
        row["z1_cystine_ssr_2mM"] = exchange_share(r, 2.0, r["Cys"], k_eq)
        # the cystine arm with NO free cysteine charged (cystine as disulfide only): what free thiol would S2 need?
        m = r["MFT"]
        need = math.sqrt(m * r["OX"] / (k_eq * (0.35 / 0.65 * MW["MFT"] / MW["MFTD"])))
        row["cystine_arm_free_RSH_needed_for_share_0p35_mM"] = need
        out[f"K{k_eq:g}"] = row
    return out


def main() -> int:
    jobs = [(v, p) for v in _VARIANTS for p in _POTS]
    try:
        os.nice(10)
    except OSError:
        pass
    with Pool(WORKERS) as pool:
        flat = pool.map(_run_one, jobs, chunksize=8)
    res: Dict[str, Dict[str, Dict[str, float]]] = {}
    for v, p, r in flat:
        res.setdefault(v, {})[p] = r
    verdict = judge(res)
    best_ox = "V4_ox10_sup1_kcysox-1"
    report = {
        "artifact": "kinetic_core_b65_plugin_screen",
        "engine_fit_report": data_paths.rel(engine._B2_FIT_REPORT),
        "measured": {k: {str(kk): vv for kk, vv in d.items()} for k, d in MEASURED.items()},
        "variant_family": {v: _VARIANTS[v][0] for v in _VARIANTS},
        "verdict": verdict,
        "post_hoc_exchange": [post_hoc_exchange(res["V1_shipped"], "V1_shipped"),
                              post_hoc_exchange(res[best_ox], best_ox)],
        "shipped_fig1_flux": {p: {k: res["V1_shipped"][p][k] for k in ("flux_thiamine_mft", "flux_sugar_mft")}
                              for p in FIG1},
        "raw": res,
    }
    OUT.write_text(json.dumps(report, indent=1, default=float))
    print(f"{'variant':34s} {'S1 cys/gsh/cyss':>22s} {'S1':>4s} {'S2 ctl/cys/gsh/cyss':>27s} {'S2':>3s} "
          f"{'S3 rise@thi10/30/100/133':>28s} {'S3':>3s} {'S4 opt,fall':>13s} {'S5':>4s}")
    for v, j in verdict.items():
        s1 = j["S1_ratio"]
        s2 = j["S2_share"]
        s3 = j["S3"]
        print(f"{v:34s} {s1['z1_cys']:6.2f}/{s1['z1_gsh']:6.2f}/{s1['z1_cystine']:6.2f} "
              f"{'P' if j['S1_pass_0p3dex'] else ('f' if j['S1_pass_floor5'] else 'F'):>4s} "
              f"{s2['z1_control']:.3f}/{s2['z1_cys']:.3f}/{s2['z1_gsh']:.3f}/{s2['z1_cystine']:.3f} "
              f"{'P' if j['S2_pass'] else 'F':>3s} "
              + "/".join(f"{s3[t]['rise_62p5_over_0p5']:.2f}" for t in s3)
              + f" {'P' if j['S3_pass_all_thi'] else ('p' if j['S3_pass_any_thi'] else 'F'):>3s} "
              f"{j['S4']['optimum_mM']:g},{j['S4']['fall_5p44_over_0']:.2f} {'P' if j['S5']['pass'] else 'F':>4s}")
    print(json.dumps(report["post_hoc_exchange"], indent=1))
    print(json.dumps(report["shipped_fig1_flux"], indent=1))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
