"""CF-076: the cysteine + xylose <-> TTCA exchange across dilute (Zhai 2021) and concentrated (Zhai 2019, Cao 2017)
pots. Rivals E0-E4 and negative controls NC1 / NC2 (prereg.yaml, written before this file was run), scored by
feasibility sampling inside declared bands -- no fitting. Standalone reduced module in the A2_MODULE_CONTRACT sec. 3
style: CONFIG_KEYS (bands + basis), rate_constants(params, ctx), OBSERVABLES; forced inputs are piecewise-linear.

Reuses the thiol-rivals cys_fate.py forced-TTCA skeleton for the dilute pots (docs/history/thiol_rivals_2026-10-31/).

Run (from the repo root):
  ./scripts/docker_maillard.sh run "OMP_NUM_THREADS=1 PYTHONPATH=$PWD nice -n 10 python <scratch>/cys_exchange.py [n_score] [n_truth]"
"""
from __future__ import annotations

import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

from src import data_paths

OUT = Path(__file__).resolve().parent
ROWS_ROUTE = data_paths.LIT_DIR / "fit_rows" / "route_submodel.yml"
ROWS_PH = data_paths.LIT_DIR / "fit_rows" / "ph_acid_submodel.yml"
R_KJ = 8.314e-3
PH = (5.5, 6.0, 7.0, 8.0)
T_GRID = np.array([0.0, 20, 40, 60, 80, 100, 120])
RIVALS = ("E0", "E1", "E2", "E3", "E4", "NC1", "NC2")
PRE = ("E0", "E1", "E2", "E3", "E4")


def arr(t_c: float, ea: float, t_ref_c: float = 100.0) -> float:
    return math.exp(-ea / R_KJ * (1.0 / (t_c + 273.15) - 1.0 / (t_ref_c + 273.15)))


# ---- forced inputs: Zhai 2021 Fig. 4b / Fig. 3a TTCA (mmol/L; "0.00" plotted = < LOD -> 0.005) -------------------
TTCA21 = {  # (pH, xyl, T_C)
    (5.5, 0, 100): [10, 9.81, 9.34, 9.10, 8.19, 7.68, 6.69], (6.0, 0, 100): [10, 9.78, 9.18, 8.79, 7.97, 5.98, 4.37],
    (7.0, 0, 100): [10, 9.14, 8.21, 7.02, 5.52, 2.08, 0.89], (8.0, 0, 100): [10, 7.79, 5.18, 3.79, 1.47, 0.005, 0.005],
    (5.5, 1, 100): [10, 9.79, 9.19, 8.79, 7.99, 7.34, 6.59], (6.0, 1, 100): [10, 9.68, 8.97, 8.36, 7.29, 5.12, 3.89],
    (7.0, 1, 100): [10, 8.57, 7.93, 6.97, 4.29, 0.005, 0.005], (8.0, 1, 100): [10, 8.12, 6.07, 2.16, 0.005, 0.005, 0.005],
    (5.5, 0, 120): [10, 8.77, 7.27, 6.18, 4.59, 3.17, 2.47], (5.5, 0, 140): [10, 7.59, 5.87, 3.98, 2.29, 1.20, 0.27]}
POTS21 = list(TTCA21)
# Zhai 2021 Fig. 5 dicarbonyl sum (3-DX + 1-DX + MGO; GO < 5 uM omitted), umol/L -> mM; E3's adduct sink only.
_DX3 = {(5.5, 0): [0, 280, 370, 520, 650, 640, 620], (6.0, 0): [0, 330, 480, 590, 750, 890, 930],
        (7.0, 0): [0, 150, 210, 380, 470, 510, 420], (8.0, 0): [0, 62, 72, 86, 98, 107, 75],
        (5.5, 1): [0, 290, 400, 580, 700, 750, 770], (6.0, 1): [0, 380, 510, 630, 810, 1030, 1090],
        (7.0, 1): [0, 170, 230, 360, 570, 770, 930], (8.0, 1): [0, 48, 76, 107, 147, 258, 359]}
_DX1 = {(5.5, 0): [0, 10, 13, 25, 37, 36, 32], (6.0, 0): [0, 17.2, 21.5, 33.7, 49.8, 62.8, 61.9],
        (7.0, 0): [0, 21, 29, 48, 61, 79, 69], (8.0, 0): [0, 32, 52, 74, 99, 117, 106],
        (5.5, 1): [0, 12, 13, 28, 40, 43, 47], (6.0, 1): [0, 17.9, 23, 34, 51, 68.8, 79],
        (7.0, 1): [0, 26, 31, 53, 75, 94, 125], (8.0, 1): [0, 31, 53, 82, 119, 172, 187]}
_MGO = {(5.5, 0): [0, 8.1, 10.7, 13.2, 13.9, 13.7, 11.4], (6.0, 0): [0, 12.9, 17.8, 24.7, 25.2, 24.8, 20.1],
        (7.0, 0): [0, 21.9, 27.9, 34.8, 31.1, 29.9, 26.9], (8.0, 0): [0, 30.2, 45.2, 57.4, 49.9, 34.7, 27.9],
        (5.5, 1): [0, 10.2, 12.7, 14.4, 19.9, 21.3, 24.7], (6.0, 1): [0, 16.9, 19.9, 21.4, 24.9, 27.8, 30.0],
        (7.0, 1): [0, 19.8, 23.1, 27.9, 31.8, 37.4, 42.5], (8.0, 1): [0, 37.4, 55.3, 63.8, 67.2, 59.7, 56.9]}
DC = {k: (np.array(_DX3[k]) + np.array(_DX1[k]) + np.array(_MGO[k])) / 1000.0 for k in _DX3}
DRIFT_KANG = {5.5: 0.4, 6.0: 1.0, 7.0: 2.1, 8.0: 3.3}   # measured pH drop of 10 mM TTCA, 120 min (Kang 2026 p. 3242)

# Zhai 2019 Fig. 3a (dossier digits; t = 0 TTCA 14.5 = the RT-formed initial condition), CAMA forced
Z19_T = np.array([0.0, 20, 40, 60, 80, 100, 120, 140, 160])
Z19_CAMA = {70: [0, 0.2, 0.5, 1.0, 1.3, 1.6, 1.9, 2.1, 2.8], 90: [0, 0.9, 1.4, 1.7, 1.9, 2.6, 3.0, 3.4, 3.7]}
Z19_CYS0, Z19_XYL0, Z19_TTCA0 = 82.6, 165.2, 14.5


def _rows(path):
    raw = yaml.safe_load(Path(path).read_text())
    rows = raw.get("rows") if isinstance(raw, dict) and "rows" in raw else None
    if rows is None:
        rows = [r for v in raw.values() if isinstance(v, list) for r in v] if isinstance(raw, dict) else raw
    return rows


def load_data():
    route = _rows(ROWS_ROUTE)
    phr = _rows(ROWS_PH)
    d = {"cys21": {}, "xyl21": {}, "z19": {70: {}, 90: {}}, "cao": {"TTCA": {}, "CAMA": {}, "Cys": {}, "PENT": {}},
         "cao_ph": {}}
    tagph = {"ph55": 5.5, "ph6": 6.0, "ph7": 7.0, "ph8": 8.0}
    for r in route:
        pot, sp = r["pot"], r["species"]
        if pot.startswith("zhai2021_") and pot.endswith("_fig6"):
            tag, ph = pot.split("_")[1], tagph[pot.split("_")[2]]
            if sp == "Cys":
                d["cys21"].setdefault((ph, 1 if tag == "ttcaxyl" else 0), {})[int(r["t_min"])] = float(r["value"])
            elif sp == "PENT":
                d["xyl21"].setdefault((ph, tag), {})[int(r["t_min"])] = float(r["value"])
        elif pot in ("zhai2019_70c", "zhai2019_90c") and sp == "TTCA":
            d["z19"][int(pot[9:11])][int(r["t_min"])] = float(r["value"])
        elif pot in ("cao2017_blank_fig4", "cao2017_blank_fig5") and sp in d["cao"]:
            d["cao"][sp][int(r["t_min"])] = float(r["value"])
    for r in phr:
        if r["pot"] == "cao2017_blank" and r["species"] == "pH":
            d["cao_ph"][int(r["t_min"])] = float(r["value"])
    return d


DATA = load_data()
K_XYL = {ph: -math.log(DATA["xyl21"][(ph, "xyl")][120] / 10.0) / 120.0 for ph in PH}
CAO_T = np.array(sorted(DATA["cao"]["TTCA"]))
CAO_CAMA_T = np.array([0.0] + list(CAO_T))
CAO_CAMA = np.array([0.0] + [DATA["cao"]["CAMA"][int(t)] for t in CAO_T])
CAO_PH_T = np.array(sorted(DATA["cao_ph"]))
CAO_PH = np.array([DATA["cao_ph"][int(t)] for t in CAO_PH_T])

# ---- CONFIG_KEYS: name -> (lo, hi, basis). Log-uniform where hi / lo > 3 ---------------------------------------
COMMON = {
    "k_th_100": (7.4e-4, 2.9e-3, "Kang 2026 Fig. S4 free Cys 10 mM pH 7, 1.47e-3 /min at 100 C x[0.5, 2]; Ea 55.1"),
    "k_cd_55": (0.0035, 0.014, "CAMA decay pH 5.5 ~100 C, de Roos 2005 x[0.5, 2]; Ea 85.7"),
    "cd_ph_slope": (0.31, 0.62, "log10 / pH above 5.5 (de Roos 5.5 -> 6.0 slope; half as lower edge)"),
    "release": (0.5, 1.0, "share of CAMA decay that frees Cys"),
    "w_cama": (0.5, 1.0, "share of TTCA's non-exchange loss through CAMA (engine 0.92)"),
    "K_ref": (0.045, 0.075, "Cao 2017 K_app 58.5-64.6 L/mol at 90 C, pH 4.2-4.7 (L/mmol)"),
}
CATAL = {"c_aa": (0.3, 30.0, "Fife 1991 general-base k_B/k0 4-83 /M, Bronsted 0.29 transfer"),
         "c_P": (4.0, 300.0, "Fife 1991 phosphate on I 83 /M of HPO4 2-; II ~20 by Bronsted")}
PHLAW = {"k_r45_100": (1e-5, 1.0, "declared, as E1 k0"), "a": (0.0, 1.0, "Fife 1991 H3O+ limb vs flat"),
         "pKa_T": (5.6, 6.6, "Kallen TC ring-N 6.24 +-0.4"), "dpKa_CT": (1.6, 2.6, "Kallen 8.38 - 6.24 = 2.14 +-0.5"),
         "Ea": (40.0, 100.0, "thiazolidine / Schiff-base class; Fife IV 58"),
         "s21": (0.5, 1.5, "x Kang 2026 measured drift"), "g21": (0.5, 2.0, "drift shape"),
         "d21": (0.0, 0.4, "+Xyl extra drift (xylose acids)"),
         "ph0_19": (6.0, 7.4, "Zhai 2019 pH after RT formation"), "phinf_19": (5.0, 6.5, "proton bookkeeping, TTCA pKa"),
         "tau_19": (5.0, 40.0, "formation timescale, min")}
CONFIG_KEYS = {
    "E0": {"k_r_100": (0.0143, 0.0572, "engine 0.0286 x[0.5, 2]")},
    "E1": {"k0_100": (1e-5, 1.0, "declared 5 decades (Kallen / Schubert fast end, Fife slow end)"),
           "Ea": (40.0, 100.0, "class"), **CATAL},
    "E2": dict(PHLAW),
    "E3": {"k_r_100": (0.0143, 0.0572, "engine"), "k_ox": (1e-4, 5e-3, "declared"), "k_s": (1e-4, 1e-1, "declared")},
    "E4": {**PHLAW, **CATAL},
    "NC1": {},
    "NC2": {},
}


def draw(rival: str, rng) -> dict:
    out = {}
    for name, (lo, hi, _) in {**COMMON, **CONFIG_KEYS[rival]}.items():
        out[name] = 10 ** rng.uniform(math.log10(lo), math.log10(hi)) if (lo > 0 and hi / lo > 3) else rng.uniform(lo, hi)
    return out


def from_unit(rival: str, u) -> dict:
    keys = {**COMMON, **CONFIG_KEYS[rival]}
    out = {}
    for (name, (lo, hi, _)), x in zip(keys.items(), u):
        out[name] = 10 ** (math.log10(lo) + x * (math.log10(hi) - math.log10(lo))) if (lo > 0 and hi / lo > 3) else lo + x * (hi - lo)
    return out


# ---- rate_constants(params, ctx): ctx = (rival, T_C, pH, AA_M, HPO4_M) ---------------------------------------------
def k_app(p, ph):
    """E2 / E4: Kallen species bookkeeping, anchored K(4.5) = K_ref."""
    pkc = p["pKa_T"] + p["dpKa_CT"]
    g = lambda x: (1 + 10 ** (x - p["pKa_T"])) / (1 + 10 ** (x - pkc))
    return p["K_ref"] * g(ph) / g(4.5)


def rate_constants(p: dict, rival: str, t_c: float, ph: float, aa_m: float, hpo4_m: float):
    """(k_f L/mmol/min, k_r /min) of the exchange."""
    if rival in ("E0", "E3"):
        kr = p["k_r_100"] * arr(t_c, 64.1)
        return p["K_ref"] * kr, kr
    if rival == "E1":
        kr = p["k0_100"] * arr(t_c, p["Ea"]) * (1 + p["c_aa"] * aa_m + p["c_P"] * hpo4_m)
        return p["K_ref"] * kr, kr
    if rival in ("E2", "E4"):
        kr = p["k_r45_100"] * arr(t_c, p["Ea"]) * 10 ** (-p["a"] * (ph - 4.5))
        if rival == "E4":
            kr *= 1 + p["c_aa"] * aa_m + p["c_P"] * hpo4_m
        return k_app(p, ph) * kr, kr
    if rival == "NC1":
        return 0.06 * 3.0, 3.0
    return 0.0, 0.0   # NC2


def insitu21(p, rival, ph_label, xyl, t):
    if rival not in ("E2", "E4"):
        return ph_label
    return ph_label - p["s21"] * DRIFT_KANG[ph_label] * (t / 120.0) ** p["g21"] - p["d21"] * xyl * t / 120.0


def insitu19(p, rival, t):
    if rival not in ("E2", "E4"):
        return 7.4
    return p["phinf_19"] + (p["ph0_19"] - p["phinf_19"]) * math.exp(-t / p["tau_19"])


def hpo4(ph, p_tot_m):
    return p_tot_m / (1 + 10 ** (7.1 - ph))


# ---- pot simulators -------------------------------------------------------------------------------------------
def sim21(p, rival, pot):
    """Dilute fed-TTCA pot, TTCA forced; returns Cys, Xyl on T_GRID and the gate violation (fraction of charge)."""
    ph, xyl, t_c = pot
    tt = np.array(TTCA21[pot], float)
    slopes = -np.diff(tt) / 20.0
    dc = DC.get((ph, xyl)) if t_c == 100 else None
    kth = p["k_th_100"] * arr(t_c, 55.1)
    kcd = p["k_cd_55"] * 10 ** (p["cd_ph_slope"] * (ph - 5.5)) * arr(t_c, 85.7)
    kx = K_XYL[ph] * arr(t_c, 100.0)
    ksink = p.get("k_ox", 0.0)
    ks = p.get("k_s", 0.0)

    def rhs(t, y):
        cys, x, cama = max(y[0], 0.0), max(y[1], 0.0), max(y[2], 0.0)
        T = np.interp(t, T_GRID, tt)
        L = slopes[min(int(t // 20.0), 5)]
        kf, kr = rate_constants(p, rival, t_c, insitu21(p, rival, ph, xyl, t), (T + cys + cama) / 1000.0, 0.0)
        r_open, r_close = kr * T, kf * cys * x
        other = r_close - r_open + L
        v = 0.0
        if other < 0.0:
            v, r_open, other = -other, r_open + other, 0.0
        r_cd = kcd * cama
        sink = (kth + ksink + (ks * np.interp(t, T_GRID, dc) if dc is not None else 0.0)) * cys
        return [r_open - r_close + p["release"] * r_cd - sink, r_open - r_close - kx * x, p["w_cama"] * other - r_cd, v]

    sol = solve_ivp(rhs, (0, 120), [0.0, 10.0 * xyl, 0.0, 0.0], t_eval=T_GRID, method="LSODA", rtol=1e-6,
                    atol=1e-9, max_step=2.0)
    return np.maximum(sol.y[0], 1e-6), np.maximum(sol.y[1], 1e-6), min(sol.y[3][-1] / 10.0, 1.0)


def sim19(p, rival, t_c):
    cama = np.array(Z19_CAMA[t_c], float)
    cslope = np.maximum(np.diff(cama) / 20.0, 0.0)
    kth = p["k_th_100"] * arr(t_c, 55.1) + p.get("k_ox", 0.0)

    def rhs(t, y):
        T, C, X = max(y[0], 0.0), max(y[1], 0.0), max(y[2], 0.0)
        A = np.interp(t, Z19_T, cama)
        kf, kr = rate_constants(p, rival, t_c, insitu19(p, rival, t), (T + C + A) / 1000.0, 0.0)
        r = kf * C * X - kr * T
        drain = cslope[min(int(t // 20.0), 7)]
        fT = T / (T + 5.0)
        return [r - fT * drain, -r - kth * C - (1 - fT) * drain, -r - (1 - fT) * drain]

    sol = solve_ivp(rhs, (0, 160), [Z19_TTCA0, Z19_CYS0 - Z19_TTCA0, Z19_XYL0 - Z19_TTCA0], t_eval=Z19_T,
                    method="LSODA", rtol=1e-6, atol=1e-8, max_step=2.0)
    return sol.y[0]


def sim_cao(p, rival):
    cslope = np.maximum(np.diff(CAO_CAMA) / np.diff(CAO_CAMA_T), 0.0)
    kth = p["k_th_100"] * arr(90.0, 55.1) + p.get("k_ox", 0.0)

    def rhs(t, y):
        T, C, X = max(y[0], 0.0), max(y[1], 0.0), max(y[2], 0.0)
        ph = float(np.interp(t, CAO_PH_T, CAO_PH))
        A = np.interp(t, CAO_CAMA_T, CAO_CAMA)
        kf, kr = rate_constants(p, rival, 90.0, ph, (T + C + A) / 1000.0, hpo4(ph, 0.2))
        r = kf * C * X - kr * T
        i = int(np.searchsorted(CAO_CAMA_T, t, side="right") - 1)
        drain = cslope[min(max(i, 0), len(cslope) - 1)]
        fT = T / (T + 10.0)
        return [r - fT * drain, -r - kth * C - (1 - fT) * drain, -r - (1 - fT) * drain]

    sol = solve_ivp(rhs, (0, 300), [0.0, 400.0, 400.0], t_eval=CAO_T, method="LSODA", rtol=1e-6, atol=1e-8,
                    max_step=2.0)
    return np.maximum(sol.y[0], 1e-3), np.maximum(sol.y[1], 1e-3), np.maximum(sol.y[2], 1e-3)


# ---- OBSERVABLES ----------------------------------------------------------------------------------------------
TESTS = ("T1", "T2", "T3", "T4", "T5", "T6", "T7")
T3_PH = (5.5, 6.0, 7.0)


def obs_from(cys21, xyl21_end, ctrl_end, gate, z19, cao):
    o = {}
    o["T1"] = np.array([np.mean([math.log10(cys21[(ph, 1)][t] / cys21[(ph, 0)][t]) for t in (60, 80, 100, 120)]) for ph in PH])
    o["T2"] = np.array([math.log10(cys21[(ph, 0)][120]) for ph in PH])
    o["T3"] = np.array([math.log10(xyl21_end[ph] / ctrl_end[ph]) for ph in T3_PH])
    o["T4"] = np.array(gate)
    f, pl = [], []
    for tc in (70, 90):
        s = z19[tc]
        den = s[60] - Z19_TTCA0
        f.append((s[20] - Z19_TTCA0) / den if den > 1.0 else 0.0)
        pl.append(math.log10(max(np.mean([s[t] for t in (60, 80, 100, 120, 140, 160)]), 1e-3)))
    o["T5"], o["T6"] = np.array(f), np.array(pl)
    T, C, X = cao
    q = {t: T[t] / (C[t] * X[t]) for t in T}
    o["T7"] = np.array([math.log10(q[30] / np.mean([q[t] for t in (60, 120, 180, 240, 300)]))])
    return o


def observed():
    d = DATA
    xt = {ph: d["xyl21"][(ph, "ttcaxyl")][120] for ph in PH}
    xc = {ph: d["xyl21"][(ph, "xyl")][120] for ph in PH}
    z19 = {tc: {0: Z19_TTCA0, **d["z19"][tc]} for tc in (70, 90)}
    cao = ({t: d["cao"]["TTCA"][t] for t in CAO_T}, {t: d["cao"]["Cys"][t] for t in CAO_T},
           {t: d["cao"]["PENT"][t] for t in CAO_T})
    return obs_from(d["cys21"], xt, xc, [0.0] * len(POTS21), z19, cao)


def predict(args):
    p, rival = args
    try:
        cys, xt, gate = {}, {}, []
        for pot in POTS21:
            c, x, v = sim21(p, rival, pot)
            gate.append(v)
            if pot[2] == 100:
                cys[(pot[0], pot[1])] = dict(zip(T_GRID.astype(int).tolist(), c.tolist()))
                if pot[1]:
                    xt[pot[0]] = x[-1]
        ctrl = {ph: 10.0 * math.exp(-K_XYL[ph] * 120) for ph in PH}
        z19 = {tc: dict(zip(Z19_T.astype(int).tolist(), sim19(p, rival, tc).tolist())) for tc in (70, 90)}
        T, C, X = sim_cao(p, rival)
        cao = tuple(dict(zip(CAO_T.astype(int).tolist(), a.tolist())) for a in (T, C, X))
        return obs_from(cys, xt, ctrl, gate, z19, cao)
    except Exception as exc:  # a failed integration is recorded, not hidden
        return {"error": repr(exc)}


# ---- likelihoods (prereg sigmas) -----------------------------------------------------------------------------
SIG = {"T1": 0.07, "T3": 0.05, "T4": 0.03, "T5": 0.06, "T6": 0.04, "T7": 0.06}
TAU2, PT2 = 0.20, math.hypot(0.035, 0.05)
_COV2 = np.full((4, 4), TAU2 ** 2) + np.eye(4) * PT2 ** 2
_CI2, _LD2 = np.linalg.inv(_COV2), math.log(np.linalg.det(_COV2))


def loglik(pred, y, test):
    d = pred - y
    if test == "T2":
        return -0.5 * np.einsum("...i,ij,...j->...", d, _CI2, d) - 0.5 * _LD2
    s = SIG[test]
    return -0.5 * np.sum((d / s) ** 2, axis=-1) - d.shape[-1] * math.log(s)


def logmean(ll):
    m = np.max(ll)
    return float(m + math.log(np.mean(np.exp(ll - m))))


def zmax(pred, y, test):
    """max |z| per draw (T2: after the best figure offset clipped to 2 tau)."""
    d = pred - y
    if test == "T2":
        off = np.clip(d.mean(axis=-1, keepdims=True), -2 * TAU2, 2 * TAU2)
        return np.max(np.abs(d - off), axis=-1) / PT2
    return np.max(np.abs(d), axis=-1) / SIG[test]


def noisy(pred, test, rng):
    if test == "T2":
        return pred + rng.normal(0, TAU2) + rng.normal(0, PT2, pred.shape)
    return pred + rng.normal(0, SIG[test], pred.shape)


def run(n_score, n_truth, workers=2):
    rng = np.random.default_rng(20261031)
    obs = observed()
    draws = {rv: [draw(rv, rng) for _ in range(n_score + n_truth)] for rv in RIVALS}
    t0 = time.time()
    preds = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for rv in RIVALS:
            res = list(ex.map(predict, [(p, rv) for p in draws[rv]], chunksize=8))
            bad = [i for i, r in enumerate(res) if "error" in r]
            for i in bad:   # replace a failed integration by a fresh draw of the same rival (recorded)
                for _ in range(5):
                    p = draw(rv, rng)
                    r = predict((p, rv))
                    if "error" not in r:
                        draws[rv][i], res[i] = p, r
                        break
            preds[rv] = {t: np.array([r[t] for r in res]) for t in TESTS}
            preds[rv]["_failed"] = len(bad)
            print(f"{rv}: {len(res)} draws, {len(bad)} failed integrations, {time.time() - t0:.0f} s", flush=True)
    return obs, draws, preds


def analyse(obs, draws, preds, n_score, n_truth):
    rng = np.random.default_rng(7)
    S = slice(0, n_score)
    out = {"observed": {t: obs[t].round(4).tolist() for t in TESTS}}
    # prediction summary
    out["pred"] = {rv: {t: {"p05": np.percentile(preds[rv][t][S], 5, axis=0).round(3).tolist(),
                            "median": np.median(preds[rv][t][S], axis=0).round(3).tolist(),
                            "p95": np.percentile(preds[rv][t][S], 95, axis=0).round(3).tolist()} for t in TESTS}
                   for rv in RIVALS}
    # power: pairwise, per test
    power = {}
    for t in TESTS:
        power[t] = {}
        for i, a in enumerate(PRE + ("NC1", "NC2")):
            for b in (PRE + ("NC1", "NC2"))[i + 1:]:
                ok = []
                for truth, other in ((a, b), (b, a)):
                    for j in range(n_truth):
                        y = noisy(preds[truth][t][n_score + j], t, rng)
                        ok.append(logmean(loglik(preds[truth][t][S], y, t)) > logmean(loglik(preds[other][t][S], y, t)))
                power[t][f"{a}|{b}"] = float(np.mean(ok))
    out["power_pairwise"] = power
    scored = [t for t in TESTS if any(v >= 0.8 for k, v in power[t].items() if "NC" not in k)]
    out["scored_tests"] = scored
    out["dropped_tests"] = [t for t in TESTS if t not in scored]
    # scores
    out["log_marginal"] = {t: {rv: logmean(loglik(preds[rv][t][S], obs[t], t)) for rv in RIVALS} for t in TESTS}
    out["feasibility"] = {t: {rv: float(np.mean(zmax(preds[rv][t][S], obs[t], t) <= 2.0)) for rv in RIVALS} for t in TESTS}

    def joint(rv, tests):
        ll = sum(loglik(preds[rv][t][S], obs[t], t) for t in tests)
        z = np.max(np.stack([zmax(preds[rv][t][S], obs[t], t) for t in tests]), axis=0)
        return logmean(ll), float(np.mean(z <= 2.0)), float(np.min(z)), int(np.argmin(z))

    lab = {"jiangnan": [t for t in scored if t != "T7"], "btbu": [t for t in scored if t == "T7"], "all": scored}
    out["joint"] = {}
    for name, tests in lab.items():
        if not tests:
            continue
        out["joint"][name] = {}
        for rv in RIVALS:
            lm, fs, zb, ib = joint(rv, tests)
            out["joint"][name][rv] = {"log_marginal": lm, "feasible_share": fs, "best_draw_max_z": zb}
            if name == "all":
                out["joint"][name][rv]["best_draw_params"] = {k: float(f"{v:.4g}") for k, v in draws[rv][ib].items()}
                out["joint"][name][rv]["best_draw_z_by_test"] = {
                    t: float(zmax(preds[rv][t][ib][None], obs[t], t)[0]) for t in tests}
    out["failed_integrations"] = {rv: preds[rv]["_failed"] for rv in RIVALS}
    return out


def search(rival, tests, obs, maxiter=12, popsize=12, seed=3):
    """Feasibility search (not a fit): min over the bands of max |z| on the scored tests."""
    from scipy.optimize import differential_evolution
    dim = len({**COMMON, **CONFIG_KEYS[rival]})

    def f(u):
        r = predict((from_unit(rival, u), rival))
        if "error" in r:
            return 1e3
        return float(max(zmax(np.asarray(r[t])[None], obs[t], t)[0] for t in tests))

    res = differential_evolution(f, [(0, 1)] * dim, maxiter=maxiter, popsize=popsize, seed=seed, polish=False,
                                 workers=2, updating="deferred", tol=1e-3)
    p = from_unit(rival, res.x)
    r = predict((p, rival))
    return {"best_max_z": float(res.fun), "evaluations": int(res.nfev),
            "z_by_test": {t: float(zmax(np.asarray(r[t])[None], obs[t], t)[0]) for t in tests},
            "params": {k: float(f"{v:.4g}") for k, v in p.items()},
            "pred": {t: np.asarray(r[t]).round(3).tolist() for t in tests}}


def main():
    n_score = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    n_truth = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    obs, draws, preds = run(n_score, n_truth)
    out = analyse(obs, draws, preds, n_score, n_truth)
    out["n_score"], out["n_truth"] = n_score, n_truth
    np.savez_compressed(OUT / "preds.npz", **{f"{rv}__{t}": preds[rv][t] for rv in RIVALS for t in TESTS})
    (OUT / "draws.json").write_text(json.dumps(draws))
    (OUT / "results.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("observed", "scored_tests", "dropped_tests", "joint", "feasibility")}, indent=1))


if __name__ == "__main__":
    main()
