#!/usr/bin/env python
"""
pyrazine_rivals.py -- power checks and feasibility scoring for the pyrazine rivals P0 / P1 / P1x / P2 / NC.

Design: prereg.yaml in this directory (written before this file; amendment_1 before any run).
Closed-form, constant-free signatures (Cauchy-Schwarz pairing bounds, a 15N label identity, label shares, mass
balances, signs); no ODE, no fit. A2-contract shape: CONFIG_KEYS (bands + basis), rival observables, positive /
negative controls, an INCONCLUSIVE outcome. Reads nothing from src/ or data/: every number below is re-typed from a
dossier (file:line in the comments) so the script is self-contained.

    python pyrazine_rivals.py --out <dir> [--draws 4000] [--seed 20261031]
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

MW = {"PZ": 80.088, "MPZ": 94.115, "DMP": 108.141}
LOG10 = math.log10


@dataclass(frozen=True)
class ConfigKey:
    name: str
    default: Any
    band: Tuple
    basis: str


CONFIG_KEYS = (
    ConfigKey("pyr.h_linear", 1.0, (0.5, 1.0), "linear product courses (Leahy ch.7; Zhou 2024 T2)"),
    ConfigKey("pyr.h_any", 1.0, (0.1, 1.0), "no time course printed"),
    ConfigKey("pyr.log10_u", 0.0, (0.0, 3.0), "aldehyde/ketone imine selectivity, class (declared)"),
    ConfigKey("pyr.log10_e", -9.0, (-2.0, math.log10(9.0)), "P2 extra/pairing MPZ flux (uninformed)"),
    ConfigKey("op.rf_calibrated_dex", 0.05, (0.0, 0.08), "per-compound RF sd, calibrated"),
    ConfigKey("op.rf_area_dex", 0.0, (-0.7, 0.7), "area-R bias, uncalibrated (Ho lab slopes give x3.4-4.2)"),
    ConfigKey("op.near_loq_dex", 0.30, (0.15, 0.40), "near-intercept compound sd"),
    ConfigKey("op.iso_sd_pt", 1.0, (0.5, 2.0), "isotopomer fraction sd, points"),
    ConfigKey("op.mh_leak_mpz", 0.0, (0.0, 0.30), "uncorrected [M-H]+ leakage, PRIMARY band"),
    ConfigKey("op.mh_leak_mpz_corrected", 0.0, (0.0, 0.05), "SECONDARY band (correction assumed)"),
    ConfigKey("op.mh_leak_pz", 0.0, (0.0, 0.05), "parent pyrazine"),
    ConfigKey("pyr.p2_g", 0.3, (0.0, 1.0), "P2b extra-route N label"),
    ConfigKey("pyr.p2_rho", 0.0, (0.0, 1.0), "P2b within-molecule N correlation"),
    ConfigKey("pyr.fA_sd_time", 0.0, (0.0, 0.03), "A-pool label drift allowed under P1"),
    ConfigKey("pyr.fM_mean", 0.45, (0.1, 0.7), "M-pool label (free)"),
    ConfigKey("pyr.fM_sd", 0.0, (0.0, 0.2), "M-pool label heterogeneity"),
    ConfigKey("strecker.log10_c_mgo", -1.15, (-2.5, 0.5), "MGO lifetime, kappa D0 /min (forced-input layer)"),
    ConfigKey("strecker.transfer_aa_dex", 0.0, (-0.5, 0.5), "Phe <-> Ala transfer"),
    ConfigKey("lab.between_dex", 0.5, (0.3, 0.7), "between-lab sd of a Strecker-type yield (CF-073)"),
)
BAND = {c.name: c.band for c in CONFIG_KEYS}
DEFAULT = {c.name: c.default for c in CONFIG_KEYS}


def R_molar(pz, mpz, dmp, unit="mass"):
    """R = MPZ^2 / (PZ x DMP), molar. unit 'mass' converts mass concentrations / mass % by MW."""
    if unit == "mass":
        pz, mpz, dmp = pz / MW["PZ"], mpz / MW["MPZ"], dmp / MW["DMP"]
    return mpz * mpz / (pz * dmp)


def r_max_from_m_prime(x: float, dmp_includes_26: bool) -> float:
    """Shared-constant bound. A, M, M' pools; DMP25 = M^2 + M'^2, DMP26 = 2 M M' (x = M'/M)."""
    if dmp_includes_26:
        return 4.0
    return 4.0 * (1.0 + x) ** 2 / (1.0 + x * x)


def x_from_26_25(ratio: float) -> float:
    """Solve 2x/(1+x^2) = ratio for x <= 1."""
    if ratio <= 0:
        return 0.0
    return (1.0 - math.sqrt(max(0.0, 1.0 - ratio * ratio))) / ratio


def slope_rel_se(r2: float, n: int) -> float:
    return math.sqrt((1.0 - r2) / (r2 * (n - 2)))


# ============================================================================================ S1 / S2 datasets
# Each row: lab, id, label, R_obs, R_max, sigma_pot (dex, per-pot noise excluding the lab-shared RF term),
# shared (dex sd of the lab-shared RF term on log R, or a symmetric bias half-width for uncalibrated sets), kind.
def leahy_rows() -> List[Dict[str, Any]]:
    rows = []
    # leahy1989_extraction.md sec. 3 Table I, umol L-1 min-1 column; r2, n per row.
    T = (95, 85, 75)
    systems = {
        # PZ, MPZ, DMP25, (DMP26 or None); each (k, r2, n) per T
        "Lys-Glc": ([(.748, .994, 22), (.102, .960, 22), (.0445, .965, 20)],
                    [(.502, .995, 22), (.0747, .967, 22), (.0282, .941, 20)],
                    [(.0287, .995, 20), (.00381, .985, 22), (.00103, .942, 16)], None),
        "Lys-Fru": ([(.283, .995, 20), (.0822, .995, 22), (.0279, .963, 21)],
                    [(.196, .995, 20), (.0687, .945, 22), (.0205, .968, 21)],
                    [(.0270, .945, 20), (.00580, .951, 20), (.00120, .935, 15)], None),
        "Lys-Rib": ([(.726, .974, 21), (.206, .916, 22), (.0645, .882, 22)],
                    [(.950, .957, 19), (.313, .929, 22), (.0779, .952, 22)],
                    [(.0234, .975, 18), (.00678, .933, 17), (.00197, .984, 18)], "x_endpoint"),
        "Asn-Glc": ([(.0214, .987, 22), (.00878, .941, 22), (.00204, .947, 22)],
                    [(.0783, .992, 22), (.0317, .926, 22), (.00607, .942, 22)],
                    [(.0134, .982, 22), (.00311, .968, 22), (.000701, .863, 14)],
                    [(.0154, .979, 21), (.00290, .821, 11), None]),
        "Asn-Fru": ([(.00472, .904, 20), (.00136, .881, 20), (.000554, .947, 20)],
                    [(.112, .927, 22), (.0289, .926, 22), (.00662, .896, 22)],
                    [(.0751, .983, 22), (.0188, .995, 22), (.00402, .940, 20)],
                    [(.0834, .960, 20), (.0234, .980, 20), (.00549, .921, 18)]),
    }
    # Lys systems: DMP26 not detected; LOD <= 3 ppb (text) against DMP25 >= ~0.16 ppm over the run -> x <= 0.02.
    # Lys-Rib: Table III end point DMP26 0.3 % vs DMP25 1.6 % -> x from 2x/(1+x^2) = 0.1875.
    for sys_, (pz, mpz, d25, d26) in systems.items():
        for i, t in enumerate(T):
            if d26 is None or d26 == "x_endpoint":
                x = 0.02 if d26 is None else x_from_26_25(0.3 / 1.6)
                dmp, rmax, d26_se = d25[i][0], r_max_from_m_prime(x, False), 0.0
            else:
                if d26[i] is None:
                    continue           # Asn-Glc 75 C: no DMP26 regression printed (not a non-detect statement)
                dmp, rmax = d25[i][0] + d26[i][0], 4.0
                d26_se = slope_rel_se(d26[i][1], d26[i][2]) * d26[i][0] / dmp
            R = mpz[i][0] ** 2 / (pz[i][0] * dmp)
            se = math.sqrt(slope_rel_se(*pz[i][1:]) ** 2 + 4 * slope_rel_se(*mpz[i][1:]) ** 2
                           + (slope_rel_se(*d25[i][1:]) * d25[i][0] / dmp) ** 2 + d26_se ** 2) / math.log(10)
            rows.append(dict(lab="Reineccius", id="leahy1989 T-I", label=f"{sys_} {t} C rate", R=R, Rmax=rmax,
                             sig_pot=se, kind="calibrated"))
    # Table III end points, 2 h 95 C (mass %), leahy1989_extraction.md sec. 3; leahy1989a Table III pH 7.
    endp = {"Lys-Glc": (55.8, 41.5, 2.4, 0.0, 0.02), "Lys-Fru": (42.3, 48.6, 7.3, 0.0, 0.02),
            "Lys-Rib": (36.8, 61.1, 1.6, 0.3, None), "Asn-Glc": (18.7, 58.8, 17.3, 4.8, None),
            "Asn-Fru": (1.0, 29.1, 37.6, 32.3, None)}
    for sys_, (pz, mpz, d25, d26, x) in endp.items():
        inc = d26 > 0
        R = R_molar(pz, mpz, d25 + d26)
        rmax = 4.0 if inc else r_max_from_m_prime(x, False)
        rd = 0.05 / 2.0  # rounding of the last printed digit (0.1 %), as relative error on the smallest share
        sig = math.sqrt(sum((0.05 / v) ** 2 for v in (pz, mpz, d25 + d26) if v > 0)) / math.log(10)
        rows.append(dict(lab="Reineccius", id="leahy1989 T-III", label=f"{sys_} 2 h 95 C end point", R=R,
                         Rmax=rmax, sig_pot=sig, kind="calibrated", dependent=True))
    rows.append(dict(lab="Reineccius", id="leahy1989a T-III", label="Lys-Glc pH 7, 2 h 95 C end point",
                     R=R_molar(58.7, 39.6, 1.7), Rmax=r_max_from_m_prime(0.02, False),
                     sig_pot=math.sqrt((0.05 / 58.7) ** 2 + 4 * (0.05 / 39.6) ** 2 + (0.05 / 1.7) ** 2) / math.log(10),
                     kind="calibrated"))
    return rows


def _sig_from_sd(pairs):
    """pairs = [(value, sd, power)] -> sd of log10 R from replicate SDs (relative, first order)."""
    return math.sqrt(sum((p * sd / v) ** 2 for v, sd, p in pairs)) / math.log(10)


def ho_rows() -> List[Dict[str, Any]]:
    rows = []
    loq = DEFAULT["op.near_loq_dex"]
    # zhou2023 Table 1 (zhou2023_extraction.md:157-160), 120 C 60 min, ug/L; DMP26 ND except pH 8 + Cys.
    x_cys = x_from_26_25(0.794 / 21.44)
    rows.append(dict(lab="Ho", id="zhou2023 T1", label="ARP pH 8 (DMP25 near intercept)",
                     R=R_molar(3.99, 3.97, 0.272), Rmax=r_max_from_m_prime(0.05, False),
                     sig_pot=math.sqrt(_sig_from_sd([(3.99, .76, 1), (3.97, .15, 2)]) ** 2 + loq ** 2), kind="calibrated"))
    rows.append(dict(lab="Ho", id="zhou2023 T1", label="ARP + Cys pH 8", R=R_molar(177.05, 99.00, 21.44 + 0.794),
                     Rmax=4.0, sig_pot=_sig_from_sd([(177.05, 8.08, 1), (99.0, 1.96, 2), (22.23, 2.0, 1)]),
                     kind="calibrated"))
    # zhou2022 Table 3, pH 8, 60 min (zhou2022_extraction.md:111-114); DMP26 not reported -> x as zhou2023 + Cys.
    for lab_, pz, mpz, dmp, near in (("ARP", (5.02, .67), (3.66, .59), (0.34, .11), True),
                                     ("ARP + 6 Cys", (19.50, 2.76), (18.30, 2.72), (9.05, .43), False),
                                     ("ARP + 14 Cys", (144.45, 9.74), (77.87, 2.52), (17.19, 1.60), False),
                                     ("ARP + 20 Cys", (169.13, 12.59), (96.30, 5.41), (23.18, 2.10), False)):
        s = _sig_from_sd([(pz[0], pz[1], 1), (mpz[0], mpz[1], 2), (dmp[0], dmp[1], 1)])
        if near:
            s = math.sqrt(_sig_from_sd([(pz[0], pz[1], 1), (mpz[0], mpz[1], 2)]) ** 2 + loq ** 2)
        rows.append(dict(lab="Ho", id="zhou2022 T3", label=f"{lab_} pH 8" + (" (DMP25 near intercept)" if near else ""),
                         R=R_molar(pz[0], mpz[0], dmp[0]), Rmax=r_max_from_m_prime(0.05 if near else x_cys, False),
                         sig_pot=s, kind="calibrated"))
    # zhou2023 Fig. 3 Cys + MGO 20 + 20 mM, figure-read (zhou2023_extraction.md:284-298, 369-383); DMP25 + DMP26.
    fig = {6: ((37.0, 46.3, 47.0), (95, 146, 160), (43, 67, 66), (5.2, 7.0, 7.9)),
           7: ((20.7, 19.4, 20.8), (66, 86, 97), (59, 94, 106), (7.4, 10.5, 12.6)),
           8: ((12.7, 14.4, 15.4), (65, 106, 117), (103, 205, 237), (10.8, 24.6, 25.9))}
    for ph, (pz, mpz, d25, d26) in fig.items():
        for j, t in enumerate((30, 60, 90)):
            dm = d25[j] + d26[j]
            s = math.sqrt(_sig_from_sd([(pz[j], 0.8, 1), (mpz[j], 3, 2), (dm, math.hypot(3, .4), 1)]) ** 2
                          + 6 * 0.05 ** 2 / 2)   # read error + half the calibration term (fig. digitisation)
            rows.append(dict(lab="Ho", id="zhou2023 F3", label=f"Cys+MGO pH {ph} {t} min", R=R_molar(pz[j], mpz[j], dm),
                             Rmax=4.0, sig_pot=s, kind="calibrated", dependent=(t != 90)))
    # zhou2023c Table 1 (zhou2023c_extraction.md:119-122), ARP + GSH / Glu, 60 min; DMP26 not reported.
    for lab_, pz, mpz, dmp in (("ARP + GSH", (8.33, 1.06), (13.1, 2.72), (3.75, .738)),
                               ("ARP + Glu", (8.33, 1.00), (29.5, 2.30), (15.1, 2.48))):
        rows.append(dict(lab="Ho", id="zhou2023c T1", label=lab_, R=R_molar(pz[0], mpz[0], dmp[0]),
                         Rmax=r_max_from_m_prime(0.1, False),
                         sig_pot=_sig_from_sd([(pz[0], pz[1], 1), (mpz[0], mpz[1], 2), (dmp[0], dmp[1], 1)]),
                         kind="calibrated"))
    return rows


def uncalibrated_rows() -> List[Dict[str, Any]]:
    rows = []
    # hwang1995 Table 1 (hwang1995_extraction.md:101-103), ug/g glucose, DMP25(6) combined -> Rmax 4.
    hw = {"ctrl": (6.71, 146.77, 132.31), "Gln": (41.99, 279.85, 115.47), "Lys": (99.94, 1253.43, 162.71),
          "Asn": (78.85, 226.53, 260.33), "Phe": (31.65, 366.44, 571.07), "Glu": (75.09, 253.66, 249.77),
          "Asp": (25.91, 147.92, 112.30), "Ile": (8.43, 148.95, 140.69), "Arg": (11.49, 61.69, 99.47)}
    for k, (a, b, c) in hw.items():
        rows.append(dict(lab="Ho(Rutgers) uncal", id="hwang1995 T1", label=f"Glc + 15N-Gly + {k}, starch 180 C",
                         R=R_molar(a, b, c), Rmax=4.0, sig_pot=0.10 * math.sqrt(6), kind="uncalibrated"))
    # martin2001 (martin2001_extraction.md:297-306), peak areas, DMP25(6) combined.
    for k, R in {"U": 19.6, "G": 6.0, "GA": 9.27, "GGn": 17.2, "GL": 4.83, "GI": 2.78, "GP": 5.59, "GM": 3.71,
                 "X": 8.72}.items():
        rows.append(dict(lab="Reading uncal", id="martin2001", label=f"potato {k}", R=R, Rmax=4.0,
                         sig_pot=0.10 * math.sqrt(6), kind="uncalibrated"))
    # deng2022 (deng2022_extraction.md:95-142), RF = 1, DMP25 only (x <= 0.1 declared).
    for k, R in {"Met/Glc 120 C 60": 0.52, "Met/Glc 120 C 240": 1.53, "ARP 120 C max": 5.71, "ARP 130 C max": 8.78}.items():
        rows.append(dict(lab="Ho uncal", id="deng2022", label=k, R=R, Rmax=r_max_from_m_prime(0.1, False),
                         sig_pot=0.10 * math.sqrt(6), kind="uncalibrated"))
    for k, R in {"vanlancker2012 max (Gly+Lys Glc)": 0.77, "shu1999 120 C": 3.13}.items():
        rows.append(dict(lab="other uncal", id=k.split()[0], label=k, R=R, Rmax=4.0, sig_pot=0.10 * math.sqrt(6),
                         kind="uncalibrated"))
    return rows


def score_bound_rows(rows, shared_sd_cal: float, bias_half_uncal: float) -> List[Dict[str, Any]]:
    """z per row: calibrated -> (log R/Rmax) / sqrt(sig_pot^2 + shared^2); uncalibrated -> (log R/Rmax - bias_max)/sig_pot."""
    out = []
    for r in rows:
        ex = LOG10(r["R"] / r["Rmax"])
        if r["kind"] == "calibrated":
            z = ex / math.sqrt(r["sig_pot"] ** 2 + shared_sd_cal ** 2)
        else:
            z = (ex - bias_half_uncal) / r["sig_pot"]
        out.append({**r, "excess_dex": ex, "z": z, "rejects_shared": z > 2.0})
    return out


def lab_level(scored, lab: str) -> Dict[str, Any]:
    """RF errors are shared inside a lab: the null survives only if ONE shared log-R bias rescues every row.
    Statistic: min over independent rows of the excess, against sqrt(shared^2 + the row's own sigma^2)."""
    rows = [r for r in scored if r["lab"] == lab and not r.get("dependent") and r["excess_dex"] > 0]
    if not rows:
        return {"lab": lab, "n_rows_above_bound": 0}
    ex_sorted = sorted(rows, key=lambda r: r["z"])
    return {"lab": lab, "n_rows_above_bound": len(rows), "n_rows_z_gt_2": sum(r["z"] > 2 for r in rows),
            "weakest_rejecting_row": ex_sorted[-1]["label"] if ex_sorted[-1]["z"] > 2 else None,
            "max_z": max(r["z"] for r in rows)}


# ============================================================================================ power, S1 / S2
def power_bound(rng, rmax: float, sig_pot: float, shared: float, n: int, kind: str, bias_half: float):
    """Null at its least favourable point (h = 1, R = Rmax); alternatives P1x / P2 (band-averaged and material)."""
    sig = math.sqrt(sig_pot ** 2 + shared ** 2) if kind == "calibrated" else sig_pot

    def observe(logR_true):
        obs = logR_true + rng.normal(0, sig, logR_true.shape)
        if kind == "uncalibrated":
            obs = obs + rng.uniform(-bias_half, bias_half, logR_true.shape)
            return (obs - LOG10(rmax) - bias_half) / sig_pot
        return (obs - LOG10(rmax)) / sig

    h = rng.uniform(*BAND["pyr.h_linear"], n)
    null = observe(np.full(n, LOG10(rmax)))
    u = 10 ** rng.uniform(*BAND["pyr.log10_u"], n)
    p1x = observe(np.log10(h * (1 + u) ** 2 / u))
    e = 10 ** rng.uniform(*BAND["pyr.log10_e"], n)
    p2 = observe(np.log10(4 * h * (1 + e) ** 2))
    # material: rival doubles R relative to the bound (prereg amendment_1)
    u_m = 10 ** rng.uniform(LOG10(5.8), BAND["pyr.log10_u"][1], n)
    e_m = 10 ** rng.uniform(LOG10(0.41), BAND["pyr.log10_e"][1], n)
    p1x_m = observe(np.log10(h * (1 + u_m) ** 2 / u_m))
    p2_m = observe(np.log10(4 * h * (1 + e_m) ** 2))
    f = lambda z: float(np.mean(z > 2.0))
    return {"size_null_h1": f(null), "power_P1x_band": f(p1x), "power_P2_band": f(p2),
            "power_P1x_material": f(p1x_m), "power_P2_material": f(p2_m)}


# ============================================================================================ S3 15N identity
# zhou2024 Table 1 (zhou2024_extraction.md sec. 3): isotopomer % with 0 / 1 / 2 15N.
Z24 = {"PZ": {30: (20.26, 50.54, 29.20), 60: (20.97, 50.06, 28.97), 90: (21.26, 50.03, 28.71), 120: (21.18, 48.97, 29.85)},
       "MPZ": {30: (39.73, 46.29, 13.98), 60: (37.75, 47.53, 14.72), 90: (35.62, 47.95, 16.43), 120: (33.94, 48.00, 18.06)},
       "DMP": {30: (35.97, 41.56, 22.47), 60: (33.34, 45.51, 21.15), 90: (32.11, 44.48, 23.41), 120: (31.60, 44.37, 24.03)}}


def D_stat(pz, mpz):
    pz = np.asarray(pz, float); mpz = np.asarray(mpz, float)
    pz = pz / pz.sum(-1, keepdims=True); mpz = mpz / mpz.sum(-1, keepdims=True)
    fA = (pz[..., 1] + 2 * pz[..., 2]) / 2
    return mpz[..., 2] / fA + mpz[..., 0] / (1 - fA) - 1.0, fA


def leak(t, r):
    """[M-H]+ leakage: each isotopologue puts a fraction r/(1+r) of its intensity one mass lower; M-1 of the
    unlabelled species is dropped (outside the cluster); renormalise."""
    t0, t1, t2 = t[..., 0], t[..., 1], t[..., 2]
    a = 1.0 / (1.0 + r); b = r / (1.0 + r)
    o = np.stack([t0 * a + t1 * b, t1 * a + t2 * b, t2 * a], -1)
    return o / o.sum(-1, keepdims=True)


def synth_S3(rng, n, rival, leak_band, material=False):
    fA = rng.uniform(0.45, 0.65, n)                       # free-amine label (any; the statistic is invariant)
    sdA = rng.uniform(*BAND["pyr.fA_sd_time"], n)
    fA_m = np.clip(fA + rng.normal(0, 1, n) * sdA, 0.01, 0.99)
    fM = rng.uniform(*BAND["pyr.fM_mean"], n)
    sdM = rng.uniform(*BAND["pyr.fM_sd"], n)
    c = rng.uniform(-1, 1, n) * sdA * sdM
    pz = np.stack([(1 - fA) ** 2, 2 * fA * (1 - fA), fA ** 2], -1)
    pairing = np.stack([(1 - fA_m) * (1 - fM) + c, fA_m * (1 - fM) + (1 - fA_m) * fM - 2 * c, fA_m * fM + c], -1)
    mpz = pairing
    if rival == "P2":
        e = 10 ** rng.uniform(*BAND["pyr.log10_e"], n)
        w = e / (1 + e)
        g = rng.uniform(*BAND["pyr.p2_g"], n)
        rho = rng.uniform(*BAND["pyr.p2_rho"], n)
        if material:
            w = rng.uniform(0.3, 0.9, n)
            # |g - pairing label| >= 0.2 or rho >= 0.5
            pair_lab = (fA_m + fM) / 2
            ok = (np.abs(g - pair_lab) >= 0.2) | (rho >= 0.5)
            while not ok.all():
                k = ~ok
                g[k] = rng.uniform(0, 1, k.sum()); rho[k] = rng.uniform(0, 1, k.sum())
                ok = (np.abs(g - pair_lab) >= 0.2) | (rho >= 0.5)
        extra = np.stack([(1 - g) ** 2 + rho * g * (1 - g), 2 * g * (1 - g) * (1 - rho), g ** 2 + rho * g * (1 - g)], -1)
        mpz = (1 - w)[:, None] * pairing + w[:, None] * extra
    r_m = rng.uniform(*leak_band, n)
    r_p = rng.uniform(*BAND["op.mh_leak_pz"], n)
    sd = rng.uniform(*BAND["op.iso_sd_pt"], n)[:, None] / 100.0
    pz_o = np.clip(leak(pz, r_p) + rng.normal(0, 1, pz.shape) * sd, 1e-4, None)
    mpz_o = np.clip(leak(mpz, r_m) + rng.normal(0, 1, mpz.shape) * sd, 1e-4, None)
    D, _ = D_stat(pz_o, mpz_o)
    return D


def s3(rng, n):
    out = {}
    obs = {t: float(D_stat(Z24["PZ"][t], Z24["MPZ"][t])[0]) for t in (30, 60, 90, 120)}
    out["D_obs"] = obs
    # P1 must also pass its own PZ binomial check (positive control on the PZ rows): report the PZ binomial residual
    pz30 = np.array(Z24["PZ"][30]) / 100; fA30 = (pz30[1] + 2 * pz30[2]) / 2
    out["PZ_binomial_residual_pt_30"] = list(np.round(100 * (pz30 - np.array([(1 - fA30) ** 2, 2 * fA30 * (1 - fA30), fA30 ** 2])), 2))
    d30 = np.array(Z24["DMP"][30]) / 100; fM30 = (d30[1] + 2 * d30[2]) / 2
    out["DMP_binomial_residual_pt_30"] = list(np.round(100 * (d30 - np.array([(1 - fM30) ** 2, 2 * fM30 * (1 - fM30), fM30 ** 2])), 2))
    for tag, band in (("primary", BAND["op.mh_leak_mpz"]), ("secondary", BAND["op.mh_leak_mpz_corrected"])):
        Dn = synth_S3(rng, n, "P1", band)
        crit = float(np.quantile(Dn, 0.95))
        crit2 = float(np.quantile(Dn, 0.05))   # two-sided report (P2 can push D either way)
        Dp = synth_S3(rng, n, "P2", band)
        Dpm = synth_S3(rng, n, "P2", band, material=True)
        rej = lambda D: (D > crit) | (D < crit2)
        out[tag] = {"leak_band": band, "D_crit_95": crit, "D_crit_05": crit2,
                    "size_two_sided": float(np.mean(rej(Dn))),
                    "power_P2_band": float(np.mean(rej(Dp))), "power_P2_material": float(np.mean(rej(Dpm))),
                    "obs_rejects_single_pool": {t: bool(obs[t] > crit or obs[t] < crit2) for t in obs},
                    "obs_pctile_under_P1_30": float(np.mean(Dn < obs[30]))}
        # worst case inside the band: leakage at the band top, sd at the band bottom
        r = band[1]
        fA = fA30; Dw = []
        for fM in np.linspace(0.1, 0.7, 61):
            t = np.array([(1 - fA) * (1 - fM), fA * (1 - fM) + (1 - fA) * fM, fA * fM])
            Dw.append(float(D_stat(pz30 * 100, leak(t, r) * 100)[0]))
        out[tag]["D_null_max_at_band_top_leak"] = max(Dw)
    # P2 extra share needed to make D = D_obs(30) with rho = 0 / rho = 1 (diagnostic)
    return out


# ============================================================================================ S4 labels
def s4(rng, n):
    obs = {"amranihemaimi1995 MPZ (Ala, Gly)": 0.00, "amranihemaimi1995 DMP25": 0.00, "low2007 MPZ (Gly)": 0.07,
           "low2007 DMP25": 0.063}
    sd = 0.04   # label-fraction sd: Low's +-15 points is on derived shares; a direct label fraction ~ +-0.03-0.05
    pred = {"P0": (0.0, 0.0), "P1": (0.0, 0.0), "P1x": (0.0, 0.0), "P2a_sugar": (0.0, 0.0), "P2b": (0.0, 0.0),
            "NC": (0.4, 1.0)}
    res = {}
    for rv, (lo, hi) in pred.items():
        p = rng.uniform(lo, hi, n)
        res[rv] = {"pred_band": [lo, hi], "passes_obs": bool(lo - 2 * sd <= 0.07)}
    # P2a_aa at material w >= 0.3 with phi_HCHO in [0.5, 1] (Low: 2,6-DMP ~45 % HCHO-labelled in the same pot)
    w = rng.uniform(0.3, 0.9, n); phi = rng.uniform(0.5, 1.0, n)
    lab = w * phi + rng.normal(0, sd, n)
    res["P2a_aa"] = {"pred_band_material": [0.15, 0.9], "power_vs_obs_0.07": float(np.mean(lab > 0.07 + 2 * sd)),
                     "passes_obs": False, "max_w_allowed_phi_0.5": (0.07 + 2 * sd) / 0.5,
                     "max_w_allowed_phi_1": 0.07 + 2 * sd}
    nc = rng.uniform(0.4, 1.0, n) + rng.normal(0, sd, n)
    res["NC"]["power"] = float(np.mean(nc > 0.07 + 2 * sd))
    res["obs"] = obs
    return res


# ============================================================================================ S5 Adams
ADAMS = [  # adams2008 Table 3: (label, pyrrole mol %, sd, DMP mol %, sd)
    ("130 C 1:1 phosphate pH 7", 0.020, 0.0073, 0.02, None),
    ("130 C 1:1 bicarbonate pH 7", 0.055, 0.0094, 0.011, 0.00064),
    ("130 C 1:1 acetate pH 4", 1.6, 0.040, 0.042, 0.005),
    ("100 C 1:1 acetate pH 4 30 min", 1.1, 0.087, 0.020, 0.00025),
    ("100 C 1:1 acetate pH 4 60 min", 1.3, 0.045, 0.026, 0.0030),
    ("100 C 2:1 acetate pH 4 30 min", 0.50, 0.019, 0.013, 0.0029),
    ("100 C 1:3 acetate pH 4 30 min", 3.7, 0.074, 0.016, 0.0040),
]


def s5(rng, n):
    op = 0.15
    rows = []
    for lab, py, spy, dm, sdm in ADAMS:
        B = LOG10(py / (2 * dm))
        sig = math.sqrt(op ** 2 * 2 + (spy / py / math.log(10)) ** 2 + ((sdm or 0.2 * dm) / dm / math.log(10)) ** 2)
        rows.append({"row": lab, "B_dex": B, "sigma": sig, "rejects_P0": B - 2 * sig > 0})
    # power: P1 with sink >= 3x pyrazine flux (material), synthetic B_obs
    sig = 0.25
    Bt = rng.uniform(0.48, 3.0, n)
    power = float(np.mean(Bt + rng.normal(0, sig, n) - 2 * sig > 0))
    # size: P0 cannot make pyrrole from its Strecker flux; at best a blank of <= 10 % of DMP (B <= -1.3)
    size = float(np.mean(-1.3 + rng.normal(0, sig, n) - 2 * sig > 0))
    ratio = {r[0]: r[1] / r[3] for r in ADAMS[5:7]}
    ratio["100 C 1:1 acetate pH 4 30 min"] = ADAMS[3][1] / ADAMS[3][3]
    order = [ratio["100 C 2:1 acetate pH 4 30 min"], ratio["100 C 1:1 acetate pH 4 30 min"],
             ratio["100 C 1:3 acetate pH 4 30 min"]]
    return {"rows": rows, "power_material": power, "size": size,
            "S5b_pyrrole_over_DMP_2:1,1:1,1:3": order, "S5b_monotone_rising": bool(order[0] < order[1] < order[2])}


# ============================================================================================ S6 Strecker vs P0
def p0_log10k(donor: str, T_c: float, ph: float) -> float:
    """B18 frozen constants (parameters_pyrazine.py FROZEN_B18), its Arrhenius and pH law (knot 7, ref 6.8)."""
    lk, ea = {"MGO": (-7.529714282148093, 114.9), "GO": (-6.5415452891122525, 103.1)}[donor]
    s_hi, s_lo = 0.19744840125387295, 0.5795142879227316
    shape = lambda x: s_hi * (x - 7) if x >= 7 else s_lo * (x - 7)
    R = 8.314e-3
    return lk - ea / (R * math.log(10)) * (1 / (T_c + 273.15) - 1 / 373.15) + shape(ph) - shape(6.8)


def s6(rng, n):
    out = {}
    obs = {"MGO": (11.1 + 4.7) / 10.0, "GO": (9.2 + 2.2) / 10.0}   # hofmann2000 Table 4 argon, mol % of Phe
    off = {"MGO": 0.0, "GO": (-0.5, 0.4)}
    for d in ("MGO", "GO"):
        lk = p0_log10k(d, 98.0, 7.0)
        lo_c, hi_c = BAND["strecker.log10_c_mgo"]
        # integral of D over 60 min under the second-order donor sink (D0 = 100 mM); max yield at the slowest sink
        def integ(log10_c):
            c = 10 ** log10_c
            return np.log1p(60 * c) / c          # in units of D0 x min
        c_grid = np.linspace(lo_c + (off[d][0] if d == "GO" else 0), hi_c + (off[d][1] if d == "GO" else 0), 400)
        y_max = 10 ** lk * 100.0 * integ(c_grid).max() * 10 ** BAND["strecker.transfer_aa_dex"][1] * 100  # mol %
        y_def = 10 ** lk * 100.0 * float(integ(DEFAULT["strecker.log10_c_mgo"])) * 100
        G_min = LOG10(obs[d] / y_max)
        thr = {k: 2 * v * math.sqrt(2) for k, v in (("lab_0.3", 0.3), ("lab_0.5", 0.5), ("lab_0.7", 0.7))}
        out[d] = {"obs_mol_pct": obs[d], "P0_yield_max_mol_pct": y_max, "P0_yield_default_mol_pct": y_def,
                  "G_min_dex": G_min, "G_default_dex": LOG10(obs[d] / y_def), "thresholds": thr,
                  "rejects_P0": {k: G_min > v for k, v in thr.items()}}
    # power at default between-lab sd: P1 with Y <= 0.1 (G_true >= 1) vs the rule at G_min; size for P0
    sig = 0.5 * math.sqrt(2)
    slack = out["MGO"]["G_default_dex"] - out["MGO"]["G_min_dex"]      # band slack the rule gives P0
    Gt = rng.uniform(1.0, 3.0, n)
    out["power_P1_Yle0.1_lab0.5"] = float(np.mean(Gt + rng.normal(0, sig, n) - slack > 2 * sig))
    out["size_P0_lab0.5"] = float(np.mean(rng.normal(0, sig, n) + rng.uniform(-0.5, 0.5, n) - slack > 2 * sig))
    return out


# ============================================================================================ S7 Adams pH
def s7(rng, n):
    obs = LOG10(0.042 / 0.02)
    pred = p0_log10k("MGO", 130, 4.0) - p0_log10k("MGO", 130, 7.0)
    sig = math.sqrt(2 * 0.15 ** 2 + 0.3 ** 2)
    p0 = pred + rng.normal(0, sig, n)
    p1 = rng.uniform(-1, 1, n) + rng.normal(0, sig, n)
    rej = lambda x: np.abs(x - pred) > 2 * sig
    return {"obs_dex": obs, "P0_pred_dex": pred, "sigma": sig, "rejects_P0": bool(abs(obs - pred) > 2 * sig),
            "size": float(np.mean(rej(p0))), "power_P1_material": float(np.mean(rej(p1)))}


# ============================================================================================ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20261031)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out_dir = Path(a.out); out_dir.mkdir(parents=True, exist_ok=True)
    shared_cal = DEFAULT["op.rf_calibrated_dex"] * math.sqrt(6)       # 2 x MPZ + PZ + DMP, per-compound sd 0.05
    shared_cal_top = BAND["op.rf_calibrated_dex"][1] * math.sqrt(6)
    bias_half = BAND["op.rf_area_dex"][1]
    rows = leahy_rows() + ho_rows() + uncalibrated_rows()
    res: Dict[str, Any] = {"config_keys": [[c.name, c.default, list(c.band), c.basis] for c in CONFIG_KEYS]}
    for tag, sh in (("rf_default", shared_cal), ("rf_band_top", shared_cal_top)):
        sc = score_bound_rows(rows, sh, bias_half)
        for r in sc:
            r["power"] = power_bound(rng, r["Rmax"], r["sig_pot"], sh, a.draws, r["kind"], bias_half)
            r["powered"] = r["power"]["size_null_h1"] <= 0.1 and min(r["power"]["power_P1x_material"],
                                                                    r["power"]["power_P2_material"]) >= 0.6
        res[f"S1S2_{tag}"] = sc
        res[f"S1S2_{tag}_lab"] = [lab_level([r for r in sc if r["powered"]], lab) for lab in
                                  ("Reineccius", "Ho", "Ho(Rutgers) uncal", "Reading uncal", "Ho uncal", "other uncal")]
    res["S3"] = s3(rng, a.draws * 5)
    res["S4"] = s4(rng, a.draws)
    res["S5"] = s5(rng, a.draws)
    res["S6"] = s6(rng, a.draws)
    res["S7"] = s7(rng, a.draws)
    (out_dir / "pyrazine_rivals_result.json").write_text(json.dumps(res, indent=1, default=float))
    write_md(res, out_dir / "pyrazine_rivals_result.md")
    print((out_dir / "pyrazine_rivals_result.md").read_text())


def write_md(res, path: Path):
    L = ["# Pyrazine rivals: power and feasibility (generated by pyrazine_rivals.py; design prereg.yaml)", ""]
    for tag in ("rf_default", "rf_band_top"):
        L += [f"## S1 / S2 mixing bound, shared-RF {tag}", "",
              "| lab | source | pot | R | R_max | excess dex | z | size | pow P1x mat | pow P2 mat | pow P2 band | powered | rejects shared |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in res[f"S1S2_{tag}"]:
            p = r["power"]
            L.append(f"| {r['lab']} | {r['id']} | {r['label']}{' (dep.)' if r.get('dependent') else ''} | {r['R']:.2f} | {r['Rmax']:.2f} | "
                     f"{r['excess_dex']:+.2f} | {r['z']:.1f} | {p['size_null_h1']:.3f} | {p['power_P1x_material']:.2f} | "
                     f"{p['power_P2_material']:.2f} | {p['power_P2_band']:.2f} | {r['powered']} | {r['rejects_shared'] and r['powered']} |")
        L += ["", "Lab level (powered rows only, independent rows): " + json.dumps(res[f"S1S2_{tag}_lab"], default=float), ""]
    L += ["## S3 15N identity (zhou2024 Table 1)", "", "```", json.dumps(res["S3"], indent=1, default=float), "```", ""]
    for k in ("S4", "S5", "S6", "S7"):
        L += [f"## {k}", "", "```", json.dumps(res[k], indent=1, default=float), "```", ""]
    path.write_text("\n".join(L))


if __name__ == "__main__":
    main()
