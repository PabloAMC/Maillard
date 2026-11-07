#!/usr/bin/env python
"""
m5_power_check.py -- the M5 case's reduced power check (case sec. 5.3; coordinator rulings 2026-11-06).

SYNTHETIC DATA ONLY. No measured value of any unit is read and no reserved data is touched: every truth is drawn from
the case's declared bands, every dataset is simulated, and the script imports nothing from the engine or from data/.
It sizes the case's DESIGNED deciders (sec. 4.4), one reduced forced-input mini-model per part:

  B  buffer catalysis: step rates read directly (forced input), 5 phosphate levels (inert zero, 5, 50, 200, 500 mM)
     x 2 pH (5.5, 7.0) x 4 readouts (Amadori formation F, Amadori decay D, 2,3-enolisation E, furfural branch U).
  F  lipid x thiol: a two-phase production-loss pool model for MFT and FFT over one cook, six arms (buffer; fresh oil
     under air; fresh oil under N2; enal dosed without oil; PV-spiked oil under N2; enal dosed into oil under N2),
     readouts per phase plus the two product classes, plus a direct K_oil/water measurement.
  W  water activity: acrylamide formation constant k_F at 5 moistures (10 - 70 % wb), moisture measured directly.

The decision rule is the case's (sec. 4.0, coordinator ruling 1): a rival is FALSIFIED only if its best fit is
inadequate (p < 0.05) under BOTH measurement-model forms:
  P  studentized per point: chi2 of t-to-normal-score residuals, each residual scaled by its own replicate SE
     (n - 1 dof per point), fit weighted by those SEs;
  Q  pooled variance: lack-of-fit F against the pooled within-cell variance (m (n - 1) dof), equal weights.
A pair (A, B) is SEPARATED on a dataset simulated from A when B is falsified and A is not. One designed dataset per
simulation, so Holm has m = 1 here.

Usage (through the wrapper; one core):
    python scripts/diagnostics/m5_power_check.py                  # all parts, writes results/validation/m5_power_check.json
    python scripts/diagnostics/m5_power_check.py --part B --nsim 40 --out /tmp/x.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
from scipy import stats
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "results" / "validation" / "m5_power_check.json"
ALPHA = 0.05


# ----------------------------------------------------------------------------- the two-form test

def _fit(resid_fn: Callable[[np.ndarray], np.ndarray], starts: Sequence[np.ndarray], lo, hi) -> np.ndarray:
    best, cost = None, np.inf
    for x0 in starts:
        try:
            r = least_squares(resid_fn, np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), method="trf", max_nfev=1500)
        except ValueError:
            continue
        if r.cost < cost:
            best, cost = r.x, r.cost
    return best


def two_form_p(Y: np.ndarray, predict: Callable[[np.ndarray], np.ndarray], starts, lo, hi) -> Tuple[float, float]:
    """Y: (m cells, n replicates) on the log scale. Returns (p_P, p_Q) for the rival's best fit under each form."""
    m, n = Y.shape
    q = len(lo)
    dof = m - q
    if dof <= 0:
        return 1.0, 1.0
    ybar = Y.mean(1)
    s = np.maximum(Y.std(1, ddof=1), 1e-6)
    se = s / np.sqrt(n)
    # form P: studentized per point
    xP = _fit(lambda th: (ybar - predict(th)) / se, starts, lo, hi)
    t = (ybar - predict(xP)) / se
    z = stats.norm.ppf(np.clip(stats.t.cdf(t, n - 1), 1e-15, 1 - 1e-15))
    pP = float(stats.chi2.sf(np.sum(z ** 2), dof))
    # form Q: pooled variance
    sp2 = np.mean(s ** 2)
    xQ = _fit(lambda th: (ybar - predict(th)), starts, lo, hi)
    F = (np.sum((ybar - predict(xQ)) ** 2) / dof) / (sp2 / n)
    pQ = float(stats.f.sf(F, dof, m * (n - 1)))
    return pP, pQ


def falsified(pP: float, pQ: float) -> bool:
    return pP < ALPHA and pQ < ALPHA


def starts_from(lo, hi, rng, k=3, centre=None):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    out = [centre] if centre is not None else []
    out += [lo + (hi - lo) * rng.uniform(0.2, 0.8, size=len(lo)) for _ in range(k)]
    return out


# ----------------------------------------------------------------------------- the generic loop

def run_part(truths, rivals, simulate, fit_flags, nsim, ns, cvs, rng, tag=""):
    """For every truth, n and CV: simulate nsim designed datasets, fit EVERY rival under both forms, and record
    P(rival falsified | truth). The truth's own entry is its false-falsification rate (type I under the rule)."""
    out = {}
    for cv in cvs:
        for n in ns:
            for truth in truths:
                flags = np.zeros((nsim, len(rivals)), bool)
                for i in range(nsim):
                    Y, ctx = simulate(truth, n, cv, rng)
                    flags[i] = [fit_flags(r, Y, ctx, rng) for r in rivals]
                key = f"{tag}cv{cv}|n{n}|truth {truth}"
                ti = rivals.index(truth)
                out[key] = {"P(falsified)": {r: float(flags[:, j].mean()) for j, r in enumerate(rivals)},
                            "P(separated: rival falsified and truth not)": {
                                r: float((flags[:, j] & ~flags[:, ti]).mean()) for j, r in enumerate(rivals) if r != truth}}
    return out


# ----------------------------------------------------------------------------- part B

B_C = np.array([0.0, 0.005, 0.05, 0.2, 0.5])
B_PH = np.array([5.5, 7.0])
B_STEPS = ("F", "D", "E", "U")
PKA2 = 7.0


def h2po4(ph):
    return 1.0 / (1.0 + 10 ** (ph - PKA2))


def b_cells():
    return [(c, ph, s) for ph in B_PH for c in B_C for s in range(4)]


CELLS_B = b_cells()
CB = np.array([x[0] for x in CELLS_B]); PHB = np.array([x[1] for x in CELLS_B]); SB = np.array([x[2] for x in CELLS_B])
IPH = (PHB == 7.0).astype(int)


def b_predict(rival: str, th: np.ndarray) -> np.ndarray:
    """log rate per cell. th[0:8] = log k_u per (pH, step); the rest rival-specific."""
    lku = th[:8][IPH * 4 + SB]
    ceff = CB * h2po4(PHB) if rival != "R1total" else CB
    if rival == "R0":
        f = np.ones_like(CB)
    elif rival in ("R1", "R1total"):
        f = 1 + th[8:12][SB] * ceff
    elif rival == "R4p":                       # no factor on the furfural branch U
        kap = np.r_[th[8:11], 0.0]
        f = 1 + kap[SB] * ceff
    elif rival == "R2c":                       # Michaelis at K = 0.5 M, initial slope kappa
        f = 1 + th[8:12][SB] * 0.5 * ceff / (0.5 + ceff)
    elif rival == "R7":                        # jump a at K7 (mM) plus linear kappa
        a, kap, K7 = th[8:12][SB], th[12:16][SB], th[16] / 1000.0
        f = 1 + a * ceff / (K7 + ceff) + kap * ceff
    else:
        raise KeyError(rival)
    return lku + np.log(f)


B_BOUNDS = {
    "R0": ([-3] * 8, [3] * 8),
    "R1": ([-3] * 8 + [0] * 4, [3] * 8 + [200] * 4),
    "R1total": ([-3] * 8 + [0] * 4, [3] * 8 + [200] * 4),
    "R4p": ([-3] * 8 + [0] * 3, [3] * 8 + [200] * 3),
    "R2c": ([-3] * 8 + [0] * 4, [3] * 8 + [400] * 4),
    "R7": ([-3] * 8 + [0] * 8 + [1.0], [3] * 8 + [10] * 4 + [200] * 4 + [5.0]),
}


def b_truth(rival: str, rng) -> np.ndarray:
    lku = rng.uniform(-1, 1, 8)
    kap = rng.uniform(5, 40, 4)
    if rival in ("R1", "R1total", "R2c"):
        return np.r_[lku, kap]
    if rival == "R4p":
        return np.r_[lku, kap[:3]]
    if rival == "R7":
        return np.r_[lku, rng.uniform(0.5, 3.0, 4), kap, rng.uniform(1, 5)]
    return lku


def b_simulate(rival, th, n, cv, rng):
    mu = b_predict(rival, th)
    return mu[:, None] + cv * rng.standard_normal((len(mu), n))


B_RIVALS = ["R0", "R1", "R1total", "R4p", "R2c", "R7"]
B_TRUTHS = ["R1", "R1total", "R4p", "R2c", "R7"]


B_DESIGNS = {
    "full": np.ones(len(CB), bool),
    "no 0.5 M arm": CB != 0.5,
    "no 5 mM arm and no inert zero": (CB != 0.005) & (CB != 0.0),
    "one pH (7.0)": PHB == 7.0,
    "no furfural-branch readout": SB != 3,
}


def run_B(nsim, ns, cvs, rng, designs=("full",)):
    out = {}
    for dname in designs:
        mask = B_DESIGNS[dname]

        def sim(truth, n, cv, rng, mask=mask):
            return b_simulate(truth, b_truth(truth, rng), n, cv, rng)[mask], mask

        def flags(r, Y, mask, rng):
            lo, hi = map(np.array, B_BOUNDS[r])
            centre = np.clip(np.r_[np.zeros(8), np.full(len(lo) - 8, 10.0)], lo, hi)
            if r == "R7":
                centre[-1] = 2.0
            if dname == "one pH (7.0)":               # the pH 5.5 rate constants are not in the data: pin them
                lo[:4] = hi[:4] - 1e-6
            return falsified(*two_form_p(Y, lambda x: b_predict(r, x)[mask], starts_from(lo, hi, rng, 2, centre), lo, hi))

        out.update(run_part(B_TRUTHS, B_RIVALS, sim, flags, nsim, ns, cvs, rng, tag=f"{dname}|"))
    return out


# ----------------------------------------------------------------------------- part F

ARMS = ("buffer", "oil_air", "oil_N2", "enal_noOil_N2", "PVoil_N2", "enal_inOil_N2")
PHI = 0.10            # oil volume fraction (Yiltirak's 10 %)
K_ENAL = 30.0         # declared enal oil/water partition for the enal-in-oil arm (only shares the enal between phases)


def f_cells(with_products=True, with_K=True, arms=ARMS):
    cells = []
    for a in arms:
        has_oil = a in ("oil_air", "oil_N2", "PVoil_N2", "enal_inOil_N2")
        for sp in ("M", "F"):
            cells.append((a, sp, "water"))
            if has_oil:
                cells.append((a, sp, "oil"))
        if with_products:
            cells.append((a, "M", "thioether"))
            cells.append((a, "M", "oxidised"))
    if with_K:
        cells += [("K", "M", "K"), ("K", "F", "K")]
    return cells


# parameter vector (13): 0 log S_M, 1 log S_F, 2 log ks_M, 3 log ks_F, 4 log10 K_M, 5 log10 K_F, 6 log cap (capture
# rate x in-situ enal under air), 7 rho_c (FFT / MFT capture), 8 log ox (oil oxidation rate x oxidant under air),
# 9 rho_o (FFT / MFT oxidation), 10 log dose ratio (dosed enal / in-situ enal), 11 log PV ratio (spiked / air oxidant),
# 12 delta (share of spiked peroxide that decomposes to enal). Each rival frees only its own subset; the rest are held at
# BASE ("off"). R2 is capture with partition negligible (K held at 10^-3); R1 partition only; R3 partition + capture;
# R6 partition + oil-phase oxidation.
NEG = -12.0
BASE = np.array([0, 0, 0, 0, -3.0, -3.0, NEG, 0, NEG, 0, 0, 0, 0], float)
F_LO = np.array([-1, -1, -0.5, -0.5, 1.0, 1.0, NEG, 0.0, NEG, 0.0, 0.0, 0.0, 0.0])
F_HI = np.array([1, 1, 1.5, 1.5, 2.5, 2.5, 4.0, 0.2, 5.5, 0.3, 1.2, 1.2, 0.3])
F_FREE = {"R1": [0, 1, 2, 3, 4, 5], "R2": [0, 1, 2, 3, 6, 7, 10, 11, 12],
          "R3": [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12], "R6": [0, 1, 2, 3, 4, 5, 8, 9, 11]}
F_CENTRE = np.array([0, 0, 0.7, 0.7, 1.75, 1.75, 1.5, 0.1, 2.5, 0.1, 0.7, 0.7, 0.1])
PHI_A = np.array([0, PHI, PHI, 0, PHI, PHI])


def f_index(cells):
    """Map each cell to a slot of the value vector built in f_values."""
    idx = []
    for c in cells:
        if c[0] == "K":
            idx.append(28 + (0 if c[1] == "M" else 1)); continue
        ai, si = ARMS.index(c[0]), (0 if c[1] == "M" else 1)
        idx.append({"water": ai * 2 + si, "oil": 12 + ai * 2 + si, "thioether": 30 + ai, "oxidised": 36 + ai}[c[2]])
    return np.array(idx)


def f_values(th, loq=0.01):
    SM, SF, ksM, ksF = np.exp(th[0:4])
    KM, KF = 10 ** th[4], 10 ** th[5]
    cap, rc, ox, ro = np.exp(th[6]), th[7], np.exp(th[8]), th[9]
    dose, pv, delta = np.exp(th[10]), np.exp(th[11]), th[12]
    E = np.array([0.0, 1.0, 0.0, dose, delta * pv, dose / (1 + PHI * K_ENAL)])
    Ox = np.array([0.0, 1.0, 0.0, 0.0, pv, 0.0])
    S, ks, K = np.array([SM, SF]), np.array([ksM, ksF]), np.array([KM, KF])
    cf, of = np.array([1.0, rc]), np.array([1.0, ro])
    w = 1.0 / (1.0 + PHI_A[:, None] * K[None, :])
    lam_c = cap * cf[None, :] * E[:, None] * w
    lam_o = ox * of[None, :] * Ox[:, None] * (1 - w)
    lam = ks[None, :] * w + lam_c + lam_o
    em = np.exp(-lam)
    P1 = S[None, :] / lam * (1 - em)
    Pint = S[None, :] / lam * (1 - (1 - em) / lam)
    ref = P1[0, 0]
    v = np.empty(42)
    v[0:12] = np.log((P1 * w).ravel() + loq * ref)
    v[12:24] = np.log((P1 * (1 - w)).ravel() + loq * ref)
    v[24:28] = 0.0
    v[28], v[29] = np.log(KM), np.log(KF)
    v[30:36] = np.log(lam_c[:, 0] * Pint[:, 0] + loq * ref)
    v[36:42] = np.log(lam_o[:, 0] * Pint[:, 0] + loq * ref)
    return v


def f_truth(rival, rng):
    th = BASE.copy()
    th[0:2] = rng.uniform(-0.5, 0.5, 2)
    th[2:4] = np.log(rng.uniform(1, 3, 2))
    th[10] = np.log(rng.uniform(1, 3))
    th[11] = np.log(rng.uniform(1, 3))
    fall = np.log(rng.uniform(3, 15))          # the fresh-oil-air MFT fall the chemistry makes (x3 - x15, Yiltirak's x10)
    if rival in ("R1", "R3", "R6"):
        th[4:6] = rng.uniform(1.0, 2.5, 2)
    if rival in ("R2", "R3"):
        th[6] = fall + th[2]; th[7] = rng.uniform(0, 0.2); th[12] = rng.uniform(0, 0.3)
    if rival == "R6":
        w = 1 / (1 + PHI * 10 ** th[4])
        th[8] = fall + th[2] - np.log(1 - w); th[9] = rng.uniform(0, 0.3)
    return th


def f_simulate(th, idx, n, cv, rng, kcv=0.2):
    mu = f_values(th)[idx]
    sd = np.where((idx == 28) | (idx == 29), kcv, cv)
    return mu[:, None] + sd[:, None] * rng.standard_normal((len(mu), n))


F_RIVALS = ["R1", "R2", "R3", "R6"]


def f_flags(r, Y, idx, rng):
    free = F_FREE[r]
    lo, hi = F_LO[free], F_HI[free]

    def pred(x):
        th = BASE.copy(); th[free] = x
        return f_values(th)[idx]

    return falsified(*two_form_p(Y, pred, starts_from(lo, hi, rng, 2, np.clip(F_CENTRE[free], lo, hi)), lo, hi))


def run_F(nsim, ns, cvs, rng, designs):
    out = {}
    for dname, (wp, wk, arms) in designs.items():
        idx = f_index(f_cells(wp, wk, arms))
        out.update(run_part(F_RIVALS, F_RIVALS, lambda truth, n, cv, rng, idx=idx: (f_simulate(f_truth(truth, rng), idx, n, cv, rng), idx),
                            f_flags, nsim, ns, cvs, rng, tag=f"{dname}|"))
    return out


# ----------------------------------------------------------------------------- part W

W_M = np.array([0.10, 0.25, 0.40, 0.55, 0.70])
X_SUGAR = 0.25        # sugar mass fraction of the dry matter (declared; the design prints it)


def w_predict(rival, th, wv=W_M):
    if rival == "R2":            # whole-sample (wet-basis) mass action: the shipped lane's basis
        g = (1 - wv)
    elif rival == "PW":          # per-water basis, no cap
        g = X_SUGAR * (1 - wv) / wv
    elif rival == "R4":          # dissolved sugar capped at solubility s (g / g water)
        g = np.minimum(X_SUGAR * (1 - wv) / wv, th[1])
    elif rival == "BAS":         # wet basis with a Bassama-like high-moisture slow-down (x0.38 at 70 % wb)
        g = (1 - wv) * np.exp(-th[1] * np.clip(wv - 0.40, 0, None) / 0.30)
    else:
        raise KeyError(rival)
    return th[0] + np.log(g)


W_BOUNDS = {"R2": ([-5], [5]), "PW": ([-5], [5]), "R4": ([-5, 0.05], [5, 3.0]), "BAS": ([-5, 0.0], [5, 3.0])}


def w_truth(rival, rng):
    if rival == "R4":
        return np.array([0.0, rng.uniform(0.15, 1.0)])     # a cap that binds inside the design's range
    if rival == "BAS":
        return np.array([0.0, np.log(2.6)])                 # Bassama's x2.6 at the wet end
    return np.array([0.0])


W_RIVALS = ["R2", "PW", "R4", "BAS"]


W_DESIGNS = {"5 moistures 10 - 70 % wb": W_M, "3 moistures (10, 40, 70 % wb)": np.array([0.10, 0.40, 0.70])}


def run_W(nsim, ns, cvs, rng):
    out = {}
    for dname, wv in W_DESIGNS.items():
        def sim(truth, n, cv, rng, wv=wv):
            mu = w_predict(truth, w_truth(truth, rng), wv)
            return mu[:, None] + cv * rng.standard_normal((len(mu), n)), wv

        def flags(r, Y, wv, rng):
            lo, hi = map(np.array, W_BOUNDS[r])
            return falsified(*two_form_p(Y, lambda x: w_predict(r, x, wv), starts_from(lo, hi, rng, 3), lo, hi))

        out.update(run_part(W_RIVALS, W_RIVALS, sim, flags, nsim, ns, cvs, rng, tag=f"{dname}|"))
    return out


# ----------------------------------------------------------------------------- main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--part", choices=("B", "F", "W", "all"), default="all")
    ap.add_argument("--nsim", type=int, default=100)
    ap.add_argument("--fnsim", type=int, default=30, help="simulations per cell for part F (its fits are the slow ones)")
    ap.add_argument("--wnsim", type=int, default=200, help="simulations per cell for part W (cheap)")
    ap.add_argument("--seed", type=int, default=20261106)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    res: Dict[str, object] = {"generated_by": "scripts/diagnostics/m5_power_check.py", "seed": a.seed, "nsim": {"B": a.nsim, "B variants": max(10, a.nsim * 2 // 3), "F": a.fnsim, "W": a.wnsim},
                              "alpha": ALPHA, "rule": "falsified iff p < alpha under BOTH forms (P studentized per point, Q pooled variance)",
                              "synthetic_only": True, "parts": {}}
    t0 = time.time()
    if a.part in ("B", "all"):
        res["parts"]["B"] = run_B(a.nsim, (2, 3, 6), (0.08, 0.15), rng)
        res["parts"]["B"].update(run_B(max(10, a.nsim * 2 // 3), (3,), (0.08,), rng,
                                       designs=[d for d in B_DESIGNS if d != "full"]))
        print(f"B done {time.time()-t0:.0f} s", flush=True)
    if a.part in ("F", "all"):
        designs = {"full (a)-(e)": (True, True, ARMS),
                   "no products (d)": (False, True, ARMS),
                   "no K (e), no products": (False, False, ARMS),
                   "arm (a) only: buffer, fresh oil air / N2, no products, no K": (False, False, ARMS[:3]),
                   "arms (a) + (b): + enal without oil, no products, no K": (False, False, ARMS[:4])}
        res["parts"]["F"] = run_F(a.fnsim, (2, 3, 6), (0.15, 0.25), rng, designs); print(f"F done {time.time()-t0:.0f} s", flush=True)
    if a.part in ("W", "all"):
        res["parts"]["W"] = run_W(a.wnsim, (2, 3, 4, 6), (0.20, 0.30), rng); print(f"W done {time.time()-t0:.0f} s", flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
