"""MFT loss law: rival structures A0 / A1 / A2 / A3 and the negative control NC-a (prereg.yaml), scored by feasibility
sampling inside declared bands -- no fitting. The formation supply is a NUISANCE drawn from a declared family per pot
so that no loss law can win by a tuned supply. Standalone (A2 contract sec. 3 style).

M' = s(t) - lambda(t, pot) M ; exponential integrator on a 0.5-min grid, vectorised over draws.
Run: ./scripts/docker_maillard.sh run "PYTHONPATH=$PWD python <scratch>/mft_loss.py [n_draws]"
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import norm

OUT = Path(__file__).resolve().parent
R = 8.314e-3
RIVALS = ("A0", "A1", "A2", "A3", "NCa")


def arr(ea, t_c, tref_c):
    return np.exp(-ea / R * (1.0 / (t_c + 273.15) - 1.0 / (tref_c + 273.15)))


def logu(rng, lo, hi, n):
    return 10 ** rng.uniform(math.log10(lo), math.log10(hi), n)


# ---- CONFIG_KEYS (bands and basis) --------------------------------------------------------------------------------
BANDS = {
    "A0": "engine k_mft_decay at R3En: 1.384 /min at 145 C (k_ref 1.362 x 10^0.007), Ea 101.25 kJ/mol, fit SD 0.172 dex",
    "A1": "k50 [5e-4, 9e-4] /min (van Seeventer entries 1 / 4, Cys-containing matrices, 59-69 %/day as first order); Ea [40, 133]",
    "A2": "lambda_E(100 C) [1e-3, 1e-1] /min DECLARED; E builds with tau_E(100 C) [60, 1000] min; Ea [40, 133]; E = 0 in clean buffer, full in a pre-reacted flavour",
    "A3": "lambda0(50 C) [1.6e-3, 1e-2] (clean buffer > 90 %/day); K_c from entry 4 (lambda4 8.13e-4 x10^+-0.1); Ea [40, 133]; lambda = lambda0 / (1 + Cys/K_c)",
    "NCa": "lambda = 0",
}


def draw(rival, rng, n):
    d = {"rival": rival, "n": n}
    if rival == "A0":
        d["k145"] = 1.384 * 10 ** rng.normal(0, 0.172, n)
        d["ea"] = np.full(n, 101.25)
    elif rival == "A1":
        d["k50"] = rng.uniform(5e-4, 9e-4, n)
        d["ea"] = rng.uniform(40, 133, n)
    elif rival == "A2":
        d["lamE100"] = logu(rng, 1e-3, 1e-1, n)
        d["tauE100"] = logu(rng, 60, 1000, n)
        d["ea"] = rng.uniform(40, 133, n)
    elif rival == "A3":
        d["lam0"] = logu(rng, 1.6e-3, 1e-2, n)
        lam4 = 8.13e-4 * 10 ** rng.uniform(-0.1, 0.1, n)
        d["Kc"] = 150.0 / np.maximum(d["lam0"] / lam4 - 1.0, 1e-3)
        d["ea"] = rng.uniform(40, 133, n)
    return d


def lam(d, t_c, t, cys=None, matrix="pot"):
    """lambda(t) /min, shape (n, len(t)). cys: (n, len(t)) free cysteine mM for A3. matrix: 'pot' (heated Maillard
    pot, E builds from zero), 'flavour' (pre-reacted: E full), 'buffer' (no electrophile)."""
    rv = d["rival"]
    n = d["n"]
    tt = np.broadcast_to(t, (n, len(t)))
    if rv == "NCa":
        return np.zeros((n, len(t)))
    if rv == "A0":
        return (d["k145"] * arr(d["ea"], t_c, 145.0))[:, None] * np.ones_like(tt)
    if rv == "A1":
        return (d["k50"] * arr(d["ea"], t_c, 50.0))[:, None] * np.ones_like(tt)
    if rv == "A2":
        a = arr(d["ea"], t_c, 100.0)
        if matrix == "buffer":
            return np.zeros((n, len(t)))
        if matrix == "flavour":
            return (d["lamE100"] * a)[:, None] * np.ones_like(tt)
        tau = (d["tauE100"] / a)[:, None]
        return (d["lamE100"] * a)[:, None] * (1.0 - np.exp(-tt / tau))
    if rv == "A3":
        return (d["lam0"] * arr(d["ea"], t_c, 50.0))[:, None] / (1.0 + cys / d["Kc"][:, None])


def cys_series(rng, n, t, c0, tau):
    c = rng.uniform(*c0, n) if c0[0] != c0[1] else np.full(n, c0[0])
    ta = logu(rng, *tau, n)
    return c[:, None] * np.exp(-t[None, :] / ta[:, None])


def supply(rng, n, t, tau, p, tex):
    tl = rng.uniform(*tau, n)[:, None]
    pw = rng.uniform(*p, n)[:, None]
    te = logu(rng, *tex, n)[:, None]
    x = np.maximum(t[None, :] - tl, 0.0)
    return np.where(x > 0, (x + 1e-9) ** pw * np.exp(-x / te), 0.0)


def integrate(s, lm, dt):
    """M_{i+1} = M_i e^{-lam dt} + s dt (1 - e^{-lam dt}) / (lam dt); returns M (n, nt)."""
    m = np.zeros_like(s)
    e = np.exp(-lm * dt)
    g = np.where(lm * dt > 1e-9, (1.0 - e) / np.maximum(lm * dt, 1e-300), 1.0)
    for i in range(s.shape[1] - 1):
        m[:, i + 1] = m[:, i] * e[:, i] + 0.5 * (s[:, i] + s[:, i + 1]) * dt * g[:, i]
    return m


def at(t, m, times):
    return np.stack([m[:, int(round(x / (t[1] - t[0])))] for x in times], axis=1)


# ---- the tests ------------------------------------------------------------------------------------------------------
def ta1_spike(d, rng, n):
    """US 10,076,126: [2H3]-MFT charged at t = 0, 4 h at 110 C, Xyl + Cys 0.25 mol/kg; quantified => survival >= thr,
    thr log-uniform [1e-3, 0.1]. Returns P(pass | draw)."""
    t = np.arange(0, 240.5, 0.5)
    cys = cys_series(rng, n, t, (30.0, 250.0), (100.0, 2000.0))
    lm = lam(d, 110.0, t, cys)
    log10_s = -np.trapezoid(lm, t, axis=1) / math.log(10)
    return np.clip((log10_s + 3.0) / 2.0, 0.0, 1.0)


def ta2_vanseeventer(d, rng, n):
    """Entries 2 (buffer), 3 (+ ribose), 4 (+150 mM Cys), 5 (+Cys + ribose; free Cys a nuisance [20, 150]), 50 C,
    1 day: fraction lost."""
    t = np.array([0.0, 1440.0])
    out = []
    for c in (0.0, 0.0, 150.0, None):
        cys = (np.full((n, 2), c) if c is not None else np.repeat(rng.uniform(20, 150, n)[:, None], 2, axis=1))
        lm = lam(d, 50.0, t, cys, matrix="buffer")[:, 0]
        out.append(1.0 - np.exp(-lm * 1440.0))
    return np.stack(out, axis=1)


def ta3_schieberle(d, rng, n):
    t = np.arange(0, 720.5, 0.5)
    s = supply(rng, n, t, (5, 20), (0.0, 1.5), (120, 1800))
    cys = cys_series(rng, n, t, (3.0, 33.0), (100.0, 2000.0))
    m = integrate(s, lam(d, 100.0, t, cys), 0.5)
    v = np.maximum(at(t, m, (30, 60, 360, 720)), 1e-300)
    return np.log10(np.stack([v[:, 1] / v[:, 0], v[:, 2] / v[:, 1], v[:, 3] / v[:, 2]], axis=1))


def ta4_cerny(d, rng, n):
    """Pot D, 95 C, 4 h: NF-route supply ~ t exp(-t / t_N) (H2S builds, fed NF slowly used), ribose route lagged;
    the routes' magnitude ratio rho log-uniform [0.003, 1] (the route question, not the loss law)."""
    t = np.arange(0, 240.5, 0.5)
    sr = supply(rng, n, t, (5, 20), (0.0, 1.5), (120, 1800))
    tn = logu(rng, 60, 2000, n)[:, None]
    sn = t[None, :] * np.exp(-t[None, :] / tn)
    rho = logu(rng, 0.003, 1.0, n)[:, None]
    sr = sr / np.trapezoid(sr, t, axis=1)[:, None]
    sn = rho * sn / np.trapezoid(sn, t, axis=1)[:, None]
    cys = cys_series(rng, n, t, (10.0, 97.0), (100.0, 2000.0))
    lm = lam(d, 95.0, t, cys)
    mn, mr = integrate(sn, lm, 0.5)[:, -1], integrate(sr, lm, 0.5)[:, -1]
    return (mn / np.maximum(mn + mr, 1e-300))[:, None]


def ta5_liu(d, rng, n):
    t = np.arange(0, 60.5, 0.5)
    s = supply(rng, n, t, (2, 8), (0.0, 1.5), (5, 120))
    cys = cys_series(rng, n, t, (80.0, 330.0), (20.0, 400.0))
    m = integrate(s, lam(d, 168.0, t, cys), 0.5)
    v = np.maximum(at(t, m, (20, 40, 60)), 1e-300)
    return np.log10(np.stack([v[:, 1] / v[:, 0], v[:, 2] / v[:, 1]], axis=1))


TESTS = {"Ta1": ta1_spike, "Ta2": ta2_vanseeventer, "Ta3": ta3_schieberle, "Ta4": ta4_cerny, "Ta5": ta5_liu}
OBS = {"Ta1": None, "Ta2": np.array([0.95, 0.95, 0.69, 0.42]),
       "Ta3": np.log10([13.8 / 4.5, 156 / 13.8, 179 / 156]), "Ta4": np.array([0.07]),
       "Ta5": np.log10([760.33 / 1252.67, 628.00 / 760.33])}


def loglik(test, pred, y, noise, real=False):
    """pred (n, k); y the datum. Returns (n,) log likelihoods."""
    if test == "Ta1":            # pred is P(pass); y = 1 (pass) or 0
        p = np.clip(pred[:, 0] if pred.ndim > 1 else pred, 1e-6, 1 - 1e-6)
        return np.log(p) if y else np.log(1 - p)
    if test == "Ta2":
        s = 0.07
        if real:                 # entries 2 and 3 censored "> 90 %"
            return (norm.logcdf((pred[:, 0] - 0.90) / s) + norm.logcdf((pred[:, 1] - 0.90) / s)
                    + norm.logpdf(y[2], pred[:, 2], s) + norm.logpdf(y[3], pred[:, 3], s))
        return np.sum(norm.logpdf(y[None, :], pred, s), axis=1)
    s = noise["cerny_share_abs"] if test == "Ta4" else noise["schieberle_ratio_dex"]
    return np.sum(norm.logpdf(y[None, :], pred, s), axis=1)


def log_marginal(test, pred, y, noise, real=False):
    ll = loglik(test, pred, y, noise, real)
    ll = ll[np.isfinite(ll)]
    if ll.size == 0:
        return -np.inf
    m = ll.max()
    return float(m + math.log(np.mean(np.exp(ll - m))))


def synth(test, pred_row, rng, noise):
    if test == "Ta1":
        return int(rng.random() < pred_row[0])
    s = 0.07 if test == "Ta2" else (noise["cerny_share_abs"] if test == "Ta4" else noise["schieberle_ratio_dex"])
    return pred_row + rng.normal(0, s, pred_row.shape)


def main():
    noise = yaml.safe_load((OUT / "prereg.yaml").read_text())["noise"]
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    n_truth = 200
    rng = np.random.default_rng(20261031)
    pred = {}
    for rv in RIVALS:
        pred[rv] = {}
        d = draw(rv, rng, n + n_truth)          # ONE constant set per draw, shared by every test (joint feasibility)
        for test, fn in TESTS.items():
            p = fn(d, rng, n + n_truth)
            pred[rv][test] = p if p.ndim > 1 else p[:, None]
    res = {"power_all": {}, "power_pairwise": {}, "score_real": {}, "feasibility": {}, "summary": {}}
    for test in TESTS:
        res["power_all"][test], res["power_pairwise"][test] = {}, {}
        lm_cache = {}
        for truth in RIVALS:
            wins = []
            for j in range(n_truth):
                y = synth(test, pred[truth][test][n + j], rng, noise)
                lm = {rv: log_marginal(test, pred[rv][test][:n], y, noise) for rv in RIVALS}
                lm_cache[(truth, j)] = lm
                wins.append(max(lm, key=lm.get) == truth)
            res["power_all"][test][truth] = float(np.mean(wins))
        for i, a in enumerate(RIVALS):
            for b in RIVALS[i + 1:]:
                ok = [lm_cache[(a, j)][a] > lm_cache[(a, j)][b] for j in range(n_truth)]
                ok += [lm_cache[(b, j)][b] > lm_cache[(b, j)][a] for j in range(n_truth)]
                res["power_pairwise"][test][f"{a}|{b}"] = float(np.mean(ok))
        y_real = 1 if test == "Ta1" else OBS[test]
        res["score_real"][test] = {rv: log_marginal(test, pred[rv][test][:n], y_real, noise, real=True) for rv in RIVALS}
        for rv in RIVALS:
            p = pred[rv][test][:n]
            res["summary"].setdefault(rv, {})[test] = {"median": np.median(p, axis=0).round(4).tolist(),
                                                      "p05": np.percentile(p, 5, axis=0).round(4).tolist(),
                                                      "p95": np.percentile(p, 95, axis=0).round(4).tolist()}
            if test == "Ta1":
                f = float(np.mean(p[:, 0] >= 0.5))
            elif test == "Ta2":
                f = float(np.mean((p[:, 0] > 0.9 - 0.14) & (p[:, 1] > 0.9 - 0.14) & (np.abs(p[:, 2] - 0.69) <= 0.14)
                                  & (np.abs(p[:, 3] - 0.42) <= 0.14)))
            else:
                s = noise["cerny_share_abs"] if test == "Ta4" else noise["schieberle_ratio_dex"]
                f = float(np.mean(np.all(np.abs(p - OBS[test]) <= 2 * s, axis=1)))
            res["feasibility"].setdefault(rv, {})[test] = f
    # joint feasibility over a set of tests (same constants; nuisances independent per pot)
    for combo in (("Ta1", "Ta2", "Ta5"), ("Ta1", "Ta5"), ("Ta1", "Ta3", "Ta5"), ("Ta1", "Ta2", "Ta3", "Ta5")):
        for rv in RIVALS:
            ok = np.ones(n, bool)
            for test in combo:
                p = pred[rv][test][:n]
                if test == "Ta1":
                    ok &= p[:, 0] >= 0.5
                elif test == "Ta2":
                    ok &= ((p[:, 0] > 0.76) & (p[:, 1] > 0.76) & (np.abs(p[:, 2] - 0.69) <= 0.14) & (np.abs(p[:, 3] - 0.42) <= 0.14))
                else:
                    sg = noise["cerny_share_abs"] if test == "Ta4" else noise["schieberle_ratio_dex"]
                    ok &= np.all(np.abs(p - OBS[test]) <= 2 * sg, axis=1)
            res["feasibility"][rv]["+".join(combo)] = float(np.mean(ok))
    res["bands"] = BANDS
    res["n_draws"], res["n_truth"] = n, n_truth
    (OUT / "mft_loss_results.json").write_text(json.dumps(res, indent=1))
    for k in ("power_all", "score_real", "feasibility"):
        print(k, json.dumps(res[k], indent=0))
    for rv in RIVALS:
        print(rv, {t: res["summary"][rv][t]["median"] for t in TESTS})


if __name__ == "__main__":
    main()
