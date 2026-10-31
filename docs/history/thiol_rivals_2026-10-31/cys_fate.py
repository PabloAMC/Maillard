"""Free-cysteine fate in Zhai 2021's fed-TTCA pots (100 C, pH 5.5 / 6 / 7 / 8, +- 10 mM xylose): rival structures
B0 / B1 / B2 / B3 / B4 and the negative control NC-b (prereg.yaml), scored by feasibility sampling inside declared
bands -- no fitting. Standalone (A2_MODULE_CONTRACT sec. 3 style): CONFIG_KEYS with bands, rate_constants(),
OBSERVABLES; TTCA and the dicarbonyls are FORCED inputs (piecewise-linear interpolants of the measured series).

Run: ./scripts/docker_maillard.sh run "PYTHONPATH=$PWD python <scratch>/cys_fate.py"
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parent
ROWS = Path("data/lit/fit_rows/route_submodel.yml")
PH = (5.5, 6.0, 7.0, 8.0)
T_GRID = np.array([0.0, 20, 40, 60, 80, 100, 120])

# ---- forced inputs, digitised Zhai 2021 Fig. 4b (TTCA) and Fig. 5 (dicarbonyl sum, umol/L -> mM) -----------------
TTCA = {  # (pH, xyl): series at T_GRID, mmol/L; "0.00" plotted = < LOD, set to 0.005
    (5.5, 0): [10, 9.81, 9.34, 9.10, 8.19, 7.68, 6.69], (6.0, 0): [10, 9.78, 9.18, 8.79, 7.97, 5.98, 4.37],
    (7.0, 0): [10, 9.14, 8.21, 7.02, 5.52, 2.08, 0.89], (8.0, 0): [10, 7.79, 5.18, 3.79, 1.47, 0.005, 0.005],
    (5.5, 1): [10, 9.79, 9.19, 8.79, 7.99, 7.34, 6.59], (6.0, 1): [10, 9.68, 8.97, 8.36, 7.29, 5.12, 3.89],
    (7.0, 1): [10, 8.57, 7.93, 6.97, 4.29, 0.005, 0.005], (8.0, 1): [10, 8.12, 6.07, 2.16, 0.005, 0.005, 0.005]}
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
DC = {k: (np.array(_DX3[k]) + np.array(_DX1[k]) + np.array(_MGO[k])) / 1000.0 for k in TTCA}   # GO < 5 uM: omitted

# xylose's own loss, first order from the no-TTCA control (Fig. 6b): -ln(X120 / 10) / 120
K_XYL = {5.5: -math.log(9.531 / 10) / 120, 6.0: -math.log(9.401 / 10) / 120, 7.0: -math.log(9.401 / 10) / 120,
         8.0: -math.log(9.359 / 10) / 120}

# ---- CONFIG_KEYS: (name, band, basis); log-uniform where band spans > x3 --------------------------------------------
# common to every rival
COMMON = {
    "k_cd_55": ((0.0035, 0.014), "CAMA decay at pH 5.5, 100 C: de Roos 2005 table (engine 0.007) x[0.5, 2]"),
    "cd_ph_slope": ((0.31, 0.62), "log10 per pH of CAMA decay above pH 5.5 (engine extends de Roos' 5.5 -> 6.0 slope 0.62; half as the lower edge)"),
    "release": ((0.5, 1.0), "fraction of CAMA decay that frees Cys (enolisation frees the amine; Strecker / adduct losses take the rest)"),
    "w_cama": ((0.5, 1.0), "share of TTCA's non-releasing loss that goes through CAMA (engine 0.92; rest = direct breakdown)"),
    "k_th": ((5.5e-4, 2.2e-3), "free-Cys thermal loss /min at 100 C (engine k_cys_thermal 1.1e-3, Kang 2026; x[0.5, 2])"),
}
RIVALS = {
    "B0": {"k_r": ((0.0143, 0.0572), "engine k_ttca_cys at 100 C 0.0286 x[0.5, 2]"),
           "K": ((0.045, 0.075), "Cao 2017 K (L/mmol), pH- and T-flat")},
    "B1": {"k_f": ((1e-6, 1.7e-4), "L/mmol/min; <= 0.17 L/mol/min from the xylose contrast (cysteine review)"),
           "K": ((0.045, 0.075), "Cao 2017")},
    "B2": {"k_r": ((0.03, 0.3), "fast exchange /min"), "K45": ((0.045, 0.075), "Cao's K at pH ~4.5"),
           "kallen_slope": ((0.55, 0.83), "Kallen 1971 dlog K_app / dpH")},
    "B3": {"k_f": ((1e-6, 1.7e-4), "as B1"), "K": ((0.045, 0.075), "Cao 2017"),
           "k_s": ((1e-4, 1e-1), "Cys + dicarbonyl, L/mmol/min: DECLARED (no measured constant on disk)")},
    "B4": {"k_f": ((1e-6, 1.7e-4), "as B1"), "K": ((0.045, 0.075), "Cao 2017"),
           "k_ox": ((1e-4, 5e-3), "/min, DECLARED; floor ~ Zheng & Ho 1994 thermal H2S"),
           "q": ((0.0, 1.0), "/mM, xylose-reductone protection, DECLARED")},
    "NCb": {"k_r": ((3.0, 3.0), "instantaneous"), "K": ((0.06, 0.06), "Cao 2017")},
    # POST HOC (written after Tb3 was read; exploratory, never admitted from this run): B1 with a base-accelerated
    # capture k_f(pH) = k_f55 10^(s_f (pH - 5.5)) (Kallen 1971: the cysteine monoanion condenses ~217x faster), K
    # pH-flat, so Cys recycled through TTCA -> CAMA -> Cys consumes xylose catalytically at pH 7-8.
    "B5posthoc": {"k_f": ((1e-6, 1.7e-4), "at pH 5.5, as B1"), "K": ((0.045, 0.075), "Cao 2017"),
                  "s_f": ((0.5, 1.2), "dlog k_f / dpH, Kallen-type monoanion reactivity (DECLARED post hoc)")},
}


def draw(rival: str, rng) -> dict:
    out = {}
    for name, (band, _) in {**COMMON, **RIVALS[rival]}.items():
        lo, hi = band
        if lo == hi:
            out[name] = lo
        elif lo > 0 and hi / lo > 3:
            out[name] = 10 ** rng.uniform(math.log10(lo), math.log10(hi))
        else:
            out[name] = rng.uniform(lo, hi)
    return out


def rate_constants(p: dict, ph: float) -> dict:
    """k at 100 C and the pot's label pH (unbuffered water: label pH only, drift not modelled -- declared)."""
    if "k_r" in p:
        k_r = p["k_r"]
        K = p["K45"] * 10 ** (p["kallen_slope"] * (ph - 4.5)) if "K45" in p else p["K"]
        k_f = K * k_r
    else:
        k_f, K = p["k_f"] * 10 ** (p.get("s_f", 0.0) * (ph - 5.5)), p["K"]
        k_r = k_f / K
    k_cd = p["k_cd_55"] * 10 ** (p["cd_ph_slope"] * (ph - 5.5))
    return {"k_r": k_r, "k_f": k_f, "k_cd": k_cd, "k_x": K_XYL[ph], "k_th": p["k_th"],
            "k_s": p.get("k_s", 0.0), "k_ox": p.get("k_ox", 0.0), "q": p.get("q", 0.0),
            "release": p["release"], "w_cama": p["w_cama"]}


def simulate(p: dict, ph: float, xyl: int) -> dict:
    k = rate_constants(p, ph)
    tt = np.array(TTCA[(ph, xyl)], float)
    slopes = -np.diff(tt) / 20.0                      # measured net TTCA loss rate per 20-min segment
    dc = DC[(ph, xyl)]

    def L(t):
        return slopes[min(int(t // 20.0), 5)]

    def rhs(t, y):
        cys, x, cama, viol = np.maximum(y, 0.0)[0], max(y[1], 0.0), max(y[2], 0.0), y[3]
        T = np.interp(t, T_GRID, tt)
        r_open = k["k_r"] * T
        r_close = k["k_f"] * cys * x
        other = r_close - r_open + L(t)               # TTCA's non-releasing loss implied by the forced series
        v = 0.0
        if other < 0.0:                               # the structure frees more than the measured TTCA loss allows:
            v = -other                                # cap the opening (record the violation) so mass stays honest
            r_open += other
            other = 0.0
        r_cd = k["k_cd"] * cama
        sink = (k["k_th"] + k["k_s"] * np.interp(t, T_GRID, dc) + k["k_ox"] / (1.0 + k["q"] * x)) * cys
        return [r_open - r_close + k["release"] * r_cd - sink,
                r_open - r_close - k["k_x"] * x,
                k["w_cama"] * other - r_cd,
                v]

    sol = solve_ivp(rhs, (0, 120), [0.0, 10.0 * xyl, 0.0, 0.0], t_eval=T_GRID, method="LSODA", rtol=1e-6, atol=1e-9,
                    max_step=2.0)
    return {"Cys": np.maximum(sol.y[0], 1e-6), "Xyl": np.maximum(sol.y[1], 1e-6), "viol": sol.y[3][-1]}


# ---- observations (keyed rows, route_submodel.yml) --------------------------------------------------------------
def load_obs():
    rows = yaml.safe_load(ROWS.read_text())["rows"]
    cys, xyl = {}, {}
    for r in rows:
        pot = r["pot"]
        if not pot.startswith("zhai2021_") or not pot.endswith("_fig6"):
            continue
        tag = pot.split("_")[1]
        ph = {"ph55": 5.5, "ph6": 6.0, "ph7": 7.0, "ph8": 8.0}[pot.split("_")[2]]
        if r["species"] == "Cys":
            cys.setdefault((ph, 1 if tag == "ttcaxyl" else 0), {})[r["t_min"]] = r["value"]
        elif r["species"] == "PENT":
            xyl.setdefault((ph, tag), {})[r["t_min"]] = r["value"]
    return cys, xyl


def observables(cys_by_pot, xyl_ttca_pot=None, xyl_ctrl=None):
    """Tb1 (per pH): mean log10 Cys(+X)/Cys(-X) over 60-120 min; Tb2 (per pH): log10 Cys(120) in the TTCA pot;
    Tb3 (per pH): log10 Xyl(+TTCA, 120) / Xyl(control, 120)."""
    tb1, tb2, tb3 = [], [], []
    for ph in PH:
        a, b = cys_by_pot[(ph, 1)], cys_by_pot[(ph, 0)]
        tb1.append(np.mean([math.log10(a[t] / b[t]) for t in (60, 80, 100, 120)]))
        tb2.append(math.log10(b[120]))
        if xyl_ttca_pot is not None:
            tb3.append(math.log10(xyl_ttca_pot[ph] / xyl_ctrl[ph]))
    return np.array(tb1), np.array(tb2), np.array(tb3)


def predict(p):
    cys, xt, viol = {}, {}, 0.0
    for ph in PH:
        for xyl in (0, 1):
            s = simulate(p, ph, xyl)
            cys[(ph, xyl)] = dict(zip(T_GRID.astype(int).tolist(), s["Cys"].tolist()))
            viol = max(viol, s["viol"] / max(10.0 - TTCA[(ph, xyl)][-1], 0.5))
            if xyl:
                xt[ph] = s["Xyl"][-1]
    ctrl = {ph: 10.0 * math.exp(-K_XYL[ph] * 120) for ph in PH}
    tb1, tb2, tb3 = observables(cys, xt, ctrl)
    return tb1, tb2, tb3, viol


# ---- likelihoods (noise from prereg.yaml) -----------------------------------------------------------------------
def loglik(tb, obs, kind, noise):
    if kind == "Tb1":
        s = noise["contrast_cys_dex"]
        return -0.5 * np.sum(((tb - obs) / s) ** 2, axis=-1) - len(PH) * math.log(s)
    if kind == "Tb3":
        s = noise["contrast_xyl_dex"]
        return -0.5 * np.sum(((tb - obs) / s) ** 2, axis=-1) - len(PH) * math.log(s)
    # Tb2: shared figure offset ~N(0, tau) + independent point noise -> multivariate normal
    tau, e = noise["level_figure_offset_dex"], math.hypot(noise["level_point_dex"], 0.05)
    n = len(PH)
    cov = np.full((n, n), tau ** 2) + np.eye(n) * e ** 2
    ci = np.linalg.inv(cov)
    d = tb - obs
    return -0.5 * np.einsum("...i,ij,...j->...", d, ci, d) - 0.5 * math.log(np.linalg.det(cov))


def log_marginal(preds, obs, kind, noise):
    ll = loglik(preds, obs, kind, noise)
    m = ll.max()
    return m + math.log(np.mean(np.exp(ll - m)))


def main():
    prereg = yaml.safe_load((OUT / "prereg.yaml").read_text())
    noise = prereg["noise"]
    rng = np.random.default_rng(20261031)
    n_score, n_truth = int(sys.argv[1]) if len(sys.argv) > 1 else 300, 60
    cys_obs, xyl_obs = load_obs()
    xt = {ph: xyl_obs[(ph, "ttcaxyl")][120] for ph in PH}
    xc = {ph: xyl_obs[(ph, "xyl")][120] for ph in PH}
    obs = dict(zip(("Tb1", "Tb2", "Tb3"), observables(cys_obs, xt, xc)))
    rivals = list(RIVALS)
    pred = {}
    for rv in rivals:
        rows = [predict(draw(rv, rng)) for _ in range(n_score + n_truth)]
        pred[rv] = {k: np.array([r[i] for r in rows]) for i, k in enumerate(("Tb1", "Tb2", "Tb3", "viol"))}
    PRE = [r for r in rivals if r != "B5posthoc"]
    # power: synthetic data from each rival's TRUTH draws (held apart from the scoring draws), realistic noise
    power = {}
    for test in ("Tb1", "Tb2", "Tb3"):
        power[test] = {}
        for truth in PRE:
            wins = []
            for j in range(n_truth):
                y = pred[truth][test][n_score + j].copy()
                if test == "Tb2":
                    y = y + rng.normal(0, noise["level_figure_offset_dex"]) + rng.normal(0, math.hypot(noise["level_point_dex"], 0.05), len(PH))
                else:
                    y = y + rng.normal(0, noise["contrast_cys_dex" if test == "Tb1" else "contrast_xyl_dex"], len(PH))
                lm = {rv: log_marginal(pred[rv][test][:n_score], y, test, noise) for rv in PRE}
                wins.append(max(lm, key=lm.get) == truth)
            power[test][truth] = float(np.mean(wins))
    # pairwise power (truth vs each alternative) -- the decision is usually binary
    pair = {}
    for test in ("Tb1", "Tb2", "Tb3"):
        pair[test] = {}
        for a in PRE:
            for b in PRE:
                if a >= b:
                    continue
                ok = []
                for truth, other in ((a, b), (b, a)):
                    for j in range(n_truth):
                        y = pred[truth][test][n_score + j].copy()
                        if test == "Tb2":
                            y = y + rng.normal(0, noise["level_figure_offset_dex"]) + rng.normal(0, math.hypot(noise["level_point_dex"], 0.05), len(PH))
                        else:
                            y = y + rng.normal(0, noise["contrast_cys_dex" if test == "Tb1" else "contrast_xyl_dex"], len(PH))
                        ok.append(log_marginal(pred[truth][test][:n_score], y, test, noise)
                                  > log_marginal(pred[other][test][:n_score], y, test, noise))
                pair[test][f"{a}|{b}"] = float(np.mean(ok))
    score = {}
    for test in ("Tb1", "Tb2", "Tb3"):
        score[test] = {rv: log_marginal(pred[rv][test][:n_score], obs[test], test, noise) for rv in rivals}
    summary = {rv: {t: {"median": np.median(pred[rv][t][:n_score], axis=0).round(3).tolist(),
                        "p05": np.percentile(pred[rv][t][:n_score], 5, axis=0).round(3).tolist(),
                        "p95": np.percentile(pred[rv][t][:n_score], 95, axis=0).round(3).tolist()}
                    for t in ("Tb1", "Tb2", "Tb3")} | {"ttca_gate_violation_median": float(np.median(pred[rv]["viol"][:n_score])),
                                                     "ttca_gate_pass_share": float(np.mean(pred[rv]["viol"][:n_score] < 0.1))}
               for rv in rivals}
    # feasibility: share of draws within 2 sigma at every pH (contrasts) / after the best figure offset (levels)
    feas = {}
    for rv in rivals:
        f = {}
        f["Tb1"] = float(np.mean(np.all(np.abs(pred[rv]["Tb1"][:n_score] - obs["Tb1"]) <= 2 * noise["contrast_cys_dex"], axis=1)))
        d2 = pred[rv]["Tb2"][:n_score] - obs["Tb2"]
        off = np.clip(d2.mean(axis=1, keepdims=True), -2 * noise["level_figure_offset_dex"], 2 * noise["level_figure_offset_dex"])
        f["Tb2"] = float(np.mean(np.all(np.abs(d2 - off) <= 2 * math.hypot(noise["level_point_dex"], 0.05), axis=1)))
        f["Tb3"] = float(np.mean(np.all(np.abs(pred[rv]["Tb3"][:n_score] - obs["Tb3"]) <= 2 * noise["contrast_xyl_dex"], axis=1)))
        ok12 = (np.all(np.abs(pred[rv]["Tb1"][:n_score] - obs["Tb1"]) <= 2 * noise["contrast_cys_dex"], axis=1)
                & np.all(np.abs(d2 - off) <= 2 * math.hypot(noise["level_point_dex"], 0.05), axis=1))
        f["Tb1_and_Tb2"] = float(np.mean(ok12))
        f["Tb1_and_Tb2_and_gate"] = float(np.mean(ok12 & (pred[rv]["viol"][:n_score] < 0.1)))
        feas[rv] = f
    # JOINT score (one constant set per draw predicts all three tests): log mean_draws prod_tests L
    def joint(rv, tests):
        ll = sum(loglik(pred[rv][t][:n_score], obs[t], t, noise) for t in tests)
        m = ll.max()
        return float(m + math.log(np.mean(np.exp(ll - m))))
    joint_scores = {"+".join(ts): {rv: joint(rv, ts) for rv in rivals}
                    for ts in (("Tb1", "Tb2"), ("Tb1", "Tb2", "Tb3"))}
    pre = [r for r in rivals if r != "B5posthoc"]
    family_power = {}
    for test in ("Tb1", "Tb2", "Tb3"):
        fam = {"fast": ("B0", "B2", "NCb"), "slow": ("B1", "B3", "B4")}
        ok = []
        for fname, members in fam.items():
            for truth in members:
                for j in range(n_truth):
                    y = pred[truth][test][n_score + j].copy()
                    if test == "Tb2":
                        y = y + rng.normal(0, noise["level_figure_offset_dex"]) + rng.normal(0, math.hypot(noise["level_point_dex"], 0.05), len(PH))
                    else:
                        y = y + rng.normal(0, noise["contrast_cys_dex" if test == "Tb1" else "contrast_xyl_dex"], len(PH))
                    lm = {rv: log_marginal(pred[rv][test][:n_score], y, test, noise) for rv in pre}
                    ok.append(max(lm, key=lm.get) in members)
        family_power[test] = float(np.mean(ok))
    out = {"observed": {k: v.round(3).tolist() for k, v in obs.items()}, "power_all_rivals": power,
           "family_power_fast_vs_slow": family_power, "joint_log_marginal": joint_scores,
           "power_pairwise": pair, "log_marginal": score, "feasibility": feas, "prediction_summary": summary,
           "n_score_draws": n_score, "n_truth_draws": n_truth, "ph": PH}
    (OUT / "cys_fate_results.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("observed", "power_all_rivals", "family_power_fast_vs_slow", "log_marginal", "joint_log_marginal", "feasibility")}, indent=1))
    for rv in rivals:
        print(rv, {t: summary[rv][t]["median"] for t in ("Tb1", "Tb2", "Tb3")}, "gate", summary[rv]["ttca_gate_pass_share"])


if __name__ == "__main__":
    main()
