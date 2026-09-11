#!/usr/bin/env python
"""
Build Wave B39 -- THE FED 3-DEOXYGLUCOSONE FIT (2026-09-11).

Pre-registered in results/validation/kinetic_core_b39_prereg.md BEFORE this ran. Five coordinates
(log10 k at 100 C): the existing k_tdg_ddg and k_ddg_hmf, and the three steps B39 added to the
trunk (r_ddg_tdg, r_ddg_dgal, r_dgal_ddg; one declared barrier). Rows: Mittelmaier 2011's six
printed fed maxima and shares plus two times of maximum, and Zhang 2021's within-study
3,4-DDG/3-DG ratios at 95-110 C. The Leitzen 2021 hold-out is never read here.

Run inside docker:  PYTHONPATH=/workspace python scripts/generators/generate_kinetic_core_b39_fit.py [--quick]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import parameters_dicarbonyl as PD  # noqa: E402
from src.kinetic_core import trunk_conditions  # noqa: E402
from src.kinetic_core.engine import TRUNK, core_parameters  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_furanic import FURANIC_PARAMETERS  # noqa: E402

WAVE = "B39"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b39_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b39_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b39_fit_report.md"
CELSIUS = 273.15
SEED = 20260911
KEYS: Tuple[str, ...] = PD.FED_3DEOXY_COORDINATES
MW = {"TDG": 162.14, "DDG": 144.13, "DGAL": 162.14}
FED_UM = 200.0        # Mittelmaier: "~200 uM each"
FED_T_C, FED_PH = 120.0, 5.0
GRID_FED = np.concatenate([np.arange(0.0, 60.0, 0.5), np.arange(60.0, 120.5, 2.0)])
ZHANG_MM, ZHANG_PH, ZHANG_MIN = 300.0, 6.5, 360.0

MITT = "Mittelmaier et al. 2011, Anal. Bioanal. Chem. 399:1689, Results 'Reaction pathways of 3-DGal formation' (mittelmaier2010_extraction.md sec. 3)"
ZHANG = "Zhang et al. 2021, Food Sci. Nutr. 9:290, Table 1 glucose-only slopes at 6 h, within-study mass ratio 3,4-DDG/3-DG (zhang2020_extraction.md sec. 2)"

# id, pot, observable, target, sigma_log, decisive, anchor
ROWS_RAW = [
    ("fedTDG_ddg_max_uM", "TDG", "max_DDG", 26.7, 0.15, True, MITT),
    ("fedTDG_ddg_tmax_min", "TDG", "tmax_DDG", 30.0, 0.25, True, MITT),
    ("fedTDG_dgal_share_60min", "TDG", "share_DGAL_60", 0.26, 0.5, True, MITT),
    ("fedDGAL_ddg_max_uM", "DGAL", "max_DDG", 46.2, 0.15, True, MITT),
    ("fedDGAL_ddg_tmax_min", "DGAL", "tmax_DDG", 20.0, 0.25, True, MITT),
    ("fedDGAL_tdg_share_60min", "DGAL", "share_TDG_60", 0.48, 0.5, True, MITT),
    ("fedDDG_tdg_max_uM", "DDG", "max_TDG", 26.9, 0.15, True, MITT),
    ("fedDDG_dgal_max_uM", "DDG", "max_DGAL", 37.9, 0.15, True, MITT),
    ("zhang_ratio_90C", "ZHANG", 90.0, (0.0005 * 6 + 0.0016) / (0.0028 * 6 + 0.0049), 0.2, False, ZHANG),
    ("zhang_ratio_95C", "ZHANG", 95.0, (0.0006 * 6 + 0.0016) / (0.0027 * 6 + 0.0063), 0.2, True, ZHANG),
    ("zhang_ratio_100C", "ZHANG", 100.0, (0.0010 * 6 + 0.0028) / (0.0041 * 6 + 0.0119), 0.2, True, ZHANG),
    ("zhang_ratio_105C", "ZHANG", 105.0, (0.0018 * 6 + 0.0024) / (0.0056 * 6 + 0.0173), 0.2, True, ZHANG),
    ("zhang_ratio_110C", "ZHANG", 110.0, (0.0015 * 6 + 0.0047) / (0.0054 * 6 + 0.017), 0.2, True, ZHANG),
]
ROWS = [dict(id=i, pot=p, observable=o, target=float(t), sigma_log=s, decisive=d, anchor=a) for i, p, o, t, s, d, a in ROWS_RAW]
DECISIVE = [r for r in ROWS if r["decisive"]]

_shipped_log10 = {
    "log10_k_tdg_ddg_100C": math.log10(FURANIC_PARAMETERS["k_tdg_ddg"].k_ref),
    "log10_k_ddg_hmf_100C": math.log10(FURANIC_PARAMETERS["k_ddg_hmf"].k_ref),
}
START = np.array([_shipped_log10.get(k, _shipped_log10["log10_k_tdg_ddg_100C"]) for k in KEYS], dtype=float)
LOWER, UPPER = START - 2.0, START + 2.0


def _params(x: np.ndarray, ph: float) -> Dict[str, Any]:
    override = {"fed_3deoxy": dict(zip(KEYS, (float(v) for v in x)))}
    params = core_parameters(TRUNK, frozen=override)
    out, _ = trunk_conditions.apply(dict(params), SimpleNamespace(ph=float(ph), water_activity=None))
    return out


def fed_run(x: np.ndarray, fed: str, quick: bool):
    params = _params(x, FED_PH)
    grid = GRID_FED[::2] if quick else GRID_FED
    return grid, integrate(params, FED_T_C + CELSIUS, {fed: FED_UM / 1000.0}, grid, rtol=1e-8, atol=1e-14)


def zhang_ratio(x: np.ndarray, t_c: float) -> float:
    params = _params(x, ZHANG_PH)
    run = integrate(params, t_c + CELSIUS, {"Glc": ZHANG_MM}, [0.0, ZHANG_MIN], rtol=1e-8, atol=1e-14)
    tdg = float(run.series("TDG")[-1]) * MW["TDG"]; ddg = float(run.series("DDG")[-1]) * MW["DDG"]
    return ddg / max(tdg, 1e-30)


def observables(x: np.ndarray, quick: bool) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for fed in ("TDG", "DDG", "DGAL"):
        grid, run = fed_run(x, fed, quick)
        series = {k: run.series(k) * 1000.0 for k in ("TDG", "DDG", "DGAL")}   # uM
        for k in ("TDG", "DDG", "DGAL"):
            if k == fed:
                continue
            i = int(np.argmax(series[k]))
            out[f"{fed}:max_{k}"] = float(series[k][i]); out[f"{fed}:tmax_{k}"] = float(grid[i])
        j = int(np.argmin(np.abs(grid - 60.0)))
        pair = series["TDG"][j] + series["DGAL"][j]
        out[f"{fed}:share_DGAL_60"] = float(series["DGAL"][j] / max(pair, 1e-30))
        out[f"{fed}:share_TDG_60"] = float(series["TDG"][j] / max(pair, 1e-30))
        out[f"{fed}:pool_60_uM"] = float(series["TDG"][j] + series["DDG"][j] + series["DGAL"][j])
    for r in ROWS:
        if r["pot"] == "ZHANG":
            out[f"ZHANG:{r['observable']}"] = zhang_ratio(x, float(r["observable"]))
    return out


def _pred(obs: Mapping[str, float], r: Mapping[str, Any]) -> float:
    key = f"{r['pot']}:{r['observable']}"
    return max(float(obs[key]), 1e-12)


def residuals(x: np.ndarray, quick: bool = True) -> np.ndarray:
    obs = observables(x, quick)
    return np.array([math.log10(_pred(obs, r) / r["target"]) / r["sigma_log"] for r in DECISIVE])


def fit_member(start: int, quick: bool, max_nfev: int) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    if start == 0:
        x0 = START.copy()
    else:
        rng = np.random.default_rng(SEED + start)
        x0 = np.clip(START + rng.normal(0.0, 0.5, size=START.shape), LOWER, UPPER)
    t0 = time.time()
    sol = least_squares(lambda v: residuals(v, quick), x0, bounds=(LOWER, UPPER), method="trf", max_nfev=max_nfev, xtol=1e-9, ftol=1e-9, diff_step=1e-3)
    r = residuals(sol.x, quick)
    return {"start": start, "x0": x0.tolist(), "x": sol.x.tolist(), "cost": float(np.sum(r * r)), "nfev": int(sol.nfev),
            "status": int(sol.status), "seconds": round(time.time() - t0, 1),
            "residuals_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(DECISIVE, r)}}


def laplace(x: np.ndarray, quick: bool) -> Dict[str, Any]:
    h = 1e-3
    jac = np.zeros((len(DECISIVE), len(KEYS)))
    for j in range(len(KEYS)):
        xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
        jac[:, j] = (residuals(xp, quick) - residuals(xm, quick)) / (2 * h)
    r = residuals(x, quick)
    dof = max(len(DECISIVE) - len(KEYS), 1)
    chi2 = float(r @ r) / dof
    fim = jac.T @ jac
    cov = np.linalg.pinv(fim) * max(chi2, 1.0)
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    diag = np.diag(fim)
    cond = np.where(diag > 0, np.sqrt(max(chi2, 1.0) / np.where(diag > 0, diag, 1.0)), np.inf)
    on_bound = [bool(abs(x[i] - LOWER[i]) < 1e-3 or abs(UPPER[i] - x[i]) < 1e-3) for i in range(len(KEYS))]
    sd = np.sqrt(np.clip(np.diag(cov), 1e-300, None)); corr = cov / np.outer(sd, sd)
    pairs = [{"a": KEYS[i], "b": KEYS[j], "corr": float(corr[i, j])} for i in range(len(KEYS)) for j in range(i + 1, len(KEYS)) if abs(corr[i, j]) >= 0.9]
    verdict = {}
    for i, k in enumerate(KEYS):
        hw = 1.96 * float(sig[i])
        verdict[k] = "AT_BOUND" if on_bound[i] else ("PINNED" if hw <= 0.5 else ("WEAK" if hw <= 1.5 else "UNIDENTIFIED"))
    return {"sigma": dict(zip(KEYS, sig.tolist())), "conditional_sigma": dict(zip(KEYS, [float(c) for c in cond])),
            "chi2_reduced": chi2, "dof": dof, "on_bound": dict(zip(KEYS, on_bound)),
            "identified": {k: bool(s < 1.0 and not b) for k, s, b in zip(KEYS, sig, on_bound)},
            "b38_verdict": verdict, "collinear_pairs": pairs, "n_pinned": sum(v == "PINNED" for v in verdict.values())}


def constants_at(x: np.ndarray, t_c: float) -> Dict[str, float]:
    params = core_parameters(TRUNK, frozen={"fed_3deoxy": dict(zip(KEYS, (float(v) for v in x)))})
    out = {}
    for key in ("k_tdg_ddg", "k_ddg_tdg", "k_ddg_dgal", "k_dgal_ddg", "k_ddg_hmf", "k_tdg_fa", "k_tdg_mgo"):
        q = params[key]
        out[key] = float(q.k_ref) * math.exp(-float(q.ea_kj_mol) * 1000.0 / 8.314462618 * (1.0 / (t_c + CELSIUS) - 1.0 / 373.15))
    return out


def build(quick: bool, max_nfev: int, starts: int) -> Dict[str, Any]:
    assert PREREG.exists(), "B39 is pre-registered; write the prereg before running the fit"
    # B39 ran BEFORE the 3-DG exits had a pH term (B40/B41 added it and B41 shipped it, so the module
    # default is now True). Reproducing B39 means running without it; the B40/B41 wrappers set their own.
    if not getattr(build, "_term_set_by_wrapper", False):
        trunk_conditions.THREE_DEOXY_EXIT_PH_TERM = False
    before_x = np.array([_shipped_log10.get(k, -30.0) for k in KEYS])   # the inert 'before'
    before_obs = observables(before_x, quick)
    members = [fit_member(s, quick, max_nfev) for s in range(starts)]
    best = min(members, key=lambda m: m["cost"])
    x = np.array(best["x"])
    after_obs = observables(x, quick)
    lap = laplace(x, quick)
    rows_out = []
    for r in ROWS:
        rows_out.append({**r, "before": _pred(before_obs, r), "after": _pred(after_obs, r),
                         "before_dex": math.log10(_pred(before_obs, r) / r["target"]), "after_dex": math.log10(_pred(after_obs, r) / r["target"])})
    dec = [r for r in rows_out if r["decisive"]]
    maxima_within_0p3 = all(abs(r["after_dex"]) <= 0.3 for r in dec if "max_uM" in r["id"])
    tmax_in_bracket = (20.0 <= after_obs["TDG:tmax_DDG"] <= 60.0) and (10.0 <= after_obs["DGAL:tmax_DDG"] <= 30.0)
    p1 = maxima_within_0p3 and tmax_in_bracket and lap["chi2_reduced"] < 3.0
    d_tdg_ddg = float(x[KEYS.index("log10_k_tdg_ddg_100C")] - _shipped_log10["log10_k_tdg_ddg_100C"])
    d_ddg_hmf = float(x[KEYS.index("log10_k_ddg_hmf_100C")] - _shipped_log10["log10_k_ddg_hmf_100C"])
    payload = {
        "artifact": "kinetic_core_b39_fit_report", "wave": WAVE,
        "provenance": provenance.provenance_block("kinetic_core_b39_fit_report", generated_by="scripts/generators/generate_kinetic_core_b39_fit.py",
                                                  wave=WAVE, inputs=[PREREG, ROOT / "data/lit/extraction_dossiers/mittelmaier2010_extraction.md",
                                                                     ROOT / "data/lit/extraction_dossiers/zhang2020_extraction.md"]),
        "prereg": data_paths.rel(PREREG), "quick": quick,
        "objective": {"form": "weighted least squares on log10(pred/target)/sigma over the decisive rows; fed pots integrated to 120 min on a 0.5-min grid, "
                              "Zhang pots to 360 min", "n_rows": len(DECISIVE), "n_free": len(KEYS), "sigma_log": {r["id"]: r["sigma_log"] for r in ROWS}},
        "declaration": {"fed_pot": {"charge_uM": FED_UM, "T_C": FED_T_C, "pH": FED_PH, "medium": "water (the PD salts are not printed and are not invented)"},
                        "zhang_pot": {"glucose_mM": ZHANG_MM, "pH": ZHANG_PH, "minutes": ZHANG_MIN},
                        "barrier_kj_mol_declared_for_new_steps": PD.EA_FED_3DEOXY_KJ_MOL,
                        "not_freed": ["k_tdg_fa (Martins, a B1 fit row)", "k_tdg_mgo (Kocadagli)", "k_glc_tdg (the glucose entry)"],
                        "holdout_never_read": "Leitzen 2021 (mp_holdout_glucose_only_autoclave_121C_Steinhagen2021)"},
        "rows": rows_out, "members": members, "best_start": best["start"],
        "frozen_parameters": {"fed_3deoxy": dict(zip(KEYS, x.tolist())), "start_log10": dict(zip(KEYS, START.tolist()))},
        "bounds": {k: [float(LOWER[i]), float(UPPER[i])] for i, k in enumerate(KEYS)},
        "laplace": lap,
        "constants_per_min": {"before_120C": constants_at(before_x, 120.0), "after_120C": constants_at(x, 120.0), "after_100C": constants_at(x, 100.0)},
        "fed_pot_before_after": {"before": {k: v for k, v in before_obs.items()}, "after": {k: v for k, v in after_obs.items()}},
        "predictions": {
            "P1_rows_fit": {"maxima_within_0p3_dex": maxima_within_0p3, "tmax_in_brackets": tmax_in_bracket, "chi2_reduced": lap["chi2_reduced"], "held": bool(p1)},
            "P2_k_tdg_ddg_up_0p5_to_1p2_dex": {"delta_dex": d_tdg_ddg, "held": bool(0.5 <= d_tdg_ddg <= 1.2)},
            "P3_k_ddg_hmf_down_at_least_0p5_dex": {"delta_dex": d_ddg_hmf, "held": bool(d_ddg_hmf <= -0.5)},
            "P5_three_pinned_reverse_pair_collinear": {"n_pinned": lap["n_pinned"], "collinear_pairs": lap["collinear_pairs"],
                                                       "held": bool(lap["n_pinned"] >= 3)},
        },
    }
    return payload


def render(p: Dict[str, Any]) -> str:
    out = [f"# Wave {p['wave']} — the fed 3-deoxyglucosone fit", "", f"_Pre-registered in `{p['prereg']}`. Generated by `{p['provenance']['generated_by']}`._", "",
           f"Objective: {p['objective']['form']}. {p['objective']['n_rows']} decisive rows, {p['objective']['n_free']} free coordinates; "
           f"best start {p['best_start']}, cost {min(m['cost'] for m in p['members']):.2f}, χ²_red {p['laplace']['chi2_reduced']:.2f}.", "",
           "## Rows", "", "| row | target | before | after | before dex | after dex | decisive |", "|---|---:|---:|---:|---:|---:|---|"]
    for r in p["rows"]:
        out.append(f"| {r['id']} | {r['target']:.4g} | {r['before']:.4g} | {r['after']:.4g} | {r['before_dex']:+.2f} | {r['after_dex']:+.2f} | {'yes' if r['decisive'] else 'reported'} |")
    out += ["", "## Coordinates (log10 k at 100 °C)", "", "| coordinate | start | fit | Δ dex | σ (Laplace) | conditional σ | on bound | B38 verdict |", "|---|---:|---:|---:|---:|---:|---|---|"]
    fp = p["frozen_parameters"]; lap = p["laplace"]
    for k in KEYS:
        out.append(f"| {k} | {fp['start_log10'][k]:.3f} | {fp['fed_3deoxy'][k]:.3f} | {fp['fed_3deoxy'][k]-fp['start_log10'][k]:+.2f} | {lap['sigma'][k]:.3f} | {lap['conditional_sigma'][k]:.3f} | {'yes' if lap['on_bound'][k] else ''} | {lap['b38_verdict'][k]} |")
    if lap["collinear_pairs"]:
        out += ["", "Collinear pairs (|corr| ≥ 0.9): " + "; ".join(f"{q['a']} ~ {q['b']} ({q['corr']:+.2f})" for q in lap["collinear_pairs"])]
    out += ["", "## Constants at 120 °C (/min)", "", "| constant | before | after |", "|---|---:|---:|"]
    for k in p["constants_per_min"]["after_120C"]:
        out.append(f"| {k} | {p['constants_per_min']['before_120C'][k]:.3e} | {p['constants_per_min']['after_120C'][k]:.3e} |")
    out += ["", "## Predictions", ""]
    for k, v in p["predictions"].items():
        out.append(f"- **{k}**: {'HELD' if v['held'] else 'REFUTED'} — " + json.dumps({kk: vv for kk, vv in v.items() if kk != 'held'}, default=str)[:240])
    out += ["", "_P4 (the Leitzen hold-out) and P6 (the panel) are judged by `generate_kinetic_core_b39_ship_rule.py` on frozen before/after panels._", ""]
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--max-nfev", type=int, default=250)
    ap.add_argument("--starts", type=int, default=2)
    a = ap.parse_args(argv)
    t0 = time.time()
    p = build(a.quick, a.max_nfev, a.starts)
    p["wall_seconds"] = round(time.time() - t0, 1)
    OUT_JSON.write_text(json.dumps(p, indent=2, default=str) + "\n"); OUT_MD.write_text(render(p))
    print(f"wrote {data_paths.rel(OUT_JSON)}: cost {min(m['cost'] for m in p['members']):.2f} | chi2_red {p['laplace']['chi2_reduced']:.2f} | "
          f"x = {json.dumps({k: round(v, 3) for k, v in p['frozen_parameters']['fed_3deoxy'].items()})} | {p['wall_seconds']} s")
    for k, v in p["predictions"].items():
        print(f"  {k}: {'HELD' if v['held'] else 'REFUTED'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
