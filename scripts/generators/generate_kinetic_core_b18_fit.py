#!/usr/bin/env python
"""
Build Wave B18 -- A PYRAZINE STEP ON THE TRUNK (2026-09-08).

Pre-registered in ``results/validation/kinetic_core_b18_prereg.md`` before any B18 number existed.
Five steps entered the trunk network (network.PYRAZINE_REACTIONS): two Strecker deaminations
(glyoxal or methylglyoxal + glycine -> aminoketone), rate-determining and FITTED here, and three
aminoketone condensations on one constant DECLARED fast (parameters_pyrazine.K_COND). The pyrazine
formation rate is then the Strecker rate over two, the "zero-order" rate Zhou 2024 printed.

WHAT IS FITTED: six coordinates.
  log10 k_go_ak (100 C, pH 6.8) and its barrier; log10 k_mgo_ak and its barrier  -- on Zhou 2024's
  six formation rates (fed 20 mM glyoxal or methylglyoxal + 20 mM alanine, water, initial pH 8,
  100 / 110 / 120 C, 0-120 min; zhou2024_extraction.md Table 2), alanine -> glycine DECLARED;
  the two slopes of the pyrazine pH term (knot 7, reference 6.8) -- on Leahy & Reineccius 1989's
  four within-study pH ratios at 95 C (leahy1989a_extraction.md: pyrazine and methylpyrazine at
  pH 7 and 5 over pH 9, lysine + glucose 100 + 100 mM; lysine -> glycine DECLARED).
Ten rows, six free. The barriers' bands are the printed value to the dossier's three-point refit
(100.59-103.1 and 111.66-114.9 kJ/mol). Two starts (the prior centre; a seeded perturbation),
scipy least_squares (trf) on log10 residuals, then a Laplace covariance at the optimum.

The observable is the MEAN formation rate over the source's regression window (0-120 min for
Zhou; Leahy's 2 h treatment), in umol L-1 min-1, so a pot whose dicarbonyl the trunk's own sinks
consume during the window is modelled as it would be measured; the glyoxal half-life under the
declared dry-glass sink (B13, k_go_sink, Ea 0) is reported beside the fit as a diagnostic.

Usage:
    python scripts/generators/generate_kinetic_core_b18_fit.py            # both starts, Laplace, report
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters, trunk_conditions  # noqa: E402
from src.kinetic_core.engine import b1_fitted  # noqa: E402
from src.kinetic_core.integrate import integrate  # noqa: E402
from src.kinetic_core.parameters_pyrazine import (  # noqa: E402
    FROZEN_B18, K_COND_DECLARED_L_PER_MMOL_MIN, PYRAZINE_FIT_PH, with_fitted_pyrazine,
)

WAVE = "B18"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b18_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b18_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b18_fit_report.md"
CELSIUS = 273.15
SEED = 20260908
WINDOW_MIN = 120.0
GRID = np.linspace(0.0, WINDOW_MIN, 13)

# ===========================================================================
# 1. THE SYSTEMS AND ROWS
# ===========================================================================
ZHOU_RATES_UMOL_L_MIN: Dict[str, Dict[float, float]] = {
    # zhou2024_extraction.md sec. 4, Table 2 (printed): zero-order slopes, 0-120 min
    "PZ": {100.0: 0.0279, 110.0: 0.0791, 120.0: 0.1507},
    "DMP": {100.0: 0.0035, 110.0: 0.0100, 120.0: 0.0230},
}
ZHOU_ANCHOR = ("Zhou et al. 2024 JAFC 72:18630 Table 2: [Ala] = [dicarbonyl] = 20 mmol/L, water, initial pH 8.0 (NaOH, "
               "unbuffered), stirred sealed vessel, 0-120 min, triplicates, HS-SPME-GC/MS with external calibration; "
               "zhou2024_extraction.md sec. 4")
LEAHY_ANCHOR = ("Leahy & Reineccius 1989 ACS Symp. Ser. 409 ch. 18 Table I: 0.1 M lysine + 0.1 M glucose, 95 C, pH 9.0 "
                "(0.1 M borate) / 7.0 and 5.0 (citrate-phosphate), pseudo-zero-order k in ppm/h; leahy1989a_extraction.md "
                "sec. 4 (within-study ratios at 95 C)")
#: Leahy's within-study ratios at 95 C (the dossier's table): k(pH x) / k(pH 9)
LEAHY_RATIOS: Dict[str, Dict[float, float]] = {
    "PZ": {7.0: 1.0 / 2.67, 5.0: 1.0 / 38.3},
    "MPZ": {7.0: 1.0 / 2.08, 5.0: 1.0 / 446.0},
}
SIGMA_ZHOU = 0.10       # the printed regressions have R2 >= 0.99; the vessel's come-up time is unstated
SIGMA_LEAHY = 0.15      # the buffer changes between arms (borate at 9, citrate-phosphate at 7 and 5)
SIGMA_LEAHY_WEAK = 0.30  # the methylpyrazine pH-5 arm has r2 0.890


def systems() -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for t_c in (100.0, 110.0, 120.0):
        out[f"zhou_go_{int(t_c)}"] = dict(initial={"GO": 20.0, "Gly": 20.0}, t_c=t_c, ph=PYRAZINE_FIT_PH, anchor=ZHOU_ANCHOR)
        out[f"zhou_mgo_{int(t_c)}"] = dict(initial={"MGO": 20.0, "Gly": 20.0}, t_c=t_c, ph=PYRAZINE_FIT_PH, anchor=ZHOU_ANCHOR)
    for ph in (9.0, 7.0, 5.0):
        out[f"leahy_95C_ph{int(ph)}"] = dict(initial={"Glc": 100.0, "Gly": 100.0}, t_c=95.0, ph=ph, anchor=LEAHY_ANCHOR)
    return out


SYSTEMS = systems()


def rows() -> Tuple[Dict[str, Any], ...]:
    out: List[Dict[str, Any]] = []
    for species, route in (("PZ", "go"), ("DMP", "mgo")):
        for t_c, k in ZHOU_RATES_UMOL_L_MIN[species].items():
            out.append(dict(id=f"zhou_{species}_rate_{int(t_c)}C", kind="mean_rate_umol_l_min", system=f"zhou_{route}_{int(t_c)}",
                            species=species, target=k, sigma_log=SIGMA_ZHOU, anchor=ZHOU_ANCHOR,
                            note="MEASURED RATE (zero-order slope over 0-120 min); alanine -> glycine declared."))
    for species in ("PZ", "MPZ"):
        for ph, ratio in LEAHY_RATIOS[species].items():
            out.append(dict(id=f"leahy_{species}_k_ph{int(ph)}_over_ph9", kind="rate_ratio", system=f"leahy_95C_ph{int(ph)}",
                            system_b="leahy_95C_ph9", species=species, target=ratio,
                            sigma_log=SIGMA_LEAHY_WEAK if (species == "MPZ" and ph == 5.0) else SIGMA_LEAHY,
                            anchor=LEAHY_ANCHOR, note="WITHIN-STUDY RATIO at 95 C; lysine -> glycine declared; the buffer changes between arms."))
    return tuple(out)


ROWS = rows()
assert len(ROWS) == 10

# ===========================================================================
# 2. THE VECTOR
# ===========================================================================
KEYS: Tuple[str, ...] = ("log10_k_go_ak_100C", "ea_go_ak_kj_mol", "log10_k_mgo_ak_100C", "ea_mgo_ak_kj_mol",
                         "ph_slope_above_7_decades_per_unit", "ph_slope_below_7_decades_per_unit")
PRIOR = np.array([FROZEN_B18[k] for k in KEYS], dtype=float)
#: bands: the log10 constants two decades either side of the prior centre (the transfer band and the
#: Strecker-over-two convention are inside it); the barriers from the printed value to the dossier's
#: three-point refit; the slopes non-negative and below 1.5 decades per unit.
LOWER = np.array([PRIOR[0] - 2.0, 100.59, PRIOR[2] - 2.0, 111.66, 0.0, 0.0])
UPPER = np.array([PRIOR[0] + 2.0, 103.1, PRIOR[2] + 2.0, 114.9, 1.5, 1.5])
BASE_PARAMETERS = dict(operative_parameters(b1_fitted()))
#: `--variant nosink` (INFORMATION ONLY, cannot ship): the B13 dry-glass glyoxal sink (k_go_sink, measured at
#: 180 C with its barrier fixed to zero by the authors) set to zero, to size how much of the fitted glyoxal
#: Strecker constant is the sink's doing. Zhou's fed-glyoxal pot shows linear pyrazine growth over 120 min
#: (figure-only), which a sink that removes 98 % of the glyoxal in two hours cannot produce.
VARIANT = "b18"


def parameters_for(x: np.ndarray) -> Dict[str, Any]:
    p = dict(BASE_PARAMETERS)
    p.update(with_fitted_pyrazine(float(x[0]), float(x[1]), float(x[2]), float(x[3])))
    if VARIANT == "nosink":
        from dataclasses import replace as _replace
        p["k_go_sink"] = _replace(p["k_go_sink"], k_ref=0.0)
    return p


def _mean_rate(run, species: str) -> float:
    """umol L-1 min-1 over the window: the endpoint over the window (the source's regression slope)."""
    c = run.series(species)
    return float(c[-1] - c[0]) / WINDOW_MIN * 1000.0


def simulate(x: np.ndarray, rtol: float = 1e-8) -> Dict[str, Any]:
    p = parameters_for(x)
    out: Dict[str, Any] = {}
    slopes = (float(x[4]), float(x[5]))
    for name, spec in SYSTEMS.items():
        proc = SimpleNamespace(ph=float(spec["ph"]), water_activity=None)
        params, notes = trunk_conditions.apply(p, proc, pyrazine_slopes=slopes)
        run = integrate(params, float(spec["t_c"]) + CELSIUS, spec["initial"], GRID, rtol=rtol, atol=1e-14)
        out[name] = run
    return out


def predictions(x: np.ndarray, runs=None) -> Dict[str, float]:
    runs = simulate(x) if runs is None else runs
    pred: Dict[str, float] = {}
    for r in ROWS:
        if r["kind"] == "mean_rate_umol_l_min":
            pred[r["id"]] = _mean_rate(runs[r["system"]], r["species"])
        else:
            a = _mean_rate(runs[r["system"]], r["species"])
            b = _mean_rate(runs[r["system_b"]], r["species"])
            pred[r["id"]] = a / b if b > 0 else float("nan")
    return pred


def residuals(x: np.ndarray) -> np.ndarray:
    pred = predictions(x)
    out = np.empty(len(ROWS))
    for i, r in enumerate(ROWS):
        p, t = pred[r["id"]], float(r["target"])
        v = math.log10((p + 1e-15) / (t + 1e-15)) / float(r["sigma_log"]) if (math.isfinite(p) and p > 0) else 25.0
        out[i] = max(-25.0, min(25.0, v))
    return out


def cost(x: np.ndarray) -> float:
    r = residuals(x)
    return float(np.sum(r * r))


# ===========================================================================
# 3. THE FIT
# ===========================================================================
def start_vector(start: int) -> np.ndarray:
    if start == 0:
        return PRIOR.copy()
    rng = np.random.default_rng(SEED + start)
    x = PRIOR.copy()
    x[0] += rng.uniform(-1.0, 1.0)
    x[2] += rng.uniform(-1.0, 1.0)
    x[1] = rng.uniform(LOWER[1], UPPER[1])
    x[3] = rng.uniform(LOWER[3], UPPER[3])
    x[4] = rng.uniform(0.05, 0.6)
    x[5] = rng.uniform(0.2, 1.2)
    return np.clip(x, LOWER, UPPER)


def fit_member(start: int, max_nfev: int = 200) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    x0 = start_vector(start)
    evals = {"n": 0}

    def f(x):
        evals["n"] += 1
        return residuals(x)

    res = least_squares(f, x0, bounds=(LOWER, UPPER), method="trf", x_scale="jac", max_nfev=max_nfev,
                        diff_step=1e-3, ftol=1e-10, xtol=1e-10)
    r = residuals(res.x)
    return {"start": start, "x0": [float(v) for v in x0], "x": [float(v) for v in res.x], "cost": float(np.sum(r * r)),
            "evals": evals["n"], "status": int(res.status), "message": str(res.message),
            "residual_by_row": {row["id"]: float(v) for row, v in zip(ROWS, r)}}


def laplace(x: np.ndarray, r: np.ndarray) -> Dict[str, Any]:
    steps = np.array([0.05, 1.0, 0.05, 1.0, 0.05, 0.05])
    J = np.zeros((len(ROWS), len(KEYS)))
    for j in range(len(KEYS)):
        h = steps[j]
        lo_ok, hi_ok = x[j] - h >= LOWER[j], x[j] + h <= UPPER[j]
        if lo_ok and hi_ok:
            xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
            J[:, j] = (residuals(xp) - residuals(xm)) / (2 * h)
        elif hi_ok:
            xp = x.copy(); xp[j] += h
            J[:, j] = (residuals(xp) - r) / h
        else:
            xm = x.copy(); xm[j] -= h
            J[:, j] = (r - residuals(xm)) / h
    dof = max(len(ROWS) - len(KEYS), 1)
    chi2_red = float(np.sum(r * r)) / dof
    jtj = J.T @ J
    sigma2 = np.linalg.pinv(jtj) * chi2_red
    sigma = np.sqrt(np.clip(np.diag(sigma2), 0.0, None))
    thresholds = np.array([3.0, 60.0, 3.0, 60.0, 1.5, 1.5])
    width = UPPER - LOWER
    on_bound = [(bool(x[j] - LOWER[j] <= 1e-3 * width[j]) or bool(UPPER[j] - x[j] <= 1e-3 * width[j])) for j in range(len(KEYS))]
    identified = [bool(sigma[j] <= thresholds[j]) and not (on_bound[j] and sigma[j] > 0.1 * width[j]) for j in range(len(KEYS))]
    return {"method": "Laplace at the optimum: J by central differences on the residual vector, Sigma = pinv(J^T J) chi2_red",
            "chi2_reduced": chi2_red, "dof": dof, "sigma": [float(v) for v in sigma], "on_bound": on_bound,
            "identified": identified, "thresholds": [float(v) for v in thresholds],
            "jtj_rank": int(np.linalg.matrix_rank(jtj)), "covariance": [[float(v) for v in row] for row in sigma2]}


# ===========================================================================
# 4. DIAGNOSTICS THE PREREG ASKS FOR
# ===========================================================================
def diagnostics(x: np.ndarray) -> Dict[str, Any]:
    runs = simulate(x)
    out: Dict[str, Any] = {}
    # the glyoxal and methylglyoxal that remain at the end of Zhou's window under the trunk's own sinks
    for name in ("zhou_go_100", "zhou_go_120", "zhou_mgo_100", "zhou_mgo_120"):
        run = runs[name]
        key = "GO" if "_go_" in name else "MGO"
        c = run.series(key)
        out[f"{name}_{key}_remaining_fraction_at_120min"] = float(c[-1] / c[0])
        ak = run.series("AKG" if key == "GO" else "AKM")
        out[f"{name}_aminoketone_mmol_per_l_at_120min"] = float(ak[-1])
    # the product's linearity in the window: the 60-min rate over the 120-min rate (1 = zero order)
    for name, species in (("zhou_go_120", "PZ"), ("zhou_mgo_120", "DMP")):
        c = runs[name].series(species)
        i60 = int(np.argmin(np.abs(GRID - 60.0)))
        out[f"{name}_{species}_rate_0_60_over_0_120"] = float((c[i60] - c[0]) / 60.0 / ((c[-1] - c[0]) / WINDOW_MIN))
    # sensitivity to the declared condensation constant (x10 and /10)
    base = predictions(x, runs)
    from src.kinetic_core import parameters_pyrazine as PP
    from dataclasses import replace as _replace
    sens = {}
    for factor in (0.1, 10.0):
        p = parameters_for(x)
        p["k_cond"] = _replace(PP.K_COND, k_ref=K_COND_DECLARED_L_PER_MMOL_MIN * factor)
        slopes = (float(x[4]), float(x[5]))
        pred = {}
        for rid, sysname, species in (("zhou_PZ_rate_120C", "zhou_go_120", "PZ"), ("zhou_DMP_rate_120C", "zhou_mgo_120", "DMP")):
            spec = SYSTEMS[sysname]
            params, _ = trunk_conditions.apply(p, SimpleNamespace(ph=float(spec["ph"]), water_activity=None), pyrazine_slopes=slopes)
            run = integrate(params, float(spec["t_c"]) + CELSIUS, spec["initial"], GRID, rtol=1e-8, atol=1e-14)
            pred[rid] = _mean_rate(run, species)
        sens[f"k_cond_x{factor:g}"] = {rid: {"rate": v, "log10_change": math.log10(v / base[rid]) if v > 0 and base[rid] > 0 else None}
                                        for rid, v in pred.items()}
    out["k_cond_sensitivity"] = sens
    return out


def _git_head() -> Dict[str, str]:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        return {"head": sha}
    except Exception as exc:  # pragma: no cover
        return {"head": f"unavailable: {exc}"}


# ===========================================================================
# 5. THE REPORT
# ===========================================================================
def build(max_nfev: int) -> Dict[str, Any]:
    members = [fit_member(0, max_nfev), fit_member(1, max_nfev)]
    best = min(members, key=lambda m: m["cost"])
    x = np.array(best["x"], dtype=float)
    r = residuals(x)
    pred = predictions(x)
    lap = laplace(x, r)
    width = UPPER - LOWER
    active = [{"key": KEYS[j], "bound": "lower" if x[j] - LOWER[j] <= 1e-3 * width[j] else "upper", "value": float(x[j])}
              for j in range(len(KEYS)) if (x[j] - LOWER[j] <= 1e-3 * width[j] or UPPER[j] - x[j] <= 1e-3 * width[j])]
    frozen = {k: float(v) for k, v in zip(KEYS, x)}
    return {
        "wave": f"{WAVE} -- a pyrazine step on the trunk (Zhou 2024's rates, Leahy 1989's pH ladder)",
        "artifact": "kinetic_core_b18_fit_report",
        "generated_by": "scripts/generators/generate_kinetic_core_b18_fit.py",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "prereg": data_paths.rel(PREREG),
        "declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 27",
        "objective": {
            "form": "10 rows (6 Zhou 2024 mean formation rates + 4 Leahy 1989 within-study pH ratios), 6 free coordinates, "
                    "log10 residuals over the rows' sigma; least_squares (trf), two starts",
            "n_rows": len(ROWS), "n_free_parameters": len(KEYS), "final_cost": best["cost"],
            "observable": "mean formation rate over the source's regression window (0-120 min), umol L-1 min-1",
            "rows": [dict(r) for r in ROWS],
        },
        "declared": {
            "alanine_to_glycine": "Zhou 2024 fed alanine; the trunk holds glycine; the products are the same and the rate is not measured for the pair: +/- 0.5 dex on every pyrazine answer (parameters_pyrazine.PYRAZINE_TRANSFER_BAND_DECADES)",
            "lysine_to_glycine": "Leahy 1989 used lysine; only the within-study pH RATIOS are read, not the levels",
            "condensation": f"k_cond = {K_COND_DECLARED_L_PER_MMOL_MIN:g} L/(mmol*min), no barrier, declared fast (Jousse 2002 R10 'fast'); sensitivity reported",
            "mixed_route": "2-methylpyrazine follows the two aminoketone pools statistically through the shared k_cond (2 sqrt of the two homo rates)",
            "ph_term": "two slopes, knot at pH 7, exactly 1 at the trunk's reference pH 6.8; the fitted constants are stored at 6.8 and Zhou's pH-8 pots carry the factor",
        },
        "bounds": {k: [float(LOWER[i]), float(UPPER[i])] for i, k in enumerate(KEYS)},
        "prior_centre": {k: float(PRIOR[i]) for i, k in enumerate(KEYS)},
        "members": members,
        "best_start": best["start"],
        "start_agreement": {k: abs(float(members[0]["x"][i]) - float(members[1]["x"][i])) for i, k in enumerate(KEYS)},
        "frozen_parameters": {"pyrazine": frozen},
        "predicted_by_row": pred,
        "residual_by_row_dex": {row["id"]: float(v) * float(row["sigma_log"]) for row, v in zip(ROWS, r)},
        "laplace": {**lap, "coordinates": list(KEYS)},
        "active_bounds": active,
        "diagnostics": diagnostics(x),
        "reference_temperature_K": 373.15,
        "reference_ph": 6.8,
    }


def render(payload: Dict[str, Any]) -> str:
    fr = payload["frozen_parameters"]["pyrazine"]
    lap = payload["laplace"]
    lines = [f"# {payload['wave']}", "",
             f"*Generated by `{payload['generated_by']}` on {payload['generated_on']}; pre-registration `{payload['prereg']}`; "
             f"{payload['declaration']}.*", "",
             "## The fitted coordinates", "",
             "| coordinate | value | Laplace sigma | identified | bound |", "|---|---:|---:|---|---|"]
    for i, k in enumerate(lap["coordinates"]):
        lines.append(f"| {k} | {fr[k]:.4f} | {lap['sigma'][i]:.3g} | {'yes' if lap['identified'][i] else 'no'} | "
                     f"{'on' if lap['on_bound'][i] else 'off'} |")
    lines += ["", f"Cost {payload['objective']['final_cost']:.3f} on {payload['objective']['n_rows']} rows, "
              f"{payload['objective']['n_free_parameters']} free; reduced chi2 {lap['chi2_reduced']:.2f}; best start {payload['best_start']}; "
              f"start agreement (max) {max(payload['start_agreement'].values()):.3g}.", "",
              "## The rows", "", "| row | target | predicted | residual (dex) |", "|---|---:|---:|---:|"]
    for row in payload["objective"]["rows"]:
        rid = row["id"]
        lines.append(f"| {rid} | {row['target']:.4g} | {payload['predicted_by_row'][rid]:.4g} | {payload['residual_by_row_dex'][rid]:+.3f} |")
    d = payload["diagnostics"]
    lines += ["", "## Diagnostics", ""]
    for k, v in d.items():
        if k != "k_cond_sensitivity":
            lines.append(f"- {k}: {v:.4g}" if isinstance(v, float) else f"- {k}: {v}")
    for tag, block in d["k_cond_sensitivity"].items():
        lines.append(f"- {tag}: " + "; ".join(f"{rid} {b['log10_change']:+.3f} dex" if b['log10_change'] is not None else f"{rid} n/a" for rid, b in block.items()))
    lines += ["", "## Declared", ""] + [f"- **{k}**: {v}" for k, v in payload["declared"].items()]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=200)
    parser.add_argument("--variant", choices=("b18", "nosink"), default="b18")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B18 is pre-registered; write the prereg before running the fit"
    global VARIANT, OUT_JSON, OUT_MD
    VARIANT = args.variant
    if VARIANT == "nosink":
        OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b18_nosink_fit_report.json"
        OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b18_nosink_fit_report.md"
    payload = build(args.max_nfev)
    payload["variant"] = VARIANT
    if VARIANT == "nosink":
        payload["wave"] += " (VARIANT nosink: the B13 dry-glass glyoxal sink zeroed; INFORMATION ONLY, cannot ship)"
        payload["generated_by"] += " --variant nosink"
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    OUT_MD.write_text(render(payload))
    fr = payload["frozen_parameters"]["pyrazine"]
    print(f"{WAVE}: cost {payload['objective']['final_cost']:.3f}; " + ", ".join(f"{k}={v:.4f}" for k, v in fr.items()))
    print(f"wrote {OUT_JSON} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
