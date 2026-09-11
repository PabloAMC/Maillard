"""
WAVE B38 (2026-09-11): the identifiability audit.

Pre-registration: results/validation/kinetic_core_b38_prereg.md (method and verdict
thresholds fixed there BEFORE this ran). Moves no constant. For every fit whose ship rule
says SHIP, at the optimum that shipped, using the fit's OWN residual vector:

  * a central-difference Jacobian J of the sigma-weighted residuals;
  * FIM = J^T J, chi2_red = 2 cost / dof, Sigma = pinv(FIM) chi2_red;
  * per coordinate: marginal sigma (everything else free), conditional sigma
    (everything else fixed), at-bound, verdict;
  * the sloppy directions: eigenvectors of the correlation-normalised FIM with
    eigenvalue < 1e-3 of the largest, named by their two heaviest loadings;
  * pairwise |corr| > 0.9 from Sigma;
  * the cross with the envelope's priors: prior-dominated vs data-dominated.

Run inside docker:  PYTHONPATH=/workspace python scripts/generators/generate_kinetic_core_b38_identifiability.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from src import data_paths, provenance
from src.kinetic_core import uncertainty as U

ROOT = data_paths.REPO_ROOT
VAL = data_paths.VALIDATION_DIR
OUT_JSON = VAL / "kinetic_core_b38_identifiability.json"
OUT_MD = VAL / "kinetic_core_b38_identifiability.md"
PREREG = VAL / "kinetic_core_b38_prereg.md"

# --- the pre-registered verdict thresholds (95 % half-width) ---------------------------
PINNED = {"log10k": 0.5, "ea": 30.0, "yield": 0.15, "pka": 0.5, "other": 0.5}
WEAK = {"log10k": 1.5, "ea": 90.0, "yield": 0.5, "pka": 1.5, "other": 1.5}
STEP = {"log10k": 1e-3, "ea": 0.5, "yield": 1e-3, "pka": 1e-3, "other": 1e-3}
NULL_EIG = 1e-3
COLLINEAR = 0.9


def _load(name: str):
    sys.argv = ["x"]
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "scripts" / "generators" / f"generate_kinetic_core_{name}_fit.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _kind(key: str) -> str:
    k = key.lower()
    if k.startswith("ea_") or "_ea_" in k or k.endswith("_ea") or "ea_kj" in k:
        return "ea"
    if "yield" in k:
        return "yield"
    if "pka" in k:
        return "pka"
    if k.startswith("k_") or k.startswith("log10_k") or "log10" in k:
        return "log10k"
    return "other"


def _jacobian(fn, x: np.ndarray, kinds: Sequence[str]) -> np.ndarray:
    r0 = np.asarray(fn(x), dtype=float)
    J = np.zeros((r0.size, x.size))
    for i in range(x.size):
        h = STEP[kinds[i]]
        xp, xm = x.copy(), x.copy()
        xp[i] += h; xm[i] -= h
        J[:, i] = (np.asarray(fn(xp), dtype=float) - np.asarray(fn(xm), dtype=float)) / (2 * h)
    return J, r0


def _audit(name: str, keys: Sequence[str], x: np.ndarray, fn, lower: Optional[np.ndarray],
           upper: Optional[np.ndarray], note: str) -> Dict[str, Any]:
    t0 = time.time()
    kinds = [_kind(k) for k in keys]
    J, r0 = _jacobian(fn, x, kinds)
    n, p = J.shape
    cost = 0.5 * float(r0 @ r0)
    dof = max(n - p, 1)
    chi2 = 2.0 * cost / dof
    fim = J.T @ J
    sigma = np.linalg.pinv(fim) * chi2
    marg = np.sqrt(np.clip(np.diag(sigma), 0.0, None))
    diag = np.diag(fim)
    cond = np.where(diag > 0, np.sqrt(chi2 / np.where(diag > 0, diag, 1.0)), np.inf)
    # correlation-normalised FIM for the eigen-analysis
    d = np.where(diag > 0, np.sqrt(diag), 1.0)
    fim_n = fim / np.outer(d, d)
    w, v = np.linalg.eigh(fim_n)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    wmax = float(w[0]) if w.size else 1.0
    null_dirs = []
    for j in range(w.size):
        if w[j] < NULL_EIG * wmax:
            load = np.argsort(np.abs(v[:, j]))[::-1][:3]
            null_dirs.append({
                "eigenvalue_over_max": float(w[j] / wmax) if wmax else None,
                "loadings": [{"key": keys[i], "weight": float(v[i, j])} for i in load],
            })
    # pairwise collinearity
    sd = np.sqrt(np.clip(np.diag(sigma), 1e-300, None))
    corr = sigma / np.outer(sd, sd)
    pairs = []
    for i in range(p):
        for j in range(i + 1, p):
            c = float(corr[i, j])
            if np.isfinite(c) and abs(c) >= COLLINEAR and np.isfinite(marg[i]) and np.isfinite(marg[j]):
                pairs.append({"a": keys[i], "b": keys[j], "corr": c})
    rows = []
    counts = {"AT_BOUND": 0, "PINNED": 0, "WEAK": 0, "UNIDENTIFIED": 0}
    collinear_not_insensitive = 0
    not_pinned = 0
    for i, k in enumerate(keys):
        kind = kinds[i]
        hw_m = 1.96 * float(marg[i])
        hw_c = 1.96 * float(cond[i])
        at_bound = False
        if lower is not None and upper is not None:
            span = max(float(upper[i] - lower[i]), 1e-12)
            at_bound = bool(min(abs(x[i] - lower[i]), abs(x[i] - upper[i])) <= 1e-6 * span)
        if at_bound:
            verdict = "AT_BOUND"
        elif hw_m <= PINNED[kind]:
            verdict = "PINNED"
        elif hw_m <= WEAK[kind]:
            verdict = "WEAK"
        else:
            verdict = "UNIDENTIFIED"
        counts[verdict] += 1
        cni = verdict in ("WEAK", "UNIDENTIFIED") and hw_c <= PINNED[kind]
        if verdict in ("WEAK", "UNIDENTIFIED"):
            not_pinned += 1
            collinear_not_insensitive += int(cni)
        insensitive = not np.isfinite(hw_c) or hw_c > WEAK[kind]
        rows.append({
            "key": k, "kind": kind, "value": float(x[i]),
            "ci95_halfwidth_marginal": hw_m if np.isfinite(hw_m) else None,
            "ci95_halfwidth_conditional": hw_c if np.isfinite(hw_c) else None,
            "sensitivity_norm": float(np.sqrt(diag[i])),
            "at_bound": at_bound, "verdict": verdict,
            "collinear_not_insensitive": bool(cni),
            "insensitive": bool(insensitive),
        })
    return {
        "fit": name, "note": note, "n_rows": n, "n_free": p, "dof": dof,
        "cost_at_optimum": cost, "reduced_chi_square": chi2,
        "verdict_counts": counts,
        "not_pinned": not_pinned, "collinear_not_insensitive": collinear_not_insensitive,
        "eigenvalues_over_max": [float(x_ / wmax) for x_ in w] if wmax else [],
        "null_directions": null_dirs, "collinear_pairs": pairs,
        "coordinates": rows, "wall_seconds": round(time.time() - t0, 1),
    }


# ---------------------------------------------------------------------------------------
def audit_b8() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    m = _load("b8")
    lap = json.loads((VAL / "kinetic_core_b8_laplace_covariance.json").read_text())
    free_keys = list(m.FREE_KEYS)
    assert list(lap["keys"]) == free_keys, "B8 laplace keys differ from FREE_KEYS"
    x_full = np.asarray(m.incumbent_vector(), dtype=float)
    x_free = np.asarray(lap["optimum"], dtype=float)
    idx = list(m.FREE_INDEX)
    x_full[idx] = x_free
    lower, upper = m.full_bounds()
    lower_f, upper_f = np.asarray(lower)[idx], np.asarray(upper)[idx]

    def fn_free(xf):
        xx = x_full.copy(); xx[idx] = xf
        return m.residual_vector(xx, True)

    free = _audit("B8 (sulfur, 23 free)", free_keys, x_free, fn_free, lower_f, upper_f,
                  "the shipped sulfur vector: B8's 23 free coordinates at its Laplace optimum on top of B2.4-half")
    # the 25 FROZEN coordinates: sensitivity only (are any of them data-visible?)
    all_keys = list(m.ALL_KEYS)
    frozen_idx = [i for i in range(len(all_keys)) if i not in idx]
    fkeys = [all_keys[i] for i in frozen_idx]
    kinds = [_kind(k) for k in fkeys]
    J, r0 = _jacobian(lambda xf: (lambda xx: m.residual_vector(xx, True))(_with(x_full, frozen_idx, xf)),
                      x_full[frozen_idx], kinds)
    chi2 = free["reduced_chi_square"]
    diag = np.diag(J.T @ J)
    frozen_rows = []
    for i, k in enumerate(fkeys):
        hw_c = 1.96 * (np.sqrt(chi2 / diag[i]) if diag[i] > 0 else np.inf)
        frozen_rows.append({"key": k, "kind": kinds[i], "value": float(x_full[frozen_idx[i]]),
                            "ci95_halfwidth_conditional": hw_c if np.isfinite(hw_c) else None,
                            "sensitivity_norm": float(np.sqrt(diag[i])),
                            "data_visible_if_freed": bool(np.isfinite(hw_c) and hw_c <= PINNED[kinds[i]])})
    return free, {"fit": "B8 frozen coordinates (25)", "coordinates": frozen_rows,
                  "n_data_visible_if_freed": sum(r["data_visible_if_freed"] for r in frozen_rows)}


def _with(x_full, idx, xf):
    xx = x_full.copy(); xx[idx] = xf; return xx


def audit_b3() -> Dict[str, Any]:
    m = _load("b3")
    rep = json.loads((VAL / "kinetic_core_b3_fit_report.json").read_text())
    keys = list(m.PARAM_ORDER) + list(m.EA_ORDER)
    x = np.array([rep["parameter_intervals"][k]["value"] for k in keys], dtype=float)
    from src.kinetic_core.parameters_acrylamide import FITTED_ACRYLAMIDE_BOUNDS_LOG10K, FITTED_ACRYLAMIDE_EA_BOUNDS
    lower = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][0] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[0]] * len(m.EA_ORDER))
    upper = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][1] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[1]] * len(m.EA_ORDER))
    return _audit("B3 (acrylamide, 11 free)", keys, x, lambda xx: m.residuals(xx, quick=True), lower, upper,
                  "the shipped acrylamide vector at the B3 report's optimum")


def audit_simple(name: str, label: str, note: str) -> Dict[str, Any]:
    m = _load(name)
    rep = json.loads((VAL / f"kinetic_core_{name}_fit_report.json").read_text())
    best = min(rep["members"], key=lambda d: d["cost"])
    x = np.asarray(best["x"], dtype=float)
    keys = list(m.KEYS)
    lower = np.asarray(m.LOWER, dtype=float); upper = np.asarray(m.UPPER, dtype=float)
    return _audit(label, keys, x, m.residuals, lower, upper, note)


def read_b1() -> Dict[str, Any]:
    rep = json.loads((VAL / "kinetic_core_b1_fit_report.json").read_text())
    found: List[Dict[str, Any]] = []

    def walk(o, path=""):
        if isinstance(o, dict):
            if any(k in o for k in ("se_log10", "standard_error", "verdict")) and ("value" in o or "estimate" in o or "log10_k" in o):
                found.append({"path": path, **{k: o[k] for k in o if k in ("value", "estimate", "log10_k", "se_log10", "standard_error", "verdict", "unit")}})
            for k, v in o.items():
                walk(v, f"{path}/{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]")
    walk(rep)
    counts: Dict[str, int] = {}
    for f in found:
        v = str(f.get("verdict", "?")); counts[v] = counts.get(v, 0) + 1
    return {"fit": "B1 (trunk) -- read from its own report, not recomputed",
            "n_entries_with_standard_errors": len(found), "verdict_counts": counts, "entries": found[:60]}


def cross_with_priors(audits: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    by_key: Dict[str, Dict[str, Any]] = {}
    for a in audits:
        for r in a.get("coordinates", []):
            if "verdict" in r:
                by_key[r["key"].lower()] = {**r, "fit": a["fit"]}
    rows = []
    for p in U.core_priors():
        segs = p.key.lower().split(".")
        hit = None
        for n in range(1, len(segs)):
            hit = by_key.get(".".join(segs[n:]))
            if hit is not None:
                break
        if hit is None:
            # b8 priors are "b8.<coordinate>.<block>": the coordinate is the MIDDLE segment; a barrier
            # prior is "b8.<family>.ea_kj_mol" against a coordinate named Ea_decay_<family>
            for seg in segs[1:]:
                hit = by_key.get(seg) or (by_key.get(f"ea_decay_{seg}") if segs[-1].startswith("ea") else None)
                if hit is not None:
                    break
        if hit is None:
            continue
        band = getattr(p, "band", None); sig = getattr(p, "sigma", None)
        if band is not None:
            half = 0.5 * abs(float(band[1]) - float(band[0]))
        elif sig is not None:
            half = 1.96 * float(sig)
        else:
            continue
        hw = hit["ci95_halfwidth_marginal"]
        dom = "prior-dominated" if (hw is None or half < hw) else "data-dominated"
        rows.append({"prior": p.key, "coordinate": hit["key"], "fit": hit["fit"], "prior_halfwidth": half,
                     "data_ci95_halfwidth": hw, "verdict": hit["verdict"], "dominance": dom,
                     "reason": str(getattr(p, "reason", ""))[:60]})
    n = len(rows); pd = sum(r["dominance"] == "prior-dominated" for r in rows)
    return {"n_priors_on_fitted_coordinates": n, "prior_dominated": pd, "data_dominated": n - pd, "rows": rows}


def main() -> int:
    t0 = time.time()
    b8_free, b8_frozen = audit_b8()
    b3 = audit_b3()
    b18 = audit_simple("b18", "B18 (dicarbonyl sinks, 6 free)", "the shipped B18 vector at its best member")
    b20 = audit_simple("b20", "B20 (glycation, 5 free)", "the shipped B20 vector at its best member")
    b21 = audit_simple("b21", "B21 (aqueous glucosone -> glyoxal, 2 free)", "the shipped B21 vector at its best member")
    b1 = read_b1()
    audits = [b8_free, b3, b18, b20, b21]
    cross = cross_with_priors(audits)
    total = sum(a["n_free"] for a in audits)
    pinned = sum(a["verdict_counts"]["PINNED"] for a in audits)
    at_bound = sum(a["verdict_counts"]["AT_BOUND"] for a in audits)
    not_pinned = sum(a["not_pinned"] for a in audits)
    cni = sum(a["collinear_not_insensitive"] for a in audits)
    b3_null = len(b3["null_directions"])
    b3_comp = {r["key"]: r for r in b3["coordinates"] if r["key"] in ("k_gln_glc", "k_ala_glc", "k_acr_ala")}
    b8_thiol = next(r for r in b8_free["coordinates"] if r["key"] == "Ea_decay_thiol_sink")
    b8_bound_ea = [r["key"] for r in b8_free["coordinates"] if r["kind"] == "ea" and r["at_bound"]]
    b18_ea = [r for r in b18["coordinates"] if r["kind"] == "ea"]
    predictions = {
        "P1_fewer_than_40pct_pinned": {"pinned": pinned, "total_free": total, "fraction": pinned / total, "held": pinned / total < 0.40},
        "P2_b3_three_null_and_competitors_insensitive": {
            "null_directions": b3_null,
            "competitors": {k: {"ci95_conditional": v["ci95_halfwidth_conditional"], "insensitive": v["insensitive"]} for k, v in b3_comp.items()},
            "held": b3_null >= 3 and all(v["ci95_halfwidth_conditional"] is None or v["ci95_halfwidth_conditional"] > 1.0 for v in b3_comp.values())},
        "P3_b8_thiol_sink_at_bound_plus_one_more_barrier": {"thiol_sink": {"value": b8_thiol["value"], "at_bound": b8_thiol["at_bound"]},
                                                            "barriers_at_bound": b8_bound_ea, "held": b8_thiol["at_bound"] and len(b8_bound_ea) >= 2},
        "P4_collinearity_dominant": {"not_pinned": not_pinned, "collinear_not_insensitive": cni,
                                     "held": not_pinned > 0 and cni >= 0.5 * not_pinned},
        "P5_priors_mostly_prior_dominated": {"n": cross["n_priors_on_fitted_coordinates"], "prior_dominated": cross["prior_dominated"],
                                             "held": cross["n_priors_on_fitted_coordinates"] > 0 and cross["prior_dominated"] > 0.5 * cross["n_priors_on_fitted_coordinates"]},
        "P6_b18_barriers_pinned_by_band": {"barriers": [{k: r[k] for k in ("key", "value", "verdict", "ci95_halfwidth_marginal", "at_bound")} for r in b18_ea],
                                           "held": all(r["at_bound"] or (r["ci95_halfwidth_marginal"] or 1e9) > 3.2 for r in b18_ea)},
        "P7_nothing_moves": {"held": True, "how": "this generator writes one artifact and reads the rest"},
    }
    payload = {
        "provenance": provenance.provenance_block(
            "kinetic_core_b38_identifiability", generated_by="scripts/generators/generate_kinetic_core_b38_identifiability.py",
            wave="B38", inputs=[PREREG,
                                VAL / "kinetic_core_b8_laplace_covariance.json", VAL / "kinetic_core_b3_fit_report.json",
                                VAL / "kinetic_core_b18_fit_report.json", VAL / "kinetic_core_b20_fit_report.json",
                                VAL / "kinetic_core_b21_fit_report.json", VAL / "kinetic_core_b1_fit_report.json"]),
        "method": "central-difference Jacobian of each shipped fit's own sigma-weighted residuals at its shipped optimum; "
                  "FIM = J^T J; Sigma = pinv(FIM) * chi2_red; marginal vs conditional 95 % half-widths; eigen-analysis of the "
                  "correlation-normalised FIM; |corr| >= 0.9 pairs; cross with uncertainty.core_priors()",
        "thresholds": {"pinned_ci95_halfwidth": PINNED, "weak_ci95_halfwidth": WEAK, "null_eigenvalue_fraction": NULL_EIG, "collinear_abs_corr": COLLINEAR},
        "summary": {"total_free_coordinates": total, "pinned": pinned, "at_bound": at_bound, "not_pinned": not_pinned,
                    "collinear_not_insensitive": cni, "fits": [a["fit"] for a in audits]},
        "predictions": predictions,
        "fits": audits, "b8_frozen": b8_frozen, "b1_from_report": b1, "priors_cross": cross,
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    OUT_MD.write_text(render_md(payload))
    print(f"wrote {data_paths.rel(OUT_JSON)}: {total} free coordinates | pinned {pinned} | at bound {at_bound} | "
          f"collinear-not-insensitive {cni}/{not_pinned} | priors on fitted coords {cross['n_priors_on_fitted_coordinates']} "
          f"(prior-dominated {cross['prior_dominated']}) | {payload['wall_seconds']} s")
    for k, v in predictions.items():
        print(f"  {k}: {'HELD' if v['held'] else 'REFUTED'}")
    return 0


def render_md(p: Mapping[str, Any]) -> str:
    s = p["summary"]
    out = ["# Wave B38 — the identifiability audit", "",
           f"_Generated by `{p['provenance']['generated_by']}`; pre-registered in `kinetic_core_b38_prereg.md`. Moves no constant._", "",
           f"**{s['total_free_coordinates']} free coordinates** across the shipped fits: **{s['pinned']} pinned** by the data, "
           f"{s['at_bound']} sitting on a declared bound, {s['not_pinned']} weak or unidentified — of which "
           f"**{s['collinear_not_insensitive']} are collinear rather than insensitive** (the data would pin them if their neighbours were held).", "",
           "## Predictions", "", "| prediction | result |", "|---|---|"]
    for k, v in p["predictions"].items():
        out.append(f"| {k} | **{'HELD' if v['held'] else 'REFUTED'}** — " + json.dumps({kk: vv for kk, vv in v.items() if kk != 'held'}, default=str)[:300].replace("|", "/") + " |")
    for a in p["fits"]:
        out += ["", f"## {a['fit']}", "", f"_{a['note']}_ — {a['n_rows']} rows, {a['n_free']} free, χ²_red {a['reduced_chi_square']:.2f}, "
                f"verdicts {a['verdict_counts']}", "",
                "| coordinate | kind | value | 95 % half-width, marginal | conditional | at bound | verdict |", "|---|---|---:|---:|---:|---|---|"]
        for r in a["coordinates"]:
            hm = "∞" if r["ci95_halfwidth_marginal"] is None else f"{r['ci95_halfwidth_marginal']:.3g}"
            hc = "∞" if r["ci95_halfwidth_conditional"] is None else f"{r['ci95_halfwidth_conditional']:.3g}"
            flag = " (collinear)" if r["collinear_not_insensitive"] else (" (insensitive)" if r["insensitive"] and r["verdict"] != "PINNED" else "")
            out.append(f"| {r['key']} | {r['kind']} | {r['value']:.3g} | {hm} | {hc} | {'yes' if r['at_bound'] else ''} | {r['verdict']}{flag} |")
        if a["null_directions"]:
            out += ["", "Directions the data cannot see (eigenvalue < 10⁻³ of the largest), by heaviest loadings:", ""]
            for d in a["null_directions"]:
                out.append("- " + ", ".join(f"{l['key']} ({l['weight']:+.2f})" for l in d["loadings"]) + f" — λ/λmax = {d['eigenvalue_over_max']:.1e}")
        if a["collinear_pairs"]:
            out += ["", "Pairs with |corr| ≥ 0.9:", ""] + [f"- {q['a']} ~ {q['b']}: {q['corr']:+.2f}" for q in a["collinear_pairs"][:25]]
    f = p["b8_frozen"]
    out += ["", f"## {f['fit']}", "", f"Frozen in B8 but data-visible if freed (conditional 95 % half-width inside the PINNED threshold): "
            f"**{f['n_data_visible_if_freed']} of {len(f['coordinates'])}**.", "",
            "| coordinate | value | conditional 95 % half-width | visible if freed |", "|---|---:|---:|---|"]
    for r in f["coordinates"]:
        hc = "∞" if r["ci95_halfwidth_conditional"] is None else f"{r['ci95_halfwidth_conditional']:.3g}"
        out.append(f"| {r['key']} | {r['value']:.3g} | {hc} | {'yes' if r['data_visible_if_freed'] else ''} |")
    b1 = p["b1_from_report"]
    out += ["", f"## {b1['fit']}", "", f"Entries carrying standard errors in the B1 report: {b1['n_entries_with_standard_errors']}; verdicts {b1['verdict_counts']}."]
    c = p["priors_cross"]
    out += ["", "## The envelope's priors against the data", "",
            f"{c['n_priors_on_fitted_coordinates']} priors sit on fitted coordinates: **{c['prior_dominated']} prior-dominated**, {c['data_dominated']} data-dominated.", "",
            "| prior | coordinate | fit | prior half-width | data 95 % half-width | verdict | who dominates |", "|---|---|---|---:|---:|---|---|"]
    for r in c["rows"]:
        hw = "∞" if r["data_ci95_halfwidth"] is None else f"{r['data_ci95_halfwidth']:.3g}"
        out.append(f"| {r['prior']} | {r['coordinate']} | {r['fit'].split(' ')[0]} | {r['prior_halfwidth']:.3g} | {hw} | {r['verdict']} | {r['dominance']} |")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.exit(main())
