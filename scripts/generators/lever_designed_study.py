"""
The designed-study path of the lever table (PREREG amendments A4 and A7, rules R0 - R5).

A designed factorial run for the table (docs/guides/FACTORIAL_LEVER_STUDY.md) arrives as pots, not as keyed
literature contrasts. This module turns pots into the rows ``build_lever_table.cell_table`` fits:
- R0 the pot SD (fixed from the declared CV, never estimated);
- R1 censored pots and the rule for which contrasts enter;
- R2 the pH band at set levels and the holding check;
- R3 the corner sign check;
- R4 the shared-protocol cluster (rows carry ``lab`` = the cluster, ``design_lab`` = the laboratory, so the cluster
  counts as one laboratory and leave-one-laboratory-out leaves it out whole; the within-cluster term w lives in
  ``build_lever_table``'s fits);
- R5 the buffered pair and its entry into the buffer-type cell.

A pot is a mapping with: ``lab``, ``cluster``, ``study``, ``odorant``, ``sugar``, ``sulfur``, ``ph_set``,
``buffer`` ("matrix" or "phosphate"), ``final_ph`` (or None), ``value`` (amount, or None when not detected),
``lod`` (that run's detection limit, or None), ``isotope_dilution`` (bool).

Implementation reading of A7 for n pots per corner arm (written before any P3 value; see the A7 implementation note
in PREREG): "one pot not detected" is "some but not all", "both pots" is "all", and sigma_pot^2 / 2 is
sigma_pot^2 / n.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy import optimize, stats

#: the factor levels of the P3 design (FACTORIAL_LEVER_STUDY.md); the first level of each is the reference
FACTORS = {"sugar identity": ("sugar", ("glucose", "ribose")),
           "sulfur source": ("sulfur", ("cysteine", "thiamine")),
           "pH": ("ph_set", (6.0, 7.0))}
PH_BAND = "pH 6-7.5, initial or label"          # R2: by the nominal set levels
HELD_MIN_DIFF = 0.3                              # R2
Z90, Z975 = stats.norm.ppf(0.95), stats.norm.ppf(0.9875)   # 90 % two-sided; 97.5 % two-sided (Bonferroni over 4)
BUFFER_PAIR_CORNER = {"sugar": "glucose", "sulfur": "cysteine"}   # R5
BUFFER_MARGIN = 0.2                              # R5, dex per pH unit
MATRIX_BUFFERED = "matrix-buffered"              # R5: the buffer-type cell's level for the matrix-only side
BOUNDED_MAX = 2                                  # R1: more than 2 of 4 bounded corners and the contrast is out


# ------------------------------------------------------------------------------------------------ R0
def sigma_pot(isotope_dilution: bool) -> float:
    cv = 0.10 if isotope_dilution else 0.20
    return math.sqrt(0.1 ** 2 + math.log10(1.0 + cv) ** 2)


# ------------------------------------------------------------------------------------------------ R1
def corner_value(pots: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """The log10 value of one corner arm (its pots), its variance and its state:
    'measured' (all detected), 'censored' (some not detected: left-censored normal ML at sigma_pot),
    'bounded' (none detected: at the limit), 'excluded' (a not-detected pot without a limit, or no pots)."""
    if not pots:
        return {"state": "excluded", "reason": "no pots"}
    sig = sigma_pot(all(bool(p.get("isotope_dilution")) for p in pots))
    det = [math.log10(p["value"]) for p in pots if p.get("value") is not None and p["value"] > 0]
    nd = [p for p in pots if not (p.get("value") is not None and p["value"] > 0)]
    if any(p.get("lod") is None or p["lod"] <= 0 for p in nd):
        return {"state": "excluded", "reason": "a not-detected pot has no detection limit", "n": len(pots)}
    lims = [math.log10(p["lod"]) for p in nd]
    n = len(pots)
    if not nd:
        return {"state": "measured", "y": float(np.mean(det)), "v": sig ** 2 / n, "n": n, "censored_pots": 0}
    if not det:
        return {"state": "bounded", "y": max(lims), "v": sig ** 2 / n, "n": n, "censored_pots": n,
                "limits": sorted(set(round(x, 9) for x in lims))}

    def nll(m):
        return -(stats.norm.logpdf(det, m, sig).sum() + stats.norm.logcdf(lims, m, sig).sum())

    lo, hi = min(det + lims) - 10 * sig, max(det + lims) + 10 * sig
    m = float(optimize.minimize_scalar(nll, bounds=(lo, hi), method="bounded",
                                       options={"xatol": 1e-10}).x)
    h = 1e-4
    info = (nll(m + h) - 2 * nll(m) + nll(m - h)) / h ** 2       # observed Fisher information
    return {"state": "censored", "y": m, "v": 1.0 / info, "n": n, "censored_pots": len(nd)}


def _corners(lever: str) -> List[Tuple[Tuple[str, Any], ...]]:
    """The 4 settings of the other two factors."""
    others = [FACTORS[k] for k in FACTORS if k != lever]
    return [((others[0][0], a), (others[1][0], b)) for a in others[0][1] for b in others[1][1]]


def _arm_pots(pots, lever, level, corner):
    key = FACTORS[lever][0]
    return [p for p in pots if p[key] == level and all(p[k] == v for k, v in corner)]


def a4_contrast(pots: Sequence[Mapping[str, Any]], lever: str) -> Dict[str, Any]:
    """One laboratory's A4 contrast for one odorant and lever (level 2 minus level 1, the mean over the 4 corners),
    with R1's entry rules, and its corner contrasts for R3. ``pots`` are that laboratory's matrix-only factorial pots
    for that odorant."""
    lo_lv, hi_lv = FACTORS[lever][1]
    corners = []
    for corner in _corners(lever):
        corners.append((corner, corner_value(_arm_pots(pots, lever, lo_lv, corner)),
                        corner_value(_arm_pots(pots, lever, hi_lv, corner))))
    counts = {"censored_pots": sum(c.get("censored_pots", 0) for _, a, b in corners for c in (a, b)),
              "bounded_corners": sum(c["state"] == "bounded" for _, a, b in corners for c in (a, b))}
    out = {"lever": lever, "levels": [lo_lv, hi_lv], **counts, "corners": []}
    if any(c["state"] == "excluded" for _, a, b in corners for c in (a, b)):
        return {**out, "enters": False, "reason": "a corner arm enters no estimate (no pots or a not-detected pot "
                                                  "without a limit)"}
    for corner, a, b in corners:
        both = a["state"] == "bounded" and b["state"] == "bounded"
        if both and a["limits"] != b["limits"]:
            return {**out, "enters": False, "reason": "a corner bounded in both arms at different limits"}
        out["corners"].append({"corner": dict(corner), "both_bounded": both,
                               "y": None if both else b["y"] - a["y"], "v": None if both else a["v"] + b["v"]})
    ylo = float(np.mean([a["y"] for _, a, _ in corners]))
    yhi = float(np.mean([b["y"] for _, _, b in corners]))
    v = (sum(a["v"] for _, a, _ in corners) + sum(b["v"] for _, _, b in corners)) / 16.0
    n_lo = sum(a["state"] == "bounded" for _, a, _ in corners)
    n_hi = sum(b["state"] == "bounded" for _, _, b in corners)
    higher_bounded, lower_bounded = (n_hi, n_lo) if yhi >= ylo else (n_lo, n_hi)
    out.update({"y": yhi - ylo, "v": v})
    if higher_bounded > BOUNDED_MAX:
        return {**out, "enters": False, "reason": "the higher arm has more than 2 of 4 corners bounded"}
    if n_lo > BOUNDED_MAX and n_hi > BOUNDED_MAX:
        return {**out, "enters": False, "reason": "both arms have more than 2 of 4 corners bounded"}
    out["enters"] = True
    out["lower_bound"] = lower_bounded > 0       # biased toward zero: the size is a lower bound
    return out


# ------------------------------------------------------------------------------------------------ R2
def ph_held(pots: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """One laboratory's holding check over its matrix-only factorial pots (any one odorant's pots: pH is per pot)."""
    lo, hi = FACTORS["pH"][1]
    sel = {lv: [p.get("final_ph") for p in pots if p["buffer"] == "matrix" and p["ph_set"] == lv] for lv in (lo, hi)}
    if not sel[lo] or not sel[hi] or any(x is None for lv in sel for x in sel[lv]):
        return {"held": False, "reason": "a factorial pot has no final pH"}
    diff = float(np.mean(sel[hi]) - np.mean(sel[lo]))
    return {"held": diff >= HELD_MIN_DIFF, "final_ph_difference": diff,
            "reason": None if diff >= HELD_MIN_DIFF else f"final pH difference {diff:.2f} < {HELD_MIN_DIFF}"}


# ------------------------------------------------------------------------------------------------ R3
def corner_check(contrasts: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """R3 over the design laboratories' A4 contrasts of one cell (those that enter)."""
    cs = [c for c in contrasts if c.get("enters")]
    if not cs:
        return {"flag": False, "reason": "no contrast enters"}
    w = np.array([1.0 / c["v"] for c in cs])
    mu = float(np.sum(w * [c["y"] for c in cs]) / w.sum())
    se = math.sqrt(1.0 / w.sum())
    a4_excl = abs(mu) > Z90 * se
    pooled = []
    for k in range(len(cs[0]["corners"])):
        rows = [c["corners"][k] for c in cs if c["corners"][k]["y"] is not None]
        if not rows:
            pooled.append({"corner": cs[0]["corners"][k]["corner"], "y": None})
            continue
        wk = np.array([1.0 / r["v"] for r in rows])
        yk = float(np.sum(wk * [r["y"] for r in rows]) / wk.sum())
        sk = math.sqrt(1.0 / wk.sum())
        pooled.append({"corner": rows[0]["corner"], "y": yk, "ci975": [yk - Z975 * sk, yk + Z975 * sk]})
    opposite = [p for p in pooled if p["y"] is not None and a4_excl
                and ((mu > 0 and p["ci975"][1] < 0) or (mu < 0 and p["ci975"][0] > 0))]
    return {"flag": bool(a4_excl and opposite), "a4_mean": mu, "a4_ci90": [mu - Z90 * se, mu + Z90 * se],
            "corners": pooled, "label": "interaction-dominated" if a4_excl and opposite else None}


# ------------------------------------------------------------------------------------------------ R5
def buffered_pair(pots: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """One laboratory's buffered pair for one odorant: D = phosphate pH slope minus matrix-only pH slope at the
    glucose + cysteine corner (dex per pH unit), and the buffer-type contrast (matrix-buffered minus phosphate, the
    mean over the two pH levels)."""
    lo, hi = FACTORS["pH"][1]
    at = [p for p in pots if all(p[k] == v for k, v in BUFFER_PAIR_CORNER.items())]
    cv = {(b, lv): corner_value([p for p in at if p["buffer"] == b and p["ph_set"] == lv])
          for b in ("matrix", "phosphate") for lv in (lo, hi)}
    bad = [k for k, c in cv.items() if c["state"] in ("excluded", "bounded")]
    if bad:
        return {"enters": False, "reason": f"{len(bad)} of the 4 arms are bounded or enter no estimate",
                "arms": {f"{b} pH {lv:g}": c["state"] for (b, lv), c in cv.items()}}
    step = hi - lo
    slope = {b: (cv[(b, hi)]["y"] - cv[(b, lo)]["y"]) / step for b in ("matrix", "phosphate")}
    vs = {b: (cv[(b, hi)]["v"] + cv[(b, lo)]["v"]) / step ** 2 for b in ("matrix", "phosphate")}
    by = [cv[("matrix", lv)]["y"] - cv[("phosphate", lv)]["y"] for lv in (lo, hi)]
    bv = sum(c["v"] for c in cv.values()) / 4.0
    return {"enters": True, "D": slope["phosphate"] - slope["matrix"], "v": vs["phosphate"] + vs["matrix"],
            "slopes": slope, "buffer_contrast": {"y": float(np.mean(by)), "v": bv},
            "censored_pots": sum(c.get("censored_pots", 0) for c in cv.values())}


def buffered_pair_reading(per_lab: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    """R5's pooled reading over laboratories (inverse variance), with each laboratory shown alone."""
    ok = {k: d for k, d in per_lab.items() if d.get("enters")}
    if not ok:
        return {"reading": "no estimate", "labs": {k: d for k, d in per_lab.items()}}
    w = np.array([1.0 / d["v"] for d in ok.values()])
    D = float(np.sum(w * [d["D"] for d in ok.values()]) / w.sum())
    se = math.sqrt(1.0 / w.sum())
    ci = [D - Z90 * se, D + Z90 * se]
    if -BUFFER_MARGIN < ci[0] and ci[1] < BUFFER_MARGIN:
        r = "same"
    elif ci[0] > BUFFER_MARGIN or ci[1] < -BUFFER_MARGIN:
        r = "different"
    else:
        r = "unresolved"
    return {"reading": r, "D": D, "ci90": ci, "margin": BUFFER_MARGIN,
            "scope": "phosphate buffering together with phosphate catalysis; not separated",
            "labs": {k: {"D": d.get("D"), "se": math.sqrt(d["v"]) if d.get("enters") else None,
                         "enters": d.get("enters"), "reason": d.get("reason")} for k, d in per_lab.items()}}


# ------------------------------------------------------------------------------------------------ rows
def _base_row(compound, lever, lab, cluster, study, matrix):
    return {"source": f"designed:{study}", "contrast": f"{study}:{lab}:{lever}:{compound}", "compound": compound,
            "lever_class": lever, "lab": cluster, "design_lab": lab, "design_cluster": cluster,
            "partition": "design", "focal": True, "peak_area": False, "matrix": [matrix], "temps": [],
            "censored_arms": 0, "non_monotone": False, "turning_point": False, "range": None, "alt": {},
            "buffer_state": "unbuffered"}


def designed_rows(pots: Sequence[Mapping[str, Any]], resolve, matrix: str = "pea isolate"
                  ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Rows for ``cell_table`` from a designed study's pots, and the study report (R1 - R5 per cell and laboratory).
    ``resolve`` maps an odorant name to its display name (compound_keys)."""
    by = defaultdict(list)
    for p in pots:
        by[(resolve(p["odorant"]), p["lab"], p["cluster"], p["study"])].append(p)
    rows, report = [], {"contrasts": [], "ph_held": {}, "buffered_pair": defaultdict(dict), "corner_check": {},
                        "protocol": {}}
    a4 = defaultdict(list)
    for (comp, lab, cluster, study), ps in sorted(by.items()):
        matrix_pots = [p for p in ps if p["buffer"] == "matrix"]
        held = ph_held(matrix_pots)
        report["ph_held"][lab] = held
        for lever in FACTORS:
            c = a4_contrast(matrix_pots, lever)
            enters = c["enters"] and (lever != "pH" or held["held"])
            reason = c.get("reason") if not c["enters"] else (None if enters else f"pH not held: {held['reason']}")
            report["contrasts"].append({"compound": comp, "lab": lab, "lever": lever, "enters": enters,
                                        "reason": reason, "y": c.get("y"), "v": c.get("v"),
                                        "lower_bound": c.get("lower_bound", False),
                                        "censored_pots": c["censored_pots"], "bounded_corners": c["bounded_corners"]})
            a4[(comp, lever)].append({**c, "enters": enters})
            if not enters:
                continue
            r = _base_row(comp, lever, lab, cluster, study, matrix)
            r["censored_arms"] = c["bounded_corners"]
            r["lower_bound"] = c["lower_bound"]
            if lever == "pH":
                step = FACTORS["pH"][1][1] - FACTORS["pH"][1][0]
                r.update({"kind": "numeric", "y": c["y"] / step, "v": c["v"] / step ** 2, "band": PH_BAND,
                          "range": list(FACTORS["pH"][1])})
            else:
                lo_lv, hi_lv = FACTORS[lever][1]
                r.update({"kind": "categorical", "band": "all", "shared_var": 0.0, "arms_used": 2,
                          "pairs": [{"level": hi_lv, "against": lo_lv, "y": c["y"], "v": c["v"],
                                     "censored": c["bounded_corners"]}]})
            rows.append(r)
        bp = buffered_pair(ps)
        report["buffered_pair"][comp][lab] = bp
        if bp["enters"]:
            r = _base_row(comp, "buffer", lab, cluster, study, matrix)
            r.update({"kind": "categorical", "band": "all", "shared_var": 0.0, "arms_used": 2, "buffer_state":
                      "buffered", "pairs": [{"level": MATRIX_BUFFERED, "against": "phosphate",
                                             "y": bp["buffer_contrast"]["y"], "v": bp["buffer_contrast"]["v"],
                                             "censored": 0}]})
            rows.append(r)
    report["corner_check"] = {f"{comp} | {lever}": corner_check(cs) for (comp, lever), cs in a4.items()}
    report["protocol"] = {f"{comp} | {lever}": protocol_mean(cs) for (comp, lever), cs in a4.items()}
    report["buffered_pair"] = {comp: {"per_lab": d, "pooled": buffered_pair_reading(d)}
                               for comp, d in report["buffered_pair"].items()}
    return rows, report


def protocol_mean(contrasts: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """R4's pea-matrix report: the design laboratories' inverse-variance mean with a 90 % interval for THIS
    protocol. No new-laboratory interval."""
    cs = [c for c in contrasts if c.get("enters")]
    if not cs:
        return {"estimate": None}
    w = np.array([1.0 / c["v"] for c in cs])
    m = float(np.sum(w * [c["y"] for c in cs]) / w.sum())
    se = math.sqrt(1.0 / w.sum())
    return {"estimate": m, "ci90_this_protocol": [m - Z90 * se, m + Z90 * se], "labs": len(cs),
            "lower_bound": any(c.get("lower_bound") for c in cs), "new_lab_interval": None}
