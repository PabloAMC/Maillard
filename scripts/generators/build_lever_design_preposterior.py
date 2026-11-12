#!/usr/bin/env python
"""
Preposterior check of candidate lever studies (STATUS_AND_PLAN sec. 0, owner row): before any request is sent, how
many lever-table cells would each design answer under the PREREG sec. 5 rule?

For every design, many simulated worlds:
  1. each cell the design touches gets a true effect for a new laboratory drawn from the lever table's own
     posterior predictive for that cell (its 90 % interval of mu and its 90 % new-laboratory interval, read as
     normals); a cell or level with no keyed estimate draws from the PREREG prior scale (mu ~ N(0, 0.5), a
     laboratory offset with SD |N(0, 0.5)|);
  2. each new laboratory shares one offset per cell across its contrasts; each contrast adds a within-laboratory
     residual (interactions, 0.1 dex) and the measurement error of its arms (10 % CV per pot by SIDA for the thiols,
     20 % otherwise, the design's replicates);
  3. optionally, a matrix shift shared by every new laboratory (they all cook in pea isolate; the keyed record does
     not), SD 0.3 dex, so independent new laboratories are not over-counted;
  4. the synthetic contrasts are keyed as PREREG sec. 3 keys a real one (one contrast per pair of arms that differ in
     one lever) and the cell is refitted with the lever table's own model and reading rule
     (``build_lever_table.numeric_fit`` / ``categorical_fit`` / ``reading``).

Reported per design: the cells newly answered (mean, 90 % range, P(>= 1), P(>= 3)), the cells that newly reach two
TEST laboratories (STATUS_AND_PLAN step 5's precondition; deterministic), the experiment-map family x lever cells
newly answered, and each cell's probability of being answered. Levers the table has no class for (yeast extract)
are counted as feeding no cell.

Development and design aid only: it moves no shipped number, reads no engine prediction, and writes only
results/validation/lever_design_preposterior.{json,md}.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import build_lever_table as LT  # noqa: E402
from scripts.generators.build_experiment_map import FAMILIES  # noqa: E402
from src import artifact_io, data_access, data_paths, provenance  # noqa: E402
from src.kinetic_core import parallel  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "lever_design_preposterior.json"

N_SIMS = 400
WITHIN_LAB_SD = 0.10          # interactions / pot-to-pot structure inside one laboratory, dex
MATRIX_SDS = (0.0, 0.3)       # a pea-matrix shift shared by every new laboratory (sensitivity)
#: (matrix shift SD, keying) runs per design: the main-effect keying is the honest one; per-pair keying at the
#: larger matrix shift shows what the table's categorical model would claim if a factorial were keyed pair by pair
CONFIGS = ((0.0, "main_effect"), (0.3, "main_effect"), (0.3, "per_pair"))
PRIOR_SD = LT.PRIOR_SCALE     # mu and the laboratory-offset scale where a cell has no keyed estimate
SIDA_CV, OTHER_CV = 0.10, 0.20

THIOLS = ["2-Methyl-3-furanthiol (MFT)", "2-Furfurylthiol (FFT)", "Bis(2-methyl-3-furyl) disulfide",
          "Hydrogen Sulfide"]
OTHERS = ["2-Methylbutanal", "3-Methylbutanal", "Phenylacetaldehyde", "Methional (3-(methylthio)propanal)",
          "2,5-Dimethylpyrazine", "Trimethylpyrazine", "2-Ethyl-3,5-dimethylpyrazine", "Pyrazines",
          "Hexanal", "Nonanal", "2-Pentylfuran", "1-Octen-3-ol", "DMHF", "Norfuraneol (NF)", "Furfural"]
ANALYTES = THIOLS + OTHERS          # the FACTORIAL_LEVER_STUDY.md analyte list, as the lever table names them


# A design: laboratories, each a list of keyed contrast blocks. A block: (lever_class, levels in keyed order
# [reference first], number of such contrasts in the laboratory, replicate pots per arm). pH blocks are numeric
# (levels are pH values); every other block is categorical. "yeast extract" has no lever class (PREREG sec. 4).
def _factorial(levers: Sequence[str], reps: int, *, ph=(5.5, 6.5)) -> List[Tuple[str, List[Any], int, int]]:
    spec = {"sugar": ("sugar identity", ["glucose", "ribose"]),
            "sulfur": ("sulfur source", ["cysteine", "thiamine"]),
            "yeast": ("yeast extract", ["none", "1 %"]),
            "ph": ("pH", list(ph))}
    k = len(levers)
    return [(spec[f][0], spec[f][1], 2 ** (k - 1), reps) for f in levers]


DESIGNS: Dict[str, Dict[str, Any]] = {
    "P3: one laboratory, 2^4 x 2 (as drafted)": {
        "pots": 36, "labs": [_factorial(["sugar", "sulfur", "yeast", "ph"], 2)]},
    "two laboratories, 2^3 x 2 each (no yeast extract)": {
        "pots": 36, "labs": [_factorial(["sugar", "sulfur", "ph"], 2)] * 2},
    "three laboratories, 2^3 x 1 each (no yeast extract)": {
        "pots": 36, "labs": [_factorial(["sugar", "sulfur", "ph"], 1)] * 3},
    "three laboratories, targeted replication (sugar x sulfur source 2 x 3, x 2)": {
        "pots": 36, "labs": [[("sugar identity", ["glucose", "ribose"], 3, 2),
                              ("sulfur source", ["cysteine", "glutathione", "none"], 2, 2)]] * 3},
    "reference: P3 in two laboratories (72 pots)": {
        "pots": 72, "labs": [_factorial(["sugar", "sulfur", "yeast", "ph"], 2)] * 2},
}


def _cv(compound: str) -> float:
    return SIDA_CV if compound in THIOLS else OTHER_CV


def _arm_var(compound: str, reps: int) -> float:
    return math.log10(1.0 + _cv(compound)) ** 2 / reps


def _ph_band(levels) -> str:
    return f"pH {LT._band(LT.PH_BANDS, (min(levels) + max(levels)) / 2.0)}, initial or label"


# --------------------------------------------------------------------------------------------- cells and fits
def _key(o: Mapping[str, Any]) -> Tuple[str, str, str]:
    return (o["compound"], o["lever_class"], o["band"])


def readings(rows: Sequence[Mapping[str, Any]], lever_class: str) -> Dict[str, Dict[str, Any]]:
    """{level or "" (numeric): {reading, labs, test_labs, estimate, ci90, new_lab_pi90}} with the table's model."""
    labs = {r["lab"] for r in rows}
    if lever_class in LT.NUMERIC:
        test = {r["lab"] for r in rows if r["partition"] == "test"}
        s = LT.numeric_fit(list(rows), LT.PRIOR_SCALE)
        return {"": {"reading": LT.reading(len(labs), s), "labs": len(labs), "test_labs": len(test),
                     "estimate": s["estimate"], "ci90": s["ci90"], "new_lab_pi90": s["new_lab_pi90"],
                     "p_same_sign_new_lab": s["p_same_sign_new_lab"]}}
    f = LT.categorical_fit(list(rows), LT.CATEGORICAL_REFERENCE[lever_class], LT.PRIOR_SCALE)
    ref = f.get("reference", LT.CATEGORICAL_REFERENCE[lever_class])
    out = {}
    for lv, s in f["levels"].items():
        if "estimate" not in s:
            continue
        test = {r["lab"] for r in rows if r["partition"] == "test"
                for p in r["pairs"] if lv in (p["level"], p["against"])}
        out[f"{lv} vs {ref}"] = {"reading": LT.reading(s["labs"], s), "labs": s["labs"], "test_labs": len(test),
                                 "estimate": s["estimate"], "ci90": s["ci90"], "new_lab_pi90": s["new_lab_pi90"],
                                 "p_same_sign_new_lab": s["p_same_sign_new_lab"]}
    return out


def _truth_params(base: Mapping[str, Dict[str, Any]], label: str) -> Optional[Tuple[float, float, float]]:
    """(estimate, se of mu, new-laboratory SD beyond mu) from a keyed level, or None (use the prior)."""
    s = base.get(label)
    if not s:
        return None
    se = (s["ci90"][1] - s["ci90"][0]) / 3.29
    pw = (s["new_lab_pi90"][1] - s["new_lab_pi90"][0]) / 3.29
    return s["estimate"], se, math.sqrt(max(pw ** 2 - se ** 2, 0.0))


# --------------------------------------------------------------------------------------------- one simulation
def _synthetic(design: Mapping[str, Any], truths, matrix_sd: float, keying: str, rng,
               tau_factor: float = 1.0) -> List[Dict[str, Any]]:
    """The keyed contrasts of one simulated run of the design. ``keying``: "main_effect" keys a laboratory's
    2^(k-1) settings of one lever as ONE contrast (their mean: arm variance and within-laboratory residual divided by
    the number of settings), so one laboratory's pots do not pose as independent heterogeneity draws; "per_pair"
    keys every pair of arms as its own contrast (PREREG sec. 3 read literally)."""
    obs = []
    drawn_shift = defaultdict(lambda: float(rng.normal(0.0, matrix_sd)) if matrix_sd else 0.0)
    matrix_shift = _ShiftLookup(truths, drawn_shift)
    for li, lab_blocks in enumerate(design["labs"]):
        lab = f"design_lab_{li + 1}"
        lab_offset = {}
        for lc, levels, n_contrasts, reps in lab_blocks:
            if lc not in LT.NUMERIC and lc not in LT.CATEGORICAL_REFERENCE:
                continue                       # no lever class: keyed, enters no cell
            k_settings = n_contrasts if keying == "main_effect" else 1
            n_keyed = 1 if keying == "main_effect" else n_contrasts
            w_sd = WITHIN_LAB_SD / math.sqrt(k_settings)
            for compound in ANALYTES:
                av = _arm_var(compound, reps) / k_settings
                for ci in range(n_keyed):
                    common = {"source": f"{lab}_study", "contrast": f"{lab}-{lc}-{compound}-{ci}",
                              "compound": compound, "lever_class": lc, "lab": lab, "partition": "test",
                              "focal": True, "peak_area": compound not in THIOLS, "matrix": ["pea isolate"],
                              "temps": [140.0], "censored_arms": 0, "arms_used": len(levels)}
                    if lc in LT.NUMERIC:
                        band = _ph_band(levels)
                        key = (compound, lc, band, "")
                        mu = truths[key]
                        off = lab_offset.setdefault(key, float(rng.normal(0.0, tau_factor * truths[("tau",) + key])))
                        true = mu + off + matrix_shift[key]
                        step = max(levels) - min(levels)
                        y = (true * step + rng.normal(0.0, w_sd) + rng.normal(0.0, math.sqrt(2 * av))) / step
                        obs.append({**common, "kind": "numeric", "y": float(y), "v": 2 * av / step ** 2,
                                    "band": band, "non_monotone": False, "range": [min(levels), max(levels)]})
                    else:
                        ref = levels[0]
                        noise_ref = rng.normal(0.0, math.sqrt(av))
                        pairs = []
                        for lv in levels[1:]:
                            key = (compound, lc, "all", f"{lv} vs {ref}")
                            mu = truths[key]
                            off = lab_offset.setdefault(key, float(rng.normal(0.0, tau_factor * truths[("tau",) + key])))
                            y = (mu + off + matrix_shift[key] + rng.normal(0.0, w_sd)
                                 + rng.normal(0.0, math.sqrt(av)) - noise_ref)
                            pairs.append({"level": lv, "against": ref, "y": float(y), "v": 2 * av})
                        obs.append({**common, "kind": "categorical", "pairs": pairs, "shared_var": av,
                                    "band": "all"})
    return obs


class _ShiftLookup:
    """The matrix shift of a cell: the truth's own (predictive mode: shared with the new kitchen) or a fresh draw."""

    def __init__(self, truths, drawn):
        self.truths, self.drawn = truths, drawn

    def __getitem__(self, key):
        k = ("shift",) + key
        return self.truths[k] if k in self.truths else self.drawn[key]


def _targets(design) -> List[Tuple[str, str, str, str]]:
    """(compound, lever class, band, level label) of every cell / level the design keys."""
    out = []
    for lab_blocks in design["labs"]:
        for lc, levels, _n, _r in lab_blocks:
            if lc in LT.NUMERIC:
                out += [(c, lc, _ph_band(levels), "") for c in ANALYTES]
            elif lc in LT.CATEGORICAL_REFERENCE:
                out += [(c, lc, "all", f"{lv} vs {levels[0]}") for c in ANALYTES for lv in levels[1:]]
    return sorted(set(out))


def run_one(job) -> Dict[str, Any]:
    design_name, matrix_sd, keying, seed = job
    design = DESIGNS[design_name]
    rng = np.random.default_rng(seed)
    cells, base = _STATE["cells"], _STATE["base"]
    truths = {}
    for t in _targets(design):
        compound, lc, band, label = t
        p = _truth_params(base.get((compound, lc, band), {}), label)
        if p is None:
            truths[t] = float(rng.normal(0.0, PRIOR_SD))
            truths[("tau",) + t] = abs(float(rng.normal(0.0, PRIOR_SD)))
        else:
            est, se, tau = p
            truths[t] = float(rng.normal(est, se))
            truths[("tau",) + t] = tau
    synth = _synthetic(design, truths, matrix_sd, keying, rng)
    by_cell = defaultdict(list)
    for o in synth:
        by_cell[_key(o)].append(o)
    result = {}
    for key, new_rows in by_cell.items():
        rows = cells.get(key, []) + new_rows
        if len({r["lab"] for r in rows}) < LT.MIN_LABS:
            continue                                   # cannot be answered: no fit needed
        for label, s in readings(rows, key[1]).items():
            result["|".join(key + (label,))] = s["reading"]
    return result


_STATE: Dict[str, Any] = {}


def _load_state() -> None:
    sources = LT.load_sources()
    labs, _meta = LT.lab_partition()
    obs, _skipped = LT.observations(sources, labs)
    cells = defaultdict(list)
    for o in obs:
        cells[_key(o)].append(o)
    _STATE["cells"] = dict(cells)
    _STATE["base"] = {k: readings(v, k[1]) for k, v in cells.items()}


def _family(compound: str) -> Optional[str]:
    return next((f for f, names in FAMILIES.items() if compound in names), None)


def _map_lever(lc: str, label: str) -> List[str]:
    lv = {"buffer concentration": "buffer"}.get(lc, lc)
    return [lv, "thiamine"] if "thiamine" in label else [lv]


def summarise(design_name: str, matrix_sd: float, keying: str, sims: Sequence[Dict[str, str]]) -> Dict[str, Any]:
    base = _STATE["base"]
    design = DESIGNS[design_name]
    targets = _targets(design)
    # every level of every cell already answered (a refitted cell reports all its levels, not only the design's)
    answered_before = {"|".join(cell + (label,)) for cell, levels in base.items() for label, s in levels.items()
                       if s["reading"] == "answered by the literature"}
    map_before = {(f, lv) for k in answered_before for t in [k.split("|")] for f in [_family(t[0])] if f
                  for lv in _map_lever(t[1], t[3])}
    keyed_before = {"|".join(cell + (label,)) for cell, levels in base.items() for label in levels}
    new_counts, lit_counts, map_counts, lost = [], [], [], []
    per_cell = defaultdict(int)
    for sim in sims:
        ans = {k for k, r in sim.items() if r == "answered by the literature"}
        newly = ans - answered_before
        lost.append(len({k for k in answered_before if sim.get(k, "answered by the literature")
                         != "answered by the literature"}))
        new_counts.append(len(newly))
        lit_counts.append(len(newly & keyed_before))   # cells with a keyed literature estimate today
        for k in newly:
            per_cell[k] += 1
        mp = {(f, lv) for k in ans for t in [k.split("|")] for f in [_family(t[0])] if f
              for lv in _map_lever(t[1], t[3])}
        map_counts.append(len(mp - map_before))
    n = len(sims)
    nc = np.array(new_counts)
    mc = np.array(map_counts)
    # step 5's precondition: TEST contrasts from >= 2 laboratories (deterministic: the design's labs are TEST)
    n_new_labs = len(design["labs"])
    step5_before = sum(1 for t in targets if (base.get(t[:3], {}).get(t[3]) or {}).get("test_labs", 0) >= 2)
    step5_after = sum(1 for t in targets
                      if (base.get(t[:3], {}).get(t[3]) or {}).get("test_labs", 0) + n_new_labs >= 2)
    skipped_levers = sorted({lc for lab in design["labs"] for lc, *_ in lab
                             if lc not in LT.NUMERIC and lc not in LT.CATEGORICAL_REFERENCE})
    top = sorted(per_cell.items(), key=lambda kv: -kv[1])[:12]
    return {"design": design_name, "pots": design["pots"], "new_labs": n_new_labs, "matrix_shift_sd": matrix_sd,
            "keying": keying,
            "sims": n, "cells_keyed": len(targets),
            "newly_answered": {"mean": float(nc.mean()), "q05": float(np.quantile(nc, 0.05)),
                               "q95": float(np.quantile(nc, 0.95)), "p_ge_1": float(np.mean(nc >= 1)),
                               "p_ge_3": float(np.mean(nc >= 3))},
            "newly_answered_with_literature": {"mean": float(np.mean(lit_counts)),
                                               "p_ge_1": float(np.mean(np.array(lit_counts) >= 1))},
            "map_cells_newly_answered": {"mean": float(mc.mean()), "p_ge_1": float(np.mean(mc >= 1))},
            "answered_lost_mean": float(np.mean(lost)),
            "step5_cells_two_test_labs": {"before": step5_before, "after": step5_after},
            "levers_feeding_no_cell": skipped_levers,
            "top_cells": [{"cell": k.replace("|", ", ").rstrip(", "), "p_answered": v / n} for k, v in top]}


def render(p: Mapping[str, Any]) -> str:
    L = ["# Preposterior check of lever-study designs", "",
         "*Generated by `scripts/generators/build_lever_design_preposterior.py` (STATUS_AND_PLAN sec. 0, owner row). "
         "For each candidate study, simulated outcomes drawn from the lever table's own posterior predictive, keyed "
         "as PREREG sec. 3 keys a real contrast, and refitted with the table's model and reading rule. A design aid: "
         "it moves no shipped number and reads no engine prediction.*", "",
         f"Settings: {p['settings']['sims']} simulations per design and matrix setting; within-laboratory residual "
         f"{p['settings']['within_lab_sd']} dex; SIDA CV {p['settings']['sida_cv']}, other analytes "
         f"{p['settings']['other_cv']}; a cell with no keyed estimate draws mu and the laboratory-offset scale from "
         f"N(0, {p['settings']['prior_sd']}). Matrix shift: an offset shared by every new laboratory (all cook in pea "
         "isolate), SD 0 or 0.3 dex.", "",
         "| design | pots | new labs | matrix shift SD | keying | cells newly answered: mean [90 %] | P(>= 1) | P(>= 3) "
         "| of which with a literature estimate today (mean) "
         "| map cells newly answered (mean) | cells with 2 TEST labs (step 5): before -> after | levers feeding no "
         "cell |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in p["results"]:
        na = r["newly_answered"]
        L.append(f"| {r['design']} | {r['pots']} | {r['new_labs']} | {r['matrix_shift_sd']} | {r['keying']} | {na['mean']:.2f} "
                 f"[{na['q05']:.0f}, {na['q95']:.0f}] | {na['p_ge_1']:.2f} | {na['p_ge_3']:.2f} | "
                 f"{r['newly_answered_with_literature']['mean']:.2f} | "
                 f"{r['map_cells_newly_answered']['mean']:.2f} | {r['step5_cells_two_test_labs']['before']} -> "
                 f"{r['step5_cells_two_test_labs']['after']} | {', '.join(r['levers_feeding_no_cell']) or '-'} |")
    L += ["", "## The cells most likely to be answered, per design (matrix shift 0.3 dex)", ""]
    for r in p["results"]:
        if (r["matrix_shift_sd"], r["keying"]) != (max(MATRIX_SDS), "main_effect"):
            continue
        L.append(f"**{r['design']}:** " + ("; ".join(f"{c['cell']} ({c['p_answered']:.2f})" for c in r["top_cells"])
                                           or "none") + ".")
        L.append("")
    L += ["Reading: a cell is one odorant x lever (x pH band, x level for categorical levers) of `lever_table.md`; "
          "'answered' is PREREG sec. 5 (>= 3 laboratories, the 90 % interval of mu excludes zero, P(a new laboratory "
          "sees the same sign) >= 0.8). 'Map cells' are the experiment map's family x lever cells. Yeast extract has "
          "no lever class in PREREG sec. 4, so a yeast-extract factor feeds no cell until an amendment adds one.", ""]
    if p.get("predictive"):
        L += render_predictive(p["predictive"])
    return "\n".join(L)


# --------------------------------------------------------------------------------------------- predictive mode
# The question the "answered" count cannot settle without circularity (it counts our own three-laboratory rule): does
# the study make the table's prediction for a NEW kitchen cooking in pea isolate better? Same 2^3, same 36 pots, split
# over 1, 2 or 3 laboratories, against the literature alone.
_F3 = ["sugar", "sulfur", "ph"]
PRED_DESIGNS: Dict[str, Dict[str, Any]] = {
    "no study (the literature alone)": {"pots": 0, "labs": []},
    "one laboratory, 2^3 x 4": {"pots": 36, "labs": [_factorial(_F3, 4)]},
    "two laboratories, 2^3 x 2 each": {"pots": 36, "labs": [_factorial(_F3, 2)] * 2},
    "three laboratories, 2^3 x 1 each": {"pots": 36, "labs": [_factorial(_F3, 1)] * 3},
}
#: (matrix shift SD, protocol factor): the design laboratories share our protocol, isolate lot and matrix, so their
#: offsets may scatter less than the literature's (factor 0.5); the new kitchen always scatters with the full tau
PRED_CONFIGS = ((0.3, 1.0), (0.3, 0.5))
_TAU_FLOOR = 1e-6


def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _prior_prediction() -> Dict[str, Any]:
    """No keyed cell and no new data: the prior's prediction (mu ~ N(0, 0.5), tau ~ |N(0, 0.5)|)."""
    sd = math.sqrt(PRIOR_SD ** 2 + PRIOR_SD ** 2)        # Var(mu) + E[tau^2]; E[tau^2] of |N(0, s)| is s^2
    return {"estimate": 0.0, "new_lab_pi90": [-1.645 * sd, 1.645 * sd], "p_same_sign_new_lab": 0.5}


def _score(pred: Mapping[str, Any], m: float, tau: float) -> Dict[str, float]:
    """Expected scores of a prediction against a new kitchen's effect y ~ N(m, tau^2)."""
    tau = max(tau, _TAU_FLOOR)
    e = pred["estimate"]
    lo, hi = pred["new_lab_pi90"]
    sign_err = 0.5 if e == 0.0 else _phi(-math.copysign(1.0, e) * m / tau)
    return {"sign_error": sign_err, "oracle_sign_error": _phi(-abs(m) / tau), "abs_error": abs(e - m),
            "coverage90": _phi((hi - m) / tau) - _phi((lo - m) / tau), "width90": hi - lo,
            "confident": float(pred["p_same_sign_new_lab"] >= LT.P_SAME_SIGN),
            "confident_wrong": float(pred["p_same_sign_new_lab"] >= LT.P_SAME_SIGN) * sign_err}


def run_predictive(job) -> List[Dict[str, Any]]:
    design_name, matrix_sd, tau_factor, seed = job
    cells, base = _STATE["cells"], _STATE["base"]
    targets = _targets(PRED_DESIGNS["three laboratories, 2^3 x 1 each"])
    rng_truth = np.random.default_rng(seed)               # common random numbers: the same worlds for every design
    truths = {}
    for t in targets:
        compound, lc, band, label = t
        p = _truth_params(base.get((compound, lc, band), {}), label)
        if p is None:
            truths[t] = float(rng_truth.normal(0.0, PRIOR_SD))
            truths[("tau",) + t] = abs(float(rng_truth.normal(0.0, PRIOR_SD)))
        else:
            est, se, tau = p
            truths[t] = float(rng_truth.normal(est, se))
            truths[("tau",) + t] = tau
        truths[("shift",) + t[:3]] = float(rng_truth.normal(0.0, matrix_sd)) if matrix_sd else 0.0
    design = PRED_DESIGNS[design_name]
    rng = np.random.default_rng([seed, len(design["labs"])])
    synth = _synthetic(design, truths, matrix_sd, "main_effect", rng, tau_factor=tau_factor)
    by_cell = defaultdict(list)
    for o in synth:
        by_cell[_key(o)].append(o)
    fitted = {}
    for key, new_rows in by_cell.items():
        fitted[key] = readings(cells.get(key, []) + new_rows, key[1])
    out = []
    for t in targets:
        cell, label = t[:3], t[3]
        pred = (fitted.get(cell) if cell in fitted else base.get(cell, {})).get(label) or _prior_prediction()
        m = truths[t] + truths[("shift",) + cell]
        row = _score(pred, m, truths[("tau",) + t])
        row["literature"] = bool((base.get(cell) or {}).get(label))
        out.append(row)
    return out


def summarise_predictive(rows_by_sim: Sequence[List[Dict[str, Any]]]) -> Dict[str, Any]:
    out = {}
    for part, keep in (("all", lambda r: True), ("with a literature estimate", lambda r: r["literature"]),
                       ("prior only", lambda r: not r["literature"])):
        rows = [r for sim in rows_by_sim for r in sim if keep(r)]
        if not rows:
            continue
        n_conf = sum(r["confident"] for r in rows)
        out[part] = {"cells": len(rows) // max(len(rows_by_sim), 1),
                     **{k: float(np.mean([r[k] for r in rows]))
                        for k in ("sign_error", "oracle_sign_error", "abs_error", "coverage90", "width90")},
                     "confident_per_run": n_conf / len(rows_by_sim),
                     "confident_sign_error": (sum(r["confident_wrong"] for r in rows) / n_conf) if n_conf else None}
    return out


def render_predictive(pr: Mapping[str, Any]) -> List[str]:
    L = ["## Does the study improve the prediction for a new kitchen? (predictive check)", "",
         "*Not the 'answered' count (that counts our own three-laboratory rule). For each simulated world, the table's "
         "prediction for one new kitchen cooking in pea isolate (its effect: mu + the matrix shift + a new laboratory "
         "offset drawn with the full tau) is scored: the chance its sign is wrong, the error of its estimate against "
         "the pea-matrix mean, and its 90 % new-laboratory interval's coverage and width. 'Oracle' is the sign error "
         "of a predictor that knows the pea-matrix mean exactly (the floor set by between-laboratory scatter). "
         "'Confident' is P(same sign in a new laboratory) >= 0.8 (the PREREG sec. 5 sign criterion). Protocol factor "
         "0.5: the design laboratories, sharing our protocol and lot, scatter half as much as the literature's, so "
         "the fit sees less between-laboratory spread than a new kitchen has. Same 36 pots; same simulated worlds "
         f"for every design ({pr['sims']} per setting).*", "",
         "| cells | protocol factor | design | sign error | oracle | abs. error of estimate (dex) | 90 % coverage | "
         "90 % width (dex) | confident calls per run | sign error of confident calls |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for part in ("with a literature estimate", "prior only", "all"):
        for r in pr["results"]:
            s = r["scores"].get(part)
            if not s:
                continue
            cse = "-" if s["confident_sign_error"] is None else f"{s['confident_sign_error']:.3f}"
            L.append(f"| {part} ({s['cells']}) | {r['protocol_factor']} | {r['design']} | {s['sign_error']:.3f} | "
                     f"{s['oracle_sign_error']:.3f} | {s['abs_error']:.3f} | {s['coverage90']:.3f} | "
                     f"{s['width90']:.2f} | {s['confident_per_run']:.1f} | {cse} |")
    L.append("")
    return L


def main_predictive(sims: int) -> int:
    _load_state()
    jobs = [(d, m, f, 7_000_000 + 100_000 * ci + s) for ci, (m, f) in enumerate(PRED_CONFIGS)
            for d in PRED_DESIGNS for s in range(sims)]
    res = parallel.map_ordered(run_predictive, jobs)
    grouped = defaultdict(list)
    for (d, m, f, _s), r in zip(jobs, res):
        grouped[(d, m, f)].append(r)
    record = data_access.load_json(OUT)
    record["predictive"] = {"sims": sims, "configs": [list(c) for c in PRED_CONFIGS],
                            "results": [{"design": d, "pots": PRED_DESIGNS[d]["pots"],
                                         "new_labs": len(PRED_DESIGNS[d]["labs"]), "matrix_shift_sd": m,
                                         "protocol_factor": f, "scores": summarise_predictive(grouped[(d, m, f)])}
                                        for m, f in PRED_CONFIGS for d in PRED_DESIGNS]}
    artifact_io.write_artifact(record, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)}")
    print("\n".join(render_predictive(record["predictive"])))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--sims", type=int, default=N_SIMS)
    ap.add_argument("--predictive", action="store_true",
                    help="add the predictive check (new-kitchen prediction error) to the existing record")
    args = ap.parse_args(argv)
    if args.predictive:
        return main_predictive(args.sims)
    _load_state()
    jobs = [(d, m, k, 1000 * di + 10_000 * ci + s) for di, d in enumerate(DESIGNS)
            for ci, (m, k) in enumerate(CONFIGS) for s in range(args.sims)]
    sims = parallel.map_ordered(run_one, jobs)
    grouped = defaultdict(list)
    for (d, m, k, _s), r in zip(jobs, sims):
        grouped[(d, m, k)].append(r)
    results = [summarise(d, m, k, grouped[(d, m, k)]) for d in DESIGNS for m, k in CONFIGS]
    payload = {"artifact": "lever_design_preposterior",
               "settings": {"sims": args.sims, "within_lab_sd": WITHIN_LAB_SD, "sida_cv": SIDA_CV,
                            "other_cv": OTHER_CV, "prior_sd": PRIOR_SD, "analytes": ANALYTES},
               "results": results}
    inputs = sorted(p for p in data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml")
                    if p.name != "TEMPLATE.yml" and ".reader_" not in p.name) + [data_paths.CONTRAST_SPLIT]
    payload["provenance"] = provenance.provenance_block(
        "lever_design_preposterior", generated_by="scripts/generators/build_lever_design_preposterior.py",
        inputs=inputs)
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)}")
    for r in results:
        print(r["design"], r["matrix_shift_sd"], r["keying"], r["newly_answered"],
              r["newly_answered_with_literature"], r["map_cells_newly_answered"], r["step5_cells_two_test_labs"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
