#!/usr/bin/env python
"""Why do laboratories disagree? An audit of the lever table's discordant cells (STATUS_AND_PLAN sec. 0, step 5b).

A DIAGNOSTIC on the keyed within-study contrasts; literature only, no engine output. Exploratory: the rules below
were fixed in this file before its first run, but after the all-data lever table had been read, so nothing here is
a pre-registered test. Three parts.

(a) Every cell with contrasts from two or more laboratories is classified, first match:
    1. "conflict between laboratories": two laboratories' mean effects have opposite signs, each more than two
       standard errors from zero;
    2. "sign changes inside one laboratory": one laboratory has contrasts of both signs, each more than two standard
       deviations from zero (the effect depends on that study's other settings);
    3. "direction agreed": every laboratory's mean effect has the same sign;
    4. "no clear effect": mixed signs, none separable from zero.
    A laboratory's mean is the inverse-variance mean of its contrasts, each variance floored by POT_SD^2.
(b) Precursor-dose cells pooled over every precursor (PREREG sec. 4 as first written) beside the cells banded by the
    kind of precursor varied (amendment A5).
(c) Which recorded conditions go with disagreement. For every pair of contrasts in one cell from two different
    laboratories, the discordance is |y_i - y_j| (dex per lever step). For each recorded condition the pairs that
    MATCH on it are compared with the pairs that do not, inside each cell, and the per-cell differences are averaged
    (weight: the smaller of the two pair counts). Interval: a bootstrap over laboratories. Also: pairs from one
    laboratory against pairs from two (what a laboratory itself is worth), and a one-way variance split.

Writes results/validation/lever_disagreement_audit.{json,md}.
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import build_lever_table as LT  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "lever_disagreement_audit.json"

POT_SD = 0.10          # dex; floor on a contrast's standard deviation (unmodelled pot-to-pot noise)
Z = 2.0                # standard errors for "separable from zero"
N_BOOT = 1000
SEED = 20261112
MODEL_LEAD = ("buffer", "water")                    # the keyed matrix field's leading category
FOOD_WORDS = ("protein", "meat", "hydrolys", "extract", "flour", "mealworm", "chicken", "pork", "beef", "yeast",
              "oil", "lard", "fish", "seed", "mushroom", "milk", "egg", "soy", "pea ")
SIZE_MIN = 0.05      # dex; smaller effects have no meaningful size ratio
CLASSES = ("conflict between laboratories", "sign changes inside one laboratory", "direction agreed",
           "no clear effect")


# ------------------------------------------------------------------------------------------------- the items
def _median(xs):
    xs = [x for x in xs if x is not None]
    return float(np.median(xs)) if xs else None


def moderators(c: Mapping[str, Any]) -> Dict[str, Any]:
    arms = c.get("arms") or []
    cond = [(a.get("conditions") or {}) for a in arms]
    matrices = {str(k.get("matrix") or "unstated").lower() for k in cond}
    stated = matrices - {"unstated"}
    model = bool(stated) and all(m.split(" (")[0].strip().startswith(MODEL_LEAD) or
                                 (m.startswith("other") and not any(w in m for w in FOOD_WORDS)) for m in stated)
    buffers = {str(k.get("buffer") or "unstated").lower() for k in cond}
    basis = {str(k.get("ph_basis") or "unstated").lower() for k in cond}
    return {
        "quantification": "absolute" if c.get("quant_type") == "absolute" else "peak area or relative",
        "matrix": None if not stated else ("model system" if model else "food or protein matrix"),
        "buffered": None if all(b in ("unstated", "") or b.startswith(("unstated", "none stated")) for b in buffers)
        else not all(b.startswith(("none", "unstated", "water")) or b == "" for b in buffers),
        "ph_basis": basis.pop() if len(basis) == 1 and basis != {"unstated"} else None,
        "digitised": bool(c.get("digitised")),
        "temp_C": _median([LT._num(k.get("temp_C")) for k in cond]),
        "time_min": _median([LT._num(k.get("time_min")) for k in cond]),
        "ph": _median([LT._num(k.get("ph")) for k in cond]),
    }


def items(sources, labs, pooled_precursor: bool = False) -> List[Dict[str, Any]]:
    """One item per numeric contrast, and one per categorical pair (levels in name order, sign flipped to match)."""
    out = []
    for doc in sources:
        for c in doc["contrasts"] or []:
            obs, _ = LT.observations([{"source": doc["source"], "contrasts": [c]}], labs)
            if not obs:
                continue
            o = obs[0]
            base = {"source": o["source"], "lab": o["lab"], "compound": o["compound"],
                    "lever_class": o["lever_class"], **moderators(c)}
            if o["kind"] == "numeric":
                band = "all" if (pooled_precursor and o["lever_class"] == "precursor dose") else o["band"]
                out.append({**base, "band": band, "y": float(o["y"]), "v": float(o["v"]),
                            "lever_mid": float(np.mean(o["range"]))})
            else:
                for p in o["pairs"]:
                    a, b = sorted((p["level"], p["against"]))
                    y = p["y"] if (p["against"], p["level"]) == (a, b) else -p["y"]
                    out.append({**base, "band": f"{b} vs {a}", "y": float(y), "v": float(p["v"]), "lever_mid": None})
    return out


def cells_of(its: Sequence[Mapping[str, Any]], min_labs: int = 2) -> Dict[Tuple[str, str, str], List[Mapping]]:
    g = defaultdict(list)
    for it in its:
        g[(it["compound"], it["lever_class"], it["band"])].append(it)
    return {k: v for k, v in sorted(g.items()) if len({r["lab"] for r in v}) >= min_labs}


# ------------------------------------------------------------------------------------------------- (a) classes
def lab_means(rows: Sequence[Mapping[str, Any]]) -> Dict[str, Tuple[float, float, int]]:
    by = defaultdict(list)
    for r in rows:
        by[r["lab"]].append(r)
    out = {}
    for lab, rs in sorted(by.items()):
        w = np.array([1.0 / (r["v"] + POT_SD ** 2) for r in rs])
        y = np.array([r["y"] for r in rs])
        out[lab] = (float(np.sum(w * y) / np.sum(w)), float(math.sqrt(1.0 / np.sum(w))), len(rs))
    return out


def classify(rows: Sequence[Mapping[str, Any]]) -> str:
    means = lab_means(rows)
    pos = [l for l, (m, se, _) in means.items() if m > Z * se]
    neg = [l for l, (m, se, _) in means.items() if m < -Z * se]
    if pos and neg:
        return CLASSES[0]
    by = defaultdict(list)
    for r in rows:
        by[r["lab"]].append(r["y"] / math.sqrt(r["v"] + POT_SD ** 2))
    if any(max(z) > Z and min(z) < -Z for z in by.values()):
        return CLASSES[1]
    signs = {np.sign(m) for m, _, _ in means.values()}
    if len(signs) == 1 and 0.0 not in signs:
        return CLASSES[2]
    return CLASSES[3]


def classify_cells(cells) -> List[Dict[str, Any]]:
    out = []
    for (compound, lc, band), rows in cells.items():
        means = lab_means(rows)
        out.append({"compound": compound, "lever_class": lc, "band": band, "labs": len(means),
                    "contrasts": len(rows), "class": classify(rows),
                    "lab_means": {l: {"mean": round(m, 3), "se": round(se, 3), "contrasts": n}
                                  for l, (m, se, n) in means.items()}})
    return out


# ------------------------------------------------------------------------------------------------- (c) pairs
def _match(key: str, a: Mapping, b: Mapping) -> Optional[bool]:
    """True / False: the pair matches / differs on this condition; None: not recorded for both."""
    if key in ("quantification", "matrix", "buffered", "ph_basis", "digitised"):
        if a[key] is None or b[key] is None:
            return None
        return a[key] == b[key]
    field = {"temperature within 10 C": "temp_C", "time within 2x": "time_min", "pH within 0.5": "ph",
             "lever midpoint within 1 step": "lever_mid"}[key]
    x, y = a[field], b[field]
    if x is None or y is None:
        return None
    if key == "temperature within 10 C":
        return abs(x - y) <= 10.0
    if key == "time within 2x":
        return x > 0 and y > 0 and abs(math.log2(x / y)) <= 1.0
    if key == "pH within 0.5":
        return abs(x - y) <= 0.5
    step = {"pH": 1.0, "temperature": 10.0}.get(a["lever_class"])
    if step is not None:
        return abs(x - y) <= step
    return x > 0 and y > 0 and abs(math.log2(x / y)) <= 1.0


CONDITIONS = ("quantification", "matrix", "buffered", "ph_basis", "digitised", "temperature within 10 C",
              "time within 2x", "pH within 0.5", "lever midpoint within 1 step")


def stratified_gap(cells, split, dist=None) -> Optional[Dict[str, Any]]:
    """Mean over cells of (mean discordance of group False) - (mean discordance of group True), weighted by the
    smaller pair count. `split(a, b)` returns True / False / None."""
    num = den = 0.0
    n_cells = n_true = n_false = 0
    d_true, d_false = [], []
    for rows in cells.values():
        t, f = [], []
        for a, b in itertools.combinations(rows, 2):
            s = split(a, b)
            if s is None:
                continue
            dd = abs(a["y"] - b["y"]) if dist is None else dist(a, b)
            if dd is None:
                continue
            (t if s else f).append(dd)
        if t and f:
            w = min(len(t), len(f))
            num += w * (float(np.mean(f)) - float(np.mean(t)))
            den += w
            n_cells += 1
            n_true += len(t)
            n_false += len(f)
            d_true += t
            d_false += f
    if den == 0:
        return None
    return {"gap": num / den, "cells": n_cells, "pairs_true": n_true, "pairs_false": n_false,
            "median_true": float(np.median(d_true)), "median_false": float(np.median(d_false))}


def bootstrap_gap(its: Sequence[Mapping[str, Any]], split, rng, dist=None, keep=None) -> Optional[Dict[str, Any]]:
    cells = {k: v for k, v in cells_of(its).items() if keep is None or k in keep}
    point = stratified_gap(cells, split, dist)
    if point is None:
        return None
    labs = sorted({it["lab"] for it in its})
    by_lab = defaultdict(list)
    for it in its:
        by_lab[it["lab"]].append(it)
    draws = []
    for _ in range(N_BOOT):
        pick = rng.choice(len(labs), size=len(labs), replace=True)
        sample = []
        for k, i in enumerate(pick):               # a resampled laboratory is a new laboratory (its own id)
            sample += [{**it, "lab": f"{labs[i]}#{k}", "_orig": labs[i]} for it in by_lab[labs[i]]]
        g = stratified_gap({k: v for k, v in cells_of(sample).items() if keep is None or k in keep},
                           lambda a, b: None if a["_orig"] == b["_orig"] and a["lab"] != b["lab"] else split(a, b),
                           dist)
        if g is not None:
            draws.append(g["gap"])
    if len(draws) >= 0.5 * N_BOOT:
        lo, hi = np.quantile(draws, [0.05, 0.95])
        point.update({"ci90": [float(lo), float(hi)], "p_gap_positive": float(np.mean(np.array(draws) > 0))})
    else:
        point.update({"ci90": None, "p_gap_positive": None})
    return point


def variance_split(cells) -> Dict[str, Any]:
    """One-way random-effects split, pooled over cells: within-laboratory and between-laboratory SD of an effect."""
    ssw = ssb = 0.0
    n_tot = n_labs = n_cells = 0
    n0_num = 0.0
    for rows in cells.values():
        by = defaultdict(list)
        for r in rows:
            by[r["lab"]].append(r["y"])
        grand = float(np.mean([r["y"] for r in rows]))
        for ys in by.values():
            m = float(np.mean(ys))
            ssw += float(np.sum((np.array(ys) - m) ** 2))
            ssb += len(ys) * (m - grand) ** 2
        n = len(rows)
        n_tot += n
        n_labs += len(by)
        n_cells += 1
        n0_num += n - sum(len(ys) ** 2 for ys in by.values()) / n
    dfw, dfb = n_tot - n_labs, n_labs - n_cells
    if dfw <= 0 or dfb <= 0:
        return {"cells": n_cells}
    msw, msb = ssw / dfw, ssb / dfb
    n0 = n0_num / dfb
    var_lab = max(0.0, (msb - msw) / n0)
    return {"cells": n_cells, "contrasts": n_tot, "sd_within_lab": math.sqrt(msw), "sd_between_lab": math.sqrt(var_lab),
            "share_between_lab": var_lab / (var_lab + msw) if (var_lab + msw) > 0 else None}


# ------------------------------------------------------------------------------------------------- build
def build(sources, labs) -> Dict[str, Any]:
    rng = np.random.default_rng(SEED)
    its = items(sources, labs)
    cells = cells_of(its)
    classified = classify_cells(cells)
    pooled = classify_cells({k: v for k, v in cells_of(items(sources, labs, pooled_precursor=True)).items()
                             if k[1] == "precursor dose"})
    banded = [c for c in classified if c["lever_class"] == "precursor dose"]

    def tally(cs, min_labs=2):
        t = Counter(c["class"] for c in cs if c["labs"] >= min_labs)
        return {k: t.get(k, 0) for k in CLASSES}

    by_lever = {}
    for lc in sorted({c["lever_class"] for c in classified}):
        cs = [c for c in classified if c["lever_class"] == lc]
        by_lever[lc] = {"cells": len(cs), **tally(cs)}

    conds = {}
    for key in CONDITIONS:
        sub = its if key != "ph_basis" else [it for it in its if it["lever_class"] == "pH"]
        g = bootstrap_gap(sub, lambda a, b, key=key: None if a["lab"] == b["lab"] else _match(key, a, b), rng)
        if g is not None:
            conds[key] = g
    lab_effect = bootstrap_gap(its, lambda a, b: a["lab"] == b["lab"], rng)

    # (c') the size of an agreed effect: for cells whose laboratories agree in direction, the distance between two
    # same-sign effects is |log10(|y_a| / |y_b|)| (a ratio of sizes), for effects of at least SIZE_MIN dex
    agreed = {(c["compound"], c["lever_class"], c["band"]) for c in classified if c["class"] == CLASSES[2]}

    def size_dist(a, b):
        if a["y"] * b["y"] <= 0 or min(abs(a["y"]), abs(b["y"])) < SIZE_MIN:
            return None
        return abs(math.log10(abs(a["y"]) / abs(b["y"])))

    size = {}
    for key in CONDITIONS:
        g = bootstrap_gap(its, lambda a, b, key=key: None if a["lab"] == b["lab"] else _match(key, a, b), rng,
                          dist=size_dist, keep=agreed)
        if g is not None:
            size[key] = g
    size_lab = bootstrap_gap(its, lambda a, b: a["lab"] == b["lab"], rng, dist=size_dist, keep=agreed)

    numeric_cells = {k: v for k, v in cells.items() if k[1] in LT.NUMERIC}
    var = {"all cells": variance_split(cells), "numeric levers": variance_split(numeric_cells),
           "categorical pairs": variance_split({k: v for k, v in cells.items() if k[1] not in LT.NUMERIC})}
    for lc in sorted({k[1] for k in cells}):
        var[lc] = variance_split({k: v for k, v in cells.items() if k[1] == lc})

    return {"artifact": "lever_disagreement_audit",
            "status": "DIAGNOSTIC, exploratory; literature contrasts only; no engine output",
            "rules": {"pot_sd_floor_dex": POT_SD, "z": Z, "bootstrap": N_BOOT, "seed": SEED, "classes": list(CLASSES)},
            "counts": {"items": len(its), "laboratories": len({it["lab"] for it in its}),
                       "cells_two_or_more_labs": len(cells),
                       "cells_three_or_more_labs": sum(1 for c in classified if c["labs"] >= 3)},
            "classes_two_or_more_labs": tally(classified), "classes_three_or_more_labs": tally(classified, 3),
            "classes_by_lever": by_lever,
            "precursor_dose": {"pooled": {"cells": len(pooled), **tally(pooled)},
                               "banded_A5": {"cells": len(banded), **tally(banded)},
                               "pooled_cells": pooled},
            "conditions": conds, "same_laboratory_vs_two": lab_effect,
            "size_in_agreed_cells": {"conditions": size, "same_laboratory_vs_two": size_lab, "cells": len(agreed)}, "variance_split": var,
            "cells": classified}


def _g(x, nd=2):
    return "-" if x is None else f"{x:+.{nd}f}"


def render(p: Mapping[str, Any]) -> str:
    c2, c3 = p["classes_two_or_more_labs"], p["classes_three_or_more_labs"]
    L = ["# Why laboratories disagree: an audit of the lever table", "",
         "*Generated by `scripts/generators/build_lever_disagreement_audit.py` (STATUS_AND_PLAN sec. 0, step 5b). "
         "A diagnostic on the keyed literature contrasts; no engine output. Exploratory: the rules were fixed in the "
         "generator before its first run, but after the lever table had been read.*", "",
         f"{p['counts']['items']} effects (one per numeric contrast, one per pair of levels of a categorical "
         f"contrast) from {p['counts']['laboratories']} laboratories; {p['counts']['cells_two_or_more_labs']} cells "
         f"hold two or more laboratories, {p['counts']['cells_three_or_more_labs']} hold three or more. A cell here is "
         "an odorant x lever x band, or an odorant x pair of levels. A pair cell 'B vs A' holds log10(B / A).", "",
         "## 1. What kind of disagreement each cell shows", "",
         "| class | cells with 2+ laboratories | cells with 3+ |", "|---|---|---|"]
    for k in CLASSES:
        L.append(f"| {k} | {c2[k]} | {c3[k]} |")
    L += ["", "Rules, first match: **conflict** = two laboratories' mean effects have opposite signs, each more than "
          f"{p['rules']['z']:.0f} standard errors from zero; **sign changes inside one laboratory** = one laboratory "
          "has contrasts of both signs, each separable from zero (the effect depends on that study's other settings); "
          "**direction agreed** = every laboratory's mean has the same sign; **no clear effect** = mixed signs, none "
          f"separable from zero. Each contrast's SD is floored at {p['rules']['pot_sd_floor_dex']} dex.", "",
          "| lever | cells | conflict | sign changes inside a lab | direction agreed | no clear effect |",
          "|---|---|---|---|---|---|"]
    for lc, t in p["classes_by_lever"].items():
        L.append(f"| {lc} | {t['cells']} | {t[CLASSES[0]]} | {t[CLASSES[1]]} | {t[CLASSES[2]]} | {t[CLASSES[3]]} |")
    L += ["", "### The cells in conflict", "",
          "| odorant | lever | band or pair | laboratory means (dex per step; contrasts) |", "|---|---|---|---|"]
    for c in p["cells"]:
        if c["class"] == CLASSES[0]:
            means = "; ".join(f"{l} {_g(m['mean'])} ({m['contrasts']})" for l, m in c["lab_means"].items())
            L.append(f"| {c['compound']} | {c['lever_class']} | {c['band']} | {means} |")
    pd = p["precursor_dose"]
    L += ["", "## 2. Precursor dose: pooled over every precursor, and banded by the precursor varied (amendment A5)",
          "", "| | cells | conflict | sign changes inside a lab | direction agreed | no clear effect |",
          "|---|---|---|---|---|---|"]
    for name, key in (("pooled (as first written)", "pooled"), ("banded (A5)", "banded_A5")):
        t = pd[key]
        L.append(f"| {name} | {t['cells']} | {t[CLASSES[0]]} | {t[CLASSES[1]]} | {t[CLASSES[2]]} | {t[CLASSES[3]]} |")
    L += ["", "## 3. Which recorded conditions go with disagreement", "",
          "For pairs of effects in one cell from two different laboratories: the mean |difference| of the pairs that "
          "DIFFER on a condition minus that of the pairs that MATCH on it, inside each cell, averaged over cells. A "
          "positive gap means that matching on the condition goes with closer agreement. Interval: bootstrap over "
          "laboratories.", "",
          "| condition | gap, dex (90 % interval) | P(gap > 0) | median, matching pairs | median, differing pairs | "
          "cells | pairs (match / differ) |", "|---|---|---|---|---|---|---|"]
    for k, g in p["conditions"].items():
        ci = "-" if not g.get("ci90") else f"{_g(g['ci90'][0])} to {_g(g['ci90'][1])}"
        pg = "-" if g.get("p_gap_positive") is None else f"{g['p_gap_positive']:.2f}"
        L.append(f"| {k} | {_g(g['gap'])} ({ci}) | {pg} | {g['median_true']:.2f} | {g['median_false']:.2f} | "
                 f"{g['cells']} | {g['pairs_true']} / {g['pairs_false']} |")
    g = p["same_laboratory_vs_two"]
    if g:
        ci = "-" if not g.get("ci90") else f"{_g(g['ci90'][0])} to {_g(g['ci90'][1])}"
        L += ["", f"**What a laboratory itself is worth.** Pairs from two laboratories differ by {_g(g['gap'])} dex "
              f"more than pairs from one laboratory in the same cell ({ci}; median {g['median_true']:.2f} within a "
              f"laboratory, {g['median_false']:.2f} between; {g['cells']} cells)."]
    sz = p["size_in_agreed_cells"]
    L += ["", f"### 3b. The SIZE of agreed effects ({sz['cells']} cells whose laboratories agree in direction)", "",
          "Same method, but the distance between two same-sign effects is |log10 of their size ratio| (effects of at "
          f"least {SIZE_MIN} dex). A positive gap: laboratories that match on the condition agree more closely on size.",
          "", "| condition | gap, dex (90 % interval) | P(gap > 0) | median, matching | median, differing | cells | "
          "pairs (match / differ) |", "|---|---|---|---|---|---|---|"]
    for k, g in sz["conditions"].items():
        ci = "-" if not g.get("ci90") else f"{_g(g['ci90'][0])} to {_g(g['ci90'][1])}"
        pg = "-" if g.get("p_gap_positive") is None else f"{g['p_gap_positive']:.2f}"
        L.append(f"| {k} | {_g(g['gap'])} ({ci}) | {pg} | {g['median_true']:.2f} | {g['median_false']:.2f} | "
                 f"{g['cells']} | {g['pairs_true']} / {g['pairs_false']} |")
    g = sz["same_laboratory_vs_two"]
    if g:
        ci = "-" if not g.get("ci90") else f"{_g(g['ci90'][0])} to {_g(g['ci90'][1])}"
        L += ["", f"Two effects from one laboratory differ in size by a median factor {10 ** g['median_true']:.1f}; "
              f"from two laboratories, {10 ** g['median_false']:.1f} (gap {_g(g['gap'])} dex, {ci})."]
    L += ["", "## 4. Variance split (one-way, pooled over cells)", "",
          "| cells | n cells | n effects | SD within a laboratory | SD between laboratories | share between |",
          "|---|---|---|---|---|---|"]
    for k, v in p["variance_split"].items():
        if "sd_within_lab" in v:
            sh = "-" if v["share_between_lab"] is None else f"{100 * v['share_between_lab']:.0f} %"
            L.append(f"| {k} | {v['cells']} | {v['contrasts']} | {v['sd_within_lab']:.2f} | "
                     f"{v['sd_between_lab']:.2f} | {sh} |")
    L += ["", "## What this cannot say", "",
          "- The conditions are those the papers print and the keying recorded. The true pH during heating, the "
          "extraction recovery, free against bound thiol and the headspace are mostly not recorded, so they cannot be "
          "tested here.",
          "- Conditions travel together: a laboratory that uses peak areas also tends to use a food matrix. A gap "
          "is an association inside cells, not a cause.",
          "- Pairs share laboratories, so the bootstrap resamples laboratories; with few laboratories per cell the "
          "intervals are wide and a null result is weak evidence of no effect.", ""]
    return "\n".join(L)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)
    sources = LT.load_sources()
    labs, _ = LT.lab_partition()
    payload = build(sources, labs)
    inputs = sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml")) + [data_paths.CONTRAST_SPLIT,
                                                                           data_paths.LEVER_TABLE_NEW_LABS]
    payload["provenance"] = provenance.provenance_block(
        OUT.stem, generated_by="scripts/generators/build_lever_disagreement_audit.py",
        inputs=[q for q in inputs if q.name != "TEMPLATE.yml" and ".reader_" not in q.name])
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)}: {payload['classes_two_or_more_labs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
