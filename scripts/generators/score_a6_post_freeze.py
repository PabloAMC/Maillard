#!/usr/bin/env python
"""
Score PREREG amendment A6's candidate rule R1 (pH banded by keyed control) on a batch of sources keyed after the A6
freeze (STATUS_AND_PLAN sec. 0, step 10). Protocol: docs/reference/lever_table/A6_SCORING_PREREG.md (revision 2); a
rule here that differs from that file is a bug here. RUN ONCE PER BATCH: each output is a frozen record, pinned by its
commit in the freshness gate.

Eligibility is decided from row counts BEFORE any fit; an ineligible batch gets a record with no prediction in it.

Writes results/validation/a6_post_freeze_score_batch<n>.{json,md}.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import build_lever_field_readings as FR  # noqa: E402
from scripts.generators import build_lever_table as LT  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

PREREG = data_paths.LEVER_TABLE_DIR / "A6_SCORING_PREREG.md"
#: post-freeze batches (new_labs.yml): sources keyed after A6 (3b2f80ad). Every batch's sources leave the training set.
BATCHES = {1: ("yang2011", "raza2020", "madruga1995", "madruga1998", "whitfield1999", "whitfield2001", "chen2000",
               "kang2026", "wang2026arp", "desousafontes2024")}
ALPHA, Z90 = 0.10, 1.6448536269514722
MIN_ROWS, MIN_LABS = 2, 2
MIN_SCORED, MIN_SOURCES, MIN_LABS_SCORED, MAX_SOURCE_SHARE = 20, 4, 3, 0.40
N_DRAWS, N_BOOT, SEED, MIN_LABELLINGS = 1000, 1000, 20261113, 20
COVERAGE_FLOOR = 0.80


def out_path(batch: int) -> Path:
    return data_paths.VALIDATION_DIR / f"a6_post_freeze_score_batch{batch}.json"


# ------------------------------------------------------------------------------------------------- pieces
def interval_score(lo: float, hi: float, y: float, alpha: float = ALPHA) -> float:
    """Gneiting & Raftery (2007) interval score of a central (1 - alpha) interval; lower is better."""
    return (hi - lo) + (2.0 / alpha) * max(lo - y, 0.0) + (2.0 / alpha) * max(y - hi, 0.0)


def control_class(c: Optional[str]) -> Optional[str]:
    return "buffered" if c == "buffered" else ("unbuffered" if c in FR.UNBUFFERED else None)


def fittable(rows: Sequence[Mapping[str, Any]]) -> bool:
    return len(rows) >= MIN_ROWS and len({r["lab"] for r in rows}) >= MIN_LABS


def default_fit(rows):
    s = LT.numeric_fit(list(rows), LT.PRIOR_SCALE)
    return float(s["new_lab_pi90"][0]), float(s["new_lab_pi90"][1])


def score_item(pi: Tuple[float, float], y: float, v: float) -> Dict[str, float]:
    h = Z90 * math.sqrt(max(v, 0.0))
    lo, hi = pi[0] - h, pi[1] + h
    return {"is": interval_score(lo, hi, y), "covered": float(lo <= y <= hi), "width": hi - lo}


def source_weighted(items: Sequence[Mapping[str, Any]], arm: str) -> float:
    by = defaultdict(list)
    for i in items:
        by[i["source"]].append(i[arm]["is"])
    return float(np.mean([np.mean(v) for v in by.values()]))


# ------------------------------------------------------------------------------------------------- eligibility
def candidates(train: Sequence[Dict[str, Any]], test: Sequence[Dict[str, Any]]):
    """Per test contrast: its cell's training rows without its own laboratory, and the R1 subset; kept only where
    both arms can fit and R1's rows differ from the cell's (A6_SCORING_PREREG sec. 3). Row counts only, no fit."""
    cells = defaultdict(list)
    for r in train:
        cells[(r["compound"], r["lever_class"], r["band"])].append(r)
    out, dropped_same = [], 0
    for t in test:
        key = (t["compound"], t["lever_class"], t["band"])
        cls = control_class(t.get("ph_control"))
        rows = [r for r in cells.get(key, []) if r["lab"] != t["lab"]]
        if cls is None or not fittable(rows):
            continue
        sub = [r for r in rows if control_class(r.get("ph_control")) == cls]
        if not fittable(sub):
            continue
        if len(sub) == len(rows):
            dropped_same += 1
            continue
        out.append({"t": t, "key": key, "cls": cls})
    return out, cells, dropped_same


def eligibility(cands) -> Dict[str, Any]:
    n = len(cands)
    srcs = defaultdict(int)
    for c in cands:
        srcs[c["t"]["source"]] += 1
    labs = {c["t"]["lab"] for c in cands}
    share = max(srcs.values()) / n if n else 0.0
    ok = n >= MIN_SCORED and len(srcs) >= MIN_SOURCES and len(labs) >= MIN_LABS_SCORED and share <= MAX_SOURCE_SHARE
    return {"eligible": ok, "contrasts": n, "sources": dict(sorted(srcs.items())), "labs": len(labs),
            "largest_source_share": share,
            "rule": f">= {MIN_SCORED} contrasts, >= {MIN_SOURCES} sources, >= {MIN_LABS_SCORED} laboratories, "
                    f"no source above {MAX_SOURCE_SHARE:.0%}"}


# ------------------------------------------------------------------------------------------------- placebo
def labellings(cells, keys) -> Tuple[int, List[List[Tuple[int, ...]]]]:
    """Per scored cell, the positions of the classed training rows and how many are buffered; the number of distinct
    joint labellings (product of binomials)."""
    total, spec = 1, []
    for k in keys:
        cls = [control_class(r.get("ph_control")) for r in cells[k]]
        pos = [i for i, c in enumerate(cls) if c]
        nb = sum(1 for c in cls if c == "buffered")
        total *= math.comb(len(pos), nb)
        spec.append((pos, nb))
    return total, spec


def placebo_draws(cells, keys, rng):
    total, spec = labellings(cells, keys)
    per_cell = []
    for pos, nb in spec:
        per_cell.append([set(c) for c in itertools.combinations(pos, nb)])
    if total <= N_DRAWS:
        for combo in itertools.product(*per_cell):
            yield total, combo
    else:
        for _ in range(N_DRAWS):
            yield total, tuple(set(rng.choice(pos, size=nb, replace=False).tolist()) if pos else set()
                               for pos, nb in spec)


# ------------------------------------------------------------------------------------------------- the score
def score(train, test, fit: Callable = default_fit, seed: int = SEED) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    cands, cells, dropped_same = candidates(train, test)
    elig = eligibility(cands)
    out: Dict[str, Any] = {"eligibility": elig, "dropped_cells_where_R1_equals_baseline": dropped_same}
    if not elig["eligible"]:
        out["verdict"] = "too few to score"
        return out
    cache: Dict[Any, Tuple[float, float]] = {}

    def pi(rows):
        k = tuple(sorted(id(r) for r in rows))
        if k not in cache:
            cache[k] = fit(rows)
        return cache[k]

    items = []
    for c in cands:
        t = c["t"]
        rows = [r for r in cells[c["key"]] if r["lab"] != t["lab"]]
        sub = [r for r in rows if control_class(r.get("ph_control")) == c["cls"]]
        items.append({"source": t["source"], "contrast": t["contrast"], "lab": t["lab"], "cell": list(c["key"]),
                      "control": c["cls"], "y": t["y"], "baseline": score_item(pi(rows), t["y"], t["v"]),
                      "rule": score_item(pi(sub), t["y"], t["v"])})
    base, rule = source_weighted(items, "baseline"), source_weighted(items, "rule")
    diff = rule - base
    keys = sorted({tuple(i["cell"]) for i in items})
    plac, total = [], None
    for total, combo in placebo_draws(cells, keys, rng):
        buffered = {k: s for k, s in zip(keys, combo)}
        draw = []
        for c in cands:
            t, k = c["t"], c["key"]
            rows = [r for r in cells[k] if r["lab"] != t["lab"]]
            lab = {id(r): ("buffered" if i in buffered[k] else "unbuffered")
                   for i, r in enumerate(cells[k]) if control_class(r.get("ph_control"))}
            sub = [r for r in rows if lab.get(id(r)) == c["cls"]]
            draw.append({"source": t["source"], "rule": score_item(pi(sub), t["y"], t["v"]) if fittable(sub)
                         else score_item(pi(rows), t["y"], t["v"])})
        plac.append(source_weighted(draw, "rule") - base)
    by_src = defaultdict(list)
    for i in items:
        by_src[i["source"]].append(i["rule"]["is"] - i["baseline"]["is"])
    srcs = sorted(by_src)
    boots = []
    for _ in range(N_BOOT):
        pick = rng.choice(len(srcs), size=len(srcs), replace=True)
        boots.append(float(np.mean([np.mean(by_src[srcs[j]]) for j in pick])))
    cov = float(np.mean([i["rule"]["covered"] for i in items]))
    out.update({
        "items": items, "labellings": int(total),
        "baseline": {"source_weighted_mean_is": base, "mean_is": float(np.mean([i["baseline"]["is"] for i in items])),
                     "coverage": float(np.mean([i["baseline"]["covered"] for i in items])),
                     "median_width": float(np.median([i["baseline"]["width"] for i in items]))},
        "rule": {"source_weighted_mean_is": rule, "mean_is": float(np.mean([i["rule"]["is"] for i in items])),
                 "coverage": cov, "median_width": float(np.median([i["rule"]["width"] for i in items]))},
        "diff": diff, "boot90_over_sources": [float(np.quantile(boots, 0.05)), float(np.quantile(boots, 0.95))],
        "placebo_p05": float(np.quantile(plac, 0.05)), "placebo_median": float(np.median(plac)),
        "placebo_share_at_or_below": float(np.mean([p <= diff for p in plac])),
        "per_source": {s: {"n": len(v), "mean_is_diff": float(np.mean(v))} for s, v in sorted(by_src.items())}})
    if total < MIN_LABELLINGS:
        out["verdict"] = f"placebo uninformative ({total} distinct labellings < {MIN_LABELLINGS})"
    elif diff < 0 and diff < out["placebo_p05"] and cov >= COVERAGE_FLOOR:
        out["verdict"] = "R1 supported"
    else:
        out["verdict"] = "R1 not supported" + (" (worse than baseline)" if diff > 0 else "")
    return out


# ------------------------------------------------------------------------------------------------- data, output
def split_observations(batch: int):
    sources = LT.load_sources()
    labs, _ = LT.lab_partition()
    obs, _ = FR.observations_with_contrasts(sources, labs)
    FR.annotate(obs, FR.load_checks())
    for o in obs:
        o.pop("_contrast", None)
    post_all = {s for b in BATCHES.values() for s in b}
    train = [o for o in obs if o["source"] not in post_all]
    test = [o for o in obs if o["source"] in set(BATCHES[batch]) and o["lever_class"] == "pH"]
    return train, test


def _f(x):
    return "-" if x is None else f"{x:+.2f}"


def render(p: Mapping[str, Any]) -> str:
    e = p["eligibility"]
    L = [f"# PREREG A6 (R1) scored on post-freeze batch {p['batch']} (STATUS_AND_PLAN sec. 0, step 10)", "",
         f"*Run once under [A6_SCORING_PREREG.md](../../docs/reference/lever_table/A6_SCORING_PREREG.md), revision 2. "
         f"Training: {p['n_train']} observations from sources in no post-freeze batch; test: {p['n_test']} pH "
         f"observations from batch {p['batch']}'s sources.*", "", f"## Verdict: {p['verdict']}", "",
         f"Eligibility ({e['rule']}): {e['contrasts']} contrasts, {len(e['sources'])} sources, {e['labs']} "
         f"laboratories, largest source share {e['largest_source_share']:.0%}. By source: "
         + (", ".join(f"{s} {n}" for s, n in e["sources"].items()) or "none") + ".",
         f"Contrasts dropped because every training row of their cell has one control class (R1 = baseline): "
         f"{p['dropped_cells_where_R1_equals_baseline']}.", ""]
    if "items" in p:
        b, r = p["baseline"], p["rule"]
        L += ["| arm | source-weighted mean interval score | contrast mean | coverage | median width |",
              "|---|---|---|---|---|",
              f"| baseline | {b['source_weighted_mean_is']:.2f} | {b['mean_is']:.2f} | {b['coverage']:.2f} | "
              f"{b['median_width']:.2f} |",
              f"| R1 | {r['source_weighted_mean_is']:.2f} | {r['mean_is']:.2f} | {r['coverage']:.2f} | "
              f"{r['median_width']:.2f} |", "",
              f"Difference (R1 - baseline) {_f(p['diff'])}; 90 % over sources {_f(p['boot90_over_sources'][0])} to "
              f"{_f(p['boot90_over_sources'][1])}; placebo ({p['labellings']} distinct labellings): 5th percentile "
              f"{_f(p['placebo_p05'])}, median {_f(p['placebo_median'])}.", ""]
    else:
        L += ["No prediction was computed: the protocol decides eligibility from row counts first.", ""]
    L += ["R2 (final pH) is not scored as a predictive rule, and R3 cannot be scored on this batch "
          "(A6_SCORING_PREREG sec. 1). A verdict changes nothing in the table."]
    return "\n".join(L) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--batch", type=int, required=True, choices=sorted(BATCHES))
    args = ap.parse_args(argv)
    out = out_path(args.batch)
    if out.exists():
        raise SystemExit(f"{out.relative_to(ROOT)} exists: each batch is scored once (A6_SCORING_PREREG.md)")
    train, test = split_observations(args.batch)
    payload = {"artifact": out.stem, "batch": args.batch, "prereg": str(PREREG.relative_to(ROOT)),
               "post_freeze_sources": list(BATCHES[args.batch]), "n_train": len(train), "n_test": len(test),
               **score(train, test)}
    inputs = ([p for p in sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml"))
               if p.name != "TEMPLATE.yml" and ".reader_" not in p.name]
              + [data_paths.CONTRAST_SPLIT, data_paths.LEVER_TABLE_NEW_LABS, PREREG])
    payload["provenance"] = provenance.provenance_block(
        out.stem, generated_by="scripts/generators/score_a6_post_freeze.py", inputs=inputs)
    artifact_io.write_artifact(payload, out, render=render)
    print(f"wrote {out.relative_to(ROOT)}: {payload['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
