#!/usr/bin/env python
"""
Score PREREG amendment A6's candidate rule R1 prequentially (STATUS_AND_PLAN sec. 0, step 12). Protocol:
docs/reference/lever_table/A6_PREQUENTIAL_PREREG.md (revision 3); a rule here that differs from that file is a bug
here. A contrast from post-freeze batch k is predicted from batches 0 .. k-1 only. Two sets, each scored ONCE:
  --set first_look    batches 1 and 2 pooled
  --set confirmation  batches 3, 4, ... pooled; a counts-only record per batch until the set is eligible, then the
                      one scored record, never again.
Eligibility is decided from row counts BEFORE any fit. The revision-2 scorer (score_a6_post_freeze.py) generates the
batch-1 and batch-2 records and is not edited; its pieces are reused here.

Writes results/validation/a6_prequential_score_first_look.{json,md} or
results/validation/a6_prequential_confirmation_through_batch<k>.{json,md}.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import score_a6_post_freeze as S  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

PREREG = data_paths.LEVER_TABLE_DIR / "A6_PREQUENTIAL_PREREG.md"
#: batches 1 and 2 are revision 2's (score_a6_post_freeze.BATCHES); batch 3 onward are added here, each registered in
#: new_labs.yml before any of its values is read.
BATCHES: Dict[int, Sequence[str]] = {
    **S.BATCHES,
    # batch 3 (new_labs.yml, fifth registration, 2026-11-14): STATUS_AND_PLAN step 15's sources, keyed for the
    # temperature / time cells; they join A6's confirmation set because every post-freeze source does (sec. 2).
    3: ("pastoriza2016", "lee2024oven", "wang2022lipid", "seo2015", "tan2011", "chansataporn2018", "goncuoglu2023")}
FIRST_LOOK = (1, 2)


def confirmation_batches() -> tuple:
    return tuple(sorted(k for k in BATCHES if k >= 3))


def out_path(set_name: str) -> Path:
    if set_name == "first_look":
        return data_paths.VALIDATION_DIR / "a6_prequential_score_first_look.json"
    return data_paths.VALIDATION_DIR / f"a6_prequential_confirmation_through_batch{max(confirmation_batches())}.json"


def batch_map(batches: Mapping[int, Sequence[str]]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for k, srcs in batches.items():
        for s in srcs:
            if s in out:
                raise ValueError(f"{s} is registered in batches {out[s]} and {k}")
            out[s] = k
    return out


# ------------------------------------------------------------------------------------------------- candidates
MAX_LAB_SHARE = 0.40
MIN_USABLE_LABELLINGS = 20
DUPLICATE_TOL = 5e-4


def prequential_candidates(obs: Sequence[Dict[str, Any]], bmap: Mapping[str, int], batches: Sequence[int]):
    """Per batch k of the set: test = its pH observations, training = observations of batches 0 .. k-1 (sec. 2).
    For each test contrast: its cell's training rows without its own laboratory and without any row of another
    source carrying the same slope (a duplicated dataset, sec. 2); kept where both arms can fit and R1's rows differ
    from the cell's. Returns (candidates, dropped because R1 = baseline, duplicate rows removed). Row counts only."""
    pooled, dropped, dups = [], 0, 0
    for k in batches:
        cells = defaultdict(list)
        for o in obs:
            if bmap.get(o["source"], 0) < k:
                cells[(o["compound"], o["lever_class"], o["band"])].append(o)
        for t in obs:
            if bmap.get(t["source"], 0) != k or t["lever_class"] != "pH":
                continue
            key = (t["compound"], t["lever_class"], t["band"])
            cls = S.control_class(t.get("ph_control"))
            if cls is None:
                continue
            same = [r for r in cells.get(key, [])
                    if r["source"] != t["source"] and abs(r["y"] - t["y"]) < DUPLICATE_TOL]
            rows = [r for r in cells.get(key, []) if r["lab"] != t["lab"] and not any(r is d for d in same)]
            sub = [r for r in rows if S.control_class(r.get("ph_control")) == cls]
            if not (S.fittable(rows) and S.fittable(sub)):
                continue
            if len(sub) == len(rows):
                dropped += 1
                continue
            dups += len(same)
            pooled.append({"t": t, "key": key, "cls": cls, "batch": k, "rows": rows})
    return pooled, dropped, dups


def eligibility(cands) -> Dict[str, Any]:
    """Revision 2's thresholds, plus no laboratory above 40 % (sec. 3)."""
    e = S.eligibility(cands)
    labs = defaultdict(int)
    for c in cands:
        labs[c["t"]["lab"]] += 1
    lab_share = max(labs.values()) / len(cands) if cands else 0.0
    e.update({"eligible": e["eligible"] and lab_share <= MAX_LAB_SHARE, "largest_lab_share": lab_share,
              "by_lab": dict(sorted(labs.items())),
              "rule": e["rule"] + f", no laboratory above {MAX_LAB_SHARE:.0%}, "
                                  f">= {MIN_USABLE_LABELLINGS} distinct usable placebo labellings"})
    return e


# ------------------------------------------------------------------------------------------------- placebo
def placebo_strata(cands, bmap):
    """Per (scored cell, batch of origin): the classed training rows any contrast of the set could use, grouped into
    blocks of (source, true class) -- buffering belongs to a source's design, not to a row (sec. 5)."""
    strata = defaultdict(lambda: defaultdict(list))
    seen = set()
    for c in cands:
        for r in c["rows"]:
            cl = S.control_class(r.get("ph_control"))
            if cl is None or id(r) in seen:
                continue
            seen.add(id(r))
            strata[(c["key"], bmap.get(r["source"], 0))][(r["source"], cl)].append(r)
    spec = []
    for skey in sorted(strata, key=str):
        blocks = [strata[skey][b] for b in sorted(strata[skey])]
        nb = sum(1 for b in sorted(strata[skey]) if b[1] == "buffered")
        spec.append((blocks, nb))
    return spec


def placebo_labellings(cands, bmap, rng):
    """Every distinct labelling (at most N_DRAWS of them; otherwise N_DRAWS random ones) and which are usable: a
    labelling is usable when every scored contrast's relabelled R1 subset can fit. Counts only, no fit."""
    spec = placebo_strata(cands, bmap)
    total = 1
    for blocks, nb in spec:
        total *= math.comb(len(blocks), nb)
    if total <= S.N_DRAWS:
        combos = list(itertools.product(*[list(itertools.combinations(range(len(b)), nb)) for b, nb in spec]))
    else:
        combos = [tuple(tuple(sorted(rng.choice(len(b), size=nb, replace=False).tolist())) for b, nb in spec)
                  for _ in range(S.N_DRAWS)]
    usable, discarded = [], 0
    for combo in combos:
        label = {}
        for (blocks, _), chosen in zip(spec, combo):
            for j, blk in enumerate(blocks):
                for r in blk:
                    label[id(r)] = "buffered" if j in chosen else "unbuffered"
        subs = [[r for r in c["rows"] if label.get(id(r)) == c["cls"]] for c in cands]
        if all(S.fittable(s) for s in subs):
            usable.append((combo, subs))
        else:
            discarded += 1
    return total, usable, discarded, len({u[0] for u in usable})


# ------------------------------------------------------------------------------------------------- the score
def weighted(items, value: Callable, by: str) -> float:
    g = defaultdict(list)
    for i in items:
        g[i[by]].append(value(i))
    return float(np.mean([np.mean(v) for v in g.values()]))


def score(obs, bmap, batches, fit: Callable = S.default_fit, seed: int = S.SEED) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    cands, dropped, dups = prequential_candidates(obs, bmap, batches)
    elig = eligibility(cands)
    out: Dict[str, Any] = {"eligibility": elig, "dropped_cells_where_R1_equals_baseline": dropped,
                           "duplicate_training_rows_removed": dups,
                           "by_batch": {str(k): sum(1 for c in cands if c["batch"] == k) for k in batches}}
    if elig["eligible"]:
        total, usable, discarded, n_distinct = placebo_labellings(cands, bmap, rng)
        out["placebo_labellings"] = {"total": int(total), "drawn": len(usable) + discarded, "usable": len(usable),
                                     "discarded_unfittable": discarded, "distinct_usable": n_distinct}
        if n_distinct < MIN_USABLE_LABELLINGS:
            elig["eligible"] = False
    if not elig["eligible"]:
        out["verdict"] = "too few to score"
        return out
    cache: Dict[Any, Any] = {}

    def pi(rows):
        k = tuple(sorted(id(r) for r in rows))
        if k not in cache:
            cache[k] = fit(rows)
        return cache[k]

    items = []
    for c in cands:
        t = c["t"]
        sub = [r for r in c["rows"] if S.control_class(r.get("ph_control")) == c["cls"]]
        items.append({"batch": c["batch"], "source": t["source"], "contrast": t["contrast"], "lab": t["lab"],
                      "cell": list(c["key"]), "control": c["cls"], "y": t["y"],
                      "baseline": S.score_item(pi(c["rows"]), t["y"], t["v"]),
                      "rule": S.score_item(pi(sub), t["y"], t["v"])})

    def d(i):
        return i["rule"]["is"] - i["baseline"]["is"]

    diff = weighted(items, d, "lab")
    plac = []
    for _, subs in usable:
        draw = [{"lab": c["t"]["lab"], "d": S.score_item(pi(s), c["t"]["y"], c["t"]["v"])["is"] - i["baseline"]["is"]}
                for c, s, i in zip(cands, subs, items)]
        plac.append(weighted(draw, lambda x: x["d"], "lab"))
    by_lab = defaultdict(list)
    for i in items:
        by_lab[i["lab"]].append(d(i))
    labs = sorted(by_lab)
    boots = [float(np.mean([np.mean(by_lab[labs[j]]) for j in rng.choice(len(labs), size=len(labs), replace=True)]))
             for _ in range(S.N_BOOT)]

    def arm(a):
        return {"lab_weighted_mean_is": weighted(items, lambda i: i[a]["is"], "lab"),
                "source_weighted_mean_is": weighted(items, lambda i: i[a]["is"], "source"),
                "mean_is": float(np.mean([i[a]["is"] for i in items])),
                "coverage": float(np.mean([i[a]["covered"] for i in items])),
                "median_width": float(np.median([i[a]["width"] for i in items]))}

    cells = sorted({tuple(i["cell"]) for i in items})
    loco = {" | ".join(k): weighted([i for i in items if tuple(i["cell"]) != k], d, "lab") for k in cells
            if any(tuple(i["cell"]) != k for i in items)}
    out.update({
        "items": items, "baseline": arm("baseline"), "rule": arm("rule"), "diff": diff,
        "diff_source_weighted": weighted(items, d, "source"),
        "boot90_over_labs": [float(np.quantile(boots, 0.05)), float(np.quantile(boots, 0.95))],
        "placebo_p05": float(np.quantile(plac, 0.05)), "placebo_median": float(np.median(plac)),
        "placebo_share_at_or_below": float(np.mean([p <= diff for p in plac])),
        "diagnostics": {
            "per_batch": {str(k): weighted([i for i in items if i["batch"] == k], d, "lab")
                          for k in batches if any(i["batch"] == k for i in items)},
            "per_cell": {" | ".join(k): {"n": sum(1 for i in items if tuple(i["cell"]) == k),
                                         "mean_is_diff": float(np.mean([d(i) for i in items if tuple(i["cell"]) == k]))}
                         for k in cells},
            "leave_one_cell_out": loco,
            "leave_one_cell_out_same_sign": sum(1 for v in loco.values() if np.sign(v) == np.sign(diff)),
            "per_lab": {l: {"n": len(v), "mean_is_diff": float(np.mean(v))} for l, v in sorted(by_lab.items())}}})
    cov = out["rule"]["coverage"]
    out["coverage_flag"] = "R1 coverage below 0.80" if cov < S.COVERAGE_FLOOR else "none"
    if diff < 0 and diff < out["placebo_p05"]:
        out["verdict"] = "R1 supported"
    else:
        out["verdict"] = "R1 not supported" + (" (worse than baseline)" if diff > 0 else "")
    return out


# ------------------------------------------------------------------------------------------------- data, output
def load_observations():
    sources = S.LT.load_sources()
    labs, _ = S.LT.lab_partition()
    obs, _ = S.FR.observations_with_contrasts(sources, labs)
    S.FR.annotate(obs, S.FR.load_checks())
    for o in obs:
        o.pop("_contrast", None)
    return obs


def _scored_confirmation_exists() -> Optional[Path]:
    for p in sorted(data_paths.VALIDATION_DIR.glob("a6_prequential_confirmation_through_batch*.json")):
        if "items" in json.loads(p.read_text()):
            return p
    return None


def _f(x):
    return f"{x:+.2f}"


FIRST_LOOK_NOTE = ("Revision 3's first draft proposed this pooled set because a source-only count called it "
                   "eligible; the premortem's laboratory cap was adopted before the count was rerun. No prediction "
                   "had been computed. This look commits nothing: the confirmation set decides "
                   "(A6_PREQUENTIAL_PREREG secs. 0, 5).")


def render(p: Mapping[str, Any]) -> str:
    e = p["eligibility"]
    name = "first look (batches 1 - 2)" if p["set"] == "first_look" else \
        f"confirmation (batches {p['batches'][0]} - {p['batches'][-1]})"
    L = [f"# PREREG A6 (R1) scored prequentially: {name}", "",
         f"*Run once under [A6_PREQUENTIAL_PREREG.md](../../docs/reference/lever_table/A6_PREQUENTIAL_PREREG.md), "
         f"revision 3. Each batch's pH contrasts are predicted from the batches keyed before it; new-laboratory "
         f"predictions only.*", ""]
    if p.get("note"):
        L += [f"**Note.** {p['note']}", ""]
    L += [f"## Verdict: {p['verdict']}", "",
          f"Eligibility ({e['rule']}): {e['contrasts']} contrasts, {len(e['sources'])} sources, {e['labs']} "
          f"laboratories, largest source share {e['largest_source_share']:.0%}, largest laboratory share "
          f"{e['largest_lab_share']:.0%}. By source: "
          + (", ".join(f"{s} {n}" for s, n in e["sources"].items()) or "none") + ". By batch: "
          + ", ".join(f"{k}: {n}" for k, n in p["by_batch"].items()) + ".",
          f"Dropped because R1 = baseline in their cell: {p['dropped_cells_where_R1_equals_baseline']}; duplicate "
          f"training rows removed: {p['duplicate_training_rows_removed']}."]
    if "placebo_labellings" in p:
        pl = p["placebo_labellings"]
        L.append(f"Placebo labellings: {pl['total']} possible, {pl['drawn']} drawn, {pl['usable']} usable "
                 f"({pl['distinct_usable']} distinct), {pl['discarded_unfittable']} discarded as unfittable.")
    L.append("")
    if "items" in p:
        L += ["| arm | laboratory-weighted mean IS | source-weighted | contrast mean | coverage | median width |",
              "|---|---|---|---|---|---|"]
        for a, lab in (("baseline", "baseline"), ("rule", "R1")):
            r = p[a]
            L.append(f"| {lab} | {r['lab_weighted_mean_is']:.2f} | {r['source_weighted_mean_is']:.2f} | "
                     f"{r['mean_is']:.2f} | {r['coverage']:.2f} | {r['median_width']:.2f} |")
        g = p["diagnostics"]
        L += ["", f"Difference (R1 - baseline, laboratory-weighted) {_f(p['diff'])} (source-weighted "
              f"{_f(p['diff_source_weighted'])}); 90 % over laboratories {_f(p['boot90_over_labs'][0])} to "
              f"{_f(p['boot90_over_labs'][1])}; placebo 5th percentile {_f(p['placebo_p05'])}, median "
              f"{_f(p['placebo_median'])}, share at or below R1's {p['placebo_share_at_or_below']:.2f}. "
              f"Coverage flag: {p['coverage_flag']}.", "",
              "Diagnostics (they do not change the verdict): per batch "
              + ", ".join(f"{k} {_f(v)}" for k, v in g["per_batch"].items())
              + f"; leave one cell out keeps the sign in {g['leave_one_cell_out_same_sign']} of "
              f"{len(g['leave_one_cell_out'])}.", "", "| cell | contrasts | mean IS difference |", "|---|---|---|"]
        L += [f"| {k} | {v['n']} | {_f(v['mean_is_diff'])} |" for k, v in g["per_cell"].items()]
        L += ["", "| laboratory | contrasts | mean IS difference |", "|---|---|---|"]
        L += [f"| {k} | {v['n']} | {_f(v['mean_is_diff'])} |" for k, v in g["per_lab"].items()]
        L.append("")
    else:
        L += ["No prediction was computed: the protocol decides eligibility from row counts first.", ""]
    L.append("What a verdict decides is fixed in A6_PREQUENTIAL_PREREG sec. 5; it changes nothing in the table.")
    return "\n".join(L) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--set", required=True, choices=("first_look", "confirmation"))
    args = ap.parse_args(argv)
    batches = FIRST_LOOK if args.set == "first_look" else confirmation_batches()
    if not batches:
        raise SystemExit("no confirmation batch is registered yet (BATCHES has nothing from 3 on)")
    out = out_path(args.set)
    if out.exists():
        raise SystemExit(f"{out.relative_to(ROOT)} exists: each record is run once (A6_PREQUENTIAL_PREREG.md)")
    if args.set == "confirmation" and (done := _scored_confirmation_exists()):
        raise SystemExit(f"{done.relative_to(ROOT)} already scored the confirmation set; it is never reopened")
    obs = load_observations()
    payload = {"artifact": out.stem, "set": args.set, "batches": list(batches), "prereg": str(PREREG.relative_to(ROOT)),
               "post_freeze_sources": {str(k): list(BATCHES[k]) for k in batches},
               **({"note": FIRST_LOOK_NOTE} if args.set == "first_look" else {}),
               **score(obs, batch_map(BATCHES), batches)}
    inputs = ([p for p in sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml"))
               if p.name != "TEMPLATE.yml" and ".reader_" not in p.name]
              + [data_paths.CONTRAST_SPLIT, data_paths.LEVER_TABLE_NEW_LABS, PREREG])
    payload["provenance"] = provenance.provenance_block(
        out.stem, generated_by="scripts/generators/score_a6_prequential.py", inputs=inputs)
    artifact_io.write_artifact(payload, out, render=render)
    print(f"wrote {out.relative_to(ROOT)}: {payload['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
