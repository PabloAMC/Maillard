#!/usr/bin/env python
"""
WAVE B46 -- THE PER-LANE OFFSET DIAGNOSTIC.
Pre-registration: results/validation/kinetic_core_b46_prereg.md (written before this ran).

WHY THIS EXISTS
---------------
`core_prediction_uncertainty.json` reports that widening every parameter prior to uncapped moves
coverage from 19 % to 21 % against a nominal 90 %, and concludes that what remains is "model-structure
error -- systematic per-lane offsets a draw around a wrong centre cannot reproduce". That sentence
names a signature. Nothing measured it. This does.

WHAT IT MEASURES, AND WHY SIGNED
--------------------------------
Every score published on this branch is an UNSIGNED fold error. Unsigned errors cannot tell a lane
that is randomly wrong from a lane that is consistently wrong, and only the second kind points at a
missing process. So this reads the signed offset

    dex = log10(predicted / measured)        positive = the model reads high

and asks, per lane: is the sign consistent, and does the offset track any condition the bundles state
(temperature, time, pH, water activity)?

WHAT IT CANNOT DO
-----------------
A correlation localises; it does not identify. An offset that tracks temperature is equally
consistent with a wrong barrier, a wrong Q10, a missing temperature-dependent channel, or a
measurement whose efficiency changes with temperature -- wave B45 met that last case, where a purge
strips a hot cell better than a cool one. Sample sizes are small (8 fat rows, 6 trunk), so every
correlation is printed with its n and the reader is expected to use it.

NOTHING HERE IS A FIT AND NOTHING IS SCORED. No constant moves.

    python scripts/generators/generate_lane_offset_diagnostic.py
"""
from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SCORES = ROOT / "results" / "validation" / "core_panel_scores.json"
OUT_JSON = ROOT / "results" / "validation" / "lane_offset_diagnostic.json"
OUT_MD = ROOT / "results" / "validation" / "lane_offset_diagnostic.md"

#: The four condition fields the bundles state often enough to correlate against.
COVARIATES = ("temp_C", "time_min", "ph", "water_activity")

#: Pre-registered thresholds (sec. 3 of the prereg). Named here so the artifact can only be read
#: against the thresholds that were written down first.
SIGN_CONSISTENCY_SYSTEMATIC = 0.75
RHO_TRACKS = 0.6

#: ADDED AFTER THE FIRST RUN, AND LABELLED AS SUCH (prereg amendment 1, sec. 6).
#: The pre-registered rule called a lane/covariate pair a "track" on |rho| alone. That is not enough
#: and the first run showed why: the fat lane returned rho = -0.87 against temperature from 8 rows
#: that are 3 pots at 2 temperatures, and the trunk returned 83 % sign consistency from 6 rows that
#: are ONE pot. A rank correlation across two levels is a two-group comparison, not a trend, and no
#: |rho| threshold can tell the difference. So the artifact now also reports how many DISTINCT pots
#: and DISTINCT covariate levels each number rests on, and carries a second, stricter verdict.
#: The pre-registered verdict is kept beside it, unchanged, so the amendment cannot hide a result.
MIN_LEVELS_FOR_A_TREND = 3


def _rank(xs: Sequence[float]) -> List[float]:
    """Average ranks, so ties do not bias the correlation."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def spearman(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    """Spearman rho, or None when it is not defined (n < 3, or either side constant)."""
    if len(xs) < 3 or len(set(xs)) < 2 or len(set(ys)) < 2:
        return None
    rx, ry = _rank(xs), _rank(ys)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return None if den == 0 else num / den


def collect() -> List[Dict[str, Any]]:
    """One record per scored row that carries both a prediction and a measurement."""
    payload = json.loads(SCORES.read_text())
    rows: List[Dict[str, Any]] = []
    for bench in payload["benchmarks"]:
        cond = bench.get("conditions") or {}
        for comp in bench.get("compounds", []):
            predicted, measured = comp.get("predicted"), comp.get("measured_ppb")
            if not predicted or not measured or predicted <= 0 or measured <= 0:
                continue
            rows.append({
                "benchmark_id": bench["benchmark_id"],
                "compound": comp["compound"],
                "lane": comp.get("lane"),
                "in_core_fit": bool(comp.get("in_core_fit")),
                "dex": math.log10(float(predicted) / float(measured)),
                "declared_share_of_prediction": comp.get("declared_share_of_prediction"),
                **{c: cond.get(c) for c in COVARIATES},
            })
    return rows


def summarise_lane(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    dex = [r["dex"] for r in rows]
    median = statistics.median(dex)
    high = sum(1 for d in dex if d > 0)
    majority = max(high, len(dex) - high)
    correlations = {}
    for cov in COVARIATES:
        pairs = [(r[cov], r["dex"]) for r in rows if r.get(cov) is not None]
        rho = spearman([p[0] for p in pairs], [p[1] for p in pairs]) if pairs else None
        levels = len({p[0] for p in pairs})
        tracks = bool(rho is not None and abs(rho) >= RHO_TRACKS)
        correlations[cov] = {
            "rho": rho, "n": len(pairs),
            "distinct_levels": levels,
            "tracks": tracks,                                   # the PRE-REGISTERED verdict
            "tracks_with_design_support": bool(tracks and levels >= MIN_LEVELS_FOR_A_TREND),
            "why_not": None if levels >= MIN_LEVELS_FOR_A_TREND else
                       f"only {levels} distinct value(s) of {cov} in this lane: a rank correlation "
                       f"here is a {levels}-group comparison, not a trend",
        }
    return {
        "n": len(dex),
        "distinct_pots": len({r["benchmark_id"] for r in rows}),
        "median_signed_dex": median,
        "median_absolute_dex": statistics.median(abs(d) for d in dex),
        "reads": "high" if median > 0 else "low",
        "rows_reading_high": high,
        "sign_consistency": majority / len(dex),
        "systematic": bool(majority / len(dex) >= SIGN_CONSISTENCY_SYSTEMATIC),
        "dex_spread": (max(dex) - min(dex)),
        "correlations": correlations,
        "worst_row": max(rows, key=lambda r: abs(r["dex"]))["benchmark_id"],
        "worst_row_dex": max(abs(r["dex"]) for r in rows),
    }


def main() -> int:
    rows = collect()
    by_lane: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by_lane[r["lane"] or "unassigned"].append(r)

    lanes = {name: summarise_lane(rs) for name, rs in sorted(by_lane.items())}
    systematic = [n for n, s in lanes.items() if s["systematic"]]
    tracking: List[Tuple[str, str, float, int]] = [
        (name, cov, c["rho"], c["n"])
        for name, s in lanes.items() for cov, c in s["correlations"].items() if c["tracks"]
    ]
    supported: List[Tuple[str, str, float, int, int]] = [
        (name, cov, c["rho"], c["n"], c["distinct_levels"])
        for name, s in lanes.items() for cov, c in s["correlations"].items()
        if c["tracks_with_design_support"]
    ]

    payload = {
        "artifact": "lane_offset_diagnostic",
        "prereg": "results/validation/kinetic_core_b46_prereg.md",
        "what_this_is": (
            "The SIGNED offset log10(predicted/measured) per lane, its sign consistency, and its rank "
            "correlation against each condition the bundles state. A diagnostic, not a fit: nothing "
            "here is scored and no constant moves. A correlation LOCALISES a structural error; it does "
            "not identify one, and it cannot separate a wrong rate from a measurement whose efficiency "
            "changes with the same covariate."
        ),
        "thresholds": {"sign_consistency_systematic": SIGN_CONSISTENCY_SYSTEMATIC,
                       "rho_tracks": RHO_TRACKS},
        "summary": {
            "rows": len(rows),
            "lanes_systematic": sorted(systematic),
            "lane_covariate_pairs_tracking": [
                {"lane": l, "covariate": c, "rho": r, "n": n} for l, c, r, n in
                sorted(tracking, key=lambda t: -abs(t[2]))
            ],
            "lane_covariate_pairs_tracking_with_design_support": [
                {"lane": l, "covariate": c, "rho": r, "n": n, "distinct_levels": lv}
                for l, c, r, n, lv in sorted(supported, key=lambda t: -abs(t[2]))
            ],
            "lanes_whose_rows_are_one_pot": sorted(
                n for n, s in lanes.items() if s["distinct_pots"] == 1),
        },
        "lanes": lanes,
        "rows": sorted(rows, key=lambda r: (r["lane"] or "", -abs(r["dex"]))),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    # --- the readable half -------------------------------------------------
    lines = ["# Per-lane offset diagnostic (wave B46)", "",
             "Signed offset `dex = log10(predicted / measured)`. Positive means the model reads high.",
             "A diagnostic, not a fit: nothing here is scored and no constant moves.", "",
             "| lane | n | pots | median signed dex | reads | sign consistency | systematic | median abs dex |",
             "|---|---:|---:|---:|---|---:|---|---:|"]
    for name, s in sorted(lanes.items(), key=lambda kv: -kv[1]["median_absolute_dex"]):
        lines.append(f"| {name} | {s['n']} | {s['distinct_pots']} | {s['median_signed_dex']:+.2f} | {s['reads']} | "
                     f"{s['sign_consistency']:.0%} | {'YES' if s['systematic'] else 'no'} | "
                     f"{s['median_absolute_dex']:.2f} |")
    lines += ["", "## What each lane's offset tracks", "",
              "Spearman rank correlation of the signed offset against each stated condition. "
              f"|rho| >= {RHO_TRACKS} is the pre-registered threshold for 'tracks'. "
              "**n is small everywhere; read the correlations with it.**", "",
              "| lane | covariate | rho | n | distinct levels | tracks (pre-registered) | trend supported by the design |",
              "|---|---|---:|---:|---:|---|---|"]
    for name, s in sorted(lanes.items()):
        for cov, c in s["correlations"].items():
            rho = "n/a" if c["rho"] is None else f"{c['rho']:+.2f}"
            lines.append(f"| {name} | {cov} | {rho} | {c['n']} | {c['distinct_levels']} | "
                         f"{'YES' if c['tracks'] else ''} | "
                         f"{'YES' if c['tracks_with_design_support'] else ('' if not c['tracks'] else 'NO -- ' + (c['why_not'] or ''))} |")
    lines += ["", "## How to read a hit", "",
              "A lane that is systematic AND tracks a covariate is where a missing process is most "
              "likely to live, because parameter uncertainty has already been ruled out globally: "
              "`core_prediction_uncertainty.json` reports that uncapping every prior moves coverage "
              "from 19 % to 21 % against a nominal 90 %.", "",
              "It does NOT say which process. An offset that grows with temperature fits a wrong "
              "barrier, a wrong Q10, a missing temperature-dependent channel, or a measurement whose "
              "efficiency changes with temperature. Wave B45 met the last of those.", ""]
    OUT_MD.write_text("\n".join(lines) + "\n")

    print(f"wrote {OUT_JSON.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}: "
          f"{len(rows)} rows | systematic lanes: {sorted(systematic) or 'none'} | "
          f"tracking (pre-registered): {[(l, c, round(r, 2)) for l, c, r, n in tracking] or 'none'} | "
          f"WITH DESIGN SUPPORT: {[(l, c, round(r, 2), f'{lv} levels') for l, c, r, n, lv in supported] or 'none'}")
    for name, s in sorted(lanes.items()):
        print(f"  {name:11s} n={s['n']:2d} ({s['distinct_pots']} pots)  median {s['median_signed_dex']:+.2f} dex  "
              f"reads {s['reads']:4s}  sign consistency {s['sign_consistency']:.0%}"
              f"{'  SYSTEMATIC' if s['systematic'] else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
