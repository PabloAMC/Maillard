"""Re-run the carried rival points on a subset of the S3 sensitivity recipes and compare with the tracked detail.

Writes nothing tracked: rows go to OUT (a JSON-lines checkpoint outside the repo's data/ and results/).
Env: SUBSET_EVERY (every Nth grid recipe; default 6), OUT, WORKERS.
"""
import gzip
import json
import math
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "scripts", "generators"))
import build_ph_balance_s3_sensitivity as g  # noqa: E402

every = int(os.environ.get("SUBSET_EVERY", "6"))
out = Path(os.environ["OUT"])
workers = int(os.environ.get("WORKERS", "4"))
grid = g.grid_recipes()
pan = g.panel_recipes()
dirr, claims = g.directional_recipes()
recipes = pan + dirr + grid[::every]
rivals = {k: v for k, v in g.point_configs().items() if k not in ("B0", "B0_inline")}
print(f"subset: panel {len(pan)}, directional {len(dirr)}, grid {len(grid[::every])} of {len(grid)}; "
      f"{len(rivals)} rival points", flush=True)

_eval = g.evaluate_recipe


def rivals_only(recipe):
    return _eval(recipe, points=rivals)


g.evaluate_recipe = rivals_only
t0 = time.time()
rows = g.run(recipes, workers, log=lambda m: print(m, flush=True), checkpoint=out)
print(f"wall {time.time() - t0:.0f} s", flush=True)

ref = {r["rid"]: r for r in json.load(gzip.open("results/validation/ph_balance_s3_sensitivity_detail.json.gz"))["recipes"]}
n_vals = n_moved6 = n_moved3 = 0
worst = []
newly_ok, newly_failed, still_failed = [], [], []
ph_moves = []
th = g.thresholds()
verdict_changes = []
for r in rows:
    old = ref[r["rid"]]
    for p, v in r["points"].items():
        o = old["points"].get(p, {})
        if v.get("not_applicable"):
            continue
        if v.get("failed") and not o.get("failed"):
            newly_failed.append((r["rid"], p, v["failed"][:80]))
        elif o.get("failed") and not v.get("failed"):
            newly_ok.append((r["rid"], p))
        elif o.get("failed") and v.get("failed"):
            still_failed.append((r["rid"], p))
        if v.get("failed") or o.get("failed"):
            continue
        for c, x in v["concentrations_ug_per_l"].items():
            y = o["concentrations_ug_per_l"].get(c)
            if y is None:
                continue
            n_vals += 1
            if x > g.ZERO_FLOOR_UG_PER_L and y > g.ZERO_FLOOR_UG_PER_L:
                d = abs(math.log10(x / y))
            elif (x > g.ZERO_FLOOR_UG_PER_L) != (y > g.ZERO_FLOOR_UG_PER_L):
                d = math.inf
            else:
                d = 0.0
            n_moved6 += d > 1e-6
            n_moved3 += d > 1e-3
            worst.append((d, r["rid"], p, c, y, x))
        for k in ("ph_cooled_end", "ph_in_situ_end"):
            if v.get(k) is not None and o.get(k) is not None:
                ph_moves.append(abs(v[k] - o[k]))
    # level verdicts: spread over B0 + all rivals, old vs new
    if all(not v.get("failed") for v in r["points"].values()) and not old.get("clamped") and r["points"]:
        for c in (old["points"].get("B0", {}).get("concentrations_ug_per_l") or {}):
            def spread(src):
                vals = [old["points"]["B0"]["concentrations_ug_per_l"][c]] + [
                    src["points"][p]["concentrations_ug_per_l"].get(c) for p in rivals if p in src["points"]]
                if any(v is None for v in vals):
                    return None
                return g.log_spread(vals)
            if any(old["points"].get(p, {}).get("failed") for p in rivals):
                continue
            so, sn = spread(old), spread(r)
            if so is None or sn is None:
                continue
            if (so > th["tau_dex"]) != (sn > th["tau_dex"]):
                verdict_changes.append((r["rid"], c, so, sn))
worst.sort(key=lambda t: -t[0] if math.isfinite(t[0]) else -1e9)
print(f"values compared {n_vals}; moved > 1e-6 dex: {n_moved6}; > 1e-3 dex: {n_moved3}")
print("largest moves (dex, rid, point, compound, old, new):")
for w in worst[:12]:
    print("  ", w)
print("pH max |move|:", max(ph_moves) if ph_moves else None, "n", len(ph_moves))
print("newly finishing:", len(newly_ok), sorted({x[0] for x in newly_ok}))
print("newly failing:", newly_failed[:10])
print("still failing:", len(still_failed), sorted({x[0] for x in still_failed}))
print("level verdict flips (tau):", verdict_changes[:20])
