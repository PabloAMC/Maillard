#!/usr/bin/env python
"""
check_feasibility_must_not_break.py -- must-not-break checks at a feasibility test's joint passes (2026-10-31).

The partition signatures cannot tell candidate structures apart (any fast 3-DX sink faster at high pH passes them;
dpo partition prereg sec. 21, amendment), so the discrimination rests on rows the structure was not designed on: the
calibrated Zhai MGO rows, free cysteine and TTCA in Zhai 2021 and Cao 2017's blanks (the route table), the thiol
signature groups and Cerny 2003's pot D. Prints them for the reference optimum and each joint pass (no fit).

    python scripts/generators/check_feasibility_must_not_break.py --arm B7ASST2C --feasibility b7asst2c
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402
import score_mechanism_signatures as sc  # noqa: E402


def _rows_of(table, keep):
    pots, rows = fit_rows.load(table)
    sel = [r for r in rows if keep(r)]
    plans = {j["pot"]: j for j in probe.plan(pots, rows) if j["pot"] in {r.pot for r in sel} and "skip" not in j}
    return sel, plans


def _resid(trial, sel, plans, group):
    out = {}
    for pot, job in plans.items():
        j = dict(job, **trial, strecker_share=0.042)
        res = joint.run_task(j)
        if "failed" in res:
            continue
        for r in sel:
            if r.pot == pot:
                p = probe.predict(r, j, res)
                if p > 0 and r.value > 0:
                    out.setdefault(group(r), []).append(math.log10(p / r.value))
    return {k: round(float(np.mean(v)), 2) for k, v in sorted(out.items())}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--arm", required=True)
    parser.add_argument("--feasibility", required=True, help="the feasibility record's stem suffix (b7asst2c, ...)")
    parser.add_argument("--reference", default="s042_b7ass_w02")
    parser.add_argument("--max", type=int, default=8)
    args = parser.parse_args(argv)
    feas = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_deoxyosone_feasibility_{args.feasibility}.json").read_text())
    ref = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.reference}.json").read_text())
    sc._configure(args.arm, 0.042)
    keys, _, _, _, start = joint.coordinates(False, joint.ARMS[args.arm])
    tables, measured, jobs = sc._measured()
    mgo = _rows_of("deoxyosone_submodel", lambda r: r.species == "MGO" and r.ledger_id.startswith("zhai")
                   and r.method not in ("semi_quantitative", "relative"))
    cys = _rows_of("route_submodel", lambda r: r.species in ("Cys", "TTCA") and r.pot.startswith(("zhai2021", "cao2017")))
    ref_theta = {k: ref["fit"]["theta"].get(k, s) for k, s in zip(keys, start)}
    cands = [("reference", ref_theta, None)] + [(f"joint {i}", jp["theta"], jp["tdp_exit"])
                                                 for i, jp in enumerate(feas["joint_passes"][:args.max])]
    for name, th, ex in cands:
        joint.TDP_EXIT_STATIC = (joint.tdp_exit_static(joint.ARMS[args.arm], **ex) if ex
                                 else joint.tdp_exit_static(joint.ARMS[args.arm], k_ref_per_min=1e-6)
                                 if "tdp_exit" in joint.ARMS[args.arm] else {})
        trial = joint.trial([th[k] for k in keys], keys)
        m = _resid(trial, *mgo, lambda r: "MGO")
        c = _resid(trial, *cys, lambda r: ("cao" if r.pot.startswith("cao") else "zhai21") + "_" + r.species)
        _, groups, cerny = sc.score(th, 0.042, tables, measured, jobs, strict=False)
        print(f"{name:9} trap {th.get('log10_k_tdp_cys_trap', float('nan')):+.2f} | MGO {m.get('MGO')} | {c} | "
              f"thiol groups {sum(groups.values())}/11 | potD {cerny['cerny2003_potD_nf_ribose_cys']['fed_nf_mft_share']:.3f}",
              flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
