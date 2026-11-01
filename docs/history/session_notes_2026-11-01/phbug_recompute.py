"""Scratch audit (read-only): re-predict the entry pH rows (Cao 2017 hold-out, Zhou 2023 score-only) of each tracked
joint-fit record at its recorded optimum, serially, with the trunk solutes OFF (what the forked parent used) and ON
(the task's own state). Writes only to the path given as argv[2]."""
import json
import sys
from pathlib import Path

sys.path.insert(0, "scripts/generators")
import numpy as np  # noqa: E402

import fit_joint_deoxyosone_thiol as joint  # noqa: E402
from fit_ph_acid_submodel import raw_table  # noqa: E402
from src.kinetic_core import fit_rows, ph_state, sulfur  # noqa: E402

ARM_OF_OLD = {"s042_e": "E", "s042_r3e": "R3E", "s042_r3en": "R3En", "s042_r3es": "R3Es"}
names = sys.argv[1].split(",")
out_path = Path(sys.argv[2])
VAL = Path("results/validation")

pots, rows = fit_rows.load(joint.ENTRY_PH_TABLE)
jobs = {}
for j in joint.probe.plan(pots, rows):
    if "skip" in j or any(pots[j["pot"]]["initial_mmol_L"].get(h) for h in joint.HEXOSES):
        continue
    if j["pot"] in joint.ENTRY_PH_HOLDOUT_POTS or j["pot"].startswith(joint.ENTRY_PH_SCORE_PREFIX):
        jobs[j["pot"]] = j
ph_rows = [r for r in joint.dedupe(r for r in rows if r.species == "pH" and r.pot in jobs)]
hold = [r for r in ph_rows if r.pot in joint.ENTRY_PH_HOLDOUT_POTS]
score = [r for r in ph_rows if r.pot not in joint.ENTRY_PH_HOLDOUT_POTS]

results = {}
for name in names:
    rec = json.loads((VAL / f"kinetic_core_joint_deoxyosone_thiol_fit_{name}.json").read_text())
    arm = rec.get("arm") or ARM_OF_OLD[name]
    arms = joint.ARMS[arm]
    joint.STRECKER_SHARE = float(rec.get("strecker_share", 0.042))
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    joint.TDP_EXIT_STATIC = dict(rec.get("tdp_exit") or {})
    joint.EA_STATIC = dict(rec.get("ea_static") or {})
    sulfur.PH_COUPLING = "inline"
    theta = rec["fit"]["theta"]
    keys = tuple(theta)
    tr = joint.trial([theta[k] for k in keys], keys)
    res = {p: joint.run_task(dict(j, **tr)) for p, j in jobs.items()}
    entry = {}
    for flag in (False, True):
        joint.install(dict(next(iter(jobs.values())), **tr))
        ph_state.TRUNK_SOLUTES_IN_BALANCE = flag
        pred = {id(r): joint.probe.predict(r, dict(jobs[r.pot], **tr), res[r.pot]) for r in ph_rows}
        entry["solutes_on" if flag else "solutes_off"] = {
            "holdout": raw_table(hold, [fit_rows.residual(r, pred[id(r)]) for r in hold]),
            "score_only": raw_table(score, [fit_rows.residual(r, pred[id(r)]) for r in score]),
            "max_abs_row_shift": None}
        entry.setdefault("_pred", {})[flag] = [pred[id(r)] for r in ph_rows]
    shifts = np.abs(np.array(entry["_pred"][True]) - np.array(entry["_pred"][False]))
    entry["max_abs_ph_shift"] = float(shifts.max())
    entry["rows_shift_gt_1e-3"] = int((shifts > 1e-3).sum())
    del entry["_pred"]
    entry["recorded"] = {"holdout": {k: v for k, v in rec["holdout_raw_at_optimum"].items() if "pH" in k},
                         "score_only": {k: v for k, v in rec["score_only_raw_at_optimum"].items() if "pH" in k}}
    results[name] = entry
    print(name, json.dumps({k: entry[k] for k in ("max_abs_ph_shift", "rows_shift_gt_1e-3")}), flush=True)
    out_path.write_text(json.dumps(results, indent=1))
