#!/usr/bin/env python
"""Re-score the entry pH rows of the tracked joint-fit records in the state each pot was integrated in.

Why (docs/history/session_notes_2026-11-01/forked_ph_bug_audit.md): before A3 phase 2 (PR #36), a forked pool
(MAILLARD_WORKERS > 1) predicted every row in the PARENT's engine state. For an entry arm, the integration had the
trunk solutes (glycine, the Amadori and Schiff-base species) in its charge balance, but the parent's cooled-sample
pH balance did not. Every fitted row, theta, cost and Laplace is unaffected (no fitted row is a pH row). The
hold-out Cao 2017 pH rows (`holdout_raw_at_optimum`, `holdout_raw_at_start`, `holdout_after_scale_at_optimum`:
the class is absolute, so after-scale equals raw) and the score-only Zhou 2023 pH rows were computed wrongly.

What this does, per record whose arm carries the pentose + amine entry (the only arms with trunk solutes):
- rebuilds the entry pH jobs and the arm's statics from the record (strecker share, DDP law, TDP exit, Ea statics);
- integrates each pot once per point (the recorded optimum and the recorded start), serially, in its own Config;
- predicts every pH row twice: in the pot's Config (CORRECTED) and with the trunk solutes out of the cooled
  balance (AS RECORDED, the forked parent's state);
- checks that AS RECORDED reproduces the record's own numbers (|diff| <= 1e-6 on every mean and RMS). A record that
  does not reproduce is reported as such, and its corrected numbers are not used.

Nothing is refitted; the records are not edited. Writes
``results/validation/kinetic_core_joint_fit_entry_ph_rescore.{json,md}``.

    ./scripts/docker_maillard.sh run "python scripts/generators/rescore_joint_fit_entry_ph_holdout.py"
    ./scripts/docker_maillard.sh run "python scripts/generators/rescore_joint_fit_entry_ph_holdout.py --only s042_b7st_w02"
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
_GEN = str(ROOT / "scripts" / "generators")
if _GEN not in sys.path:
    sys.path.insert(0, _GEN)

import fit_joint_deoxyosone_thiol as joint  # noqa: E402
from fit_ph_acid_submodel import raw_table  # noqa: E402

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import config, fit_rows, parallel  # noqa: E402

RECORD_PREFIX = "kinetic_core_joint_deoxyosone_thiol_fit_"
OUT = data_paths.VALIDATION_DIR / "kinetic_core_joint_fit_entry_ph_rescore"
# records written before the "arm" field existed (the audit's ARM_OF_OLD)
ARM_OF_OLD = {"s042_e": "E", "s042_r3e": "R3E", "s042_r3en": "R3En", "s042_r3es": "R3Es"}
CAO_KEY = "cao2017_fig4c_ph | pH"
REPRODUCE_TOL = 1e-6


def record_names() -> List[str]:
    """Every tracked joint-fit record whose arm has the entry (trunk solutes in the integration's balance)."""
    out = []
    for path in sorted(data_paths.VALIDATION_DIR.glob(f"{RECORD_PREFIX}*.json")):
        name = path.stem[len(RECORD_PREFIX):]
        rec = json.loads(path.read_text(encoding="utf-8"))
        arm = rec.get("arm") or ARM_OF_OLD.get(name, "none")
        if "entry" in joint.ARMS.get(arm, frozenset()):
            out.append(name)
    return out


def set_statics(rec: Dict[str, Any], arms) -> None:
    """The module statics ``trial`` / ``job_config`` read, as the record's run had them (fit_joint main())."""
    joint.STRECKER_SHARE = float(rec.get("strecker_share", 0.042))
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    joint.TDP_EXIT_STATIC = dict(rec.get("tdp_exit") or {})
    joint.EA_STATIC = dict(rec.get("ea_static") or {})


def entry_ph_jobs():
    pots, rows = fit_rows.load(joint.ENTRY_PH_TABLE)
    jobs = {}
    for j in joint.probe.plan(pots, rows):
        if "skip" in j or any(pots[j["pot"]]["initial_mmol_L"].get(h) for h in joint.HEXOSES):
            continue
        if j["pot"] in joint.ENTRY_PH_HOLDOUT_POTS or j["pot"].startswith(joint.ENTRY_PH_SCORE_PREFIX):
            jobs[j["pot"]] = j
    ph_rows = list(joint.dedupe(r for r in rows if r.species == "pH" and r.pot in jobs))
    return jobs, ph_rows


def score_point(theta: Dict[str, float], jobs, ph_rows) -> Dict[str, Any]:
    keys = tuple(theta)
    tr = joint.trial([theta[k] for k in keys], keys)
    tasks = {p: parallel.Configured(dict(j, **tr), joint.job_config(dict(j, **tr))) for p, j in jobs.items()}
    runs = {p: parallel.run_under(t, lambda t=t: joint.run_pot(t.item)) for p, t in tasks.items()}
    out: Dict[str, Any] = {}
    preds: Dict[str, List[float]] = {}
    for mode, solutes in (("corrected", None), ("as_recorded", False)):
        vals = []
        for r in ph_rows:
            task = tasks[r.pot]
            if solutes is not None:
                task = parallel.Configured(task.item, task.config.replace({"ph.trunk_solutes_in_balance": solutes}))
            got = runs[r.pot]
            if "failed" in got:
                vals.append(float("nan"))
                continue
            vals.append(parallel.run_under(task, lambda r=r, t=task, g=got: joint.probe.predict(r, t.item, g)))
        preds[mode] = vals
        hold = [(r, v) for r, v in zip(ph_rows, vals) if r.pot in joint.ENTRY_PH_HOLDOUT_POTS]
        score = [(r, v) for r, v in zip(ph_rows, vals) if r.pot not in joint.ENTRY_PH_HOLDOUT_POTS]
        out[mode] = {
            "holdout": raw_table([r for r, _ in hold], [fit_rows.residual(r, v) for r, v in hold]),
            "score_only": raw_table([r for r, _ in score], [fit_rows.residual(r, v) for r, v in score]),
        }
    shifts = [abs(a - b) for a, b in zip(preds["corrected"], preds["as_recorded"])]
    out["max_abs_ph_shift"] = max(shifts) if shifts else 0.0
    out["rows_shift_gt_1e-3"] = sum(s > 1e-3 for s in shifts)
    return out


def reproduces(recomputed: Dict[str, Any], recorded: Dict[str, Any]) -> bool:
    if set(recomputed) != set(recorded):
        return False
    return all(abs(recomputed[k][f] - recorded[k][f]) <= REPRODUCE_TOL for k in recorded for f in ("mean", "rms"))


def rescore(name: str, jobs, ph_rows) -> Dict[str, Any]:
    rec = json.loads((data_paths.VALIDATION_DIR / f"{RECORD_PREFIX}{name}.json").read_text(encoding="utf-8"))
    arm = rec.get("arm") or ARM_OF_OLD[name]
    set_statics(rec, joint.ARMS[arm])
    ph_keys = lambda table: {k: v for k, v in (table or {}).items() if k.endswith("| pH")}  # noqa: E731
    recorded = {"optimum": {"holdout": ph_keys(rec.get("holdout_raw_at_optimum")),
                            "score_only": ph_keys(rec.get("score_only_raw_at_optimum"))},
                "start": {"holdout": ph_keys(rec.get("holdout_raw_at_start"))}}
    entry: Dict[str, Any] = {"arm": arm, "recorded": recorded}
    for point, theta in (("optimum", rec["fit"]["theta"]), ("start", rec["start"])):
        got = score_point(dict(theta), jobs, ph_rows)
        checks = {part: reproduces(got["as_recorded"][part], recorded[point][part]) for part in recorded[point]}
        got["reproduces_record"] = all(checks.values())
        got["reproduce_checks"] = checks
        entry[point] = got
    return entry


def verdict_line(name: str, e: Dict[str, Any]) -> str:
    def cell(point, mode):
        v = e[point][mode]["holdout"].get(CAO_KEY)
        return "-" if v is None else f"{v['mean']:+.3f} / {v['rms']:.3f} ({'within' if v['within'] else 'FAIL'})"

    rep = "yes" if e["optimum"]["reproduces_record"] and e["start"]["reproduces_record"] else "**NO**"
    return (f"| {name} | {e['arm']} | {rep} | {cell('optimum', 'as_recorded')} | {cell('optimum', 'corrected')} | "
            f"{cell('start', 'as_recorded')} | {cell('start', 'corrected')} | {e['optimum']['max_abs_ph_shift']:.3f} |")


def render(payload: Dict[str, Any]) -> str:
    lines = ["# Joint-fit entry pH rows re-scored in each pot's own state (forked-pool pH bug)", "",
             "Generated by `scripts/generators/rescore_joint_fit_entry_ph_holdout.py`. Audit: "
             "`docs/history/session_notes_2026-11-01/forked_ph_bug_audit.md`. Fix: PR #36 (A3 phase 2).", "",
             "Before PR #36 a forked pool predicted each row in the parent's state, so an entry arm's pH rows lost "
             "the trunk solutes from the cooled-sample charge balance. Fitted rows, theta, cost and Laplace are not "
             "affected. Only the pH rows below are. Nothing is refitted and no record is edited.", "",
             f"Cao 2017 pH hold-out (`{CAO_KEY}`, n 7, tolerance 0.3, absolute class so after-scale = raw): "
             "mean / RMS (verdict). \"Reproduces\" = the AS RECORDED mode matches the record's own numbers to "
             f"{REPRODUCE_TOL:g} at both points.", "",
             "| record | arm | reproduces | optimum, as recorded | optimum, corrected | start, as recorded | "
             "start, corrected | max pH shift at optimum |",
             "|---|---|---|---|---|---|---|---|"]
    lines += [verdict_line(n, e) for n, e in payload["records"].items()]
    flips = [n for n, e in payload["records"].items()
             if (e["optimum"]["as_recorded"]["holdout"].get(CAO_KEY) or {}).get("within")
             != (e["optimum"]["corrected"]["holdout"].get(CAO_KEY) or {}).get("within")]
    zhou = {n: e["optimum"]["max_abs_ph_shift"] for n, e in payload["records"].items()}
    zhou_rms_moves = [abs(e[p]["corrected"]["score_only"][k]["rms"] - e[p]["as_recorded"]["score_only"][k]["rms"])
                      for e in payload["records"].values() for p in ("optimum", "start")
                      for k in e[p]["corrected"]["score_only"]]
    lines += ["", "## Summary", "",
              f"- Records re-scored: {len(payload['records'])}; reproducing their record: "
              f"{sum(e['optimum']['reproduces_record'] and e['start']['reproduces_record'] for e in payload['records'].values())}.",
              f"- Cao verdict at the optimum changes in: {', '.join(flips) or 'none'}.",
              f"- Largest pH shift of any row at an optimum: {max(zhou.values()):.3f}." if zhou else "",
              "- Zhou 2023 score-only pH (no trunk solutes in those pots): largest RMS change over every record and "
              f"point {max(zhou_rms_moves, default=0.0):.1e} (the audit expected none).",
              f"- Wall time {payload['wall_seconds']:.0f} s, serial."]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", help="comma-separated record suffixes (default: every entry-arm record)")
    args = parser.parse_args(argv)
    t0 = time.time()
    config.update({"ph.coupling": "inline"})   # fit_joint main()'s base for every task
    names = args.only.split(",") if args.only else record_names()
    jobs, ph_rows = entry_ph_jobs()
    records = {}
    for name in names:
        records[name] = rescore(name, jobs, ph_rows)
        print(name, "reproduces:", records[name]["optimum"]["reproduces_record"],
              "max shift:", round(records[name]["optimum"]["max_abs_ph_shift"], 4), flush=True)
    payload = {
        "artifact": OUT.name,
        "provenance": provenance.provenance_block(
            OUT.name, generated_by="scripts/generators/rescore_joint_fit_entry_ph_holdout.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{joint.ENTRY_PH_TABLE}.yml",
                    ROOT / "scripts/generators/fit_joint_deoxyosone_thiol.py",
                    ROOT / "scripts/generators/build_submodel_probe.py"]
                   + [data_paths.VALIDATION_DIR / f"{RECORD_PREFIX}{n}.json" for n in names]),
        "audit": "docs/history/session_notes_2026-11-01/forked_ph_bug_audit.md",
        "reproduce_tolerance": REPRODUCE_TOL,
        "records": records,
        "wall_seconds": time.time() - t0,
    }
    if args.only:   # a partial run is a check, never the tracked artifact
        print(json.dumps({n: {p: records[n][p]["corrected"]["holdout"].get(CAO_KEY) for p in ("optimum", "start")}
                          for n in names}, indent=1))
        return 0
    OUT.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=float) + "\n", encoding="utf-8")
    OUT.with_suffix(".md").write_text(render(payload) + "\n", encoding="utf-8")
    print(f"wrote {OUT}.{{json,md}} in {payload['wall_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
