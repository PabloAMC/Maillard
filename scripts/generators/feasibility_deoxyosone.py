#!/usr/bin/env python
"""
feasibility_deoxyosone.py -- can a STRUCTURE reproduce the deoxyosone partition's signatures anywhere in its bands?

Owner, 2026-10-31: "let's focus on getting the mechanisms right, before getting the fitting constants right"
(STATUS_AND_PLAN sec. 5, mechanisms before constants). A structure is tested for feasibility before any fit: the
arm's structural constants are drawn uniformly within declared bands (the rest held at a reference optimum) and every
draw is scored on the scale-free partition signatures of score_deoxyosone_signatures.py. A structure no draw can pass
is rejected without fitting. Design: results/validation/kinetic_core_dpo_partition_prereg.md sec. 16.

    python scripts/generators/feasibility_deoxyosone.py --arm B7ASST --reference s042_b7ass_w02 --draws 300
Writes results/validation/kinetic_core_deoxyosone_feasibility_<arm>.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402
import score_deoxyosone_signatures as sd  # noqa: E402
import score_mechanism_signatures as sc  # noqa: E402

#: the declared bands (prereg sec. 16): (theta key, low, high) -- absolute, or "+-" relative to the reference
BANDS: Tuple[Tuple[str, str, float, float], ...] = (
    ("dpo_share7", "abs", 0.02, 0.98),
    ("split_ph_c", "abs", 5.0, 10.0),
    ("aa_pk", "abs", 5.0, 7.5),
    ("log10_k_dpo_aa", "rel", -1.0, 1.0),
    ("offset_k_dpo_c2c3", "abs", 0.0, 3.0),
    ("offset_k_dpo_ddp", "rel", -1.0, 1.0),
    # v2 (prereg sec. 17): the sources -- applied only when the arm carries them
    ("offset_k_pent_caramel", "abs", -3.0, 0.5),
    ("offset_k_ttca_deg", "abs", -3.0, 0.5),
    ("offset_k_dpo_nf", "rel", -1.0, 1.0),
    # v4 (prereg sec. 19): the base-favoured 3-DX exit to MGO (log10 k at 100 C, pH 5; the base factor rises with pH)
    ("log10_k_tdp_mgo", "abs", -5.0, -1.0),
    # v5 (prereg sec. 20): the cysteine trap of 3-DX, L/(mmol min) at 100 C / pH 5
    ("log10_k_tdp_cys_trap", "abs", -5.0, 0.0),
)
TDP_EXIT_BAND = {"log10_k_ref_per_min": (-3.5, -1.0), "ea_kj_mol": (30.0, 100.0)}
THRESHOLDS = {"time_slope": 5, "ph_segments": 4, "xylose_double_ratio": 3}


def _synthetic_rows(args, rows, plans, share):
    """The power test (prereg sec. 23): the TDP / DPO rows replaced by a known structure's predictions with log-normal
    noise -- the truth is a joint pass of another feasibility record (its arm, constants and TDP exit)."""
    import dataclasses
    stem, index = args.synthetic_from.rsplit(":", 1)
    rec = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_deoxyosone_feasibility_{stem}.json").read_text())
    truth = rec["joint_passes"][int(index)]
    saved = (dict(joint.TDP_EXIT_STATIC), dict(joint.EA_STATIC), dict(joint.DDP_LAW_STATIC))
    sc._configure(rec["arm"], share)
    joint.TDP_EXIT_STATIC = joint.tdp_exit_static(joint.ARMS[rec["arm"]], **truth["tdp_exit"]) if truth["tdp_exit"] else {}
    tkeys, *_ = joint.coordinates(False, joint.ARMS[rec["arm"]])
    trial = joint.trial([truth["theta"][k] for k in tkeys], tkeys)
    rng = np.random.default_rng(args.seed + 7)
    out = []
    for r in rows:
        if r.pot in plans and r.species in ("TDP", "DPO"):
            j = dict(plans[r.pot], **trial, strecker_share=share)
            res = joint.run_task(j)
            p = probe.predict(r, j, res) if "failed" not in res else float("nan")
            out.append(dataclasses.replace(r, value=float(p) * 10.0 ** (args.synthetic_sigma * rng.standard_normal())))
        else:
            out.append(r)
    joint.TDP_EXIT_STATIC, joint.EA_STATIC, joint.DDP_LAW_STATIC = saved
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--arm", required=True, choices=tuple(joint.ARMS))
    parser.add_argument("--reference", required=True, help="the fit whose optimum holds the unsampled constants")
    parser.add_argument("--draws", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--tag", default="", help="suffix for the output stem (a rerun with another seed or size)")
    parser.add_argument("--synthetic-from", metavar="RECORD:INDEX", help="power test (prereg sec. 23): replace the "
                        "TDP / DPO rows by joint pass INDEX of feasibility record RECORD (e.g. b7asst2m:0) plus noise")
    parser.add_argument("--synthetic-sigma", type=float, default=0.15, help="log10 noise per synthetic row")
    args = parser.parse_args(argv)
    from src.kinetic_core import parallel

    ref_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.reference}.json"
    ref = json.loads(ref_path.read_text())
    sc._configure(args.arm, float(ref["strecker_share"]))
    keys, _kinds, _lo, _hi, start = joint.coordinates(False, joint.ARMS[args.arm])
    base = {k: float(ref["fit"]["theta"].get(k, s0)) for k, s0 in zip(keys, start)}
    rng = np.random.default_rng(args.seed)
    draws: List[Dict[str, Any]] = []
    for _ in range(args.draws):
        theta = dict(base)
        for key, mode, lo, hi in BANDS:
            if key in theta:
                u = lo + (hi - lo) * rng.random()
                theta[key] = u if mode == "abs" else base[key] + u
        exit_ = None
        if "tdp_exit" in joint.ARMS[args.arm]:
            lk = TDP_EXIT_BAND["log10_k_ref_per_min"]
            ea = TDP_EXIT_BAND["ea_kj_mol"]
            exit_ = {"k_ref_per_min": 10.0 ** (lk[0] + (lk[1] - lk[0]) * rng.random()),
                     "ea_kj_mol": ea[0] + (ea[1] - ea[0]) * rng.random()}
        draws.append({"theta": theta, "tdp_exit": exit_})
    pots, rows = fit_rows.load(sd.TABLE)
    plans = {j["pot"]: j for j in probe.plan(pots, rows) if j["pot"] in sd._wanted()}
    if args.synthetic_from:
        rows = _synthetic_rows(args, rows, plans, float(ref["strecker_share"]))

    def one(draw):
        try:
            joint.TDP_EXIT_STATIC = (joint.tdp_exit_static(joint.ARMS[args.arm], **draw["tdp_exit"])
                                     if draw["tdp_exit"] else {})
            th = draw["theta"]
            _, passed = sd.score(joint.trial([th[k] for k in keys], keys), float(ref["strecker_share"]), plans, rows,
                                 strict=False)
            return {"passed": passed} if passed is not None else {"error": "a pot failed"}
        except Exception as exc:  # noqa: BLE001 - a draw the engine cannot run is counted, not hidden
            return {"error": repr(exc)[:200]}

    results = parallel.map_ordered(one, draws)
    ok = [(d, r["passed"]) for d, r in zip(draws, results) if "error" not in r]
    reach = {fam: sum(1 for _, p in ok if p[fam] >= t) for fam, t in THRESHOLDS.items()}
    best = {fam: max((p[fam] for _, p in ok), default=None) for fam in THRESHOLDS}
    joint_ok = [(d, p) for d, p in ok if all(p[f] >= t for f, t in THRESHOLDS.items())]
    pairs = {f"{a}+{b}": sum(1 for _, p in ok if p[a] >= THRESHOLDS[a] and p[b] >= THRESHOLDS[b])
             for a, b in (("time_slope", "ph_segments"), ("time_slope", "xylose_double_ratio"),
                          ("ph_segments", "xylose_double_ratio"))}
    stem = data_paths.VALIDATION_DIR / (f"kinetic_core_deoxyosone_feasibility_{args.arm.lower()}"
                                        + (f"_{args.tag}" if args.tag else ""))
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/feasibility_deoxyosone.py",
                   inputs=[ref_path, ROOT / "src/kinetic_core/sulfur.py", fit_rows.FIT_ROWS_DIR / f"{sd.TABLE}.yml"]),
               "arm": args.arm, "reference": args.reference, "draws": args.draws, "seed": args.seed,
               "runnable_draws": len(ok), "bands": [list(b) for b in BANDS], "tdp_exit_band": TDP_EXIT_BAND,
               "thresholds": THRESHOLDS, "draws_reaching": reach, "best": best, "pairs_reaching": pairs,
               "feasible": bool(joint_ok), "joint_passes": [{"theta": d["theta"], "tdp_exit": d["tdp_exit"],
                                                            "passed": p} for d, p in joint_ok[:20]],
               # v3 (prereg sec. 18): where each family is reached, so a tension between families can be read
               "reaching_draws": {fam: [{"sampled": {k: d["theta"][k] for k, *_ in BANDS if k in d["theta"]},
                                         "tdp_exit": d["tdp_exit"], "passed": p}
                                        for d, p in ok if p[fam] >= t][:40] for fam, t in THRESHOLDS.items()},
               "errors_sample": [r["error"] for r in results if "error" in r][:5]}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload, default=float) + "\n")
    lines = [f"# Deoxyosone partition feasibility: arm {args.arm}", "",
             f"*Generated by `scripts/generators/feasibility_deoxyosone.py --arm {args.arm} --reference {args.reference} "
             f"--draws {args.draws}` (seed {args.seed}); design: kinetic_core_dpo_partition_prereg.md sec. 16. The arm's "
             f"structural constants drawn uniformly in declared bands, the rest at {args.reference}'s optimum; each draw "
             f"scored on the scale-free partition signatures. {len(ok)} of {args.draws} draws ran.*", "",
             f"**Feasible: {'YES' if joint_ok else 'NO'}** ({len(joint_ok)} draw(s) pass all three families).", "",
             "| family | threshold | draws reaching it | best draw |", "|---|---|---|---|"]
    lines += [f"| {fam} | >= {t} | {reach[fam]} | {best[fam]} |" for fam, t in THRESHOLDS.items()]
    lines += ["", "| pair | draws reaching both |", "|---|---|"] + [f"| {k} | {v} |" for k, v in pairs.items()]
    if joint_ok:
        lines += ["", "## Joint passes (first 20): the constants", ""]
        lines += [f"- {json.dumps({k: round(v, 3) for k, v in d['theta'].items() if k in dict((b[0], 1) for b in BANDS)})}"
                  f" tdp_exit {d['tdp_exit']}: {p}" for d, p in joint_ok[:20]]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: feasible={bool(joint_ok)} reach={reach} best={best}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
