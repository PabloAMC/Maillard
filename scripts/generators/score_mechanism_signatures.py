#!/usr/bin/env python
"""
score_mechanism_signatures.py -- the mechanism signatures of a joint-fit arm (STRATEGY sec. 9, 2026-10-30).

Owner: "we should make sure we are figuring out the right mechanisms, not the specific numbers". A signature is a
within-paper log ratio that a wrong mechanism cannot reproduce whatever its constants -- a route share, a lag, a pH or
temperature direction, a plateau at a precursor's exhaustion, a precursor-addition contrast. Each is computed from
the two-reader keyed tables (measured) and from the engine at a fit's optimum (predicted), and judged against a
tolerance declared here, before any arm it scores was read (2026-10-30: arm E running, round 3 not yet run).

    python scripts/generators/score_mechanism_signatures.py --fit s042_r2b
Writes results/validation/kinetic_core_mechanism_signatures_<fit>.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402

Point = Tuple[str, str, str, float]   # (table, pot, species, t_min)
T = "thiol_submodel"
Z = "thiol_zhou2023_ph_series"
D = "ph_acid_davidek_acids"

#: (name, mechanism tested, numerator points, denominator points, tolerance in dex, sign must match)
SIGNATURES: Tuple[Tuple[str, str, List[Point], List[Point], float, bool], ...] = (
    ("schieberle MFT lag 60 / 720 min", "an early MFT source (not only accumulating norfuraneol)",
     [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 60.0)], [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 720.0)],
     0.3, True),
    ("schieberle FFT lag 60 / 720 min", "FFT's early supply",
     [(T, "schieberle2000_rib_cys_ph5_100C", "FFT", 60.0)], [(T, "schieberle2000_rib_cys_ph5_100C", "FFT", 720.0)],
     0.3, True),
    # v2 (2026-10-31, before any converged arm was read): the SHAPE of Schieberle's MFT -- near-linear after a ~9 min
    # dead time, then a stop near 6 h (720 / 360 = 1.15) -- which an accumulating-norfuraneol route (>= 3x over the
    # same interval in the s042 engine) and a draining reservoir predict differently (literature pass, prereg sec. 14)
    ("schieberle MFT plateau 720 / 360 min", "a reservoir that runs dry near 6 h, not accumulating norfuraneol",
     [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 720.0)], [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 360.0)],
     0.25, False),
    ("schieberle MFT early rise 60 / 30 min", "near-linear formation after a short dead time (no slow intermediate)",
     [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 60.0)], [(T, "schieberle2000_rib_cys_ph5_100C", "MFT", 30.0)],
     0.25, True),
    ("hofmann MFT pH 7 / pH 3 (hold-out)", "the pH direction of MFT formation",
     [(T, "hofmann1998_cys_ribose_ph7", "MFT", 20.0)], [(T, "hofmann1998_cys_ribose_ph3", "MFT", 20.0)], 0.5, True),
    ("hofmann FFT pH 7 / pH 3 (hold-out)", "the pH direction of FFT formation",
     [(T, "hofmann1998_cys_ribose_ph7", "FFT", 20.0)], [(T, "hofmann1998_cys_ribose_ph3", "FFT", 20.0)], 0.5, True),
    ("zhou2023 furfural pH 8 / 6", "furfural supply vs pH (reads the pH model too)",
     [(Z, "zhou2023_arp_cys_120c_ph8", "FUR", 60.0)], [(Z, "zhou2023_arp_cys_120c_ph6", "FUR", 60.0)], 0.3, True),
    ("zhou2023 FFT pH 8 / 6", "FFT vs pH (reads the pH model too)",
     [(Z, "zhou2023_arp_cys_120c_ph8", "FFT", 60.0)], [(Z, "zhou2023_arp_cys_120c_ph6", "FFT", 60.0)], 0.3, True),
    ("zhou2023 MFT pH 7 / 6", "MFT's rising side vs pH",
     [(Z, "zhou2023_arp_cys_120c_ph7", "MFT", 60.0)], [(Z, "zhou2023_arp_cys_120c_ph6", "MFT", 60.0)], 0.3, True),
    ("zhou2023 MFT pH 8 / 7", "MFT's falling side vs pH",
     [(Z, "zhou2023_arp_cys_120c_ph8", "MFT", 60.0)], [(Z, "zhou2023_arp_cys_120c_ph7", "MFT", 60.0)], 0.3, True),
    ("davidek acetic 240 / 60 min, pH 8", "no acetic-acid sink (plateau at xylose exhaustion)",
     [(D, "davidek2008_xyl_gly_120c_ph80", "AA", 240.0)], [(D, "davidek2008_xyl_gly_120c_ph80", "AA", 60.0)],
     0.1, False),
    ("davidek acetic 240 / 60 min, pH 6", "no acetic-acid sink",
     [(D, "davidek2008_xyl_gly_120c_ph60", "AA", 240.0)], [(D, "davidek2008_xyl_gly_120c_ph60", "AA", 60.0)],
     0.1, False),
    ("davidek acetic pH 8 / 6 at 240 min", "base-favoured beta-cleavage of the 1-deoxypentosone",
     [(D, "davidek2008_xyl_gly_120c_ph80", "AA", 240.0)], [(D, "davidek2008_xyl_gly_120c_ph60", "AA", 240.0)],
     0.15, True),
    ("davidek xylose pH 8 / 6 at 30 min", "the pentose + amine entry's pH law",
     [(D, "davidek2008_xyl_gly_120c_ph80", "PENT", 30.0)], [(D, "davidek2008_xyl_gly_120c_ph60", "PENT", 30.0)],
     0.3, True),
    ("zhai MFT +Cys / TTCA, 120 C", "cysteine's stoichiometric role in MFT formation",
     [(T, "zhai2023jafc_ttcacys_ph7_120C", "MFT", t) for t in (40.0, 80.0, 120.0, 180.0)],
     [(T, "zhai2023jafc_ttca_ph7_120C", "MFT", t) for t in (40.0, 80.0, 120.0, 180.0)], 0.3, True),
    ("zhai MFT 140 / 120 C, +Cys (hold-out)", "the temperature dependence of MFT formation",
     [(T, "zhai2023jafc_ttcacys_ph7_140C", "MFT", t) for t in (40.0, 80.0)],
     [(T, "zhai2023jafc_ttcacys_ph7_120C", "MFT", t) for t in (40.0, 80.0)], 0.3, True),
    ("zhai FFT 140 / 120 C, +Cys (hold-out)", "the temperature dependence of FFT formation",
     [(T, "zhai2023jafc_ttcacys_ph7_140C", "FFT", t) for t in (40.0, 80.0)],
     [(T, "zhai2023jafc_ttcacys_ph7_120C", "FFT", t) for t in (40.0, 80.0)], 0.3, True),
)

#: v1 frozen 2026-10-31 (premortem #4; v2 the same night, two MFT-shape signatures added before any converged arm): each signature's mechanism group and its role in the joint fit -- "fitted"
#: (its rows are fit rows; for the Davidek rows, in every arm with the pentose entry), "held_out" (declared hold-out,
#: FIT_HOLDOUT_DECLARATION.md), "unfitted" (in no fit table, but read since round 2). Signatures are correlated
#: within a group, so progress is counted in groups: a group passes when all of its signatures pass. The isotope
#: route tests (Cerny 2003) are one more group, read from the fit's own check.
GROUP_OF: Dict[str, Tuple[str, str]] = {
    "schieberle MFT lag 60 / 720 min": ("MFT timing (an early source)", "fitted"),
    "schieberle FFT lag 60 / 720 min": ("FFT timing", "fitted"),
    "schieberle MFT plateau 720 / 360 min": ("MFT supply exhaustion (a reservoir)", "fitted"),
    "schieberle MFT early rise 60 / 30 min": ("MFT timing (an early source)", "fitted"),
    "hofmann MFT pH 7 / pH 3 (hold-out)": ("MFT pH direction", "held_out"),
    "zhou2023 MFT pH 7 / 6": ("MFT pH direction", "unfitted"),
    "zhou2023 MFT pH 8 / 7": ("MFT pH direction", "unfitted"),
    "hofmann FFT pH 7 / pH 3 (hold-out)": ("FFT and furfural pH direction", "held_out"),
    "zhou2023 FFT pH 8 / 6": ("FFT and furfural pH direction", "unfitted"),
    "zhou2023 furfural pH 8 / 6": ("FFT and furfural pH direction", "unfitted"),
    "davidek acetic 240 / 60 min, pH 8": ("no acetic-acid sink", "fitted"),
    "davidek acetic 240 / 60 min, pH 6": ("no acetic-acid sink", "fitted"),
    "davidek acetic pH 8 / 6 at 240 min": ("beta-cleavage pH law", "fitted"),
    "davidek xylose pH 8 / 6 at 30 min": ("pentose + amine entry pH law", "fitted"),
    "zhai MFT +Cys / TTCA, 120 C": ("cysteine's role in MFT formation", "fitted"),
    "zhai MFT 140 / 120 C, +Cys (hold-out)": ("temperature dependence", "held_out"),
    "zhai FFT 140 / 120 C, +Cys (hold-out)": ("temperature dependence", "held_out"),
}
CERNY_GROUP = "isotope routes (Cerny 2003)"
SIGNATURES_VERSION, SIGNATURES_SHA = "v2", "0ecbfb24f608"


def signatures_sha() -> str:
    """The list's fingerprint: an edit to any signature, tolerance, group or role changes it, so the list cannot move
    under a comparison without a new version."""
    import hashlib
    blob = json.dumps([list(s) for s in SIGNATURES] + [sorted(GROUP_OF.items())], sort_keys=True, default=list)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def _arm_of(stem: str) -> str:
    suffix = stem.split("_")[-1]
    by_lower = {k.lower(): k for k in joint.ARMS}
    return by_lower.get(suffix, "none")


def _configure(arm: str, strecker_share: float) -> None:
    joint.STRECKER_SHARE = float(strecker_share)
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in joint.ARMS[arm] else {}
    joint.EA_STATIC = {"k_dpo_nf": joint.LAU_NF_EA_KJ_MOL} if "dpo_nf_lau" in joint.ARMS[arm] else {}
    joint.TDP_EXIT_STATIC = joint.tdp_exit_static(joint.ARMS[arm])   # B7a R-T (the split and B7D's k are in theta)
    sulfur.PH_COUPLING = "inline"


def _measured() -> Tuple[Dict[str, Any], Dict[Point, float], Dict[Tuple[str, str], Dict[str, Any]]]:
    """The two-reader tables, every signature point's measured value, and one engine job per pot."""
    tables: Dict[str, Any] = {}
    measured: Dict[Point, float] = {}
    jobs: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for table, pot, species, t in sorted({p for _, _, num, den, _, _ in SIGNATURES for p in num + den}):
        if table not in tables:
            tables[table] = fit_rows.load(table)
        pots, rows = tables[table]
        hit = [r for r in rows if r.pot == pot and r.species == species and abs(r.t_min - t) < 1e-6]
        if not hit:
            raise SystemExit(f"no measured row for {(table, pot, species, t)}")
        measured[(table, pot, species, t)] = hit[0].value
        if (table, pot) not in jobs:
            jobs[(table, pot)] = next(j for j in probe.plan(pots, rows) if j["pot"] == pot)
    return tables, measured, jobs


def score(theta: Dict[str, float], strecker_share: float, tables, measured, jobs, strict: bool = True
          ) -> Tuple[List[Dict[str, Any]], Dict[str, bool], Dict[str, Any]]:
    """Every signature, every group's pass and the Cerny route checks at one theta (the arm already configured).
    ``strict``: a pot the engine cannot run stops the scoring (a fit's optimum); otherwise its signatures fail."""
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    state: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for (table, pot), job in jobs.items():
        j = dict(job, **trial, strecker_share=float(strecker_share))
        res = joint.run_task(j)
        if strict and "failed" in res:
            raise SystemExit(f"{pot}: {res['failed']}")
        state[(table, pot)] = {"job": j, "res": res}

    def predicted(point: Point) -> float:
        table, pot, species, t = point
        st = state[(table, pot)]
        if "failed" in st["res"]:
            return float("nan")
        row = next(r for r in tables[table][1] if r.pot == pot and r.species == species and abs(r.t_min - t) < 1e-6)
        return probe.predict(row, st["job"], st["res"])

    out: List[Dict[str, Any]] = []
    for name, mechanism, num, den, tol, sign in SIGNATURES:
        m = sum(math.log10(measured[a] / measured[b]) for a, b in zip(num, den)) / len(num)
        p_vals = [(predicted(a), predicted(b)) for a, b in zip(num, den)]
        p = (sum(math.log10(x / y) for x, y in p_vals) / len(p_vals)
             if all(x > 0 and y > 0 for x, y in p_vals) else float("nan"))
        ok = math.isfinite(p) and abs(p - m) <= tol and (not sign or (p == 0 or m == 0 or (p > 0) == (m > 0)))
        group, role = GROUP_OF[name]
        out.append({"signature": name, "mechanism": mechanism, "group": group, "role": role, "measured_log10": m,
                    "predicted_log10": p, "tolerance_dex": tol, "sign_must_match": sign, "pass": bool(ok)})
    cerny = joint.cerny_route_checks([theta[k] for k in keys], keys)
    groups: Dict[str, bool] = {}
    for s in out:
        groups[s["group"]] = groups.get(s["group"], True) and s["pass"]
    groups[CERNY_GROUP] = bool(cerny.get("verdict")) and "FAIL" not in str(cerny["verdict"])
    return out, groups, cerny


def chance(arm: str, n: int, seed: int, strecker_share: float = 0.042) -> int:
    """The chance baseline (premortem #4): pass rates per signature and per group when every coordinate of the arm is
    drawn uniformly within its declared box -- what a constant set chosen without the data passes. A signature or
    group that passes often by chance is weak evidence when it passes at a fit."""
    import numpy as np
    from src.kinetic_core import parallel
    _configure(arm, strecker_share)
    keys, _kinds, lower, upper, _start = joint.coordinates(False, joint.ARMS[arm])
    rng = np.random.default_rng(seed)
    draws = [dict(zip(keys, map(float, lower + (upper - lower) * rng.random(len(keys))))) for _ in range(n)]
    tables, measured, jobs = _measured()

    def one(theta):
        try:
            out, groups, _ = score(theta, strecker_share, tables, measured, jobs, strict=False)
            return {"signatures": {s["signature"]: s["pass"] for s in out}, "groups": groups}
        except Exception as exc:   # a draw the engine cannot run counts as failing every signature
            return {"error": repr(exc)[:200]}

    results = parallel.map_ordered(one, draws)
    ok = [r for r in results if "error" not in r]
    rate = {name: sum(r["signatures"][name] for r in ok) / max(1, len(ok)) for name, *_ in SIGNATURES}
    groups = sorted({g for g, _ in GROUP_OF.values()} | {CERNY_GROUP})
    g_rate = {g: sum(r["groups"].get(g, False) for r in ok) / max(1, len(ok)) for g in groups}
    g_count = [sum(r["groups"].values()) for r in ok]
    stem = data_paths.VALIDATION_DIR / f"kinetic_core_mechanism_signatures_chance_{arm.lower()}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/score_mechanism_signatures.py",
                   inputs=[ROOT / "src/kinetic_core/sulfur.py"] + [fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in sorted(tables)]),
               "arm": arm, "draws": n, "seed": seed, "runnable_draws": len(ok),
               "signatures_version": SIGNATURES_VERSION, "signatures_sha": SIGNATURES_SHA,
               "signature_pass_rate": rate, "group_pass_rate": g_rate,
               "groups_passed_distribution": {str(k): g_count.count(k) for k in sorted(set(g_count))},
               "errors_sample": [r["error"] for r in results if "error" in r][:5]}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    lines = [f"# Mechanism signatures: chance baseline, arm {arm}", "",
             f"*Generated by `scripts/generators/score_mechanism_signatures.py --chance {n} --arm {arm}` (seed {seed}). "
             "Every coordinate drawn uniformly within the arm's declared box; signature list "
             f"{SIGNATURES_VERSION} ({SIGNATURES_SHA}). {len(ok)} of {n} draws ran.*", "",
             "| group | chance pass rate |", "|---|---|"]
    lines += [f"| {g} | {r:.2f} |" for g, r in g_rate.items()]
    lines += ["", "| signature | role | chance pass rate |", "|---|---|---|"]
    lines += [f"| {name} | {GROUP_OF[name][1]} | {rate[name]:.2f} |" for name, *_ in SIGNATURES]
    lines += ["", "Groups passed per draw: " + ", ".join(f"{k}: {v}" for k, v in payload["groups_passed_distribution"].items())]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", help="the joint fit's stem suffix (s042, s042_r2b, s042_e, ...)")
    parser.add_argument("--chance", type=int, metavar="N", help="the chance baseline from N uniform draws (with --arm)")
    parser.add_argument("--arm", choices=tuple(joint.ARMS), help="the arm whose box the chance baseline draws from")
    parser.add_argument("--seed", type=int, default=20261031)
    args = parser.parse_args(argv)
    if signatures_sha() != SIGNATURES_SHA or set(GROUP_OF) != {s[0] for s in SIGNATURES}:
        raise SystemExit(f"the signature list changed (sha {signatures_sha()}, pinned {SIGNATURES_SHA}): "
                         "bump SIGNATURES_VERSION and re-pin, and never compare arms across versions")
    if args.chance:
        if not args.arm:
            raise SystemExit("--chance needs --arm")
        return chance(args.arm, args.chance, args.seed)
    if not args.fit:
        raise SystemExit("--fit or --chance")
    fit_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json"
    fit = json.loads(fit_path.read_text())
    # the arm and the Strecker share from the fit's own record (premortem 2026-10-31: a tagged re-fit's stem no longer
    # ends in its arm); the stem only for records written before the fit stored its arm
    theta, arm = fit["fit"]["theta"], fit.get("arm") or _arm_of(args.fit)
    _configure(arm, float(fit["strecker_share"]))
    if fit.get("ea_static") is not None:   # the barriers the fit held fixed (a profile point), from its record
        joint.EA_STATIC = dict(fit["ea_static"])
    if fit.get("tdp_exit") is not None:    # B7a R-T's constants, from the fit's record (before B7a: from the arm)
        joint.TDP_EXIT_STATIC = dict(fit["tdp_exit"])
    tables, measured, jobs = _measured()
    out, groups, cerny = score(theta, float(fit["strecker_share"]), tables, measured, jobs)
    stem = data_paths.VALIDATION_DIR / f"kinetic_core_mechanism_signatures_{args.fit}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/score_mechanism_signatures.py",
                   inputs=[fit_path, ROOT / "src/kinetic_core/sulfur.py"]
                   + [fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in sorted(tables)]),
               "fit": data_paths.rel(fit_path), "arm": arm,
               "signatures_version": SIGNATURES_VERSION, "signatures_sha": SIGNATURES_SHA,
               "fit_converged": not fit["fit"].get("budget_exhausted", False),
               "groups_passed": groups, "signatures": out, "cerny_route_checks": cerny}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    passed = sum(s["pass"] for s in out)
    g_passed = sum(groups.values())
    lines = [f"# Mechanism signatures: {args.fit} (arm {arm})", "",
             f"*Generated by `scripts/generators/score_mechanism_signatures.py` from `{data_paths.rel(fit_path)}`. "
             "Within-paper log10 ratios, measured (two-reader tables) vs predicted at the fit's optimum; tolerances "
             "declared before the arms were read (STRATEGY sec. 9).*", "",
             f"**{g_passed} of {len(groups)} mechanism groups pass; {passed} of {len(out)} signatures within "
             f"tolerance** (signature list {SIGNATURES_VERSION}, {SIGNATURES_SHA}). Cerny 2003 routes: "
             f"{cerny.get('verdict')}", "",
             ("" if payload["fit_converged"] else "**The fit stopped on its evaluation budget, not at an optimum: no "
              "verdict is read from this record (STATUS_AND_PLAN sec. 5, rule a).**"), "",
             "| group | passes |", "|---|---|"]
    lines += [f"| {g} | {'yes' if ok else 'NO'} |" for g, ok in groups.items()]
    lines += ["", "| signature | group | role | measured | predicted | tol | pass |", "|---|---|---|---|---|---|---|"]
    lines += [f"| {s['signature']} | {s['group']} | {s['role']} | {s['measured_log10']:+.2f} | "
              f"{s['predicted_log10']:+.2f} | {s['tolerance_dex']:.2f}{' (sign)' if s['sign_must_match'] else ''} | "
              f"{'yes' if s['pass'] else 'NO'} |" for s in out]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: {passed}/{len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
