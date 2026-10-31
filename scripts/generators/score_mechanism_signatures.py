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


def _arm_of(stem: str) -> str:
    suffix = stem.split("_")[-1]
    by_lower = {k.lower(): k for k in joint.ARMS}
    return by_lower.get(suffix, "none")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", required=True, help="the joint fit's stem suffix (s042, s042_r2b, s042_e, ...)")
    args = parser.parse_args(argv)
    fit_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json"
    fit = json.loads(fit_path.read_text())
    # the arm and the Strecker share from the fit's own record (premortem 2026-10-31: a tagged re-fit's stem no longer
    # ends in its arm); the stem only for records written before the fit stored its arm
    theta, arm = fit["fit"]["theta"], fit.get("arm") or _arm_of(args.fit)
    joint.STRECKER_SHARE = float(fit["strecker_share"])
    keys = tuple(theta)
    joint.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in joint.ARMS[arm] else {}
    joint.EA_STATIC = {"k_dpo_nf": joint.LAU_NF_EA_KJ_MOL} if "dpo_nf_lau" in joint.ARMS[arm] else {}
    trial = joint.trial([theta[k] for k in keys], keys)
    sulfur.PH_COUPLING = "inline"

    tables: Dict[str, Any] = {}
    measured: Dict[Point, float] = {}
    needed = {p for _, _, num, den, _, _ in SIGNATURES for p in num + den}
    jobs: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for table, pot, species, t in needed:
        if table not in tables:
            tables[table] = fit_rows.load(table)
        pots, rows = tables[table]
        hit = [r for r in rows if r.pot == pot and r.species == species and abs(r.t_min - t) < 1e-6]
        if not hit:
            raise SystemExit(f"no measured row for {(table, pot, species, t)}")
        measured[(table, pot, species, t)] = hit[0].value
        if (table, pot) not in jobs:
            jobs[(table, pot)] = next(j for j in probe.plan(pots, rows) if j["pot"] == pot)

    predicted_state: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for (table, pot), job in jobs.items():
        j = dict(job, **trial, strecker_share=float(fit["strecker_share"]))
        res = joint.run_task(j)
        if "failed" in res:
            raise SystemExit(f"{pot}: {res['failed']}")
        predicted_state[(table, pot)] = {"job": j, "res": res}

    def predicted(point: Point) -> float:
        table, pot, species, t = point
        st = predicted_state[(table, pot)]
        row = next(r for r in tables[table][1] if r.pot == pot and r.species == species and abs(r.t_min - t) < 1e-6)
        return probe.predict(row, st["job"], st["res"])

    out: List[Dict[str, Any]] = []
    for name, mechanism, num, den, tol, sign in SIGNATURES:
        m = sum(math.log10(measured[a] / measured[b]) for a, b in zip(num, den)) / len(num)
        p_vals = [(predicted(a), predicted(b)) for a, b in zip(num, den)]
        p = (sum(math.log10(x / y) for x, y in p_vals) / len(p_vals)
             if all(x > 0 and y > 0 for x, y in p_vals) else float("nan"))
        ok = math.isfinite(p) and abs(p - m) <= tol and (not sign or (p == 0 or m == 0 or (p > 0) == (m > 0)))
        out.append({"signature": name, "mechanism": mechanism, "measured_log10": m, "predicted_log10": p,
                    "tolerance_dex": tol, "sign_must_match": sign, "pass": bool(ok)})

    cerny = joint.cerny_route_checks([theta[k] for k in keys], keys)
    stem = data_paths.VALIDATION_DIR / f"kinetic_core_mechanism_signatures_{args.fit}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/score_mechanism_signatures.py",
                   inputs=[fit_path, ROOT / "src/kinetic_core/sulfur.py"]
                   + [fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in sorted(tables)]),
               "fit": data_paths.rel(fit_path), "arm": arm, "signatures": out, "cerny_route_checks": cerny}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    passed = sum(s["pass"] for s in out)
    lines = [f"# Mechanism signatures: {args.fit} (arm {arm})", "",
             f"*Generated by `scripts/generators/score_mechanism_signatures.py` from `{data_paths.rel(fit_path)}`. "
             "Within-paper log10 ratios, measured (two-reader tables) vs predicted at the fit's optimum; tolerances "
             "declared before the arms were read (STRATEGY sec. 9).*", "",
             f"**{passed} of {len(out)} signatures within tolerance.** Cerny 2003 routes: {cerny.get('verdict')}", "",
             "| signature | mechanism | measured | predicted | tol | pass |", "|---|---|---|---|---|---|"]
    lines += [f"| {s['signature']} | {s['mechanism']} | {s['measured_log10']:+.2f} | {s['predicted_log10']:+.2f} | "
              f"{s['tolerance_dex']:.2f}{' (sign)' if s['sign_must_match'] else ''} | {'yes' if s['pass'] else 'NO'} |"
              for s in out]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: {passed}/{len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
