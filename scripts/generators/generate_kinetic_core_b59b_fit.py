"""
WAVE B59b -- de Roos's route (TTCA reservoir, cysteine Amadori compound) refitted on each B64 reading (2026-10-29).

B60's objective on one B64 thermal reading, with ``CYS_AMADORI_ROUTE`` on (its constants measured or anchored on Cao's
pot; the direct cysteine-catalysed pentose enolisations retired) and, optionally, B59a-r's Strecker sulfur at one of
its two transferred shares. Start 0 is the reading's best B64 member (isothermal: B60's optimum). Nothing else
changes. Pre-registration: results/validation/kinetic_core_b59b_prereg.md.

Free set (Amendment 2): a kind-2 refit refits only the coordinates the data identify (WAVES.md, "Evidence decides,
scores diagnose"), here the 16 that B60's Laplace marks identified; the other 7 stay at the start vector's values.
``--free-set all`` reproduces the withdrawn first runs (all 23 free) as a diagnosis only.

    python scripts/generators/generate_kinetic_core_b59b_fit.py --baseline isothermal
    python scripts/generators/generate_kinetic_core_b59b_fit.py --baseline isothermal --strecker-share 0.042
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b64_fit as B64  # noqa: E402
from generate_kinetic_core_b59a_fit import best_member  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402

B8 = B64.B8
LAPLACE = data_paths.VALIDATION_DIR / "kinetic_core_b60_laplace_covariance.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b59b_prereg.md"


def install(strecker_share=None) -> dict:
    ps.CYS_AMADORI_ROUTE = dataclasses.replace(ps.CYS_AMADORI_ROUTE, enabled=True)
    if strecker_share is not None:
        ps.STRECKER_SULFUR = dataclasses.replace(ps.STRECKER_SULFUR, enabled=True, h2s_share=float(strecker_share))
    return {"cys_amadori_route": dataclasses.asdict(ps.CYS_AMADORI_ROUTE),
            "strecker_sulfur": dataclasses.asdict(ps.STRECKER_SULFUR) if strecker_share is not None else None}


def identified_keys() -> tuple:
    """The coordinates B60's Laplace identifies (its own rule: sigma within the kind's threshold, not bound-pinned)."""
    laplace = json.loads(LAPLACE.read_text(encoding="utf-8"))
    keys = tuple(k for k, ok in zip(laplace["keys"], laplace["identified"]) if ok)
    assert set(keys) <= set(B8.FREE_KEYS), set(keys) - set(B8.FREE_KEYS)
    return keys


def restrict_free_set(keys: tuple) -> None:
    B8.FREE_KEYS = tuple(k for k in B8.FREE_KEYS if k in keys)
    B8.FREE_INDEX = tuple(B8.ALL_KEYS.index(k) for k in B8.FREE_KEYS)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", choices=sorted(B64.READINGS), required=True)
    parser.add_argument("--strecker-share", dest="share", type=float, choices=ps.STRECKER_SULFUR.alternatives)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--budget", type=int, default=1200)
    parser.add_argument("--cost-only", action="store_true", help="the baseline's vector with the route on")
    parser.add_argument("--free-set", dest="free_set", choices=("identified", "all"), default="identified")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B59b is pre-registered; write the prereg before running the fit"
    B64.install_reading(args.baseline)
    record = install(args.share)
    origin_path = best_member(args.baseline)
    payload = json.loads(origin_path.read_text())
    lower, upper = B64.full_bounds()
    origin = np.clip(np.array(payload["x_full"] if "x_full" in payload else B64.frozen_vector(payload), dtype=float),
                     lower, upper)
    tag = (f"{args.baseline}" + (f"_strecker{args.share:g}" if args.share is not None else "")
           + ("_all_free" if args.free_set == "all" else ""))
    if args.free_set == "identified":
        restrict_free_set(identified_keys())
    record.update({"baseline_reading": args.baseline, "start_from": data_paths.rel(origin_path),
                   "free_set": args.free_set, "free_keys": list(B8.FREE_KEYS)})
    if args.cost_only:
        r = B64.B60.residual_vector(origin, True)
        print(json.dumps({"tag": tag, "cost_baseline_vector_route_on": float(0.5 * np.sum(r ** 2))}))
        return 0
    member_dir = data_paths.VALIDATION_DIR / f"kinetic_core_b59b_{tag}_members"
    member_dir.mkdir(parents=True, exist_ok=True)
    B8.MEMBER_DIR = member_dir
    B8.OUT_FIT_REPORT = data_paths.VALIDATION_DIR / f"kinetic_core_b59b_{tag}_fit_report.json"
    B8.incumbent_vector = lambda: origin
    member = B8.fit_member(0, args.max_nfev, True, args.budget)
    member["wave"] = "B59b"
    member["b59b"] = record
    dest = member_dir / f"b59b_{tag}_s0.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
