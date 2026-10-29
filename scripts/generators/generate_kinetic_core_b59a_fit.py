"""
WAVE B59a -- the Strecker sulfur released as H2S, refitted on each B64 baseline (2026-10-28).

B64's objective, on one of its thermal readings (the baseline), with ``STRECKER_SULFUR`` on at one of its two
measured shares (Kobayasi & Fujimaki 1965: 0.042, 0.080; both lower bounds). Start 0 is that reading's best B64
member. Nothing else changes: rows, bands, free set, numerics (inline pH coupling, compiled charge balance).
Pre-registration: results/validation/kinetic_core_b59a_prereg.md.

    python scripts/generators/generate_kinetic_core_b59a_fit.py --baseline ramp_plus_hold --share 0.042
    python scripts/generators/generate_kinetic_core_b59a_fit.py --baseline ramp_plus_hold --share 0.042 --cost-only
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

import generate_kinetic_core_b64_fit as B64  # noqa: E402  (B60's objective, inline pH, the thermal readings)
from src import data_paths  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402

B8 = B64.B8
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b59a_prereg.md"


def best_member(reading: str) -> Path:
    if reading == "isothermal":  # the shipped B60 optimum (a fit report, which frozen_vector reads)
        return data_paths.VALIDATION_DIR / "kinetic_core_b60_fit_report.json"
    files = sorted((data_paths.VALIDATION_DIR / f"kinetic_core_b64_{reading}_members").glob("*.json"))
    if not files:
        raise SystemExit(f"no B64 member for reading {reading!r}: run B64 first")
    return min(files, key=lambda f: json.loads(f.read_text())["cost"])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", choices=sorted(B64.READINGS), required=True)
    parser.add_argument("--share", type=float, choices=ps.STRECKER_SULFUR.alternatives, required=True)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--start-from", dest="start_from", type=Path,
                        help="continue from this member's x_full instead of the baseline's best B64 member")
    parser.add_argument("--cost-only", action="store_true", help="the baseline's vector with the switch on (candidate B)")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B59a is pre-registered; write the prereg before running the fit"
    B64.install_reading(args.baseline)
    ps.STRECKER_SULFUR = dataclasses.replace(ps.STRECKER_SULFUR, enabled=True, h2s_share=float(args.share))
    origin_path = args.start_from or best_member(args.baseline)
    lower, upper = B64.full_bounds()
    payload = json.loads(origin_path.read_text())
    x0 = payload["x_full"] if "x_full" in payload else B64.frozen_vector(payload)
    origin = np.clip(np.array(x0, dtype=float), lower, upper)
    tag = f"{args.baseline}_share{args.share:g}"
    record = {"baseline_reading": args.baseline, "h2s_share": args.share, "start_from": data_paths.rel(origin_path),
              "strecker_sulfur_source": ps.STRECKER_SULFUR.source}
    if args.cost_only:
        r = B64.B60.residual_vector(origin, True)
        print(json.dumps({**record, "cost_baseline_vector_switch_on": float(0.5 * np.sum(r ** 2))}))
        return 0
    member_dir = data_paths.VALIDATION_DIR / f"kinetic_core_b59a_{tag}_members"
    member_dir.mkdir(parents=True, exist_ok=True)
    B8.MEMBER_DIR = member_dir
    B8.OUT_FIT_REPORT = data_paths.VALIDATION_DIR / f"kinetic_core_b59a_{tag}_fit_report.json"
    B8.incumbent_vector = lambda: origin
    member = B8.fit_member(0, args.max_nfev, True, args.budget)
    member["wave"] = "B59a"
    member["b59a"] = record
    member["thermal_history_b64"] = {"reading": args.baseline}
    suffix = "_cont" if args.start_from is not None else ""
    dest = member_dir / f"b59a_{tag}_s0{suffix}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
