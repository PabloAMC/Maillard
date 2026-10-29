#!/usr/bin/env python
"""
Build Wave B61 -- NORFURANEOL'S HYDROXIDE-CATALYSED DECAY, MEASURED, REFITTED (2026-10-28).

Pre-registered in ``results/validation/kinetic_core_b61_prereg.md`` before any B61 number existed. Candidate
(A) of that file: B60's objective, physics and numerics (importing the B60 generator installs B51's rows and
pH-state corrections and B60's quick setting) with the B61 step NF -> fragments, hydroxide-catalysed, switched
on (``parameters_sulfur.NF_HYDROXIDE_DECAY``; Mikami 2017's 0.032 1/min at pH 8 / 100 C, the class-median
barrier). The step's constant is measured and NOT fitted. Same free set, bands, sigmas, weighting, optimiser,
budget and two-start protocol; start 0 is B60's optimum.

Usage (run the starts in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b61_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b61_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b61_fit.py --consolidate
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b60_fit as B60  # noqa: E402  (B51's rows and physics, B60's numerics and quick setting)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

B51, B9, B23, B8 = B60.B51, B60.B9, B60.B23, B60.B8

#: Candidate (A): B60 with the B61 step on, for the whole process.
ps.NF_HYDROXIDE_DECAY = dataclasses.replace(ps.NF_HYDROXIDE_DECAY, enabled=True)

WAVE = "B61"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b61_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b61_fit_report.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b61_prereg.md"
B60_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b60_fit_report.json"

QUICK_SETTING = B60.QUICK_SETTING
REMOVED_ROW_IDS = B60.REMOVED_ROW_IDS
ROW_CORRECTIONS = B60.ROW_CORRECTIONS
SYSTEM_CORRECTIONS = B60.SYSTEM_CORRECTIONS
FREE_KEYS = B60.FREE_KEYS
FREE_INDEX = B60.FREE_INDEX
ALL_KEYS = B60.ALL_KEYS
full_bounds = B60.full_bounds
residual_vector = B60.residual_vector


def incumbent_vector() -> np.ndarray:
    """B60's frozen optimum, clipped into the same bands."""
    lower, upper = full_bounds()
    return np.clip(frozen_vector(json.loads(B60_REPORT.read_text(encoding="utf-8"))), lower, upper)


def _step_record() -> dict:
    d = ps.NF_HYDROXIDE_DECAY
    return {"enabled": d.enabled, "k_ph8_per_min": d.k_ph8_per_min, "t_anchor_c": d.t_anchor_c,
            "ea_kj_mol": d.ea_kj_mol, "log10_k_anchor_ph5": d.log10_k_anchor_ph5, "source": d.source}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B61 is pre-registered; write the prereg before running the fit"
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)

    # Rebind B8's module-level paths and incumbent so its runner writes B61's artifacts.
    B8.MEMBER_DIR = MEMBER_DIR
    B8.OUT_FIT_REPORT = OUT_FIT_REPORT
    B8.incumbent_vector = incumbent_vector

    if args.consolidate:
        payload = B8.consolidate()
        payload["wave"] = "B61 -- norfuraneol's hydroxide-catalysed decay (Mikami 2017), refitted on B60's objective (2026-10-28)"
        payload["generated_by"] = "scripts/generators/generate_kinetic_core_b61_fit.py"
        payload["prereg"] = data_paths.rel(PREREG)
        payload["objective"]["form"] = (
            "B60's objective (B47's 52 corrected rows, B51's pH-state corrections, B60's numerics) with the "
            "B61 hydroxide-catalysed norfuraneol decay on (Mikami 2017, measured; not fitted)"
        )
        payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
        payload["objective"]["b47_removed_row_ids"] = list(REMOVED_ROW_IDS)
        payload["ph_state_corrections_b51"] = ph_state.B51.as_dict()
        payload["ph_numerics_b60"] = {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT,
                                      "relaxation": sulfur.PH_FIXED_POINT_RELAXATION,
                                      "quick_setting": QUICK_SETTING}
        payload["nf_hydroxide_decay_b61"] = _step_record()
        payload["objective"]["b47_corrected_row_ids"] = sorted(ROW_CORRECTIONS)
        payload["objective"]["b47_corrected_systems"] = sorted(SYSTEM_CORRECTIONS)
        payload["objective"]["not_comparable_note"] = (
            "TOTAL COST IS COMPARABLE TO B60's (same 52 rows, same numerics) but NOT TO B9's; the ship rule "
            "scores every arm on these 52 rows for the like-for-like number."
        )
        payload["objective"]["n_new_rows"] = 0
        payload["objective"]["new_row_ids"] = []
        payload["start_vector"] = "B60's frozen optimum (start 0); B8's perturbation protocol (start 1)"
        lower, upper = full_bounds()
        fr = payload["frozen_parameters"]
        x = frozen_vector(payload)
        active = []
        for i, key in enumerate(ALL_KEYS):
            width = float(upper[i] - lower[i])
            if (x[i] - lower[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "lower", "value": float(x[i])})
            elif (upper[i] - x[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "upper", "value": float(x[i])})
        payload["active_bounds"] = active
        assert fr["log10_k_ref_at_145C"]["k_nf_decay"] == json.loads(B60_REPORT.read_text())["frozen_parameters"][
            "log10_k_ref_at_145C"]["k_nf_decay"], "the neutral norfuraneol decay is frozen in B61"
        OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
        print(f"wrote {OUT_FIT_REPORT}; active bounds: {[a['key'] for a in active]}")
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["ph_state_corrections_b51"] = ph_state.B51.as_dict()
    member["ph_numerics_b60"] = {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT,
                                 "quick_setting": QUICK_SETTING}
    member["nf_hydroxide_decay_b61"] = _step_record()
    member["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
    dest = MEMBER_DIR / f"b61_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
