#!/usr/bin/env python
"""
Build Wave B60 -- B51 REFITTED ON A CONVERGED pH TRAJECTORY (2026-10-28).

Pre-registered in ``results/validation/kinetic_core_b60_prereg.md`` before any B60 number existed. Candidate
(A) of that file, the kind-2 refit of a kind-1 correction: B51's objective and physics (B47's 52 corrected
rows; every ``ph_state.B51`` correction on; importing the B51 generator installs both) integrated on B60's
log-and-linear pH-node grid (``sulfur.PH_NODE_GRID``, the default since B60) with B60's flux budget. Same
network, free set, bands, sigmas, weighting, optimiser, budget and two-start protocol; start 0 is B51's optimum.

Usage (each start ~30-40 min; run them in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b60_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b60_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b60_fit.py --consolidate
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b51_fit as B51  # noqa: E402  (installs B47's rows; B51's physics on)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import ph_state, sulfur  # noqa: E402

B9, B23, B8, B47 = B51.B9, B51.B23, B51.B8, B51.B47

#: Candidate (A): B51's corrections on, the converged pH-node grid, for the whole process.
ph_state.B51 = ph_state.B51_ALL_ON
assert sulfur.PH_NODE_GRID == "log_and_linear" and sulfur.PH_FIXED_POINT == "relaxed_converged"

WAVE = "B60"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b60_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b60_fit_report.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b60_prereg.md"
B51_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b51_fit_report.json"

VERIFICATION = data_paths.VALIDATION_DIR / "kinetic_core_b60_verification.json"


def install_quick_setting() -> dict:
    """Amendment 1: the optimiser's quick setting is the one the verification measured converged
    (cheapest (ph_nodes, rtol) whose 52 observables lie within 0.03 dex of the careful setting)."""
    if getattr(B23.integrate_sulfur, "_b60_quick", None):
        return B23.integrate_sulfur._b60_quick
    choice = json.loads(VERIFICATION.read_text(encoding="utf-8"))["quick_setting"]["choice"]
    assert choice is not None, "the B60 verification found no converged quick setting; do not fit"
    B23.PH_NODES_QUICK = int(choice["ph_nodes"])
    original = B23.integrate_sulfur
    quick_rtol = float(choice["rtol"])

    def _quick_rtol(*args, **kwargs):
        if kwargs.get("rtol") == 1e-5:  # run_systems' quick tolerance
            kwargs["rtol"] = quick_rtol
        return original(*args, **kwargs)

    _quick_rtol._b60_quick = choice
    B23.integrate_sulfur = _quick_rtol
    return choice


#: Installed at import, so every consumer of this module (the Laplace, the profile, the fit-target ledger)
#: differentiates the objective the optimiser minimised, not the 9-node default.
QUICK_SETTING = install_quick_setting()


REMOVED_ROW_IDS = B51.REMOVED_ROW_IDS
ROW_CORRECTIONS = B51.ROW_CORRECTIONS
SYSTEM_CORRECTIONS = B51.SYSTEM_CORRECTIONS
FREE_KEYS = B51.FREE_KEYS
FREE_INDEX = B51.FREE_INDEX
ALL_KEYS = B51.ALL_KEYS
full_bounds = B51.full_bounds
residual_vector = B51.residual_vector


def incumbent_vector() -> np.ndarray:
    """B51's frozen optimum, clipped into the same bands."""
    lower, upper = full_bounds()
    return np.clip(frozen_vector(json.loads(B51_REPORT.read_text(encoding="utf-8"))), lower, upper)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B60 is pre-registered; write the prereg before running the fit"
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    quick_setting = QUICK_SETTING

    # Rebind B8's module-level paths and incumbent so its runner writes B47's artifacts.
    B8.MEMBER_DIR = MEMBER_DIR
    B8.OUT_FIT_REPORT = OUT_FIT_REPORT
    B8.incumbent_vector = incumbent_vector

    if args.consolidate:
        payload = B8.consolidate()
        payload["wave"] = "B60 -- B51 refitted on a converged pH trajectory (2026-10-28)"
        payload["generated_by"] = "scripts/generators/generate_kinetic_core_b60_fit.py"
        payload["prereg"] = data_paths.rel(PREREG)
        payload["objective"]["form"] = (
            "B51's objective (B47's 52 corrected rows, B51's pH-state corrections) on B60's "
            "log-and-linear pH-node grid and closed flux budget"
        )
        payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
        payload["objective"]["b47_removed_row_ids"] = list(REMOVED_ROW_IDS)
        payload["ph_state_corrections_b51"] = ph_state.B51.as_dict()
        payload["ph_numerics_b60"] = {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT,
                                      "relaxation": sulfur.PH_FIXED_POINT_RELAXATION,
                                      "quick_setting": quick_setting}
        payload["objective"]["b47_corrected_row_ids"] = sorted(ROW_CORRECTIONS)
        payload["objective"]["b47_corrected_systems"] = sorted(SYSTEM_CORRECTIONS)
        payload["objective"]["not_comparable_note"] = (
            "TOTAL COST IS COMPARABLE TO B47's (same 52 rows) but NOT TO B9's; the ship rule scores "
            "every arm on these 52 rows for the like-for-like number."
        )
        payload["objective"]["n_new_rows"] = 0
        payload["objective"]["new_row_ids"] = []
        payload["start_vector"] = "B51's frozen optimum (start 0); B8's perturbation protocol (start 1)"
        lower, upper = full_bounds()
        fr = payload["frozen_parameters"]
        x = np.array([fr["log10_k_ref_at_145C"][k] for k in B23.PARAM_ORDER]
                     + [fr["lumped_formation_Ea_kJ_mol"]]
                     + [fr["decay_Ea_kJ_mol"][f] for f in B23.DECAY_FAMILY_ORDER]
                     + [fr["ph_drift"]["acid_yield_per_sink_event"],
                        fr["ph_drift"]["arp_secondary_ammonium_pKa"]], dtype=float)
        active = []
        for i, key in enumerate(ALL_KEYS):
            width = float(upper[i] - lower[i])
            if (x[i] - lower[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "lower", "value": float(x[i])})
            elif (upper[i] - x[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "upper", "value": float(x[i])})
        payload["active_bounds"] = active
        OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
        print(f"wrote {OUT_FIT_REPORT}; active bounds: {[a['key'] for a in active]}")
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["ph_state_corrections_b51"] = ph_state.B51.as_dict()
    member["ph_numerics_b60"] = {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT,
                                 "quick_setting": quick_setting}
    member["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
    dest = MEMBER_DIR / f"b60_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
