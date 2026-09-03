#!/usr/bin/env python
"""
Build Wave B9 -- THE FIT / VALIDATE SPLIT (2026-09-03).

Pre-registered in ``results/validation/kinetic_core_b9_prereg.md`` before any B9
number existed. Owner's rule (2026-09-03): PRIMARY evidence -- rate constants,
activation energies, fed-intermediate yields, conversions, within-study ratios --
is FIT evidence; end-to-end concentrations measured in full precursor systems are
VALIDATION evidence. B8's objective violated that on eight rows: the Hofmann &
Schieberle 1998 Table 1 levels of FFT and MFT in the ribose / xylose / glucose /
fructose + cysteine pots at pH 5, which are ALSO four scored panel bundles
(``hofmann1998_{ribose,xylose,glucose,fructose}_cysteine_145C_20min_pH5``).

B9 is B8 with exactly those eight rows removed from the objective and NOTHING
else changed:

  * same network, same T-structure, same 23-coordinate free set, same declared
    bands (B8's ``full_bounds``), same W-HALF pH weighting, same optimiser and
    budget, same two-start protocol;
  * the objective drops from 62 to 54 rows: 46 rows shared with B8 plus the
    in-situ intermediate levels (NF, furan-2-aldehyde in the ribose pot, Hofmann
    T5) and every step-level row (fed intermediates T3 / T4 / T10, Kang / Zhai,
    Zhou, Whitfield, Cerny, van Seeventer, Feng) -- none of which a bundle scores;
  * start 0 is B8's frozen optimum (the incumbent), start 1 a perturbed vector.

What B9 gives back: the four Hofmann pH-5 bundles become OUT-OF-SAMPLE for the
kinetic core (the xylose bundle returns to the hold-out panel), and the sulfur
lane's absolute scale is anchored by yields per mmol fed rather than by the
levels it is then scored against. What it costs, if anything, the prereg says
how to read.

Usage (each start is ~25 min; run them in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b9_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b9_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b9_fit.py --consolidate
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
import generate_kinetic_core_b8_fit as B8  # noqa: E402  (installs B8's four rows into B23)
from src import data_paths  # noqa: E402

WAVE = "B9"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b9_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b9_prereg.md"
B8_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b8_fit_report.json"

#: The eight end-to-end LEVEL rows that are also scored panel bundles. Removed.
VALIDATION_ROW_IDS: Tuple[str, ...] = (
    "hofmann_ribose_FFT", "hofmann_ribose_MFT",
    "hofmann_xylose_FFT", "hofmann_xylose_MFT",
    "hofmann_glucose_FFT", "hofmann_glucose_MFT",
    "hofmann_fructose_FFT", "hofmann_fructose_MFT",
)

# Same interface as B8 for the derived-artifact generators (laplace, profile, fit targets).
FREE_KEYS = B8.FREE_KEYS
FREE_INDEX = B8.FREE_INDEX
ALL_KEYS = B8.ALL_KEYS
full_bounds = B8.full_bounds
residual_vector = B8.residual_vector


def install_b9_rows() -> None:
    """Remove the eight validation rows from the imported B8 objective, once."""
    before = len(B23.ACTIVE_FIT_ROWS)
    B23.ACTIVE_FIT_ROWS = tuple(r for r in B23.ACTIVE_FIT_ROWS if r["id"] not in VALIDATION_ROW_IDS)
    B23.FIT_ROWS = tuple(r for r in B23.FIT_ROWS if r["id"] not in VALIDATION_ROW_IDS)
    removed = before - len(B23.ACTIVE_FIT_ROWS)
    assert removed in (0, len(VALIDATION_ROW_IDS)), removed
    for row in B23.ACTIVE_FIT_ROWS:
        assert row.get("benchmark_id") not in {
            "hofmann1998_ribose_cysteine_145C_20min_pH5", "hofmann1998_xylose_cysteine_145C_20min_pH5",
            "hofmann1998_glucose_cysteine_145C_20min_pH5", "hofmann1998_fructose_cysteine_145C_20min_pH5",
        }, row["id"]


install_b9_rows()
assert len(B23.ACTIVE_FIT_ROWS) == 54, len(B23.ACTIVE_FIT_ROWS)


def incumbent_vector() -> np.ndarray:
    """B8's frozen optimum, clipped into the same bands (it already is)."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    report = json.loads(B8_FIT_REPORT.read_text(encoding="utf-8"))
    lower, upper = full_bounds()
    return np.clip(frozen_vector(report), lower, upper)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B9 is pre-registered; write the prereg before running the fit"
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)

    # Rebind B8's module-level paths and incumbent so its runner writes B9's artifacts.
    B8.MEMBER_DIR = MEMBER_DIR
    B8.OUT_FIT_REPORT = OUT_FIT_REPORT
    B8.incumbent_vector = incumbent_vector

    if args.consolidate:
        payload = B8.consolidate()
        payload["wave"] = "B9 -- the fit / validate split (owner rule 2026-09-03)"
        payload["generated_by"] = "scripts/generators/generate_kinetic_core_b9_fit.py"
        payload["prereg"] = data_paths.rel(PREREG)
        payload["objective"]["form"] = (
            "B8's objective with the eight Hofmann 1998 Table 1 LEVEL rows removed "
            "(ribose / xylose / glucose / fructose FFT and MFT at pH 5): 54 declared FIT rows"
        )
        payload["objective"]["removed_row_ids"] = list(VALIDATION_ROW_IDS)
        payload["objective"]["not_comparable_note"] = (
            "TOTAL COST IS NOT COMPARABLE TO B8's: eight rows fewer. "
            "`sum_r2_level_shared_with_b2_4` still sums the non-pH rows both waves scored, "
            "minus the eight removed, so compare it to B8's value recomputed over the same 46."
        )
        payload["objective"]["n_new_rows"] = 0
        payload["objective"]["new_row_ids"] = []
        payload["start_vector"] = "B8's frozen optimum (start 0); B8's perturbation protocol (start 1)"
        lower, upper = full_bounds()
        x = np.array([payload["frozen_parameters"]["log10_k_ref_at_145C"][k] for k in B23.PARAM_ORDER]
                     + [payload["frozen_parameters"]["lumped_formation_Ea_kJ_mol"]]
                     + [payload["frozen_parameters"]["decay_Ea_kJ_mol"][f] for f in B23.DECAY_FAMILY_ORDER]
                     + [payload["frozen_parameters"]["ph_drift"]["acid_yield_per_sink_event"],
                        payload["frozen_parameters"]["ph_drift"]["arp_secondary_ammonium_pKa"]], dtype=float)
        active = []
        for i, key in enumerate(ALL_KEYS):
            width = float(upper[i] - lower[i])
            if (x[i] - lower[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "lower", "value": float(x[i])})
            elif (upper[i] - x[i]) <= 1e-3 * width:
                active.append({"key": key, "bound": "upper", "value": float(x[i])})
        payload["active_bounds"] = active
        OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
        print(f"rewrote {OUT_FIT_REPORT} as B9; active bounds: {[a['key'] for a in active]}")
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["removed_row_ids"] = list(VALIDATION_ROW_IDS)
    dest = MEMBER_DIR / f"b9_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
