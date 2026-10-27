#!/usr/bin/env python
"""
Build Wave B47 -- THE SULFUR OBJECTIVE, CORRECTED AGAINST ITS PDFS (2026-10-27).

Pre-registered in ``results/validation/kinetic_core_b47_prereg.md`` before any B47
number existed. On 2026-10-27 the 54 rows of B9's objective were checked against the
source PDFs for the first time; six rows and two systems did not say what their papers
say. B47 is B9 with exactly those corrected and NOTHING else changed:

  * targets: the Whitfield mercaptoketone ratio (16.3 summed three compounds by mass;
    2-mercapto-3-pentanone alone is 4.99 molar), the Kang free-cysteine ceiling (the
    120 C curve, 12.1 %, not the 140 C 16.3 %), and five Hofmann 1998 mol % rows taken
    from the printed micrograms instead of the one-figure percentages;
  * removed: the Cerny isomer split (Cerny 2007 prints 0 for the numerator at pH 5) and
    the Yaghmur FFT-share ceiling (a 65 C / 15 h measurement applied to a 145 C pot);
  * systems: Whitfield's fed-norfuraneol pots at the paper's 50 mM (and ~97 mM H2S),
    with the mol % basis following; the Whitfield, van Seeventer and Cerny 2007b
    buffers declared from their sources (numerically unchanged, 0.5 M phosphate);
  * same network, T-structure, 23-coordinate free set, bands, sigmas, W-HALF pH
    weighting, optimiser, budget and two-start protocol; start 0 is B9's optimum.

Usage (each start is ~25 min; run them in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b47_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b47_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b47_fit.py --consolidate
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b9_fit as B9  # noqa: E402  (installs B9's 54 rows into B23)
from src import data_paths  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

B23 = B9.B23
B8 = B9.B8

WAVE = "B47"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b47_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b47_fit_report.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b47_prereg.md"

#: Rows whose condition or value has no counterpart in the source. Removed.
REMOVED_ROW_IDS: Tuple[str, ...] = ("cerny_isomer_split", "yaghmur_fft_share_ceiling")

#: Row -> (field, new value, anchor appended). The evidence is the prereg's table.
ROW_CORRECTIONS: Dict[str, Dict[str, Any]] = {
    "whitfield_mercaptoketone_over_MFT": dict(
        target=4.99,
        anchor="Whitfield & Mottram 1999 Table 1, cysteine system: 2-mercapto-3-pentanone 77.5 against "
               "MFT 15 ug per 10 mg norfuraneol = (77.5/118.20)/(15/114.16) = 4.99 molar. B47: B9's 16.3 "
               "was (92 + 74.5 + 77.5)/15, three mercaptoketones by mass"),
    "kang_free_cys_yield_ceiling": dict(
        target=0.121,
        anchor="Kang 2026 SI Fig. S3 (digitised 2026-10-27): the 120 C curve peaks at ~1.21 mM of 10 mM "
               "TTCA (12.1 %) at 60 min. B47: B9's 16.3 % was the 140 C curve"),
    "whitfield_nf_cys_MFT": dict(basis=50.0),
    "whitfield_nf_h2s_MFT": dict(basis=50.0),
    "fed_ribose_h2s_MFT": dict(target=0.0132, anchor="Hofmann 1998 Table 4: 15.1 ug -> 0.0132 mol% (printed 0.01)"),
    "fed_c2c3_MFT_pH3": dict(target=0.0136, anchor="Hofmann 1998 Table 8, pH 3.0: 15.5 ug -> 0.0136 mol% (printed 0.01)"),
    "fed_thiamine_MFT": dict(target=0.0072, anchor="Hofmann 1998 Table 10: 8.2 ug -> 0.0072 mol% (printed 0.01)"),
    "fed_tdp_h2s_FFT": dict(target=0.069, anchor="Hofmann 1998 Table 3: 78.6 ug -> 0.069 mol% (printed 0.08)"),
    "fed_nf_cys_MFT": dict(target=0.0445, anchor="Hofmann 1998 Table 10: 50.8 ug -> 0.0445 mol% (printed 0.05)"),
}


def _source_buffer(source: str) -> BufferSpec:
    return BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source=source)


SYSTEM_CORRECTIONS: Dict[str, Dict[str, Any]] = {
    "whitfield_nf_cys": dict(
        initial={"NF": 50.0, "Cys": 50.0},
        buffer=_source_buffer("Whitfield & Mottram 1999 Methods: 0.5 M phosphate, pH 4.5; "
                              "1 mL 0.1 M norfuraneol + 1 mL 0.1 M cysteine")),
    "whitfield_nf_h2s": dict(
        initial={"NF": 50.0, "H2S": 97.0},
        buffer=_source_buffer("Whitfield & Mottram 1999 Methods: 0.5 M phosphate, pH 4.5; "
                              "1 mL 0.1 M norfuraneol + 1 mL ~0.2 M H2S")),
    "vanseeventer_130C": dict(buffer=_source_buffer("van Seeventer 2001 Methods: 0.5 M phosphate, pH 5")),
    "cerny_ternary": dict(buffer=_source_buffer("Cerny & Briffod 2007b Methods: 0.5 M phosphate, pH 5.00")),
}

# Same interface as B9 for the derived-artifact generators (laplace, profile, fit targets).
FREE_KEYS = B9.FREE_KEYS
FREE_INDEX = B9.FREE_INDEX
ALL_KEYS = B9.ALL_KEYS
full_bounds = B9.full_bounds
residual_vector = B9.residual_vector


def _correct(row: Dict[str, Any]) -> Dict[str, Any]:
    fix = ROW_CORRECTIONS.get(row["id"])
    if not fix:
        return row
    new = dict(row)
    for key, value in fix.items():
        new[key] = f"{value} || B9: {row.get('anchor', '')}" if key == "anchor" else value
    new["b47_corrected_from"] = {k: row.get(k) for k in fix if k != "anchor"}
    return new


def install_b47_rows() -> None:
    """Patch B9's 54 rows and four systems in place, once."""
    if getattr(B23, "_B47_INSTALLED", False):
        return
    B23.ACTIVE_FIT_ROWS = tuple(_correct(r) for r in B23.ACTIVE_FIT_ROWS if r["id"] not in REMOVED_ROW_IDS)
    B23.FIT_ROWS = tuple(_correct(r) for r in B23.FIT_ROWS if r["id"] not in REMOVED_ROW_IDS)
    for name, fix in SYSTEM_CORRECTIONS.items():
        B23.SYSTEMS[name] = {**B23.SYSTEMS[name], **fix}
    ids = {r["id"] for r in B23.ACTIVE_FIT_ROWS}
    assert set(ROW_CORRECTIONS) <= ids, set(ROW_CORRECTIONS) - ids
    B23._B47_INSTALLED = True


install_b47_rows()
assert len(B23.ACTIVE_FIT_ROWS) == 52, len(B23.ACTIVE_FIT_ROWS)


def incumbent_vector() -> np.ndarray:
    """B9's frozen optimum, clipped into the same bands (it already is)."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    report = json.loads(B9.OUT_FIT_REPORT.read_text(encoding="utf-8"))
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
    assert PREREG.exists(), "B47 is pre-registered; write the prereg before running the fit"
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)

    # Rebind B8's module-level paths and incumbent so its runner writes B47's artifacts.
    B8.MEMBER_DIR = MEMBER_DIR
    B8.OUT_FIT_REPORT = OUT_FIT_REPORT
    B8.incumbent_vector = incumbent_vector

    if args.consolidate:
        payload = B8.consolidate()
        payload["wave"] = "B47 -- the sulfur objective, corrected against its PDFs (2026-10-27)"
        payload["generated_by"] = "scripts/generators/generate_kinetic_core_b47_fit.py"
        payload["prereg"] = data_paths.rel(PREREG)
        payload["objective"]["form"] = (
            "B9's objective with six rows and four systems corrected against the source PDFs "
            "and two rows removed: 52 declared FIT rows"
        )
        payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
        payload["objective"]["b47_removed_row_ids"] = list(REMOVED_ROW_IDS)
        payload["objective"]["b47_corrected_row_ids"] = sorted(ROW_CORRECTIONS)
        payload["objective"]["b47_corrected_systems"] = sorted(SYSTEM_CORRECTIONS)
        payload["objective"]["not_comparable_note"] = (
            "TOTAL COST IS NOT COMPARABLE TO B9's: two rows fewer and seven targets changed. "
            "The ship rule scores B9's vector on B47's 52 rows for the like-for-like number."
        )
        payload["objective"]["n_new_rows"] = 0
        payload["objective"]["new_row_ids"] = []
        payload["start_vector"] = "B9's frozen optimum (start 0); B8's perturbation protocol (start 1)"
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
    member["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
    dest = MEMBER_DIR / f"b47_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
