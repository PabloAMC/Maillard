#!/usr/bin/env python
"""
Build Wave B27 -- THE OXIDANT THE THREE REFUSED SINK WAVES WERE NEVER GIVEN (2026-09-11).

Pre-registered in ``results/validation/kinetic_core_b27_prereg.md`` (sections 1-8 on 2026-09-09; sections
9 and 10 on 2026-09-10/11, before this fit was started). B16's 64-row objective, plus ONE new row (Whitfield
1999's MFT disulfide share, a lower bound), the three charge corrections of prereg sec. 4 installed at
import and restored on exit, and ONE new coordinate appended to the vector the way B11 and B25 appended
theirs: ``log10_phi_redox``, the fraction of mercaptoketone-forming events that deliver one oxidant
equivalent (parameters_sulfur.apply_dicarbonyl_redox; sulfur.py ch_redox_mp3p). 24 free.

Start 0 is B9's optimum with phi at the centre of its band; start 1 is B8's perturbation protocol.

Usage (each start ~30 min; run both in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b27_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b27_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b27_fit.py --consolidate
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

import generate_kinetic_core_b16_fit as B16  # noqa: E402  (installs B16's ten rows into B23 at import)
import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
import generate_kinetic_core_b8_fit as B8  # noqa: E402
import generate_kinetic_core_b9_fit as B9  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    DICARBONYL_REDOX_BOUNDS_LOG10_YIELD, MEASURED_SULFUR, apply_dicarbonyl_redox, with_fitted_sulfur,
)
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

WAVE = "B27"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b27_prereg.md"
B9_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b27_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b27_fit_report.json"

# B16's bindings, captured before configure() rebinds B8's names to B27's.
B16.configure(False)
_B16_FULL_BOUNDS, _B16_INCUMBENT_VECTOR = B8.full_bounds, B8.incumbent_vector
_B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS = B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS
_B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT = B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT
_B23_BUILD_PARAMETERS_ORIGINAL = B23.build_parameters

# ===========================================================================
# THE VECTOR: 48 + 1
# ===========================================================================
ALL_KEYS: Tuple[str, ...] = tuple(_B8_ALL_KEYS) + ("log10_phi_redox",)
K_SLOT = len(ALL_KEYS) - 1
FREE_KEYS: Tuple[str, ...] = tuple(_B8_FREE_KEYS) + ("log10_phi_redox",)
FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(_B8_FREE_CLAUSE_OF),
    "log10_phi_redox": ("B27 (the dicarbonyl redox couple): log10 of the fraction of mercaptoketone-forming events that "
                        "deliver one oxidant equivalent; a branching ratio, no barrier of its own; ceiling 1 is physical"),
}
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))
COORDINATE_OF_OVERRIDES: Dict[str, Dict[str, str]] = {
    "log10_phi_redox": {"block": "dicarbonyl_redox", "key": "log10_ox_yield_per_mercaptoketone", "kind": "log10k"},
}
assert len(FREE_KEYS) == 24 and len(ALL_KEYS) == 49

# ===========================================================================
# THE THREE CHARGE CORRECTIONS AND THE ONE NEW ROW (prereg sec. 4), installed at import
# ===========================================================================
#: The printed buffer, replacing the ASSUMED Hofmann spec the frozen generator carried for this pot.
BUFFER_WHITFIELD = BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=True,
    source=("Whitfield & Mottram 1999 Methods: 'cysteine in 0.5 M phosphate buffer at pH 4.5'; a saturated H2S "
            "solution in the same buffer; whitfield1999_extraction.md sec. 3"),
)
#: mmol/L in the final 2 mL, from the printed masses (11.4 mg NF / 114.10; 12.1 mg Cys / 121.16; 6.6 mg H2S / 34.08).
NF_MMOL_L, CYS_MMOL_L, H2S_MMOL_L = 50.0, 50.0, 97.0
CORRECTED_SYSTEMS: Dict[str, Dict[str, Any]] = {
    "whitfield_nf_cys": dict(
        initial={"NF": NF_MMOL_L, "Cys": CYS_MMOL_L}, t_c=140.0, minutes=60.0, ph=4.5, buffer=BUFFER_WHITFIELD,
        anchor="Whitfield & Mottram 1999 Table 1, NF 50 + cysteine 50 mmol/L, 0.5 M phosphate pH 4.5 (B27 charge correction)"),
    "whitfield_nf_h2s": dict(
        initial={"NF": NF_MMOL_L, "H2S": H2S_MMOL_L}, t_c=140.0, minutes=60.0, ph=4.5, buffer=BUFFER_WHITFIELD,
        anchor="Whitfield & Mottram 1999 Table 1, NF 50 + H2S ~97 mmol/L, 0.5 M phosphate pH 4.5 (B27 charge correction)"),
}
CORRECTED_ROWS: Dict[str, Dict[str, Any]] = {
    # free MFT 15 ug per 10 mg NF = 0.150 mol %; the six MFT-bearing disulfides carry a further ~8.0 ug of MFT
    # equivalent, so the MFT actually MADE is 0.230 mol % -- and the model species MFT is the FREE pool now that
    # the dimer is carried separately, so the row must score the total against MFT + 2 x dimer.
    "whitfield_nf_cys_MFT": dict(
        id="whitfield_nf_cys_MFT", system="whitfield_nf_cys", kind="molpct_total",
        species_terms={"MFT": 1, "MFTD": 2}, basis=NF_MMOL_L, target=0.230, sigma_log=0.5,
        anchor="Whitfield & Mottram 1999 T1: NF + cysteine, pH 4.5, 140 C, free MFT 0.150 mol % + ~0.080 mol % "
               "as six MFT-bearing disulfides = 0.230 mol % TOTAL (whitfield1999_extraction.md sec. 2; DHS, "
               "response factors ASSUMED 1)"),
    "whitfield_nf_h2s_MFT": dict(
        id="whitfield_nf_h2s_MFT", system="whitfield_nf_h2s", kind="molpct",
        species="MFT", basis=NF_MMOL_L, target=0.120, sigma_log=0.5,
        anchor="Whitfield & Mottram 1999 T1: NF + H2S 1:2, pH 4.5, 0.120 mol % free MFT over the printed 50 mmol/L basis",
        note="two labs, two methods, agreement within ~1.6x on the H2S channel (against Hofmann's 0.19 mol% at 145 C)"),
}
NEW_ROW: Dict[str, Any] = dict(
    id="whitfield_mft_disulfide_share_floor", system="whitfield_nf_cys", kind="floor",
    species_terms={"numerator": {"MFTD": 2}, "denominator": {"MFT": 1, "MFTD": 2}},
    target=0.35, sigma_log=0.3,
    anchor="Whitfield & Mottram 1999 T1 (cysteine column): six MFT-bearing disulfides (47, 48, 54, 55, 62, 64) carry "
           "~8.0 ug of MFT equivalent against 15 ug free MFT, so ~35 % of the MFT made is in a disulfide. A LOWER BOUND: "
           "disulfides are far less volatile than thiols under a 60 C dynamic headspace and unit response factors are "
           "assumed, so over-prediction is not a failure (kinetic_core_b27_prereg.md sec. 4)",
)

_ORIGINAL_SYSTEMS = {name: dict(B23.SYSTEMS[name]) for name in CORRECTED_SYSTEMS}
_ORIGINAL_ROWS = {r["id"]: r for r in B23.FIT_ROWS if r["id"] in CORRECTED_ROWS}


def _swap_rows(rows, replacements: Dict[str, Dict[str, Any]]):
    return tuple(replacements.get(r["id"], r) for r in rows)


def install_b27_corrections(on: bool) -> None:
    if on:
        for name, spec in CORRECTED_SYSTEMS.items():
            B23.SYSTEMS[name] = dict(spec)
        B23.ACTIVE_FIT_ROWS = _swap_rows(B23.ACTIVE_FIT_ROWS, CORRECTED_ROWS)
        B23.FIT_ROWS = _swap_rows(B23.FIT_ROWS, CORRECTED_ROWS)
        if NEW_ROW["id"] not in {r["id"] for r in B23.ACTIVE_FIT_ROWS}:
            B23.ACTIVE_FIT_ROWS = tuple(B23.ACTIVE_FIT_ROWS) + (NEW_ROW,)
            B23.FIT_ROWS = tuple(B23.FIT_ROWS) + (NEW_ROW,)
    else:
        for name, spec in _ORIGINAL_SYSTEMS.items():
            B23.SYSTEMS[name] = dict(spec)
        B23.ACTIVE_FIT_ROWS = tuple(r for r in _swap_rows(B23.ACTIVE_FIT_ROWS, _ORIGINAL_ROWS) if r["id"] != NEW_ROW["id"])
        B23.FIT_ROWS = tuple(r for r in _swap_rows(B23.FIT_ROWS, _ORIGINAL_ROWS) if r["id"] != NEW_ROW["id"])


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B16_FULL_BOUNDS()
    return (np.append(np.array(lower, dtype=float), [DICARBONYL_REDOX_BOUNDS_LOG10_YIELD[0]]),
            np.append(np.array(upper, dtype=float), [DICARBONYL_REDOX_BOUNDS_LOG10_YIELD[1]]))


def phi_from_vector(x: np.ndarray) -> float:
    return 10.0 ** float(x[K_SLOT])


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    fitted, formation_ea, decay_ea, _drift = B23.unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B23.B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, formation_ea, decay_ea))
    apply_dicarbonyl_redox(parameters, float(x[K_SLOT]))
    return parameters


def incumbent_vector() -> np.ndarray:
    """B9's frozen optimum with phi at the centre of its band."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    x9 = frozen_vector(json.loads(B9_FIT_REPORT.read_text(encoding="utf-8")))
    assert len(x9) == B23.N_K + B23.N_EXTRA, len(x9)
    centre = [0.5 * (DICARBONYL_REDOX_BOUNDS_LOG10_YIELD[0] + DICARBONYL_REDOX_BOUNDS_LOG10_YIELD[1])]
    lower, upper = full_bounds()
    return np.clip(np.append(x9, centre), lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


def restore() -> None:
    B23.build_parameters = _B23_BUILD_PARAMETERS_ORIGINAL
    install_b27_corrections(False)
    B8.full_bounds, B8.incumbent_vector = _B16_FULL_BOUNDS, _B16_INCUMBENT_VECTOR
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = _B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS
    B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT = _B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT
    B16.restore()


def configure() -> None:
    B16.configure(False)
    install_b27_corrections(True)
    B23.build_parameters = build_parameters
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = ALL_KEYS, FREE_KEYS, FREE_INDEX, FROZEN_KEYS
    B8.FREE_CLAUSE_OF = FREE_CLAUSE_OF
    B8.full_bounds, B8.incumbent_vector = full_bounds, incumbent_vector
    B8.MEMBER_DIR, B8.OUT_FIT_REPORT = MEMBER_DIR, OUT_FIT_REPORT


configure()


def consolidate() -> Dict[str, Any]:
    payload = B8.consolidate()
    best = min(B8.load_members(), key=lambda m: m["cost"])
    x = np.array(best["x_full"], dtype=float)
    assert len(x) == len(ALL_KEYS), len(x)
    fr = payload["frozen_parameters"]
    fr["dicarbonyl_redox"] = {"log10_ox_yield_per_mercaptoketone": float(x[K_SLOT])}
    fr["dicarbonyl_redox_values"] = {"phi_ox_per_mercaptoketone": phi_from_vector(x)}
    payload["wave"] = f"{WAVE} -- the dicarbonyl redox couple: oxidant as a co-product of the mercaptoketone flux"
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b27_fit.py"
    payload["prereg"] = data_paths.rel(PREREG)
    payload["declaration"] = "docs/reference/FIT_HOLDOUT_DECLARATION.md (amendment on shipping)"
    payload["objective"]["form"] = ("B16's 64 rows with the three Whitfield charge corrections + Whitfield 1999's MFT disulfide "
                                    "share as a lower bound (65 rows); phi appended (24 free)")
    payload["objective"]["n_new_rows"] = 1
    payload["objective"]["new_row_ids"] = [NEW_ROW["id"]]
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    payload["objective"]["n_free_parameters"] = len(FREE_KEYS)
    payload["objective"]["n_frozen"] = len(FROZEN_KEYS)
    payload["free_set"] = {"n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS), "keys": list(FREE_KEYS),
                           "clause": FREE_CLAUSE_OF, "frozen_keys": list(FROZEN_KEYS)}
    payload["redox_band_log10_yield"] = list(DICARBONYL_REDOX_BOUNDS_LOG10_YIELD)
    payload["charge_corrections"] = {
        "whitfield_nf_cys": {"before": _ORIGINAL_SYSTEMS["whitfield_nf_cys"]["initial"], "after": CORRECTED_SYSTEMS["whitfield_nf_cys"]["initial"]},
        "whitfield_nf_h2s": {"before": _ORIGINAL_SYSTEMS["whitfield_nf_h2s"]["initial"], "after": CORRECTED_SYSTEMS["whitfield_nf_h2s"]["initial"]},
        "whitfield_nf_cys_MFT": {"before": {"target": _ORIGINAL_ROWS["whitfield_nf_cys_MFT"]["target"], "basis": _ORIGINAL_ROWS["whitfield_nf_cys_MFT"]["basis"], "kind": "molpct (free MFT)"},
                                 "after": {"target": 0.230, "basis": NF_MMOL_L, "kind": "molpct_total (MFT + 2 x dimer)"}},
        "buffer": "0.5 M phosphate pH 4.5, printed (was the ASSUMED Hofmann spec)",
    }
    payload["start_vector"] = "B9's frozen optimum with phi at the centre of its band (start 0); B8's perturbation protocol (start 1)"
    lower, upper = full_bounds()
    active = []
    for i, key in enumerate(ALL_KEYS):
        width = float(upper[i] - lower[i])
        if (x[i] - lower[i]) <= 1e-3 * width:
            active.append({"key": key, "bound": "lower", "value": float(x[i])})
        elif (upper[i] - x[i]) <= 1e-3 * width:
            active.append({"key": key, "bound": "upper", "value": float(x[i])})
    payload["active_bounds"] = active
    OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
    print(f"rewrote {OUT_FIT_REPORT} as {WAVE}: log10 phi {x[K_SLOT]:.3f} (phi {phi_from_vector(x):.3g}); active bounds: {[a['key'] for a in active]}")
    return payload


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B27 is pre-registered; write the prereg before running the fit"
    configure()
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    assert len(B23.ACTIVE_FIT_ROWS) == 65, len(B23.ACTIVE_FIT_ROWS)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["log10_phi_redox"] = float(member["x_full"][K_SLOT])
    dest = MEMBER_DIR / f"b27_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
