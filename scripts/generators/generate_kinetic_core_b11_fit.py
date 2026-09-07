#!/usr/bin/env python
"""
Build Wave B11 -- OXYGEN AS AN INPUT (2026-09-07).

Pre-registered in ``results/validation/kinetic_core_b11_prereg.md`` (sec. 9 amendments).
B11 is B9's objective (54 rows) and free set (23) plus the TWO oxygen consumers of the
two-pool oxygen state (`sulfur.py`: ox_supply, ch_cys_ox, ch_red_ox_*), appended to the
vector as log10 coordinates 48 and 49 with the declared bands
`parameters_sulfur.OXYGEN_BOUNDS_LOG10K`. Every fit system is charged with a headspace
reservoir `OXR` from the vessel table below (Hofmann 1998's volumes are stated; every
other system gets the declared default and is marked as identifying neither consumer).

Usage (two starts in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b11_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b11_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b11_fit.py --consolidate
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

import generate_kinetic_core_b2_3_fit as B23  # noqa: E402
import generate_kinetic_core_b8_fit as B8  # noqa: E402
import generate_kinetic_core_b9_fit as B9  # noqa: E402  (removes the eight Hofmann level rows)
from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    K_OX_SUPPLY_PER_UNIT_MIN, MEASURED_SULFUR, OX_RESERVOIR_DEFAULT_UNITS, OX_SAT_MMOL_L,
    OXYGEN_BOUNDS_LOG10K, OXYGEN_FITTED_KEYS, oxygen_parameters, with_fitted_sulfur,
)

WAVE = "B11"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b11_prereg.md"
B9_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b11_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b11_fit_report.json"

# ===========================================================================
# 1. THE VESSEL TABLE OF THE FIT SYSTEMS (prereg sec. 9.3), in reservoir units
# ===========================================================================
_HOF_100 = (0.8706 + 0.027) / 0.100 / OX_SAT_MMOL_L      # 100 mL in a 200 mL autoclave, air
_HOF_50 = (1.306 + 0.0135) / 0.050 / OX_SAT_MMOL_L       # 50 mL fed pots, same autoclave
VESSEL_TABLE: Dict[str, Tuple[float, str]] = {
    "hofmann_pentose_pH5": (_HOF_100, "Hofmann 1998: 100 mL in a 200 mL autoclave (stated)"),
    "hofmann_glucose_pH5": (_HOF_100, "Hofmann 1998: 100 mL in a 200 mL autoclave (stated)"),
    "hofmann_fructose_pH5": (_HOF_100, "Hofmann 1998: 100 mL in a 200 mL autoclave (stated)"),
}
_FED = ("fed_ribose_h2s", "fed_tdp_h2s", "fed_furfural_h2s", "fed_nf_h2s", "fed_nf_cys",
        "fed_c2c3", "fed_c2c3_pH3", "fed_c2c3_pH7", "fed_thiamine", "fed_mgo_h2s_1to1", "fed_mgo_h2s_1to2")
for _name in _FED:
    VESSEL_TABLE[_name] = (_HOF_50, "Hofmann 1998: 50 mL fed pot in the 200 mL autoclave (stated)")
DEFAULT_BASIS = "vessel volumes unstated in the source: the declared default reservoir"


def reservoir_for(system: str) -> Tuple[float, str, bool]:
    if system in VESSEL_TABLE:
        units, basis = VESSEL_TABLE[system]
        return units, basis, True
    return OX_RESERVOIR_DEFAULT_UNITS, DEFAULT_BASIS, False


# ===========================================================================
# 2. THE VECTOR: 48 + 2
# ===========================================================================
_B8_FULL_BOUNDS = B8.full_bounds
_B8_INCUMBENT_VECTOR = B8.incumbent_vector
_B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS = B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS
_B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT = B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT
_B23_BUILD_PARAMETERS_ORIGINAL = B23.build_parameters
_B23_INITIALS_ORIGINAL = {name: dict(spec["initial"]) for name, spec in B23.SYSTEMS.items()}

ALL_KEYS: Tuple[str, ...] = tuple(B8.ALL_KEYS) + ("log10_k_cys_ox", "log10_k_red_ox")
CYS_SLOT, RED_SLOT = len(ALL_KEYS) - 2, len(ALL_KEYS) - 1
FREE_KEYS: Tuple[str, ...] = tuple(B8.FREE_KEYS) + ("log10_k_cys_ox", "log10_k_red_ox")
FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(B8.FREE_CLAUSE_OF),
    "log10_k_cys_ox": "R7 oxygen consumer (B11): cysteine autoxidation; expected unidentified",
    "log10_k_red_ox": "R7 oxygen consumer (B11): the reductone pool; expected unidentified",
}
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))
COORDINATE_OF_OVERRIDES: Dict[str, Dict[str, str]] = {
    "log10_k_cys_ox": {"block": "oxygen_log10_k", "key": "k_cys_ox", "kind": "log10k"},
    "log10_k_red_ox": {"block": "oxygen_log10_k", "key": "k_red_ox", "kind": "log10k"},
}
assert len(FREE_KEYS) == 25 and len(ALL_KEYS) == 50


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B8_FULL_BOUNDS()
    lo = list(np.array(lower, dtype=float)) + [OXYGEN_BOUNDS_LOG10K[k][0] for k in OXYGEN_FITTED_KEYS]
    hi = list(np.array(upper, dtype=float)) + [OXYGEN_BOUNDS_LOG10K[k][1] for k in OXYGEN_FITTED_KEYS]
    return np.array(lo), np.array(hi)


def oxygen_from_vector(x: np.ndarray) -> Dict[str, float]:
    return {"k_cys_ox": 10.0 ** float(x[CYS_SLOT]), "k_red_ox": 10.0 ** float(x[RED_SLOT])}


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    fitted, formation_ea, decay_ea, _drift = B23.unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B23.B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, formation_ea, decay_ea))
    parameters.update(oxygen_parameters(**oxygen_from_vector(x)))
    return parameters


def incumbent_vector() -> np.ndarray:
    """B9's frozen optimum with both consumers at the geometric centre of their bands."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    x9 = frozen_vector(json.loads(B9_FIT_REPORT.read_text(encoding="utf-8")))
    assert len(x9) == B23.N_K + B23.N_EXTRA, len(x9)
    centres = [0.5 * (OXYGEN_BOUNDS_LOG10K[k][0] + OXYGEN_BOUNDS_LOG10K[k][1]) for k in OXYGEN_FITTED_KEYS]
    lower, upper = full_bounds()
    return np.clip(np.append(x9, centres), lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


def install_reservoirs() -> None:
    for name, spec in B23.SYSTEMS.items():
        units, _basis, _stated = reservoir_for(name)
        spec["initial"] = {**_B23_INITIALS_ORIGINAL.get(name, spec["initial"]), "OXR": units, "OXV": 0.0}


def restore() -> None:
    B23.build_parameters = _B23_BUILD_PARAMETERS_ORIGINAL
    for name, spec in B23.SYSTEMS.items():
        if name in _B23_INITIALS_ORIGINAL:
            spec["initial"] = dict(_B23_INITIALS_ORIGINAL[name])
    B8.full_bounds, B8.incumbent_vector = _B8_FULL_BOUNDS, _B8_INCUMBENT_VECTOR
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = _B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS
    B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT = _B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT


def configure() -> None:
    install_reservoirs()
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
    oxy = oxygen_from_vector(x)
    fr = payload["frozen_parameters"]
    fr["oxygen"] = {**oxy, "k_ox_supply": K_OX_SUPPLY_PER_UNIT_MIN, "ox_sat_mmol_l": OX_SAT_MMOL_L,
                    "reservoir_default_units": OX_RESERVOIR_DEFAULT_UNITS}
    fr["oxygen_log10_k"] = {"k_cys_ox": float(x[CYS_SLOT]), "k_red_ox": float(x[RED_SLOT])}
    payload["wave"] = f"{WAVE} -- oxygen as an input (two-pool state; consumers fitted, expected unidentified)"
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b11_fit.py"
    payload["prereg"] = data_paths.rel(PREREG)
    payload["declaration"] = "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 23"
    payload["objective"]["form"] = "B9's 54 rows, every system charged with its headspace reservoir; two oxygen consumers appended (25 free)"
    payload["objective"]["n_new_rows"] = 0
    payload["objective"]["new_row_ids"] = []
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    payload["objective"]["n_free_parameters"] = len(FREE_KEYS)
    payload["objective"]["n_frozen"] = len(FROZEN_KEYS)
    payload["free_set"] = {"n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS), "keys": list(FREE_KEYS),
                           "clause": FREE_CLAUSE_OF, "frozen_keys": list(FROZEN_KEYS)}
    payload["vessel_table"] = {
        name: {"reservoir_units": reservoir_for(name)[0], "basis": reservoir_for(name)[1], "stated": reservoir_for(name)[2]}
        for name in B23.SYSTEMS
    }
    payload["oxygen_bands_log10k"] = {k: list(v) for k, v in OXYGEN_BOUNDS_LOG10K.items()}
    payload["start_vector"] = "B9's frozen optimum with both consumers at the centre of their log10 bands (start 0); the perturbation protocol (start 1)"
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
    print(f"rewrote {OUT_FIT_REPORT} as {WAVE}: oxygen {oxy}; active bounds: {[a['key'] for a in active]}")
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
    assert PREREG.exists(), "B11 is pre-registered; write the prereg before running the fit"
    configure()
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    assert len(B23.ACTIVE_FIT_ROWS) == 54, len(B23.ACTIVE_FIT_ROWS)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["oxygen"] = oxygen_from_vector(np.array(member["x_full"], dtype=float))
    dest = MEMBER_DIR / f"b11_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
