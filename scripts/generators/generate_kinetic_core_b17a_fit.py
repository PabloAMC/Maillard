#!/usr/bin/env python
"""
Build Wave B17 -- THE SINK STRUCTURE, variant (a): A SATURABLE THIOETHER SINK ON A POOL THE POT MAKES (2026-09-09).

Pre-registered in ``results/validation/kinetic_core_b17_prereg.md`` sec. 2 (a); run after variant (b)
failed T1, in the pre-registration's own order. B16's objective unchanged (B9's 54 rows + Schieberle
2000's seven within-study ratios at 100 C + Zhai 2021's three TTCA rows = 64 rows; every B9 band
kept, the thiol-sink ceiling at 102 kJ/mol) and ONE new coordinate appended to the vector, the way
B11 appended its oxygen consumers: ``log10_mele_site_yield``, the electrophile sites made per
deoxyosone decayed (``ch_mele_from_dpo`` / ``_tdp`` / ``_ddp`` in sulfur.py, rate = yield x
k_osone_decay with the carbonyl-sink family's barrier), which feed the MEASURED thioether channel
(Hofmann 2002's k_thioether, Stack 2018's K(T)) with a pool the pot itself makes and exhausts. The
release constant of variant (b) stays at its inert zero. 24 free. Start 0 is B9's optimum with the
yield at the centre of its band; start 1 is B8's perturbation protocol. Two starts, the
600-evaluation budget, quick mode for the search.

Usage (each start ~30 min; run both in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b17a_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b17a_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b17a_fit.py --consolidate
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
    MEASURED_SULFUR, MELE_SITE_YIELD_BOUNDS_LOG10, mele_site_parameters, with_fitted_sulfur,
)

WAVE = "B17a"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b17_prereg.md"
B9_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b17a_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b17a_fit_report.json"

# B16's bindings, captured before configure() rebinds B8's names to B17a's.
B16.configure(False)
_B16_FULL_BOUNDS, _B16_INCUMBENT_VECTOR = B8.full_bounds, B8.incumbent_vector   # B16's (bound by B16.configure)
_B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS = B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS
_B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT = B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT
_B23_BUILD_PARAMETERS_ORIGINAL = B23.build_parameters

# ===========================================================================
# THE VECTOR: 48 + 1
# ===========================================================================
ALL_KEYS: Tuple[str, ...] = tuple(_B8_ALL_KEYS) + ("log10_mele_site_yield",)
YIELD_SLOT = len(ALL_KEYS) - 1
FREE_KEYS: Tuple[str, ...] = tuple(_B8_FREE_KEYS) + ("log10_mele_site_yield",)
FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(_B8_FREE_CLAUSE_OF),
    "log10_mele_site_yield": "B17 variant (a): electrophile sites per deoxyosone decayed, feeding the measured thioether channel; barrier is the carbonyl-sink family's",
}
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))
COORDINATE_OF_OVERRIDES: Dict[str, Dict[str, str]] = {
    "log10_mele_site_yield": {"block": "mele_site_log10_yield", "key": "mele_site_yield", "kind": "log10k"},
}
assert len(FREE_KEYS) == 24 and len(ALL_KEYS) == 49


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B16_FULL_BOUNDS()
    return (np.append(np.array(lower, dtype=float), MELE_SITE_YIELD_BOUNDS_LOG10[0]),
            np.append(np.array(upper, dtype=float), MELE_SITE_YIELD_BOUNDS_LOG10[1]))


def yield_from_vector(x: np.ndarray) -> float:
    return 10.0 ** float(x[YIELD_SLOT])


def k_mele_site_from_vector(x: np.ndarray) -> Tuple[float, float]:
    """(k_mele_site at 145 C, its barrier): the yield times the fitted k_osone_decay, the carbonyl-sink Ea."""
    fitted, _formation_ea, decay_ea, _drift = B23.unpack(x)
    return yield_from_vector(x) * 10.0 ** float(fitted["k_osone_decay"]), float(decay_ea["carbonyl_sink"])


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    fitted, formation_ea, decay_ea, _drift = B23.unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B23.B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, formation_ea, decay_ea))
    k_site, ea_site = k_mele_site_from_vector(x)
    parameters.update(mele_site_parameters(k_mele_site=k_site, ea_kj_mol=ea_site))
    return parameters


def incumbent_vector() -> np.ndarray:
    """B9's frozen optimum with the site yield at the centre of its band."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    x9 = frozen_vector(json.loads(B9_FIT_REPORT.read_text(encoding="utf-8")))
    assert len(x9) == B23.N_K + B23.N_EXTRA, len(x9)
    centre = 0.5 * (MELE_SITE_YIELD_BOUNDS_LOG10[0] + MELE_SITE_YIELD_BOUNDS_LOG10[1])
    lower, upper = full_bounds()
    return np.clip(np.append(x9, centre), lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


def restore() -> None:
    B23.build_parameters = _B23_BUILD_PARAMETERS_ORIGINAL
    B8.full_bounds, B8.incumbent_vector = _B16_FULL_BOUNDS, _B16_INCUMBENT_VECTOR
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = _B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS
    B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT = _B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT
    B16.restore()


def configure() -> None:
    B16.configure(False)
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
    k_site, ea_site = k_mele_site_from_vector(x)
    fr["mele_site_log10_yield"] = {"mele_site_yield": float(x[YIELD_SLOT])}
    fr["mele_site"] = {"yield_sites_per_osone_decayed": yield_from_vector(x), "k_mele_site_per_min_at_145C": k_site,
                       "ea_kj_mol_carbonyl_sink_family": ea_site}
    payload["wave"] = f"{WAVE} -- the sink structure, variant (a): a saturable thioether sink on a pool the pot makes"
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b17a_fit.py"
    payload["prereg"] = data_paths.rel(PREREG)
    payload["declaration"] = "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 29"
    payload["objective"]["form"] = "B16's 64 rows (B9's 54 + 7 Schieberle 2000 ratios + 3 Zhai 2021 TTCA rows); one site yield appended (24 free)"
    payload["objective"]["n_new_rows"] = 0
    payload["objective"]["new_row_ids"] = []
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    payload["objective"]["n_free_parameters"] = len(FREE_KEYS)
    payload["objective"]["n_frozen"] = len(FROZEN_KEYS)
    payload["free_set"] = {"n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS), "keys": list(FREE_KEYS),
                           "clause": FREE_CLAUSE_OF, "frozen_keys": list(FROZEN_KEYS)}
    payload["site_yield_band_log10"] = list(MELE_SITE_YIELD_BOUNDS_LOG10)
    payload["thiol_sink_ceiling_kj_mol"] = float(full_bounds()[1][ALL_KEYS.index("Ea_decay_thiol_sink")])
    payload["start_vector"] = "B9's frozen optimum with the site yield at the centre of its band (start 0); B8's perturbation protocol (start 1)"
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
    print(f"rewrote {OUT_FIT_REPORT} as {WAVE}: log10 site yield {x[YIELD_SLOT]:.3f}; active bounds: {[a['key'] for a in active]}")
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
    assert PREREG.exists(), "B17 (both variants) is pre-registered; write the prereg before running the fit"
    configure()
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    assert len(B23.ACTIVE_FIT_ROWS) == 64, len(B23.ACTIVE_FIT_ROWS)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["log10_mele_site_yield"] = float(member["x_full"][YIELD_SLOT])
    dest = MEMBER_DIR / f"b17a_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
