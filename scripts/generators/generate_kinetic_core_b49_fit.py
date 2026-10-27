#!/usr/bin/env python
"""
Build Wave B49 -- NORFURANEOL'S OTHER SULFIDE PRODUCTS (2026-10-27).

Pre-registered in ``results/validation/kinetic_core_b49_prereg.md`` before any B49 number existed.
B47's 52 corrected rows, plus two FIT rows from Whitfield & Mottram 1999 Table 2 at pH 4.5 (the
sulfur-heterocycle classes against MFT, fed norfuraneol + H2S or cysteine), plus ONE appended
coordinate: ``log10_k_nf_shet``, the constant of ``r_nf_shet`` (NF + HS- -> SHET; sulfur.py),
which proceeds through the hydrosulfide anion so its pH response is declared, not fitted. Its
barrier is the lumped formation barrier the vector already carries. Everything at pH 6.5
(Whitfield & Mottram 2001) is held out for the ship rule.

Usage (each start ~15-25 min; run them in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b49_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b49_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b49_fit.py --consolidate
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

import generate_kinetic_core_b47_fit as B47  # noqa: E402  (installs B47's 52 rows into B23)
from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    MEASURED_SULFUR, NF_SHET_BOUNDS_LOG10K, nf_shet_parameters, with_fitted_sulfur,
)

B23, B8, B9 = B47.B23, B47.B8, B47.B9

WAVE = "B49"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b49_prereg.md"
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b49_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b49_fit_report.json"

#: Declared mean molar masses of the Whitfield classes (g/mol), prereg sec. 3.
MW_DITHIOLANONES, MW_THIOPHENES, MW_THIOPHENONES, MW_MFT = 141.0, 112.0, 114.0, 114.16


def _shet_over_mft(dithiolanones: float, thiophenes: float, thiophenones: float, mft: float) -> float:
    shet = dithiolanones / MW_DITHIOLANONES + thiophenes / MW_THIOPHENES + thiophenones / MW_THIOPHENONES
    return shet / (mft / MW_MFT)


NEW_ROWS: Tuple[Dict[str, Any], ...] = (
    dict(id="whitfield_h2s_shet_over_MFT", system="whitfield_nf_h2s", kind="ratio",
         species="SHET", species_b="MFT", target=_shet_over_mft(443.0, 21.0, 39.0, 12.0), sigma_log=0.4,
         anchor="Whitfield & Mottram 1999 Table 2 (reprinted in Whitfield & Mottram 2001 Table 2), NF + H2S, pH 4.5: "
                "dithiolanones/dithianones 443, thiophenes 21, thiophenones 39, MFT 12 ug per 10 mg NF; molar with "
                "declared class masses 141/112/114"),
    dict(id="whitfield_cys_shet_over_MFT", system="whitfield_nf_cys", kind="ratio",
         species="SHET", species_b="MFT", target=_shet_over_mft(41.0, 94.0, 94.0, 15.0), sigma_log=0.4,
         anchor="Whitfield & Mottram 1999 Table 2, NF + cysteine, pH 4.5: 41, 94, 94; MFT 15 ug per 10 mg NF"),
)

_B8_ALL_KEYS, _B8_FREE_KEYS = B8.ALL_KEYS, B8.FREE_KEYS
_B23_BUILD_ORIGINAL = B23.build_parameters   # B9's builder, for the ship rule's B9 comparison
ALL_KEYS: Tuple[str, ...] = tuple(_B8_ALL_KEYS) + ("log10_k_nf_shet",)
SHET_SLOT = len(ALL_KEYS) - 1
FREE_KEYS: Tuple[str, ...] = tuple(_B8_FREE_KEYS) + ("log10_k_nf_shet",)
FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(B8.FREE_CLAUSE_OF),
    "log10_k_nf_shet": "B49: norfuraneol + HS- -> sulfur heterocycles; barrier shared with the lumped formation barrier",
}
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))
COORDINATE_OF_OVERRIDES: Dict[str, Dict[str, str]] = {
    "log10_k_nf_shet": {"block": "nf_shet_log10_k", "key": "k_nf_shet", "kind": "log10k"},
}
_B47_FULL_BOUNDS = B47.full_bounds


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B47_FULL_BOUNDS()
    return (np.append(np.array(lower, dtype=float), NF_SHET_BOUNDS_LOG10K[0]),
            np.append(np.array(upper, dtype=float), NF_SHET_BOUNDS_LOG10K[1]))


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    fitted, formation_ea, decay_ea, _drift = B23.unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B23.B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, formation_ea, decay_ea))
    parameters.update(nf_shet_parameters(k_nf_shet=10.0 ** float(x[SHET_SLOT]), ea_kj_mol=float(formation_ea)))
    return parameters


def incumbent_vector() -> np.ndarray:
    """B47's optimum with the new constant at the centre of its band."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    x47 = frozen_vector(json.loads(B47.OUT_FIT_REPORT.read_text(encoding="utf-8")))
    centre = 0.5 * (NF_SHET_BOUNDS_LOG10K[0] + NF_SHET_BOUNDS_LOG10K[1])
    lower, upper = full_bounds()
    return np.clip(np.append(x47, centre), lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


def install_b49_rows() -> None:
    if getattr(B23, "_B49_INSTALLED", False):
        return
    B23.ACTIVE_FIT_ROWS = tuple(B23.ACTIVE_FIT_ROWS) + NEW_ROWS
    B23.FIT_ROWS = tuple(B23.FIT_ROWS) + NEW_ROWS
    B23._B49_INSTALLED = True


def configure() -> None:
    install_b49_rows()
    B23.build_parameters = build_parameters
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = ALL_KEYS, FREE_KEYS, FREE_INDEX, FROZEN_KEYS
    B8.FREE_CLAUSE_OF = FREE_CLAUSE_OF
    B8.full_bounds, B8.incumbent_vector = full_bounds, incumbent_vector
    B8.MEMBER_DIR, B8.OUT_FIT_REPORT = MEMBER_DIR, OUT_FIT_REPORT


configure()
assert len(B23.ACTIVE_FIT_ROWS) == 54, len(B23.ACTIVE_FIT_ROWS)


def consolidate() -> Dict[str, Any]:
    payload = B8.consolidate()
    best = min(B8.load_members(), key=lambda m: m["cost"])
    x = np.array(best["x_full"], dtype=float)
    assert len(x) == len(ALL_KEYS), len(x)
    fr = payload["frozen_parameters"]
    fr["nf_shet_log10_k"] = {"k_nf_shet": float(x[SHET_SLOT])}
    payload["wave"] = f"{WAVE} -- norfuraneol's other sulfide products (2026-10-27)"
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b49_fit.py"
    payload["prereg"] = data_paths.rel(PREREG)
    payload["objective"]["form"] = ("B47's 52 corrected rows + 2 Whitfield 1999 heterocycle ratios at pH 4.5 (54 rows); "
                                    "one branch constant appended (24 free)")
    payload["objective"]["n_new_rows"] = len(NEW_ROWS)
    payload["objective"]["new_row_ids"] = [r["id"] for r in NEW_ROWS]
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(B47.REMOVED_ROW_IDS)
    payload["objective"]["b47_corrected_row_ids"] = sorted(B47.ROW_CORRECTIONS)
    payload["objective"]["n_free_parameters"] = len(FREE_KEYS)
    payload["free_set"] = {"n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS), "keys": list(FREE_KEYS),
                           "clause": FREE_CLAUSE_OF, "frozen_keys": list(FROZEN_KEYS)}
    payload["nf_shet_band_log10k"] = list(NF_SHET_BOUNDS_LOG10K)
    payload["start_vector"] = "B47's optimum with log10 k_nf_shet at its band centre (start 0); B8's perturbation protocol (start 1)"
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
    print(f"wrote {OUT_FIT_REPORT}: log10 k_nf_shet {x[SHET_SLOT]:.3f}; cost {best['cost']:.3f}; "
          f"active bounds: {[a['key'] for a in active]}")
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
    assert PREREG.exists(), "B49 is pre-registered; write the prereg before running the fit"
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["log10_k_nf_shet"] = float(member["x_full"][SHET_SLOT])
    dest = MEMBER_DIR / f"b49_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
