#!/usr/bin/env python
"""
Build Wave B10 -- THE TEMPERATURE STRUCTURE OF THE SULFUR LANE (2026-09-06).

Pre-registered in ``results/validation/kinetic_core_b10_prereg.md`` before any B10
number existed. Owner's instruction (2026-09-06): "take the most sensible decision
long term." What B10 changes, and nothing else:

  * ONE formation barrier becomes TWO, by route (`parameters_sulfur.FORMATION_ROUTE_OF`):
    `Ea_sugar_trunk` takes the lumped barrier's slot in the vector, `Ea_thiol_assembly`
    is appended as the 49th coordinate. Both are FREE, with bands narrowed by the
    prefactor rule around sourced centres (Zhang 2026 k16; Chan & Reineccius 1994).
  * SIX within-study FOLD rows from Yiltirak 2026's time-compensated ladder enter the
    objective (owner's rule: within-study ratios are primary evidence; the eight LEVELS
    stay on the hold-out panel and become fit-adjacent). Four new systems carry the pot.
  * The ambient oxidant is charged consistently in the engine (B11 prereg sec. 2.1);
    the fit systems already carried it, so the objective is unchanged by that.
  * Everything else is B9: same 54 rows, same 23 free coordinates (now 25), same
    weighting, optimiser, budget and two-start protocol; start 0 is B9's optimum with
    the two route barriers at their band centres.

`--without-yiltirak` runs the pre-registered leave-Yiltirak-out variant (T5): same
free set, the six fold rows and four systems NOT installed, artifacts under their own
names, so the reader can see what Yiltirak alone taught.

Usage (each start is ~25-40 min; run the four in parallel, then consolidate):
    python scripts/generators/generate_kinetic_core_b10_fit.py --start 0
    python scripts/generators/generate_kinetic_core_b10_fit.py --start 1
    python scripts/generators/generate_kinetic_core_b10_fit.py --without-yiltirak --start 0
    python scripts/generators/generate_kinetic_core_b10_fit.py --without-yiltirak --start 1
    python scripts/generators/generate_kinetic_core_b10_fit.py --consolidate
    python scripts/generators/generate_kinetic_core_b10_fit.py --without-yiltirak --consolidate
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
import generate_kinetic_core_b8_fit as B8  # noqa: E402  (installs B8's four rows into B23)
import generate_kinetic_core_b9_fit as B9  # noqa: E402  (removes the eight Hofmann level rows)
from src import data_paths  # noqa: E402
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    FORMATION_EA_BOUNDS_BY_ROUTE,
    FORMATION_EA_PRIOR_CENTRE,
    FORMATION_EA_PRIOR_SOURCE,
    FORMATION_ROUTE_OF,
    FORMATION_ROUTES,
    MEASURED_SULFUR,
    OX_AMBIENT_MMOL_L,
    with_fitted_sulfur,
)
from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

WAVE = "B10"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b10_prereg.md"
B9_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"

# Rebound by `configure()` for the leave-Yiltirak-out variant.
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b10_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b10_fit_report.json"
WITH_YILTIRAK = True

assert OX_AMBIENT_MMOL_L == B23.OX_AMBIENT_MMOL_L, "engine and fit must charge the same ambient oxidant"

# ===========================================================================
# 1. THE FOUR YILTIRAK SYSTEMS AND THE SIX FOLD ROWS
# ===========================================================================
# data/lit/extraction_dossiers/yiltirak2026_extraction.md sec. 2 (Methods, verbatim)
# and sec. 3 (Table S3, re-typed from the docx on disk). 25 mM ribose + 25 mM
# cysteine in 0.5 M potassium phosphate pH 5.5 (bench, tap water); 3 mL in a 20 mL
# PTFE-capped tube under air, stirred; four "equivalent" cooks. The vessel is
# recorded on the bundles (R1) and NOT modelled here: the ambient oxidant is
# what every other fit system carries, and B11 is the wave that reads the vessel.
BUFFER_YILTIRAK = BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=True,
    source=("Yiltirak et al. 2026 sec. 2.1 / 2.3, read from data/articles/Yiltirak2026.pdf: "
            "'Potassium phosphate buffer (0.5 M, pH 5.5) was made ... in tap water'; 'ribose "
            "(25 mM) and cysteine (25 mM) in potassium phosphate buffer (0.5 M, pH 5.5)'"))

YILTIRAK_RUNGS: Tuple[Tuple[str, float, float], ...] = (
    ("yiltirak_100C_4h", 100.0, 240.0),
    ("yiltirak_110C_2h", 110.0, 120.0),
    ("yiltirak_120C_1h", 120.0, 60.0),
    ("yiltirak_130C_30min", 130.0, 30.0),
)

B10_SYSTEMS: Dict[str, Dict[str, Any]] = {
    name: dict(
        initial={"PENT": 25.0, "Cys": 25.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=t_c, minutes=minutes, ph=5.5, buffer=BUFFER_YILTIRAK,
        anchor=(f"Yiltirak 2026 sec. 2.4: {t_c:.0f} C / {minutes:.0f} min, model (iv) buffer arm; "
                "yiltirak2026_extraction.md sec. 2"),
    )
    for name, t_c, minutes in YILTIRAK_RUNGS
}

#: Table S3, Buffer rows, ug/L (mean +/- SD, n = 3): MFT 6.88/3.29/2.4/1.71, FFT 1.28/1.46/1.68/1.62.
TABLE_S3_BUFFER: Dict[str, Dict[str, float]] = {
    "MFT": {"yiltirak_100C_4h": 6.88, "yiltirak_110C_2h": 3.29, "yiltirak_120C_1h": 2.4, "yiltirak_130C_30min": 1.71},
    "FFT": {"yiltirak_100C_4h": 1.28, "yiltirak_110C_2h": 1.46, "yiltirak_120C_1h": 1.68, "yiltirak_130C_30min": 1.62},
}
_S3_QUOTE = ("Table S3 (mmc1.docx), Buffer rows verbatim: '6.88+/-0.98aA 1.28+/-0.08cB' (100 C), "
             "'3.29+/-0.06bB 1.46+/-0.19abB' (110 C), '2.4+/-0.16bcB 1.68+/-0.21aA' (120 C), "
             "'1.71+/-0.12cB 1.62+/-0.11abB' (130 C); columns MFT | FFT, ug/L")
_BUNDLE = {
    "yiltirak_100C_4h": "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026",
    "yiltirak_110C_2h": "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026",
    "yiltirak_120C_1h": "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026",
    "yiltirak_130C_30min": "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026",
}
_COMPOUND = {"MFT": "2-Methyl-3-furanthiol (MFT)", "FFT": "2-Furfurylthiol (FFT)"}

#: sigma_log 0.10: the printed SDs give 0.02-0.06 dex per rung and the come-up time of a
#: 3 mL tube in a heating block is not reported; 0.10 leaves room for it (prereg sec. 4).
SIGMA_LOG_FOLD = 0.10


def _fold_rows() -> Tuple[Dict[str, Any], ...]:
    rows = []
    names = [n for n, _t, _m in YILTIRAK_RUNGS]
    for species in ("MFT", "FFT"):
        for lo, hi in zip(names[:-1], names[1:]):
            target = TABLE_S3_BUFFER[species][hi] / TABLE_S3_BUFFER[species][lo]
            rows.append(dict(
                id=f"yiltirak_{species}_fold_{hi.split('_')[1]}_over_{lo.split('_')[1]}",
                system=hi, system_b=lo, kind="cross_system_ratio", species=species,
                target=target, sigma_log=SIGMA_LOG_FOLD,
                anchor=(f"Yiltirak 2026 {_S3_QUOTE}; {species} {TABLE_S3_BUFFER[species][lo]} -> "
                        f"{TABLE_S3_BUFFER[species][hi]} ug/L over {lo} -> {hi}"),
                # the fit-target index lists BOTH bundles the fold reads, so their
                # levels leave the strict out-of-sample count (prereg sec. 2, 8)
                benchmark_id=_BUNDLE[hi], benchmark_compound=_COMPOUND[species],
                benchmark_id_b=_BUNDLE[lo],
                note=("WITHIN-STUDY FOLD: the response factor, the lab and the vessel cancel; "
                      "what is left is the cook's temperature-time response. The LEVELS are not "
                      "read (they stay validation)."),
            ))
    return tuple(rows)


B10_FIT_ROWS: Tuple[Dict[str, Any], ...] = _fold_rows()
assert len(B10_FIT_ROWS) == 6


def install_b10_rows(with_yiltirak: bool) -> None:
    """Splice (or leave out) the four systems and six rows; idempotent."""
    present = {r["id"] for r in B23.ACTIVE_FIT_ROWS}
    if with_yiltirak:
        for name, spec in B10_SYSTEMS.items():
            B23.SYSTEMS.setdefault(name, spec)
        new = tuple(r for r in B10_FIT_ROWS if r["id"] not in present)
        B23.ACTIVE_FIT_ROWS = tuple(B23.ACTIVE_FIT_ROWS) + new
        B23.FIT_ROWS = tuple(B23.FIT_ROWS) + new
    else:
        ids = {r["id"] for r in B10_FIT_ROWS}
        B23.ACTIVE_FIT_ROWS = tuple(r for r in B23.ACTIVE_FIT_ROWS if r["id"] not in ids)
        B23.FIT_ROWS = tuple(r for r in B23.FIT_ROWS if r["id"] not in ids)


# B8's own functions, captured BEFORE `configure` rebinds the module's names to
# B10's, so B10's versions can build on them without calling themselves.
_B8_FULL_BOUNDS = B8.full_bounds
_B8_INCUMBENT_VECTOR = B8.incumbent_vector
_B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS = B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS
_B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT = B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT

# ===========================================================================
# 2. THE VECTOR: 48 + 1
# ===========================================================================
#: B9's 48 coordinates in their order (the lumped slot now carries Ea_sugar_trunk)
#: plus Ea_thiol_assembly appended, so every B2.3 index below N_K + 5 keeps meaning.
ALL_KEYS: Tuple[str, ...] = tuple(B8.ALL_KEYS) + ("Ea_thiol_assembly",)
SUGAR_SLOT = B23.N_K            # the coordinate B2.3 calls "Ea_lumped_formation"
THIOL_SLOT = len(ALL_KEYS) - 1

FREE_KEYS: Tuple[str, ...] = tuple(B8.FREE_KEYS) + ("Ea_lumped_formation", "Ea_thiol_assembly")
FREE_CLAUSE_OF: Dict[str, str] = {
    **dict(B8.FREE_CLAUSE_OF),
    "Ea_lumped_formation": "R6 route barrier, sugar trunk (B10; the slot keeps its B2.3 name)",
    "Ea_thiol_assembly": "R6 route barrier, thiol assembly (B10)",
}
FREE_INDEX: Tuple[int, ...] = tuple(ALL_KEYS.index(k) for k in FREE_KEYS)
#: For the Laplace / profile derivatives: the B2.3-named slot is the sugar-trunk route here.
COORDINATE_OF_OVERRIDES: Dict[str, Dict[str, str]] = {
    "Ea_lumped_formation": {"block": "formation_Ea_by_route_kJ_mol", "key": "sugar_trunk", "kind": "Ea"},
    "Ea_thiol_assembly": {"block": "formation_Ea_by_route_kJ_mol", "key": "thiol_assembly", "kind": "Ea"},
}
FROZEN_KEYS: Tuple[str, ...] = tuple(k for k in ALL_KEYS if k not in set(FREE_KEYS))
assert len(FREE_KEYS) == 25 and len(ALL_KEYS) == 49


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B8_FULL_BOUNDS()
    lower, upper = np.array(lower, dtype=float), np.array(upper, dtype=float)
    lower[SUGAR_SLOT], upper[SUGAR_SLOT] = FORMATION_EA_BOUNDS_BY_ROUTE["sugar_trunk"]
    lo_t, hi_t = FORMATION_EA_BOUNDS_BY_ROUTE["thiol_assembly"]
    return np.append(lower, lo_t), np.append(upper, hi_t)


def route_ea_from_vector(x: np.ndarray) -> Dict[str, float]:
    return {"sugar_trunk": float(x[SUGAR_SLOT]), "thiol_assembly": float(x[THIOL_SLOT])}


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    """B2.3's builder with the two route barriers in place of the lumped one."""
    fitted, _lumped, decay_ea, _drift = B23.unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B23.B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, route_ea_from_vector(x), decay_ea))
    return parameters


def incumbent_vector() -> np.ndarray:
    """B9's frozen optimum with both route barriers at their band centres, clipped."""
    from generate_kinetic_core_b8_laplace import frozen_vector

    report = json.loads(B9_FIT_REPORT.read_text(encoding="utf-8"))
    x9 = frozen_vector(report)
    assert len(x9) == B23.N_K + B23.N_EXTRA, len(x9)
    x = np.append(x9, FORMATION_EA_PRIOR_CENTRE["thiol_assembly"])
    x[SUGAR_SLOT] = FORMATION_EA_PRIOR_CENTRE["sugar_trunk"]
    lower, upper = full_bounds()
    return np.clip(x, lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


_B23_BUILD_PARAMETERS_ORIGINAL = B23.build_parameters


def restore() -> None:
    """Undo `configure`: B2.3's own builder back, the six rows and four systems out."""
    B23.build_parameters = _B23_BUILD_PARAMETERS_ORIGINAL
    install_b10_rows(False)
    for name in B10_SYSTEMS:
        B23.SYSTEMS.pop(name, None)
    B8.full_bounds, B8.incumbent_vector = _B8_FULL_BOUNDS, _B8_INCUMBENT_VECTOR
    B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS = _B8_ALL_KEYS, _B8_FREE_KEYS, _B8_FREE_INDEX, _B8_FROZEN_KEYS
    B8.FREE_CLAUSE_OF, B8.MEMBER_DIR, B8.OUT_FIT_REPORT = _B8_FREE_CLAUSE_OF, _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT


def configure(with_yiltirak: bool) -> None:
    """Bind the B8 runner to B10's vector, bounds, rows and artifact paths."""
    global MEMBER_DIR, OUT_FIT_REPORT, WITH_YILTIRAK
    WITH_YILTIRAK = with_yiltirak
    tag = "b10" if with_yiltirak else "b10_noyil"
    MEMBER_DIR = data_paths.VALIDATION_DIR / f"kinetic_core_{tag}_members"
    OUT_FIT_REPORT = data_paths.VALIDATION_DIR / f"kinetic_core_{tag}_fit_report.json"
    install_b10_rows(with_yiltirak)
    B23.build_parameters = build_parameters
    B8.ALL_KEYS = ALL_KEYS
    B8.FREE_KEYS = FREE_KEYS
    B8.FREE_INDEX = FREE_INDEX
    B8.FROZEN_KEYS = FROZEN_KEYS
    B8.FREE_CLAUSE_OF = FREE_CLAUSE_OF
    B8.full_bounds = full_bounds
    B8.incumbent_vector = incumbent_vector
    B8.MEMBER_DIR = MEMBER_DIR
    B8.OUT_FIT_REPORT = OUT_FIT_REPORT


configure(True)   # the default binding; `--without-yiltirak` rebinds before running


def _active_bounds(x: np.ndarray):
    lower, upper = full_bounds()
    out = []
    for i, key in enumerate(ALL_KEYS):
        width = float(upper[i] - lower[i])
        if (x[i] - lower[i]) <= 1e-3 * width:
            out.append({"key": key, "bound": "lower", "value": float(x[i])})
        elif (upper[i] - x[i]) <= 1e-3 * width:
            out.append({"key": key, "bound": "upper", "value": float(x[i])})
    return out


def consolidate() -> Dict[str, Any]:
    payload = B8.consolidate()          # writes OUT_FIT_REPORT in B2.3's schema
    best = min(B8.load_members(), key=lambda m: m["cost"])
    x = np.array(best["x_full"], dtype=float)
    assert len(x) == len(ALL_KEYS), len(x)
    routes = route_ea_from_vector(x)
    fr = payload["frozen_parameters"]
    fr["formation_Ea_by_route_kJ_mol"] = routes
    fr["lumped_formation_Ea_kJ_mol"] = routes["sugar_trunk"]   # the slot's value, for pre-B10 readers
    fr["formation_route_of"] = dict(FORMATION_ROUTE_OF)
    payload["wave"] = (f"{WAVE} -- the temperature structure: one formation barrier -> two by route"
                       + ("" if WITH_YILTIRAK else " (LEAVE-YILTIRAK-OUT variant, prereg T5)"))
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b10_fit.py"
    payload["prereg"] = data_paths.rel(PREREG)
    payload["declaration"] = "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 20"
    payload["objective"]["form"] = (
        "B9's 54 rows" + (" + 6 within-study Yiltirak 2026 folds (60 rows)" if WITH_YILTIRAK else " (the six Yiltirak folds NOT installed)")
        + "; two FREE route barriers replace the frozen lumped one"
    )
    payload["objective"]["n_new_rows"] = len(B10_FIT_ROWS) if WITH_YILTIRAK else 0
    payload["objective"]["new_row_ids"] = [r["id"] for r in B10_FIT_ROWS] if WITH_YILTIRAK else []
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    payload["objective"]["not_comparable_note"] = (
        "TOTAL COST IS NOT COMPARABLE TO B9's when the six folds are installed. "
        "`sum_r2_level_shared_with_b2_4` is the like-for-like comparator over the rows both waves scored."
    )
    payload["free_set"] = {
        "n_free": len(FREE_KEYS), "n_frozen": len(FROZEN_KEYS),
        "keys": list(FREE_KEYS), "clause": FREE_CLAUSE_OF, "frozen_keys": list(FROZEN_KEYS),
    }
    payload["objective"]["n_free_parameters"] = len(FREE_KEYS)
    payload["objective"]["n_frozen"] = len(FROZEN_KEYS)
    payload["t_structure"]["formation_routes"] = {
        route: {
            "fitted_kj_mol": routes[route],
            "band_kj_mol": list(FORMATION_EA_BOUNDS_BY_ROUTE[route]),
            "prior_centre_kj_mol": FORMATION_EA_PRIOR_CENTRE[route],
            "prior_source": FORMATION_EA_PRIOR_SOURCE[route],
            "n_steps": sum(1 for k, r in FORMATION_ROUTE_OF.items() if r == route),
        }
        for route in FORMATION_ROUTES
    }
    payload["start_vector"] = ("B9's frozen optimum with both route barriers at their prior centres "
                               "(start 0); B8's perturbation protocol around it (start 1)")
    payload["active_bounds"] = _active_bounds(x)
    OUT_FIT_REPORT.write_text(json.dumps(payload, indent=2, default=str))
    print(f"rewrote {OUT_FIT_REPORT} as {WAVE}: routes {routes}; active bounds: {[a['key'] for a in payload['active_bounds']]}")
    return payload


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    parser.add_argument("--without-yiltirak", dest="with_yiltirak", action="store_false", default=True)
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B10 is pre-registered; write the prereg before running the fit"
    configure(args.with_yiltirak)
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    expected_rows = 60 if args.with_yiltirak else 54
    assert len(B23.ACTIVE_FIT_ROWS) == expected_rows, len(B23.ACTIVE_FIT_ROWS)

    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE + ("" if args.with_yiltirak else "_noyil")
    member["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    member["route_barriers"] = route_ea_from_vector(np.array(member["x_full"], dtype=float))
    dest = MEMBER_DIR / f"{'b10' if args.with_yiltirak else 'b10_noyil'}_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
