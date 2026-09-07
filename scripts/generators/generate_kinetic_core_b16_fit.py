#!/usr/bin/env python
"""
Build Wave B16 -- THE THIOL SINKS AGAINST A 100 C TIME SERIES AND MEASURED TTCA DECAY (2026-09-07).

Pre-registered in ``results/validation/kinetic_core_b16_prereg.md``. Base wave B9 (54 rows, 23 free).
Ten rows added: seven within-study ratios from Schieberle, Hofmann & Muench 2000 Table IV (the
Hofmann 1998 pentose pot held at 100 C for 30 / 60 / 360 / 720 min, stable isotope dilution) and
three TTCA-remaining rows from Zhai et al. 2021 (TTCA 10 mM, pH 7, 100 / 120 / 140 C, 60 min,
zero-order fits). No new coordinate. Two variants: ``b16`` keeps every B9 band (the thiol-sink barrier
ceiling 102 kJ/mol, the owner's decision); ``b16_lift`` raises that ceiling to 160 kJ/mol as
INFORMATION and cannot ship under the prereg.

Usage:
    python scripts/generators/generate_kinetic_core_b16_fit.py --start 0 [--lift-ceiling]
    python scripts/generators/generate_kinetic_core_b16_fit.py --consolidate [--lift-ceiling]
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
from src.kinetic_core.parameters_sulfur import OX_AMBIENT_MMOL_L  # noqa: E402

WAVE = "B16"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b16_prereg.md"
B9_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json"
LIFTED_THIOL_SINK_CEILING_KJ_MOL = 160.0
MEMBER_DIR = data_paths.VALIDATION_DIR / "kinetic_core_b16_members"
OUT_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b16_fit_report.json"
LIFT = False

# ===========================================================================
# 1. THE NEW SYSTEMS AND ROWS
# ===========================================================================
_HOF = B23.SYSTEMS["hofmann_pentose_pH5"]
SCHIEBERLE_MINUTES: Tuple[float, ...] = (30.0, 60.0, 360.0, 720.0)
B16_SYSTEMS: Dict[str, Dict[str, Any]] = {
    f"schieberle_100C_{int(m)}": dict(
        initial=dict(_HOF["initial"]), t_c=100.0, minutes=m, ph=5.0, buffer=_HOF["buffer"], grid_points=25,
        anchor=(f"Schieberle, Hofmann & Muench 2000 Table IV footnote: ribose 10 mmol + cysteine 3.3 mmol in "
                f"phosphate buffer (100 mL; 0.5 mol/L; pH 5.0) at 100 C in a laboratory autoclave, {m:.0f} min; "
                "schieberle2000_extraction.md sec. 2"),
    )
    for m in SCHIEBERLE_MINUTES
}
ZHAI_TTCA_K: Dict[float, Tuple[float, float]] = {100.0: (10.331, 0.0271), 120.0: (9.9718, 0.0651), 140.0: (9.3375, 0.0813)}
for t_c in ZHAI_TTCA_K:
    B16_SYSTEMS[f"zhai2021_ttca_{int(t_c)}"] = dict(
        initial={"TTCA": 10.0, "OX": OX_AMBIENT_MMOL_L}, t_c=t_c, minutes=60.0, ph=7.0, buffer=B23.BUFFER_NONE,
        grid_points=25,
        anchor=(f"Zhai et al. 2021 JAFC 69:10648, Methods + Fig. 3a: aqueous TTCA 10 mmol/L, pH 7 (NaOH), oil bath "
                f"{t_c:.0f} C; zero-order fit c = c0 - k t; zhai2021_extraction.md sec. 2"),
    )

#: Table IV, ug per 100 mL pot.
TABLE_IV: Dict[str, Dict[float, float]] = {
    "MFT": {30.0: 4.5, 60.0: 13.8, 360.0: 156.0, 720.0: 179.0},
    "FFT": {30.0: 2.0, 60.0: 3.1, 360.0: 110.0, 720.0: 132.0},
}
SIGMA_LOG_FOLD = 0.10      # the source prints no SD; B10's fold sigma
SIGMA_LOG_TTCA = 0.05      # the zero-order fits' R2 0.95-0.99 on six points
FACTOR_13_SIGMA = 0.15     # "by a factor of 13", a rounded statement


def _rows() -> Tuple[Dict[str, Any], ...]:
    rows = []
    for species in ("MFT", "FFT"):
        for m in SCHIEBERLE_MINUTES[1:]:
            rows.append(dict(
                id=f"schieberle_{species}_fold_{int(m)}_over_30", system=f"schieberle_100C_{int(m)}",
                system_b="schieberle_100C_30", kind="cross_system_ratio", species=species,
                target=TABLE_IV[species][m] / TABLE_IV[species][30.0], sigma_log=SIGMA_LOG_FOLD,
                anchor=(f"Schieberle 2000 Table IV: {species} {TABLE_IV[species][30.0]} -> {TABLE_IV[species][m]} ug "
                        f"over 30 -> {m:.0f} min at 100 C (SIDA)"),
                note="WITHIN-STUDY FOLD in the fit's own Hofmann 1998 pot: the sink/formation balance at 100 C.",
            ))
    rows.append(dict(
        id="schieberle_MFT_145C20min_over_100C360min", system="hofmann_pentose_pH5", system_b="schieberle_100C_360",
        kind="cross_system_ratio", species="MFT", target=1.0 / 13.0, sigma_log=FACTOR_13_SIGMA,
        anchor=("Schieberle 2000 p. 141: 'by a factor of 13 higher yields of e.g. 2-methyl-3-furanthiol were "
                "obtained [100 C / 6 h] than in the mixture reacted for 20 min at 145 C (cf. Tables III and IV)'"),
        note="Links the 100 C series to the objective's 145 C level: the same pot, two cooks.",
    ))
    for t_c, (c0, k) in ZHAI_TTCA_K.items():
        rows.append(dict(
            id=f"zhai2021_ttca_remaining_{int(t_c)}C_60min", system=f"zhai2021_ttca_{int(t_c)}", kind="conc",
            species="TTCA", target=c0 - k * 60.0, sigma_log=SIGMA_LOG_TTCA,
            anchor=(f"Zhai 2021 Fig. 3a fit at {t_c:.0f} C: y = -{k} x + {c0} (mmol/L vs min); "
                    f"{c0 - k * 60.0:.2f} mmol/L remaining at 60 min"),
            note="A MEASURED RATE on a species the lane carries (the first temperature series on k_ttca_deg).",
        ))
    return tuple(rows)


B16_FIT_ROWS: Tuple[Dict[str, Any], ...] = _rows()
assert len(B16_FIT_ROWS) == 10


def install_b16_rows(on: bool) -> None:
    present = {r["id"] for r in B23.ACTIVE_FIT_ROWS}
    if on:
        for name, spec in B16_SYSTEMS.items():
            B23.SYSTEMS.setdefault(name, spec)
        new = tuple(r for r in B16_FIT_ROWS if r["id"] not in present)
        B23.ACTIVE_FIT_ROWS = tuple(B23.ACTIVE_FIT_ROWS) + new
        B23.FIT_ROWS = tuple(B23.FIT_ROWS) + new
    else:
        ids = {r["id"] for r in B16_FIT_ROWS}
        B23.ACTIVE_FIT_ROWS = tuple(r for r in B23.ACTIVE_FIT_ROWS if r["id"] not in ids)
        B23.FIT_ROWS = tuple(r for r in B23.FIT_ROWS if r["id"] not in ids)
        for name in B16_SYSTEMS:
            B23.SYSTEMS.pop(name, None)


# ===========================================================================
# 2. THE VECTOR: B9's 48, B9's 23 free; the lift variant widens one band
# ===========================================================================
_B8_FULL_BOUNDS, _B8_INCUMBENT_VECTOR = B8.full_bounds, B8.incumbent_vector
_B8_MEMBER_DIR, _B8_OUT_FIT_REPORT = B8.MEMBER_DIR, B8.OUT_FIT_REPORT
ALL_KEYS, FREE_KEYS, FREE_INDEX, FROZEN_KEYS = B8.ALL_KEYS, B8.FREE_KEYS, B8.FREE_INDEX, B8.FROZEN_KEYS
FREE_CLAUSE_OF = B8.FREE_CLAUSE_OF
THIOL_SINK_SLOT = ALL_KEYS.index("Ea_decay_thiol_sink")
assert len(FREE_KEYS) == 23


def full_bounds() -> Tuple[np.ndarray, np.ndarray]:
    lower, upper = _B8_FULL_BOUNDS()
    lower, upper = np.array(lower, dtype=float), np.array(upper, dtype=float)
    if LIFT:
        upper[THIOL_SINK_SLOT] = LIFTED_THIOL_SINK_CEILING_KJ_MOL
    return lower, upper


def incumbent_vector() -> np.ndarray:
    from generate_kinetic_core_b8_laplace import frozen_vector

    x9 = frozen_vector(json.loads(B9_FIT_REPORT.read_text(encoding="utf-8")))
    assert len(x9) == B23.N_K + B23.N_EXTRA, len(x9)
    lower, upper = full_bounds()
    return np.clip(x9, lower, upper)


def residual_vector(x_full: np.ndarray, quick: bool) -> np.ndarray:
    return B8.residual_vector(x_full, quick)


def restore() -> None:
    install_b16_rows(False)
    B8.full_bounds, B8.incumbent_vector = _B8_FULL_BOUNDS, _B8_INCUMBENT_VECTOR
    B8.MEMBER_DIR, B8.OUT_FIT_REPORT = _B8_MEMBER_DIR, _B8_OUT_FIT_REPORT


def configure(lift: bool = False) -> None:
    global LIFT, MEMBER_DIR, OUT_FIT_REPORT
    LIFT = bool(lift)
    tag = "b16_lift" if LIFT else "b16"
    MEMBER_DIR = data_paths.VALIDATION_DIR / f"kinetic_core_{tag}_members"
    OUT_FIT_REPORT = data_paths.VALIDATION_DIR / f"kinetic_core_{tag}_fit_report.json"
    install_b16_rows(True)
    B8.full_bounds, B8.incumbent_vector = full_bounds, incumbent_vector
    B8.MEMBER_DIR, B8.OUT_FIT_REPORT = MEMBER_DIR, OUT_FIT_REPORT


configure(False)


def consolidate() -> Dict[str, Any]:
    payload = B8.consolidate()
    best = min(B8.load_members(), key=lambda m: m["cost"])
    x = np.array(best["x_full"], dtype=float)
    tag = "b16_lift" if LIFT else "b16"
    payload["wave"] = (f"{WAVE}{' (lift variant: thiol-sink ceiling 160, information only)' if LIFT else ''} -- "
                       "the thiol sinks against Schieberle 2000's 100 C series and Zhai 2021's TTCA decay")
    payload["generated_by"] = "scripts/generators/generate_kinetic_core_b16_fit.py" + (" --lift-ceiling" if LIFT else "")
    payload["prereg"] = data_paths.rel(PREREG)
    payload["declaration"] = "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendment 26"
    payload["objective"]["form"] = "B9's 54 rows + 7 Schieberle 2000 within-study ratios + 3 Zhai 2021 TTCA-remaining rows (64 rows, 23 free)"
    payload["objective"]["n_new_rows"] = len(B16_FIT_ROWS)
    payload["objective"]["new_row_ids"] = [r["id"] for r in B16_FIT_ROWS]
    payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS)
    payload["thiol_sink_ceiling_kj_mol"] = float(full_bounds()[1][THIOL_SINK_SLOT])
    payload["variant"] = tag
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
    print(f"rewrote {OUT_FIT_REPORT} as {tag}: thiol-sink Ea {x[THIOL_SINK_SLOT]:.1f} (ceiling {upper[THIOL_SINK_SLOT]:.0f}); "
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
    parser.add_argument("--lift-ceiling", dest="lift", action="store_true")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B16 is pre-registered; write the prereg before running the fit"
    configure(args.lift)
    MEMBER_DIR.mkdir(parents=True, exist_ok=True)
    assert len(B23.ACTIVE_FIT_ROWS) == 64, len(B23.ACTIVE_FIT_ROWS)
    if args.consolidate:
        consolidate()
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE + ("_lift" if LIFT else "")
    dest = MEMBER_DIR / f"{'b16_lift' if LIFT else 'b16'}_s{args.start}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
