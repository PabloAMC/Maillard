"""
WAVE B64 -- B60's objective refitted on the thermal history each source records (2026-10-28).

Every sulfur fit system was integrated isothermal from t = 0. Two sources say otherwise:
- Hofmann & Schieberle 1998 (14 systems): the same laboratory, pot and 200 mL Roth autoclave are run
  "by raising the temperature within 20 min from 20 to 145 °C" in Hofmann 1995 and 1997 (verbatim; the
  1997 paper also calls that pot "thermally treated for 20 min at 145 °C", the 1998 wording). Reading
  ``ramp_only`` (the point, p ~ 0.8) integrates 20 -> 145 C over 20 min; ``ramp_plus_hold`` (the declared
  alternative, p ~ 0.12) adds a 20 min hold at 145 C.
- van Seeventer 2001 (1 system): "from room temperature to 130 °C in 10 min and was kept at 130 °C for
  20 min", stated, so both readings carry it.
Every other system's source states no come-up and stays isothermal (declared in the prereg). Nothing else
changes: B60's rows, bands, free set, weighting and numerics. Pre-registration:
results/validation/kinetic_core_b64_prereg.md.

Usage (starts in parallel, niced, then consolidate):
    python scripts/generators/generate_kinetic_core_b64_fit.py --reading ramp_only --start 0
    python scripts/generators/generate_kinetic_core_b64_fit.py --reading ramp_only --start 1
    python scripts/generators/generate_kinetic_core_b64_fit.py --reading ramp_only --consolidate
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b60_fit as B60  # noqa: E402  (B51's rows and physics, B60's numerics and quick setting)
from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import parallel, ph_state, sulfur  # noqa: E402

B51, B9, B23, B8 = B60.B51, B60.B9, B60.B23, B60.B8
CELSIUS = 273.15

#: B64 runs on the inline pH coupling (C5): the limit B60's node fixed point converges to, within
#: 0.002 dex on the probe pots, without the node grid or the sweeps (kinetic_core_b64_prereg.md sec. 2).
sulfur.PH_COUPLING = "inline"

WAVE = "B64"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b64_prereg.md"
B60_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b60_fit_report.json"

#: The Hofmann 1998 systems: every aqueous 145 C / 20 min pot of that paper in the fit.
HOFMANN_1998_SYSTEMS = (
    "hofmann_pentose_pH5", "hofmann_glucose_pH5", "hofmann_fructose_pH5",
    "fed_ribose_h2s", "fed_tdp_h2s", "fed_furfural_h2s", "fed_nf_h2s", "fed_nf_cys",
    "fed_c2c3", "fed_c2c3_pH3", "fed_c2c3_pH7", "fed_thiamine", "fed_mgo_h2s_1to1", "fed_mgo_h2s_1to2",
)
#: The recorded profiles, ((t_min, T_C), ...) from t = 0, per reading.
READINGS: Dict[str, Dict[str, tuple]] = {
    "ramp_only": {
        **{name: ((0.0, 20.0), (20.0, 145.0)) for name in HOFMANN_1998_SYSTEMS},
        "vanseeventer_130C": ((0.0, 20.0), (10.0, 130.0), (30.0, 130.0)),
    },
    "ramp_plus_hold": {
        **{name: ((0.0, 20.0), (20.0, 145.0), (40.0, 145.0)) for name in HOFMANN_1998_SYSTEMS},
        "vanseeventer_130C": ((0.0, 20.0), (10.0, 130.0), (30.0, 130.0)),
    },
}



def _newtonian(tau: float, n: int = 60) -> tuple:
    """Amendment 2: a bath-heated pot's Newtonian approach to 145 C over the stated 20 min."""
    import math

    return tuple((20.0 * i / n, 145.0 - 125.0 * math.exp(-(20.0 * i / n) / tau)) for i in range(n + 1))


#: Amendment 2 (2026-10-28, before any such fit): the ramp's SHAPE moved the ramp-only cost 45.4 -> 33.4 with
#: no refit (generate_kinetic_core_b64_shape_sensitivity.py), so the Newtonian shapes are refitted as readings
#: of their own. A swept range (within 6 K and 1 K of 145 C at 20 min), reported side by side; none is chosen.
for _tau, _name in ((20.0 / 3.0, "ramp_only_newtonian_tau6.7"), (20.0 / 5.0, "ramp_only_newtonian_tau4")):
    READINGS[_name] = {**{name: _newtonian(_tau) for name in HOFMANN_1998_SYSTEMS},
                       "vanseeventer_130C": READINGS["ramp_only"]["vanseeventer_130C"]}

#: The shipped reading while B64 chooses none (its verdict): no profile on any system, i.e. B60's objective exactly.
READINGS["isothermal"] = {}

QUICK_SETTING = B60.QUICK_SETTING
REMOVED_ROW_IDS = B60.REMOVED_ROW_IDS
FREE_KEYS = B60.FREE_KEYS
ALL_KEYS = B60.ALL_KEYS
full_bounds = B60.full_bounds


def install_reading(reading: str) -> Dict[str, tuple]:
    """Put each recorded profile on its system and route ``run_systems`` through it."""
    profiles = READINGS[reading]
    for name, points in profiles.items():
        spec = B23.SYSTEMS[name]
        assert max(c for _, c in points) <= float(spec["t_c"]) + 1e-9, name  # a Newtonian approach never reaches it
        B23.SYSTEMS[name] = {**spec, "minutes": float(points[-1][0]), "profile_c": tuple(points)}
    B23.run_systems = _run_systems_with_profiles
    return profiles


def _run_systems_with_profiles(parameters: Mapping[str, Any], quick: bool, drift=None) -> Dict[str, Any]:
    """``B23.run_systems`` with each system's recorded profile passed to the integrator; the
    systems run through ``parallel.map_ordered`` (serial unless ``MAILLARD_WORKERS`` > 1)."""

    def integrate_one(name: str):
        spec = B23.SYSTEMS[name]
        minutes = float(spec["minutes"])
        points = int(spec.get("grid_points", 0))
        grid = np.linspace(0.0, minutes, points) if points >= 3 else np.array([0.0, minutes])
        profile = spec.get("profile_c")
        return B23.integrate_sulfur(
            parameters, float(spec["t_c"]) + CELSIUS, spec["initial"], grid,
            ph=float(spec["ph"]), buffer_spec=spec.get("buffer"), ph_drift=drift,
            ph_nodes=B23.PH_NODES_QUICK if quick else B23.PH_NODES_CAREFUL,
            ph_iterations=B23.PH_ITERS_QUICK if quick else B23.PH_ITERS_CAREFUL,
            rtol=1e-5 if quick else 1e-8, atol=1e-14,
            temperature_profile=(
                None if profile is None
                else (np.array([t for t, _ in profile]), np.array([c for _, c in profile]) + CELSIUS)
            ),
        )

    names = list(B23.SYSTEMS)
    runs = parallel.map_ordered(integrate_one, names)
    return {name: {"run": run, "spec": B23.SYSTEMS[name]} for name, run in zip(names, runs)}


def incumbent_vector() -> np.ndarray:
    """B60's frozen optimum, clipped into the same bands."""
    lower, upper = full_bounds()
    return np.clip(frozen_vector(json.loads(B60_REPORT.read_text(encoding="utf-8"))), lower, upper)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reading", choices=sorted(READINGS), required=True)
    parser.add_argument("--start", type=int)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    parser.add_argument("--budget", type=int, default=B8.EVAL_BUDGET)
    parser.add_argument("--consolidate", action="store_true")
    parser.add_argument("--start-from", dest="start_from", type=Path,
                        help="a member JSON whose x_full is start 0's vector instead of B60's optimum "
                             "(the cross-restart check against local minima, owner 2026-10-28)")
    parser.add_argument("--cost-only", action="store_true",
                        help="print B60's vector's cost under the reading (candidate C) and exit")
    args = parser.parse_args(argv)
    assert PREREG.exists(), "B64 is pre-registered; write the prereg before running the fit"
    profiles = install_reading(args.reading)
    member_dir = data_paths.VALIDATION_DIR / f"kinetic_core_b64_{args.reading}_members"
    out_report = data_paths.VALIDATION_DIR / f"kinetic_core_b64_{args.reading}_fit_report.json"
    record = {"reading": args.reading, "profiles_c": {k: [list(p) for p in v] for k, v in profiles.items()},
              "profile_node_step_k": sulfur.PROFILE_NODE_STEP_K, "ph_coupling": sulfur.PH_COUPLING}

    if args.cost_only:
        r = B60.residual_vector(incumbent_vector(), args.quick)
        print(json.dumps({"reading": args.reading, "cost_b60_vector": float(0.5 * np.sum(r ** 2)),
                          "n_rows": int(r.size)}))
        return 0

    member_dir.mkdir(parents=True, exist_ok=True)
    B8.MEMBER_DIR = member_dir
    B8.OUT_FIT_REPORT = out_report
    B8.incumbent_vector = incumbent_vector
    if args.start_from is not None:
        lower, upper = full_bounds()
        origin = np.clip(np.array(json.loads(args.start_from.read_text())["x_full"], dtype=float), lower, upper)
        B8.incumbent_vector = lambda: origin
        record["start_from"] = data_paths.rel(args.start_from)

    if args.consolidate:
        payload = B8.consolidate()
        payload["wave"] = f"B64 -- B60's objective on the recorded thermal history, reading {args.reading} (2026-10-28)"
        payload["generated_by"] = "scripts/generators/generate_kinetic_core_b64_fit.py"
        payload["prereg"] = data_paths.rel(PREREG)
        payload["thermal_history_b64"] = record
        payload["objective"]["form"] = (
            "B60's objective (B47's 52 corrected rows, B51's pH state, B60's numerics) with each fit "
            "system integrated along the thermal history its source records")
        payload["objective"]["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
        payload["ph_state_corrections_b51"] = ph_state.B51.as_dict()
        payload["ph_numerics_b60"] = {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT,
                                      "relaxation": sulfur.PH_FIXED_POINT_RELAXATION,
                                      "quick_setting": QUICK_SETTING}
        payload["start_vector"] = "B60's frozen optimum (start 0); B8's perturbation protocol (start 1)"
        lower, upper = full_bounds()
        x = frozen_vector(payload)
        payload["active_bounds"] = [
            {"key": key, "bound": side, "value": float(x[i])}
            for i, key in enumerate(ALL_KEYS)
            for side, hit in (("lower", (x[i] - lower[i]) <= 1e-3 * (upper[i] - lower[i])),
                              ("upper", (upper[i] - x[i]) <= 1e-3 * (upper[i] - lower[i])))
            if hit
        ]
        out_report.write_text(json.dumps(payload, indent=2, default=str))
        print(f"wrote {out_report}; active bounds: {[a['key'] for a in payload['active_bounds']]}")
        return 0
    if args.start is None:
        parser.error("--start is required unless --consolidate or --cost-only")
    member = B8.fit_member(args.start, args.max_nfev, args.quick, args.budget)
    member["wave"] = WAVE
    member["thermal_history_b64"] = record
    member["ph_state_corrections_b51"] = ph_state.B51.as_dict()
    member["removed_row_ids"] = list(B9.VALIDATION_ROW_IDS) + list(REMOVED_ROW_IDS)
    suffix = f"_from_{args.start_from.stem}" if args.start_from is not None else ""
    dest = member_dir / f"b64_{args.reading}_s{args.start}{suffix}.json"
    dest.write_text(json.dumps(member, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
