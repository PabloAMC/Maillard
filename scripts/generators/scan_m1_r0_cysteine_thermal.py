#!/usr/bin/env python
"""The engine scan behind M1.1's exact-density regression: R0's decisive levels along k_cys_thermal (2026-11-05).

WHY
===
M1's results review (results/validation/kinetic_core_m1_results_review.md) found R0's RES_Cys predictive under the
alternative floor to be about -746 by an exact 1-D integral along ``k_cys_thermal`` (R0's one identified coordinate
that moves its cysteine; Laplace sd 0.279), against the linearised -175. M1's records keep R0's levels and Jacobian at
the carried optimum only, so the integral needs the engine at shifted points. This script runs that scan once
(25 points, delta = -1.2 ... +1.2 dex in steps of 0.1, the reviewer's grid) at R0's carried M1 optimum, with M1's
scoring numerics (inline, rtol 1e-8), isothermal reading, and writes the levels as a small test fixture,
``tests/fixtures/m1_r0_cysteine_thermal_scan.json``. ``tests/unit/test_predictive_densities.py`` integrates it with
``src/kinetic_core/predictive_densities.py``. A reproduction of the review's check, not a decisive re-run (approved by
the coordinator, 2026-11-05): it changes no M1 record and decides nothing.

    ./scripts/docker_maillard.sh run "MAILLARD_WORKERS=4 nice -n 10 python scripts/generators/scan_m1_r0_cysteine_thermal.py"
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "fixtures" / "m1_r0_cysteine_thermal_scan.json"
GENERATED_BY = "scripts/generators/scan_m1_r0_cysteine_thermal.py"
COORDINATE = "k_cys_thermal"
DELTAS = tuple(round(-1.2 + 0.1 * i, 1) for i in range(25))
READING = "isothermal"


def main() -> int:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    import fit_m1_supply_combinations as G  # the M1 generator: its forms, its numerics, its decisive levels
    from src import provenance

    forms, _ = G.forms_for_score(G._read(G.CHECKS_OUT))
    fm = forms["R0"]
    i = list(fm.obj.keys).index(COORDINATE)
    points = []
    for d in DELTAS:
        th = fm.theta.copy()
        th[i] += d
        points.append(fm.obj.point(th))
    G.FAILURES.clear()
    levels = G.decisive_many(points, READING, fm.numerics)
    failures = list(G.FAILURES)
    dec = json.loads((G.V / "kinetic_core_m1_decisive_comparison.json").read_text(encoding="utf-8"))
    recorded = np.asarray(dec["levels"]["R0"][READING]["levels"], dtype=float)
    at_zero = np.asarray(levels[DELTAS.index(0.0)], dtype=float)
    carried = G.carried_path("R0")
    payload = {
        "artifact": "m1_r0_cysteine_thermal_scan",
        "what": "R0's log10 decisive levels (mmol/L reservoir, ug/L D1; unfloored) at k_cys_thermal + delta dex, the "
                "other coordinates at R0's carried M1 optimum",
        "rival": "R0", "coordinate": COORDINATE, "carried_theta": float(fm.theta[i]),
        "laplace_sd": float(fm.lap["sigma"][COORDINATE]),
        "carried_from_start": int(G._read(carried)["carried_from_start"]),
        "reading": READING, "numerics": G.dataclasses.asdict(fm.numerics),
        "point_labels": G.point_labels(), "deltas": list(DELTAS),
        "levels": [np.asarray(lv, dtype=float).tolist() for lv in levels],
        "failed_pot_runs": len(failures),
        "check_delta_zero_vs_decisive_record_max_abs_dex": float(np.nanmax(np.abs(at_zero - recorded))),
        "provenance": provenance.provenance_block(
            "m1_r0_cysteine_thermal_scan", generated_by=GENERATED_BY,
            inputs=[carried, G.CHECKS_OUT, G.V / "kinetic_core_m1_decisive_comparison.json",
                    ROOT / G.GENERATED_BY, ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/parameters_sulfur.py", ROOT / "src/kinetic_core/config.py",
                    ROOT / "scripts/generators/build_submodel_probe.py", Path(__file__)]),
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(DELTAS)} points, {len(failures)} failed pot runs, delta 0 vs the "
          f"decisive record max |diff| {payload['check_delta_zero_vs_decisive_record_max_abs_dex']:.2e} dex")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
