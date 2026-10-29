"""
WAVE B64 verification V5 and candidate C (2026-10-28; kinetic_core_b64_prereg.md sec. 3-4).

V5: under each reading's profiles, B60's quick setting (the optimiser's) agrees with the careful setting
(81 pH nodes, rtol 1e-8) within 0.03 dex on all 52 observables at B60's optimum, as B60's verification
required of the isothermal objective. Candidate C: B60's vector's cost on the 52 rows under each reading,
beside the isothermal cost. Writes results/validation/kinetic_core_b64_verification.json.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b64_fit as B64  # noqa: E402
import generate_kinetic_core_b60_verification as V60  # noqa: E402
from src import data_paths  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b64_verification.json"
V5_TOL_DEX = 0.03


def _cost(x: np.ndarray) -> float:
    r = B64.B60.residual_vector(x, True)
    return float(0.5 * np.sum(r ** 2))


def main() -> int:
    x = B64.incumbent_vector()
    quick = B64.QUICK_SETTING
    report = {"artifact": "kinetic_core_b64_verification", "prereg": data_paths.rel(B64.PREREG),
              "quick_setting": quick, "cost_isothermal_b60_vector": _cost(x), "readings": {}}
    systems_before = {k: dict(v) for k, v in B64.B23.SYSTEMS.items()}
    run_systems_before = B64.B23.run_systems
    for reading in sorted(B64.READINGS):
        B64.B23.SYSTEMS.clear()
        B64.B23.SYSTEMS.update({k: dict(v) for k, v in systems_before.items()})
        B64.B23.run_systems = run_systems_before
        B64.install_reading(reading)
        t0 = time.perf_counter()
        fast = V60._observables(x, True, int(quick["ph_nodes"]), "log_and_linear", rtol=float(quick["rtol"]))
        wall_quick = time.perf_counter() - t0
        careful = V60._observables(x, False, 81, "log_and_linear")
        dex = {k: V60._dex(fast[k], careful[k], "ph_endpoint" if "ph" in k.lower() and "endpoint" in k.lower()
                           else "level") for k in careful}
        worst = max(dex, key=dex.get)
        report["readings"][reading] = {
            "v5_max_dex_quick_vs_careful": dex[worst], "v5_worst_row": worst,
            "v5_rows_over_tol": sorted(k for k, v in dex.items() if v > V5_TOL_DEX),
            "v5_pass": dex[worst] <= V5_TOL_DEX, "wall_seconds_quick_objective": wall_quick,
            "cost_b60_vector": _cost(x),
        }
        print(reading, json.dumps(report["readings"][reading]), flush=True)
    OUT.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
