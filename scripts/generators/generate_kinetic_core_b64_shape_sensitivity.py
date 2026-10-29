"""
WAVE B64, amendment 1 item 2: how much the ramp's SHAPE matters (choices register CR-04; 2026-10-28).

Not a refit. At each reading's best member the 52 fit rows are rescored with the Hofmann systems' heat-up
replaced by (i) a Newtonian approach to 145 C, T(t) = 145 - 125 exp(-t / tau), for tau = 20/3 and 20/5 min
(the pot within 6 K and within 1 K of 145 C at the stated 20 min), and (ii) the linear ramp followed by a
linear cool-down to 20 C over 10 or 20 min. Every value is reported; none is chosen. A swept range, not a
number: the sources give the ramp's ends and duration, not its shape or the cooling rate.
Writes results/validation/kinetic_core_b64_shape_sensitivity.json.

    python scripts/generators/generate_kinetic_core_b64_shape_sensitivity.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import generate_kinetic_core_b64_fit as B64  # noqa: E402
from generate_kinetic_core_b64_discriminator import _members  # noqa: E402
from src import data_paths  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "kinetic_core_b64_shape_sensitivity.json"


def newtonian(tau: float, until: float = 20.0, hold: float = 0.0, n: int = 60) -> tuple:
    times = np.linspace(0.0, until, n + 1)
    points = [(float(t), 145.0 - 125.0 * math.exp(-t / tau)) for t in times]
    if hold > 0.0:
        points.append((until + hold, points[-1][1]))
    return tuple(points)


def shapes(reading: str) -> dict:
    hold = 20.0 if reading == "ramp_plus_hold" else 0.0
    end = 20.0 + hold
    base = ((0.0, 20.0), (20.0, 145.0)) + (((end, 145.0),) if hold else ())
    return {
        "linear (fitted)": base,
        "newtonian tau=20/3": newtonian(20.0 / 3.0, hold=hold),
        "newtonian tau=20/5": newtonian(20.0 / 5.0, hold=hold),
        "linear + cool 10 min": base + ((end + 10.0, 20.0),),
        "linear + cool 20 min": base + ((end + 20.0, 20.0),),
    }


def main() -> int:
    members = _members()
    systems = {k: dict(v) for k, v in B64.B23.SYSTEMS.items()}
    report = {"artifact": "kinetic_core_b64_shape_sensitivity", "prereg": data_paths.rel(B64.PREREG), "by_reading": {}}
    for reading, member in members.items():
        x = np.array(member["x_full"], dtype=float)
        rows = {}
        for label, profile in shapes(reading).items():
            B64.B23.SYSTEMS.clear()
            B64.B23.SYSTEMS.update({k: dict(v) for k, v in systems.items()})
            B64.READINGS["_sweep"] = {**{n: profile for n in B64.HOFMANN_1998_SYSTEMS},
                                      "vanseeventer_130C": B64.READINGS[reading]["vanseeventer_130C"]}
            B64.install_reading("_sweep")
            r = B64.B60.residual_vector(x, True)
            rows[label] = float(0.5 * np.sum(r ** 2))
            print(reading, label, round(rows[label], 2), flush=True)
        report["by_reading"][reading] = {"member_cost": member["cost"], "cost_by_shape": rows}
    OUT.write_text(json.dumps(report, indent=2))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
