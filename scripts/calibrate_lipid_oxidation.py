"""Calibration aid for the lipid-oxidation hydroperoxide kinetics (S27 Workstream A).

This does NOT auto-fit against the external hold-out (that set is frozen). It:

  1. Exercises the kinetic core (predict_hexanal_generation) directly at the real
     in-panel anchor conditions, sweeping `max_conversion_fraction`, so a human can
     confirm that the saturation cap (a) preserves the calibrated low-temperature
     load and (b) bounds the runaway high-temperature x time extrapolation.
  2. Runs the full evaluator on the in-panel anchors as a PIN-PRESERVATION guard
     (these must stay ~1.0x — the per-process-state registry pins them) and on the
     frozen hold-out for context only.

Real measured anchors (data/benchmarks/*.json measured_volatiles, conc_ppb):
  pea_iso 40C hexanal 260 | soy_iso 40C hexanal 380 | pea UHT 140C hexanal 782.
No numbers are invented here.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.lipid_oxidation import (  # noqa: E402
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    _saturated_extent,
    predict_hexanal_generation,
)

# (label, profile, temp_C, time_min) at real in-panel anchor conditions.
ANCHOR_CONDITIONS = [
    ("pea_iso 40C/10min", PEA_LIPID_PROFILE, 40.0, 10.0),
    ("soy_iso 40C/10min", SOY_LIPID_PROFILE, 40.0, 10.0),
    ("pea_iso UHT 140C/0.1min", PEA_LIPID_PROFILE, 140.0, 0.1),
    ("pea_iso roasted 160C/30min", PEA_LIPID_PROFILE, 160.0, 30.0),
    ("spi/wg HME ~150C/3min", SOY_LIPID_PROFILE, 150.0, 3.0),
]


def core_shape_sweep() -> None:
    print("=== kinetic-core total_hydroperoxide vs max_conversion_fraction ===")
    print("(low-T anchors should be ~flat across the sweep; high-T should saturate)\n")
    caps = [0.25, 0.5, 1.0, 2.0]
    header = f"{'condition':30s}" + "".join(f"  cap={c:<6}" for c in caps)
    print(header)
    for label, profile, temp_c, t_min in ANCHOR_CONDITIONS:
        loads = []
        for cap in caps:
            # Reproduce the core extent at this cap (independent of the JSON value).
            from src.lipid_oxidation import _kinetics, _oxidation_rate_per_min
            rate = _oxidation_rate_per_min(temp_c, iron_ppm=profile.pro_oxidant_iron_ppm)
            progress = rate * t_min
            extent = _saturated_extent(progress, cap)
            lin = profile.linoleic_acid_pct / 100.0
            tl = profile.total_lipid_pct / 100.0
            loads.append(extent * lin * tl * _kinetics()["hydroperoxide_scale"])
        row = f"{label:30s}" + "".join(f"  {v:10.3g}" for v in loads)
        print(row)
    print()


def evaluator_guard() -> None:
    try:
        from scripts.diagnose_lipid_bias import main as diag_main
    except Exception as exc:  # pragma: no cover - diagnostic convenience
        print(f"(skipping evaluator guard: {exc})")
        return
    diag_main()


if __name__ == "__main__":
    core_shape_sweep()
    evaluator_guard()
