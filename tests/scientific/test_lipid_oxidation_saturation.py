"""S27 Workstream A — substrate-limited saturation of the lipid-oxidation kinetics.

Verifies that the first-order conversion cap (a) preserves the calibrated
low-temperature / short-time regime (linear limit) and (b) bounds the runaway
high temperature x time extrapolation that previously produced absurd hold-out
over-predictions, without breaking the matrix-benchmark pins.
"""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.lipid_oxidation import (
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    _saturated_extent,
    lipid_profile_for,
    predict_hexanal_generation,
)


def test_saturated_extent_linear_in_low_progress_regime():
    # progress << max_conversion -> extent ~= progress (<0.5% deviation),
    # so previously calibrated low-T anchors are preserved.
    for progress in (1e-5, 1e-4, 1e-3):
        extent = _saturated_extent(progress, max_conversion=1.0)
        assert math.isclose(extent, progress, rel_tol=5e-3)


def test_saturated_extent_caps_at_max_conversion():
    # Deep extrapolation cannot exceed the physical 100% conversion bound.
    assert _saturated_extent(1e6, max_conversion=1.0) <= 1.0
    assert _saturated_extent(50.0, max_conversion=0.5) <= 0.5
    # Monotonic but bounded.
    assert _saturated_extent(5.0, 1.0) > _saturated_extent(1.0, 1.0)
    assert _saturated_extent(5.0, 1.0) < 1.0


def test_saturated_extent_handles_degenerate_inputs():
    assert _saturated_extent(0.0, 1.0) == 0.0
    assert _saturated_extent(-1.0, 1.0) == 0.0
    # max_conversion <= 0 falls back to the raw (clamped) progress.
    assert _saturated_extent(3.0, 0.0) == 3.0


def test_default_path_is_byte_identical_to_linear_model():
    # The cap ships DISABLED (max_conversion_fraction=null), so the production
    # hydroperoxide term must EXACTLY equal the pre-S27 linear computation, in the
    # same multiplication order.
    #
    # DRIFT GUARD, not a coverage target. The trace lipid compounds carry near-zero-width
    # CIs, so they flip inside/outside on any perturbation of this term -- including a
    # bit-level float reassociation that changes no physics. This test exists to catch that
    # silent drift, so that a change in the in-panel coverage number always corresponds to a
    # real modelling change and never to an accidental refactor of the arithmetic here.
    #
    # It deliberately does NOT pin a coverage figure. The number moves for legitimate
    # reasons (the panel shrank from 19 to 16 benchmarks on 2026-08-26 when three benchmarks
    # with unlocatable sources were quarantined; in-panel coverage went 37/48 -> 35/41), and
    # pinning it here would turn an observation into a target to be defended. Regenerate
    # results/validation/prediction_uncertainty.md for the current value.
    from src.lipid_oxidation import _kinetics, _oxidation_rate_per_min

    assert _kinetics()["max_conversion_fraction"] is None, "cap must be disabled by default"
    for temp_c, t_min in ((40.0, 10.0), (100.0, 45.0), (160.0, 30.0)):
        out = predict_hexanal_generation(PEA_LIPID_PROFILE, temp_C=temp_c, time_min=t_min)
        rate = _oxidation_rate_per_min(temp_c, iron_ppm=PEA_LIPID_PROFILE.pro_oxidant_iron_ppm)
        linear = rate * (PEA_LIPID_PROFILE.linoleic_acid_pct / 100.0) * (
            PEA_LIPID_PROFILE.total_lipid_pct / 100.0
        ) * t_min * _kinetics()["hydroperoxide_scale"]
        assert out["total_hydroperoxide"] == linear


def test_conversion_factor_is_identity_when_disabled_and_bounds_when_enabled():
    from src.lipid_oxidation import _conversion_factor

    # Disabled (None / <= 0) => exactly 1.0, regardless of magnitude.
    assert _conversion_factor(1.0, 30.0, None) == 1.0
    assert _conversion_factor(1.0, 30.0, 0.0) == 1.0
    # Enabled => < 1.0 for a large progress, and ~1.0 for a tiny progress.
    assert _conversion_factor(0.05, 30.0, 1.0) < 1.0          # progress 1.5, capped
    assert math.isclose(_conversion_factor(1e-6, 10.0, 1.0), 1.0, rel_tol=1e-3)  # tiny progress


def test_lipid_profile_registry_round_trips_known_matrices():
    pea = lipid_profile_for("pea_iso")
    assert pea is not None and pea.linoleic_acid_pct == PEA_LIPID_PROFILE.linoleic_acid_pct
    soy = lipid_profile_for("soy_iso")
    assert soy is not None and soy.pro_oxidant_iron_ppm == SOY_LIPID_PROFILE.pro_oxidant_iron_ppm
    assert lipid_profile_for("unregistered_matrix") is None
    assert lipid_profile_for(None) is None
