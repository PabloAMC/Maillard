"""The route calibration's plumbing and de Roos's measured decay profile (prereg V2)."""
from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

import pytest

from src.kinetic_core import parameters_sulfur as ps

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

DE_ROOS_FIG4 = {3.0: 5.6e-3, 3.5: 5.3e-3, 4.0: 5.1e-3, 4.5: 4.8e-3, 5.0: 4.5e-3, 5.25: 5.1e-3, 5.5: 7.0e-3, 6.0: 1.43e-2}


@pytest.mark.parametrize("ph,k", sorted(DE_ROOS_FIG4.items()))
def test_decay_reproduces_de_roos_fig4_at_100C(ph, k):
    assert ps.CYS_AMADORI_ROUTE.k_decay(373.15, ph) == pytest.approx(k, rel=1e-9)


def test_decay_profile_is_u_shaped_and_held_below_ph3():
    route = ps.CYS_AMADORI_ROUTE
    assert route.k_decay(373.15, 4.0) > route.k_decay(373.15, 5.0) < route.k_decay(373.15, 6.0)
    assert route.k_decay(373.15, 2.0) == pytest.approx(route.k_decay(373.15, 3.0))


def test_above_measured_ph_alternatives():
    extend = ps.CYS_AMADORI_ROUTE
    flat = dataclasses.replace(extend, above_measured_ph="flat")
    assert flat.k_decay(373.15, 7.4) == pytest.approx(1.43e-2)
    slope = (math.log10(1.43e-2) - math.log10(7.0e-3)) / 0.5
    assert math.log10(extend.k_decay(373.15, 7.4)) == pytest.approx(math.log10(1.43e-2) + slope * 1.4)


@pytest.fixture(scope="module")
def calib():
    # Through wave_generator: importing it pulls in B64 -> B60 -> B2.3, which splice B2.3's tables;
    # imported bare, they stayed spliced for every later test (B2.4's control cost moved 22.56 -> 26.13).
    from tests.support import wave_generator
    with wave_generator("generate_kinetic_core_cysteine_route_calibration") as C:
        yield C


def test_rows_are_unique_and_counted(calib):
    keys = [(r["pot"], r["species"], r["minutes"]) for r in calib.ROWS]
    assert len(keys) == len(set(keys))
    counts = {s: sum(r["source"] == s for r in calib.ROWS) for s in calib.SOURCES}
    assert counts == {"cao2017": 24, "kang_s3": 18, "kang_s4": 12, "zhai2019": 64, "zhai2021": 3}


def test_start_vector_inside_bounds(calib):
    assert all(lo <= x <= hi for x, lo, hi in zip(calib.START, calib.LOWER, calib.UPPER))


def test_overrides_reach_the_parameters_and_switches_restore(calib):
    import numpy as np
    theta = calib.START.copy()
    theta[0] -= 1.0
    p = calib.overridden(calib._B60_PARAMS, theta)
    assert p["k_ttca_cys"].k_ref == pytest.approx(calib._B60_PARAMS["k_ttca_cys"].k_ref / 10.0)
    assert p["k_ttca_deg"].ea_kj_mol == pytest.approx(theta[1])
    before = (ps.CYS_AMADORI_ROUTE, ps.TTCA_FORMATION)
    with calib.installed(theta, "flat"):
        assert ps.CYS_AMADORI_ROUTE.enabled and ps.CYS_AMADORI_ROUTE.above_measured_ph == "flat"
    assert (ps.CYS_AMADORI_ROUTE, ps.TTCA_FORMATION) == before
    assert not ps.CYS_AMADORI_ROUTE.enabled
    assert np.isfinite(calib.START).all()
