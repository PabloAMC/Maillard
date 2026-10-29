"""src/kinetic_core/parallel.py (strategy M2, 2026-10-29): faster, and not one number moved."""
from __future__ import annotations

import numpy as np
import pytest

from src.kinetic_core import parallel

SWITCH = {"value": "default"}


def _square(x):
    return x * x


def _fail_on_three(x):
    if x == 3:
        raise ValueError("three")
    return x


def test_map_ordered_keeps_order_and_matches_serial():
    items = list(range(11))
    assert parallel.map_ordered(_square, items, n_workers=4) == [x * x for x in items]


def test_a_closure_and_the_callers_state_reach_the_workers():
    SWITCH["value"] = "set by the caller"
    offset = 7
    try:
        out = parallel.map_ordered(lambda x: (x + offset, SWITCH["value"]), [1, 2, 3], n_workers=3)
    finally:
        SWITCH["value"] = "default"
    assert out == [(8, "set by the caller"), (9, "set by the caller"), (10, "set by the caller")]


def test_a_worker_exception_is_raised_in_the_caller():
    with pytest.raises(ValueError, match="three"):
        parallel.map_ordered(_fail_on_three, [1, 2, 3, 4], n_workers=2)


def test_workers_default_to_one_and_read_the_environment(monkeypatch):
    monkeypatch.delenv("MAILLARD_WORKERS", raising=False)
    assert parallel.workers() == 1
    monkeypatch.setenv("MAILLARD_WORKERS", "5")
    assert parallel.workers() == 5
    monkeypatch.setenv("MAILLARD_WORKERS", "not a number")
    assert parallel.workers() == 1


def test_b23_systems_are_bit_identical_in_parallel(wave_generator, monkeypatch):
    """Three of B2.3's systems at its own vector, serial and on three workers: every array equal."""
    with wave_generator("generate_kinetic_core_b2_3_fit") as b23:
        keep = list(b23.SYSTEMS)[:3]
        monkeypatch.setattr(b23, "SYSTEMS", {k: b23.SYSTEMS[k] for k in keep})
        import generate_kinetic_core_b2_4_fit as b24
        x = b24.b23_vector()
        parameters = b23.build_parameters(x)
        drift = b23.unpack(x)[3]
        monkeypatch.setenv("MAILLARD_WORKERS", "1")
        serial = b23.run_systems(parameters, True, drift)
        monkeypatch.setenv("MAILLARD_WORKERS", "3")
        forked = b23.run_systems(parameters, True, drift)
    assert list(serial) == list(forked) == keep
    for name in keep:
        a, b = serial[name]["run"], forked[name]["run"]
        assert np.array_equal(a.concentrations, b.concentrations), name
        assert np.array_equal(a.times_min, b.times_min), name
        assert (a.ph_series is None) == (b.ph_series is None)
        if a.ph_series is not None:
            assert np.array_equal(a.ph_series, b.ph_series), name


def test_a_call_inside_a_worker_runs_serially_and_still_answers():
    """A worker is a daemon and cannot fork; the nested call (a Jacobian point integrating its pots) runs serially."""
    def inner(x):
        return sum(parallel.map_ordered(_square, [x, x + 1], n_workers=4)), parallel._IN_WORKER
    out = parallel.map_ordered(inner, [1, 2, 3], n_workers=3)
    assert out == [(1 + 4, True), (4 + 9, True), (9 + 16, True)]
    assert parallel._IN_WORKER is False
