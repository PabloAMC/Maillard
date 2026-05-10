"""Unit tests for validate_ts_mode (no xTB calls; pure linear algebra)."""
from __future__ import annotations

import numpy as np
import pytest

from src.xtb_backend import validate_ts_mode


def test_validate_ts_mode_concentrated_on_expected():
    """A mode whose motion is fully on atoms 0,1 passes for expected=[0,1]."""
    v = np.zeros(9)
    v[0:3] = [1.0, 0.0, 0.0]
    v[3:6] = [-1.0, 0.0, 0.0]
    out = validate_ts_mode(v, expected_atoms=[0, 1])
    assert out["pass"] is True
    assert abs(out["concentration"] - 1.0) < 1e-9
    # Top atoms should rank atom 0 and 1 first
    top_idx = [a for a, _ in out["top_atoms"][:2]]
    assert set(top_idx) == {0, 1}


def test_validate_ts_mode_dispersed_fails():
    """A mode evenly spread across 10 atoms fails when only 1 is expected."""
    n = 10
    v = np.ones(3 * n) / np.sqrt(3 * n)
    out = validate_ts_mode(v, expected_atoms=[3], threshold=0.3)
    assert out["pass"] is False
    # 1 atom out of 10 evenly spread → ~0.10
    assert abs(out["concentration"] - 0.10) < 1e-2


def test_validate_ts_mode_threshold_boundary():
    n = 5
    v = np.zeros(3 * n)
    v[0:3] = [0.5, 0.5, 0.0]   # atom 0 carries half the squared-norm mass
    v[3:6] = [0.5, 0.5, 0.0]   # atom 1 carries the other half
    out = validate_ts_mode(v, expected_atoms=[0], threshold=0.4)
    assert abs(out["concentration"] - 0.5) < 1e-9
    assert out["pass"] is True


def test_validate_ts_mode_zero_vector():
    out = validate_ts_mode(np.zeros(9), expected_atoms=[0])
    assert out["pass"] is False
    assert out["concentration"] == 0.0


def test_validate_ts_mode_invalid_shape():
    with pytest.raises(ValueError):
        validate_ts_mode(np.zeros(7), expected_atoms=[0])


def test_validate_ts_mode_ignores_out_of_range_atoms():
    v = np.zeros(9)
    v[0:3] = [1.0, 0.0, 0.0]
    out = validate_ts_mode(v, expected_atoms=[0, 99])  # 99 silently ignored
    assert abs(out["concentration"] - 1.0) < 1e-9
