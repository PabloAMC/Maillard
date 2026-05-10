"""Phase-2 hardening: TS imaginary-frequency classifier."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.dft_refiner import (
    TS_IMAGINARY_FREQ_FLOOR_CM1,
    _classify_ts_frequencies,
)


def test_true_ts_one_significant_imaginary():
    v = _classify_ts_frequencies([-450.2, 120.0, 250.5, 800.1])
    assert v["is_true_saddle"] is True
    assert v["verdict"] == "true_ts"
    assert v["n_imaginary"] == 1
    assert v["n_significant_imaginary"] == 1
    assert v["most_negative_cm1"] == pytest.approx(-450.2)


def test_minimum_no_imaginary():
    v = _classify_ts_frequencies([12.3, 45.0, 120.0, 800.0])
    assert v["is_true_saddle"] is False
    assert v["verdict"] == "minimum"
    assert v["n_imaginary"] == 0
    assert v["n_significant_imaginary"] == 0
    assert v["most_negative_cm1"] is None


def test_weak_imaginary_below_floor_treated_as_invalid():
    # |f|=12 cm^-1 imaginary mode is residual translation/rotation noise.
    v = _classify_ts_frequencies([-12.0, 100.0, 200.0])
    assert v["is_true_saddle"] is False
    assert v["verdict"] == "weak_imaginary"
    assert v["n_imaginary"] == 1
    assert v["n_significant_imaginary"] == 0


def test_high_order_saddle_two_significant_imaginary():
    v = _classify_ts_frequencies([-600.0, -300.0, 100.0, 800.0])
    assert v["is_true_saddle"] is False
    assert v["verdict"] == "high_order_saddle"
    assert v["n_imaginary"] == 2
    assert v["n_significant_imaginary"] == 2


def test_no_frequencies():
    v = _classify_ts_frequencies(None)
    assert v["is_true_saddle"] is False
    assert v["verdict"] == "no_frequencies"
    v2 = _classify_ts_frequencies([])
    assert v2["verdict"] == "no_frequencies"


def test_floor_threshold_value():
    # Mode at exactly the floor is NOT counted as significant (strict >).
    v = _classify_ts_frequencies([-TS_IMAGINARY_FREQ_FLOOR_CM1, 100.0])
    assert v["n_imaginary"] == 1
    assert v["n_significant_imaginary"] == 0
    assert v["verdict"] == "weak_imaginary"
    # Just below floor (more negative) IS significant.
    v2 = _classify_ts_frequencies([-(TS_IMAGINARY_FREQ_FLOOR_CM1 + 1.0), 100.0])
    assert v2["n_significant_imaginary"] == 1
    assert v2["verdict"] == "true_ts"


def test_imaginary_freqs_sorted_ascending():
    v = _classify_ts_frequencies([-100.0, -500.0, -200.0, 300.0])
    assert v["imaginary_freqs_cm1"] == [-500.0, -200.0, -100.0]
    assert v["most_negative_cm1"] == -500.0


def test_serializable_to_json():
    v = _classify_ts_frequencies([-450.0, 200.0])
    s = json.dumps(v)
    rt = json.loads(s)
    assert rt["verdict"] == "true_ts"
    assert rt["is_true_saddle"] is True


def test_promotion_gate_uses_ts_validation_file(tmp_path: Path):
    """Smoke test: ts_validation.json with verdict != true_ts must be readable
    and recognizable as invalid by the script's gate logic."""
    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    v = _classify_ts_frequencies([12.0, 100.0])  # minimum
    (ckpt / "ts_validation.json").write_text(json.dumps(v))
    loaded = json.loads((ckpt / "ts_validation.json").read_text())
    assert loaded["is_true_saddle"] is False
    assert loaded["verdict"] == "minimum"
