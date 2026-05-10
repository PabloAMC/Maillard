"""Tests for S22.3 gap heatmap renderer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def test_build_grid_orders_by_max_voi(tmp_path: Path):
    from scripts.generators import generate_gap_heatmap

    payload = {
        "candidates": [
            {"benchmark_id": "B1", "compound": "MFT", "voi_score": 1.0, "inside_ci": True},
            {"benchmark_id": "B1", "compound": "Hexanal", "voi_score": 0.5, "inside_ci": False},
            {"benchmark_id": "B2", "compound": "MFT", "voi_score": 7.0, "inside_ci": False},
        ]
    }
    bench_labels, comp_labels, voi, outside = generate_gap_heatmap.build_grid(payload)
    assert bench_labels[0] == "B2"  # higher max VoI first
    assert comp_labels[0] == "MFT"
    # Outside-CI marker preserved.
    assert bool(outside[0, 0]) is True


def test_render_writes_png(tmp_path: Path):
    from scripts.generators import generate_gap_heatmap

    payload = {
        "candidates": [
            {"benchmark_id": "B1", "compound": "MFT", "voi_score": 1.0, "inside_ci": True},
            {"benchmark_id": "B2", "compound": "MFT", "voi_score": 7.0, "inside_ci": False},
        ]
    }
    bench_labels, comp_labels, voi, outside = generate_gap_heatmap.build_grid(payload)
    out = tmp_path / "heatmap.png"
    generate_gap_heatmap.render(bench_labels, comp_labels, voi, outside, out)
    assert out.exists()
    assert out.stat().st_size > 1000  # non-trivial PNG
