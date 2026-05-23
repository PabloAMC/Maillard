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
    bench_labels, comp_labels, voi, outside, planned_mask = generate_gap_heatmap.build_grid(payload)
    # Since PLANNED_BENCHMARKS has 5 elements, B2 should be at index 5 (first after the planned ones)
    assert bench_labels[5] == "B2"
    assert comp_labels[0] == "MFT"
    # Outside-CI marker preserved.
    assert bool(outside[5, 0]) is True


def test_render_writes_png(tmp_path: Path):
    from scripts.generators import generate_gap_heatmap

    payload = {
        "candidates": [
            {"benchmark_id": "B1", "compound": "MFT", "voi_score": 1.0, "inside_ci": True},
            {"benchmark_id": "B2", "compound": "MFT", "voi_score": 7.0, "inside_ci": False},
        ]
    }
    bench_labels, comp_labels, voi, outside, planned_mask = generate_gap_heatmap.build_grid(payload)
    out = tmp_path / "heatmap.png"
    generate_gap_heatmap.render(bench_labels, comp_labels, voi, outside, planned_mask, out)
    assert out.exists()
    assert out.stat().st_size > 1000  # non-trivial PNG


def test_generate_html_briefs_writes_html(tmp_path: Path):
    from scripts.generators import generate_gap_heatmap

    payload = {
        "candidates": [
            {
                "benchmark_id": "cys_glucose_150C_Farmer1999",
                "compound": "2-methyl-3-furanthiol",
                "voi_score": 7.7,
                "inside_ci": True,
                "matrix_family": "free",
                "measured_ppb": 15.0,
                "predicted_p5": 0.004,
                "predicted_p50": 31.7,
                "predicted_p95": 586.9,
                "rank": 1,
                "rationale": "High value gap",
                "suggested_doe_template": "blocking_benchmark_gap"
            }
        ]
    }
    out = tmp_path / "briefs.html"
    generate_gap_heatmap.generate_html_briefs(payload, out)
    assert out.exists()
    content = out.read_text(encoding="utf-8")
    assert "2-methyl-3-furanthiol" in content
    assert "Sulfur-Maillard (MFT) Pathway" in content
    assert "blocking_benchmark_gap" in content
    assert "D-Ribose/D-Glucose" in content


def test_build_grid_keeps_distinct_compounds_when_short_labels_collide():
    from scripts.generators import generate_gap_heatmap

    payload = {
        "candidates": [
            {
                "benchmark_id": "cml_cel_commercial_pbma_Foods2023",
                "compound": "Nε-(Carboxymethyl)lysine (CML)",
                "voi_score": 1.5,
                "inside_ci": False,
            },
            {
                "benchmark_id": "cml_cel_commercial_pbma_Foods2023",
                "compound": "Nε-(Carboxyethyl)lysine (CEL)",
                "voi_score": 1.2,
                "inside_ci": False,
            },
        ]
    }

    bench_labels, comp_labels, voi, outside, planned_mask = generate_gap_heatmap.build_grid(payload)

    assert bench_labels[5] == "cml_cel_commercial_pbma_Foods2023"
    assert len(comp_labels) == 2
    assert comp_labels[0] == "Nε-(Carboxymethyl)lysine (CML)"
    assert comp_labels[1] == "Nε-(Carboxyethyl)lysine (CEL)"
    assert voi.shape == (6, 2)
    assert int(outside.sum()) == 2

