"""
tests/unit/test_trust_visibility.py

P5 gate tests: verifies that trust signals (prediction posture banner,
calibration source per compound, [G] accessibility warning) fire correctly.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path
from typing import Any, Dict, List
from unittest.mock import MagicMock, patch

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.usability_reports import (
    _build_compound_confidence_rows,
    ConfidenceAssessment,
)
from src.presentation import render_decision_summary_cli, render_deep_explainability_cli


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_confidence_assessment(
    score: float = 90.0,
    tier: str = "high",
    neighborhood: str = "primary_free_precursor",
) -> ConfidenceAssessment:
    return ConfidenceAssessment(
        tier=tier,
        score=score,
        benchmark_neighborhood=neighborhood,
        support_ratio=1.0,
        avg_uncertainty_kcal=3.0,
        dominant_factors=[],
        recommended_posture="Suitable for quantitative prioritization.",
    )


def _make_result(
    *,
    projection_metadata: Dict[str, Any],
    confidence_metadata: Dict[str, Any],
    target_score: float = 50.0,
    off_flavour_risk: float = 0.0,
    safety_score: float = 0.0,
    targets: List[Dict[str, Any]] | None = None,
    avg_uncertainty: float = 3.0,
    precursor_contributions: Dict[str, float] | None = None,
    bottleneck_precursor: str = "none",
    bottleneck_severity: float = 0.0,
    suppressed_compounds: List | None = None,
    predicted_ppb: Dict[str, float] | None = None,
    radar: Dict | None = None,
    flagged_toxics: List | None = None,
) -> MagicMock:
    result = MagicMock()
    result.projection_metadata = projection_metadata
    result.confidence_metadata = confidence_metadata
    result.target_score = target_score
    result.off_flavour_risk = off_flavour_risk
    result.safety_score = safety_score
    result.targets = targets or []
    result.avg_uncertainty = avg_uncertainty
    result.precursor_contributions = precursor_contributions or {}
    result.bottleneck_precursor = bottleneck_precursor
    result.bottleneck_severity = bottleneck_severity
    result.suppressed_compounds = suppressed_compounds or []
    result.predicted_ppb = predicted_ppb or {}
    result.radar = radar or {}
    result.flagged_toxics = flagged_toxics or []
    return result


_FREE_PRECURSOR_META = {
    "CCCCCC=O": {
        "compound": "Hexanal",
        "proxy_ppb": 50.0,
        "observable_ppb": 45.0,
        "proxy_to_observable_ratio": 0.9,
        "matrix_factor": 0.95,
        "headspace_factor": 0.90,
        "volatile_class": "aldehyde",
        "process_state": "heated",
        "retention_runtime_mode": "dynamic",
        "calibration_source": "literature_direct",
        "calibration_evidence_strength": "externally_benchmarked",
        "calibration_fallback_mode": "none",
        "evidence_state": "externally_benchmarked",
        "target_class": "adverse_lipid_markers",
        "browning_index": 0.1,
    }
}

_MATRIX_LOW_META = {
    "CCCCCC=O": {
        **_FREE_PRECURSOR_META["CCCCCC=O"],
        "matrix_factor": 0.20,  # below 0.35 threshold → accessibility warning
        "headspace_factor": 0.25,
        "calibration_source": "class_fallback",
        "calibration_evidence_strength": "class_anchored",
    }
}


# ---------------------------------------------------------------------------
# P5 Test 1: free-precursor → QUANTITATIVE MODE banner
# ---------------------------------------------------------------------------

def test_prediction_posture_banner_quantitative(capsys):
    confidence_pkg = {
        "tier": "high",
        "score": 88.0,
        "benchmark_neighborhood": "primary_free_precursor",
        "recommended_posture": "Suitable for quantitative prioritization.",
        "dominant_factors": [],
        "decision_mode": "quantitative_recommendation",
        "prediction_mode": "benchmark_supported_quantitative",
        "calibration_diagnostics": {"summary": "Within envelope.", "extrapolation_axes": []},
        "compound_confidence": [],
        "sensitivity_summary": {},
    }
    result = _make_result(
        projection_metadata=_FREE_PRECURSOR_META,
        confidence_metadata=confidence_pkg,
    )
    render_decision_summary_cli(result, warnings=[])
    captured = capsys.readouterr().out
    assert "QUANTITATIVE MODE" in captured, f"Expected 'QUANTITATIVE MODE' in output:\n{captured}"
    assert "✅" in captured, f"Expected ✅ in output:\n{captured}"
    assert "Decision Mode    : quantitative_recommendation" in captured


# ---------------------------------------------------------------------------
# P5 Test 2: matrix (pea_iso) → DIRECTIONAL MODE banner
# ---------------------------------------------------------------------------

def test_prediction_posture_banner_directional(capsys):
    confidence_pkg = {
        "tier": "medium",
        "score": 68.0,
        "benchmark_neighborhood": "matrix_intake_only",
        "recommended_posture": "Verify absolute concentrations experimentally.",
        "dominant_factors": ["Plant-matrix support is still intake/headspace validated."],
        "decision_mode": "directional_hypothesis",
        "prediction_mode": "ranking_supported",
        "calibration_diagnostics": {"summary": "Class-level anchor applied.", "extrapolation_axes": ["benchmark_neighborhood"]},
        "compound_confidence": [],
        "sensitivity_summary": {},
    }
    result = _make_result(
        projection_metadata=_MATRIX_LOW_META,
        confidence_metadata=confidence_pkg,
    )
    render_decision_summary_cli(result, warnings=[])
    captured = capsys.readouterr().out
    assert "DIRECTIONAL MODE" in captured, f"Expected 'DIRECTIONAL MODE' in output:\n{captured}"
    assert "⚠️" in captured, f"Expected ⚠️ in output:\n{captured}"
    assert "Decision Mode    : directional_hypothesis" in captured
    assert "Extrapolation    : benchmark_neighborhood" in captured


# ---------------------------------------------------------------------------
# P5 Test 3: compound_confidence rows include calibration_source
# ---------------------------------------------------------------------------

def test_compound_confidence_includes_calibration_source():
    assessment = _make_confidence_assessment()
    result = _make_result(projection_metadata=_FREE_PRECURSOR_META, confidence_metadata={})
    rows = _build_compound_confidence_rows(result, assessment, top_n=3)
    assert rows, "Expected at least one compound confidence row"
    for row in rows:
        assert "calibration_source" in row, "calibration_source missing from row"
        assert "calibration_evidence_strength" in row, "calibration_evidence_strength missing from row"
        assert row["evidence_state"] == "externally_benchmarked"
        assert row["target_class"] == "adverse_lipid_markers"
        assert row["support_origin"] == "standard_matrix_support"
        assert row["reachability_status"] == "chemically_reachable"
        assert row["observable_assumption_summary"] == "dynamic | none | standard_matrix_support"
        assert isinstance(row["calibration_source"], str)
        assert row["calibration_source"] != ""


# ---------------------------------------------------------------------------
# P5 Test 4: [G] panel fires when matrix_factor < 0.35
# ---------------------------------------------------------------------------

def test_accessibility_warning_surfaces_when_matrix_factor_low(capsys):
    """[G] ACCESSIBILITY WARNING must appear when top compound has matrix_factor < 0.35."""
    assessment = _make_confidence_assessment(score=60.0, tier="medium", neighborhood="matrix_intake_only")
    result_for_rows = _make_result(projection_metadata=_MATRIX_LOW_META, confidence_metadata={})
    compound_rows = _build_compound_confidence_rows(result_for_rows, assessment, top_n=3)

    # Check that the row is flagged
    dominated = [r for r in compound_rows if r.get("accessibility_dominated")]
    assert dominated, "Expected at least one row with accessibility_dominated=True for matrix_factor=0.20"

    confidence_pkg = {
        "tier": "medium",
        "score": 60.0,
        "benchmark_neighborhood": "matrix_intake_only",
        "recommended_posture": "Use directionally.",
        "dominant_factors": [],
        "prediction_mode": "ranking_supported",
        "calibration_diagnostics": {"summary": "Class anchor."},
        "compound_confidence": compound_rows,
        "sensitivity_summary": {},
    }
    result = _make_result(
        projection_metadata=_MATRIX_LOW_META,
        confidence_metadata=confidence_pkg,
        predicted_ppb={"CCCCCC=O": 45.0},
    )
    render_deep_explainability_cli(result)
    captured = capsys.readouterr().out
    assert "[G] ACCESSIBILITY WARNING" in captured, (
        f"Expected '[G] ACCESSIBILITY WARNING' in output when matrix_factor=0.20:\n{captured}"
    )
    assert "Accessibility (not chemistry)" in captured
