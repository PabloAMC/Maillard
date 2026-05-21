from __future__ import annotations

import pytest

from src.matrix_recalibration import assess_refit_candidate


def _uncertainty_payload(*, inside: bool, residual: float, trust_hits: int, trust_total: int):
    return {
        "summary": {
            "ci_coverage_hits": trust_hits,
            "matched_compound_count": trust_total,
            "ci_coverage_rate": trust_hits / trust_total,
        },
        "benchmarks": [
            {
                "benchmark_id": "bench-1",
                "compounds": [
                    {
                        "compound": "Hexanal",
                        "measured_ppb": 10.0,
                        "predicted_p50": 10.0 * (10 ** residual),
                        "inside_ci": inside,
                    }
                ],
            }
        ],
    }


def test_assess_refit_candidate_accepts_only_with_improved_residual_and_no_ci_regression():
    baseline = _uncertainty_payload(inside=True, residual=0.50, trust_hits=39, trust_total=48)
    candidate = _uncertainty_payload(inside=True, residual=0.20, trust_hits=39, trust_total=48)

    decision = assess_refit_candidate(baseline, candidate, baseline, candidate)

    assert decision["accepted"] is True
    assert decision["ci_regressions"] == 0
    assert decision["candidate_median_abs_log10_residual"] == pytest.approx(0.20)


def test_assess_refit_candidate_rejects_ci_regression_even_when_residual_improves():
    baseline = _uncertainty_payload(inside=True, residual=0.50, trust_hits=39, trust_total=48)
    candidate = _uncertainty_payload(inside=False, residual=0.20, trust_hits=39, trust_total=48)

    decision = assess_refit_candidate(baseline, candidate, baseline, candidate)

    assert decision["accepted"] is False
    assert decision["ci_regressions"] == 1