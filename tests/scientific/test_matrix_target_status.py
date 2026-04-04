import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_target_status_artifact
from src.presentation import render_matrix_target_status_markdown


def test_matrix_target_status_distinguishes_quantitative_from_internal_and_directional_support(matrix_target_status_payload):
    payload = matrix_target_status_payload
    by_id = {row["benchmark_id"]: row for row in payload["benchmarks"]}

    pea_off = by_id["pea_isolate_40C_PratapSingh2021"]
    assert pea_off["support_counts"]["quantitative_closed"] >= 1
    assert pea_off["target_profile"] == "adverse_only"

    pea_meaty = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_meaty["support_counts"]["internal_candidate"] >= 1
    assert pea_meaty["target_profile"] == "mixed"
    assert pea_meaty["promotion_ready"] is False
    assert pea_meaty["evidence_or_calibration_priority_ready"] is True
    assert pea_meaty["mechanistic_priority_ready"] is False
    assert pea_meaty["blocker_class"] == "observable_or_calibration_blocker"
    assert pea_meaty["promotion_claim_posture"] == "not_a_promotion_lane"
    assert pea_meaty["next_best_action"] == "improve_observable_or_calibration"
    assert pea_meaty["best_computational_action"] == "improve_observable_or_calibration_before_qm"

    summary = payload["summary"]
    assert summary["quantitative_closed"] >= 1
    assert summary["internal_candidate"] >= 1
    assert summary["evidence_or_calibration_priority_ready"] >= 1

    markdown = render_matrix_target_status_markdown(payload)
    assert "Matrix Target Status" in markdown
    assert "quantitative_closed" in markdown
    assert "internal_candidate" in markdown
    assert "observable_or_calibration_blocker" in markdown
    assert "Evidence Origin" in markdown
    assert "Claim Posture" in markdown
    assert "Evidence/calibration-priority benchmarks" in markdown
    assert "improve_observable_or_calibration" in markdown
    assert "improve_observable_or_calibration_before_qm" in markdown
