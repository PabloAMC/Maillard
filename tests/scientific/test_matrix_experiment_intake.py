import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_experiment_intake import (
    build_matrix_experiment_support_delta_artifact,
    load_matrix_experiment_intake,
    load_matrix_experiment_intake_schema,
)
from src.presentation import render_matrix_experiment_support_delta_markdown


def test_matrix_experiment_intake_support_delta_strengthens_external_matrix_evidence():
    schema = load_matrix_experiment_intake_schema()
    assert schema["contract_id"] == "matrix_experiment_intake_v1"

    payload = build_matrix_experiment_support_delta_artifact(
        ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml"
    )

    assert payload["aligned_benchmark"]["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"
    assert payload["promotion_assessment"]["promotion_ready_after"] is False
    assert payload["promotion_assessment"]["promotion_claim_allowed"] is True
    assert payload["promotion_assessment"]["readiness_change"] == "evidence_strengthened_not_yet_promotable"
    assert payload["promotion_assessment"]["promotion_blocker_after"] == "depends on internal or transferred support"
    assert payload["promotion_assessment"]["landing_recommendation"] == "land_in_benchmark_candidate_or_blocker_registry"
    assert payload["support_delta"]["delta_counts"]["strengthened"] >= 1

    compounds = {row["compound"]: row for row in payload["compounds"]}
    assert compounds["2-furfurylthiol"]["support_status"] == "quantitative_closed"
    assert compounds["bis(2-methyl-3-furyl) disulfide"]["support_status"] == "directional_support"

    markdown = render_matrix_experiment_support_delta_markdown(payload)
    assert "Matrix Experiment Support Delta" in markdown
    assert "Promotion claim allowed: yes" in markdown
    assert "evidence_strengthened_not_yet_promotable" in markdown


def test_internal_matrix_experiment_support_delta_stays_outside_promotion_claims():
    payload = load_matrix_experiment_intake(ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml")
    payload["experiment_id"] = "pea_iso_ribose_cysteine_internal_protocol_pilot_example"
    payload["source_kind"] = "internal_experiment"
    payload["provenance"]["origin"] = "internal_experiment"
    payload["provenance"]["source_reference"] = "internal:protocol_pilot_example"

    artifact = build_matrix_experiment_support_delta_artifact(payload)

    assert artifact["promotion_assessment"]["promotion_ready_after"] is False
    assert artifact["promotion_assessment"]["promotion_claim_allowed"] is False
    assert artifact["promotion_assessment"]["landing_recommendation"] == "land_in_internal_candidate_or_calibration_payload_not_promotion_claim"
    assert (
        artifact["promotion_assessment"]["promotion_claim_policy"]
        == "internal_measurements_strengthen_calibration_and_comparator_evidence_but_do_not_unlock_external_decision_ready_claims"
    )

    markdown = render_matrix_experiment_support_delta_markdown(artifact)
    assert "Promotion claim allowed: no" in markdown
    assert "internal_measurements_strengthen_calibration_and_comparator_evidence_but_do_not_unlock_external_decision_ready_claims" in markdown