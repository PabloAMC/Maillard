import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_experiment_intake import build_matrix_experiment_support_delta_artifact, load_matrix_experiment_intake_schema
from src.presentation import render_matrix_experiment_support_delta_markdown


def test_matrix_experiment_intake_support_delta_strengthens_external_matrix_evidence():
    schema = load_matrix_experiment_intake_schema()
    assert schema["contract_id"] == "matrix_experiment_intake_v1"

    payload = build_matrix_experiment_support_delta_artifact(
        ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml"
    )

    assert payload["aligned_benchmark"]["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"
    assert payload["promotion_assessment"]["promotion_ready_after"] is False
    assert payload["promotion_assessment"]["readiness_change"] == "evidence_strengthened_not_yet_promotable"
    assert payload["promotion_assessment"]["promotion_blocker_after"] == "depends on internal or transferred support"
    assert payload["support_delta"]["delta_counts"]["strengthened"] >= 1

    compounds = {row["compound"]: row for row in payload["compounds"]}
    assert compounds["2-furfurylthiol"]["support_status"] == "quantitative_closed"
    assert compounds["bis(2-methyl-3-furyl) disulfide"]["support_status"] == "directional_support"

    markdown = render_matrix_experiment_support_delta_markdown(payload)
    assert "Matrix Experiment Support Delta" in markdown
    assert "evidence_strengthened_not_yet_promotable" in markdown