import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_promotion_contract_artifact
from src.presentation import render_matrix_promotion_contract_markdown


def test_matrix_promotion_contract_exposes_explicit_rule_and_selects_primary_target_lane(matrix_promotion_contract_payload):
    payload = matrix_promotion_contract_payload

    assert payload["promotion_rule"]["contract_id"] == "matrix_external_decision_ready_v1"
    assert payload["promotion_rule"]["minimum_quantitative_closed_targets"] == 2
    assert payload["selected_promotion_target"]["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"

    by_id = {row["benchmark_id"]: row for row in payload["benchmarks"]}
    pea_candidate = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    assert pea_candidate["promotion_ready"] is False
    assert any(req["key"] == "external_quantitative_origin" and req["passed"] is False for req in pea_candidate["requirements"])

    markdown = render_matrix_promotion_contract_markdown(payload)
    assert "Matrix Promotion Contract" in markdown
    assert "matrix_external_decision_ready_v1" in markdown
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in markdown