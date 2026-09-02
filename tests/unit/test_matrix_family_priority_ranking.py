import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_family_priority_ranking import (  # noqa: E402
    build_matrix_family_priority_ranking_artifact,
    render_matrix_family_priority_ranking_markdown,
)


def test_matrix_family_priority_ranking_orders_next_scope_choices():
    payload = build_matrix_family_priority_ranking_artifact()

    by_family = {row["matrix_family"]: row for row in payload["families"]}
    assert by_family["pea_isolate"]["scope_priority"] == "active_matrix_priority"
    assert by_family["soy_isolate"]["scope_priority"] == "active_matrix_priority"
    assert by_family["mycoprotein"]["scope_priority"] == "bounded_next_candidate"
    assert by_family["coconut_oil_co_matrix"]["scope_priority"] == "scope_gap_to_rank"
    assert by_family["other_plant_proteins"]["scope_priority"] == "scope_gap_to_rank"

    families_in_order = [row["matrix_family"] for row in payload["families"]]
    assert families_in_order.index("mycoprotein") < families_in_order.index("coconut_oil_co_matrix")
    assert families_in_order.index("coconut_oil_co_matrix") < families_in_order.index("other_plant_proteins")
    assert "mycoprotein" in payload["summary"]["bounded_next_candidates"]
    assert "coconut_oil_co_matrix" in payload["summary"]["scope_gap_priorities"]
    assert by_family["pea_isolate"]["evidence_landing"] == "external_benchmark_plus_calibration_review"
    assert by_family["mycoprotein"]["evidence_landing"] == "bounded_calibration_prior"
    assert by_family["coconut_oil_co_matrix"]["evidence_landing"] == "family_specific_calibration_or_tradeoff_benchmark"


