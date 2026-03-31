import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_primary_benchmark_campaign import (  # noqa: E402
    build_matrix_primary_benchmark_campaign_artifact,
    render_matrix_primary_benchmark_campaign_markdown,
)


def test_matrix_primary_benchmark_campaign_targets_selected_mixed_lane():
    payload = build_matrix_primary_benchmark_campaign_artifact()

    assert payload["summary"]["selected_matrix"] == "pea_iso"
    assert payload["summary"]["primary_data_blocker"] == "external_quantitative_measured_volatiles_for_mixed_matrix_lane"

    pea_row = next(row for row in payload["arms"] if row["matrix"] == "pea_iso")
    assert "comparator_is_measured_volatiles" in pea_row["would_close_requirements"]
    assert "external_quantitative_origin" in pea_row["would_close_requirements"]
    assert "Hexanal" in pea_row["evidence_or_calibration_blockers"]
    assert pea_row["mechanistic_blockers"] == []
    assert pea_row["calibration_closure_action"] == "calibration_closed"
    assert pea_row["calibration_passed"] is True
    assert pea_row["hexanal_ratio"] is not None
    assert pea_row["nonanal_ratio"] is not None
    assert len(pea_row["required_desirable_targets"]) == 4


def test_matrix_primary_benchmark_campaign_markdown_surfaces_promotion_delta():
    markdown = render_matrix_primary_benchmark_campaign_markdown(build_matrix_primary_benchmark_campaign_artifact())

    assert "Matrix Primary Benchmark Campaign" in markdown
    assert "Promotion Delta" in markdown
    assert "Evidence/Calibration Blockers" in markdown
    assert "Calibration Route" in markdown
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in markdown