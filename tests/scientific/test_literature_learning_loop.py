import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_learning_loop import (
    build_literature_learning_loop_payload,
    render_literature_learning_loop_markdown,
)


def test_learning_loop_payload_links_ready_references_to_runtime_artifacts():
    payload = build_literature_learning_loop_payload(ROOT)
    rows = {row["id"]: row for row in payload["ready_reference_rows"]}

    assert payload["summary"]["ready_reference_count"] == 5
    assert rows["trikusuma_2019"]["encoding_status"] == "encoded_runtime_artifact"
    assert rows["lincoln_2025"]["encoding_status"] == "encoded_runtime_artifact"
    assert any(item["artifact_id"] == "pea_isolate_uht_140C_Trikusuma2019" for item in rows["trikusuma_2019"]["runtime_artifacts"])
    assert any(item["artifact_id"] == "asen_2022_pea_denaturation" for item in rows["asen_2022"]["runtime_artifacts"])


def test_learning_loop_reviews_matrix_priors_and_structural_gaps():
    payload = build_literature_learning_loop_payload(ROOT)
    prior_rows = {row["protein_type"]: row for row in payload["matrix_prior_review"]}

    assert prior_rows["myco"]["has_accessibility_window"] is True
    assert "directional_only" in prior_rows["myco"]["uncertainty_postures"]
    assert "extrusion_structured" in prior_rows["myco"]["process_state_applicability"]
    assert any(row["gap_id"] == "ppi_meaty_positive_matrix_benchmark" for row in payload["intake_structural_gap_review"])
    assert any(row["gap_id"] == "intact_spi_ppi_quantified_mft_fft" for row in payload["process_gap_review"])

    markdown = render_literature_learning_loop_markdown(payload)
    assert "Literature Learning Loop" in markdown
    assert "Matrix Prior Review" in markdown
    assert "myco" in markdown