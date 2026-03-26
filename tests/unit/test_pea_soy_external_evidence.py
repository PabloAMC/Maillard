import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pea_soy_external_evidence import (  # noqa: E402
    build_pea_soy_external_evidence_artifact,
    render_pea_soy_external_evidence_markdown,
)


def test_pea_soy_external_evidence_marks_current_external_surface_as_insufficient_for_promotion():
    payload = build_pea_soy_external_evidence_artifact()

    assert payload["summary"]["lane_count"] == 2
    assert payload["summary"]["external_mixed_meaty_positive_ready"] == 0
    assert payload["summary"]["promotion_ready_today"] == 0
    by_protein = {row["protein_type"]: row for row in payload["lanes"]}
    assert by_protein["pea_iso"]["external_target_profile"] == "adverse_only"
    assert by_protein["soy_iso"]["external_target_profile"] == "adverse_only"
    assert by_protein["pea_iso"]["required_external_package"]["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"


def test_pea_soy_external_evidence_markdown_surfaces_missing_mixed_external_package():
    markdown = render_pea_soy_external_evidence_markdown(build_pea_soy_external_evidence_artifact())

    assert "Pea Soy External Evidence" in markdown
    assert "Mixed Meaty-Positive External Present" in markdown
    assert "land_external_quantitative_mixed_matrix_benchmark_package" in markdown