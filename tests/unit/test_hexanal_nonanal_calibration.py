import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.hexanal_nonanal_calibration import (  # noqa: E402
    build_hexanal_nonanal_calibration_artifact,
    render_hexanal_nonanal_calibration_markdown,
)


def test_hexanal_nonanal_calibration_surfaces_prediction_change_cascade():
    payload = build_hexanal_nonanal_calibration_artifact()

    assert payload["summary"]["closed_marker_count"] == 4
    assert payload["summary"]["marker_count"] == 4
    assert len(payload["prediction_change_cascade"]) == 4
    first = payload["prediction_change_cascade"][0]
    assert first["closure_state"] == "calibration_closed"
    assert first["next_decision_gate"] == "retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence"


def test_hexanal_nonanal_markdown_includes_prediction_validation_chain():
    markdown = render_hexanal_nonanal_calibration_markdown(build_hexanal_nonanal_calibration_artifact())

    assert "Prediction Validation Chain" in markdown
    assert "Closed markers: 4 / 4" in markdown