import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.hexanal_nonanal_calibration import (  # noqa: E402
    build_hexanal_nonanal_calibration_artifact,
    render_hexanal_nonanal_calibration_markdown,
)


def test_hexanal_nonanal_calibration_closure_passes_for_protocol_pilot_lanes():
    payload = build_hexanal_nonanal_calibration_artifact()

    assert payload["summary"]["lane_count"] == 2
    assert payload["summary"]["closed_lane_count"] == 2
    assert payload["summary"]["hazard_lane_count"] == 0

    by_protein = {row["protein_type"]: row for row in payload["lanes"]}
    pea_row = by_protein["pea_iso"]
    soy_row = by_protein["soy_iso"]

    assert pea_row["closure_action"] == "calibration_closed"
    assert soy_row["closure_action"] == "calibration_closed"
    assert pea_row["compounds"]["Hexanal"]["in_bounds"] is True
    assert pea_row["compounds"]["Nonanal"]["in_bounds"] is True
    assert soy_row["compounds"]["Hexanal"]["ratio"] == 1.0
    assert soy_row["compounds"]["Nonanal"]["ratio"] == 1.0


def test_hexanal_nonanal_calibration_markdown_surfaces_internal_closure_policy():
    markdown = render_hexanal_nonanal_calibration_markdown(build_hexanal_nonanal_calibration_artifact())

    assert "Hexanal Nonanal Calibration Closure" in markdown
    assert "calibration_closed" in markdown
    assert "Closed lanes: 2" in markdown
    assert "still do not count as external promotion claims" not in markdown
    assert "retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence" in markdown