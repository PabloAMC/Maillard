import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pea_soy_mixed_external_package import (  # noqa: E402
    build_pea_soy_mixed_external_package_artifact,
    render_pea_soy_mixed_external_package_markdown,
)


def test_pea_soy_mixed_external_package_tracks_both_matrix_lanes():
    payload = build_pea_soy_mixed_external_package_artifact()

    assert payload["summary"]["matrix_count"] == 2
    assert payload["summary"]["package_status"] == "specified_not_yet_measured"
    matrices = {row["matrix"]: row for row in payload["matrices"]}
    assert set(matrices) == {"pea_iso", "soy_iso"}
    assert "external_quantitative_origin" in matrices["pea_iso"]["would_close_requirements"]
    assert matrices["soy_iso"]["current_external_anchor"] == "soy_isolate_40C_PratapSingh2021"


def test_pea_soy_mixed_external_package_markdown_surfaces_dual_bundle():
    markdown = render_pea_soy_mixed_external_package_markdown(build_pea_soy_mixed_external_package_artifact())

    assert "Pea Soy Mixed External Package" in markdown
    assert "Pea isolate mixed meaty-positive package" in markdown
    assert "Soy isolate mixed meaty-positive package" in markdown