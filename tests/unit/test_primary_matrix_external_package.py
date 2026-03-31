import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.primary_matrix_external_package import (  # noqa: E402
    build_primary_matrix_external_package_artifact,
    render_primary_matrix_external_package_markdown,
)


def test_primary_matrix_external_package_selects_one_priority_lane():
    payload = build_primary_matrix_external_package_artifact()

    assert payload["summary"]["selected_matrix"] == "pea_iso"
    assert payload["package"]["status"] == "specified_not_yet_measured"
    assert "external_quantitative_origin" in payload["package"]["would_close_requirements"]


def test_primary_matrix_external_package_markdown_surfaces_decision_delta():
    markdown = render_primary_matrix_external_package_markdown(build_primary_matrix_external_package_artifact())

    assert "Primary Matrix External Package" in markdown
    assert "Would close requirements" in markdown