import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dha_lysinoalanine_external_package import (  # noqa: E402
    build_dha_lysinoalanine_external_package_artifact,
    render_dha_lysinoalanine_external_package_markdown,
)


def test_dha_lysinoalanine_external_package_specifies_two_matrix_bundle():
    payload = build_dha_lysinoalanine_external_package_artifact()

    assert payload["summary"]["package_status"] == "specified_not_yet_measured"
    assert payload["summary"]["sprint_name"] == "dha_lysinoalanine_external_package"
    assert payload["summary"]["matrix_count"] == 2
    by_matrix = {row["matrix"]: row for row in payload["matrices"]}
    assert "Lysinoalanine (LAL)" in by_matrix["pea_iso"]["direct_damage_targets"]
    assert "direct_crosslink_marker_external_quantified" in by_matrix["soy_iso"]["would_close_requirements"]
    assert by_matrix["soy_iso"]["supportive_anchor_ids"] == ["de_leyn_2019_thiamine_retention"]


def test_dha_lysinoalanine_external_package_markdown_surfaces_bundle_and_policy():
    markdown = render_dha_lysinoalanine_external_package_markdown(build_dha_lysinoalanine_external_package_artifact())

    assert "DHA Lysinoalanine External Package" in markdown
    assert "Global Measurement Bundle" in markdown
    assert "Reactive lysine fraction" in markdown
    assert "Lysinoalanine (LAL)" in markdown
    assert "specified_not_yet_measured" in markdown