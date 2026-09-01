import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.objective_progress import (  # noqa: E402
    build_objective_progress_artifact,
    render_objective_progress_markdown,
)


def test_objective_progress_summarizes_external_gaps():
    payload = build_objective_progress_artifact()

    assert payload["summary"]["objective_count"] == 2
    by_id = {row["objective_id"]: row for row in payload["objectives"]}
    assert by_id["external_mixed_meaty_positive_package"]["closed_count"] == 0
    assert by_id["extrusion_direct_damage_package"]["closed_count"] == 0
    assert payload["selected_next_family_sprint"] == "dha_lysinoalanine_external_package"


def test_objective_progress_markdown_surfaces_prediction_effects():
    markdown = render_objective_progress_markdown(build_objective_progress_artifact())

    assert "Objective Progress" in markdown
    assert "External mixed meaty-positive package" in markdown
    assert "Pea Soy Mixed External Package Delta" in markdown
    assert "Extrusion Direct Closure Delta" in markdown
    assert "DHA Lysinoalanine External Package Delta" in markdown
    assert "dha_lysinoalanine_external_package" in markdown