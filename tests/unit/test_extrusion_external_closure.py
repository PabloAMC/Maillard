import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extrusion_external_closure import (  # noqa: E402
    build_extrusion_external_closure_artifact,
    render_extrusion_external_closure_markdown,
)


def test_extrusion_external_closure_artifact_tracks_root_blocker_and_next_sprint():
    payload = build_extrusion_external_closure_artifact()

    assert payload["summary"]["matrix_count"] == 2
    assert payload["summary"]["root_blocker"] == "external_quantitative_direct_damage_markers_under_extrusion"
    assert payload["summary"]["selected_next_family_sprint"] == "dha_lysinoalanine_external_package"
    by_matrix = {row["matrix"]: row for row in payload["matrices"]}
    assert by_matrix["pea_iso"]["direct_closure_ready"] is False
    assert "Reactive lysine fraction" in by_matrix["pea_iso"]["direct_markers_missing"]
    assert by_matrix["soy_iso"]["supportive_anchor_ids"] == ["de_leyn_2019_thiamine_retention"]
    assert by_matrix["soy_iso"]["contextual_anchor_ids"] == ["troise_2018_soy_thermal_history"]


