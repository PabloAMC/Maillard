import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_family_next_action import (  # noqa: E402
    build_matrix_family_next_action_artifact,
    render_matrix_family_next_action_markdown,
)


def test_matrix_family_next_action_advances_only_one_bounded_family():
    payload = build_matrix_family_next_action_artifact()

    assert payload["summary"]["chosen_family"] == "mycoprotein"
    assert payload["summary"]["chosen_action"] == "advance_bounded_next_family"

    by_family = {row["matrix_family"]: row for row in payload["decisions"]}
    assert by_family["mycoprotein"]["decision"] == "advance_now"
    assert by_family["coconut_oil_co_matrix"]["decision"] == "defer_as_scope_gap"
    assert by_family["other_plant_proteins"]["decision"] == "defer_as_scope_gap"
    assert by_family["pea_isolate"]["decision"] == "defer_until_primary_matrix_lane_moves"


