import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.scope_gap_guard import (  # noqa: E402
    build_scope_gap_guard_artifact,
    render_scope_gap_guard_markdown,
)


def test_scope_gap_guard_blocks_scope_gap_families():
    payload = build_scope_gap_guard_artifact()

    assert payload["summary"]["blocked_family_count"] == 2
    assert payload["summary"]["blocked_families"] == ["coconut_oil_co_matrix", "other_plant_proteins"]


def test_scope_gap_guard_markdown_surfaces_guard_decision():
    markdown = render_scope_gap_guard_markdown(build_scope_gap_guard_artifact())

    assert "Scope Gap Guard" in markdown
    assert "blocked_from_active_expansion" in markdown