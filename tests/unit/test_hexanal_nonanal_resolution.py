import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.hexanal_nonanal_resolution import (  # noqa: E402
    build_hexanal_nonanal_resolution_artifact,
    render_hexanal_nonanal_resolution_markdown,
)


def test_hexanal_nonanal_resolution_keeps_internal_closure_separate_from_promotion():
    payload = build_hexanal_nonanal_resolution_artifact()

    assert payload["summary"]["marker_count"] == 4
    assert payload["summary"]["reduced_marker_count"] == 4
    assert all(row["promotion_effect"] == "no_external_promotion_upgrade" for row in payload["markers"])


def test_hexanal_nonanal_resolution_markdown_surfaces_remaining_external_blocker():
    markdown = render_hexanal_nonanal_resolution_markdown(build_hexanal_nonanal_resolution_artifact())

    assert "Hexanal Nonanal Resolution" in markdown
    assert "Remaining External Blocker" in markdown