import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mycoprotein_reference import (  # noqa: E402
    build_mycoprotein_reference_artifact,
    render_mycoprotein_reference_markdown,
)


def test_mycoprotein_reference_uses_bounded_priors_and_next_family_decision():
    payload = build_mycoprotein_reference_artifact()

    assert payload["summary"]["matrix_family"] == "mycoprotein"
    assert payload["summary"]["decision"] == "advance_now"
    assert payload["summary"]["evidence_surface"] == "bounded_calibration_prior"
    assert payload["reference_windows"]["denaturation"]["midpoint_celsius"] == 78.0


def test_mycoprotein_reference_markdown_surfaces_bounded_policy():
    markdown = render_mycoprotein_reference_markdown(build_mycoprotein_reference_artifact())

    assert "Mycoprotein Reference" in markdown
    assert "bounded_calibration_prior" in markdown
    assert "Decision: advance_now" in markdown