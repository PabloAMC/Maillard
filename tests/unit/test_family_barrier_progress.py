import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_barrier_progress import (  # noqa: E402
    build_family_barrier_progress_artifact,
    render_family_barrier_progress_markdown,
)


def test_family_barrier_progress_distinguishes_explicit_and_non_explicit_lanes():
    payload = build_family_barrier_progress_artifact()

    assert payload["summary"]["family_count"] == 10
    by_family = {row["chemistry_family"]: row for row in payload["families"]}
    assert by_family["amino_acid_sugar_core"]["explicit_fast_barrier_count"] > 0
    assert by_family["amino_acid_sugar_core"]["arrhenius_anchor_count"] > 0
    assert by_family["alternative_protein_matrix_scope"]["barrier_lane_stage"] == "no_explicit_barrier_lane"
    assert by_family["fermentation_pretreatment"]["barrier_lane_stage"] == "no_explicit_barrier_lane"


def test_family_barrier_progress_markdown_surfaces_barrier_and_prediction_readout():
    markdown = render_family_barrier_progress_markdown(build_family_barrier_progress_artifact())

    assert "Family Barrier Progress" in markdown
    assert "Barrier Stage" in markdown
    assert "Prediction Readout" in markdown