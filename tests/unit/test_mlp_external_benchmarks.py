import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mlp_external_benchmarks import (  # noqa: E402
    build_external_mlp_landscape_payload,
    build_external_mlp_evidence_index,
    render_external_mlp_landscape_markdown,
)


def test_external_mlp_landscape_tracks_research_priority_without_granting_adoption():
    payload = build_external_mlp_landscape_payload()
    by_candidate = {row["candidate_id"]: row for row in payload["candidates"]}

    assert by_candidate["aimnet2_shortlist"]["prior_strength"] == "strong"
    assert by_candidate["aimnet2_shortlist"]["selection_priority"] == "high"
    assert "ts_initialization" in by_candidate["aimnet2_shortlist"]["supported_roles"]
    assert by_candidate["materials_foundation_example"]["domain_relevance"] == "low"


def test_external_mlp_evidence_index_and_rendering_cover_new_sota_candidates():
    index = build_external_mlp_evidence_index()
    markdown = render_external_mlp_landscape_markdown(build_external_mlp_landscape_payload())

    assert "mace_omol_external_2026q1" in index
    assert index["orbmol_external_2026q1"].prior_strength == "strong"
    assert "External ML Accelerator Landscape" in markdown
    assert "mace_omol_shortlist" in markdown