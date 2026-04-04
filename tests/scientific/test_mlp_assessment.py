import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chemistry_benchmark_validator import build_mlp_assessment_artifact  # noqa: E402
from src.reaction_benchmark import build_reaction_benchmark_artifact  # noqa: E402
from src.reporting import _build_scientific_surface  # noqa: E402


pytestmark = [pytest.mark.slow]


def test_mlp_artifacts_and_registry_are_exposed_in_reporting_surface():
    surface = _build_scientific_surface(ROOT)

    assert surface["reaction_benchmark_set"] == "data/lit/reaction_benchmark_set.json"
    assert surface["mlp_candidate_registry"] == "data/lit/mlp_candidate_registry.json"
    assert surface["mlp_external_benchmark_evidence"] == "data/lit/mlp_external_benchmark_evidence.json"
    assert surface["mlp_geometry_benchmark_set"] == "data/lit/geometry_benchmark_set.json"
    assert surface["mlp_ts_seed_benchmark_set"] == "data/lit/ts_seed_benchmark_set.json"


def test_mlp_assessment_quarantines_current_barrier_surrogate_and_keeps_no_default_policy():
    benchmark_payload = build_reaction_benchmark_artifact()
    assessment_payload = build_mlp_assessment_artifact()

    assert benchmark_payload["summary"]["entry_count"] >= 6
    by_candidate = {row["candidate_id"]: row for row in assessment_payload["candidates"]}
    assert by_candidate["mace_off_medium"]["decision"] == "quarantine"
    assert "geometry_backend_unavailable" in by_candidate["mace_omol_shortlist"]["stop_reasons"]
    assert "ts_seed_backend_unavailable" in by_candidate["aimnet2_shortlist"]["stop_reasons"]
    assert by_candidate["aimnet2_shortlist"]["external_prior_strength"] == "strong"
    assert by_candidate["orbmol_shortlist"]["external_selection_priority"] == "high"
    assert assessment_payload["summary"]["default_policy"] == "no_default_mlp_until_reaction_benchmark_passes"