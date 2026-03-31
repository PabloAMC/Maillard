import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_family_coverage import (  # noqa: E402
    build_matrix_family_coverage_artifact,
    render_matrix_family_coverage_markdown,
)


def test_matrix_family_coverage_artifact_distinguishes_explicit_from_indirect_support():
    payload = build_matrix_family_coverage_artifact()
    by_family = {row["matrix_family"]: row for row in payload["families"]}

    assert by_family["pea_isolate"]["runtime_posture"] == "directional_matrix"
    assert by_family["pea_isolate"]["support_class"] == "explicit_supported"
    assert by_family["soy_isolate"]["runtime_posture"] == "directional_matrix"
    assert by_family["coconut_oil_co_matrix"]["runtime_posture"] == "indirect_generic_support"
    assert by_family["coconut_oil_co_matrix"]["expansion_status"] == "blocked_on_family_specific_evidence"
    assert by_family["coconut_oil_co_matrix"]["scope_priority"] == "scope_gap_to_rank"
    assert "coconut-specific lipid profile" in by_family["coconut_oil_co_matrix"]["what_is_not_supported"]
    assert by_family["mycoprotein"]["expansion_status"] == "bounded_expansion_candidate"
    assert by_family["mycoprotein"]["scope_priority"] == "bounded_next_candidate"
    assert "other_plant_proteins" in payload["summary"]["open_gap_families"]
    assert payload["summary"]["bounded_expansion_candidates"] == ["mycoprotein"]
    assert "pea_isolate" in payload["summary"]["active_scope_priorities"]
    assert "soy_isolate" in payload["summary"]["active_scope_priorities"]
    assert "coconut_oil_co_matrix" in payload["summary"]["scope_gap_priorities"]


def test_matrix_family_coverage_markdown_mentions_coconut_gap_policy():
    markdown = render_matrix_family_coverage_markdown(build_matrix_family_coverage_artifact())

    assert "Matrix Family Coverage" in markdown
    assert "Expansion Gates" in markdown
    assert "Scope Priority" in markdown
    assert "coconut_oil_co_matrix" in markdown
    assert "Scope-gap priorities" in markdown
    assert "Indirect-only families" in markdown