import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_ingestion_plan import (  # noqa: E402
    build_family_ingestion_plan_artifact,
    render_family_ingestion_plan_markdown,
)


def test_family_ingestion_plan_prioritizes_first_wave_for_product_value():
    payload = build_family_ingestion_plan_artifact()
    by_family = {row["slr_family"]: row for row in payload["families"]}

    assert payload["summary"]["recommended_first_wave"][:4] == ["02", "07", "10", "08"]
    assert {"11", "12", "13"}.issubset(set(payload["summary"]["recommended_first_wave"]))
    assert by_family["02"]["runtime_concept"] == "lipid_crosstalk_lane"
    assert by_family["07"]["runtime_concept"] == "carbonyl_donor_hierarchy"
    assert by_family["10"]["runtime_concept"] == "fermentation_pretreatment_node"
    assert by_family["08"]["strategic_posture"] == "guardrail_lane"


def test_family_ingestion_plan_markdown_mentions_machine_readable_policy():
    markdown = render_family_ingestion_plan_markdown(build_family_ingestion_plan_artifact())

    assert "Family Ingestion Plan" in markdown
    assert "Recommended first wave: 02, 07, 10, 08" in markdown
    assert "Deep Research Priority Surface" in markdown
    assert "Policy: extend_the_quantitative_core_by_explicit_family_lanes_with_machine_readable_payloads_not_narrative_only_docs" in markdown
    assert "Identifier policy: scope_family_id_uses_the_same_canonical_chemistry_family_id_as_payloads_and_validation_outputs" in markdown


def test_family_ingestion_plan_surfaces_deep_research_priority_queue():
    payload = build_family_ingestion_plan_artifact()
    summary = payload["summary"]
    by_family = {row["slr_family"]: row for row in payload["families"]}

    assert summary["backlog_family_count"] >= 5
    assert summary["backlog_citation_count"] >= 20
    assert summary["recommended_runtime_queue"][:3] == ["02", "07", "11"]
    assert summary["recommended_next_family"]["slr_family"] == "02"
    assert "adverse lipid markers" in summary["recommended_next_family"]["recommended_slice"]["focus"]
    assert by_family["11"]["deep_research_backlog"]["high_confidence_count"] >= 4
    assert by_family["03"]["deep_research_backlog"]["high_confidence_count"] >= 3