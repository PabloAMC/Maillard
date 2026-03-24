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

    assert payload["summary"]["recommended_first_wave"] == ["02", "07", "10", "08"]
    assert by_family["02"]["runtime_concept"] == "lipid_crosstalk_lane"
    assert by_family["07"]["runtime_concept"] == "carbonyl_donor_hierarchy"
    assert by_family["10"]["runtime_concept"] == "fermentation_pretreatment_node"
    assert by_family["08"]["strategic_posture"] == "guardrail_lane"


def test_family_ingestion_plan_markdown_mentions_machine_readable_policy():
    markdown = render_family_ingestion_plan_markdown(build_family_ingestion_plan_artifact())

    assert "Family Ingestion Plan" in markdown
    assert "Recommended first wave: 02, 07, 10, 08" in markdown
    assert "Policy: extend_the_quantitative_core_by_explicit_family_lanes_with_machine_readable_payloads_not_narrative_only_docs" in markdown
    assert "Identifier policy: scope_family_id_uses_the_same_canonical_chemistry_family_id_as_payloads_and_validation_outputs" in markdown