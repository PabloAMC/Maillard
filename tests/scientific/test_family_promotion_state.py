import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_promotion_state import build_family_promotion_state_artifact, render_family_promotion_state_markdown


def test_family_promotion_state_promotes_family_07_as_benchmark_linked_support():
    payload = build_family_promotion_state_artifact()
    by_id = {row["chemistry_family"]: row for row in payload["families"]}

    assert by_id["carbonyl_donor_hierarchy"]["promotion_state"] == "benchmark_linked_support"
    assert by_id["carbonyl_donor_hierarchy"]["support_benchmark_count"] >= 1
    assert by_id["carbonyl_donor_hierarchy"]["explicit_uncertainty_bounds"] is True
    assert "carbonyl_donor_hierarchy" in payload["summary"]["promoted_non_quantitative_families"]

    markdown = render_family_promotion_state_markdown(payload)
    assert "Family Promotion State" in markdown
    assert "benchmark_linked_support" in markdown
    assert "carbonyl_donor_hierarchy" in markdown