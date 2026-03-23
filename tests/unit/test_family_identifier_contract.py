import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.family_identifier_contract import (  # noqa: E402
    build_family_identifier_contract_artifact,
    render_family_identifier_contract_markdown,
)


def test_family_identifier_contract_aligns_scope_plan_and_payloads():
    payload = build_family_identifier_contract_artifact()

    assert payload["summary"]["scope_registry_covers_plan"] is True
    assert payload["summary"]["payload_surfaces_use_only_canonical_families"] is True
    assert payload["summary"]["chemistry_and_matrix_axes_are_disjoint"] is True
    assert payload["summary"]["plan_families_missing_from_scope"] == []
    assert payload["summary"]["payload_unknown_family_ids"] == []


def test_family_identifier_contract_markdown_mentions_axis_policy():
    markdown = render_family_identifier_contract_markdown(build_family_identifier_contract_artifact())

    assert "Family Identifier Contract" in markdown
    assert "Chemistry/matrix axis overlap: none" in markdown
    assert "Axis policy: chemistry_family_and_matrix_family_remain_separate_axes_even_when_a_payload_references_both" in markdown