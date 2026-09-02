import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chemistry_family_scope import (  # noqa: E402
    build_chemistry_family_scope_artifact,
    render_chemistry_family_scope_markdown,
)


def test_chemistry_family_scope_recommends_lipid_crosstalk_as_next_family():
    payload = build_chemistry_family_scope_artifact()
    by_family = {row["family_id"]: row for row in payload["families"]}

    assert payload["summary"]["recommended_next_family"] == "lipid_oxidation_and_carbonylic_crosstalk"
    assert by_family["amino_acid_sugar_core"]["current_status"] == "first_class_core"
    assert by_family["lipid_oxidation_and_carbonylic_crosstalk"]["current_status"] == "partially_encoded_high_priority"
    assert by_family["nucleotide_and_ribose_support"]["current_status"] == "bounded_lane"


def test_chemistry_family_scope_uses_canonical_numbered_family_ids():
    payload = build_chemistry_family_scope_artifact()
    family_ids = {row["family_id"] for row in payload["families"]}

    assert "carbonyl_donor_hierarchy" in family_ids
    assert "fermentation_pretreatment" in family_ids
    assert "off_note_and_maillard_suppression" in family_ids


