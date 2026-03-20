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
    assert by_family["heme_nucleotide_vitamin_catalysis"]["current_status"] == "open_gap"


def test_chemistry_family_scope_markdown_mentions_ingestion_policy():
    markdown = render_chemistry_family_scope_markdown(build_chemistry_family_scope_artifact())

    assert "Chemistry Family Scope" in markdown
    assert "lipid_oxidation_and_carbonylic_crosstalk" in markdown
    assert "Ingestion policy:" in markdown