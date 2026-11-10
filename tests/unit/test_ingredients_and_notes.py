"""Ingredients in % w/w -> precursors (src/ingredients.py) and note headings (src/notes.py); sec. 0 step 3."""
import pytest

from src import ingredients, notes
from src.comparative_cli import SpecError, validate_spec

BASE = {"temp_C": 150, "time_min": 5, "ph": 5.5, "aw": 0.95}


def test_every_mapped_value_names_its_source():
    for name, entry in ingredients.table().items():
        for pre, row in (entry.get("precursors") or {}).items():
            assert row.get("source"), (name, pre)
            lo, hi = row["band"]
            assert lo <= row["central"] <= hi, (name, pre)
        if entry.get("protein_fraction"):
            assert entry["protein_fraction"].get("source") and entry.get("matrix")


def test_yeast_extract_maps_at_one_kg_per_litre_and_names_what_it_misses():
    spec, report = ingredients.expand({"ingredients": {"yeast_extract": 1.0}, "precursors": {"D-Ribose": 10.0}})
    cys = ingredients.table()["yeast_extract"]["precursors"]["L-Cysteine"]["central"]
    assert spec["precursors"]["L-Cysteine"] == pytest.approx(10.0 * cys)       # 1 % = 10 g/L
    assert spec["precursors"]["D-Ribose"] == 10.0                               # the spec's own charge is kept
    assert "ingredients" not in spec
    assert any("ribose" in k for k in report["ingredients"]["yeast_extract"]["not_mapped"])


def test_isolate_sets_matrix_and_protein():
    spec, report = ingredients.expand({"ingredients": {"pea isolate": 10}, "precursors": {"D-Ribose": 5}})
    assert spec["matrix"] == "pea_isolate" and spec["protein_g_per_l"] == pytest.approx(84.0)


def test_unknown_ingredient_and_two_isolates_refuse():
    with pytest.raises(ingredients.IngredientError):
        ingredients.expand({"ingredients": {"onion powder": 1}})
    with pytest.raises(ingredients.IngredientError):
        ingredients.expand({"ingredients": {"pea_isolate": 5, "soy_isolate": 5}, "precursors": {"D-Ribose": 1}})


def test_validate_spec_accepts_ingredients_without_precursors_and_reports():
    spec = validate_spec({**BASE, "ingredients": {"yeast_extract": 0.5}}, label="t")
    assert spec["precursors"]["L-Cysteine"] > 0 and spec["ingredient_report"]
    with pytest.raises(SpecError):
        validate_spec({**BASE, "ingredients": {"onion powder": 1}}, label="t")


def test_note_headings():
    assert notes.note_of("MFT").startswith("meaty")
    assert notes.note_of("Hexanal").startswith("off-notes")
    assert notes.note_of("2-acetylthiazole").startswith("roasted")
    assert notes.note_of("something unknown") == "other"
