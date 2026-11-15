"""data/perception/thresholds.yml (perception step 7.1, docs/reference/perception/SCOPING_2026-11.md sec. 4).

Pinned:

1. the file validates against data/schemas/perception_thresholds.schema.json;
2. every odorant of the lever table (results/validation/lever_table.json, the compounds that enter a cell) appears
   exactly once, as a cited threshold or an explicit no_measured_threshold, and nothing else appears;
3. every record has medium, units, criterion, route and source; its paper is in data/keys/papers.yml; its
   compound_id is the lever spelling's registry entry, and a printed name that resolves elsewhere or nowhere
   carries a name_mapping;
4. no record cites the uncited target files, the Zhou 2023 SI or Thong 2024 Table 4;
5. one criterion per odorant per medium, and Czerny 2008 wins wherever it has a value.
6. the Vega 1994 gelatin records equal the ladder src/kinetic_core/matrix_oav.py already carries.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict

import jsonschema
import pytest

from src import compound_keys, data_access, data_paths, paper_keys

REQUIRED = ("medium", "units", "criterion", "route", "source")
FORBIDDEN_SOURCES = ("desirable_targets", "off_flavour_targets", "pathways.md", "zhou", "thong")


@pytest.fixture(scope="module")
def payload():
    return data_access.load_perception_thresholds()


def _lever_odorants():
    table = data_access.load_json(data_paths.VALIDATION_DIR / "lever_table.json")
    return {cell["compound"] for cell in table["cells"]}


def _records(payload):
    for entry in payload["odorants"]:
        for rec in entry.get("thresholds", []):
            yield entry, rec


def test_validates_against_its_schema(payload):
    schema = json.loads(data_paths.PERCEPTION_THRESHOLDS_SCHEMA.read_text(encoding="utf-8"))
    errors = [f"{'/'.join(map(str, e.path))}: {e.message}" for e in jsonschema.Draft202012Validator(schema).iter_errors(payload)]
    assert not errors, "\n".join(errors[:20])


def test_every_lever_odorant_appears_exactly_once(payload):
    counts = Counter(entry["compound"] for entry in payload["odorants"])
    assert not [c for c, n in counts.items() if n > 1], "odorants listed twice"
    lever = _lever_odorants()
    assert set(counts) == lever, (
        f"missing: {sorted(lever - set(counts))}; not in the lever table: {sorted(set(counts) - lever)}"
    )


def test_status_matches_content(payload):
    for entry in payload["odorants"]:
        if entry["status"] == "cited_threshold":
            assert entry.get("thresholds"), entry["compound"]
        else:
            assert entry["status"] == "no_measured_threshold" and entry.get("reason"), entry["compound"]


def test_required_keys_and_identities(payload):
    papers = {p.paper_id for p in paper_keys.all_papers()}
    for entry in payload["odorants"]:
        key = compound_keys.resolve(entry["compound"])
        assert key is not None and key.id == entry["compound_id"], entry["compound"]
    ids = Counter()
    for entry, rec in _records(payload):
        ids[rec["id"]] += 1
        for k in REQUIRED:
            assert rec.get(k), f"{rec['id']}: {k} missing"
        assert rec["source_paper_id"] in papers, f"{rec['id']}: {rec['source_paper_id']} not in data/keys/papers.yml"
        printed = compound_keys.resolve(rec["compound_as_printed"])
        if printed is None or printed.id != entry["compound_id"]:
            assert rec.get("name_mapping"), f"{rec['id']}: printed name does not resolve to {entry['compound_id']}"
        if rec["source_kind"] == "secondary_quotation":
            assert rec.get("primary_citation"), rec["id"]
        if rec["units"] == "ug/L of water":
            assert rec["value_ug_per_l"] == rec["value"], rec["id"]
    assert not [i for i, n in ids.items() if n > 1], "duplicate record ids"


def test_no_record_uses_an_uncited_source(payload):
    for _entry, rec in _records(payload):
        text = " ".join(str(rec.get(k, "")) for k in ("source", "source_paper_id", "primary_citation")).lower()
        assert not [s for s in FORBIDDEN_SOURCES if s in text], f"{rec['id']} cites an excluded source: {text}"


def test_one_criterion_per_odorant_per_medium_and_czerny_first(payload):
    czerny_compounds = {entry["compound"] for entry, rec in _records(payload)
                        if rec["source_paper_id"] == "czerny2008" and rec["source_kind"] == "primary_measurement"}
    for entry in payload["odorants"]:
        by_medium = defaultdict(set)
        for rec in entry.get("thresholds", []):
            by_medium[rec["medium"]].add(rec["criterion"])
            if entry["compound"] in czerny_compounds and rec["medium"] == "water":
                assert rec["source_paper_id"] == "czerny2008", (
                    f"{rec['id']}: a non-Czerny water value for an odorant Czerny 2008 measured")
        for medium, criteria in by_medium.items():
            assert len(criteria) == 1, f"{entry['compound']} in {medium}: {sorted(criteria)}"


def test_loader_raises_when_the_file_is_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(data_paths, "PERCEPTION_THRESHOLDS", tmp_path / "absent.yml")
    with pytest.raises(data_access.DataFileMissing):
        data_access.load_perception_thresholds()


def test_gelatin_records_reproduce_the_ladder_in_matrix_oav(payload):
    """The Vega 1994 gelatin values are keyed with their existing source; they must equal the engine's copy."""
    from src.kinetic_core.matrix_oav import MATRIX_THRESHOLDS

    engine_key = {"pentanal": "pentanal", "hexanal": "hexanal", "heptanal": "heptanal",
                  "e_2_octenal": "t_2_octenal", "e_e_2_4_decadienal": "tt_2_4_decadienal"}
    ladder = {(r.compound, r.temperature_c): r.value_ug_per_l for r in MATRIX_THRESHOLDS if r.matrix == "gelatin_3pct"}
    keyed = [(entry, rec) for entry, rec in _records(payload) if rec["medium"] == "gelatin_3pct"]
    assert len(keyed) == 20
    for entry, rec in keyed:
        assert ladder[(engine_key[entry["compound_id"]], rec["temperature_c"])] == rec["value_ug_per_l"], rec["id"]
