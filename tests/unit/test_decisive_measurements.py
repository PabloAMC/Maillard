"""The decisive-measurements record (docs/reference/decisive_measurements.yml) and the wishlist section rendered
from it (2026-11-05).

The record is hand-kept, so what the tests guard is that it cannot silently rot: every entry has the fields a reader
needs, every CF / EX / CR id and path it cites exists (or is declared pending with where it lives), no private file is
cited, pH bases are tagged, and the wishlist renders it as its own, labelled-curated last section without touching
sections 1 - 7.
"""
from __future__ import annotations

import copy

import pytest

from src import data_wishlist, decisive_measurements as dm


def _entry(**overrides):
    entry = {
        "id": "toy_question", "plan_step": "M9", "added": "2026-11-05",
        "question": "Which toy route holds?",
        "rivals_from": ["CF-001"], "rivals": ["R1: route one", "R2: route two"],
        "why_not_on_disk": [{"blocker": "No toy pot exists.", "cites": ["CF-001", "docs/toy_case.md"]}],
        "measurement": {
            "source": "toy case sec. 4", "summary": "A fed toy pot.",
            "arms": [
                {"name": "arm 1", "design_written": True, "conditions": "toy 10 mM", "species": ["toy"],
                 "ph_values": "5, 7", "ph_basis": ["cooled"], "temperatures_C": [100], "o2": "argon vs air",
                 "replicates": "n = 3"},
                {"name": "arm 2", "design_written": False, "requirement": "design not yet written: a second toy pot"},
            ],
        },
        "decides": [{"outcome": "toy rises", "means": "R1"}],
        "unlocks": ["the toy refusal"], "existing_data_route": None, "status": "open",
    }
    entry.update(overrides)
    return entry


IDS = {"CF": {"CF-001"}, "EX": {"EX-001"}, "CR": {"CR-01"}}


def _exists(path: str) -> bool:
    return path == "docs/toy_case.md"


def _problems(payload):
    return dm.validate(payload, ids=IDS, exists=_exists)[0]


# ---------------------------------------------------------------------------- the tracked record
def test_the_tracked_record_loads_and_every_citation_resolves():
    entries, _warnings = dm.load()
    ids = [e["id"] for e in entries]
    assert len(ids) == len(set(ids)) >= 7
    # the questions the owner named on 2026-11-05 all have an entry
    for expected in ("thiol_supply_structure", "pentose_split_vs_removal", "norfuraneol_decay_ph_law", "strecker_routes",
                     "glucose_mft_carbon_origin", "thiol_absolute_second_lab", "methylglyoxal_self_loss"):
        assert expected in ids
    assert all(e["status"] in dm.STATUSES for e in entries)


def test_m3_is_provisional_until_its_prereg_freezes_and_the_rest_are_open():
    entries = {e["id"]: e for e in dm.load()[0]}
    assert entries["strecker_routes"]["status"] == "provisional"
    assert "power check" in entries["strecker_routes"]["status_note"]
    assert {e["status"] for k, e in entries.items() if k != "strecker_routes"} == {"open"}


def test_the_registers_are_read_by_prefix():
    ids = dm.register_ids()
    assert "CF-107" in ids["CF"] and "EX-059" in ids["EX"] and "CR-19" in ids["CR"]


# ---------------------------------------------------------------------------- the schema check
def test_a_well_formed_entry_passes():
    assert _problems({"entries": [_entry()]}) == []


@pytest.mark.parametrize("field", ["question", "rivals", "why_not_on_disk", "measurement", "decides", "unlocks",
                                   "existing_data_route", "status", "added", "plan_step"])
def test_a_missing_required_field_is_refused(field):
    entry = _entry()
    del entry[field]
    assert any(field in p for p in _problems({"entries": [entry]}))


def test_an_unknown_register_id_is_refused_unless_declared_pending():
    entry = _entry(why_not_on_disk=[{"blocker": "b", "cites": ["CF-999"]}])
    assert any("CF-999" in p for p in _problems({"entries": [entry]}))
    pending = {"entries": [entry], "pending_records": {"CF-999": "PR #0, unmerged"}}
    problems, warnings = dm.validate(pending, ids=IDS, exists=_exists)
    assert problems == [] and warnings == []


def test_a_pending_mark_whose_record_has_landed_is_a_warning_not_a_failure():
    payload = {"entries": [_entry()], "pending_records": {"CF-001": "PR #0, unmerged"}}
    problems, warnings = dm.validate(payload, ids=IDS, exists=_exists)
    assert problems == [] and any("CF-001" in w and "now on this branch" in w for w in warnings)


def test_an_unknown_path_is_refused():
    entry = _entry(rivals_from=["docs/no_such_case.md"])
    assert any("docs/no_such_case.md" in p for p in _problems({"entries": [entry]}))


def test_a_private_article_is_never_citable():
    private = dm._PRIVATE_PREFIX + "somepaper.pdf"
    entry = _entry(why_not_on_disk=[{"blocker": "b", "cites": [private]}])
    payload = {"entries": [entry], "pending_records": {private: "laptop"}}
    assert any("private" in p for p in _problems(payload))


def test_a_ph_basis_must_be_tagged():
    entry = _entry()
    entry["measurement"]["arms"][0]["ph_basis"] = ["room temperature"]
    assert any("pH basis" in p for p in _problems({"entries": [entry]}))


def test_an_unwritten_design_must_say_so_and_a_written_one_must_be_complete():
    entry = _entry()
    entry["measurement"]["arms"][1]["requirement"] = "a second toy pot"
    assert any("design not yet written" in p for p in _problems({"entries": [entry]}))
    entry = _entry()
    del entry["measurement"]["arms"][0]["o2"]
    assert any("o2" in p for p in _problems({"entries": [entry]}))


def test_provisional_and_decided_entries_carry_their_reason():
    assert any("status_note" in p for p in _problems({"entries": [_entry(status="provisional")]}))
    assert any("decided_by" in p for p in _problems({"entries": [_entry(status="decided")]}))
    assert any("status" in p for p in _problems({"entries": [_entry(status="maybe")]}))


def test_duplicate_ids_and_single_rival_entries_are_refused():
    assert any("duplicate" in p for p in _problems({"entries": [_entry(), _entry()]}))
    assert any("rivals" in p for p in _problems({"entries": [_entry(rivals=["R1 only"])]}))


def test_load_raises_with_every_problem_named(tmp_path):
    import yaml

    path = tmp_path / "record.yml"
    path.write_text(yaml.safe_dump({"entries": [_entry(status="maybe", rivals_from=["CF-999"])]}), encoding="utf-8")
    with pytest.raises(dm.RecordError) as err:
        dm.load(path, ids=IDS, exists=_exists)
    assert "maybe" in str(err.value) and "CF-999" in str(err.value)


# ---------------------------------------------------------------------------- the wishlist section
def _payload_with(rows):
    payload = copy.deepcopy(data_wishlist.build())
    payload["structural_questions"] = rows
    return payload


def test_the_section_renders_one_row_per_entry_and_says_it_is_curated():
    rows = dm.summarise([_entry(), _entry(id="toy_two", question="Is the toy | pipe escaped?")])
    text = "\n".join(data_wishlist._render_structural_questions({"structural_questions": rows}))
    assert data_wishlist.STRUCTURAL_SECTION_TITLE in text and "Curated, not derived" in text
    table = [line for line in text.splitlines() if line.startswith("| ") and "---" not in line]
    assert len(table) == 1 + 2   # the header and one row per entry
    assert "toy / pipe escaped" in text and "design not yet written: a second toy pot" in text
    assert "| M9 |" in text and "| open |" in text


def test_the_section_is_last_and_leaves_sections_1_to_7_unchanged():
    payload = data_wishlist.build()
    with_section = data_wishlist.render_markdown(payload)
    without = data_wishlist.render_markdown(_payload_with([]))
    assert with_section.startswith(without.rstrip("\n"))
    assert data_wishlist.STRUCTURAL_SECTION_TITLE not in without
    head, tail = with_section.split(data_wishlist.STRUCTURAL_SECTION_TITLE)
    assert "## 7." in head and "## " not in tail.replace("## 8.", "")


def test_the_built_wishlist_carries_every_recorded_question():
    payload = data_wishlist.build()
    entries = dm.load()[0]
    assert [q["id"] for q in payload["structural_questions"]] == [e["id"] for e in entries]
    assert payload["summary"]["structural_questions"] == len(entries)
    assert payload["summary"]["structural_questions_open"] == sum(e["status"] == "open" for e in entries)
    text = data_wishlist.render_markdown(payload)
    for q in payload["structural_questions"]:
        assert f"`{q['id']}`" in text and q["plan_step"] in text
    assert "decisive_measurements.yml" in str(payload["provenance"]["inputs"])
