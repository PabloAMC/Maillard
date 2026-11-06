"""The S1 evidence-tier gate (scripts/ci/evidence_tier_gate.py) catches each failure mode it exists for.

Every case plants one violation in a synthetic record under tmp_path (never data/ or results/).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "ci"))
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import evidence_tier_gate as gate  # noqa: E402

from src import evidence_tiers as et, provenance  # noqa: E402

SHIPPED = ["k_a.log10_k_ref", "k_a.ea_kj_mol", "k_b.log10_k_ref", "k_b.ea_kj_mol"]
BASELINE_TIERS = {"k_a.log10_k_ref": "measured", "k_a.ea_kj_mol": "analogue",
                  "k_b.log10_k_ref": "fitted_unidentified", "k_b.ea_kj_mol": "structural"}


@pytest.fixture
def world(tmp_path, monkeypatch):
    """A synthetic record, its input and its baseline; the validation dir is tmp_path/results/validation."""
    vdir = tmp_path / "results" / "validation"
    vdir.mkdir(parents=True)
    monkeypatch.setattr(gate, "VALIDATION_DIR", vdir)
    source = tmp_path / "input.yml"
    source.write_text("x: 1\n", encoding="utf-8")
    rows = [{"id": i, "constant": i.split(".")[0], "tier": t, "basis": "test", "identification": None,
             "tier_changes": []} for i, t in BASELINE_TIERS.items()]
    record = {"provenance": {"inputs": [{"path": str(source), "sha256": provenance.sha256_of(source)}]},
              "rows": rows, "tier_changes_without_a_row": []}
    (tmp_path / "baseline.json").write_text(json.dumps({"tiers": BASELINE_TIERS}), encoding="utf-8")

    def run(mutate=None, shipped=SHIPPED):
        payload = json.loads(json.dumps(record))
        if mutate:
            mutate(payload, {r["id"]: r for r in payload["rows"]})
        (tmp_path / "record.json").write_text(json.dumps(payload), encoding="utf-8")
        problems, _ = gate.check(tmp_path / "record.json", tmp_path / "baseline.json", shipped=shipped, root=tmp_path)
        return problems

    run.tmp = tmp_path
    run.vdir = vdir
    run.source = source
    return run


def _profile(world, name="profile.json", text="k_b profiled with every unidentified coordinate over its band"):
    (world.vdir / name).write_text(text, encoding="utf-8")
    return f"results/validation/{name}"


def test_the_clean_record_passes(world):
    assert world() == []


# --- 1. a shipped coordinate without a tier
def test_a_shipped_coordinate_without_a_row_fails(world):
    problems = world(shipped=SHIPPED + ["k_c.log10_k_ref"])
    assert len(problems) == 1 and "k_c.log10_k_ref" in problems[0] and "no row" in problems[0]


@pytest.mark.parametrize("tier", ["unclassified", None, "fitted"])
def test_a_row_without_a_tier_fails(world, tier):
    problems = world(lambda p, rows: rows["k_a.ea_kj_mol"].update(tier=tier))
    assert any("k_a.ea_kj_mol: no tier" in x for x in problems)


def test_a_duplicated_row_fails(world):
    problems = world(lambda p, rows: p["rows"].append(dict(rows["k_a.ea_kj_mol"])))
    assert any("k_a.ea_kj_mol: 2 rows" in x for x in problems)


# --- 2. a stale record
def test_a_stale_input_fails(world):
    def mutate(p, rows):
        world.source.write_text("x: 2\n", encoding="utf-8")
    problems = world(mutate)
    assert len(problems) == 1 and "stale input" in problems[0]


def test_a_record_without_provenance_fails(world):
    problems = world(lambda p, rows: p.pop("provenance"))
    assert any("no provenance" in x for x in problems)


def test_a_register_entry_naming_no_row_fails(world):
    problems = world(lambda p, rows: p["tier_changes_without_a_row"].append({"coordinate": "k_gone.ea_kj_mol"}))
    assert any("k_gone.ea_kj_mol" in x for x in problems)


# --- 3. fitted_identified without the S1 profile
def _upgrade_k_b(world, profile, **entry):
    change = {"coordinate": "k_b.log10_k_ref", "from": "fitted_unidentified", "to": "fitted_identified",
              "date": "2026-11-07", "citation": "test fit report", "artifact": profile, "second_reader": "reader B"}
    change.update(entry)

    def mutate(p, rows):
        rows["k_b.log10_k_ref"].update(tier="fitted_identified", identification={"profile_artifact": profile},
                                       tier_changes=[change])
    return mutate


def test_fitted_identified_with_its_profile_and_review_passes(world):
    assert world(_upgrade_k_b(world, _profile(world))) == []


def test_fitted_identified_without_a_profile_fails(world):
    problems = world(_upgrade_k_b(world, "results/validation/no_such_profile.json"))
    assert any("without its S1 profile" in x for x in problems)


def test_a_profile_outside_results_validation_fails(world):
    (world.tmp / "elsewhere.json").write_text("k_b", encoding="utf-8")
    problems = world(_upgrade_k_b(world, "elsewhere.json"))
    assert any("without its S1 profile" in x for x in problems)


def test_a_profile_that_does_not_name_the_constant_fails(world):
    problems = world(_upgrade_k_b(world, _profile(world, text="a profile of k_other only")))
    assert any("does not name k_b" in x for x in problems)


# --- 4. an upgrade without a second reader (or without any entry)
def test_an_upgrade_without_a_second_reader_fails(world):
    problems = world(_upgrade_k_b(world, _profile(world), second_reader=""))
    assert problems == ["k_b.log10_k_ref: change fitted_unidentified -> fitted_identified: an upgrade without a second_reader"]


def test_an_upgrade_without_an_artifact_fails(world):
    profile = _profile(world)

    def mutate(p, rows):
        _upgrade_k_b(world, profile, artifact=None)(p, rows)
    problems = world(mutate)
    assert any("an upgrade without an artifact" in x for x in problems)


def test_an_upgrade_with_no_register_entry_fails(world):
    problems = world(lambda p, rows: rows["k_b.ea_kj_mol"].update(tier="measured"))
    assert problems == ["k_b.ea_kj_mol: tier structural -> measured is an upgrade with no entry in "
                        "docs/reference/evidence_tier_changes.yml"]


def test_a_downgrade_needs_only_a_citation(world):
    change = {"coordinate": "k_a.log10_k_ref", "from": "measured", "to": "conflicting_measurements",
              "date": "2026-11-07", "citation": "CF-999"}

    def ok(p, rows):
        rows["k_a.log10_k_ref"].update(tier="conflicting_measurements", tier_changes=[change])
    assert world(ok) == []

    def uncited(p, rows):
        rows["k_a.log10_k_ref"].update(tier="conflicting_measurements", tier_changes=[dict(change, citation="")])
    assert world(uncited) == ["k_a.log10_k_ref: change measured -> conflicting_measurements: no citation"]

    def unregistered(p, rows):
        rows["k_a.log10_k_ref"].update(tier="conflicting_measurements")
    assert any("is a downgrade with no entry" in x for x in world(unregistered))


def test_a_change_whose_from_is_not_the_tier_before_fails(world):
    change = {"coordinate": "k_a.ea_kj_mol", "from": "bounded", "to": "fitted_unidentified", "citation": "x"}
    problems = world(lambda p, rows: rows["k_a.ea_kj_mol"].update(tier="fitted_unidentified", tier_changes=[change]))
    assert any("is not the tier before it (analogue)" in x for x in problems)


def test_a_new_coordinate_entering_at_an_upgrade_tier_needs_review(world):
    new = {"id": "k_c.log10_k_ref", "constant": "k_c", "tier": "measured", "basis": "t", "identification": None,
           "tier_changes": []}
    problems = world(lambda p, rows: p["rows"].append(dict(new)), shipped=SHIPPED + ["k_c.log10_k_ref"])
    assert problems == ["k_c.log10_k_ref: tier None -> measured is an upgrade with no entry in "
                        "docs/reference/evidence_tier_changes.yml"]
    assert world(lambda p, rows: p["rows"].append(dict(new, tier="fitted_unidentified"))) == []


@pytest.mark.parametrize("before, after, up", [
    ("fitted_unidentified", "fitted_identified", True), ("fitted_unidentified", "measured", True),
    ("structural", "bounded", True), ("withdrawn", "analogue", True), (None, "measured", True),
    ("analogue", "measured", True), ("bounded", "fitted_identified", True),
    ("measured", "analogue", False), ("analogue", "bounded", False), ("measured", "conflicting_measurements", False),
    ("fitted_identified", "fitted_unidentified", False), ("measured", "measured", False),
    (None, "fitted_unidentified", False),
])
def test_what_counts_as_an_upgrade(before, after, up):
    assert et.is_upgrade(before, after) is up


def test_the_baseline_is_write_once(tmp_path, monkeypatch):
    import build_evidence_tiers as gen

    frozen = tmp_path / "baseline.json"
    frozen.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(gen, "BASELINE", frozen)
    monkeypatch.setattr(gen, "build", lambda: pytest.fail("must refuse before building"))
    assert gen.main(["--freeze-baseline"]) == 1
    assert frozen.read_text(encoding="utf-8") == "{}"


def test_the_gate_passes_on_the_repository():
    problems, stats = gate.check(shipped=et.shipped_coordinate_ids())
    assert problems == []
    assert stats["shipped"] > 0 and stats["by_tier"].get("fitted_identified", 0) == 0
