"""scripts/reconcile_source_readings.py: two independent readings agree, or the difference is shown (2026-10-29)."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import reconcile_source_readings as rec  # noqa: E402
from src import data_paths  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402

CHECKS = "{conditions: pass, units: pass, quantification: pass, duplicates: pass, values: pass}"


def _verification(reader, ids, second="pending"):
    return "entries:\n" + "".join(
        f"  - ledger_id: {i}\n    source_file: {i}.pdf\n    read_by: {reader}\n    second_reader: {second}\n"
        f"    date: 2026-10-29\n    quote: \"q\"\n    checks: {CHECKS}\n" for i in ids)


def _setup(tmp_path, monkeypatch, first_rows, second_rows, first_ids=("a", "b"), second_ids=("a", "b", "c")):
    vdir, rdir = tmp_path / "v", tmp_path / "r"
    (vdir / "second").mkdir(parents=True)
    (rdir / "second").mkdir(parents=True)
    (vdir / "t.yml").write_text(_verification("reader1", first_ids))
    (vdir / "second" / "t.yml").write_text(_verification("reader2", second_ids, second="self"))
    pots1 = "pots:\n  p1: {ledger_id: a, t_c: 90, minutes: 60, label_ph: 7.0, initial_mmol_L: {PENT: 400, Cys: 400}}\n" \
            "  p2: {ledger_id: b, t_c: 120, minutes: 30, label_ph: 6.5, initial_mmol_L: {Glc: 200}}\n"
    pots2 = "pots:\n  x: {ledger_id: a, t_c: 90.0, minutes: 60, label_ph: 7, initial_mmol_L: {Cys: 400.0, PENT: 400}}\n" \
            "  y: {ledger_id: b, t_c: 120, minutes: 30, label_ph: 6.5, initial_mmol_L: {Glc: 200}}\n"
    (rdir / "t.yml").write_text("status: draft\n" + pots1 + "rows:\n" + first_rows)
    (rdir / "second" / "t.yml").write_text("status: draft\n" + pots2 + "rows:\n" + second_rows)
    monkeypatch.setattr(data_paths, "SOURCE_VERIFICATION_DIR", vdir)
    monkeypatch.setattr(fit_rows, "FIT_ROWS_DIR", rdir)
    monkeypatch.setattr(data_paths, "VALIDATION_DIR", tmp_path)
    return vdir


def test_matching_readings_agree_despite_different_pot_names(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch,
           "  - {pot: p1, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, read_error: 0.05}\n"
           "  - {pot: p2, species: pH, t_min: 30, value: 5.50, unit: pH, read_error: 0.03}\n",
           "  - {pot: x, species: TTCA, t_min: 20, value: 10.5, unit: mmol_L, read_error: 0.05}\n"
           "  - {pot: y, species: pH, t_min: 30, value: 5.48, unit: pH, read_error: 0.03}\n")
    r = rec.reconcile("t")
    assert r["agree"] == ["a", "b"] and r["disagree"] == {}
    assert r["only_second"] == ["c"] and r["only_first"] == []


def test_a_value_beyond_both_errors_is_a_disagreement(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch,
           "  - {pot: p1, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, read_error: 0.05}\n",
           "  - {pot: x, species: TTCA, t_min: 20, value: 14.0, unit: mmol_L, read_error: 0.05}\n")
    r = rec.reconcile("t")
    assert "a" in r["disagree"] and any("10.0 vs 14.0" in p for p in r["disagree"]["a"])


def test_a_pot_with_different_conditions_does_not_match(tmp_path, monkeypatch):
    vdir = _setup(tmp_path, monkeypatch, "", "")
    (fit_rows.FIT_ROWS_DIR / "second" / "t.yml").write_text(
        "status: draft\npots:\n  x: {ledger_id: a, t_c: 90, minutes: 60, label_ph: 5.5, initial_mmol_L: {PENT: 400, Cys: 400}}\n"
        "  y: {ledger_id: b, t_c: 120, minutes: 30, label_ph: 6.5, initial_mmol_L: {Glc: 200}}\nrows: []\n")
    r = rec.reconcile("t")
    assert "a" in r["disagree"] and "b" in r["agree"]
    assert vdir.exists()


def test_promotion_names_the_second_reader_on_agreed_entries_only(tmp_path, monkeypatch):
    vdir = _setup(tmp_path, monkeypatch,
                  "  - {pot: p1, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, read_error: 0.05}\n",
                  "  - {pot: x, species: TTCA, t_min: 20, value: 14.0, unit: mmol_L, read_error: 0.05}\n")
    r = rec.reconcile("t")
    second = rec._entries(vdir / "second" / "t.yml")
    assert rec.promote("t", r["agree"], second) == 1
    entries = {e["ledger_id"]: e for e in yaml.safe_load((vdir / "t.yml").read_text())["entries"]}
    assert entries["b"]["second_reader"] == "reader2"
    assert entries["a"]["second_reader"] == "pending"
