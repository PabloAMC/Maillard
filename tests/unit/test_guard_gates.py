"""The 2026-10-29 guards catch what they exist to catch."""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "ci"))

import import_time_switch_gate as switch_gate  # noqa: E402
import source_verification_gate as verification_gate  # noqa: E402


@pytest.fixture(autouse=True)
def _no_real_verification_files(tmp_path, monkeypatch):
    """The gate also reads data/lit/source_verification/*.yml; the synthetic cases below must not see real entries."""
    monkeypatch.setattr(verification_gate, "VERIFICATION_DIR", tmp_path / "no_such_dir")


@pytest.mark.parametrize("code, expected", [
    ("import x\nsulfur.PH_COUPLING = 'inline'\n", ["sulfur.PH_COUPLING"]),
    ("if True:\n    ps.ROUTE = 1\n", ["ps.ROUTE"]),
    ("try:\n    ph_state.B51 = 2\nexcept Exception:\n    pass\n", ["ph_state.B51"]),
    ("a.X, b.Y = 1, 2\n", ["a.X", "b.Y"]),
    ("def f():\n    sulfur.PH_COUPLING = 'inline'\n", []),          # inside a function: fine
    ("class C:\n    def g(self):\n        ps.ROUTE = 1\n", []),
    ("sulfur.lower_case = 1\nx = 3\n", []),                         # not a switch
])
def test_switch_gate_finds_import_time_assignments(code, expected):
    hits = switch_gate._top_level_switch_assignments(ast.parse(code))
    assert [name for _, name in hits] == expected


def test_switch_gate_passes_on_the_repository():
    assert switch_gate.check() == []


def _write(tmp_path, ledger, verification, rows=None):
    ledgers = tmp_path / "ledgers"
    ledgers.mkdir()
    (ledgers / "sulfur.yml").write_text(ledger)
    (tmp_path / "verification.yml").write_text(verification)
    fit_rows = tmp_path / "fit_rows"
    fit_rows.mkdir()
    if rows is not None:
        (fit_rows / "rows.yml").write_text(rows)
    return ledgers, tmp_path / "verification.yml", fit_rows


LEDGER = "entries:\n  - {id: a, status: in_model_fit}\n  - {id: b, status: digitised_in_dossier}\n"
GOOD = """entries:
  - ledger_id: a
    source_file: paper.pdf
    read_by: agent-1
    second_reader: agent-2
    date: 2026-10-29
    quote: "heated at 90 C in 0.2 M phosphate, pH 5.5"
    checks: {conditions: pass, units: pass, quantification: corrected, duplicates: pass, values: pass}
"""


def test_an_unverified_row_cannot_enter_a_fit(tmp_path, monkeypatch):
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, GOOD, "status: ready\npots:\n  pa: {ledger_id: a}\n  pb: {ledger_id: b}\nrows: [{pot: pa}, {pot: pb}]\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, _ = verification_gate.check()
    assert problems == ["rows.yml is ready but its source 'b' is not verified by two readers"]


def test_a_failed_check_or_a_single_reader_is_refused(tmp_path, monkeypatch):
    bad = GOOD.replace("values: pass", "values: fail").replace("agent-2", "agent-1")
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, bad, "status: ready\npots:\n  pa: {ledger_id: a}\nrows: [{pot: pa}]\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, _ = verification_gate.check()
    assert any("second reader must differ" in p for p in problems)
    assert any("not verified by two readers" in p for p in problems)


def test_the_debt_is_reported_not_failed(tmp_path, monkeypatch):
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, "entries: []\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, notes = verification_gate.check()
    assert problems == [] and "1 in_model_fit entries of earlier waves unverified" in notes[0]


def test_anchor_tokens_and_resolution():
    from src import evidence_anchors as ea

    text = "martins2005_extraction.md (Table 2); data/lit/extraction_dossiers/k3_final_parameter_inventory.md sec. A.1"
    assert ea.anchor_tokens(text) == ["martins2005_extraction.md",
                                      "data/lit/extraction_dossiers/k3_final_parameter_inventory.md"]
    assert ea.resolve("martins2005_extraction.md") is not None
    assert ea.resolve("FIT_HOLDOUT_DECLARATION.md") is not None      # a uniquely named tracked file
    assert ea.resolve("no_such_paper2099" + "_extraction.md") is None   # built here so the coverage gate does not read it as a citation
    assert ea.dossier_stems(text) == ["martins2005"]


def test_every_maillard_constant_cites_a_resolvable_source():
    sys.path.insert(0, str(ROOT / "scripts" / "ci"))
    import evidence_anchor_gate

    problems, n = evidence_anchor_gate.check()
    assert problems == [] and n > 100


def test_unpinned_constants_gate_catches_a_silent_omission_and_a_stale_entry(tmp_path, monkeypatch):
    import json as _json

    import unpinned_constants_gate as gate

    table = tmp_path / "table.json"
    register = tmp_path / "register.yml"
    table.write_text(_json.dumps({"summary": {"unpinned_constants_that_matter": ["k_a", "k_b"],
                                              "unpinned_and_rate_not_sampled": ["k_a", "k_b"]}}))
    register.write_text("entries:\n  - {key: k_a, planned_test: sub-model 2}\n  - {key: k_old, planned_test: done}\n")
    monkeypatch.setattr(gate, "TABLE", table)
    monkeypatch.setattr(gate, "REGISTER", register)
    problems = gate.check()
    assert any(p.startswith("k_b:") and "not in" in p for p in problems)
    assert any(p.startswith("k_old:") and "only shrinks" in p for p in problems)
    assert gate.check.__module__ == "unpinned_constants_gate"


def test_unpinned_constants_gate_passes_on_the_repository():
    import unpinned_constants_gate as gate

    assert gate.check() == []


def test_a_pending_second_reader_is_not_a_verification(tmp_path, monkeypatch):
    pending = GOOD.replace("second_reader: agent-2", "second_reader: pending")
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, pending, "status: ready\npots:\n  pa: {ledger_id: a}\nrows: [{pot: pa}]\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "VERIFICATION_DIR", tmp_path / "none")
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, _ = verification_gate.check()
    assert problems == ["rows.yml is ready but its source 'a' is not verified by two readers"]
