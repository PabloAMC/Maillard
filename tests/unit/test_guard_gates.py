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
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, GOOD, "rows:\n  - {ledger_id: a}\n  - {ledger_id: b}\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, _ = verification_gate.check()
    assert problems == ["rows.yml: row of 'b' enters a fit but its source is not verified"]


def test_a_failed_check_or_a_single_reader_is_refused(tmp_path, monkeypatch):
    bad = GOOD.replace("values: pass", "values: fail").replace("agent-2", "agent-1")
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, bad, "rows:\n  - {ledger_id: a}\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, _ = verification_gate.check()
    assert any("second reader must differ" in p for p in problems)
    assert any("not verified" in p for p in problems)


def test_the_debt_is_reported_not_failed(tmp_path, monkeypatch):
    ledgers, verification, fit_rows = _write(tmp_path, LEDGER, "entries: []\n")
    monkeypatch.setattr(verification_gate, "LEDGERS", ledgers)
    monkeypatch.setattr(verification_gate, "VERIFICATION", verification)
    monkeypatch.setattr(verification_gate, "FIT_ROWS", fit_rows)
    problems, notes = verification_gate.check()
    assert problems == [] and "1 in_model_fit entries of earlier waves unverified" in notes[0]
