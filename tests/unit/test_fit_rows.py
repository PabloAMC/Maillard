"""src/kinetic_core/fit_rows.py (guard 1, 2026-10-29): fit data come from keyed tables, never from code."""
from __future__ import annotations

import math

import pytest

from src import data_paths
from src.kinetic_core import fit_rows, fitting

LEDGER = """entries:
  - {id: paper_a_fig1, lab: "Lab A", quant_class: calibrated, flags: []}
  - {id: paper_b_tab2, lab: "Lab B", quant_class: semi_quantitative, flags: []}
  - {id: paper_c_ph, lab: "Lab C", quant_class: other, flags: []}
  - {id: paper_b_reprint, lab: "Lab B", quant_class: semi_quantitative, flags: ["duplicate_of: paper_b_tab2"]}
"""
VERIFIED = """entries:
""" + "".join(f"""  - ledger_id: {k}
    source_file: {k}.pdf
    read_by: agent-1
    second_reader: agent-2
    date: 2026-10-29
    quote: "heated at 90 C"
    checks: {{conditions: pass, units: pass, quantification: pass, duplicates: pass, values: pass}}
""" for k in ("paper_a_fig1", "paper_b_tab2", "paper_c_ph", "paper_b_reprint"))
POT = "{{ledger_id: {lid}, t_c: 90, minutes: 60, initial_mmol_L: {{PENT: 400, Cys: 400}}, condition_quote: \"x\"}}"


def _setup(tmp_path, monkeypatch, table: str, verified: str = VERIFIED):
    (tmp_path / "labs.yml").write_text("groups:\n" + "".join(
        f"  lab_{x.lower()}: {{strings: [\"Lab {x}\"], evidence: test}}\n" for x in "ABC"))
    monkeypatch.setattr(fit_rows, "LABORATORIES", tmp_path / "labs.yml")
    monkeypatch.setattr(data_paths, "SOURCE_VERIFICATION_DIR", tmp_path / "no_dir")
    ledgers = tmp_path / "ledger"
    ledgers.mkdir()
    (ledgers / "sulfur.yml").write_text(LEDGER)
    (tmp_path / "verification.yml").write_text(verified)
    rows_dir = tmp_path / "fit_rows"
    rows_dir.mkdir()
    (rows_dir / "t.yml").write_text(table)
    monkeypatch.setattr(data_paths, "MEASUREMENT_LEDGER_DIR", ledgers)
    monkeypatch.setattr(data_paths, "SOURCE_VERIFICATION", tmp_path / "verification.yml")
    monkeypatch.setattr(fit_rows, "FIT_ROWS_DIR", rows_dir)


def _table(rows: str, pots: str = "") -> str:
    return ("status: ready\npots:\n"
            f"  a: {POT.format(lid='paper_a_fig1')}\n"
            f"  b: {POT.format(lid='paper_b_tab2')}\n"
            f"  c: {POT.format(lid='paper_c_ph')}\n" + pots + "rows:\n" + rows)


def test_laboratory_and_method_come_from_the_ledger(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table(
        "  - {pot: a, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, sd: 1.0}\n"
        "  - {pot: b, species: MFT, t_min: 20, value: 5.0, unit: ug_L, read_error: 0.1}\n"
        "  - {pot: c, species: pH, t_min: 20, value: 5.5, unit: pH, sd: 0.05}\n"))
    _, rows = fit_rows.load("t")
    by_pot = {r.pot: r for r in rows}
    assert (by_pot["a"].laboratory, by_pot["a"].method) == ("lab_a", "calibrated")
    assert (by_pot["b"].laboratory, by_pot["b"].method) == ("lab_b", "semi_quantitative")
    assert by_pot["c"].method == "calibrated" and not by_pot["c"].log_scale
    assert by_pot["a"].sigma == pytest.approx(math.log10(math.e) * 0.1)      # 10 % SD as dex
    assert by_pot["c"].sigma == pytest.approx(0.05)                            # pH units
    assert fit_rows.measurement_model().nuisance_dof(rows) == 1                # the semi-quantitative series


def test_an_unverified_source_is_refused(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table("  - {pot: a, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, sd: 1}\n"),
           verified="entries: []\n")
    with pytest.raises(PermissionError, match="not verified"):
        fit_rows.load("t")


def test_a_row_cannot_contradict_the_ledgers_quantification_class(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table(
        "  - {pot: b, species: MFT, t_min: 20, value: 5.0, unit: ug_L, sd: 1, method: calibrated}\n"))
    with pytest.raises(ValueError, match="the ledger says 'semi_quantitative'"):
        fit_rows.load("t")


def test_a_row_without_an_error_is_refused(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table("  - {pot: a, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L}\n"))
    with pytest.raises(ValueError, match="no sd and no read_error"):
        fit_rows.load("t")


def test_an_unknown_class_needs_a_declared_method_and_its_basis(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table("  - {pot: c, species: HMF, t_min: 20, value: 2.0, unit: mmol_L, sd: 0.2}\n"))
    with pytest.raises(ValueError, match="method_basis"):
        fit_rows.load("t")


def test_a_reprint_is_dropped_when_its_original_is_present(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table(
        "  - {pot: b, species: MFT, t_min: 20, value: 5.0, unit: ug_L, read_error: 0.1}\n"
        "  - {pot: r, species: MFT, t_min: 20, value: 5.0, unit: ug_L, read_error: 0.1}\n",
        pots=f"  r: {POT.format(lid='paper_b_reprint')}\n"))
    _, rows = fit_rows.load("t")
    assert [r.pot for r in rows] == ["b"]


def test_residuals_are_log_ratios_for_amounts_and_differences_for_ph():
    amount = fit_rows.FitRow("t", "a", "x", "L", "calibrated", "TTCA", 20.0, 10.0, "mmol_L", 0.04, "s")
    ph = fit_rows.FitRow("t", "c", "x", "L", "calibrated", "pH", 20.0, 5.5, "pH", 0.05, "s2")
    assert fit_rows.residual(amount, 100.0) == pytest.approx(1.0)
    assert fit_rows.residual(ph, 5.3) == pytest.approx(-0.2)
    assert isinstance(fit_rows.measurement_model(), fitting.MeasurementModel)


def test_a_draft_table_is_refused_for_fitting(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table("  - {pot: a, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, sd: 1}\n")
           .replace("status: ready", "status: draft"))
    with pytest.raises(PermissionError, match="not 'ready'"):
        fit_rows.load("t")
    assert len(fit_rows.load("t", require_verified=False)[1]) == 1


def test_an_unmapped_laboratory_is_refused(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, _table("  - {pot: a, species: TTCA, t_min: 20, value: 10.0, unit: mmol_L, sd: 1}\n"))
    (tmp_path / "labs.yml").write_text("groups: {}\n")
    with pytest.raises(KeyError, match="laboratories.yml"):
        fit_rows.load("t")
