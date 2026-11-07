"""The frozen-fit drift gate (scripts/ci/frozen_fit_drift_gate.py) catches each failure it exists for.

Toy fits run through the gate's own ``check`` / ``main`` with the evaluation replaced by fixed probes and the
allow-list, register and record written under tmp_path (never data/ or results/). One test runs the gate on the
repository for one allow-listed drifted fit (B18) and one that reproduces (B20).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "ci"))

import frozen_fit_drift_gate as gate  # noqa: E402

FITS = ("T1", "T2")


def _probe(fit, s_today, s_recorded=10.0, rows_then=("a", "b"), rows_now=("a", "b")):
    return {"fit": fit, "S_recorded": s_recorded, "S_today": s_today, "rel_tol": 0.01, "n_rows_then": len(rows_then),
            "n_rows_now": len(rows_now), "rows_then": list(rows_then), "rows_now": list(rows_now)}


def _record(fits):
    """fit -> (class, S_today, n_rows_now), as the reproducibility record keeps them."""
    return {"fits": {f: {"verdict": {"class": c}, "S_today": s, "n_rows_now": n} for f, (c, s, n) in fits.items()}}


@pytest.fixture
def world(tmp_path, monkeypatch):
    """Run the gate on toy fits: returns run(probes, entries, record) -> (exit code, problems)."""
    register = tmp_path / "conflicts_register.yml"
    register.write_text(yaml.safe_dump({"entries": [{"id": "CF-900"}]}), encoding="utf-8")
    monkeypatch.setattr(gate, "REGISTER", register)
    monkeypatch.setattr(gate, "ALLOWLIST", tmp_path / "allowlist.yml")
    monkeypatch.setattr(gate, "RECORD", tmp_path / "record.json")
    monkeypatch.setattr(gate, "_generator", lambda: SimpleNamespace(FIT_ORDER=FITS, PROBE_PREFIX="X "))

    def run(probes, entries, record):
        (tmp_path / "allowlist.yml").write_text(yaml.safe_dump({"entries": entries}), encoding="utf-8")
        (tmp_path / "record.json").write_text(json.dumps(record), encoding="utf-8")
        monkeypatch.setattr(gate, "evaluate", lambda fits, jobs=2: {f: probes[f] for f in fits})
        problems, _ = gate.check()
        return gate.main([]), problems

    return run


def _entry(s, cls="MATERIAL", cf="CF-900", n=2, precision=1e-5):
    return {"fit": "T2", "class": cls, "S": s, "S_rel_precision": precision, "n_rows": n, "cf": cf}


REPRODUCING = _probe("T1", 10.05)


def test_reproducing_fits_pass(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 9.99)}, [],
                           _record({"T1": ("REPRODUCES", 10.05, 2), "T2": ("REPRODUCES", 9.99, 2)}))
    assert (code, problems) == (0, [])


def test_drifted_fit_without_an_entry_fails(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 25.0)}, [],
                           _record({"T1": ("REPRODUCES", 10.05, 2)}))
    assert code == 1
    assert any(p.startswith("T2: DRIFTED") for p in problems)


def test_row_set_change_fails_at_the_same_s(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 10.0, rows_now=("a", "c"))}, [],
                           _record({"T1": ("REPRODUCES", 10.05, 2)}))
    assert code == 1
    assert any("T2: DRIFTED" in p and "withdrawn ['b']" in p for p in problems)


def test_allow_listed_drift_passes(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 25.0)}, [_entry(25.0)],
                           _record({"T1": ("REPRODUCES", 10.05, 2), "T2": ("DRIFTED-MATERIAL", 25.0, 2)}))
    assert (code, problems) == (0, [])


def test_stale_entry_fails_when_s_moves_again(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 26.0)}, [_entry(25.0)],
                           _record({"T1": ("REPRODUCES", 10.05, 2), "T2": ("DRIFTED-MATERIAL", 25.0, 2)}))
    assert code == 1
    assert any(p.startswith("T2: STALE allow-list entry") for p in problems)


def test_entry_for_a_fit_that_reproduces_again_fails(world):
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 10.0)}, [_entry(10.0)],
                           _record({"T1": ("REPRODUCES", 10.05, 2), "T2": ("DRIFTED-MATERIAL", 10.0, 2)}))
    assert code == 1
    assert any("reproduces again" in p for p in problems)


@pytest.mark.parametrize("entry, record_t2, needle", [
    (_entry(25.0, cf="CF-901"), ("DRIFTED-MATERIAL", 25.0, 2), "not an entry of"),
    (_entry(25.0, cls="BENIGN"), ("DRIFTED-MATERIAL", 25.0, 2), "classes it DRIFTED-MATERIAL"),
    (_entry(25.0, cls="MINOR"), ("DRIFTED-MATERIAL", 25.0, 2), "class 'MINOR'"),
    (_entry(25.0, precision=0.05), ("DRIFTED-MATERIAL", 25.0, 2), "S_rel_precision"),
    (_entry(25.0), ("DRIFTED-MATERIAL", 24.0, 2), "not judged at this S"),
    (_entry(25.0), None, "has no audit of T2"),
])
def test_malformed_or_unjudged_entries_fail(world, entry, record_t2, needle):
    rec = {"T1": ("REPRODUCES", 10.05, 2)}
    if record_t2 is not None:
        rec["T2"] = record_t2
    code, problems = world({"T1": REPRODUCING, "T2": _probe("T2", 25.0)}, [entry], _record(rec))
    assert code == 1
    assert any(needle in p for p in problems), problems


def test_record_drift_without_an_entry_fails_even_when_not_evaluated(world, monkeypatch):
    monkeypatch.setattr(gate, "SLOW_FITS", ("T2",))  # T2 is skipped without --slow
    code, problems = world({"T1": REPRODUCING}, [],
                           _record({"T1": ("REPRODUCES", 10.05, 2), "T2": ("DRIFTED-MATERIAL", 25.0, 2)}))
    assert code == 1
    assert any("no allow-list entry names it" in p for p in problems)


def test_gate_passes_on_the_repository_for_a_drifted_and_a_reproducing_fit():
    problems, info = gate.check(only=["B18", "B20"])
    assert problems == []
    assert set(info["probes"]) == {"B18", "B20"}
    assert info["allowed"] == ["B18 MATERIAL CF-136"]
