"""Leave-one-laboratory-out of the thiol lane (results/validation/kinetic_core_thiol_lolo_prereg.md sec. 3): no
laboratory's rows enter its own fold's fit, every row is held out exactly once, and a reprint counts as its
original's laboratory."""
import sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))

import fit_thiol_lolo as lolo  # noqa: E402


def _row(ledger_id, laboratory, i):
    return SimpleNamespace(ledger_id=ledger_id, laboratory=laboratory, i=i)


def _check(rows, folds):
    held_once = Counter(id(r) for _, held in folds.values() for r in held)
    assert held_once == Counter(id(r) for r in rows), "every row is held out exactly once"
    for lab, (train, held) in folds.items():
        assert held and all(lolo.lab_of(r) == lab for r in held)
        assert not any(lolo.lab_of(r) == lab for r in train), f"{lab}'s rows leak into its own fit"
        assert len(train) + len(held) == len(rows)
        assert {id(r) for r in train}.isdisjoint({id(r) for r in held})


def test_folds_hold_out_exactly_one_laboratory_and_never_fit_it():
    rows = [_row("zhai2021_fig5", "jiangnan_zhang_x", 0), _row("davidek2008_fig1", "nestle_blank", 1),
            _row("zhou2022_fig1", "jiangnan_zhang_x", 2), _row("hofmann1998_t4", "dfa_hofmann_schieberle", 3),
            _row("schieberle2000_t4", "dfa_hofmann_schieberle", 4)]
    folds = lolo.make_folds(rows)
    assert sorted(folds) == ["dfa_hofmann_schieberle", "jiangnan_zhang_x", "nestle_blank"]
    _check(rows, folds)
    assert [r.i for r in folds["jiangnan_zhang_x"][1]] == [0, 2]


def test_a_reprint_is_held_out_with_its_original_laboratory():
    # Kang 2026 (Zhejiang Gongshang) reprints Jiangnan's Zhai 2023 Food Chem data (CF-005): never a second laboratory
    rows = [_row("zhai2023foodchem_table1", "jiangnan_zhang_x", 0),
            _row("kang2026_si_table_s4_ttca_tladder", "zhejiang_gongshang_zhai_fu", 1),
            _row("desclaux2006_fig2", "reading_mottram", 2)]
    folds = lolo.make_folds(rows)
    assert "zhejiang_gongshang_zhai_fu" not in folds
    assert [r.i for r in folds["jiangnan_zhang_x"][1]] == [0, 1]
    assert [r.i for r in folds["reading_mottram"][0]] == [0, 1]
    _check(rows, folds)


def test_the_real_fit_table_has_six_laboratories_and_no_leak():
    data, _, _ = lolo.load_fit_table()   # refuses to run unless it equals R3En's 521-row table
    folds = lolo.make_folds(data)
    assert {lab: len(held) for lab, (_, held) in folds.items()} == {
        "jiangnan_zhang_x": 454, "dalian_polytechnic_chen_d": 19, "nestle_blank": 19, "dfa_hofmann_schieberle": 16,
        "reading_mottram": 8, "ochanomizu_murata": 5}
    _check(data, folds)
