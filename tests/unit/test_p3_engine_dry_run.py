"""The P3 engine dry run's design bookkeeping (no integration)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))

import build_p3_engine_dry_run as D  # noqa: E402


def test_each_lever_has_four_matched_pairs_differing_only_in_that_lever():
    for idx, lever in enumerate(("sugar identity", "sulfur source", "pH")):
        pairs = D._pairs(lever)
        assert len(pairs) == 4
        assert len({lo for lo, _ in pairs}) == 4
        for lo, hi in pairs:
            diff = [i for i in range(3) if lo[i] != hi[i]]
            assert diff == [idx]
            assert lo[idx] == (D.SUGARS, D.SULFUR, D.PH)[idx][0]


def test_isolate_reading_charges_protein_and_added_precursors():
    spec = D.corner_spec("ribose", "thiamine", 7.0, "pea_isolate")
    assert spec.precursors["ribose"] == D.SUGAR_MM and spec.precursors["thiamine"] == D.SULFUR_MM
    assert spec.process.matrix == "pea_isolate" and spec.process.protein_g_per_l > 100
    assert spec.process.buffer.kind == "none"


def test_water_reading_and_buffers():
    assert D.corner_spec("glucose", "cysteine", 6.0, "added_precursors").process.matrix == "water"
    assert D.corner_spec("glucose", "cysteine", 6.0, "pea_isolate_ph_held").process.buffer.kind == "clamped"
    phos = D.corner_spec("glucose", "cysteine", 6.0, "pea_isolate", phosphate=True).process.buffer
    assert phos.kind == "phosphate" and phos.phosphate_mol_l == D.PHOSPHATE_M
