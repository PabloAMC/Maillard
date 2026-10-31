"""B7a-2 (2026-10-31; kinetic_core_dpo_partition_prereg.md sec. 13): the DPO beta-cleavage's saturating pH law is an
option; the default slope law is unchanged, bit for bit."""
import sys
from pathlib import Path

import pytest

from src.kinetic_core import parameters_sulfur as ps

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))


def test_the_default_slope_law_is_bit_identical():
    entry = ps.PentoseAmineEntry(aa_ph_slope=0.15)
    for ph in (4.6, 5.5, 6.0, 7.43, 8.0):
        assert entry.aa_ph_factor(ph) == 10.0 ** (0.15 * (float(ph) - 6.0))


def test_the_saturating_law_is_one_at_ph6_and_flat_above_its_pk():
    entry = ps.PentoseAmineEntry(aa_law="saturating", aa_pk=6.0)
    assert entry.aa_ph_factor(6.0) == pytest.approx(1.0)
    assert entry.aa_ph_factor(5.0) == pytest.approx(2.0 / 11.0)
    assert entry.aa_ph_factor(10.0) == pytest.approx(2.0, rel=1e-3)
    assert entry.aa_ph_factor(8.0) / entry.aa_ph_factor(7.0) < 1.1


def test_an_unknown_law_is_refused():
    with pytest.raises(ValueError):
        ps.PentoseAmineEntry(aa_law="linear").aa_ph_factor(6.0)


def test_the_arms_swap_the_slope_for_the_pk():
    import fit_joint_deoxyosone_thiol as joint
    keys, _, lower, upper, start = joint.coordinates(False, joint.ARMS["B7AS"])
    assert "aa_ph_slope" not in keys and "aa_pk" in keys
    i = keys.index("aa_pk")
    assert (lower[i], upper[i], start[i]) == (5.0, 7.5, 6.1)
    entry = joint.trial(list(start), keys)["entry"]
    assert entry["aa_law"] == "saturating" and entry["aa_pk"] == pytest.approx(6.1)
    keys_n, *_ = joint.coordinates(False, joint.ARMS["R3En"])
    assert "aa_ph_slope" in keys_n and "aa_pk" not in keys_n
