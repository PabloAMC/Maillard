"""B7 hub (2026-10-31; kinetic_core_dpo_partition_prereg.md sec. 20): the cysteine trap of 3-DX is a switch -- off,
the reaction is not in the network; on, it is balanced and the fit's install() sets and resets it."""
import sys
from pathlib import Path

import pytest

from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import sulfur

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))


def _trap():
    return next(r for r in sulfur.FULL_REACTIONS if r.key == "r_tdp_cys_trap")


def test_the_trap_is_carbon_nitrogen_and_sulfur_balanced_in_the_lumped_fragments():
    r = _trap()
    assert r.reactants == {"TDP": 1, "Cys": 1}
    assert r.products == {"CBX": 1, "FRAG_C": 7, "FRAG_N": 1, "FRAG_S": 1}


def test_the_trap_is_off_by_default_and_its_rate_scales_with_the_constant():
    assert ps.TDP_CYSTEINE_TRAP.k_l_per_mmol_min_100c_ph5 == 0.0
    on = ps.TdpCysteineTrap(k_l_per_mmol_min_100c_ph5=0.1)
    assert on.k_at(373.15) == pytest.approx(0.1)
    assert on.k_at(393.15) > on.k_at(373.15)


def test_install_sets_and_resets_the_trap():
    import fit_joint_deoxyosone_thiol as joint
    keys, _, lower, upper, start = joint.coordinates(False, joint.ARMS["B7ASST2C"])
    assert "log10_k_tdp_cys_trap" in keys
    job = joint.trial(list(start), keys)
    assert job["tdp_cys_trap_k"] == pytest.approx(10.0 ** -6.0)
    joint.install(dict(job, routes=job["routes"]))
    assert ps.TDP_CYSTEINE_TRAP.k_l_per_mmol_min_100c_ph5 == pytest.approx(1e-6)
    keys_n, _, _, _, start_n = joint.coordinates(False, joint.ARMS["R3En"])
    joint.install(joint.trial(list(start_n), keys_n))
    assert ps.TDP_CYSTEINE_TRAP.k_l_per_mmol_min_100c_ph5 == 0.0
