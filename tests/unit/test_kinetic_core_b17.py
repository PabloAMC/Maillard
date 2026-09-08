"""Wave B17 (2026-09-08): the disulfide gives the thiol back -- the engine side (kinetic_core_b17_prereg.md)."""
from __future__ import annotations

import pytest

from src.kinetic_core import sulfur
from src.kinetic_core.engine import SULFUR, core_parameters, frozen_parameters
from src.kinetic_core.parameters_sulfur import (
    DIMER_RELEASE_BOUNDS_LOG10K, MEASURED_SULFUR, ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL, dimer_release_parameters,
)


def test_the_two_release_steps_exist_share_one_constant_and_balance():
    by_key = {r.key: r for r in sulfur.SULFUR_REACTIONS}
    assert by_key["ch_dimer_release_mft"].reactants == {"MFTD": 1} and by_key["ch_dimer_release_mft"].products == {"MFT": 2}
    assert by_key["ch_dimer_release_fft"].reactants == {"FFTD": 1} and by_key["ch_dimer_release_fft"].products == {"FFT": 2}
    assert by_key["ch_dimer_release_mft"].parameter_key == by_key["ch_dimer_release_fft"].parameter_key == "k_dimer_release"
    sulfur.validate_sulfur_balance(sulfur.FULL_REACTIONS)


def test_the_release_constant_is_inert_by_default_and_carries_the_dimerisation_barrier():
    p = MEASURED_SULFUR["k_dimer_release"]
    assert p.k_ref == 0.0 and p.ea_kj_mol == ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL and p.order == 1
    shipped = core_parameters(SULFUR)
    assert shipped["k_dimer_release"].k_ref == 0.0          # no shipped report carries a B17 block
    assert "dimer_release_log10_k" not in frozen_parameters(SULFUR)
    assert DIMER_RELEASE_BOUNDS_LOG10K == (-10.0, 0.5)


def test_an_override_block_switches_the_release_on_and_only_that():
    on = core_parameters(SULFUR, frozen={"dimer_release_log10_k": {"k_dimer_release": -3.0}})
    off = core_parameters(SULFUR)
    assert on["k_dimer_release"].k_ref == pytest.approx(1e-3)
    assert on["k_dimer_release"].ea_kj_mol == ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL
    for key in off:
        if key != "k_dimer_release":
            assert on[key] == off[key], key
    k = sulfur.sulfur_rate_constants_at(on, 418.15, 5.0)
    assert k["ch_dimer_release_mft"] == pytest.approx(1e-3) and k["ch_dimer_release_fft"] == pytest.approx(1e-3)
    assert dimer_release_parameters()["k_dimer_release"].k_ref == 0.0
