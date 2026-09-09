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


# ---- variant (a), 2026-09-09: the pot makes its own electrophile sites --------------------------
from src.kinetic_core.parameters_sulfur import MELE_SITE_YIELD_BOUNDS_LOG10, mele_site_parameters  # noqa: E402


def test_variant_a_the_three_site_steps_exist_share_one_constant_and_balance():
    by_key = {r.key: r for r in sulfur.SULFUR_REACTIONS}
    for osone in ("DPO", "TDP", "DDP"):
        r = by_key[f"ch_mele_from_{osone.lower()}"]
        assert r.reactants == {osone: 1} and r.products == {"FRAG_C": 5, "MELE": 1} and r.parameter_key == "k_mele_site"
    sulfur.validate_sulfur_balance(sulfur.FULL_REACTIONS)


def test_variant_a_the_site_constant_is_inert_by_default():
    p = MEASURED_SULFUR["k_mele_site"]
    assert p.k_ref == 0.0 and p.order == 1
    assert core_parameters(SULFUR)["k_mele_site"].k_ref == 0.0     # no shipped report carries a B17a block
    assert "mele_site_log10_yield" not in frozen_parameters(SULFUR)
    assert MELE_SITE_YIELD_BOUNDS_LOG10 == (-4.0, 0.2)
    assert mele_site_parameters()["k_mele_site"].k_ref == 0.0


def test_variant_a_an_override_block_switches_the_sites_on_as_yield_times_k_osone_decay():
    off = core_parameters(SULFUR)
    on = core_parameters(SULFUR, frozen={"mele_site_log10_yield": {"mele_site_yield": -1.0}})
    fr = frozen_parameters(SULFUR)
    expected = 0.1 * 10.0 ** fr["log10_k_ref_at_145C"]["k_osone_decay"]
    assert on["k_mele_site"].k_ref == pytest.approx(expected)
    assert on["k_mele_site"].ea_kj_mol == pytest.approx(fr["decay_Ea_kJ_mol"]["carbonyl_sink"])
    for key in off:
        if key != "k_mele_site":
            assert on[key] == off[key], key
    k = sulfur.sulfur_rate_constants_at(on, 418.15, 5.0)
    for osone in ("dpo", "tdp", "ddp"):
        assert k[f"ch_mele_from_{osone}"] == pytest.approx(expected)
