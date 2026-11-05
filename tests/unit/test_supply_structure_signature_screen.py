"""Fast tests of the pure functions of the supply-structure signature screen (no engine runs)."""
import math

import pytest

from scripts.generators import supply_structure_signature_screen as S


def test_cao_slope_bound_reproduces_the_prereg_numbers():
    r = S.cao_k_slope()
    assert r["slope"] == pytest.approx(0.0127, abs=5e-4)
    assert r["ci90_lo"] == pytest.approx(-0.155, abs=2e-3)
    assert r["ci90_hi"] == pytest.approx(0.180, abs=2e-3)
    assert r["dof"] == 4


def test_cao_slope_detects_a_planted_slope():
    ph = S.CAO_PH
    k = [60.0 * 10 ** (1.0 * (p - ph[-1])) for p in ph]
    r = S.cao_k_slope(k, ph, read_rel=0.0)
    assert r["slope"] == pytest.approx(1.0, abs=1e-9)
    assert r["ci90_lo"] > 0.9


def test_kallen_law_is_flat_below_pka_and_unit_slope_above():
    assert abs(S.kallen_shift_dex(6.24, 4.15, 4.67)) < 0.01
    assert S.kallen_shift_dex(5.2, 4.15, 4.67) == pytest.approx(0.0752, abs=1e-3)
    assert S.kallen_shift_dex(4.0, 7.0, 8.0) == pytest.approx(1.0, abs=1e-3)


def test_drain_ceiling_from_cao_pool_closure():
    assert S.drain_ceiling_per_min() == pytest.approx(60 / (65 * 300))
    assert S.drain_ceiling_per_min() >= 10 ** S.DRAIN_CEILING_LOG10   # the band's 3e-3 sits inside the ceiling


def test_zhai_loss_ratio_measured_and_guards():
    assert S.zhai_loss_ratio(S.ZHAI_FRAC80["8"], S.ZHAI_FRAC80["5.5"]) == pytest.approx(9.6, abs=0.05)
    assert math.isnan(S.zhai_loss_ratio(0.0, 0.5))
    assert S.zhai_loss_ratio(0.5, 0.5) == pytest.approx(1.0)


def test_feasibility_on_the_measured_values_passes_everything():
    measured = {"cao_K_mean": 60.0, "cao_K_range_dex": 0.043,
                "zhai_lossratio_8_over_55": S.zhai_loss_ratio(S.ZHAI_FRAC80["8"], S.ZHAI_FRAC80["5.5"]),
                "zh_TTCA30": 17.7, "zh_TTCA50": 18.2, "zh_CAMA50": 1.61, "zh_Cys50": 0.8016, "zh_PENT50": 1.22}
    f = S.feasibility(measured)
    assert all(f.values())
    assert S.signature_verdicts(f) == {"a": True, "b": True, "c": True}


def test_feasibility_fails_a_pH_flat_ttca_loss_and_high_free_cysteine():
    f = S.feasibility({"cao_K_mean": 41.0, "cao_K_range_dex": 0.48, "zhai_lossratio_8_over_55": 1.0,
                       "zh_TTCA30": 13.6, "zh_TTCA50": 12.3, "zh_CAMA50": 1.13, "zh_Cys50": 17.4, "zh_PENT50": 14.6})
    assert f == {"a_level": False, "a_flat": False, "b_ph_order": False, "c_ttca": False, "c_cama": True,
                 "c_free": False}


def test_modify_rates_touches_only_its_keys():
    keys = ("r_ttca_cys", "r_ttca_form", "r_cama_form", "r_ttca_deg", "r_cys_thermal", "r_pent_caramel", "r_other")
    base = {k: 1.0 for k in keys}
    assert S.modify_rates(dict(base), {}, 373.15, 6.5) == base
    mod = S.mod_of("S3X", {"K": 0.075, "log10_drain90": -3.0, "Ea_drain": S.EA_LUMPED, "exit_n": 1.0})
    k = S.modify_rates(dict(base), mod, S.T_REF_DRAIN_K, 6.5)
    assert k["r_ttca_form"] == pytest.approx(1.25)
    assert k["r_cama_form"] == pytest.approx(10.0) and k["r_ttca_deg"] == pytest.approx(10.0)
    assert k["r_cys_thermal"] == pytest.approx(1.001) and k["r_pent_caramel"] == pytest.approx(1.001)
    assert k["r_ttca_cys"] == 1.0 and k["r_other"] == 1.0
    k = S.modify_rates(dict(base), S.mod_of("S2pow_rel", {"m": 1.0}), 373.15, S.PH_REF_IN_SITU + 1.0)
    assert k["r_ttca_cys"] == pytest.approx(0.1) and k["r_ttca_form"] == 1.0
