"""Wave B27 (2026-09-11): the dicarbonyl redox couple, pre-registered, gated before the fit, NOT SHIPPED
(kinetic_core_b27_prereg.md; kinetic_core_b27_ship_rule.md). The structure stays in the code, inert."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import sulfur
from src.kinetic_core.engine import SULFUR, core_parameters, frozen_parameters
from src.kinetic_core.parameters_sulfur import (
    DICARBONYL_REDOX_BOUNDS_LOG10_YIELD, MEASURED_SULFUR, apply_dicarbonyl_redox, dicarbonyl_redox_parameters,
)

SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b27_ship_rule.json"


def test_the_redox_branch_exists_balances_and_carries_the_mercaptoketone_steps_ph_factor():
    by_key = {r.key: r for r in sulfur.SULFUR_REACTIONS}
    r = by_key["ch_redox_mp3p"]
    assert r.reactants == {"NF": 1, "H2S": 1} and r.products == {"MP3P": 1, "OX": 1} and r.parameter_key == "k_redox_mp3p"
    assert sulfur.REACTION_PH_FACTOR["ch_redox_mp3p"] == sulfur.REACTION_PH_FACTOR["r_nf_mp3p"] == "neutral_h2s"
    sulfur.validate_sulfur_balance(sulfur.FULL_REACTIONS)


def test_inert_by_default_and_absent_from_the_shipped_block():
    assert MEASURED_SULFUR["k_redox_mp3p"].k_ref == 0.0 and MEASURED_SULFUR["k_redox_mp3p"].order == 2
    assert core_parameters(SULFUR)["k_redox_mp3p"].k_ref == 0.0
    assert "dicarbonyl_redox" not in frozen_parameters(SULFUR)
    assert dicarbonyl_redox_parameters()["k_redox_mp3p"].k_ref == 0.0
    assert DICARBONYL_REDOX_BOUNDS_LOG10_YIELD == (-4.0, 0.0)


def test_phi_splits_the_mercaptoketone_flux_exactly_and_touches_nothing_else():
    off = core_parameters(SULFUR)
    on = core_parameters(SULFUR, frozen={"dicarbonyl_redox": {"log10_ox_yield_per_mercaptoketone": -0.5}})
    phi = 10.0 ** -0.5
    k = off["k_nf_mp3p"].k_ref
    assert on["k_nf_mp3p"].k_ref == pytest.approx((1.0 - phi) * k)
    assert on["k_redox_mp3p"].k_ref == pytest.approx(phi * k)
    assert on["k_nf_mp3p"].k_ref + on["k_redox_mp3p"].k_ref == pytest.approx(k, rel=1e-12)
    assert on["k_redox_mp3p"].ea_kj_mol == off["k_nf_mp3p"].ea_kj_mol
    for key in off:
        if key not in ("k_nf_mp3p", "k_redox_mp3p"):
            assert on[key] == off[key], key
    rates = sulfur.sulfur_rate_constants_at(on, 413.15, 4.5)
    assert rates["ch_redox_mp3p"] / rates["r_nf_mp3p"] == pytest.approx(phi / (1.0 - phi))
    with pytest.raises(ValueError):
        apply_dicarbonyl_redox(dict(off), 0.1)   # phi > 1 is not a fraction


def test_the_generator_installs_the_printed_charges_and_restores_them():
    import sys
    sys.path.insert(0, str(data_paths.REPO_ROOT / "scripts" / "generators"))
    import generate_kinetic_core_b2_3_fit as B23
    import generate_kinetic_core_b27_fit as B27  # configure() at import
    try:
        assert len(B23.ACTIVE_FIT_ROWS) == 65 and len(B27.ALL_KEYS) == 49 and B27.K_SLOT == 48
        assert B23.SYSTEMS["whitfield_nf_cys"]["initial"] == {"NF": 50.0, "Cys": 50.0}
        assert B23.SYSTEMS["whitfield_nf_h2s"]["initial"] == {"NF": 50.0, "H2S": 97.0}
        assert B23.SYSTEMS["whitfield_nf_cys"]["buffer"].declared and B23.SYSTEMS["whitfield_nf_cys"]["buffer"].phosphate_mol_l == 0.5
        row = next(r for r in B23.ACTIVE_FIT_ROWS if r["id"] == "whitfield_nf_cys_MFT")
        assert row["kind"] == "molpct_total" and row["species_terms"] == {"MFT": 1, "MFTD": 2} and row["target"] == 0.230
        floor = next(r for r in B23.ACTIVE_FIT_ROWS if r["id"] == "whitfield_mft_disulfide_share_floor")
        assert floor["kind"] == "floor" and floor["target"] == 0.35
    finally:
        B27.restore()
    assert len(B23.ACTIVE_FIT_ROWS) == 64 or "whitfield_mft_disulfide_share_floor" not in {r["id"] for r in B23.ACTIVE_FIT_ROWS}
    assert B23.SYSTEMS["whitfield_nf_cys"]["initial"] == {"NF": 20.0, "Cys": 20.0}
    assert next(r for r in B23.ACTIVE_FIT_ROWS if r["id"] == "whitfield_nf_cys_MFT")["target"] == 0.150


def test_the_gate_fired_before_the_fit_and_the_record_says_why():
    p = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    assert p["verdict"] == "NOT FITTED -- DO NOT SHIP" and p["fitted"] is False
    g1 = p["G1"]
    # The ambient pots are not oxidant-limited: consumers use under 1 % of the pool, and the whole
    # mercaptoketone flux at phi = 1 adds under 1 % -- against the ~9x a decade on the share needs.
    assert not g1["reachable"] and g1["best_available_at_phi_1"] < 0.01
    for pot in g1["pots"].values():
        assert pot["consumed_fraction"] < 0.01
    # The one pot the structure fixes reaches the 35 % floor only at the physical ceiling.
    assert p["G2"]["phi_first_reaching_floor"] == 1.0 and p["G2"]["at_ceiling"]
    # And section 9's targeted T2 held: the Kumazawa rows do not move, because they carry no norfuraneol.
    assert p["G3"]["kumazawa_max_abs_growth_dex"] < 0.05 and p["G3"]["n_over_0_3_dex"] == 0
