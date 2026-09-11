"""ENV-B34 (2026-09-11): prior rows for the 3-deoxyglucosone limb and the amine-free sugar entries,
banded on the source's printed 95 % HPD (kinetic_core_env_b34_prereg.md). No centre moves."""
from __future__ import annotations

import math

import numpy as np
import pytest

from src.kinetic_core import uncertainty as U
from src.kinetic_core.engine import TRUNK, core_parameters
from src.kinetic_core.parameters_dicarbonyl import AGREEING_SINK_KEYS, HPD_SINK_BANDS, HPD_SINK_KEYS, with_disputed_sinks
from src.kinetic_core.parameters_furanic import FURANIC_PARAMETERS


def test_the_four_constants_carry_the_printed_hpd_and_no_centre_moves():
    assert HPD_SINK_KEYS == ("k_glc_tdg", "k_tdg_ddg", "k_fru_int", "k_fru_odg")
    assert HPD_SINK_BANDS["k_glc_tdg"]["k_rel_hpd"] == pytest.approx(2.44 / 4.19)
    assert HPD_SINK_BANDS["k_tdg_ddg"]["ea_hpd_kj_mol"] == 6.3
    rows = {p.key: p for p in U.core_priors() if p.key.startswith("b34.")}
    assert len(rows) == 8 and all(p.sampled and p.distribution == "uniform_band" for p in rows.values())
    for key in HPD_SINK_KEYS:
        base = FURANIC_PARAMETERS[key]
        k = rows[f"b34.{key}.log10_k_100C"]; e = rows[f"b34.{key}.ea_kj_mol"]
        assert k.band[0] < math.log10(base.k_ref) < k.band[1]      # the centre is inside its own band
        assert e.band[0] < base.ea_kj_mol < e.band[1]
        assert k.centre == pytest.approx(math.log10(base.k_ref)) and e.centre == base.ea_kj_mol
    # the shipped parameters are untouched: a prior row is not a centre move
    assert core_parameters(TRUNK)["k_tdg_ddg"].k_ref == FURANIC_PARAMETERS["k_tdg_ddg"].k_ref


def test_k_tdg_ddg_left_the_agreeing_list_and_the_reason_says_why():
    from src.kinetic_core.parameters_dicarbonyl import AGREEING_SINK_REASON
    assert "k_tdg_ddg" not in AGREEING_SINK_KEYS and "k_ddg_hmf" in AGREEING_SINK_KEYS
    assert "32x" in AGREEING_SINK_REASON and "DRY" in AGREEING_SINK_REASON


def test_a_draw_reaches_the_engine_through_the_same_hook_env_b13_uses():
    d = U.draw_from_rng(np.random.default_rng(7), 0)
    block = d.core.maillard["disputed_sinks"]
    assert set(HPD_SINK_KEYS) <= set(block)
    op = core_parameters(TRUNK, frozen=d.core.maillard)
    base = core_parameters(TRUNK)
    for key in HPD_SINK_KEYS:
        rel = float(HPD_SINK_BANDS[key]["k_rel_hpd"]); hpd = float(HPD_SINK_BANDS[key]["ea_hpd_kj_mol"])
        assert base[key].k_ref * (1 - rel) * 0.999 <= op[key].k_ref <= base[key].k_ref * (1 + rel) * 1.001
        assert abs(op[key].ea_kj_mol - base[key].ea_kj_mol) <= hpd + 1e-9
    # the hook accepts both waves' keys and still refuses a stranger
    assert "k_glc_tdg" in with_disputed_sinks({"k_glc_tdg": {"ea_kj_mol": 90.0}})
    with pytest.raises(KeyError):
        with_disputed_sinks({"k_ddg_hmf": {"ea_kj_mol": 1.0}})
