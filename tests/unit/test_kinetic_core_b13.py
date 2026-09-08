"""Wave B13 (2026-09-07): glucosone, glyoxal and diacetyl on the trunk lane, trunk-only."""
from __future__ import annotations

import pytest

from src.kinetic_core import network, species
from src.kinetic_core.engine import (
    FormulationSpec, ProcessSpec, ThermalProgram, declare_envelope, predict,
)
from src.kinetic_core.parameters_dicarbonyl import DICARBONYL_KEYS, DICARBONYL_PARAMETERS


def _pot(name="pot", precursors=None, t_c=120.0, minutes=60.0):
    return FormulationSpec(name=name, precursors=precursors or {"D-Glucose": 100.0, "Glycine": 100.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes)))


def test_the_five_steps_run_on_the_trunk_only_and_balance():
    assert set(network.DICARBONYL_REACTION_KEYS) == {"r_glc_g", "r_g_go", "r_odg_da", "r_go_sink", "r_da_sink"}
    assert len(network.TRUNK_REACTIONS) == len(network.REACTIONS) + 5 + 5   # + B18's five pyrazine steps
    assert not set(network.DICARBONYL_REACTION_KEYS) & set(network.REACTION_KEYS)
    network.validate_balance(network.TRUNK_REACTIONS)          # raises on an unbalanced step
    from src.kinetic_core.sulfur import FULL_REACTION_KEYS
    assert not set(network.DICARBONYL_REACTION_KEYS) & set(FULL_REACTION_KEYS)


def test_the_species_are_appended_after_every_existing_one():
    keys = list(species.SPECIES_KEYS)
    assert keys[-8:-5] == ["G", "GO", "DA"]      # B18 appended PZ / DMP / MPZ / AKG / AKM after them
    assert species.INDEX["HMF"] < species.INDEX["G"]


def test_the_constants_are_kocadagli_re_referenced_from_180C():
    import math
    p = DICARBONYL_PARAMETERS["k_g_go"]
    # k at 180 C must be the published k_b
    assert p.k_at(453.15) == pytest.approx(737.0e-3, rel=1e-6)
    assert p.ea_kj_mol == pytest.approx(93.8)
    assert DICARBONYL_PARAMETERS["k_odg_da"].k_at(453.15) == pytest.approx(12.2e-3, rel=1e-6)
    assert DICARBONYL_PARAMETERS["k_glc_g"].k_at(453.15) == pytest.approx(0.069e-3, rel=1e-6)
    go_sink = DICARBONYL_PARAMETERS["k_go_sink"]
    assert go_sink.ea_kj_mol == 0.0 and go_sink.k_at(373.15) == go_sink.k_at(453.15)
    assert "ea_fixed_to_zero_by_authors" in go_sink.flags
    da_sink = DICARBONYL_PARAMETERS["k_da_sink"]
    assert da_sink.k_ref == 0.0 and "rate_zero_in_source" in da_sink.flags
    assert set(DICARBONYL_KEYS) == {"k_glc_g", "k_g_go", "k_odg_da", "k_go_sink", "k_da_sink"}


def test_a_trunk_pot_answers_the_three_and_the_old_answers_barely_move():
    run = predict(_pot(), ["5-HMF", "glyoxal", "diacetyl", "glucosone", "methylglyoxal"])
    assert run.answered, run.declaration.reasons
    assert run.require("glyoxal") > 0.0 and run.require("diacetyl") > 0.0
    # At 120 C the source's own constants make glucosone ACCUMULATE (its onward step is
    # 0.016 /min there against a 0.033 /min glyoxal sink at every temperature): a finding
    # of the wave, recorded in the prereg, not an assertion to force.
    assert run.require("glucosone") > 0.0
    # the reference prediction with the five steps switched off, for the 0.1 % pin
    from src.kinetic_core import engine
    from src.kinetic_core.engine import TRUNK, core_parameters
    spec = _pot()
    decl = declare_envelope(spec, ["5-HMF"])
    operative = core_parameters(TRUNK)
    off = {k: (v if k not in DICARBONYL_KEYS else v.__class__(**{**v.__dict__, "k_ref": 0.0})) for k, v in operative.items()}
    with_steps, _ = engine._integrate_program(TRUNK, operative, dict(decl.mapped_precursors), spec.process)
    without, _ = engine._integrate_program(TRUNK, off, dict(decl.mapped_precursors), spec.process)
    for key in ("HMF", "DMHF", "MGO", "TDG", "AMA"):
        assert abs(with_steps[key] / without[key] - 1.0) < 1e-3, key
    assert without["GO"] == 0.0 and with_steps["GO"] > 0.0


def test_a_sulfur_lane_request_for_a_dicarbonyl_is_refused_by_name():
    run = predict(_pot(precursors={"D-Ribose": 100.0, "L-Cysteine": 33.0}), ["2-methyl-3-furanthiol", "glyoxal"])
    assert not run.answered
    assert any("DICARBONYL" in r and "trunk" in r for r in run.declaration.reasons), run.declaration.reasons
