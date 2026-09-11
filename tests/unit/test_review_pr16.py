"""2026-09-11: the pre-merge review of PR #16, ten findings, each fixed and pinned here."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.comparative_cli import predict_core, spec_to_core
from src.kinetic_core.engine import (
    FormulationSpec, ProcessSpec, ThermalProgram, default_targets_for, predict,
)


def _pot(prec, t=140.0, m=20.0, ph=5.0, **kw):
    return FormulationSpec("p", prec, ProcessSpec(ThermalProgram.isothermal(t, m), ph=ph, **kw))


def test_1_methionine_and_proline_on_a_non_trunk_pot_no_longer_crash():
    """MET/PRO/PYRL are trunk species; on the sulfur and acrylamide lanes the declaration promised
    they were 'recorded and not charged' and the integrator then raised KeyError on them."""
    run = predict(_pot({"cysteine": 5, "ribose": 5, "methionine": 1}), ["2-methyl-3-furanthiol"])
    assert run.answered and run.concentrations_ug_per_l["2-methyl-3-furanthiol"] > 0
    assert run.run_metadata["recorded_not_charged"] == {"MET": 1.0}
    run = predict(_pot({"asparagine": 5, "glucose": 5, "proline": 1}, t=160), ["acrylamide"])
    assert run.answered and run.run_metadata["recorded_not_charged"] == {"PRO": 1.0}


@pytest.mark.parametrize("prec,target,t", [
    ({"cysteine": 5}, "2-methyl-3-furanthiol", 140.0),
    ({"thiamine": 5}, "2-furfurylthiol", 140.0),
    ({"asparagine": 5}, "acrylamide", 160.0),
    ({"glycine": 5}, "5-hydroxymethylfurfural", 120.0),
    ({"glucose": 0.0}, "5-hydroxymethylfurfural", 120.0),
])
def test_2_a_target_no_reaction_chain_can_reach_is_refused_not_answered_zero(prec, target, t):
    """The one rule the three earlier guards were projections of: a target must be REACHABLE from
    the charge in the lane's reaction set. Each of these used to answer exactly 0.0."""
    run = predict(_pot(prec, t=t, ph=6.0), [target])
    assert not run.answered
    assert "CHARGES NO PRECURSOR" in " ".join(run.declaration.reasons)


def test_2b_a_pentose_is_not_a_sulfur_source():
    run = predict(_pot({"ribose": 5}), ["2-methyl-3-furanthiol"])
    assert not run.answered and "NO sulfur source" in " ".join(run.declaration.reasons)


def test_2c_a_mixed_request_answers_what_it_can_and_names_the_rest():
    spec = _pot({"ribose": 100, "cysteine": 33}, t=145.0)
    run = predict(spec, ["2-methyl-3-furanthiol", "methanethiol"])
    assert run.answered
    assert "2-methyl-3-furanthiol" in run.concentrations_ug_per_l
    assert "methanethiol" not in run.concentrations_ug_per_l
    assert list(run.run_metadata["refused_targets"]) == ["methanethiol"]
    assert run.declaration.unreachable_targets == ("methanethiol",)
    assert any("NOT ANSWERED, BY NAME: 'methanethiol'" in w for w in run.declaration.warnings)
    # and the default target list is what the core can report for THIS charge, literally
    assert "methanethiol" not in default_targets_for({"D-Ribose": 100.0, "L-Cysteine": 33.0})
    assert "2-methyl-3-furanthiol (MFT)" in default_targets_for({"D-Ribose": 100.0, "L-Cysteine": 33.0})


def test_3_the_cli_carries_carried_volatiles_and_atmosphere_into_the_engine():
    spec = {"name": "s", "precursors": {"Pea Protein Isolate": 1000.0}, "temp_C": 140.0, "time_min": 0.1,
            "ph": 6.0, "matrix": "pea protein isolate", "carried_volatiles": {"hexanal": 331.0},
            "atmosphere": "air"}
    core = spec_to_core(spec)
    assert core.process.carried_volatiles == {"hexanal": 331.0} and core.process.atmosphere == "air"
    payload = predict_core({**spec, "targets": ["hexanal"]})
    assert payload["answered"]
    schema = json.loads(data_paths.SPEC_SCHEMA.read_text())
    assert "carried_volatiles" in schema["properties"] and "atmosphere" in schema["properties"]


def test_4_a_carried_level_declared_under_one_alias_applies_to_a_request_under_another():
    base = dict(matrix="pea protein isolate")
    spec = _pot({"Pea Protein Isolate": 1000.0}, t=140.0, m=0.1, ph=6.0,
                carried_volatiles={"2-pentylfuran": 59.4}, **base)
    a = predict(spec, ["2-pentyl furan"]).concentrations_ug_per_l["2-pentyl furan"]
    b = predict(spec, ["2-pentylfuran"]).concentrations_ug_per_l["2-pentylfuran"]
    bare = predict(_pot({"Pea Protein Isolate": 1000.0}, t=140.0, m=0.1, ph=6.0, **base),
                   ["2-pentylfuran"]).concentrations_ug_per_l["2-pentylfuran"]
    assert a == pytest.approx(b) and a > bare


def test_6_a_declared_zero_is_a_declaration():
    """A source that prints 'not detected' for the unheated control has declared a starting state."""
    kw = dict(matrix="pea protein isolate", t=40.0, m=10.0, ph=6.0,
              vessel=type("V", (), {"closure": "no cook"})())
    undeclared = predict(_pot({"Pea Protein Isolate": 1000.0}, **kw), ["hexanal"])
    declared0 = predict(_pot({"Pea Protein Isolate": 1000.0}, carried_volatiles={"hexanal": 0.0}, **kw), ["hexanal"])
    assert not undeclared.answered and "NEVER COOKED" in " ".join(undeclared.declaration.reasons)
    assert declared0.answered


def test_7_a_protein_loading_warns_that_it_moved_a_trunk_answer():
    loaded = predict(_pot({"glucose": 150, "glycine": 50}, t=120, m=30, ph=6.8, matrix="pea_isolate",
                          protein_g_per_l=30.0), ["5-hydroxymethylfurfural"])
    bare = predict(_pot({"glucose": 150, "glycine": 50}, t=120, m=30, ph=6.8), ["5-hydroxymethylfurfural"])
    assert loaded.concentrations_ug_per_l["5-hydroxymethylfurfural"] != bare.concentrations_ug_per_l["5-hydroxymethylfurfural"]
    assert any("BOUND-LYSINE POOL IS CHARGED" in w for w in loaded.declaration.warnings)
    assert not any("BOUND-LYSINE POOL IS CHARGED" in w for w in bare.declaration.warnings)


def test_8_the_carried_split_uses_the_engine_s_own_numbers():
    """No second integration: the split is a subtraction on what the engine applied and kept."""
    from src.kinetic_core.scoring import _carried_split
    spec = _pot({"Pea Protein Isolate": 1000.0}, t=140.0, m=0.1, ph=6.0, matrix="pea protein isolate",
                carried_volatiles={"hexanal": 331.0})
    run = predict(spec, ["hexanal"])
    total = run.concentrations_ug_per_l["hexanal"]
    split = _carried_split(run, 331.0, "hexanal", "ppb", total, 782.0)
    assert split["formed_predicted"] == pytest.approx(total - 331.0)
    assert split["formed_measured"] == pytest.approx(782.0 - 331.0)
    assert 0.0 < split["declared_share_of_prediction"] < 1.0


def test_10_b28_s_ship_rule_reads_two_frozen_panels_and_nothing_live():
    src = (data_paths.REPO_ROOT / "scripts" / "generators" / "generate_kinetic_core_b28_ship_rule.py").read_text()
    assert "score_panel()" not in src and '"git"' not in src
    for name in ("core_panel_scores_before_b28.json", "core_panel_scores_after_b28.json"):
        assert (data_paths.VALIDATION_DIR / "_b28_baseline" / name).exists()
    rule = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_b28_ship_rule.json").read_text())
    assert rule["T2"]["refused_before_after"] == [25, 18]   # the wave's own lift, no later wave's
