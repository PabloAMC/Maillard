"""A thiol donor's effect on thiamine's MFT route is undecided: the outputs its sensitivity rival moves are refused
(STEP5_RECORD second amendment, 2026-11-26; B66; src/kinetic_core/cysteine_thiamine_mft_refusal.py).

Refused in a pot charged with thiamine and cysteine (or glutathione, or TTCA), answered without either, bit-identical where nothing moves; the
reason names the sources and what would lift it; compare says why; the switch off reproduces origin/main on every
sulfur-lane panel pot (the thiamine rule's fixture: origin/main 50a224ab, S3 off).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.comparative_cli import spec_to_core, validate_spec
from src.kinetic_core import config, engine, ph_balance_refusal, sulfur
from src.kinetic_core import cysteine_thiamine_mft_refusal as ctr
from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import thiamine_h2s_refusal as thr

ROOT = Path(__file__).resolve().parents[2]
MFT, FFT, FUR = "2-methyl-3-furanthiol (MFT)", "2-furfurylthiol (FFT)", "furfural"
#: thiamine + cysteine + ribose, pH held at 5.7 (clamped: no S3 rivals), 100 C (below the thiamine H2S rule's 105 C)
PAIR = {"name": "thi_cys_rib_pH5.7_100C", "precursors": {"thiamine": 10.0, "cysteine": 20.0, "ribose": 30.0},
        "temp_C": 100.0, "time_min": 60.0, "ph": 5.7, "aw": 0.98, "buffer": {"kind": "clamped"}}
THI_ONLY = dict(PAIR, name="thi_rib", precursors={"thiamine": 10.0, "ribose": 30.0})
CYS_ONLY = dict(PAIR, name="cys_rib", precursors={"cysteine": 20.0, "ribose": 30.0})
GSH_PAIR = dict(PAIR, name="thi_gsh_rib", precursors={"thiamine": 10.0, "glutathione": 10.0, "ribose": 30.0})
TTCA_PAIR = dict(PAIR, name="thi_ttca", precursors={"thiamine": 10.0, "ttca": 10.0})


def _core(spec):
    return spec_to_core(validate_spec(spec, label="s"))


def _same(a, b) -> bool:
    return (a.declaration == b.declaration and dict(a.concentrations_ug_per_l) == dict(b.concentrations_ug_per_l)
            and dict(a.species_mmol_per_l) == dict(b.species_mmol_per_l))


def test_the_threshold_is_the_thiamine_rules():
    assert ctr.threshold_dex() == thr.threshold_dex()
    assert ctr.threshold_dex() == pytest.approx(0.5 * ph_balance_refusal.thresholds()["omega_dex"])


def test_the_rival_is_off_when_shipped_and_names_its_sources():
    assert not ps.CYSTEINE_THIAMINE_MFT_RIVAL.enabled
    assert config.Config.shipped()["sulfur.cysteine_thiamine_mft_rival"] == ps.CysteineThiamineMftRival()
    assert ps.CYSTEINE_THIAMINE_MFT_RIVAL.factor == pytest.approx(7.0)
    for needle in ("Grosch & Zeiler-Hilgart 1992", "Table II", "x7.0", "glutathione", "Zhang et al. 2024",
                   "10.1021/bk-1992-0490.ch014"):
        assert needle in ps.CYSTEINE_THIAMINE_MFT_RIVAL.source


def test_the_rival_moves_only_the_thiamine_mft_step():
    params = engine.core_parameters(engine.SULFUR)
    on = config.Config.from_globals().replace(
        {"sulfur.cysteine_thiamine_mft_rival": ps.CysteineThiamineMftRival(enabled=True)})
    for t_c, ph in ((100.0, 5.7), (140.0, 5.0)):
        base = sulfur.sulfur_rate_constants_at(params, t_c + 273.15, ph)
        with config.applied(on):
            rival = sulfur.sulfur_rate_constants_at(params, t_c + 273.15, ph)
        assert {k: v for k, v in rival.items() if k != "r_hmp_mft"} == {k: v for k, v in base.items()
                                                                          if k != "r_hmp_mft"}
        assert rival["r_hmp_mft"] == pytest.approx(7.0 * base["r_hmp_mft"])


def test_mft_is_refused_in_a_thiamine_plus_cysteine_pot_with_the_reason():
    spec = _core(PAIR)
    run = engine.predict(spec, [MFT])
    assert not run.answered and not run.concentrations_ug_per_l
    reason = run.declaration.reasons[0]
    assert reason.startswith(ctr.CYSTEINE_THIAMINE_REFUSED_TAG)
    for needle in ("Grosch & Zeiler-Hilgart 1992", "Zhang et al. 2024", "x7.0", "x19", "x0.38", "glutathione",
                   "WHAT WOULD LIFT IT", "cysteine_boost_thiol_yield_concentrated", MFT):
        assert needle in reason, needle
    unrefused = engine.predict(spec, [MFT], cysteine_thiamine_refusal=False)
    assert "HIGHER" not in reason and "probably" not in reason      # no direction is claimed (premortem cause 5)
    assert unrefused.answered          # the shipped number, which the rule withholds


def test_a_partly_refused_request_drops_mft_by_name_and_keeps_what_does_not_move():
    spec = _core(PAIR)
    full = engine.predict(spec, [MFT, FFT, FUR], cysteine_thiamine_refusal=False)
    run = engine.predict(spec, [MFT, FFT, FUR])
    verdict = ctr.cached_verdict(spec, [MFT, FFT, FUR])
    assert verdict is not None and verdict.refused == (MFT,)
    assert verdict.spread_dex[MFT] > 0.3 and verdict.spread_dex[FFT] < 0.05 and verdict.spread_dex[FUR] < 0.05
    assert run.answered and MFT not in run.concentrations_ug_per_l
    for c in (FFT, FUR):
        assert run.concentrations_ug_per_l[c] == full.concentrations_ug_per_l[c]
    assert run.run_metadata["refused_targets"][MFT].startswith(ctr.CYSTEINE_THIAMINE_REFUSED_TAG)
    assert run.declaration.mapped_targets[MFT] not in run.species_mmol_per_l
    assert run.run_metadata["cysteine_thiamine_mft_refusal"]["refused"] == [MFT]


@pytest.mark.parametrize("pot", [GSH_PAIR, TTCA_PAIR], ids=["glutathione", "ttca"])
def test_the_other_thiol_donors_trigger_the_rule_too(pot):
    """Grosch 1992 Table II: glutathione raises thiamine's MFT x4.4 (Zhang 2024: x15); TTCA opens to cysteine in the pot."""
    spec = _core(pot)
    run = engine.predict(spec, [MFT])
    assert ctr.applies(spec, engine.declare_envelope(spec, [MFT]))
    assert not run.answered and run.declaration.reasons[0].startswith(ctr.CYSTEINE_THIAMINE_REFUSED_TAG)


@pytest.mark.parametrize("pot", [THI_ONLY, CYS_ONLY], ids=["no_donor", "no_thiamine"])
def test_answered_and_bit_identical_without_either_partner(pot):
    spec = _core(pot)
    run = engine.predict(spec, [MFT, FFT])
    assert run.answered and set(run.concentrations_ug_per_l) == {MFT, FFT}
    assert _same(run, engine.predict(spec, [MFT, FFT], cysteine_thiamine_refusal=False))
    assert "cysteine_thiamine_mft_refusal" not in run.run_metadata
    assert not ctr.applies(spec, run.declaration)           # no rival integrated at all


def test_compare_says_why_the_pair_is_refused():
    out = engine.compare(_core(PAIR), _core(THI_ONLY), [MFT])
    assert out["comparable"] is False
    assert out["reason"].startswith(ctr.CYSTEINE_THIAMINE_REFUSED_TAG)


def test_the_switch_off_is_the_unrefused_engine():
    spec = _core(PAIR)
    with config.applied(config.Config.from_globals().replace({"sulfur.cysteine_thiamine_mft_refusal": False})):
        off = engine.predict(spec, [MFT])
    assert off.answered and _same(off, engine.predict(spec, [MFT], cysteine_thiamine_refusal=False))


def test_the_other_rules_own_integrations_never_decide_a_verdict_of_this_rule():
    """S3's rival points and the thiamine rule's two runs are the model's own numbers (cysteine_thiamine_refusal=False),
    and a verdict is keyed by the Config, so one decided under a rival's Config never stands for the shipped one."""
    spec = _core(dict(PAIR, name="thi_cys_rib_unbuffered_140C", temp_C=140.0, ph=5.0, buffer={"kind": "none"}))
    (changes, drift), = list(ph_balance_refusal.rival_points().values())[:1]
    ctr._CACHE.clear()
    with config.applied(config.Config.from_globals().replace(changes)):
        ph_balance_refusal._rival_values((spec, (MFT,), drift))
    thr._CACHE.clear()
    thr.assess(spec, [FFT])
    assert not ctr._CACHE


def test_on_the_sulfur_panel_the_rule_refuses_only_in_pots_charged_with_both():
    """Every sulfur-lane panel pot (the thiamine rule's fixture lists them; its own test shows both switches off
    reproduce origin/main): with this rule on, an output is dropped only where thiamine and a donor are both charged,
    and what is dropped there is MFT or its dimer."""
    from src.kinetic_core.panel import bundle_targets, core_spec, load_bundle

    reference = json.loads((ROOT / "tests" / "fixtures" / "thiamine_refusal_origin_main_sulfur_panel.json").read_text())
    dropped = {}
    with config.applied(config.Config.from_globals().replace({"sulfur.thiamine_h2s_refusal": False})):
        for bundle, expected in reference.items():
            bench = load_bundle(ROOT / bundle)
            spec, targets = core_spec(bench, use_buffer=True), list(bundle_targets(bench))
            run = engine.predict(spec, targets, balance_refusal=False)
            gone = set(expected) - set(run.concentrations_ug_per_l)
            both = ctr.applies(spec, engine.declare_envelope(spec, targets))
            assert both or not gone, bundle
            if gone:
                dropped[bundle] = gone
    assert dropped, "the panel holds thiamine + cysteine pots (Bolton 1994 at least)"
    assert any("Bolton1994" in b for b in dropped)
    for bundle, gone in dropped.items():
        assert all("furanthiol" in c.lower() or "furyl) disulfide" in c.lower() for c in gone), (bundle, gone)
