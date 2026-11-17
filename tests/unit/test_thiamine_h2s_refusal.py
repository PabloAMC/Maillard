"""Thiamine's H2S share where its zero is unmeasured: the outputs it moves are refused (STEP5_RECORD amendment,
2026-11-17; src/kinetic_core/thiamine_h2s_refusal.py).

Refused inside the region (thiamine, pH < 6.25, above 105 C), answered outside it (100 C; pH 7 held; no thiamine),
bit-identical where nothing moves; the reason names both sources and what would lift it; compare says why; the
switch off reproduces origin/main on every sulfur-lane panel pot (reference values computed on origin/main 50a224ab,
S3 off so that the model's own numbers are compared).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from src.comparative_cli import spec_to_core, validate_spec
from src.kinetic_core import config, engine, ph_balance_refusal, sulfur
from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import thiamine_h2s_refusal as thr

ROOT = Path(__file__).resolve().parents[2]
FFT, FUR, H2S = "2-furfurylthiol (FFT)", "furfural", "hydrogen sulfide"
#: thiamine + ribose, pH held at 5.0 (a clamped buffer: no S3 rivals to integrate), 140 C, 20 min
THI_PH5 = {"name": "thi_rib_pH5_140C", "precursors": {"thiamine": 10.0, "ribose": 30.0}, "temp_C": 140.0,
           "time_min": 20.0, "ph": 5.0, "aw": 0.98, "buffer": {"kind": "clamped"}}
THI_PH7 = dict(THI_PH5, name="thi_rib_pH7_140C", ph=7.0)
THI_100C = dict(THI_PH5, name="thi_rib_pH5_100C", temp_C=100.0)
CYS_PH5 = dict(THI_PH5, name="cys_rib_pH5_140C", precursors={"cysteine": 10.0, "ribose": 30.0})


def _core(spec):
    return spec_to_core(validate_spec(spec, label="s"))


def _same(a, b) -> bool:
    return (a.declaration == b.declaration and dict(a.concentrations_ug_per_l) == dict(b.concentrations_ug_per_l)
            and dict(a.species_mmol_per_l) == dict(b.species_mmol_per_l))


def test_the_per_arm_threshold_is_half_of_omega():
    th = thr.thresholds()
    assert th["omega_dex"] == ph_balance_refusal.thresholds()["omega_dex"]
    assert thr.threshold_dex() == pytest.approx(0.5 * th["omega_dex"])


def test_the_rival_covers_only_the_unmeasured_region_and_is_off_when_shipped():
    rival = ps.ThiamineH2sUnmeasuredRival(enabled=True)
    assert not ps.THIAMINE_H2S_UNMEASURED_RIVAL.enabled
    assert config.Config.shipped()["sulfur.thiamine_h2s_unmeasured_rival"] == ps.ThiamineH2sUnmeasuredRival()
    assert rival.covers(6.0, 120.0 + 273.15) and rival.covers(4.0, 140.0 + 273.15)
    assert not rival.covers(6.25, 140.0 + 273.15)          # the measured log-linear share starts at the floor
    assert not rival.covers(6.0, 105.0 + 273.15)           # at or below 105 C the zero is the measured one
    assert rival.share == pytest.approx(0.36)
    for needle in ("Dwivedi & Arnold 1971", "1972", "0.014"):
        assert needle in rival.source


def test_the_rate_constants_move_only_inside_the_region():
    params = engine.core_parameters(engine.SULFUR)
    on = config.Config.from_globals().replace(
        {"sulfur.thiamine_h2s_unmeasured_rival": ps.ThiamineH2sUnmeasuredRival(enabled=True)})
    for t_c, ph, moves in ((140.0, 5.0, True), (100.0, 5.0, False), (140.0, 7.0, False), (140.0, 6.3, False)):
        base = sulfur.sulfur_rate_constants_at(params, t_c + 273.15, ph)
        with config.applied(on):
            rival = sulfur.sulfur_rate_constants_at(params, t_c + 273.15, ph)
        assert {k: v for k, v in rival.items() if k != "r_thi_h2s"} == {k: v for k, v in base.items()
                                                                          if k != "r_thi_h2s"}
        if moves:
            assert base["r_thi_h2s"] == 0.0
            assert rival["r_thi_h2s"] == pytest.approx(0.36 / 0.64 * (base["r_thi_hmp"] + base["r_thi_mesh"]))
        else:
            assert rival["r_thi_h2s"] == base["r_thi_h2s"]


def test_refused_inside_the_region_with_the_reason():
    spec = _core(THI_PH5)
    run = engine.predict(spec, [FFT])
    assert not run.answered and not run.concentrations_ug_per_l
    reason = run.declaration.reasons[0]
    assert reason.startswith(thr.THIAMINE_REFUSED_TAG)
    for needle in ("Dwivedi & Arnold 1971", "Dwivedi & Arnold 1972", "0.014", "105 C", "pH 6.25", "P3",
                   "thiamine series", "WHAT WOULD LIFT IT", FFT):
        assert needle in reason, needle
    unrefused = engine.predict(spec, [FFT], thiamine_refusal=False)
    assert unrefused.answered          # the shipped zero's number, which the rule withholds


def test_a_partly_refused_request_drops_the_moved_output_by_name_and_keeps_the_rest():
    spec = _core(THI_PH5)
    full = engine.predict(spec, [FFT, FUR], thiamine_refusal=False)
    run = engine.predict(spec, [FFT, FUR])
    verdict = thr.cached_verdict(spec, [FFT, FUR])
    assert verdict is not None and FFT in verdict.refused and FUR not in verdict.refused
    assert run.answered and FFT not in run.concentrations_ug_per_l
    assert run.concentrations_ug_per_l[FUR] == full.concentrations_ug_per_l[FUR]
    assert run.run_metadata["refused_targets"][FFT].startswith(thr.THIAMINE_REFUSED_TAG)
    assert run.declaration.mapped_targets[FFT] not in run.species_mmol_per_l
    assert run.run_metadata["thiamine_h2s_refusal"]["refused"] == [FFT]


@pytest.mark.parametrize("pot", [THI_100C, THI_PH7, CYS_PH5], ids=["100C", "pH7_held", "no_thiamine"])
def test_answered_and_bit_identical_outside_the_region(pot):
    spec = _core(pot)
    run = engine.predict(spec, [FFT, H2S])
    assert run.answered and set(run.concentrations_ug_per_l) == {FFT, H2S}
    assert _same(run, engine.predict(spec, [FFT, H2S], thiamine_refusal=False))
    assert "thiamine_h2s_refusal" not in run.run_metadata
    if pot is not THI_PH7:
        assert not thr.applies(spec, run.declaration)       # no rival integrated at all


def test_compare_says_why_a_thiamine_arm_is_refused():
    out = engine.compare(_core(THI_PH5), _core(CYS_PH5), [FFT])
    assert out["comparable"] is False
    assert out["reason"].startswith(thr.THIAMINE_REFUSED_TAG)


def test_the_switch_off_is_the_unrefused_engine():
    spec = _core(THI_PH5)
    with config.applied(config.Config.from_globals().replace({"sulfur.thiamine_h2s_refusal": False})):
        off = engine.predict(spec, [FFT])
    assert off.answered and _same(off, engine.predict(spec, [FFT], thiamine_refusal=False))


def test_s3_rival_points_never_decide_or_read_a_thiamine_verdict():
    """S3's own integrations are the model's numbers (thiamine_refusal=False), and a verdict is keyed by the Config, so
    a verdict decided under an S3 rival point can never stand for the shipped one (found 2026-11-17 in the directional
    panel: CER07's spreads read 2.32 dex under the leak, 0.53 at the shipped point)."""
    spec = _core(dict(THI_PH5, name="thi_rib_pH5_unbuffered", buffer={"kind": "none"}))
    (changes, drift), = list(ph_balance_refusal.rival_points().values())[:1]
    thr._CACHE.clear()
    with config.applied(config.Config.from_globals().replace(changes)):
        ph_balance_refusal._rival_values((spec, (FFT,), drift))
        assert not thr._CACHE
        thr.assess(spec, [FFT])
    assert thr.cached_verdict(spec, [FFT]) is None          # decided under the rival's Config, not the shipped one


#: origin/main (50a224ab), every sulfur-lane pot of the core panel, all its targets, S3 off: {bundle: {target: hex}}
REFERENCE = json.loads((ROOT / "tests" / "fixtures" / "thiamine_refusal_origin_main_sulfur_panel.json").read_text())


def test_switch_off_is_bit_identical_to_origin_main_on_the_sulfur_panel_pots():
    from src.kinetic_core.panel import bundle_targets, core_spec, load_bundle

    assert len(REFERENCE) >= 20
    with config.applied(config.Config.from_globals().replace({"sulfur.thiamine_h2s_refusal": False})):
        for bundle, expected in REFERENCE.items():
            bench = load_bundle(ROOT / bundle)
            run = engine.predict(core_spec(bench, use_buffer=True), list(bundle_targets(bench)),
                                 balance_refusal=False)
            got = {c: float(v).hex() for c, v in run.concentrations_ug_per_l.items()}
            assert got == expected, bundle
