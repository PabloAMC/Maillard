"""
The contrast benchmark scorer's frozen rules (prereg v6, results/validation/kinetic_core_contrast_benchmark_score_prereg.md):
the name map, the pair rules, the censored inclusion, the zero / refusal rule and the hand table. Envelope and YAML only;
nothing here integrates the engine.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import pytest

from src import data_paths

sys.path.insert(0, str(data_paths.REPO_ROOT / "scripts" / "generators"))
import score_contrast_benchmark as S  # noqa: E402


def test_name_map_follows_the_rule():
    assert S.precursor_name("GSH") == "glutathione"
    assert S.precursor_name("LYS") == "lysine"     # mapped; the envelope refuses LYS + PENT as a lane conflict
    assert S.precursor_name("TPP") is None
    assert S.target_name("MFTD") is not None and S.target_name("ACTZ") == "2-acetylthiazole"
    assert S.target_name("V_FURFURYL_ALCOHOL") is None


def test_restating_field_is_not_a_second_factor_but_a_scalar_one_is():
    base = {"initial_mmol_L": {"Cys": 400, "PENT": 400}, "aldehyde": "none", "label_ph": 5.0}
    hexanal = {"initial_mmol_L": {"Cys": 400, "PENT": 400, "V_HEXANAL": 100}, "aldehyde": "hexanal", "label_ph": 5.0}
    assert S.factors_between(base, hexanal, broad=True) == ["V_HEXANAL"]
    lipid_and_buffer = dict(base, lipid="1 %", buffer="phosphate")
    assert len(S.factors_between(dict(base, buffer="pyro"), lipid_and_buffer, broad=True)) == 2
    swap = {"initial_mmol_L": {"Cys": 400, "GLC": 400}, "label_ph": 5.0}
    assert len(S.factors_between(base, swap, broad=False)) == 2


@pytest.mark.parametrize("ea, eb, kind", [
    ({"value": 2.0, "sd": 0.2, "n": 3}, {"value": 1.0, "sd": 0.1, "n": 3}, "detected"),
    ({"value": 1.9, "sd": 0.4, "n": 3}, {"nd": True, "bound": 1.5}, "censored"),     # weighted, no 1.4x cut
    ({"value": 1.2, "sd": 1.1, "n": 3}, {"nd": True, "bound": 1.5}, "excluded"),     # below the other arm's limit
    ({"nd": True, "bound": 0.2}, {"nd": True, "bound": 0.2}, "excluded"),
    ({"value": 1.0}, {"pnq": True}, "excluded"),
])
def test_censored_inclusion(ea, eb, kind):
    assert S._measured("heldout_grosch1992_tab2", ea, eb)["kind"] == kind


def test_sign_weight_follows_the_papers_noise():
    sure = S._measured("heldout_grosch1992_tab2", {"value": 10.0, "sd": 0.5, "n": 3}, {"value": 1.0, "sd": 0.1, "n": 3})
    noisy = S._measured("heldout_grosch1992_tab2", {"value": 2.3, "sd": 1.0, "n": 4}, {"value": 1.2, "sd": 1.1, "n": 3})
    assert sure["weight"] > 0.99 and 0.0 < noisy["weight"] < sure["weight"]
    a, b = {"value": 4.7, "sd": 1.2, "n": 3}, {"value": 5.2, "sd": 0.9, "n": 3}      # Grosch MFT no. 1 vs no. 3
    t_weight = S._measured("heldout_grosch1992_tab2", a, b)["weight"]
    normal = 2 * S._cdf(abs(math.log10(4.7 / 5.2)) / S._measured("heldout_grosch1992_tab2", a, b)["sigma_m"], math.inf) - 1
    assert t_weight < normal                                # Student-t for SDs from 3 trials
    xu = S._measured("heldout_xu2007_tab1", {"value": 2.0}, {"value": 1.0})         # stated CV < 20 %, n not stated
    assert xu["sigma_m"] == pytest.approx(math.hypot(0.2, 0.2) / math.log(10), rel=1e-3)


def test_zero_and_refusal_rule():
    assert S.p_of(None, None) == 0.0 and S.p_of(0.0, 0.0) == 0.0
    assert S.p_of(None, 1.0) is None and S.p_of(1.0, 0.0) is None
    assert S.p_of(1.0 + 1e-8, 1.0) == 0.0                  # below the cross-platform float floor
    assert S.p_of(1.01, 1.0, floor=0.01) == 0.0 and S.p_of(1.01, 1.0, floor=0.001) != 0.0   # the item's measured floor
    assert S.p_of(1.005, 1.0) != 0.0                       # a real, if small, prediction
    assert S.p_of(2.0, 1.0) == pytest.approx(math.log10(2.0))
    item = {"kind": "detected", "m": 0.5}
    assert S.s_hit(item, None) == 0.0 and S.s_hit(item, 0.0) == 0.0 and S.s_hit(item, 0.3) == 1.0
    assert math.isinf(S.e_err(item, None))
    assert S._median([1.0, math.inf, -math.inf, 2.0]) == 1.0   # rank-based: no NaN


def test_hand_table_covers_the_65_v2_test_contrasts_once():
    assert not set(S.ATLAS_COVERAGE) & set(S.ATLAS_SPECS)
    assert len(S.ATLAS_COVERAGE) + len(S.ATLAS_SPECS) == 65
    assert {c for c, _ in S.ATLAS_COVERAGE.values()} <= {
        "inputs_missing", "operator_missing", "species_absent", "precursor_unmapped", "factor_not_translated"}


def test_engine_value_reads_the_requested_alias(monkeypatch):
    """The v4 review's B8: predictions are keyed by the requested alias; a refusal is the only None."""
    class Pred:
        answered = True
        concentrations_ug_per_l = {"2-methyl-3-furanthiol": 4.0}
        def require(self, k):
            return self.concentrations_ug_per_l[k]
    monkeypatch.setattr(S.engine, "predict", lambda spec, targets: Pred())
    pot = {"initial_mmol_L": {"PENT": 35, "Cys": 50}, "t_c": 140, "minutes": 60, "label_ph": 5.0}
    assert S.engine_value("heldout_xu2007_tab1", "p", pot, "MFT", clamped=False) == (4.0, False)


def test_score_block_end_to_end_on_fabricated_predictions():
    """The score path with made-up numbers (no engine): floors, refusals, weights, intervals and readings compute."""
    items = [
        {"table": "heldout_grosch1992_tab2", "a": "a1", "b": "b1", "factor": "Cys", "species": "MFT", "ladder": "L1",
         "exposure": "clean", "methionine": False, "kind": "detected", "m": 0.5, "weight": 0.9, "sigma_m": 0.1,
         "refused_arm": None},
        {"table": "heldout_xu2007_tab1", "a": "a2", "b": "b2", "factor": "label_ph", "species": "FUR", "ladder": "L2",
         "exposure": "clean", "methionine": False, "kind": "censored", "m_sign": -1.0, "weight": 0.8, "refused_arm": None},
        {"table": "heldout_meynier1995_tab1_2_4", "a": "a3", "b": "b3", "factor": "label_ph", "species": "NF",
         "ladder": "L3", "exposure": "clean", "methionine": False, "kind": "detected", "m": -0.3, "weight": 0.5,
         "sigma_m": 0.2, "refused_arm": "a3"},
    ]
    preds = {("heldout_grosch1992_tab2", "a1", "MFT"): 2.0, ("heldout_grosch1992_tab2", "b1", "MFT"): 1.0,
             ("heldout_xu2007_tab1", "a2", "FUR"): 1.0, ("heldout_xu2007_tab1", "b2", "FUR"): 1.0 + 1e-4,
             ("heldout_meynier1995_tab1_2_4", "a3", "NF"): None, ("heldout_meynier1995_tab1_2_4", "b3", "NF"): 3.0}
    tight = dict(preds)
    tight[("heldout_xu2007_tab1", "b2", "FUR")] = 1.0 + 3e-3          # floor ~1.3e-3 dex > |p| ~4e-5: p = 0
    block = S.score_block(items, preds, tight, np.random.default_rng(0))
    assert block["S"]["items"] == 3 and block["S"]["zero_or_refused_one_arm"] == 1 and block["S"]["p_zero"] == 1
    assert block["S"]["share_weighted"] == pytest.approx(0.9 / 2.2)
    assert block["E"]["items"] == 2 and block["E"]["infinite_errors"] == 1
    reading = S.readings(block, {"curve": {"per_pair": {"0.5": 0.05}}, "false_pass_rate": 0.05})
    assert set(reading) >= {"direction", "log_ratio_error", "status"}
