"""
Wave B26 -- the plant-protein binding row, and the guard that made it safe to add.

The frozen literals here are read from `results/validation/kinetic_core_b26_ship_rule.md` and
`kinetic_core_b26_prereg.md`. If the registry moves, these fail, which is the point.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from src.kinetic_core.matrix_oav import (
    fit_class_binding_constants,
    fit_unsaturation_penalty,
    predict_matrix_shift,
)
from src.kinetic_core.parameters_matrix import (
    CHAIN_LENGTH_SLOPE_PER_CH2,
    COMPOUND_STRUCTURE,
    HOLDOUT_SEALED_BINDING,
    MATRIX_LOADING,
    REVERSIBLE_BINDING,
    UNSATURATION_OBSERVATIONS_EXCLUDED,
)

REPO = Path(__file__).resolve().parents[2]
SHIP_RULE = REPO / "results/validation/kinetic_core_b26_ship_rule.json"


def _rows():
    return {p.key: p for p in REVERSIBLE_BINDING}


def test_the_three_pea_rows_carry_the_printed_partition_pair():
    rows = _rows()
    # K_g = (matrix/gas over buffer/gas - 1) / 10 g/L, the registry's own form.
    assert rows["kg_hexanal_pea"].value == pytest.approx((116.37 / 32.90 - 1) / 10.0, rel=2e-4)
    assert rows["kg_z_2_penten_1_ol_pea"].value == pytest.approx((2065.44 / 1460.49 - 1) / 10.0, rel=2e-3)
    assert rows["kg_t_2_octenal_pea"].value == pytest.approx((2203.85 / 455.93 - 1) / 10.0, rel=2e-4)
    for key in ("kg_hexanal_pea", "kg_z_2_penten_1_ol_pea", "kg_t_2_octenal_pea"):
        assert rows[key].medium == "pea_protein_1pct"
        assert rows[key].method == "static_headspace_partition"
        assert rows[key].ph_of_measurement == 7.6
        assert rows[key].temperature_c == 37.0
        assert "bi2022" in rows[key].source_anchor


def test_the_alkenal_row_is_quarantined_and_pools_into_nothing():
    rows = _rows()
    assert rows["kg_t_2_octenal_pea"].provenance["quarantined_as_binding"] is True
    classes = fit_class_binding_constants()
    assert "alkenal" not in classes, (
        "a quarantined row must not create a class; the alkenal contrast is carried by the "
        "unsaturation penalty, which this wave deliberately does not touch")


def test_the_n_alkanal_class_is_now_two_rows_and_pools_them_geometrically():
    rows = _rows()
    classes = fit_class_binding_constants()
    expected = math.exp((math.log(rows["kg_hexanal_dairy"].value)
                         + math.log(rows["kg_hexanal_pea"].value)) / 2.0)
    assert classes["n_alkanal"]["n_fit_rows"] == 2
    assert classes["n_alkanal"]["k_g_l_per_g"] == pytest.approx(expected, rel=1e-12)
    assert classes["n_alkanal"]["k_g_l_per_g"] == pytest.approx(0.054038, rel=1e-4)
    # the members must be distinguishable: two rows on the same compound, different media
    assert classes["n_alkanal"]["members"] == ["hexanal@skim_milk", "hexanal@pea_protein_1pct"]
    # the chain-length surrogate follows, and only through the measured slope
    assert classes["branched_alkanal"]["k_g_l_per_g"] == pytest.approx(
        expected / CHAIN_LENGTH_SLOPE_PER_CH2, rel=1e-12)
    # the reference loading is still the largest member's, unchanged
    assert classes["n_alkanal"]["reference_loading_g_per_l"] == 33.9


def test_the_alcohol_gets_its_own_class_with_no_panel_consumer():
    classes = fit_class_binding_constants()
    assert classes["alkenol"]["members"] == ["z_2_penten_1_ol@pea_protein_1pct"]
    assert COMPOUND_STRUCTURE["z_2_penten_1_ol"].alpha_beta_unsaturated_carbonyl is False


def test_bi_contrast_is_excluded_from_the_penalty_and_the_penalty_did_not_move():
    penalty = fit_unsaturation_penalty()
    assert penalty["n_fit_rows"] == 2
    assert penalty["penalty_x"] == pytest.approx(math.sqrt(2.81 * 4.95), rel=1e-12)
    assert "unsat_penalty_pea" in UNSATURATION_OBSERVATIONS_EXCLUDED
    reason = UNSATURATION_OBSERVATIONS_EXCLUDED["unsat_penalty_pea"]
    assert "SAME-CARBON" in reason and "chain-length" in reason


def test_the_loading_says_it_is_inherited_across_sections():
    """
    The loading is the paper's own printed number for its binding assay, and the partition
    run these constants come from says the conditions were the same without restating it.
    So it is carried as printed composition -- the standing the skim-milk row has, whose
    composition is cited rather than measured by its own paper -- and the inheritance is
    named in the notes rather than dressed up as a band nobody measured.
    """
    loading = MATRIX_LOADING["pea_protein_1pct"]
    assert loading.protein_g_per_l == 10.0
    assert loading.ph == 7.6
    assert loading.evidence_class == "measured_ratio"
    assert "INHERITED ACROSS SECTIONS" in loading.notes
    assert loading.protein_lo_g_per_l == loading.protein_hi_g_per_l == 10.0
    # and what the inheritance is worth is computed, not asserted
    payload = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    band = [payload["T4"]["loading_double_20_g_per_L"]["hong_hexanal_prediction_x"],
            payload["T4"]["loading_half_5_g_per_L"]["hong_hexanal_prediction_x"]]
    assert band[0] < payload["T4"]["live"]["hong_hexanal_prediction_x"] < band[1]


def test_the_hold_out_rows_improve_and_the_ceiling_breaks():
    """The pre-registered decisive test, re-run from the shipped artifact."""
    payload = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    assert payload["verdict"] == "SHIP"
    assert payload["T2"]["rows_made_worse"] == []
    rows = {r["compound"]: r for r in payload["T2"]["rows"]}
    assert rows["hexanal"]["fold_before"] == pytest.approx(50.3, rel=1e-2)
    assert rows["hexanal"]["fold_after"] == pytest.approx(15.28, rel=1e-2)
    # and the finding the wave pre-registered at 75 %: the term crosses its evidence ceiling
    assert payload["T3"]["rows_now_over_the_ceiling"] == ["hexanal"]
    assert rows["hexanal"]["explained_share_after"] > 0.25
    assert rows["hexanal"]["ceiling_flag"], "the layer must still say so out loud"


def test_the_live_prediction_matches_the_shipped_artifact():
    payload = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    live = predict_matrix_shift("hexanal", "soy_paste_hong").predicted_ratio
    assert live == pytest.approx(payload["T4"]["live"]["hong_hexanal_prediction_x"], rel=1e-12)
    assert live == pytest.approx(8.673, rel=1e-3)


def test_nothing_the_wave_promised_not_to_touch_moved():
    payload = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    assert payload["T5"]["kinetic_panel_identical"] is True
    assert payload["T5"]["sealed_keys_still_valueless"] is True
    carried = {p.key for p in REVERSIBLE_BINDING}
    assert carried.isdisjoint(HOLDOUT_SEALED_BINDING)
    unmoved = set(payload["T1"]["classes_unmoved"])
    assert {"ester", "methyl_ketone", "diketone", "lactone", "furanone"} <= unmoved


def test_the_blind_prediction_file_is_refused_by_default():
    """
    The defect this wave found. `kinetic_core_b4_frozen_predictions.json` records a prediction made
    before its wave read the paired thresholds; re-running B4 after a later wave changed the
    registry would have overwritten it, date and all, with a prediction made by somebody who had
    seen the answer, and nothing afterwards could tell.
    """
    generator = REPO / "scripts/generators/generate_kinetic_core_b4_fit.py"
    source = generator.read_text(encoding="utf-8")
    assert "--refreeze" in source
    assert "if FROZEN_PREDICTIONS.exists() and not args.refreeze:" in source
    frozen = json.loads(
        (REPO / "results/validation/kinetic_core_b4_frozen_predictions.json").read_text())
    # it still holds the pre-B26 constants, which is what makes it the blind record
    assert frozen["class_binding_constants"]["n_alkanal"]["k_g_l_per_g"] == pytest.approx(0.01151)
    assert frozen["class_binding_constants"]["n_alkanal"]["n_fit_rows"] == 1
