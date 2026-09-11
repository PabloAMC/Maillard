"""Wave B31 (2026-09-10): what the isolate arrived with (kinetic_core_b31_prereg.md).

Two halves. `conditions.carried_volatiles` lets a bundle DECLARE the level a pot starts
with, so a formation model is not charged for raw material it never had to make; and a
pot that was never cooked, for which no source prints a starting state, is REFUSED
rather than answered with a formation from zero.
"""
from __future__ import annotations

import json
import math

import pytest

from src import data_paths
from src.kinetic_core import panel
from src.kinetic_core.engine import UNCOOKED_LOOH_CONVERSION_LIMIT, predict
from src.kinetic_core.parameters_lipid import Q10_ASSUMPTION, k_looh_decomp_per_min

BENCH = data_paths.BENCHMARKS_DIR
TRIKUSUMA = BENCH / "pea_isolate_uht_140C_Trikusuma2019.json"
UNCOOKED = (
    BENCH / "pea_isolate_40C_PratapSingh2021.json",
    BENCH / "soy_isolate_40C_PratapSingh2021.json",
    BENCH / "external_validation" / "external_validation_bi_2020_raw_pea_hexanal.json",
    BENCH / "external_validation" / "external_validation_liu_2023_ppi_offnote_baseline.json",
)
#: Three bundles whose vessel ALSO says "no cook", meaning "there is no vessel to record"
#: rather than "nobody heated this": a synthetic snapshot and two commercial products whose
#: conditions block is a proxy operating point. The string alone must never refuse a row.
HOT_NO_COOK = (
    BENCH / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    BENCH / "resconi_2023_pbma_beef_identity_benchmark.json",
    BENCH / "cml_cel_commercial_pbma_Foods2023.json",
)


def _extent(bench_path, q10=None):
    cond = json.loads(bench_path.read_text())["conditions"]
    k = k_looh_decomp_per_min(float(cond["temp_C"]), q10)
    return 1.0 - math.exp(-k * float(cond["time_min"]))


# --------------------------------------------------------------------------
# the declared half
# --------------------------------------------------------------------------

def test_trikusuma_declares_its_own_unheated_column_and_the_engine_reads_it():
    bench = json.loads(TRIKUSUMA.read_text())
    carried = bench["conditions"]["carried_volatiles"]
    # EXTENDED BY WAVE B35 (2026-09-11): the source prints a control column for SEVEN compounds and
    # B31 declared three of them. Declaring some and not the others was arbitrary, so all seven now
    # carry one. The three B31 shipped are asserted by value; the four B35 added come with them.
    assert carried["hexanal"] == 331.0 and carried["2-pentylfuran"] == 59.4 and carried["nonanal"] == 8.24
    assert carried == {"hexanal": 331.0, "2-pentylfuran": 59.4, "nonanal": 8.24,
                       "2,5-dimethylpyrazine": 2.46, "methional": 0.55,
                       "2-acetyl-1-pyrroline": 0.29, "(E,E)-2,4-decadienal": 0.06}
    assert bench["conditions"]["carried_volatiles_anchor"]
    spec = panel.core_spec(bench)
    assert spec.process.carried_volatiles == carried


def test_a_declared_carried_level_adds_and_only_adds():
    """Zero by default: a pot that declares nothing is bit-for-bit what it was."""
    bench = json.loads(TRIKUSUMA.read_text())
    with_declared = panel.core_spec(bench)
    bench["conditions"].pop("carried_volatiles")
    without = panel.core_spec(bench)
    assert without.process.carried_volatiles in (None, {})
    for compound, declared in with_declared.process.carried_volatiles.items():
        # B35 added four more declarations to this pot, three of which the engine REFUSES by name
        # (B18's pyrazine, B22's methional, B24's pyrroline arms all ship inert). A declaration on a
        # refused compound is still correct to record -- it is what the source printed -- and there
        # is simply no answer to compare, so the walk is over the answerable ones.
        run_a, run_b = predict(without, [compound]), predict(with_declared, [compound])
        if not (run_a.answered and run_b.answered):
            continue
        a = run_a.absolutes()[compound]
        b = run_b.absolutes()[compound]
        # The carried level is added BEFORE the matrix-binding factor, so the increase is
        # the declared amount times whatever fraction of it survives binding: strictly
        # positive, and never more than the declared level itself.
        gained = float(b.point_ug_per_l) - float(a.point_ug_per_l)
        assert 0.0 < gained <= declared + 1e-9


def test_only_the_three_declared_rows_moved_and_all_three_improved():
    scores = json.loads((data_paths.VALIDATION_DIR / "core_panel_scores.json").read_text())
    rows = {
        (b["benchmark_id"], c["compound"]): c
        for b in scores["benchmarks"] for c in b["compounds"]
    }
    for compound, fold in (("hexanal", 3.0), ("2-pentylfuran", 3.0), ("nonanal", 3.0)):
        row = rows[("pea_isolate_uht_140C_Trikusuma2019", compound)]
        # All three were outside the 3x band before the declaration and are inside it after.
        assert row["fold_error"] < fold, (compound, row["fold_error"])
        assert row["within_band"]


# --------------------------------------------------------------------------
# the refused half (T3)
# --------------------------------------------------------------------------

def test_the_threshold_sits_in_an_empty_gap_and_survives_the_q10_band():
    """Not a tuned knob: nothing in the panel lies between the two sides."""
    cold = max(_extent(p) for p in UNCOOKED)
    warm = _extent(TRIKUSUMA)  # 140 C for 6 s, the mildest real cook on the panel
    assert cold < UNCOOKED_LOOH_CONVERSION_LIMIT < warm
    assert warm / cold > 50.0
    for q10 in (Q10_ASSUMPTION.lo, Q10_ASSUMPTION.hi):
        assert max(_extent(p, q10) for p in UNCOOKED) < UNCOOKED_LOOH_CONVERSION_LIMIT
        assert _extent(TRIKUSUMA, q10) > UNCOOKED_LOOH_CONVERSION_LIMIT


def test_every_refused_pot_says_in_its_own_provenance_that_nobody_heated_it():
    """The refusal rests on a datum the bundles recorded, not on the size of the miss."""
    for path in UNCOOKED:
        vessel = json.loads(path.read_text())["conditions"]["vessel"]
        assert vessel["closure"] == "no cook"
        note = vessel["provenance_note"].lower()
        assert "never heated" in note or "unheated" in note


@pytest.mark.parametrize("path", UNCOOKED, ids=lambda p: p.stem)
def test_an_uncooked_pot_refuses_its_lipid_rows_and_names_the_cure(path):
    bench = json.loads(path.read_text())
    for compound in panel.bundle_targets(bench):
        run = predict(panel.core_spec(bench), [compound])
        if compound.lower() not in {"hexanal", "2-pentylfuran", "nonanal"}:
            continue
        assert not run.answered, f"{path.stem}/{compound} answered"
        reason = " ".join(run.declaration.reasons)
        assert "NEVER COOKED" in reason
        assert "carried_volatiles" in reason  # the refusal must name its own cure


@pytest.mark.parametrize("path", HOT_NO_COOK, ids=lambda p: p.stem)
def test_the_no_cook_string_alone_refuses_nothing(path):
    """`closure: "no cook"` is overloaded; clause 2 is what keeps these three answerable."""
    bench = json.loads(path.read_text())
    assert bench["conditions"]["vessel"]["closure"] == "no cook"
    assert _extent(path) > UNCOOKED_LOOH_CONVERSION_LIMIT
    for compound in panel.bundle_targets(bench):
        run = predict(panel.core_spec(bench), [compound])
        assert "NEVER COOKED" not in " ".join(run.declaration.reasons)


def test_the_never_cooked_rule_refuses_only_uncooked_pots_and_leaves_every_cooked_miss_standing():
    scores = json.loads((data_paths.VALIDATION_DIR / "core_panel_scores.json").read_text())
    refused = {
        (r["benchmark_id"], r["compound"]) for r in scores["refused_compounds"]
        if "NEVER COOKED" in r["reason"]
    }
    assert refused == {
        ("pea_isolate_40C_PratapSingh2021", "hexanal"),
        ("pea_isolate_40C_PratapSingh2021", "2-pentylfuran"),
        ("soy_isolate_40C_PratapSingh2021", "hexanal"),
        ("soy_isolate_40C_PratapSingh2021", "2-pentylfuran"),
        ("external_validation_bi_2020_raw_pea_hexanal", "hexanal"),
        # B35 (2026-09-11): Bi 2020 prints nonanal for the raw pea too, and the same pot is the same
        # pot -- it was never cooked, so its nonanal is refused on exactly the rule its hexanal was.
        ("external_validation_bi_2020_raw_pea_hexanal", "nonanal"),
        ("external_validation_liu_2023_ppi_offnote_baseline", "hexanal"),
        ("external_validation_liu_2023_ppi_offnote_baseline", "nonanal"),
    }
    # THE GUARD AGAINST A SELF-SERVING RULE. Refusing rows raises the headline by
    # arithmetic alone, so the rule has to be shown NOT to be reaching for the misses:
    # the panel's largest lipid miss is in a pot that WAS cooked, and it survives.
    lipid_misses = [
        (b["benchmark_id"], c["compound"], c["fold_error"])
        for b in scores["benchmarks"] for c in b["compounds"]
        if c["lane"] == "lipid" and c["fold_error"] and c["fold_error"] > 3.0
    ]
    assert lipid_misses, "the rule swept every lipid miss off the panel"
    assert max(f for _, _, f in lipid_misses) > 100.0


def test_a_declared_row_publishes_how_much_of_its_answer_was_declared():
    """The correction to this wave's own claim, frozen so it cannot be quietly dropped.

    Two of the three declared rows are inside the 3x band on about six per cent of their own
    answer. A reader who sees 2.21x must see 93.5 % in the same artifact.
    """
    scores = json.loads((data_paths.VALIDATION_DIR / "core_panel_scores.json").read_text())
    rows = {
        (b["benchmark_id"], c["compound"]): c
        for b in scores["benchmarks"] for c in b["compounds"]
    }
    expected = {  # compound: (declared share, formed-only fold), both to two figures
        "hexanal": (0.935, 19.7),
        "2-pentylfuran": (0.922, 20.6),
        "nonanal": (0.528, 2.14),
    }
    for compound, (share, formed) in expected.items():
        row = rows[("pea_isolate_uht_140C_Trikusuma2019", compound)]
        assert row["carried_declared_ug_per_l"] > 0.0
        assert row["declared_share_of_prediction"] == pytest.approx(share, abs=0.005)
        assert row["fold_error_formed_only"] == pytest.approx(formed, rel=0.02)
        # The whole point: the total flatters, the formed-only does not.
        assert row["fold_error"] < row["fold_error_formed_only"]
    # Only nonanal is a chemistry result on the formed part.
    inside_on_formed = [
        c for c in expected
        if rows[("pea_isolate_uht_140C_Trikusuma2019", c)]["fold_error_formed_only"] <= 3.0
    ]
    assert inside_on_formed == ["nonanal"]


def test_every_other_row_declares_nothing_and_says_so_with_nulls():
    scores = json.loads((data_paths.VALIDATION_DIR / "core_panel_scores.json").read_text())
    declared = [
        (b["benchmark_id"], c["compound"]) for b in scores["benchmarks"] for c in b["compounds"]
        if c.get("carried_declared_ug_per_l") is not None
    ]
    # B35: four, not three -- (E,E)-2,4-decadienal joined them and is answered at 1.62x.
    assert len(declared) == 4 and {b for b, _ in declared} == {"pea_isolate_uht_140C_Trikusuma2019"}
    for b in scores["benchmarks"]:
        for c in b["compounds"]:
            if c["carried_declared_ug_per_l"] is None:
                # Present as an explicit null, not absent: a reader scanning for the split must be
                # able to tell "declares nothing" from "this artifact predates the split".
                for key in ("declared_share_of_prediction", "formed_predicted",
                            "formed_measured", "fold_error_formed_only"):
                    assert key in c and c[key] is None
