"""
WAVE B46: the per-lane offset diagnostic. Pre-registration:
results/validation/kinetic_core_b46_prereg.md.

The artifact moves no constant, so nothing else in the suite would notice if it started lying. These
tests hold the two things the wave's conclusion rests on: that the numbers are what the outcome
section says, and that the honesty machinery added in amendment 1 stays switched on. The second
matters more. The amendment exists because the pre-registered criterion called a two-point comparison
a trend, and a future edit that quietly dropped the design-support columns would restore exactly the
overclaim it was written to stop.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "results" / "validation" / "lane_offset_diagnostic.json"
PREREG = ROOT / "results" / "validation" / "kinetic_core_b46_prereg.md"


def _payload():
    return json.loads(ART.read_text())


def test_every_lane_reports_how_many_distinct_pots_it_rests_on():
    """A lane statement built from one pot is a statement about that pot."""
    for name, lane in _payload()["lanes"].items():
        assert lane["distinct_pots"] >= 1, name
        assert lane["distinct_pots"] <= lane["n"], name


def test_the_trunk_lane_is_flagged_as_a_single_pot():
    p = _payload()
    assert "trunk" in p["summary"]["lanes_whose_rows_are_one_pot"]
    assert p["lanes"]["trunk"]["distinct_pots"] == 1


def test_every_correlation_carries_its_distinct_level_count_and_both_verdicts():
    for name, lane in _payload()["lanes"].items():
        for cov, c in lane["correlations"].items():
            assert "distinct_levels" in c, (name, cov)
            assert "tracks" in c and "tracks_with_design_support" in c, (name, cov)
            # the stricter verdict can never be true where the pre-registered one is false
            assert not (c["tracks_with_design_support"] and not c["tracks"]), (name, cov)
            if c["tracks_with_design_support"]:
                assert c["distinct_levels"] >= 3, (name, cov)


def test_the_fat_lanes_temperature_correlation_is_not_credited_as_a_trend():
    """P3 passed the pre-registered threshold and failed the stricter one. It must stay failed:
    those 8 rows span two temperatures."""
    c = _payload()["lanes"]["lipid"]["correlations"]["temp_C"]
    assert c["tracks"] is True, "the pre-registered verdict is kept, unchanged"
    assert c["tracks_with_design_support"] is False
    assert c["distinct_levels"] == 2
    assert c["why_not"]


def test_the_sulfur_lane_reads_high_which_refuted_the_prediction():
    """P2 said this lane would read LOW. It reads high, and that is the wave's main finding."""
    lane = _payload()["lanes"]["sulfur"]
    assert lane["reads"] == "high"
    assert lane["median_signed_dex"] > 0.5
    assert lane["systematic"] is True
    assert lane["distinct_pots"] >= 10


def test_the_sulfur_temperature_trend_is_the_one_claim_with_design_behind_it():
    supported = _payload()["summary"]["lane_covariate_pairs_tracking_with_design_support"]
    hit = [s for s in supported if s["lane"] == "sulfur" and s["covariate"] == "temp_C"]
    assert hit, "the sulfur temperature trend must survive the design-support filter"
    assert abs(hit[0]["rho"]) >= 0.6
    assert hit[0]["distinct_levels"] >= 6


def test_the_artifact_says_it_is_not_a_fit():
    what = _payload()["what_this_is"].lower()
    assert "not a fit" in what and "no constant moves" in what


def test_the_amendment_is_recorded_as_post_hoc_rather_than_folded_into_the_method():
    text = PREREG.read_text()
    assert "## 6. Amendment 1 — made AFTER seeing the first run" in text
    assert text.index("## 3. Predictions") < text.index("## 5. Outcome")
    assert text.index("## 5. Outcome") < text.index("## 6. Amendment 1")


def test_the_wishlist_derives_the_portfolio_gap_from_this_diagnostic():
    """B46's conclusion was that three of four lanes cannot support a trend statement at all, which
    is a gap in the benchmark portfolio rather than in the model. That must reach the wishlist as a
    DERIVED section, so it cannot go stale while the diagnostic moves."""
    wishlist = json.loads((ROOT / "results" / "validation" / "data_wishlist.json").read_text())
    gaps = wishlist["lane_design_gaps"]
    assert wishlist["summary"]["lanes_that_cannot_support_a_trend"] == len(gaps)
    by_lane = {g["lane"]: g for g in gaps}
    # the trunk's rows are one pot: every trend statement about it is blocked
    assert by_lane["trunk"]["single_pot"] is True
    assert "ONE pot" in by_lane["trunk"]["what_is_blocked"]
    # the fat lane is flat in temperature, which is what withdrew P3
    assert "temp_C" in by_lane["lipid"]["conditions_with_fewer_than_3_levels"]
    # sulfur is the one lane with a design behind it, so it must NOT appear here
    assert "sulfur" not in by_lane
    for g in gaps:
        assert g["measurement"], g["lane"]


def test_the_experiments_guide_records_that_formation_is_now_implicated():
    """Experiment 1 was written across several rounds as a story about removal. The diagnostic says
    the reacting pots read HIGH, and the sink is already pinned. If that correction is dropped, the
    guide sends someone to measure the wrong half."""
    guide = (ROOT / "docs" / "guides" / "EXPERIMENTS.md").read_text()
    exp1 = guide.split("### 1. Where the thiols go")[1].split("### 2.")[0]
    # prose is hard-wrapped, so a phrase can straddle a newline: normalise whitespace before matching
    flat = " ".join(exp1.split()).lower()
    assert "reads high on both meaty thiols, not low" in flat
    assert "the formation side is now equally implicated" in flat
    assert "99.5" in flat
