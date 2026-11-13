"""Synthetic-data tests of the prequential A6 scorer (scripts/generators/score_a6_prequential.py; protocol
docs/reference/lever_table/A6_PREQUENTIAL_PREREG.md). No keyed data is read: each A6 record is run once."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import score_a6_prequential as P  # noqa: E402


def _row(lab, y, control, src, cell=("X", "pH", "b"), v=0.01):
    return {"compound": cell[0], "lever_class": cell[1], "band": cell[2], "lab": lab, "source": src,
            "contrast": f"{src}-{lab}-{y:.3f}", "y": y, "v": v, "ph_control": control}


def _fit(rows):
    y = np.array([r["y"] for r in rows])
    sd = max(float(np.std(y)), 0.05)
    return float(np.mean(y) - 2 * sd), float(np.mean(y) + 2 * sd)


def test_a_batch_trains_on_earlier_batches_only():
    # pre-freeze: buffered only. batch 1 adds unbuffered rows from two labs; batch 2's unbuffered contrast can use
    # them, batch 1's cannot (its own batch is not training).
    obs = [_row("A", 0.1, "buffered", "a"), _row("B", 0.2, "buffered", "b"),
           _row("C", 1.0, "adjusted_unbuffered", "c1"), _row("D", 1.1, "unadjusted", "d1"),
           _row("E", 1.2, "adjusted_unbuffered", "e2")]
    bmap = {"c1": 1, "d1": 1, "e2": 2}
    cands, _, _ = P.prequential_candidates(obs, bmap, (1, 2))
    assert [(c["batch"], c["t"]["source"]) for c in cands] == [(2, "e2")]


def test_ineligible_set_runs_no_fit():
    obs = [_row("A", 0.1, "buffered", "a"), _row("B", 0.2, "adjusted_unbuffered", "b"), _row("T", 0, "buffered", "t")]

    def boom(rows):
        raise AssertionError("no fit may run on an ineligible set")

    out = P.score(obs, {"t": 1}, (1,), fit=boom)
    assert out["verdict"] == "too few to score"


def test_a_source_in_two_batches_is_refused():
    with pytest.raises(ValueError):
        P.batch_map({1: ("s",), 2: ("s",)})


def test_a_real_split_wins_prequentially_and_the_placebo_does_not():
    rng = np.random.default_rng(3)
    obs, bmap = [], {}
    for c in range(6):
        cell = (f"X{c}", "pH", "b")
        for i in range(8):                               # pre-freeze: buffered near -1, unbuffered near +1
            ctl = "buffered" if i < 4 else "adjusted_unbuffered"
            obs.append(_row(f"L{c}{i}", (-1 if i < 4 else 1) + rng.normal(0, 0.1), ctl, f"p{c}{i}", cell))
        for k in (1, 2):                                 # two batches, two test sources each
            for j in range(2):
                ctl = "buffered" if j else "unadjusted"
                src = f"s{k}{j}"
                bmap[src] = k
                obs.append(_row(f"T{k}{j}", (-1 if j else 1) + rng.normal(0, 0.1), ctl, src, cell))
    out = P.score(obs, bmap, (1, 2), fit=_fit)
    assert out["eligibility"]["eligible"] and out["by_batch"] == {"1": 12, "2": 12}
    assert out["diff"] < 0 and out["verdict"] == "R1 supported"
    assert out["placebo_median"] > out["diff"]


def test_a_duplicated_dataset_leaves_training():
    obs = [_row("A", 0.1, "buffered", "a"), _row("B", 0.2, "buffered", "b"), _row("C", 0.9, "unadjusted", "c"),
           _row("D", 1.0, "unadjusted", "d"), _row("E", 1.3, "unadjusted", "e"), _row("T", 1.3, "unadjusted", "t")]
    cands, _, dups = P.prequential_candidates(obs, {"t": 1}, (1,))
    assert dups == 1 and {r["source"] for r in cands[0]["rows"]} == {"a", "b", "c", "d"}


def test_one_laboratory_above_40_percent_is_ineligible():
    cands = [{"t": {"source": f"s{i}", "lab": "L" if i < 9 else f"M{i}"}} for i in range(20)]
    e = P.eligibility(cands)
    assert e["largest_lab_share"] == 0.45 and not e["eligible"]


def test_placebo_moves_whole_source_blocks_and_discards_unfittable_draws():
    # one cell: three buffered sources (two rows each) and two unbuffered; the test contrast is unbuffered
    obs = [_row(f"B{i}", 0.1 * i + 0.01 * j, "buffered", f"b{i}") for i in range(3) for j in range(2)]
    obs += [_row(f"U{i}", 1.0 + 0.1 * i, "unadjusted", f"u{i}") for i in range(2)]
    obs += [_row("T", 0.7, "unadjusted", "t")]
    cands, _, _ = P.prequential_candidates(obs, {"t": 1}, (1,))
    total, usable, discarded, distinct = P.placebo_labellings(cands, {"t": 1}, np.random.default_rng(0))
    assert total == 10                       # C(5 blocks, 3 buffered)
    # R1 needs 2 laboratories: a draw making only 2 blocks unbuffered passes; all 10 draws do (5 - 3 = 2)
    assert discarded == 0 and distinct == 10
    for _, subs in usable:
        srcs = [r["source"] for r in subs[0]]
        assert all(srcs.count(s) == (2 if s.startswith("b") else 1) for s in srcs)


def test_a_draw_whose_relabelled_subset_cannot_fit_is_discarded_not_replaced():
    # two buffered sources from ONE laboratory: the draw that labels exactly them unbuffered leaves R1 one laboratory
    obs = [_row("B", 0.1, "buffered", "b0"), _row("B", 0.2, "buffered", "b1"),
           _row("U0", 1.0, "unadjusted", "u0"), _row("U1", 1.1, "unadjusted", "u1"), _row("T", 0.7, "unadjusted", "t")]
    cands, _, _ = P.prequential_candidates(obs, {"t": 1}, (1,))
    total, usable, discarded, distinct = P.placebo_labellings(cands, {"t": 1}, np.random.default_rng(0))
    assert (total, discarded, distinct) == (6, 1, 5)
