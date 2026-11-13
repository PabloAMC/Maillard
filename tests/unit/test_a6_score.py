"""Synthetic-data tests of the A6 scorer (scripts/generators/score_a6_post_freeze.py; protocol
docs/reference/lever_table/A6_SCORING_PREREG.md). No keyed data is read: A6 is scored once per batch."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import score_a6_post_freeze as S  # noqa: E402


def _row(lab, y, control, src=None, cell=("X", "pH", "b"), v=0.01):
    return {"compound": cell[0], "lever_class": cell[1], "band": cell[2], "lab": lab, "source": src or lab,
            "contrast": f"{src or lab}-{y}", "y": y, "v": v, "ph_control": control}


def _fit(rows):
    """A stand-in for the table's fit: mean +/- 2 SD of the rows (fast, deterministic)."""
    y = np.array([r["y"] for r in rows])
    sd = max(float(np.std(y)), 0.05)
    return float(np.mean(y) - 2 * sd), float(np.mean(y) + 2 * sd)


def test_interval_score_rewards_coverage_and_narrowness():
    assert S.interval_score(-1, 1, 0) == 2.0
    assert S.interval_score(-1, 1, 2) == 2.0 + 20.0 * 1.0
    assert S.interval_score(-0.5, 0.5, 0) < S.interval_score(-1, 1, 0)


def test_too_few_is_decided_before_any_fit():
    train = [_row(f"L{i}", 0.1 * i, "buffered" if i % 2 else "adjusted_unbuffered") for i in range(6)]
    test = [_row("T1", 0.2, "buffered", src="s1")]

    def boom(rows):
        raise AssertionError("no fit may run on an ineligible batch")

    out = S.score(train, test, fit=boom)
    assert out["verdict"] == "too few to score"
    assert out["eligibility"]["contrasts"] == 1


def test_own_laboratory_is_removed_and_one_class_cells_are_dropped():
    train = [_row("A", 0.1, "buffered"), _row("B", 0.2, "buffered"), _row("C", 0.3, "buffered")]
    test = [_row("A", 0.1, "buffered", src="t")]
    cands, _, dropped = S.candidates(train, test)
    assert cands == [] and dropped == 1          # every training row buffered: R1 = baseline, dropped


def test_a_real_split_wins_and_a_placebo_does_not():
    rng = np.random.default_rng(1)
    train, test = [], []
    for c in range(6):                            # six cells; buffered rows sit near -1, unbuffered near +1
        cell = (f"X{c}", "pH", "b")
        for i in range(8):
            ctl = "buffered" if i < 4 else "adjusted_unbuffered"
            mu = -1.0 if ctl == "buffered" else 1.0
            train.append(_row(f"L{c}{i}", mu + rng.normal(0, 0.1), ctl, cell=cell))
        for j in range(4):                         # test contrasts from 4 sources and 4 labs per cell
            ctl = "buffered" if j % 2 else "adjusted_unbuffered"
            mu = -1.0 if ctl == "buffered" else 1.0
            test.append(_row(f"T{j}", mu + rng.normal(0, 0.1), ctl, src=f"s{j}", cell=cell))
    out = S.score(train, test, fit=_fit)
    assert out["eligibility"]["eligible"]
    assert out["diff"] < 0 and out["verdict"] == "R1 supported"
    assert out["placebo_median"] > out["diff"]


def test_labelling_count_is_the_product_of_binomials():
    rows = {("X", "pH", "b"): [_row("A", 0, "buffered"), _row("B", 0, "adjusted_unbuffered"),
                               _row("C", 0, "adjusted_unbuffered"), _row("D", 0, None)]}
    total, spec = S.labellings(rows, [("X", "pH", "b")])
    assert total == 3 and spec == [([0, 1, 2], 1)]   # the unstated row stays fixed
