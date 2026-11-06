"""The within-study contrast pairs skip a MIXED pair (one absolute row, one relative-only row) and count it.

2026-11-05, H0 follow-up: a study holding the same compound both answered (an absolute ppb level) and refused
(only a centred log level, the thiol lane's ``relative_use_only``) used to raise TypeError in ``engine_dex``
(float minus None). No such pair exists in today's scorecard; this builds one.
"""
from __future__ import annotations

import math

from scripts.generators import build_contrast_scores as C


def _row(benchmark, *, engine=None, centred=None, measured=10.0, lane="trunk"):
    return {"benchmark": benchmark, "study": "study_x", "compound": "methylglyoxal", "lane": lane,
            "measured": measured, "engine": engine, "engine_centred_log10": centred, "out_of_sample": True}


def test_mixed_pair_is_skipped_and_counted():
    rows = [_row("pot_trunk_a", engine=20.0), _row("pot_trunk_b", engine=40.0, measured=30.0),
            _row("pot_sulfur", centred=0.3, lane="sulfur")]
    pairs, mixed = C._pairs(rows)
    assert mixed == 2                                     # the sulfur pot against each trunk pot
    assert [(p["a"], p["b"]) for p in pairs] == [("pot_trunk_a", "pot_trunk_b")]
    assert math.isclose(pairs[0]["engine_dex"], math.log10(20.0 / 40.0))


def test_two_relative_only_rows_still_pair():
    rows = [_row("pot_a", centred=0.25, lane="sulfur"), _row("pot_b", centred=-0.25, lane="sulfur")]
    pairs, mixed = C._pairs(rows)
    assert mixed == 0 and len(pairs) == 1
    assert math.isclose(pairs[0]["engine_dex"], 0.5)


def test_render_reports_the_mixed_count():
    payload = {"definition": "d", "scores": {"all_pairs": {"pairs": 0}}, "by_lane": {}, "mixed_pairs_skipped": 3}
    assert "Mixed pairs skipped" in C.render(payload) and ": 3." in C.render(payload)
