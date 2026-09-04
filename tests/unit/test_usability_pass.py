"""The usability pass (2026-09-04): the verb the docs name exists, a declared-unidentified arm
never yields a magnitude, a zero A arm is undefined, and the data wishlist is built from the
tracked artifacts and says what each measurement would unlock."""
from __future__ import annotations

import subprocess
import sys

import pytest

from src import data_paths


ROOT = data_paths.REPO_ROOT


def _maillard(*args):
    return subprocess.run([sys.executable, "scripts/maillard.py", *args], cwd=str(ROOT),
                          capture_output=True, text=True, timeout=600)


def test_rank_is_accepted_as_well_as_rank_experiments_and_wishlist_is_a_verb():
    help_text = _maillard("--help").stdout
    assert "rank" in help_text and "wishlist" in help_text and "Six verbs" in help_text
    assert _maillard("rank", "--help").returncode == 0
    assert _maillard("rank-experiments", "--help").returncode == 0
    assert _maillard("wishlist", "--help").returncode == 0


def test_a_zero_a_arm_is_undefined_not_higher_in_b():
    from src.kinetic_core.matrix_oav import compare_formulations

    payload = compare_formulations({"x": 0.0, "y": 2.0}, {"x": 5.0, "y": 1.0}, "A", "B")
    rows = {r["compound"]: r for r in payload["rows"]}
    assert rows["x"]["direction"] == "undefined" and rows["x"]["ratio_a_over_b"] == 0.0
    assert rows["y"]["direction"] == "higher_in_A"
    assert payload["n_undefined"] == 1


def test_a_hexose_only_arm_yields_no_ratio_for_the_thiols():
    from src.comparative_cli import spec_to_core
    from src.kinetic_core import engine

    base = {"temp_C": 140.0, "time_min": 30.0, "ph": 5.0, "aw": 0.98}
    ribose = spec_to_core({"name": "ribose", "precursors": {"cysteine": 10.0, "ribose": 10.0}, **base})
    glucose = spec_to_core({"name": "glucose", "precursors": {"cysteine": 10.0, "glucose": 10.0}, **base})
    out = engine.compare(ribose, glucose, ["2-methyl-3-furanthiol", "2-furfurylthiol", "furfural"])
    assert out["comparable"]
    rows = {r["compound"]: r for r in out["ratios"]["rows"]}
    for thiol in ("2-methyl-3-furanthiol", "2-furfurylthiol"):
        assert rows[thiol]["direction"] == "undefined"
        assert rows[thiol]["ratio_a_over_b"] is None
        assert rows[thiol]["unidentified_arm"] == "b"
    assert "unidentified_arm" not in rows["furfural"]
    assert out["ratios"]["n_undefined"] >= 2


def test_the_data_wishlist_builds_from_the_tracked_artifacts_and_names_the_known_gaps():
    from src import data_wishlist

    payload = data_wishlist.build()
    keys = {c["key"] for c in payload["unidentified_coordinates"]}
    assert "k_glc_ha" in keys and "Ea_decay_thiol_sink" in keys
    glc = next(c for c in payload["unidentified_coordinates"] if c["key"] == "k_glc_ha")
    assert glc["reaction"] == "r_glc_c2c3" and "MFT" in glc["unlocks_observables"] and "glucose" in glc["unlocks_from_charges"]
    assert payload["summary"]["not_evaluable_rows"] == len(payload["not_evaluable_rows"]) >= 1
    assert any("2-pentylfuran" in g["what"] for g in payload["refused_targets"])
    aw = next(a for a in payload["thin_axes"] if a["axis"] == "moisture_aw")
    assert aw["structural_block"] and aw["evaluable"] == 0
    assert any(u["if_measured"] == "k_glc_ha" for u in payload["what_you_could_predict"])
    text = data_wishlist.render_markdown(payload)
    assert text.startswith("# Data wishlist") and "## 6. What you could predict" in text


def test_the_reaction_parser_reads_every_sulfur_step_and_shared_constants():
    from src import data_wishlist

    reactions = data_wishlist.sulfur_reactions()
    keys = [r["key"] for r in reactions]
    assert len(reactions) >= 45 and "k_glc_ha" in keys and "k_dimer_mft" in keys  # ch_* names included
    assert keys.count("k_thiol_decay") >= 5
    assert "FFT" in data_wishlist.downstream_observables(reactions, ["FUR"])
    assert "PENT" in data_wishlist.upstream_charges(reactions, ["MFT"])


def test_wishlist_reports_agreeing_claims_needed_to_reach_trust():
    from src.data_wishlist import _claims_to_trust
    from src.directional_reliability import verdict_for, VERDICT_TRUST

    n = _claims_to_trust(4, 5)
    assert n is not None and verdict_for(4 + n, 5 + n) == VERDICT_TRUST and verdict_for(4 + n - 1, 5 + n - 1) != VERDICT_TRUST
