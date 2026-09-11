"""
docs/guides/REACTION_TREES.md is the appendix a reader reaches for when they want to know what the
model actually integrates. Until 2026-09-11 it drew three trees and never said a fourth path existed,
so a reader finished it believing the model had three. It has four, and the fourth is the one whose
answers most need a structural caveat: its absolute rate is an assumption carried far outside the
span its own source licenses.

These tests hold the correction. They are about the DOCUMENT, not the model: they fail if the fat
path is dropped from the appendix again, or if it is quietly described as a step list.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "docs" / "guides" / "REACTION_TREES.md"
FIG = ROOT / "docs" / "assets" / "thiol_sink" / "14_fat_path.png"
GEN = ROOT / "scripts" / "generators" / "build_reaction_tree.py"


def test_the_appendix_documents_all_four_paths():
    text = DOC.read_text()
    for heading in ("## The sugar path", "## The pentose–cysteine path",
                    "## The acrylamide path", "## The fat path, which is not a tree"):
        assert heading in text, heading


def test_the_appendix_says_four_and_does_not_call_the_fat_path_a_step_list():
    text = DOC.read_text()
    assert "carries **four** paths" in text
    # the old wording claimed every figure was a step list; it must not come back
    assert "The three figures are the actual step lists" not in text
    fat = text.split("## The fat path")[1].split("## Other figures")[0]
    assert "no step list to draw" in fat.lower()
    assert "It declares no reactions" in fat


def test_the_fat_path_figure_exists_and_is_linked():
    assert FIG.exists() and FIG.stat().st_size > 10_000
    assert "14_fat_path.png" in DOC.read_text()


def test_the_two_older_fat_figures_are_no_longer_orphaned():
    """Both were built before the path had a section to sit in and neither was referenced."""
    text = DOC.read_text()
    for name in ("24_fat_path_hexanal.png", "31_lipid_slate_crosscheck.png"):
        assert name in text, name
        assert (ROOT / "docs" / "assets" / "thiol_sink" / name).exists(), name


def test_the_appendix_carries_the_rate_caveat_the_module_itself_raises():
    """The fat lane's own Q10 note warns that its rate is an assumption licensed for 15-40 C. A
    reader of the appendix must meet that fact, because every absolute number on the path rests on
    it."""
    fat = DOC.read_text().split("## The fat path")[1].split("## Other figures")[0].lower()
    assert "15 and 40" in fat
    assert "branch distribution is measured" in fat
    assert "absolute rate is not" in fat


def test_the_generator_no_longer_claims_one_figure_per_lane_and_then_lists_three():
    src = GEN.read_text()
    assert "THE ENGINE HAS FOUR LANES" in src
    assert "def draw_fat_path" in src
    assert '"14_fat_path.png"' in src
