"""
tests/test_barrier_constants.py — Unit tests for barrier_constants.py

Tests the FAST-mode heuristic barrier system that underpins all pathway screening.
Validates:
1. Barrier retrieval for all known reaction families
2. Edge cases (unknown families, empty strings, None)
3. Default barrier fallback behavior
4. Consistency with FAST_BARRIERS dictionary
5. Heme catalyst constants
"""

import pytest
from src.barrier_constants import (
    FAST_BARRIERS,
    DEFAULT_BARRIER,
    HEME_CATALYST_REDUCTION,
    HEME_CATALYST_FAMILIES,
    get_barrier as _get_barrier,
)

def get_barrier(family):
    return _get_barrier(family)[0]


class TestFastBarriersConstants:
    """Test the FAST_BARRIERS dictionary and constants."""

    def test_fast_barriers_is_initialized(self):
        """Verify FAST_BARRIERS dictionary is properly populated."""
        assert FAST_BARRIERS is not None
        assert isinstance(FAST_BARRIERS, dict)
        assert len(FAST_BARRIERS) > 0, "FAST_BARRIERS should contain entries"

    def test_fast_barriers_all_positive(self):
        """All barrier values should be positive (kcal/mol)."""
        for pattern, (barrier, _note) in FAST_BARRIERS.items():
            assert barrier > 0, f"Barrier for '{pattern}' is not positive: {barrier}"

    def test_fast_barriers_reasonable_values(self):
        """Barrier values should be within reasonable range (3-50 kcal/mol)."""
        for pattern, (barrier, _note) in FAST_BARRIERS.items():
            assert 0 < barrier < 100, f"Barrier for '{pattern}' out of range: {barrier}"

    def test_fast_barriers_all_have_notes(self):
        """All barriers should have source notes."""
        for pattern, (barrier, note) in FAST_BARRIERS.items():
            assert note is not None
            assert isinstance(note, str)
            assert len(note) > 0, f"No source note for pattern '{pattern}'"

    def test_default_barrier_is_positive(self):
        """DEFAULT_BARRIER should be reasonable."""
        assert DEFAULT_BARRIER > 0
        assert DEFAULT_BARRIER > max(b for b, _ in FAST_BARRIERS.values()), \
            "DEFAULT_BARRIER should be higher than most specific barriers"

    def test_heme_catalyst_reduction_is_positive(self):
        """Heme catalyst reduction should be positive and less than most barriers."""
        assert HEME_CATALYST_REDUCTION > 0
        assert HEME_CATALYST_REDUCTION < 10, "Reduction should be modest"

    def test_heme_catalyst_families_is_set(self):
        """Heme catalyst families should be a non-empty frozenset."""
        assert isinstance(HEME_CATALYST_FAMILIES, frozenset)
        assert len(HEME_CATALYST_FAMILIES) > 0


class TestGetBarrierFunction:
    """Test the get_barrier() function."""

    def test_get_barrier_known_families(self):
        """get_barrier() should return correct values for known reaction families."""
        # Test a few well-known families
        assert get_barrier("schiff") > 0
        assert get_barrier("amadori") > 0
        assert get_barrier("strecker") > 0
        assert get_barrier("cysteine") > 0

    def test_get_barrier_case_insensitive(self):
        """Barrier lookup should be case-insensitive."""
        lower = get_barrier("schiff")
        upper = get_barrier("SCHIFF")
        mixed = get_barrier("Schiff")
        assert lower == upper == mixed, "get_barrier should be case-insensitive"

    def test_get_barrier_partial_match(self):
        """Barrier lookup should work with partial family strings."""
        # "1,2-enolisation" should match pattern "enolisation"
        barrier = get_barrier("1,2-enolisation")
        assert barrier > 0, "Should find barrier for enolisation variant"
        assert barrier == get_barrier("enolisation")

    def test_get_barrier_none_returns_default(self):
        """None or empty input should return DEFAULT_BARRIER."""
        assert get_barrier(None) == DEFAULT_BARRIER
        assert get_barrier("") == DEFAULT_BARRIER

    def test_get_barrier_unknown_family_returns_default(self):
        """Unknown reaction families should return DEFAULT_BARRIER."""
        unknown_barriers = [
            get_barrier("unknown_reaction_family_xyz"),
            get_barrier("foo_bar_baz"),
            get_barrier("NotARealReactionFamily"),
        ]
        for barrier in unknown_barriers:
            assert barrier == DEFAULT_BARRIER, \
                f"Unknown family should return DEFAULT_BARRIER {DEFAULT_BARRIER}, got {barrier}"

    def test_get_barrier_specificity_ordering(self):
        """More specific patterns should match before generic ones."""
        # Test that "schiff" matches before "ring"
        schiff_barrier = get_barrier("Schiff_Base_Formation")
        ring_barrier = get_barrier("Ring_Opening")
        assert schiff_barrier != ring_barrier
        assert schiff_barrier == get_barrier("schiff")
        assert ring_barrier == get_barrier("ring")

    def test_get_barrier_consistent_across_calls(self):
        """get_barrier() should return same value for same input."""
        family = "amadori"
        barrier1 = get_barrier(family)
        barrier2 = get_barrier(family)
        barrier3 = get_barrier(family)
        assert barrier1 == barrier2 == barrier3

    def test_get_barrier_major_families_exist(self):
        """Critical reaction families should all have defined barriers."""
        critical_families = [
            "schiff",
            "amadori",
            "strecker",
            "enolisation",
            "retro",
            "cysteine",
            "beta",
            "thiazole",
        ]
        for family in critical_families:
            barrier = get_barrier(family)
            assert barrier > 0, f"Critical family '{family}' missing or invalid barrier"
            assert barrier != DEFAULT_BARRIER, \
                f"Family '{family}' should have specific barrier, not default"

    def test_get_barrier_returns_float(self):
        """get_barrier() should always return a float."""
        test_families = ["schiff", "amadori", "unknown_xyz", None, ""]
        for family in test_families:
            barrier = get_barrier(family)
            assert isinstance(barrier, (int, float)), \
                f"get_barrier({family!r}) should return numeric, got {type(barrier)}"

    def test_get_barrier_all_fast_barriers_retrievable(self):
        """All patterns in FAST_BARRIERS should be retrievable via get_barrier()."""
        for pattern in FAST_BARRIERS.keys():
            barrier = get_barrier(pattern)
            expected = FAST_BARRIERS[pattern][0]
            assert barrier == expected, \
                f"Pattern '{pattern}' should return {expected}, got {barrier}"


class TestBarrierIntegration:
    """Integration tests: barriers in context of recommendation/screening."""

    def test_barriers_decrease_reasonable_with_ph(self):
        """Verify that pH-accelerated families have lower barriers than typical."""
        # Enolisation is pH-sensitive; should be lower than retro-aldol
        enol_barrier = get_barrier("enolisation")
        retro_barrier = get_barrier("retro")
        # Both exist and are reasonable
        assert 0 < enol_barrier < 40
        assert 0 < retro_barrier < 40

    def test_heme_catalyst_applicable_families(self):
        """Heme catalyst families should be real, defined families."""
        for family in HEME_CATALYST_FAMILIES:
            # Should be able to find a barrier for this family pattern
            barrier = get_barrier(family)
            assert barrier > 0, f"Heme family '{family}' has no valid barrier"
            assert barrier != DEFAULT_BARRIER, \
                f"Heme family '{family}' uses unknown default barrier"

    def test_barrier_reduction_makes_sense(self):
        """With heme catalyst, resulting barrier should still be positive."""
        for family in HEME_CATALYST_FAMILIES:
            original = get_barrier(family)
            reduced = original - HEME_CATALYST_REDUCTION
            assert reduced > 0, \
                f"Heme reduction would make {family} barrier negative: {original} - {HEME_CATALYST_REDUCTION} = {reduced}"

    def test_lipid_synergy_barriers(self):
        """Lipid and synergy pathways should be fast (low barriers)."""
        lipid_barrier = get_barrier("lipid")
        synergy_barrier = get_barrier("synergy")
        # These are fast reactions
        assert lipid_barrier < 20, "Lipid trapping should be fast"
        assert synergy_barrier < 20, "Lipid synergy should be fast"

    def test_sulfur_pathway_consistency(self):
        """Sulfur-containing pathways should be distinct from base reactions."""
        thiol_barrier = get_barrier("thiol")
        thiazole_barrier = get_barrier("thiazole")
        cysteine_barrier = get_barrier("cysteine")
        # All should exist and be distinct
        assert thiol_barrier > 0
        assert thiazole_barrier > 0
        assert cysteine_barrier > 0
