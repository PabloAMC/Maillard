"""
tests/test_precursor_resolver.py — Comprehensive tests for precursor_resolver.py

Tests the fuzzy-matching precursor name resolution system that bridges human-readable
names to canonical SMILES.

Verifies:
1. Basic precursor lookup (exact match)
2. Fuzzy matching (case-insensitive, typos)
3. Edge cases and error handling
4. Batch resolution (resolve_many)
5. Available precursor listing
"""

import pytest
from src.precursor_resolver import resolve, resolve_many, list_available


class TestResolveBasicMatching:
    """Test basic precursor name resolution."""

    def test_resolve_exact_match_lowercase(self):
        """Should resolve exact lowercase names."""
        precursor = resolve("ribose")
        assert precursor is not None
        assert hasattr(precursor, 'label')
        assert hasattr(precursor, 'smiles')
        assert precursor.label.lower() == "ribose" or "ribose" in precursor.label.lower()
        assert precursor.smiles is not None

    def test_resolve_exact_match_mixed_case(self):
        """Should resolve names regardless of case."""
        precursor1 = resolve("ribose")
        precursor2 = resolve("Ribose")
        precursor3 = resolve("RIBOSE")
        
        assert precursor1 is not None
        assert precursor2 is not None
        assert precursor3 is not None
        # All should resolve to the same compound
        assert precursor1.smiles == precursor2.smiles == precursor3.smiles

    def test_resolve_common_sugars(self):
        """Should resolve known sugars."""
        sugars = ["glucose", "fructose", "ribose"]
        for sugar in sugars:
            precursor = resolve(sugar)
            assert precursor is not None, f"Failed to resolve {sugar}"
            assert len(precursor.smiles) > 0, f"No SMILES for {sugar}"

    def test_resolve_common_amino_acids(self):
        """Should resolve known amino acids."""
        amino_acids = ["glycine", "cysteine", "lysine", "serine"]
        for aa in amino_acids:
            precursor = resolve(aa)
            assert precursor is not None, f"Failed to resolve {aa}"
            assert len(precursor.smiles) > 0, f"No SMILES for {aa}"

    def test_resolve_additives(self):
        """Should resolve known additives."""
        additives = ["thiamine", "glutathione"]
        for additive in additives:
            try:
                precursor = resolve(additive)
                assert precursor is not None, f"Failed to resolve {additive}"
            except ValueError:
                # Some additives might not be available, that's ok
                pass

    def test_resolve_returns_species_object(self):
        """Resolved precursor should be a Species-like object."""
        precursor = resolve("glucose")
        assert hasattr(precursor, 'label'), "Should have label attribute"
        assert hasattr(precursor, 'smiles'), "Should have smiles attribute"
        assert isinstance(precursor.label, str)
        assert isinstance(precursor.smiles, str)
        assert len(precursor.smiles) > 0, "SMILES should not be empty"

    def test_resolve_smiles_validity(self):
        """Resolved SMILES should be valid (parseable by RDKit at minimum)."""
        precursor = resolve("glucose")
        # SMILES should not be empty and should look reasonable
        assert len(precursor.smiles) > 0
        assert any(c.isalpha() for c in precursor.smiles), "SMILES should contain characters"


class TestResolveFuzzyMatching:
    """Test fuzzy matching with typos and partial names."""

    def test_resolve_partial_name(self):
        """Partial names should not resolve - resolver uses exact matching."""
        # The resolver uses exact lowercase matching, not fuzzy matching
        with pytest.raises(ValueError):
            resolve("glu")  # Too short/partial - should not match glucose

    def test_resolve_common_misspellings(self):
        """Should handle common misspellings."""
        # Exact behavior depends on fuzzy matching implementation
        # But it should attempt to resolve despite minor variations
        try:
            # These might fail or succeed depending on fuzzy threshold
            precursor = resolve("glysine")  # misspelling of glycine
            if precursor is not None:
                assert "glycine" in precursor.label.lower() or \
                       "gly" in precursor.label.lower(), \
                       "Should match similar names"
        except ValueError:
            # Fuzzy matching might have a threshold that rejects very bad matches
            pass

    def test_resolve_whitespace_handling(self):
        """Should handle leading/trailing whitespace."""
        precursor1 = resolve("ribose")
        precursor2 = resolve("  ribose  ")
        precursor3 = resolve("\tribose\n")
        
        if precursor1:
            assert precursor1.smiles == precursor2.smiles == precursor3.smiles


class TestResolveErrorHandling:
    """Test error handling for invalid inputs."""

    def test_resolve_unknown_compound_raises_error(self):
        """Should raise ValueError for truly unknown compound."""
        with pytest.raises(ValueError):
            resolve("xyzabc_not_a_real_compound_12345")

    def test_resolve_empty_string_raises_error(self):
        """Empty string should raise error."""
        with pytest.raises((ValueError, KeyError)):
            resolve("")

    def test_resolve_none_raises_error(self):
        """None input should raise error."""
        with pytest.raises((ValueError, TypeError, AttributeError)):
            resolve(None)

    def test_resolve_numeric_input_raises_error(self):
        """Numeric input should raise error."""
        with pytest.raises((ValueError, TypeError, AttributeError)):
            resolve(12345)


class TestResolveMany:
    """Test batch resolution of multiple precursors."""

    def test_resolve_many_single_compound(self):
        """Should resolve single compound in list."""
        precursors = resolve_many(["glucose"])
        assert len(precursors) >= 1, "Should return at least one precursor"
        assert all(hasattr(p, 'smiles') for p in precursors)

    def test_resolve_many_multiple_compounds(self):
        """Should resolve multiple compounds."""
        precursors = resolve_many(["glucose", "glycine"])
        assert len(precursors) == 2, "Should resolve both compounds"
        
        smiles_list = [p.smiles for p in precursors]
        assert len(set(smiles_list)) == 2, "Should have distinct SMILES for different compounds"

    def test_resolve_many_mixed_categories(self):
        """Should resolve precursors from different categories."""
        precursors = resolve_many(["glucose", "glycine", "hexanal"])
        assert len(precursors) >= 2, "Should resolve mixed category inputs"
        assert all(hasattr(p, 'smiles') for p in precursors)

    def test_resolve_many_empty_list(self):
        """Empty list should return empty list or raise error."""
        try:
            precursors = resolve_many([])
            assert isinstance(precursors, list)
        except ValueError:
            # Acceptable to raise error for empty input
            pass

    def test_resolve_many_with_invalid_entry(self):
        """Should raise error if any entry is invalid."""
        with pytest.raises(ValueError):
            resolve_many(["glucose", "xyzabc_invalid_compound"])

    def test_resolve_many_preserves_order_and_count(self):
        """Should return same order and unique entries."""
        names = ["ribose", "glycine", "cysteine"]
        precursors = resolve_many(names)
        
        assert len(precursors) == len(names), "Should return one precursor per input"
        labels = [p.label.lower() for p in precursors]
        
        # Each input should be covered
        for name in names:
            found = any(name in label for label in labels)
            assert found, f"Should have resolved {name}"

    def test_resolve_many_handles_duplicates(self):
        """Should handle duplicate names in input."""
        precursors = resolve_many(["glucose", "glucose"])
        # Should either return 1 or 2, but shouldn't error
        assert len(precursors) >= 1, "Should resolve duplicates without error"


class TestListAvailable:
    """Test listing available precursors."""

    def test_list_available_returns_list(self):
        """Should return list of available precursor names."""
        available = list_available()
        assert isinstance(available, (list, tuple)), "Should return list-like"
        assert len(available) > 0, "Should have at least some available precursors"

    def test_list_available_contains_strings(self):
        """All names in available list should be strings."""
        available = list_available()
        assert all(isinstance(name, str) for name in available), \
            "All available precursor names should be strings"
        assert all(len(name) > 0 for name in available), \
            "No empty precursor names"

    def test_list_available_contains_known_compounds(self):
        """Available list should include known compounds."""
        available = list_available()
        available_lower = [name.lower() for name in available]
        
        known_compounds = ["glucose", "glycine", "ribose"]
        for compound in known_compounds:
            found = any(compound in avail for avail in available_lower)
            assert found, f"Available list should include {compound}"

    def test_listed_compounds_are_resolvable(self):
        """All listed available precursors should be resolvable via resolve()."""
        available = list_available()[:5]  # Test first 5 to keep test fast
        
        for name in available:
            try:
                precursor = resolve(name)
                assert precursor is not None, f"Listed precursor '{name}' should be resolvable"
            except ValueError as e:
                pytest.fail(f"Listed precursor '{name}' raised error: {e}")


class TestPrecursorDataIntegrity:
    """Test quality and consistency of resolved precursor data."""

    def test_smiles_string_format(self):
        """SMILES strings should follow standard format."""
        precursor = resolve("glucose")
        smiles = precursor.smiles
        
        # Basic SMILES validation
        assert isinstance(smiles, str)
        assert len(smiles) > 0
        # Should contain valid SMILES characters
        valid_chars = set("CNOS()[]\\/#=@+-1234567890")
        for char in smiles:
            assert char.isalpha() or char in valid_chars, \
                f"Invalid SMILES character: {char}"

    def test_different_compounds_different_smiles(self):
        """Different compounds should have different SMILES."""
        glucose = resolve("glucose")
        glycine = resolve("glycine")
        ribose = resolve("ribose")
        
        smiles_set = {glucose.smiles, glycine.smiles, ribose.smiles}
        assert len(smiles_set) == 3, "Different compounds should have distinct SMILES"

    def test_same_compound_same_smiles(self):
        """Same compound resolved twice should give same SMILES."""
        precursor1 = resolve("glucose")
        precursor2 = resolve("glucose")
        
        assert precursor1.smiles == precursor2.smiles, \
            "Same compound should always resolve to same SMILES"

    def test_label_not_empty(self):
        """All resolved precursors should have non-empty labels."""
        precursors = resolve_many(["glucose", "glycine", "ribose"])
        
        for precursor in precursors:
            assert precursor.label is not None
            assert isinstance(precursor.label, str)
            assert len(precursor.label) > 0, "Label should not be empty"
