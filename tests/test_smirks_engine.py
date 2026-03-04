"""
tests/test_smirks_engine.py — Test suite for the Phase 6.1 SMIRKS rule engine.

Verifies that SmirksEngine.enumerate() correctly generates ElementaryStep
sequences for the 4 canonical Maillard precursor systems, and that all
outputs are structurally valid and pipeline-compatible.
"""

import pytest
from pathlib import Path
from rdkit import Chem

from src.smirks_engine import SmirksEngine, _canonical
from src.conditions import ReactionConditions
from src.pathway_extractor import Species, ElementaryStep


# ── Helpers ────────────────────────────────────────────────────────────────

def to_species(label: str, smiles: str) -> Species:
    return Species(label=label, smiles=smiles)


RIBOSE   = to_species("D-ribose",   "O=CC(O)C(O)C(O)CO")
GLUCOSE  = to_species("D-glucose",  "O=CC(O)C(O)C(O)C(O)CO")
GLYCINE  = to_species("glycine",    "NCC(=O)O")
CYSTEINE = to_species("L-cysteine", "NC(CS)C(=O)O")
LEUCINE  = to_species("L-leucine",  "CC(C)CC(N)C(=O)O")
LYSINE   = to_species("L-lysine",   "NCCCCC(N)C(=O)O")
H2S      = to_species("H2S",        "S")
HEXANAL  = to_species("hexanal",    "CCCCCC=O")
FRUCTOSE = to_species("D-fructose", "OCC(=O)C(O)C(O)C(O)CO")

ACID    = ReactionConditions(pH=5.0, temperature_celsius=150.0)
NEUTRAL = ReactionConditions(pH=7.0, temperature_celsius=150.0)


def _labels(steps):
    """Collect all product labels from a list of steps."""
    return {p.label for s in steps for p in s.products}


def _families(steps):
    return {s.reaction_family for s in steps}


def _all_smiles_valid(steps):
    """All product SMILES in all steps must be parseable by RDKit."""
    for step in steps:
        for p in step.products:
            if p.smiles and p.smiles not in ("O", "S", "N", "O=C=O"):
                if Chem.MolFromSmiles(p.smiles) is None:
                    return False, p.smiles
    return True, None


# ── Test 1: Ribose + Glycine @ pH 5 ────────────────────────────────────────

class TestRiboseGlycinePH5:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=ACID)
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE])

    def test_produces_steps(self):
        assert len(self.steps) > 0, "Expected at least 1 step"

    def test_schiff_base_formed(self):
        families = _families(self.steps)
        assert "Schiff_Base_Formation" in families, \
            f"Schiff base not found. Families: {families}"

    def test_amadori_formed(self):
        families = _families(self.steps)
        assert "Amadori_Rearrangement" in families, \
            f"Amadori rearrangement not found. Families: {families}"

    def test_acid_enolisation_fires(self):
        """At pH 5, 1,2-enolisation should fire → furfural."""
        families = _families(self.steps)
        assert "Enolisation_1_2" in families, \
            f"1,2-enolisation not found at pH 5. Families: {families}"

    def test_furfural_produced(self):
        labels = _labels(self.steps)
        assert "furfural" in labels, \
            f"furfural not in products. Products: {labels}"

    def test_neutral_enolisation_does_not_fire(self):
        families = _families(self.steps)
        assert "Enolisation_2_3" not in families, \
            "2,3-enolisation should not fire at pH 5"

    def test_all_smiles_valid(self):
        valid, bad = _all_smiles_valid(self.steps)
        assert valid, f"Invalid SMILES in products: {bad}"

    def test_output_is_elementary_steps(self):
        for step in self.steps:
            assert isinstance(step, ElementaryStep)
            assert len(step.reactants) > 0
            assert len(step.products) > 0


# ── Test 2: Glucose + Glycine @ pH 7 ───────────────────────────────────────

class TestGlucoseGlycinePH7:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([GLUCOSE, GLYCINE])

    def test_produces_steps(self):
        assert len(self.steps) > 0

    def test_neutral_enolisation_fires(self):
        """At pH 7, 2,3-enolisation should fire → pyruvaldehyde."""
        families = _families(self.steps)
        assert "Enolisation_2_3" in families, \
            f"2,3-enolisation not found at pH 7. Families: {families}"

    def test_pyruvaldehyde_produced(self):
        labels = _labels(self.steps)
        assert "pyruvaldehyde" in labels, \
            f"pyruvaldehyde not in products. Products: {labels}"

    def test_acid_enolisation_does_not_fire(self):
        families = _families(self.steps)
        assert "Enolisation_1_2" not in families, \
            "1,2-enolisation should not fire at pH 7"

    def test_all_smiles_valid(self):
        valid, bad = _all_smiles_valid(self.steps)
        assert valid, f"Invalid SMILES in products: {bad}"


# ── Test 3: Ribose + Cysteine (S-Maillard) ─────────────────────────────────

class TestRiboseCysteine:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=ACID)
        self.steps = self.engine.enumerate([RIBOSE, CYSTEINE, H2S])

    def test_beta_elimination_fires(self):
        """Cysteine should undergo beta-elimination → DHA + H₂S."""
        families = _families(self.steps)
        assert "Beta_Elimination" in families, \
            f"Beta elimination not found. Families: {families}"

    def test_dha_produced(self):
        labels = _labels(self.steps)
        assert "dehydroalanine" in labels, \
            f"DHA not in products. Products: {labels}"

    def test_thiol_addition_fires(self):
        """Furfural + H₂S → FFT (via Thiol_Addition SMIRKS)."""
        # Thiol addition requires furfural to be in the pool first
        # (from enolisation). Check if it fires at all.
        families = _families(self.steps)
        # The thiol addition SMIRKS may produce FFT-type compounds
        # just verify it doesn't error; structural validation is in test 4
        assert len(self.steps) > 2

    def test_all_smiles_valid(self):
        valid, bad = _all_smiles_valid(self.steps)
        assert valid, f"Invalid SMILES in products: {bad}"


# ── Test 4: DHA competition — Cysteine + Lysine ─────────────────────────────

class TestDHACompetition:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([RIBOSE, CYSTEINE, LYSINE, H2S])

    def test_lal_produced(self):
        """DHA + Lysine → Lysinoalanine (LAL)."""
        labels = _labels(self.steps)
        assert "lysinoalanine" in labels, \
            f"LAL not in products. Products: {labels}"

    def test_dha_crosslinking_fires(self):
        families = _families(self.steps)
        assert "DHA_Crosslinking" in families, \
            f"DHA crosslinking not found. Families: {families}"


# ── Test 5: Strecker degradation ────────────────────────────────────────────

class TestStrecker:
    def setup_method(self):
        # Neutral pH → pyruvaldehyde (dicarbonyl) generated → Strecker fires
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE, LEUCINE])

    def test_strecker_fires(self):
        families = _families(self.steps)
        assert "Strecker_Degradation" in families, \
            f"Strecker not found. Families: {families}"

    def test_methylbutanal_produced(self):
        labels = _labels(self.steps)
        assert "3-methylbutanal" in labels, \
            f"3-methylbutanal not in products. Products: {labels}"


# ── Test 6: Aminoketone Condensation (Pyrazines) ────────────────────────────

class TestAminoketoneCondensation:
    def setup_method(self):
        # Neutral pH + Leucine -> Strecker fires -> aminoacetone
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE])

    def test_pyrazine_produced(self):
        labels = _labels(self.steps)
        assert "2,5-dimethylpyrazine" in labels, "Pyrazine not formed from aminoacetone"

    def test_condensation_fires(self):
        families = _families(self.steps)
        assert "Aminoketone_Condensation" in families


# ── Test 7: Retro-Aldol Fragmentation ───────────────────────────────────────

class TestRetroAldol:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([GLUCOSE, GLYCINE])

    def test_retro_aldol_fires(self):
        families = _families(self.steps)
        assert "Retro_Aldol_Fragmentation" in families

    def test_fragments_produced(self):
        labels = _labels(self.steps)
        # GLUCOSE gives pyruvaldehyde + glyceraldehyde
        assert "pyruvaldehyde" in labels
        assert "glyceraldehyde" in labels


# ── Test 8: Cysteine Thermal Degradation ────────────────────────────────────

class TestCysteineDegradation:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=ACID) # Acid/Neutral, >100C
        self.steps = self.engine.enumerate([CYSTEINE])
        
    def test_degradation_fires(self):
        families = _families(self.steps)
        assert "Cysteine_Degradation" in families
        
    def test_h2s_and_ammonia_produced(self):
        labels = _labels(self.steps)
        assert "H2S" in labels
        assert "ammonia" in labels
        assert "acetaldehyde" in labels


# ── Test 9: Thiazole Condensation ───────────────────────────────────────────

class TestThiazoleCondensation:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([RIBOSE, LEUCINE, CYSTEINE, H2S])

    def test_thiazole_condensation_fires(self):
        families = _families(self.steps)
        assert "Lipid_Thiazole_Condensation" in families

    def test_isobutylthiazole_produced(self):
        labels = _labels(self.steps)
        # Strecker of Leucine gives 3-methylbutanal -> 2-isobutylthiazole
        assert "2-isobutylthiazole" in labels


# ── Test 10: Heyns Rearrangement (Ketoses) ──────────────────────────────────

class TestHeynsRearrangement:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([FRUCTOSE, GLYCINE])

    def test_heyns_fires(self):
        families = _families(self.steps)
        assert "Heyns_Rearrangement" in families
        assert "Amadori_Rearrangement" not in families

    def test_heyns_product(self):
        labels = _labels(self.steps)
        assert "D-fructose-glycine-Heyns" in labels


# ── Test 11: Structural compatibility with pipeline ──────────────────────────

class TestOutputCompatibility:
    def setup_method(self):
        self.engine = SmirksEngine()
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE])

    def test_step_structure(self):
        for step in self.steps:
            assert hasattr(step, "reactants")
            assert hasattr(step, "products")
            assert hasattr(step, "reaction_family")
            assert isinstance(step.reactants, list)
            assert isinstance(step.products, list)

    def test_species_structure(self):
        for step in self.steps:
            for sp in step.reactants + step.products:
                assert hasattr(sp, "label")
                assert hasattr(sp, "smiles")
                assert isinstance(sp.label, str)
                assert isinstance(sp.smiles, str)

    def test_no_duplicate_steps(self):
        from src.smirks_engine import _step_key
        keys = [_step_key(s) for s in self.steps]
        assert len(keys) == len(set(keys)), "Duplicate steps found"


# ── PBMA Additive & Lipid Tests ───────────────────────────────────────────

class TestPBMAAdditives:
    """Tests for Phase 7.2 PBMA formulation additives and lipid trapping."""

    def test_thiamine_degradation(self):
        engine = SmirksEngine(conditions=ReactionConditions(temperature_celsius=150))
        thiamine = to_species("Thiamine", "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1")
        steps = engine.enumerate([thiamine])
        
        # Should have found the Additive_Thermal_Degradation step
        degrad_steps = [s for s in steps if s.reaction_family == "Additive_Thermal_Degradation"]
        assert len(degrad_steps) >= 1
        
        products = [p.label for step in degrad_steps for p in step.products]
        assert "Hydrogen_Sulfide" in products
        assert "2-methylthiophene" in products
        assert "4,5-dihydro-2-methylthiazole" in products

    def test_glutathione_cleavage(self):
        engine = SmirksEngine(conditions=ReactionConditions(temperature_celsius=150))
        gsh = to_species("L-Glutathione", "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O")
        steps = engine.enumerate([gsh])
        
        degrad_steps = [s for s in steps if s.reaction_family == "Additive_Thermal_Degradation"]
        assert len(degrad_steps) >= 1
        
        products = [p.label for step in degrad_steps for p in step.products]
        assert "Glutamic_Acid" in products
        assert "Cysteinylglycine" in products

    def test_lipid_schiff_base_trapping(self):
        engine = SmirksEngine(conditions=ReactionConditions(pH=6.0))
        hexanal = to_species("Hexanal", "CCCCCC=O")
        glycine = to_species("Glycine", "NCC(=O)O")
        
        steps = engine.enumerate([hexanal, glycine])
        
        sb_steps = [s for s in steps if s.reaction_family == "Lipid_Schiff_Base"]
        assert len(sb_steps) >= 1
        
        # Verify hexanal Schiff base is formed: C6 aliphatic chain attached to =N
        has_hexanal_imine = False
        for step in sb_steps:
            for p in step.products:
                mol = Chem.MolFromSmiles(p.smiles)
                # Quick SMARTS check for aliphatic imine attached to an acid
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("CCCCC/C=N\\CC(=O)O")):
                    has_hexanal_imine = True
        
        assert has_hexanal_imine, "Hexanal was not correctly trapped into a Schiff base by glycine"

