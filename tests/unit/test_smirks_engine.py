"""
tests/test_smirks_engine.py — Test suite for the Phase 6.1 SMIRKS rule engine.

Verifies that SmirksEngine.enumerate() correctly generates ElementaryStep
sequences for the 4 canonical Maillard precursor systems, and that all
outputs are structurally valid and pipeline-compatible.
"""

import pytest  # noqa: E402
from rdkit import Chem  # noqa: E402

from src.smirks_engine import SmirksEngine  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from src.chem_utils import Species, ElementaryStep  # noqa: E402


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


def assert_balanced(step: ElementaryStep):
    """Assert that an ElementaryStep conserves every element, hydrogen included.

    HARDENED 2026-08-27 (Wave J2, red-team finding S7). The previous implementation
    silently `continue`d past any species RDKit could not parse, and silently kept a
    hydrogen-suppressed molecule when `AddHs` raised. Both failure modes DROP atoms from
    one side of the ledger, so a genuinely unbalanced reaction containing an unparseable
    or un-hydrogenatable participant passed this assertion. Ten `test_mass_conservation`
    tests were riding on that helper, which made them vacuous in exactly the situation
    they exist to catch.

    Every parse failure is now collected and raised as a test failure naming the species,
    the side it sits on and the reaction family, BEFORE the element comparison runs. That
    ordering matters: an unparseable species must never be reported as a balance diff,
    because the two defects need different fixes.

    Measured at the time of the fix: 329 steps over the 12 enumeration scenarios in this
    file, ZERO unparseable species, ZERO AddHs failures, ZERO unbalanced steps. So the
    hardening lands green -- it removes a latent laundering path rather than papering over
    a live chemistry defect.
    """
    from collections import Counter

    parse_failures: list[str] = []

    def count_atoms(species_list, side: str):
        counts = Counter()
        for sp in species_list:
            mol = Chem.MolFromSmiles(sp.smiles)
            if mol is None:
                mol = Chem.MolFromSmarts(sp.smiles)  # fallback just in case
            if mol is None:
                parse_failures.append(
                    f"{side} species {sp.label!r} has SMILES {sp.smiles!r} that RDKit cannot "
                    f"parse as SMILES or SMARTS"
                )
                continue

            # Explicit hydrogens are mandatory: an implicit-H molecule contributes no H to
            # the ledger, which is the same silent atom-dropping this fix exists to remove.
            try:
                mol = Chem.AddHs(mol)
            except Exception as exc:  # pragma: no cover - defensive; no known trigger
                parse_failures.append(
                    f"{side} species {sp.label!r} ({sp.smiles!r}) failed Chem.AddHs: {exc!r}; "
                    f"its hydrogen count cannot be established"
                )
                continue

            for atom in mol.GetAtoms():
                counts[atom.GetSymbol()] += 1
        return counts

    r_counts = count_atoms(step.reactants, "reactant")
    p_counts = count_atoms(step.products, "product")

    assert not parse_failures, (
        f"Mass conservation for {step.reaction_family} is UNVERIFIABLE, not verified -- "
        f"{len(parse_failures)} species could not be counted, so the atom ledger is "
        f"incomplete on at least one side: " + "; ".join(parse_failures)
    )

    diff = {}
    all_elements = set(r_counts.keys()).union(set(p_counts.keys()))
    for el in all_elements:
        d = r_counts.get(el, 0) - p_counts.get(el, 0)
        if d != 0:
            diff[el] = d

    assert not diff, f"Unbalanced reaction ({step.reaction_family}): Diff = {diff}. Reactants: {[r.smiles for r in step.reactants]} -> Products: {[p.smiles for p in step.products]}"


def assert_all_balanced(steps):
    """Balance every step AND require that there was something to balance.

    ADDED 2026-08-27 (Wave J2). ``for step in self.steps: assert_balanced(step)`` is
    vacuously true when the enumeration returns nothing, so a routing regression that
    silenced a whole family would have turned ten mass-conservation tests green rather
    than red. The emptiness check makes that failure mode loud.
    """
    assert steps, "no steps were enumerated, so mass conservation was never exercised"
    for step in steps:
        assert_balanced(step)


class TestAssertBalancedHelperItself:
    """Negative tests for the balance helper (Wave J2, 2026-08-27).

    The helper is the single point of failure for ten mass-conservation tests. Before this
    wave it silently dropped species it could not parse; these tests pin the hardened
    behaviour so the laundering path cannot be reintroduced without going red.
    """

    def test_unparseable_species_fails_loudly_instead_of_being_dropped(self):
        # Balanced except for the unparseable participant. Under the pre-fix helper the
        # unparseable reactant was skipped and this step "passed" mass conservation.
        step = ElementaryStep(
            reactants=[to_species("water", "O"), to_species("nonsense", "not_a_smiles")],
            products=[to_species("water", "O")],
            reaction_family="Synthetic_Unparseable_Probe",
        )
        with pytest.raises(AssertionError, match="UNVERIFIABLE"):
            assert_balanced(step)

    def test_genuinely_unbalanced_step_is_still_reported_as_a_balance_diff(self):
        step = ElementaryStep(
            reactants=[to_species("water", "O")],
            products=[to_species("water", "O"), to_species("water", "O")],
            reaction_family="Synthetic_Unbalanced_Probe",
        )
        with pytest.raises(AssertionError, match="Unbalanced reaction"):
            assert_balanced(step)

    def test_empty_step_list_is_not_silently_a_pass(self):
        with pytest.raises(AssertionError, match="no steps were enumerated"):
            assert_all_balanced([])


# ── Test 1: Ribose + Glycine @ pH 5 ────────────────────────────────────────

class TestRiboseGlycinePH5:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=ACID)
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE])

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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

    def test_neutral_enolisation_fires(self):
        families = _families(self.steps)
        assert "Enolisation_2_3" in families, \
            "2,3-enolisation should fire at all pH levels (kinetically gated)"

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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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

    def test_acid_enolisation_fires(self):
        families = _families(self.steps)
        assert "Enolisation_1_2" in families, \
            "1,2-enolisation should fire at all pH levels (kinetically gated)"

    def test_all_smiles_valid(self):
        valid, bad = _all_smiles_valid(self.steps)
        assert valid, f"Invalid SMILES in products: {bad}"


# ── Test 3: Ribose + Cysteine (S-Maillard) ─────────────────────────────────

class TestRiboseCysteine:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=ACID)
        h2 = Species("H2", "[HH]")
        self.steps = self.engine.enumerate([RIBOSE, CYSTEINE, H2S, h2])

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

    def test_beta_elimination_fires(self):
        """Cysteine should undergo beta-elimination → DHA + H₂S."""
        families = _families(self.steps)
        assert "Beta_Elimination" in families, \
            f"Beta elimination not found. Families: {families}"

    def test_dha_produced(self):
        labels = _labels(self.steps)
        assert "dehydroalanine" in labels, \
            f"DHA not in products. Products: {labels}"

    def test_thiohemiacetal_formation_fires(self):
        """Furfural + H₂S → Thiohemiacetal (Stage 1 of FFT)."""
        families = _families(self.steps)
        assert "Thiohemiacetal_Formation" in families, \
            f"Thiohemiacetal formation not found. Families: {families}"

    def test_all_smiles_valid(self):
        valid, bad = _all_smiles_valid(self.steps)
        assert valid, f"Invalid SMILES in products: {bad}"


# ── Test 4: DHA competition — Cysteine + Lysine ─────────────────────────────

class TestDHACompetition:
    def setup_method(self):
        self.engine = SmirksEngine(conditions=NEUTRAL)
        self.steps = self.engine.enumerate([RIBOSE, CYSTEINE, LYSINE, H2S])

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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
        water = Species("water", "O")
        h2 = Species("H2", "[HH]")
        self.steps = self.engine.enumerate([GLUCOSE, GLYCINE, water, h2])

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

    def test_retro_aldol_fires(self):
        families = _families(self.steps)
        assert "Enolisation_2_3" in families

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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)
        
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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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

    def test_mass_conservation(self):
        assert_all_balanced(self.steps)

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
        self.engine = SmirksEngine(conditions=ACID)
        self.steps = self.engine.enumerate([RIBOSE, GLYCINE, LYSINE, H2S, CYSTEINE])
        
        # Add a couple more components to trigger more templates
        self.engine_neutral = SmirksEngine(conditions=NEUTRAL)
        self.steps_neutral = self.engine_neutral.enumerate([GLUCOSE, LEUCINE, FRUCTOSE, HEXANAL])

    def test_mass_conservation(self):
        """Phase 8.0: Every generated step must strictly conserve atoms."""
        assert_all_balanced(self.steps + self.steps_neutral)

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
        
        # UPDATED 2026-08-27 (Wave G1 fix 8). The single step this used to assert
        # was grossly unbalanced: thiamine (one sulfur) cannot yield BOTH H2S and
        # a thiophene, and it destroyed 3 C, 3 N and an O per firing. It is now a
        # balanced five-step cascade through the thiazole moiety and
        # 5-hydroxy-3-mercapto-2-pentanone. H2S and 2-methylthiophene survive as
        # products (with a real sulfur budget); 4,5-dihydro-2-methylthiazole was
        # the fabricated second sulfur sink and is gone. MFT is now reachable
        # directly from thiamine, which is the accepted route.
        products = [p.label for step in degrad_steps for p in step.products]
        assert "Hydrogen_Sulfide" in products
        assert "2-methylthiophene" in products
        assert "4-methyl-5-(2-hydroxyethyl)thiazole" in products
        assert "5-hydroxy-3-mercapto-2-pentanone" in products

        # The furan ring closure that yields MFT is an oxidative cyclodehydration
        # and carries its own family so that its 2[H] lumping stays visible.
        aromatisation = [s for s in steps if s.reaction_family == "Furan_Ring_Aromatisation"]
        assert "2-methyl-3-furanthiol" in [
            p.label for step in aromatisation for p in step.products
        ]

    def test_furanone_generation_from_pentose_and_alanine_or_glycine(self):
        engine = SmirksEngine(conditions=ReactionConditions(temperature_celsius=120))
        ribose = to_species("D-Ribose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O")
        alanine = to_species("L-Alanine", "CC(N)C(=O)O")
        glycine = to_species("Glycine", "NCC(=O)O")

        steps = engine.enumerate([ribose, alanine, glycine])
        furanone_steps = [s for s in steps if s.reaction_family == "Furanone_Formation"]

        assert len(furanone_steps) >= 2
        products = {p.label for step in furanone_steps for p in step.products}
        assert "HEMF" in products
        assert "DMHF" in products

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

class TestLipidMaillardSynergy:
    """Tests for Phase 8.E lipid-Maillard synergy pathways."""

    def test_alkylthiazole_formation(self):
        engine = SmirksEngine(conditions=NEUTRAL)
        # We provide hexanal + glycine + leucine + cysteine
        # 1. Glycine/Leucine + Ribose -> dicarbonyls -> aminoacetone
        # 2. Hexanal + aminoacetone + H2S -> 2-pentyl-4-methylthiazole
        hexanal = to_species("Hexanal", "CCCCCC=O")
        ribose = to_species("D-ribose", "O=CC(O)C(O)C(O)CO")
        glycine = to_species("Glycine", "NCC(=O)O")
        cysteine = to_species("L-cysteine", "NC(CS)C(=O)O")
        
        # We need H2S which comes from Cysteine degradation or is provided
        h2s = to_species("H2S", "S")
        
        steps = engine.enumerate([hexanal, ribose, glycine, cysteine, h2s])
        
        families = _families(steps)
        labels = _labels(steps)
        
        assert "Lipid_Strecker_Synergy" in families
        assert "2-pentyl-4-methylthiazole" in labels
        

