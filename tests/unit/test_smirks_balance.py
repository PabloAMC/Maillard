from collections import Counter

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from src.smirks_engine import _SMIRKS_RULES
from src.reaction_templates import _amadori_cascade, _strecker_step, _beta_elimination_steps
from src.chem_utils import Species
from src.conditions import ReactionConditions

def _check_balance(reactants, products):
    """Verify that the total atom counts match between reactants and products."""
    r_atoms = {}
    for r in reactants:
        mol = Chem.MolFromSmiles(r.smiles)
        assert mol is not None, f"unparseable reactant SMILES {r.smiles!r} ({r.label!r})"
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            r_atoms[symbol] = r_atoms.get(symbol, 0) + 1

    p_atoms = {}
    for p in products:
        mol = Chem.MolFromSmiles(p.smiles)
        assert mol is not None, f"unparseable product SMILES {p.smiles!r} ({p.label!r})"
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            p_atoms[symbol] = p_atoms.get(symbol, 0) + 1

    assert r_atoms == p_atoms, f"Balance failure: Reactants {r_atoms} != Products {p_atoms}"


def _template_atom_maps(rxn, side: str) -> Counter:
    count = rxn.GetNumReactantTemplates() if side == "reactant" else rxn.GetNumProductTemplates()
    get = rxn.GetReactantTemplate if side == "reactant" else rxn.GetProductTemplate
    maps: Counter = Counter()
    for index in range(count):
        for atom in get(index).GetAtoms():
            if atom.GetAtomMapNum():
                maps[atom.GetAtomMapNum()] += 1
    return maps


def _unmapped_heavy_atoms(rxn, side: str) -> Counter:
    count = rxn.GetNumReactantTemplates() if side == "reactant" else rxn.GetNumProductTemplates()
    get = rxn.GetReactantTemplate if side == "reactant" else rxn.GetProductTemplate
    symbols: Counter = Counter()
    for index in range(count):
        for atom in get(index).GetAtoms():
            if atom.GetAtomMapNum():
                continue
            symbol = atom.GetSymbol()
            # Explicit [H] on the product side of an H-transfer template is legitimate: the
            # hydrogen it re-materialises was implicit (and therefore unmappable) on the
            # reactant side. Heavy atoms have no such excuse.
            if symbol in ("H", "*"):
                continue
            symbols[symbol] += 1
    return symbols


@pytest.mark.parametrize("name, family, smirks, gate", _SMIRKS_RULES)
def test_smirks_rule_template_conserves_its_mapped_atoms(name, family, smirks, gate):
    """Each SMIRKS template must conserve atom maps and invent no unmapped heavy atoms.

    REWRITTEN 2026-08-27 (Wave J2, red-team finding S7: "no-op balance parametrizations").
    The previous body built ``rxn`` and then executed a bare ``pass`` -- five parametrized
    cases named "rule balance" that asserted nothing at all. Renamed to state exactly what
    is now checked, because it is NOT full stoichiometric balance: balancing a template in
    the abstract needs substrates, and the concrete, substrate-level balance check over the
    whole enumerated network lives in
    ``tests/unit/test_smirks_engine.py::assert_all_balanced`` (329 steps at the time of
    writing).

    What this does establish, per template and without substrates:
      * the SMIRKS parses and RDKit's own validator reports zero errors;
      * every atom map number on one side appears on the other -- a map dropped from the
        product side is the classic way a template deletes an atom;
      * no map number is duplicated in the products, which is how a template silently
        clones an atom;
      * the product side introduces no unmapped HEAVY atom, i.e. no mass from nowhere.
    Verified discriminating: mutating any of the five templates to drop a product map or
    add an unmapped heavy atom fails this test.
    """
    rxn = AllChem.ReactionFromSmarts(smirks)
    assert rxn is not None, f"{name}: SMIRKS does not parse"

    num_warnings, num_errors = rxn.Validate(silent=True)
    assert num_errors == 0, f"{name}: RDKit reports {num_errors} error(s) in the template"

    assert rxn.GetNumReactantTemplates() >= 1, f"{name}: no reactant template"
    assert rxn.GetNumProductTemplates() >= 1, f"{name}: no product template"

    reactant_maps = _template_atom_maps(rxn, "reactant")
    product_maps = _template_atom_maps(rxn, "product")
    assert reactant_maps, f"{name}: template has no mapped atoms, so nothing is traced"

    lost = sorted(set(reactant_maps) - set(product_maps))
    gained = sorted(set(product_maps) - set(reactant_maps))
    assert not lost, f"{name} ({family}): atom map(s) {lost} vanish between reactants and products"
    assert not gained, f"{name} ({family}): atom map(s) {gained} appear only in the products"

    duplicated = sorted(number for number, count in product_maps.items() if count > 1)
    assert not duplicated, f"{name} ({family}): atom map(s) {duplicated} are duplicated in the products"

    invented = _unmapped_heavy_atoms(rxn, "product") - _unmapped_heavy_atoms(rxn, "reactant")
    assert not invented, (
        f"{name} ({family}): product template introduces unmapped heavy atoms {dict(invented)} "
        f"-- mass appearing from nowhere"
    )

def test_amadori_cascade_balance():
    sugar = Species(label="glucose", smiles="OCC(O)C(O)C(O)C(O)C=O")
    aa = Species(label="glycine", smiles="NCC(=O)O")

    steps = _amadori_cascade(sugar, aa)
    # ADDED 2026-08-27 (Wave J2): without this the test passes when the cascade returns
    # nothing, which is precisely the regression it exists to catch.
    assert steps, "the Amadori cascade produced no steps, so nothing was balance-checked"
    for step in steps:
        _check_balance(step.reactants, step.products)

def test_strecker_step_balance():
    dic = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
    aa = Species(label="alanine", smiles="CC(N)C(=O)O")

    step = _strecker_step(dic, aa)
    # TIGHTENED 2026-08-27 (Wave J2): was `if step:`, so a Strecker step that stopped
    # firing turned this test green instead of red. Pyruvaldehyde + alanine is the
    # canonical Strecker substrate pair; if it returns None that is a real regression.
    assert step is not None, "Strecker step did not fire for pyruvaldehyde + alanine"
    _check_balance(step.reactants, step.products)

def test_beta_elimination_balance():
    cys = Species(label="cysteine", smiles="NC(CS)C(=O)O")
    steps = _beta_elimination_steps(cys, [])
    # ADDED 2026-08-27 (Wave J2): same vacuity guard as above.
    assert steps, "beta-elimination produced no steps from cysteine, so nothing was checked"
    for step in steps:
        _check_balance(step.reactants, step.products)
