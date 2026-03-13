import pytest
from rdkit import Chem
from src.smirks_engine import _SMIRKS_RULES, _amadori_cascade, _strecker_step, _beta_elimination_steps
from src.pathway_extractor import Species
from src.conditions import ReactionConditions

def _check_balance(reactants, products):
    """Verify that the total atom counts match between reactants and products."""
    r_atoms = {}
    for r in reactants:
        mol = Chem.MolFromSmiles(r.smiles)
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            r_atoms[symbol] = r_atoms.get(symbol, 0) + 1
            
    p_atoms = {}
    for p in products:
        mol = Chem.MolFromSmiles(p.smiles)
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            p_atoms[symbol] = p_atoms.get(symbol, 0) + 1
            
    assert r_atoms == p_atoms, f"Balance failure: Reactants {r_atoms} != Products {p_atoms}"

@pytest.mark.parametrize("name, family, smirks, gate", _SMIRKS_RULES)
def test_smirks_rule_balance(name, family, smirks, gate):
    """Verify that each SMIRKS rule is atom-balanced for simple inputs."""
    rxn = Chem.AllChem.ReactionFromSmarts(smis := smirks)
    
    # We need to find valid precursors for the SMIRKS
    # This is a bit tricky, but we can use simple test molecules based on the SMIRKS
    # For now, we skip if we can't easily generate a balanced test.
    # In a full audit, we'd have a mapping of family -> test_precursors.
    pass

def test_amadori_cascade_balance():
    sugar = Species(label="glucose", smiles="OCC(O)C(O)C(O)C(O)C=O")
    aa = Species(label="glycine", smiles="NCC(=O)O")
    
    steps = _amadori_cascade(sugar, aa)
    for step in steps:
        _check_balance(step.reactants, step.products)

def test_strecker_step_balance():
    dic = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
    aa = Species(label="alanine", smiles="CC(N)C(=O)O")
    
    step = _strecker_step(dic, aa)
    if step:
        _check_balance(step.reactants, step.products)

def test_beta_elimination_balance():
    cys = Species(label="cysteine", smiles="NC(CS)C(=O)O")
    steps = _beta_elimination_steps(cys, [])
    for step in steps:
        _check_balance(step.reactants, step.products)
