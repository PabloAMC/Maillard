from rdkit import Chem
from src.pathway_extractor import Species


def _parse_mol(smiles: str):
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)

# SMARTS for classifying precursors
_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[CH1]=O")          # aliphatic aldehyde
_KETONE_SMARTS = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]") # internal ketone

def is_sugar(s: Species) -> bool:
    """Heuristic: has an aldehyde OR ketone AND at least 2 hydroxyl groups."""
    m = _parse_mol(s.smiles)
    if m is None:
        return False
    has_aldehyde = m.HasSubstructMatch(_ALDEHYDE_SMARTS)
    has_ketone = m.HasSubstructMatch(_KETONE_SMARTS)
    # Count terminal OH groups
    oh_count = sum(
        1 for atom in m.GetAtoms()
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1
        and atom.GetDegree() == 1
    )
    carbon_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8)
    cyclic_sugar_like = carbon_count in {5, 6} and oxygen_count >= 4 and oh_count >= 2
    return (has_aldehyde or has_ketone or cyclic_sugar_like) and oh_count >= 2

def is_ketose(s: Species) -> bool:
    """Heuristic: has a ketone C=O and multiple OH."""
    m = _parse_mol(s.smiles)
    if not m: 
        return False
    has_ketone = m.HasSubstructMatch(_KETONE_SMARTS)
    oh_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1)
    return has_ketone and oh_count >= 2

def is_hexose(s: Species) -> bool:
    """Heuristic: 6 carbons + is a sugar."""
    m = _parse_mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return is_sugar(s) and c_count == 6

def is_pentose(s: Species) -> bool:
    """Heuristic: 5 carbons + is a sugar."""
    m = _parse_mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return is_sugar(s) and c_count == 5
