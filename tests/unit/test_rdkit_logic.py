import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from src.smirks_engine import SmirksEngine

def test_smarts_hydration_logic():
    """Verify basic SMARTS reaction for furfural hydration/reduction."""
    smarts = "[c:1][CH1:2]=[O:3].[SH2:4].[HH]>>[*:1][CH2:2][SH1:4].[O:3]"
    rxn = AllChem.ReactionFromSmarts(smarts)
    furfural = Chem.MolFromSmiles("O=Cc1ccco1")
    h2s = Chem.MolFromSmiles("S")
    h2 = Chem.MolFromSmiles("[HH]")

    prods = rxn.RunReactants((furfural, h2s, h2))
    assert len(prods) > 0
    # Check that sulfur was incorporated
    product_smiles = [Chem.MolToSmiles(p[0]) for p in prods]
    assert any("S" in s for s in product_smiles)

def test_ring_opening_smirks():
    """Verify SMIRKS for hemiacetal ring opening (Ribose/Glucose)."""
    # A hemiacetal carbon is connected to a ring oxygen, an OH, an H, and another Carbon.
    smirks = "[O:1]-[CH1:2](-[OH1:3])-[C:4]>>[OH1:1].[CH1:2](=[O:3])-[C:4]"
    rxn = AllChem.ReactionFromSmarts(smirks)
    
    ribose = Chem.MolFromSmiles('OC[C@H]1OC(O)[C@H](O)[C@@H]1O')
    glucose = Chem.MolFromSmiles('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')
    
    # Ribose
    prods_r = rxn.RunReactants((ribose,))
    assert len(prods_r) >= 1
    # Check if ANY of the product segments for the first mapping has the carbonyl
    rib_smiles = [Chem.MolToSmiles(m) for m in prods_r[0]]
    assert any("C=O" in s or "O=C" in s for s in rib_smiles)
    
    # Glucose
    prods_g = rxn.RunReactants((glucose,))
    assert len(prods_g) >= 1
    glc_smiles = [Chem.MolToSmiles(m) for m in prods_g[0]]
    assert any("C=O" in s or "O=C" in s for s in glc_smiles)
