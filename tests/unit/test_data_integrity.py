import pytest
from rdkit import Chem
from data.reactions.curated_pathways import PATHWAYS

def get_comp(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if not mol: return None
        mol = Chem.AddHs(mol)
        comp = {}
        for a in mol.GetAtoms():
            s = a.GetSymbol()
            comp[s] = comp.get(s, 0) + 1
        return comp
    except:
        return None

@pytest.mark.parametrize("pathway_name", PATHWAYS.keys())
def test_curated_pathway_atom_balance(pathway_name):
    """Verify that every step in the curated pathways is atom-balanced."""
    steps = PATHWAYS[pathway_name]
    for i, step in enumerate(steps):
        r_comp = {}
        for s in step.reactants:
            c = get_comp(s.smiles)
            assert c is not None, f"Failed to parse reactant SMILES: {s.smiles}"
            for k, v in c.items():
                r_comp[k] = r_comp.get(k, 0) + v
        
        p_comp = {}
        for s in step.products:
            c = get_comp(s.smiles)
            assert c is not None, f"Failed to parse product SMILES: {s.smiles}"
            for k, v in c.items():
                p_comp[k] = p_comp.get(k, 0) + v
        
        if r_comp != p_comp:
            all_elements = set(r_comp.keys()) | set(p_comp.keys())
            diffs = {}
            for e in all_elements:
                d = p_comp.get(e, 0) - r_comp.get(e, 0)
                if d != 0:
                    diffs[e] = d
            
            error_msg = (
                f"[{pathway_name}] Step {i+1} ({step.reaction_family}) is UNBALANCED:\n"
                f"  Missing in products: { {k: -v for k, v in diffs.items() if v < 0} }\n"
                f"  Extra in products: { {k: v for k, v in diffs.items() if v > 0} }"
            )
            assert r_comp == p_comp, error_msg
