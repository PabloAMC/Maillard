from rdkit import Chem
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
print("Reactant Cys:", get_comp("NC(CS)C(=O)O"))
print("Product Pyruv:", get_comp("CC(=O)C=O"))
print("Product S:", get_comp("S"))
print("Product N:", get_comp("N"))
