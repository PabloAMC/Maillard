from rdkit import Chem
from rdkit.Chem import AllChem

smarts = "[c:1][CH1:2]=[O:3].[SH2:4].[HH]>>[*:1][CH2:2][SH1:4].[O:3]"
rxn = AllChem.ReactionFromSmarts(smarts)
furfural = Chem.MolFromSmiles("O=Cc1ccco1")
h2s = Chem.MolFromSmiles("S")
h2 = Chem.MolFromSmiles("[HH]")

prods = rxn.RunReactants((furfural, h2s, h2))
for p_tuple in prods:
    print("Products:", [Chem.MolToSmiles(p) for p in p_tuple])
