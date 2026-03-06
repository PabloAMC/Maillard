from rdkit import Chem
from rdkit.Chem import AllChem

ribose = Chem.MolFromSmiles('OC[C@H]1OC(O)[C@H](O)[C@@H]1O')
glucose = Chem.MolFromSmiles('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')

# A hemiacetal carbon is connected to a ring oxygen, an OH, an H, and another Carbon.
# We want to break the C-O(ring) bond.
smirks = "[O:1]-[CH1:2](-[OH1:3])-[C:4]>>[OH1:1].[CH1:2](=[O:3])-[C:4]"
rxn = AllChem.ReactionFromSmarts(smirks)

print("Ribose open:")
for prod in rxn.RunReactants((ribose,)):
    print(Chem.MolToSmiles(prod[0]))

print("Glucose open:")
for prod in rxn.RunReactants((glucose,)):
    print(Chem.MolToSmiles(prod[0]))
