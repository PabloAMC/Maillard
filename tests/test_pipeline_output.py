import sys
from pathlib import Path
sys.path.insert(0, str(Path('.')))
from src.conditions import ReactionConditions
from src.smirks_engine import SmirksEngine
from src.precursor_resolver import resolve_many
from src.recommend import Recommender
from rdkit import Chem

precursors = resolve_many(["ribose", "cysteine"])
engine = SmirksEngine(ReactionConditions(pH=5.5, temperature_celsius=150))
steps = engine.enumerate(precursors)

def _canon(smi):
    try:
        can = set(Chem.MolToSmiles(Chem.MolFromSmiles(smi)).split('.'))
        return max(can, key=len)
    except:
        return smi

print("Generated products:")
for s in steps:
    for p in s.products:
        print(f"  {p.label}: {_canon(p.smiles)}")
