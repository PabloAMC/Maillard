from typing import Optional
from functools import lru_cache

try:
    from rdkit import Chem
except ImportError:
    Chem = None

@lru_cache(maxsize=4096)
def canonicalize_smiles(smi: str, fallback_to_original: bool = False, strip_salts: bool = False) -> Optional[str]:
    """
    Generates a canonical SMILES string using RDKit, stripping stereochemistry
    if necessary, and optionally isolating the largest fragment (strip_salts).
    Returns None if the SMILES is entirely invalid or RDKit is not installed,
    unless fallback_to_original is True.
    """
    if Chem is None or not smi:
        return smi if fallback_to_original else None
        
    try:
        m = Chem.MolFromSmiles(smi)
        if m:
            can_smi = Chem.MolToSmiles(m)
            if strip_salts and "." in can_smi:
                return max(set(can_smi.split(".")), key=len)
            return can_smi
    except Exception:
        pass
    return smi if fallback_to_original else None
