from typing import Optional
from functools import lru_cache

try:
    from rdkit import Chem
except ImportError:
    Chem = None

@lru_cache(maxsize=4096)
def _mol_cached(smi: str) -> Optional[Chem.Mol]:
    """Internal cached Mol parsing for reuse."""
    if not smi or Chem is None:
        return None
    return Chem.MolFromSmiles(smi)


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


def parse_mol(smi: str, cloned: bool = True) -> Optional[Chem.Mol]:
    """
    Returns an RDKit Mol object from SMILES. If cloned=True, 
    returns a copy to prevent unintended in-place mutations.
    """
    m = _mol_cached(smi)
    return Chem.Mol(m) if (m and cloned) else m


def calculate_mw(smi: str) -> float:
    """Returns the molecular weight of a SMILES string."""
    from rdkit.Chem import Descriptors
    m = parse_mol(smi, cloned=False)
    return Descriptors.MolWt(m) if m else 9999.0
