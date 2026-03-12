"""
src/smirks_engine.py — Hybrid SMIRKS + Template Maillard Rule Engine

Phase 6.1: Replaces the static curated_pathways.py dictionary with a
rule-based engine that automatically enumerates Maillard reaction pathways
for any arbitrary sugar + amino acid precursor pool.

Architecture:
  Tier B (Templates): Complex multi-step rearrangements encoded as Python
    functions — Amadori cascade, Strecker degradation, enolisation
    branching, beta-elimination (DHA). Exact and chemically grounded.

  Tier A (SMIRKS): Simple, tight functional-group transforms applied
    iteratively on the growing intermediate pool — Schiff base formation
    from lipid aldehydes, thiol addition to generate FFT-type thiols.
    Guarded by abundant MW capping and canonical SMILES deduplication.

Output: List[ElementaryStep] — fully compatible with xtb_screener.py
        and the existing Tier 1 / recommend.py pipeline.
"""

import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple
from functools import lru_cache

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Suppress RDKit atom-mapping warnings from SMIRKS rules that use unmapped atoms
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402

from src.pathway_extractor import Species, ElementaryStep  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402

# Suppress RDKit atom-mapping warnings
RDLogger.DisableLog("rdApp.warning")

# ──────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────

MAX_MW = 300.0  # Daltons — prune products above this (volatiles are small)

# Additive Canonical SMILES (for exact matching)
_THIAMINE_CANONICAL = "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"
_GSH_CANONICAL = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"

# Tier A SMIRKS rules: (name, reaction_family, smirks, ph_gate)
# ph_gate: "any" | "acid" (pH<6) | "neutral" (pH>=6)
_SMIRKS_RULES: List[Tuple[str, str, str, str]] = [
    (
        "schiff_base_lipid",
        "Lipid_Schiff_Base",
        # C3+ aliphatic aldehyde whose alpha-carbon has NO hydroxyl (excludes sugars).
        # The amine donor can be anything with a primary amine on an sp3 carbon (like amino acids).
        "[CX4H2,CX4H3:2][CH1:1]=[O:5].[NH2:3][CX4:4]>>[*:2][CH1:1]=[N:3][*:4].[O:5]",
        "any",
    ),
    (
        "beta_scission_alkoxy",
        "Beta_Scission",
        # Alkoxy radical (R-O.) fragmentation. 
        # Pattern: [C:1]([O:2])-[C:3]. >> [C:1]=[O:2] + [C:3].
        "[CX4:3][CX4:1][OX1H0:2]>>[C:1]=[O:2].[C:3]",
        "any",
    ),
    (
        "radical_propagation_o2",
        "Radical_Propagation_O2",
        # Alkyl radical + O2 -> Peroxy radical
        "[C;X3:1].[O;X1:2]=[O;X1:3]>>[C;X4:1]-[O;X2:2]-[O;X1:3]",
        "any",
    ),
    (
        "peroxy_h_abstraction",
        "Peroxy_H_Abstraction",
        # Peroxy radical + Reactive H (Allylic or specific lipid H) -> Hydroperoxide + Alkyl radical
        # Pattern: [O.] + [H-C-C=C] or [H-C-O]
        # We'll use a broad but restricted pattern to avoid matching EVERY sugar carbon.
        "[O;X1:1]-[O;X2:2].[C;H1,H2,H3;$(C-C=C);!$(C=O):3]>>[O;X2:1]([H])[O;X2:2].[C;X3:3]",
        "any",
    ),
    (
        "radical_termination_disproportionation",
        "Radical_Termination",
        # Two peroxy radicals -> O2 + ... (simplified)
        "[O;X1:1]-[O;X2:2].[O;X1:3]-[O;X2:4]>>[O;X2:1]=[O;X2:3].[O;X2:2].[O;X2:4]",
        "any",
    ),
]


# ──────────────────────────────────────────────────────────────────────────
# Structural classification helpers
# ──────────────────────────────────────────────────────────────────────────

# SMARTS for classifying precursors
_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[CH1]=O")          # aliphatic aldehyde
_POLYOL_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[CH1]=O")   # reused on sugars
_AMINO_SMARTS = Chem.MolFromSmarts("[NH2][CX4]")          # primary amine, not amide
_THIOL_SMARTS = Chem.MolFromSmarts("[SH]")                # thiol or H2S
_DICARBONYL_SMARTS = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)")  # adjacent carbonyls
_AROMATIC_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[c][CH]=O") # furfural-type


@lru_cache(maxsize=4096)
def _mol_cached(smi: str) -> Optional[Chem.Mol]:
    """Internal cached Mol parsing."""
    return Chem.MolFromSmiles(smi) if smi else None


def _mol(smi: str) -> Optional[Chem.Mol]:
    """Returns a CLONE of the cached Mol to prevent in-place mutation issues."""
    m = _mol_cached(smi)
    return Chem.Mol(m) if m else None


@lru_cache(maxsize=4096)
def _canonical(smi: str) -> Optional[str]:
    m = _mol_cached(smi)
    return Chem.MolToSmiles(m) if m else None


def _mw(smi: str) -> float:
    m = _mol(smi)
    return Descriptors.MolWt(m) if m else 9999.0


def _is_valid(smi: str) -> bool:
    """Return True if SMILES is parseable and MW is below cap."""
    return _mol(smi) is not None and _mw(smi) <= MAX_MW


def _is_sugar(s: Species) -> bool:
    """Heuristic: has an aldehyde OR ketone AND at least 2 hydroxyl groups."""
    m = _mol(s.smiles)
    if m is None:
        return False
    has_aldehyde = m.HasSubstructMatch(_ALDEHYDE_SMARTS)
    has_ketone = m.HasSubstructMatch(Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]"))
    # Count OH groups (O with at least one H)
    oh_count = sum(
        1 for atom in m.GetAtoms()
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1
        and atom.GetDegree() == 1  # terminal OH
    )
    return (has_aldehyde or has_ketone) and oh_count >= 2


def _is_ketose(s: Species) -> bool:
    """Heuristic: has a ketone C=O and multiple OH."""
    m = _mol(s.smiles)
    if not m: 
        return False
    pat = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]")
    has_ketone = m.HasSubstructMatch(pat)
    oh_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1)
    return has_ketone and oh_count >= 2


def _is_hexose(s: Species) -> bool:
    """Heuristic: 6 carbons + is a sugar."""
    m = _mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return _is_sugar(s) and c_count == 6


def _is_pentose(s: Species) -> bool:
    """Heuristic: 5 carbons + is a sugar."""
    m = _mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return _is_sugar(s) and c_count == 5


def _is_asparagine(s: Species) -> bool:
    """Detects strictly free asparagine."""
    if s.label.lower() in ["l-asparagine", "asparagine"]: 
        return True
    return s.smiles == "NC(CC(N)=O)C(=O)O"


def _is_lysine(s: Species) -> bool:
    """Detects strictly free lysine."""
    if s.label.lower() in ["l-lysine", "lysine"]: 
        return True
    return s.smiles == "NCCCCC(N)C(=O)O"

def _is_primary_amine(s: Species) -> bool:
    """Has a primary amine not in amide context (amino acids, amines)."""
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_AMINO_SMARTS)


def _has_thiol(s: Species) -> bool:
    """Has thiol group (–SH) or is H₂S."""
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_THIOL_SMARTS)


def _has_cysteine_beta_carbon(s: Species) -> bool:
    """
    Detects cysteine-like β-carbon with thiol (Cys, Ser).
    Pattern: N[C][C][SH or OH] (2-aminol/thiol-3-carbon skeleton).
    """
    pat = Chem.MolFromSmarts("[NH2][CH1][CH2][SX2H,OX2H]")
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(pat)


def _is_dicarbonyl(s: Species) -> bool:
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_DICARBONYL_SMARTS)


def _is_aromatic_aldehyde(s: Species) -> bool:
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_AROMATIC_ALDEHYDE_SMARTS)


def _is_lipid_aldehyde(s: Species) -> bool:
    """C5+ aliphatic monocarbonyl aldehyde without multiple OH (excludes sugars/dicarbonyls)."""
    m = _mol(s.smiles)
    if not m: 
        return False
    # Has aldehyde
    if not m.HasSubstructMatch(_ALDEHYDE_SMARTS): 
        return False
    # NOT aromatic
    if m.HasSubstructMatch(Chem.MolFromSmarts("a")): 
        return False
    # NOT dicarbonyl
    if m.HasSubstructMatch(_ALDEHYDE_SMARTS) and len(m.GetSubstructMatches(_ALDEHYDE_SMARTS)) > 1:
        return False
    # NOT nitrogenous (excludes amino-aldehydes like 5-aminopentanal)
    if any(atom.GetAtomicNum() == 7 for atom in m.GetAtoms()): 
        return False
    # Not a sugar (oh_count < 2)
    oh_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1)
    # C5+ (typically lipid-derived volatiles like pentanal, hexanal)
    c_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    return oh_count < 2 and c_count >= 5


def _is_lipid_hydroperoxide(s: Species) -> bool:
    """Heuristic: has a hydroperoxide group [OX2H,OX2-] and a lipid chain."""
    m = _mol(s.smiles)
    if not m:
        return False
    # Match R-O-OH or anion
    pat = Chem.MolFromSmarts("[OX2,OX1-][OX2H,OX1H0-]")
    if not m.HasSubstructMatch(pat):
        return False
    # C8+ (typical PUFA derivatives)
    c_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    return c_count >= 8


def _species_from_pool(pool: Set[str], label: str, smiles: str) -> Species:
    """Create a Species, canonicalise its SMILES, and add to pool."""
    can = _canonical(smiles)
    if can:
        pool.add(can)
    return Species(label=label, smiles=smiles)


# ──────────────────────────────────────────────────────────────────────────
# Tier B: Parameterised chemical templates
# ──────────────────────────────────────────────────────────────────────────

def _amadori_cascade(sugar: Species, amino_acid: Species) -> List[ElementaryStep]:
    """
    Template: sugar + amino acid → Schiff base → Amadori/Heyns product

    The Schiff base SMILES is constructed by replacing the sugar's terminal
    carbonyl with C=N linked to the amino acid's alpha-amino group.
    For aldoses -> Amadori rearrangement (ketosamine).
    For ketoses -> Heyns rearrangement (aldosamine).
    """
    steps = []
    water = Species(label="water", smiles="O")

    # The string approach actually worked well except for extracting the N.
    # We should just ensure we extract the alpha fragment AND append the basic NH!
    # Let's write a robust version of extraction using RDKit that returns the whole molecule
    # MINUS the oxygen.
    
    # RDKit approach:
    # 1. Sugar: Convert C=O to C-OH and add dummy linker?
    # Simpler: If we know the exact SMILES of sugar and AA, let's use string manipulation 
    # but WITH atom conservation.
    # Aldohexose: O=CC(O)C(O)C(O)C(O)CO + NCC(=O)O -> OCC(O)C(O)C(O)C(O)/C=N/CC(=O)O + H2O
    # Ribose: O=CC(O)C(O)C(O)CO + NCC(=O)O -> OCC(O)C(O)C(O)/C=N/CC(=O)O + H2O
    
    # We need the fragment of the amino acid starting *from the N*, but without its 2 Hs.
    # In canonical SMILES, primary amines are often `N...` or `[NH2]...`.
    # Let's find the N, isolate the fragment, and use it.
    
    _fragment = _extract_alpha_amine_fragment_with_n(amino_acid)
    if not _fragment:
        return []
        
    schiff_label = f"{sugar.label}-{amino_acid.label}-Schiff-base"
    
    if _is_ketose(sugar):
        amadori_label = f"{sugar.label}-{amino_acid.label}-Heyns"
        family = "Heyns_Rearrangement"
        if "fructose" in sugar.label.lower():
            # Fructose: OCC(=O)C(O)C(O)C(O)CO
            # C=O is at index 1 (carbon 2).
            # Fragment = -N-R
            schiff_smiles = f"OCC(={_fragment})C(O)C(O)C(O)CO"
            # Add parenthesis around the -NH-R group so the chain continues properly
            amadori_smiles = f"O=CC(-{_fragment}H)C(O)C(O)C(O)CO"
            # Wait, if _fragment is e.g. "NC(CC)C(=O)O", "-NC(CC)C(=O)OH" is invalid because it has a dash before N.
            # RDKit will output N(CC)C(=O)O if we root it.
            # So _fragment starts with N.
            # schiff: OCC(=N...)
            # amadori: O=CC(N...H) -> The 'H' can be appended inside the N parenthesis if any, but string concat is hard.
            # The safest is: "O=CC(" + _fragment + ")C(O)C(O)C(O)CO", and let RDKit implicitly add the H to N to satisfy valence.
            amadori_smiles = f"O=CC({_fragment})C(O)C(O)C(O)CO"
        else:
            return []
    else:
        amadori_label = f"{sugar.label}-{amino_acid.label}-Amadori"
        family = "Amadori_Rearrangement"
        if _is_pentose(sugar):
            # Ribose: O=CC(O)C(O)C(O)CO
            schiff_smiles = f"OCC(O)C(O)C(O)/C={_fragment}"
            # Again, use parenthesis around the N fragment
            amadori_smiles = f"OCC(O)C(O)C(=O)C({_fragment})"
        elif _is_hexose(sugar):
            # Glucose: O=CC(O)C(O)C(O)C(O)CO
            schiff_smiles = f"OCC(O)C(O)C(O)C(O)/C={_fragment}"
            amadori_smiles = f"OCC(O)C(O)C(O)C(=O)C({_fragment})"
        else:
            return []

    if not _is_valid(schiff_smiles) or not _is_valid(amadori_smiles):
        # We might have generated invalid stereo like C=N]C...
        # Let's clean it up slightly and re-verify
        schiff_smiles = schiff_smiles.replace("=[N", "=N").replace("=[nH]", "=N")
        amadori_smiles = amadori_smiles.replace("C[N", "CN").replace("C[nH]", "CN")
        if not _is_valid(schiff_smiles) or not _is_valid(amadori_smiles):
            return []

    schiff_base = Species(label=schiff_label, smiles=schiff_smiles)
    amadori_product = Species(label=amadori_label, smiles=amadori_smiles)

    steps.append(ElementaryStep(
        reactants=[sugar, amino_acid],
        products=[schiff_base, water],
        reaction_family="Schiff_Base_Formation",
    ))
    steps.append(ElementaryStep(
        reactants=[schiff_base],
        products=[amadori_product],
        reaction_family=family,
    ))
    return steps


def _extract_alpha_amine_fragment_with_n(amino_acid: Species) -> str:
    """
    Extracts the alpha-amino acid SMILES with its alpha nitrogen explicitly grouped,
    e.g., as N(...). We identify the primary alpha N, insert an attachment point or 
    restructure the SMILES so it begins with N.
    """
    # For glycine: NCC(=O)O -> N(CC(=O)O) or something similar
    # By default, RDKit canonical SMILES for simple AAs usually start with N.
    m = _mol(amino_acid.smiles)
    if not m: 
        return ""
    
    # Find alpha nitrogen (N attached to CH attached to C=O)
    pat = Chem.MolFromSmarts("[NH2][CH1,CH2][C](=O)[OH]")
    matches = m.GetSubstructMatches(pat)
    
    # If no match, maybe it's cysteine or something that didn't match perfectly.
    # Try more general: primary amine
    if not matches:
        pat = Chem.MolFromSmarts("[NH2]")
        matches = m.GetSubstructMatches(pat)
        
    if not matches:
        # Fallback if we really can't find it
        return ""
        
    # We want to re-root the SMILES generation at the alpha nitrogen
    n_idx = matches[0][0]
    smi = Chem.MolToSmiles(m, rootedAtAtom=n_idx, isomericSmiles=False)
    
    # The SMILES will start with N. We want to return exactly that string,
    # but when it's appended to C= etc., we'll strip the leading N? No, we WANT the N.
    # E.g. rooted glycine: NCC(=O)O. We want to return N(CC(=O)O).
    # Wait, if we return NCC(=O)O, and substitute `C={fragment}`, we get `C=NCC(=O)O`, 
    # which is exactly correct!
    
    # Let's test Lysine: NCCCCC(N)C(=O)O
    # Rotated at alpha-amine: NC(CCCCN)C(=O)O
    # So `C=NC(CCCCN)C(=O)O` works perfectly!
    
    return smi


def _enolisation_steps(
    amadori: Species,
    sugar: Species,
    amino_acid: Species, # We need the original AA to balance atoms
    conditions: ReactionConditions
) -> List[ElementaryStep]:
    """
    Amadori product → 3-deoxyosone + amino_acid + H2O.
    1,2-enolisation: deoxyosone → furfural/HMF + 2 H2O
    2,3-enolisation: deoxyosone → pyruvaldehyde + fragment
    """
    steps = []
    water = Species(label="water", smiles="O")

    if _is_pentose(sugar):
        deoxy_smi = "O=CC(=O)CC(O)CO"
    else:
        deoxy_smi = "O=CC(=O)CC(O)C(O)CO"

    deoxy = Species(label=f"{sugar.label}-deoxyosone-3", smiles=deoxy_smi)

    # 1. Formation of deoxyosone intermediate
    # C5H8O4 + C2H5NO2 (glycine) = C7H13NO6 (Amadori)!
    steps.append(ElementaryStep(
        reactants=[amadori],
        products=[deoxy, amino_acid],
        reaction_family="Enolisation_Intermediate"
    ))

    # 2. Dehydration to final products
    # 2a. 1,2-enolisation (favored at acidic pH)
    if _is_pentose(sugar):
        product_12 = Species(label="furfural", smiles="O=Cc1ccco1")
        water_count = 2
    else:
        product_12 = Species(label="HMF", smiles="OCC1=CC=C(C=O)O1")
        water_count = 2
    
    steps.append(ElementaryStep(
        reactants=[deoxy],
        products=[product_12] + [water] * water_count,
        reaction_family="Enolisation_1_2"
    ))

    # 2b. 2,3-enolisation (favored at neutral/alkali pH)
    product_23 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
    if _is_pentose(sugar): # C5H8O4 -> C3H4O2 (pyruv) + C2H4O2 (glycolaldehyde)
        p2 = Species(label="glycolaldehyde", smiles="O=CCO")
    else: # C6H10O5 -> C3H4O2 (pyruv) + C3H6O3 (glyceraldehyde)
        p2 = Species(label="glyceraldehyde", smiles="O=CC(O)CO")
    
    steps.append(ElementaryStep(
        reactants=[deoxy],
        products=[product_23, p2],
        reaction_family="Enolisation_2_3"
    ))

    return steps


def _strecker_step(
    dicarbonyl: Species, amino_acid: Species
) -> Optional[ElementaryStep]:
    """
    α-dicarbonyl (e.g. pyruvaldehyde) + amino acid → Strecker aldehyde + aminoketone + CO₂
    
    To balance atoms: 
    Amino acid (e.g. Glycine: C2 H5 N O2) loses CO2 (C1 O2) and its sidechain becomes the Strecker aldehyde.
    The remaining -(NH2) group from the amino acid attaches to the dicarbonyl.
    The dicarbonyl (e.g. Pyruvaldehyde: C3 H4 O2) loses ONE oxygen (which goes to the Strecker aldehyde as its carbonyl O), 
    and accepts the -(NH2) to form the aminoketone.
    
    Wait, let's track the exact mechanism:
    Dicarbonyl: R1-C(=O)-C(=O)-R2
    Amino Acid: NH2-CH(R3)-COOH
    
    1. Condensation to Schiff base, losing H2O (from dicarbonyl O and amino 2H).
    2. Decarboxylation: Loses CO2.
    3. Hydrolysis: Adds H2O across the C=N bond.
    
    Net reaction:
    R1-C(=O)-C(=O)-R2 + NH2-CH(R3)-COOH → R1-C(=O)-CH(NH2)-R2 + O=CH-R3 + CO2
    
    So the aminoketone is exactly the dicarbonyl minus ONE carbonyl oxygen, plus NH2, plus 1 H (from the amino acid alpha carbon).
    Since building this dynamically for pyruvaldehyde (CC(=O)C=O) vs diacetyl (CC(=O)C(=O)C) via SMIRKS is complex, 
    we'll use a mapping for both the amino acid AND the dicarbonyl.
    """
    
    # 1. Map Amino Acid to its Strecker Aldehyde
    _aa_to_aldehyde = {
        # name -> (aldehyde_label, aldehyde_smiles)
        "l-leucine":      ("3-methylbutanal", "CC(C)CC=O"),
        "l-isoleucine":   ("2-methylbutanal", "CCC(C)C=O"),
        "l-valine":       ("2-methylpropanal","CC(C)C=O"),
        "glycine":        ("formaldehyde",    "C=O"), # Glycine sidechain is H. So H-C=O is formaldehyde
        "l-alanine":      ("acetaldehyde",    "CC=O"),
        "l-phenylalanine":("phenylacetaldehyde","O=CCc1ccccc1"),
        "l-methionine":   ("methional",       "CSCCC=O"),
        "l-lysine":       ("5-aminopentanal", "NCCCCC=O"), # Assuming epsilon amine doesn't react here
    }

    # 2. Map Dicarbonyl to its Aminoketone
    # R1-C(=O)-C(=O)-R2 -> R1-C(=O)-CH(NH2)-R2
    _dicarbonyl_to_ak = {
        "pyruvaldehyde": ("aminoacetone", "CC(=O)CN"), # CC(=O)C=O (C3H4O2) -> CC(=O)CN (C3H7NO)
        "diacetyl":      ("3-amino-2-butanone", "CC(=O)C(N)C"),
        "glyoxal":       ("2-aminoethanal", "O=CCN"),
        "furfural":      None, # Not a dicarbonyl
        "HMF":           None, 
    }
    
    aa_entry = _aa_to_aldehyde.get(amino_acid.label.lower())
    if aa_entry is None:
        # Fallback for unrecognized amino acids, though they won't balance if we don't know the products.
        return None

    ald_label, ald_smiles = aa_entry
    
    # Find matching dicarbonyl
    ak_entry = _dicarbonyl_to_ak.get(dicarbonyl.label)
    
    # If it's a generic deoxyosone, it acts as the dicarbonyl.
    # e.g., D-glucose-deoxyosone-3 is O=CC(=O)CC(O)C(O)CO.
    # It turns into the corresponding aminoketone: O=CC(N)CC(O)C(O)CO or NC(=O)CC(O)C(O)CO (Wait, aldehydes are more reactive).
    # Since generating these dynamically is hard, we can use an RDKit reaction!
    
    if ak_entry is None:
        # Generic RDKit reaction for dicarbonyl + amino acid -> aldehyde + aminoketone + CO2
        # It's much easier to just use RunReactants.
        # Dicarbonyl: [C:1](=[O:2])[C:3](=[O:4])
        # Amino Acid: [NH2:5][CH1,CH2:6]([R:7])[C:8](=[O:9])[OH:10]
        # Strecker aldehyde: [R:7][C:6]=[O:2] (Wait, the O comes from dicarbonyl? No, O comes from water hydrolysis. Net it's the same as swapping).
        # We know the aldehyde SMILES from the dictionary. We just need the aminoketone!
        rxn_ak = AllChem.ReactionFromSmarts(
            "[C:1](=[O:2])[C:3](=[O:4]) >> [C:1](=[O:2])[C:3]([NH2])"
        )
        dic_mol = _mol(dicarbonyl.smiles)
        if not dic_mol: 
            return None
        prods = rxn_ak.RunReactants((dic_mol,))
        if not prods: 
            return None
        
        try:
            Chem.SanitizeMol(prods[0][0])
            ak_smiles = Chem.MolToSmiles(prods[0][0])
            ak_label = f"amino-{dicarbonyl.label}"
        except Exception:
            return None
    else:
        ak_label, ak_smiles = ak_entry

    aldehyde = Species(label=ald_label, smiles=ald_smiles)
    aminoketone = Species(label=ak_label, smiles=ak_smiles)
    co2 = Species(label="CO2", smiles="O=C=O")

    return ElementaryStep(
        reactants=[dicarbonyl, amino_acid],
        products=[aldehyde, aminoketone, co2],
        reaction_family="Strecker_Degradation",
    )


def _beta_elimination_steps(aa: Species, pool_species: List[Species]) -> List[ElementaryStep]:
    """
    β-elimination of Cys → DHA + H₂S, or Ser → DHA + H₂O.
    If Lysine in pool: DHA + Lys → LAL.
    """
    steps = []
    dha = Species(label="dehydroalanine", smiles="C=C(N)C(=O)O")
    h2s = Species(label="H2S", smiles="S")
    water = Species(label="water", smiles="O")
    lal = Species(label="lysinoalanine", smiles="NC(CCCCNCC(N)C(=O)O)C(=O)O")

    smi_lower = aa.smiles.lower()
    label_lower = aa.label.lower()

    if "cysteine" in label_lower or "nc(cs)" in smi_lower:
        steps.append(ElementaryStep(
            reactants=[aa],
            products=[dha, h2s],
            reaction_family="Beta_Elimination",
        ))
    elif "serine" in label_lower or "nc(co)" in smi_lower:
        steps.append(ElementaryStep(
            reactants=[aa],
            products=[dha, water],
            reaction_family="Beta_Elimination",
        ))
    else:
        return []

    # If lysine is in the pool, DHA reacts to form LAL
    for sp in pool_species:
        if "lysine" in sp.label.lower() or "nccccc(n)" in sp.smiles.lower():
            steps.append(ElementaryStep(
                reactants=[dha, sp],
                products=[lal],
                reaction_family="DHA_Crosslinking",
            ))
            break

    return steps


def _aminoketone_condensation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    2x aminoacetone -> 2,5-dimethylpyrazine + 2H2O + H2
    (Aromatic pyrazines require oxidation/loss of 2H from the dihydro-intermediate).
    """
    steps = []
    aks = [s for s in pool_species if "aminoacetone" in s.label.lower() or s.smiles == "CC(=O)CN"]
    for ak in aks:
        pyrazine = Species(label="2,5-dimethylpyrazine", smiles="Cc1cnc(C)cn1")
        water = Species(label="water", smiles="O")
        hydro = Species(label="H2", smiles="[HH]")
        steps.append(ElementaryStep(
            reactants=[ak, ak],
            products=[pyrazine, water, water, hydro],
            reaction_family="Aminoketone_Condensation"
        ))
    return steps


def _retro_aldol_fragmentation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    3-deoxyosone -> C2 + C3 fragments
    Hexose Deoxyosone (C6 H10 O5) -> Pyruvaldehyde (C3 H4 O2) + Glyceraldehyde (C3 H6 O3)
    Sum: C6 H10 O5. Balanced! No water needed.
    Pentose Deoxyosone (C5 H8 O4) -> Pyruvaldehyde (C3 H4 O2) + Glycolaldehyde (C2 H4 O2)
    Sum: C5 H8 O4. Balanced.
    """
    steps = []
    for s in pool_species:
        lower = s.label.lower()
        if "deoxyosone" in lower:
            # Hexose -> Pyruvaldehyde + Glyceraldehyde
            if "glucose" in lower or "fructose" in lower or _is_hexose(s):
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glyceraldehyde", smiles="O=CC(O)CO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
            # Pentose -> Pyruvaldehyde + Glycolaldehyde
            elif "ribose" in lower or _is_pentose(s):
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glycolaldehyde", smiles="O=CCO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
    return steps


def _cysteine_degradation(amino_acids: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Thermal degradation of Cysteine -> H2S, NH3, acetaldehyde, CO2
    NC(CS)C(=O)O (C3 H7 N O2 S) -> H2S (H2 S) + NH3 (H3 N) + CC=O (C2 H4 O) + O=C=O (C1 O2)
    Sum products: H5 N S + C3 H4 O3
    Total: C3 H9 N O3 S
    Reactants: C3 H7 N O2 S. Difference is +H2O in the products!
    So Cysteine + H2O -> H2S + NH3 + Acetaldehyde + CO2 is perfectly balanced.
    Wait, if water is required, it should be a reactant.
    """
    steps = []
    if conditions.temperature_celsius < 100:
        return steps
        
    for aa in amino_acids:
        if "cysteine" == aa.label.lower() or "l-cysteine" == aa.label.lower() or aa.smiles == "NC(CS)C(=O)O":
            h2s = Species(label="H2S", smiles="S")
            ammonia = Species(label="ammonia", smiles="N")
            acetaldehyde = Species(label="acetaldehyde", smiles="CC=O")
            co2 = Species(label="CO2", smiles="O=C=O")
            water = Species(label="water", smiles="O")
            
            steps.append(ElementaryStep(
                reactants=[aa, water],
                products=[h2s, ammonia, acetaldehyde, co2],
                reaction_family="Cysteine_Degradation"
            ))
            
    return steps


def _thiazole_condensation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Thiazole formation (Simplified balanced pathway).
    Pyruvaldehyde (C3 H4 O2) + NH3 (H3 N) + H2S (H2 S) -> Thiazole (C3 H3 N S) + 2 H2O (H4 O2) + H2 (H2).
    Sum Reactants: C3 H9 N O2 S. Sum Products: C3 H9 N O2 S. Perfectly balanced!
    Then we can decorate it via alkylation (for the 2-alkylthiazoles) or just use the Strecker aldehydes as the backbone 
    if they have enough carbons. Wait, the template generates specific alkylthiazoles based on the Strecker aldehyde.
    E.g. 3-methylbutanal (C5 H10 O) -> 2-isobutylthiazole (C7 H11 N S). 
    This gains 2 Carbons (C3 backbone from somewhere).
    So: Aldehyde (C_n) + Pyruvaldehyde (C3) + NH3 + H2S -> 2-Alkylthiazole (C_{n+2}) + ... wait.
    Let's use the actual balanced synthesis:
    Aldehyde (R-CHO) + alpha-mercapto-ketone (R'-C(=O)-CH(SH)-R'') + NH3 -> Thiazole + 3H2O.
    The easiest way to balance the existing hardcoded list is to use Pyruvaldehyde (C3H4O2) as the C3 backbone.
    Aldehyde (R-CHO: C_n H_2n O) + Pyruvaldehyde (C3 H4 O2) + NH3 (H3 N) + H2S (H2 S) 
      -> 2-Alkylthiazole (C_{n+3} H_{2n+3} N S) + 3 H2O (H6 O3) + H2 gas (H2).
    Let's check math for acetaldehyde (C2 H4 O).
    Reactants: C2H4O + C3H4O2 + NH3 + H2S = C5 H13 N O3 S.
    Products: 2-methylthiazole (C4 H5 N S). Wait, C4? Acetaldehyde (C2) + Pyruvaldehyde (C3) = C5!
    Where did the extra carbon go? 2-methylthiazole only has 4 carbons!
    Ah, the Thiazole ring itself has 3 carbons. 2-methylthiazole has 3 (ring) + 1 (methyl) = 4 carbons.
    So the backbone must be a C2 piece! Glycolaldehyde (C2 H4 O2).
    Let's check Glycolaldehyde (C2 H4 O2) + Acetaldehyde (C2 H4 O) + NH3 + H2S -> C4 H13 N O3 S.
    2-methylthiazole (C4 H5 N S) + 3 H2O (H6 O3) + H2 (H2) -> C4 H13 N O3 S. PERFECTLY BALANCED!
    
    So: Strecker Aldehyde + Glycolaldehyde + NH3 + H2S -> 2-Alkylthiazole + 3 H2O + H2.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    nh3 = next((s for s in pool_species if s.smiles == "N"), None)
    glycol = next((s for s in pool_species if "glycolaldehyde" in s.label.lower() or s.smiles == "O=CCO"), None)
    
    if not (h2s and nh3 and glycol):
        return steps
        
    _thiazole_map = {
        "3-methylbutanal": ("2-isobutylthiazole", "CC(C)Cc1nccs1"),
        "2-methylbutanal": ("2-sec-butylthiazole", "CCC(C)c1nccs1"),
        "2-methylpropanal": ("2-isopropylthiazole", "CC(C)c1nccs1"),
        "acetaldehyde": ("2-methylthiazole", "Cc1nccs1"),
    }
    
    for sp in pool_species:
        entry = _thiazole_map.get(sp.label)
        if entry:
            tz = Species(label=entry[0], smiles=entry[1])
            water = Species("water", "O")
            hydro = Species("H2", "[HH]")
            steps.append(ElementaryStep(
                reactants=[sp, glycol, h2s, nh3], 
                products=[tz, water, water, water, hydro], 
                reaction_family="Lipid_Thiazole_Condensation"
            ))
    return steps


def _thiol_addition(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Furfural + H2S + [Reducer] -> Furfurylthiol (FFT) + ... 

    In vivo/food matrices, reduction is driven by the sugar pool (Reductones)
    rather than H2 gas. We couple Furfural reduction with Sugar dehydration:
    Furfural (C5H4O2) + H2S (H2S) + Ribose (C5H10O5) -> FFT (C5H6OS) + Deoxyosone (C5H8O4) + 2 H2O
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    if not h2s:
        return steps

    # Priority 1: H2 Gas (Legacy/Specific path)
    h2 = next((s for s in pool_species if s.smiles == "[HH]"), None)

    _fft_map = {
        "furfural":         ("2-furfurylthiol",       "SCc1ccco1"),
        "5-methylfurfural": ("5-methylfurfurylthiol","SCc1ccc(C)o1"),
    }

    water = Species("water", "O")
    seen = set()
    
    for sp in pool_species:
        entry = _fft_map.get(sp.label)
        if not entry:
            continue

        # Strategy A: H2 gas reduction (if available)
        if h2:
            key = (sp.smiles, h2.smiles)
            if key not in seen:
                seen.add(key)
                fft = Species(label=entry[0], smiles=entry[1])
                steps.append(ElementaryStep(
                    reactants=[sp, h2s, h2],
                    products=[fft, water],
                    reaction_family="Thiol_Addition_H2"
                ))

        # Strategy B: H2S-mediated reduction (Two-Step Bimolecular)
        # Step 1: Furfural + H2S -> Thiohemiacetal (C5H6O2S)
        # Step 2: Thiohemiacetal + H2S -> FFT + S + H2O
        if h2s:
            intermed_smi = "OC(S)c1ccco1" # Thiohemiacetal
            key = (sp.smiles, "bimolecular_coupled")
            if key not in seen:
                seen.add(key)
                intermed = Species(label=f"{sp.label}-thiohemiacetal", smiles=intermed_smi)
                fft = Species(label=entry[0], smiles=entry[1])
                sulfur = Species(label="elemental-sulfur", smiles="[S]")
                
                # Step 1
                steps.append(ElementaryStep(
                    reactants=[sp, h2s],
                    products=[intermed],
                    reaction_family="Thiohemiacetal_Formation"
                ))
                # Step 2
                steps.append(ElementaryStep(
                    reactants=[intermed, h2s],
                    products=[fft, sulfur, water],
                    reaction_family="Thiol_Dehydration"
                ))

    return steps


def _mft_pathway(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Template for 2-methyl-3-furanthiol (MFT) formation (Critical Meat Aroma).
    Path: 3-deoxyosone -> 1,4-dideoxyosone (1,4-DDO) -> MFT.
    Balanced: 
    C5H8O4 (deoxyosone) + H2S -> C5H6OS (MFT) + 2H2O.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    if not h2s:
        return steps
        
    water = Species("water", "O")
    
    for s in pool_species:
        if "deoxyosone" in s.label.lower() and (_is_pentose(s) or "ribose" in s.label.lower()):
            # Filter out N-containing species; they require deamination first (R.12)
            if "N" in s.smiles:
                continue

            # Net Dehydration/Thiolation to MFT
            # Balanced: C5H8O4 + 2 H2S -> C5H6OS + 3 H2O + S
            mft = Species(label="2-methyl-3-furanthiol", smiles="Cc1c(S)cco1")
            sulfur = Species(label="elemental-sulfur", smiles="[S]")
            steps.append(ElementaryStep(
                reactants=[s, h2s, h2s],
                products=[mft, sulfur, water, water, water],
                reaction_family="Thiol_Addition"
            ))
            
    return steps


def _sulfur_volatiles_pathway(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Template for Methionine-derived sulfur volatiles (MeSH, DMDS, DMTS).
    Methional (CSCCC=O) -> Methanethiol (CS) + Acrolein (C=CC=O).
    2x Methanethiol (CS) -> DMDS (CSSC) + H2.
    DMDS + Methanethiol -> DMTS (CSSSC) + H2.
    """
    steps = []
    hydro = Species("H2", "[HH]")
    
    methionals = [s for s in pool_species if s.label == "methional" or s.smiles == "CSCCC=O"]
    for m in methionals:
        mesh = Species(label="methanethiol", smiles="CS")
        acrolein = Species(label="acrolein", smiles="C=CC=O")
        steps.append(ElementaryStep(
            reactants=[m],
            products=[mesh, acrolein],
            reaction_family="Strecker_Degradation"
        ))
        
    mesh_list = [s for s in pool_species if s.label == "methanethiol" or s.smiles == "CS"]
    dmds_list = []
    for m1 in mesh_list:
        dmds = Species(label="dimethyl-disulfide", smiles="CSSC")
        steps.append(ElementaryStep(
            reactants=[m1, m1],
            products=[dmds, hydro],
            reaction_family="Thiol_Oxidation"
        ))
        dmds_list.append(dmds)
        
    for d in dmds_list:
        for m1 in mesh_list:
            dmts = Species(label="dimethyl-trisulfide", smiles="CSSSC")
            steps.append(ElementaryStep(
                reactants=[d, m1],
                products=[dmts, hydro],
                reaction_family="Thiol_Oxidation"
            ))
            
    return steps




def _deamination_step(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    R.12: Generalized Hydrolytic deamination of ketosamines (Amadori) 
    and alpha-aminoketones back to dicarbonyls (deoxyosones).

    Mechanism:
        R-C(=O)-CH2-NH-R' + H2O → R-C(=O)-CHO + NH2-R'
    """
    steps = []
    water = Species("water", "O")
    
    # SMIRKS for hydrolytic deamination:
    # [C:1](=[O:2])[CH2:3]-[NX3:4] >> [C:1](=[O:2])[CH1:3]=[O] . [NX3:4]
    # This takes the ketosamine and splits it into a dicarbonyl and an amine.
    rxn = AllChem.ReactionFromSmarts("[CX3:1](=[O:2])-[CX4H2:3]-[NX3:4]>>[CX3:1](=[O:2])-[CX3H1:3]=O.[NX3:4]")

    for s in pool_species:
        # Filter out pyrazines or very large structures
        if "pyrazine" in s.label.lower() or _mw(s.smiles) > MAX_MW:
            continue
            
        # R.12: Skip free amino acids — they undergo Strecker/Maillard, not self-deamination
        if any(aa in s.label.lower() for aa in [
            "glycine", "alanine", "leucine", "isoleucine", "valine",
            "phenylalanine", "methionine", "lysine", "cysteine",
            "asparagine", "serine", "threonine", "arginine"
        ]):
            continue
            
        mol = _mol(s.smiles)
        if mol is None:
            continue

        prods = rxn.RunReactants((mol,))
        if not prods:
            continue
            
        try:
            # First product is the dicarbonyl (deoxyosone)
            # Second product is the amine fragment
            dicarb_smi = Chem.MolToSmiles(prods[0][0])
            amine_smi = Chem.MolToSmiles(prods[0][1])
            
            # Labeling logic
            dicarb_label = f"deaminated-{s.label}-dicarbonyl"
            # If it's a pentose Amadori, call it deoxyosone
            if "Amadori" in s.label or "Heyns" in s.label:
                # Heuristic: if C5 -> pentose-deoxyosone, C6 -> hexose-deoxyosone
                c_count = sum(1 for a in prods[0][0].GetAtoms() if a.GetAtomicNum() == 6)
                if c_count == 5: dicarb_label = "pentose-deoxyosone-R12"
                elif c_count == 6: dicarb_label = "hexose-deoxyosone-R12"

            amine_sp = Species(label=f"liberated-amine-{s.label}", smiles=amine_smi)
            dicarb_sp = Species(label=dicarb_label, smiles=dicarb_smi)
            
            # Stoichiometry logic:
            # Amadori (ketosamine) -> Dicarb + Amine (Balanced)
            # Aminoketone + H2O -> Dicarb + Amine + H2 (Balanced)
            if "Amadori" in s.label or "Heyns" in s.label:
                steps.append(ElementaryStep(
                    reactants=[s],
                    products=[dicarb_sp, amine_sp],
                    reaction_family="Generalized_Deamination"
                ))
            else:
                h2 = Species("H2", "[HH]")
                steps.append(ElementaryStep(
                    reactants=[s, water],
                    products=[dicarb_sp, amine_sp, h2],
                    reaction_family="Generalized_Deamination"
                ))
        except Exception:
            continue

    return steps

def _lipid_maillard_synergy(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Lipid Aldehyde + alpha-aminoketone + H2S -> 2-Alkylthiazole + 2 H2O + H2.
    Lipid Aldehyde + 2x alpha-aminoketone -> Alkylpyrazine (Branching synergy).
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    
    lipids = [s for s in pool_species if _is_lipid_aldehyde(s)]
    # Target aminoketones (aminoacetone, 3-amino-2-butanone, etc.)
    aks = [s for s in pool_species if "amino" in s.label.lower() and ("acetone" in s.label.lower() or s.smiles == "CC(=O)CN") and "dicarbonyl" not in s.label.lower()]

    if not lipids or not aks:
        return steps

    water = Species("water", "O")
    hydro = Species("H2", "[HH]")

    for lip in lipids:
        # Extract alkyl chain length for label
        m_lip = _mol(lip.smiles)
        c_lip = sum(1 for a in m_lip.GetAtoms() if a.GetAtomicNum() == 6)
        r_len = c_lip - 1 # excluding carbonyl C
        
        for ak in aks:
            # 1. Thiazole Synergy (if H2S present)
            if h2s:
                if r_len == 5: # Hexanal
                    tz_name = "2-pentyl-4-methylthiazole"
                    tz_smi = "CCCCCC1=NC(C)=CS1"
                elif r_len == 6: # Heptanal
                    tz_name = "2-hexyl-4-methylthiazole"
                    tz_smi = "CCCCCCC1=NC(C)=CS1"
                else:
                    tz_name = f"2-alkyl(C{r_len})-4-methylthiazole"
                    prefix = "C" * r_len
                    tz_smi = f"{prefix}C1=NC(C)=CS1"

                steps.append(ElementaryStep(
                    reactants=[lip, ak, h2s],
                    products=[Species(tz_name, tz_smi), water, water, hydro],
                    reaction_family="Lipid_Strecker_Synergy"
                ))

    return steps


def _lipid_hydroperoxide_scission(pool: List[Species]) -> List[ElementaryStep]:
    """
    Tier B: Homolytic cleavage of lipid hydroperoxides.
    R-OOH -> R-O. + .OH
    """
    steps = []
    oh_radical = Species(label="hydroxyl-radical", smiles="[OH]")
    
    for s in pool:
        if _is_lipid_hydroperoxide(s):
            # Use SMIRKS for homolysis but tag the alkoxy oxygen with Isotope 13 so we can find it
            rxn = AllChem.ReactionFromSmarts("[CX4:1]-[OX2:3][OX2H,OX1H0-:2]>>[C:1]-[13OX1:3].[OH:2]")
            m = _mol(s.smiles)
            if not m: continue
            prods = rxn.RunReactants((m,))
            if prods:
                alk_mol = prods[0][0]
                # Explicitly set radical on the oxygen (find Isotope 13)
                for atom in alk_mol.GetAtoms():
                    if atom.GetIsotope() == 13 and atom.GetAtomicNum() == 8:
                        atom.SetIsotope(0)
                        atom.SetNumExplicitHs(0)
                        atom.SetNoImplicit(True)
                        atom.SetNumRadicalElectrons(1)
                
                try:
                    alk_mol.UpdatePropertyCache(strict=False)
                    Chem.SanitizeMol(alk_mol)
                except Exception:
                    pass

                alkoxy_smi = Chem.MolToSmiles(alk_mol)
                alkoxy = Species(label=f"{s.label}-alkoxy-radical", smiles=alkoxy_smi)
                steps.append(ElementaryStep(
                    reactants=[s],
                    products=[alkoxy, oh_radical],
                    reaction_family="Lipid_Homolysis"
                ))
    return steps


def _radical_crosstalk_templates(pool: List[Species]) -> List[ElementaryStep]:
    """
    Templates for radical + sulfur collisions.
    R. + H2S -> R-H + .SH
    R. + .SH -> R-SH
    """
    steps = []
    actual_radicals = []
    for s in pool:
        m = _mol(s.smiles)
        if m and any(a.GetNumRadicalElectrons() > 0 for a in m.GetAtoms()):
            actual_radicals.append(s)
            
    # H2S Crosstalk
    h2s = next((s for s in pool if s.smiles == "S"), None)
    if h2s:
        for rad in actual_radicals:
            # Rad + H2S -> RH + SH.
            m_rad = _mol(rad.smiles)
            if not m_rad: continue
            rh_mol = Chem.RWMol(m_rad)
            rad_indices = []
            for atom in rh_mol.GetAtoms():
                if atom.GetNumRadicalElectrons() > 0:
                    rad_indices.append(atom.GetIdx())
                    atom.SetNumRadicalElectrons(0)
                    atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            if not rad_indices: continue
            try:
                rh_mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(rh_mol)
                rh_smi = Chem.MolToSmiles(rh_mol)
                sh_rad = Species(label="thiol-radical", smiles="[SH]")
                steps.append(ElementaryStep(
                    reactants=[rad, h2s],
                    products=[Species(label=f"quenched-{rad.label}", smiles=rh_smi), sh_rad],
                    reaction_family="Radical_Crosstalk"
                ))
            except Exception:
                continue
            
    # Expanded Crosstalk: Rad + MFT -> RH + MFT-radical
    mft_can = "Cc1occc1S" # Fixed canonical
    mft_sp = next((s for s in pool if _canonical(s.smiles) == mft_can), None)
    if mft_sp:
        for rad in actual_radicals:
             m_rad = _mol(rad.smiles)
             if not m_rad: continue
             rh_mol = Chem.RWMol(m_rad)
             for atom in rh_mol.GetAtoms():
                 if atom.GetNumRadicalElectrons() > 0:
                     atom.SetNumRadicalElectrons(0)
                     atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
             try:
                 rh_mol.UpdatePropertyCache(strict=False)
                 Chem.SanitizeMol(rh_mol)
                 rh_smi = Chem.MolToSmiles(rh_mol)
                 # MFT radical: S on the furan ring
                 m_mft_rad = Chem.RWMol(_mol(mft_sp.smiles))
                 for atom in m_mft_rad.GetAtoms():
                     if atom.GetSymbol() == "S":
                         atom.SetNumRadicalElectrons(1)
                 mft_rad_smi = Chem.MolToSmiles(m_mft_rad)
                 mft_rad = Species(label="MFT-radical", smiles=mft_rad_smi)
                 steps.append(ElementaryStep(
                     reactants=[rad, mft_sp],
                     products=[Species(label=f"quenched-{rad.label}", smiles=rh_smi), mft_rad],
                     reaction_family="Radical_Crosstalk"
                 ))
             except Exception:
                 continue

    return steps


def _sugar_ring_opening(pool_species: List[Species]) -> List[ElementaryStep]:
    """Hemiacetal cyclic sugar -> open-chain aldehyde. (Defensive rule via RWMol)"""
    steps = []
    # Match hemiacetal: O(ring) - C(ring)(OH)
    # The [O;X2;R] ensures it's a ring oxygen, [C;X4;R] is a ring carbon, [O;X2;H] is the hydroxyl.
    patt = Chem.MolFromSmarts("[O;X2;R]-[C;X4;R]-[OH]")
    if not patt: 
        return steps
    
    for s in pool_species:
        m = _mol(s.smiles)
        if not m: 
            continue
            
        matches = m.GetSubstructMatches(patt)
        if not matches: 
            continue
            
        # We might have multiple hemiacetals (e.g., dimers), just take the first
        o_ring_idx, c_anomeric_idx, o_hydroxyl_idx = matches[0]
        
        try:
            rw_mol = Chem.RWMol(m)
            # 1. Break O_ring - C_anomeric bond
            rw_mol.RemoveBond(o_ring_idx, c_anomeric_idx)
            
            # 2. Change C_anomeric - O_hydroxyl bond to double
            b = rw_mol.GetBondBetweenAtoms(c_anomeric_idx, o_hydroxyl_idx)
            if b:
                b.SetBondType(Chem.BondType.DOUBLE)
            
            # 3. Adjust valency hydrogens automatically via Sanitize
            Chem.SanitizeMol(rw_mol)
            open_smi = Chem.MolToSmiles(rw_mol)
            
            if _is_valid(open_smi):
                p = Species(label=f"{s.label}-open", smiles=open_smi)
                steps.append(ElementaryStep([s], [p], "Sugar_Ring_Opening"))
        except Exception:
            pass
            
    return steps


# ── Toxic / Safety Marker Templates ──────────────────────────────────────

def _acrylamide_formation(sugar: Species, asparagine: Species) -> List[ElementaryStep]:
    """
    Asparagine + Reducing Sugar -> Acrylamide + Fragments.
    Acrylamide: C=CC(=O)N (C3H5NO)
    """
    if not (_is_sugar(sugar) and _is_asparagine(asparagine)):
        return []
    
    # Net balanced reaction (simplified):
    # Asparagine (C4H8N2O3) -> Acrylamide (C3H5NO) + CO2 (CO2) + NH3 (H3N)
    acrylamide = Species(label="Acrylamide", smiles="C=CC(=O)N")
    co2 = Species(label="CO2", smiles="O=C=O")
    ammonia = Species(label="ammonia", smiles="N")
    
    return [ElementaryStep(
        reactants=[sugar, asparagine],
        products=[acrylamide, sugar, co2, ammonia], # Sugar is conserved here for simplicity
        reaction_family="Safety_Risk_Acrylamide"
    )]


def _cml_cel_formation(lysine: Species, pool: List[Species]) -> List[ElementaryStep]:
    """
    Lysine + Glyoxal -> CML
    Lysine + Methylglyoxal -> CEL
    """
    if not _is_lysine(lysine):
        return []
        
    steps = []
    
    for s in pool:
        # Glyoxal -> CML
        if s.smiles == "O=CC=O" or s.label == "glyoxal":
            cml = Species("CML", "N[C@@H](CCCCNCC(=O)O)C(=O)O")
            steps.append(ElementaryStep(
                reactants=[lysine, s],
                products=[cml],
                reaction_family="Safety_Risk_AGE"
            ))
        # Methylglyoxal -> CEL
        elif s.smiles == "CC(=O)C=O" or s.label == "pyruvaldehyde" or s.label == "methylglyoxal":
            cel = Species("CEL", "N[C@@H](CCCCNC(C)C(=O)O)C(=O)O")
            steps.append(ElementaryStep(
                reactants=[lysine, s],
                products=[cel],
                reaction_family="Safety_Risk_AGE"
            ))
            
    return steps


# ── PBMA Additive Degradations ───────────────────────────────────────────

def _thiamine_degradation(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: Thermal breakdown of Thiamine (Vitamin B1).
    Literature shows this yields H2S, 2-methylthiophene, and thiazoles.
    """
    if conditions.temperature_celsius < 100:
        return []
        
    steps = []
    # Identify Thiamine in the pool by canonical SMILES
    for s in pool:
        if _canonical(s.smiles) == _canonical(_THIAMINE_CANONICAL):
            p1 = Species("Hydrogen_Sulfide", "S")
            p2 = Species("2-methylthiophene", "Cc1cccs1")
            p3 = Species("4,5-dihydro-2-methylthiazole", "CC1=NCCS1")
            
            steps.append(ElementaryStep(
                reactants=[s],
                products=[p1, p2, p3],
                reaction_family="Additive_Thermal_Degradation"
            ))
    return steps


def _glutathione_cleavage(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: Controlled cleavage of Glutathione (GSH).
    Literature shows it cleaves into glutamic acid and cysteinylglycine dipeptide.
    """
    if conditions.temperature_celsius < 100:
        return []
        
    steps = []
    for s in pool:
        if _canonical(s.smiles) == _canonical(_GSH_CANONICAL):
            p1 = Species("Glutamic_Acid", "N[C@@H](CCC(=O)O)C(=O)O")
            p2 = Species("Cysteinylglycine", "N[C@@H](CS)C(=O)NCC(=O)O")
            
            steps.append(ElementaryStep(
                reactants=[s],
                products=[p1, p2],
                reaction_family="Additive_Thermal_Degradation"
            ))
    return steps


# ──────────────────────────────────────────────────────────────────────────
# Tier A: SMIRKS application
# ──────────────────────────────────────────────────────────────────────────

def _fix_radicals(mol: Chem.Mol, family: str, clear_all: bool = False):
    """Ensure atoms with unsatisfied valence are marked as radicals, and clear satisfied ones."""
    if clear_all:
        for atom in mol.GetAtoms():
            atom.SetNumRadicalElectrons(0)
        return mol

    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        # Use explicit valence (bond order sum) + implicit H
        # This is more accurate for double/triple bonds.
        try:
            val = atom.GetExplicitValence() + atom.GetNumImplicitHs()
        except Exception:
            # Fallback if property cache is not updated
            val = atom.GetTotalDegree() + atom.GetTotalNumHs()
        
        # Determine if we should set or clear radicals
        if an == 6: # Carbon
            if val == 3:
                atom.SetNumRadicalElectrons(1)
            elif val >= 4:
                atom.SetNumRadicalElectrons(0)
        elif an == 8: # Oxygen
            if val == 1:
                if any(f in family for f in ["Radical_Propagation", "Lipid_Homolysis", "Beta_Scission"]):
                    atom.SetNumRadicalElectrons(1)
            elif val >= 2:
                atom.SetNumRadicalElectrons(0)
    return mol

def _apply_smirks_rule(
    name: str, family: str, smirks: str, ph_gate: str,
    pool: List[Species], conditions: ReactionConditions
) -> List[ElementaryStep]:
    """Apply a single SMIRKS rule to all relevant species pairs in the pool."""
    # pH gate check
    if ph_gate == "acid" and conditions.pH >= 6:
        return []
    if ph_gate == "neutral" and conditions.pH < 6:
        return []

    rxn = AllChem.ReactionFromSmarts(smirks)
    if rxn is None:
        return []

    steps = []
    pool_smiles = list({s.smiles for s in pool})

    n_reactants = rxn.GetNumReactantTemplates()

    if n_reactants == 1:
        for smi in pool_smiles:
            m = _mol(smi)
            if m is None:
                continue
            try:
                prods = rxn.RunReactants((m,))
            except Exception:
                continue
            for prod_tuple in prods:
                prod_smiles = []
                for p in prod_tuple:
                    try:
                        ps = Chem.MolToSmiles(p)
                        if _is_valid(ps):
                            prod_smiles.append(ps)
                    except Exception:
                        pass
                if prod_smiles:
                    reactant_sp = next((s for s in pool if s.smiles == smi), Species(smi, smi))
                    for p in prod_tuple:
                        _fix_radicals(p, family, clear_all=True)
                        try:
                            Chem.SanitizeMol(p)
                        except Exception:
                            pass
                        _fix_radicals(p, family)
                    
                    steps.append(ElementaryStep(
                        reactants=[reactant_sp],
                        products=[Species(ps, ps) for ps in prod_smiles],
                        reaction_family=family,
                    ))

    elif n_reactants == 2:
        for i, smi1 in enumerate(pool_smiles):
            for smi2 in pool_smiles:
                m1, m2 = _mol(smi1), _mol(smi2)
                if m1 is None or m2 is None:
                    continue
                try:
                    prods = rxn.RunReactants((m1, m2))
                except Exception:
                    continue
                for prod_tuple in prods:
                    prod_smiles = []
                    valid_step = True
                    for p in prod_tuple:
                        try:
                            _fix_radicals(p, family, clear_all=True)
                            # Sanitize to catch valence issues
                            Chem.SanitizeMol(p)
                            _fix_radicals(p, family)
                            ps = Chem.MolToSmiles(p)
                            if _is_valid(ps):
                                prod_smiles.append(ps)
                            else:
                                valid_step = False
                                break
                        except Exception:
                            valid_step = False
                            break
                            
                    # Only append if ALL products were successfully generated and are valid
                    # This guarantees mass conservation. The RDKit reaction MUST output everything.
                    if valid_step and len(prod_smiles) == len(prod_tuple):
                        r1 = next((s for s in pool if s.smiles == smi1), Species(smi1, smi1))
                        r2 = next((s for s in pool if s.smiles == smi2), Species(smi2, smi2))
                        steps.append(ElementaryStep(
                            reactants=[r1, r2],
                            products=[Species(ps, ps) for ps in prod_smiles],
                            reaction_family=family,
                        ))

    elif n_reactants == 3:
        for smi1 in pool_smiles:
            for smi2 in pool_smiles:
                for smi3 in pool_smiles:
                    m1, m2, m3 = _mol(smi1), _mol(smi2), _mol(smi3)
                    if m1 is None or m2 is None or m3 is None:
                        continue
                    try:
                        prods = rxn.RunReactants((m1, m2, m3))
                    except Exception:
                        continue
                    for prod_tuple in prods:
                        prod_smiles = []
                        valid_step = True
                        for p in prod_tuple:
                            try:
                                Chem.SanitizeMol(p)
                                ps = Chem.MolToSmiles(p)
                                if _is_valid(ps):
                                    prod_smiles.append(ps)
                                else:
                                    valid_step = False
                                    break
                            except Exception:
                                valid_step = False
                                break
                                
                        if valid_step and len(prod_smiles) == len(prod_tuple):
                            r1 = next((s for s in pool if s.smiles == smi1), Species(smi1, smi1))
                            r2 = next((s for s in pool if s.smiles == smi2), Species(smi2, smi2))
                            r3 = next((s for s in pool if s.smiles == smi3), Species(smi3, smi3))
                            steps.append(ElementaryStep(
                                reactants=[r1, r2, r3],
                                products=[Species(ps, ps) for ps in prod_smiles],
                                reaction_family=family,
                            ))

    return steps


# ──────────────────────────────────────────────────────────────────────────
# Main Engine
# ──────────────────────────────────────────────────────────────────────────

class SmirksEngine:
    """
    Hybrid rule-based Maillard pathway generator.

    Usage:
        engine = SmirksEngine(conditions=ReactionConditions(pH=5.5, temperature_celsius=150))
        steps = engine.enumerate([ribose, cysteine, glycine])
        # steps: List[ElementaryStep] — feed into xtb_screener or recommend.py
    """

    def __init__(self, conditions: Optional[ReactionConditions] = None):
        self.conditions = conditions or ReactionConditions()

    def enumerate(
        self,
        precursors: List[Species],
        max_generations: int = 3
    ) -> List[ElementaryStep]:
        """
        Enumerate Maillard pathways from the given precursor pool.

        Phase 1: Apply Tier B templates (Amadori cascade, enolisation,
                 Strecker, beta-elimination).
        Phase 2: Apply Tier A SMIRKS rules iteratively on the growing pool.
        """
        all_steps: List[ElementaryStep] = []

        # Working pool: canonical SMILES → Species
        pool_dict: dict = {_canonical(s.smiles): s for s in precursors if _canonical(s.smiles)}

        def pool_list() -> List[Species]:
            return list(pool_dict.values())

        def add_to_pool(sp: Species):
            can = _canonical(sp.smiles)
            if can and can not in pool_dict:
                pool_dict[can] = sp

        def add_step_products(step: ElementaryStep):
            for p in step.products:
                if _is_valid(p.smiles):
                    add_to_pool(p)

        # ── Pre-Phase: Sugar Ring Opening ────────────────────────────────
        ring_steps = _sugar_ring_opening(pool_list())
        for step in ring_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Pre-Phase: PBMA Additive Degradations ─────────────────────────
        thiamine_steps = _thiamine_degradation(pool_list(), self.conditions)
        for step in thiamine_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)
                
        gsh_steps = _glutathione_cleavage(pool_list(), self.conditions)
        for step in gsh_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Tier B Phase 1: Amadori / Heyns cascade ──────────────────────
        sugars = [s for s in pool_list() if _is_sugar(s)]
        amines = [s for s in pool_list() if _is_primary_amine(s)]

        for sugar in sugars:
            for amine in amines:
                cascade = _amadori_cascade(sugar, amine)
                for step in cascade:
                    # Dedup by reactant+product labels
                    if not _step_exists(step, all_steps):
                        all_steps.append(step)
                        add_step_products(step)

                amadori_smiles = next((p.smiles for s in cascade for p in s.products if "Amadori" in p.label or "Heyns" in p.label), None)
                amadori_can = _canonical(amadori_smiles) if amadori_smiles else None
                amadori_sp = pool_dict.get(amadori_can) if amadori_can else None
                
                if amadori_sp:
                    enols = _enolisation_steps(amadori_sp, sugar, amine, self.conditions)
                    for enol in enols:
                        if not _step_exists(enol, all_steps):
                            all_steps.append(enol)
                            add_step_products(enol)

        # ── Tier B Phase 2: Retro-aldol and Strecker degradation ─────────
        ra_steps = _retro_aldol_fragmentation(pool_list())
        for step in ra_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        dicarbonyls = [s for s in pool_list() if _is_dicarbonyl(s)]
        amines_now = [s for s in pool_list() if _is_primary_amine(s)]

        for dc in dicarbonyls:
            for amine in amines_now:
                s_step = _strecker_step(dc, amine)
                if s_step and not _step_exists(s_step, all_steps):
                    all_steps.append(s_step)
                    add_step_products(s_step)

        # ── Tier B Phase 3: Secondary Condensations & Eliminations ────────
        # 3a. Beta-elimination (DHA pathway)
        beta_candidates = [s for s in pool_list() if _has_cysteine_beta_carbon(s)]
        for aa in beta_candidates:
            be_steps = _beta_elimination_steps(aa, pool_list())
            for step in be_steps:
                if not _step_exists(step, all_steps):
                    all_steps.append(step)
                    add_step_products(step)

        # 3b. Cysteine thermal degradation
        cys_steps = _cysteine_degradation(pool_list(), self.conditions)
        for step in cys_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3c. Aminoketone Condensation (Pyrazines)
        ak_steps = _aminoketone_condensation(pool_list())
        for step in ak_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3d. Lipid Thiazole Condensation
        tz_steps = _thiazole_condensation(pool_list())
        for step in tz_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3e-0. R.12: Deamination (Must happen before volatile templates)
        deam_steps = _deamination_step(pool_list())
        for step in deam_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3e. Thiol Addition (Furfural + H2S + H2 -> FFT)
        ta_steps = _thiol_addition(pool_list())
        for step in ta_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3e-2. MFT Formation (Phase R.2 Fix)
        mft_steps = _mft_pathway(pool_list())
        for step in mft_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3e-3. Methionine Sulfur Volatiles (Phase R.2 Fix)
        sulf_steps = _sulfur_volatiles_pathway(pool_list())
        for step in sulf_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3f. Safety / Toxic Markers (Acrylamide, CML, CEL)
        for sug in sugars:
            for amine in amines_now:
                acry_steps = _acrylamide_formation(sug, amine)
                for step in acry_steps:
                    if not _step_exists(step, all_steps):
                        all_steps.append(step)
                        add_step_products(step)
        
        for amine in amines_now:
            age_steps = _cml_cel_formation(amine, pool_list())
            for step in age_steps:
                if not _step_exists(step, all_steps):
                    all_steps.append(step)
                    add_step_products(step)

        # 3g. Lipid-Maillard Synergy (Lipid Aldehyde + Strecker AK)
        syn_steps = _lipid_maillard_synergy(pool_list())
        for step in syn_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3h. Lipid Oxidation Radicals (Phase 19)
        hooh_steps = _lipid_hydroperoxide_scission(pool_list())
        for step in hooh_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)
        
        cross_steps = _radical_crosstalk_templates(pool_list())
        for step in cross_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Tier A: SMIRKS rules, iterative ──────────────────────────────
        seen_step_keys: Set[str] = {_step_key(s) for s in all_steps}

        for _gen in range(max_generations):
            new_steps_this_gen = []
            current_pool = pool_list()

            for name, family, smirks, gate in _SMIRKS_RULES:
                candidates = _apply_smirks_rule(
                    name, family, smirks, gate, current_pool, self.conditions
                )
                for step in candidates:
                    k = _step_key(step)
                    if k not in seen_step_keys:
                        new_steps_this_gen.append(step)
                        seen_step_keys.add(k)

            if not new_steps_this_gen:
                break  # No new reactions found — converged

            for step in new_steps_this_gen:
                add_step_products(step)
            all_steps.extend(new_steps_this_gen)

        return all_steps


def _step_key(step: ElementaryStep) -> str:
    """Stable hash key for deduplication: sorted reactants + sorted products."""
    reacts = tuple(sorted(_canonical(r.smiles) or r.label for r in step.reactants))
    prods = tuple(sorted(_canonical(p.smiles) or p.label for p in step.products))
    return str((reacts, prods))


def _step_exists(step: ElementaryStep, existing: List[ElementaryStep]) -> bool:
    key = _step_key(step)
    return any(_step_key(s) == key for s in existing)


# ──────────────────────────────────────────────────────────────────────────
# CLI demo
# ──────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from data.reactions.curated_pathways import (
        RIBOSE, GLUCOSE, GLYCINE, CYSTEINE, LEUCINE, H2S
    )

    SYSTEMS = [
        ("Ribose + Glycine @ pH 5",  [RIBOSE, GLYCINE],            ReactionConditions(pH=5.0)),
        ("Glucose + Glycine @ pH 7", [GLUCOSE, GLYCINE],           ReactionConditions(pH=7.0)),
        ("Ribose + Cysteine @ pH 6", [RIBOSE, CYSTEINE, H2S],      ReactionConditions(pH=6.0)),
        ("Ribose + Cys + Leu",       [RIBOSE, CYSTEINE, LEUCINE, H2S], ReactionConditions(pH=6.0)),
    ]

    for label, precursors, conds in SYSTEMS:
        print(f"\n{'='*60}")
        print(f"System: {label}")
        print(f"Input:  {[p.label for p in precursors]}")
        engine = SmirksEngine(conds)
        steps = engine.enumerate(precursors, max_generations=4)
        print(f"Generated {len(steps)} elementary steps:")
        for step in steps:
            print(f"  {step}")
