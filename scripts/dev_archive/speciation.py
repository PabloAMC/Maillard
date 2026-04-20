import math
from typing import Dict, Any

# Standard pKa values for Maillard-relevant functional groups
# [alpha-amino, side-chain]
PKA_REGISTRY = {
    "lysine": {"alpha_amino": 8.95, "side_chain": 10.53},
    "arginine": {"alpha_amino": 9.04, "side_chain": 12.48},
    "histidine": {"alpha_amino": 9.17, "side_chain": 6.00},
    "cysteine": {"alpha_amino": 10.78, "side_chain": 8.33, "thiol": 8.33},
    "glycine": {"alpha_amino": 9.60, "side_chain": None},
    "alanine": {"alpha_amino": 9.69, "side_chain": None},
    "leucine": {"alpha_amino": 9.60, "side_chain": None},
    "isoleucine": {"alpha_amino": 9.68, "side_chain": None},
    "valine": {"alpha_amino": 9.62, "side_chain": None},
    "methionine": {"alpha_amino": 9.21, "side_chain": None},
    "proline": {"alpha_amino": 10.60, "side_chain": None},
    "phenylalanine": {"alpha_amino": 9.13, "side_chain": None},
    "tyrosine": {"alpha_amino": 9.11, "side_chain": 10.07}, # phenolic OH
    "tryptophan": {"alpha_amino": 9.39, "side_chain": None},
    "serine": {"alpha_amino": 9.15, "side_chain": None},
    "threonine": {"alpha_amino": 9.10, "side_chain": None},
    "glutamine": {"alpha_amino": 9.13, "side_chain": None},
    "asparagine": {"alpha_amino": 8.80, "side_chain": None},
    "glutamic_acid": {"alpha_amino": 9.47, "side_chain": 4.07}, # carboxyl
    "aspartic_acid": {"alpha_amino": 9.60, "side_chain": 3.65}, # carboxyl
    "generic_amino_acid": {"alpha_amino": 9.50, "side_chain": None},
    
    # Simple amines
    "ammonia": {"alpha_amino": 9.25, "side_chain": None},
}

def henderson_hasselbalch_unprotonated_fraction(pH: float, pKa: float) -> float:
    """
    Returns the fraction of the unprotonated (basic) form [A-] or [B].
    For amines/thiols acting as nucleophiles in Maillard, the unprotonated
    free base (R-NH2, R-S-) is the reactive species.
    
    Fraction = 10^(pH - pKa) / (1 + 10^(pH - pKa))
    """
    try:
        ratio = 10.0 ** (pH - pKa)
        return ratio / (1.0 + ratio)
    except OverflowError:
        return 1.0 if pH > pKa else 0.0

def calculate_active_fraction(
    molecule_name: str, 
    pH: float, 
    is_peptide_bound: bool = False, 
    denaturation_state: float = 1.0
) -> float:
    """
    Calculates the effective reactive fraction of an amino acid or amine
    based on Henderson-Hasselbalch speciation and structural accessibility.
    
    Args:
        molecule_name: Name of the precursor (e.g. 'lysine', 'cysteine', 'glycine').
        pH: Reaction pH.
        is_peptide_bound: If True, indicates the amino acid is intrinsic to the
                          protein matrix (not a free monomer spike).
        denaturation_state: Structural accessibility (0.0=folded/native, 1.0=fully denatured).
                            Only impacts side-chains of peptide-bound residues.
    
    Returns:
        Effective reactive fraction (0.0 to 1.0).
    """
    name_lower = molecule_name.lower().replace("-", "_").replace(" ", "_")
    
    # Strip common prefixes like 'l-' or 'd-'
    if name_lower.startswith("l_") or name_lower.startswith("d_"):
        name_lower = name_lower[2:]
        
    # If the molecule is not in the PKA_REGISTRY, it's likely a sugar, lipid, or additive
    # lacking an active amine/thiol. We do not speciate these (fraction = 1.0).
    if name_lower not in PKA_REGISTRY:
        if not any(x in name_lower for x in ["amine", "acid", "glycine", "alanine"]):
            return 1.0
        pka_data = PKA_REGISTRY["generic_amino_acid"]
    else:
        pka_data = PKA_REGISTRY[name_lower]
        
    alpha_pka = pka_data.get("alpha_amino")
    side_pka = pka_data.get("side_chain")
    
    # 1. Calculate chemical active fractions (unprotonated free base)
    alpha_active = henderson_hasselbalch_unprotonated_fraction(pH, alpha_pka) if alpha_pka else 0.0
    side_active = henderson_hasselbalch_unprotonated_fraction(pH, side_pka) if side_pka else 0.0
    
    # 2. Apply structural Priors
    if is_peptide_bound:
        # If bound in a protein chain, the alpha-amino group is tied up in a peptide bond.
        # Only the N-terminus is free. We assume N-termini represent ~1% of total residues for a large intact protein.
        n_terminus_abundance = 0.01 
        effective_alpha = alpha_active * n_terminus_abundance
        
        if side_pka:
            # Side chains (Lys, Arg, Cys, His) are present, but buried vs exposed.
            # Denaturation state scales their accessibility.
            # Base accessibility = 0.2 (even perfectly folded proteins have some surface residues)
            base_accessibility = 0.2
            accessibility = base_accessibility + denaturation_state * (1.0 - base_accessibility)
            effective_side = side_active * accessibility
            
            # The total effective fraction is dominated by the side chain for reactive residues like Lysine
            return max(effective_alpha, effective_side)
        else:
            # For aliphatic AAs lacking reactive side chains, their only contribution when bound is at the N-terminus
            return effective_alpha
    else:
        # Free amino acids are fully accessible.
        # Both alpha and side-chain (if present) can react.
        # For Lysine, both amines can form Schiff bases, but typically one reacts faster.
        # We take the maximum active fraction as the governing rate-limiting pool.
        if side_pka:
            return max(alpha_active, side_active)
        else:
            return alpha_active
