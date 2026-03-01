"""
data/reactions/curated_pathways.py — Hand-curated Maillard reaction pathways.

Defines the 5 core Maillard cascades (A–E) as explicit ElementaryStep objects
with validated SMILES, ready to feed into the xTB screening pipeline.

Sources: pathways.md, Maillard_meat.md, Maillard_Plant_based.md
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.pathway_extractor import Species, ElementaryStep


def _s(label: str, smiles: str) -> Species:
    """Shorthand to create a Species."""
    return Species(label=label, smiles=smiles)


# ──────────────────────────────────────────────────────────────────────
# Species Library (validated SMILES)
# ──────────────────────────────────────────────────────────────────────

# Sugars (open-chain aldehyde form)
RIBOSE       = _s("D-ribose",        "O=CC(O)C(O)C(O)CO")
GLUCOSE      = _s("D-glucose",       "O=CC(O)C(O)C(O)C(O)CO")

# Amino acids
GLYCINE      = _s("glycine",         "NCC(=O)O")
CYSTEINE     = _s("L-cysteine",      "NC(CS)C(=O)O")
LEUCINE      = _s("L-leucine",       "CC(C)CC(N)C(=O)O")
LYSINE       = _s("L-lysine",        "NCCCCC(N)C(=O)O")
SERINE       = _s("L-serine",        "NC(CO)C(=O)O")

# Small molecules
WATER        = _s("water",           "O")
CO2          = _s("CO2",             "O=C=O")
H2S          = _s("H2S",            "S")
NH3          = _s("ammonia",         "N")

# Key intermediates
PYRUVALDEHYDE = _s("pyruvaldehyde",   "CC(=O)C=O")
FURFURAL      = _s("furfural",       "O=Cc1ccco1")
FFT           = _s("2-furfurylthiol", "SCc1ccco1")
MFT           = _s("2-methyl-3-furanthiol", "Cc1occc1S")

# Strecker aldehydes
METHYLBUTANAL_3 = _s("3-methylbutanal", "CC(C)CC=O")

# DHA pathway
DHA           = _s("dehydroalanine",   "C=C(N)C(=O)O")
LAL           = _s("lysinoalanine",    "NC(CCCCNCC(N)C(=O)O)C(=O)O")

# Schiff bases and Amadori products
RIBOSE_GLY_SCHIFF = _s("ribose-glycine-Schiff-base",
                        "OCC(O)C(O)C(O)/C=N/CC(=O)O")
RIBOSE_GLY_AMADORI = _s("ribose-glycine-Amadori",
                         "OCC(O)C(O)C(=O)CNCC(=O)O")
DEOXYOSONE_3  = _s("3-deoxyosone",    "O=CC(=O)CC(O)CO")

# Lipid off-flavour
HEXANAL       = _s("hexanal",         "CCCCCC=O")
HEXANAL_SCHIFF = _s("hexanal-glycine-Schiff-base",
                     "CCCCC/C=N/CC(=O)O")

# Alpha-aminoketone (Strecker co-product)
AMINOACETONE  = _s("aminoacetone",    "CC(=O)CN")


# ──────────────────────────────────────────────────────────────────────
# Pathway Definitions
# ──────────────────────────────────────────────────────────────────────

PATHWAYS = {

    # ── Pathway A: Core Maillard Cascade (Ribose + Glycine) ──────────
    "A_Core_Maillard_Ribose_Gly": [
        # Step 1: Schiff base formation
        ElementaryStep(
            reactants=[RIBOSE, GLYCINE],
            products=[RIBOSE_GLY_SCHIFF, WATER],
            reaction_family="Schiff_Base_Formation",
        ),
        # Step 2: Amadori rearrangement (1,2-proton shift)
        ElementaryStep(
            reactants=[RIBOSE_GLY_SCHIFF],
            products=[RIBOSE_GLY_AMADORI],
            reaction_family="Amadori_Rearrangement",
        ),
        # Step 3: 1,2-Enolisation → 3-deoxyosone
        ElementaryStep(
            reactants=[RIBOSE_GLY_AMADORI],
            products=[DEOXYOSONE_3, GLYCINE],
            reaction_family="Enolisation",
        ),
        # Step 4: Cyclisation/dehydration → furfural
        ElementaryStep(
            reactants=[DEOXYOSONE_3],
            products=[FURFURAL, WATER, WATER],
            reaction_family="Sugar_Dehydration",
        ),
    ],

    # ── Pathway B: Strecker Degradation (Leucine → 3-methylbutanal) ──
    "B_Strecker_Leu": [
        # Step 1: α-dicarbonyl + amino acid → Strecker aldehyde + aminoketone + CO₂
        ElementaryStep(
            reactants=[PYRUVALDEHYDE, LEUCINE],
            products=[METHYLBUTANAL_3, AMINOACETONE, CO2],
            reaction_family="Strecker_Degradation",
        ),
    ],

    # ── Pathway C: S-Maillard (Ribose + Cysteine → FFT) ─────────────
    "C_S_Maillard_FFT": [
        # Step 1: Cysteine thermal degradation
        ElementaryStep(
            reactants=[CYSTEINE],
            products=[PYRUVALDEHYDE, H2S, NH3],
            reaction_family="Cysteine_Degradation",
        ),
        # Step 2: Ribose dehydration → furfural
        ElementaryStep(
            reactants=[RIBOSE],
            products=[FURFURAL, WATER, WATER, WATER],
            reaction_family="Sugar_Dehydration",
        ),
        # Step 3: H₂S + furfural → FFT + H₂O
        ElementaryStep(
            reactants=[FURFURAL, H2S],
            products=[FFT, WATER],
            reaction_family="Thiol_Addition",
        ),
    ],

    # ── Pathway D: Off-flavour trapping (Hexanal + Glycine) ──────────
    "D_Offflavour_Trapping": [
        # Step 1: Hexanal + amino acid → non-volatile Schiff base
        ElementaryStep(
            reactants=[HEXANAL, GLYCINE],
            products=[HEXANAL_SCHIFF, WATER],
            reaction_family="Lipid_Schiff_Base",
        ),
    ],

    # ── Pathway E: DHA Competition (Cysteine → DHA → LAL) ────────────
    "E_DHA_Competition": [
        # Step 1: β-elimination → DHA + H₂S
        ElementaryStep(
            reactants=[CYSTEINE],
            products=[DHA, H2S],
            reaction_family="Beta_Elimination",
        ),
        # Step 2: DHA + Lysine → Lysinoalanine (LAL)
        ElementaryStep(
            reactants=[DHA, LYSINE],
            products=[LAL],
            reaction_family="DHA_Crosslinking",
        ),
    ],
}


if __name__ == "__main__":
    print("Curated Maillard Pathways")
    print("=" * 60)
    for name, steps in PATHWAYS.items():
        print(f"\n{name} ({len(steps)} steps):")
        for i, step in enumerate(steps, 1):
            print(f"  {i}. {step}")
