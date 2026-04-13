"""Hand-curated Maillard reaction pathways used by the screening pipeline."""

from src.pathway_extractor import ElementaryStep, Species


def _species(label: str, smiles: str) -> Species:
    return Species(label=label, smiles=smiles)


RIBOSE = _species("D-ribose", "O=CC(O)C(O)C(O)CO")
GLUCOSE = _species("D-glucose", "O=CC(O)C(O)C(O)C(O)CO")
GLYCINE = _species("glycine", "NCC(=O)O")
CYSTEINE = _species("L-cysteine", "NC(CS)C(=O)O")
LEUCINE = _species("L-leucine", "CC(C)CC(N)C(=O)O")
LYSINE = _species("L-lysine", "NCCCCC(N)C(=O)O")
WATER = _species("water", "O")
CO2 = _species("CO2", "O=C=O")
H2S = _species("H2S", "S")
NH3 = _species("ammonia", "N")
PYRUVALDEHYDE = _species("pyruvaldehyde", "CC(=O)C=O")
FURFURAL = _species("furfural", "O=Cc1ccco1")
FFT = _species("2-furfurylthiol", "SCc1ccco1")
METHYLBUTANAL_3 = _species("3-methylbutanal", "CC(C)CC=O")
ACETALDEHYDE = _species("acetaldehyde", "CC=O")
DHA = _species("dehydroalanine", "C=C(N)C(=O)O")
LAL = _species("lysinoalanine", "NC(CCCCNCC(N)C(=O)O)C(=O)O")
RIBOSE_GLY_SCHIFF = _species("ribose-glycine-Schiff-base", "OCC(O)C(O)C(O)/C=N/CC(=O)O")
RIBOSE_GLY_AMADORI = _species("ribose-glycine-Amadori", "OCC(O)C(O)C(=O)CNCC(=O)O")
GLUCOSE_GLY_SCHIFF = _species("glucose-glycine-Schiff-base", "OCC(O)C(O)C(O)C(O)/C=N/CC(=O)O")
GLUCOSE_GLY_AMADORI = _species("glucose-glycine-Amadori", "OCC(O)C(O)C(O)C(=O)CNCC(=O)O")
DEOXYOSONE_3 = _species("3-deoxyosone", "O=CC(=O)CC(O)CO")
GLUCOSE_DEOXYOSONE_3 = _species("glucose-3-deoxyosone", "O=CC(=O)CC(O)C(O)CO")
HMF = _species("HMF", "OCC1=CC=C(C=O)O1")
HEXANAL = _species("hexanal", "CCCCCC=O")
HEXANAL_GLY_SCHIFF = _species("hexanal-glycine-Schiff-base", "CCCCC/C=N/CC(=O)O")
HEXANAL_LYS_SCHIFF = _species("hexanal-lysine-Schiff-base", "CCCCC/C=N/CCCCC(N)C(=O)O")
AMINOACETONE = _species("aminoacetone", "CC(=O)CN")
HYDROGEN = _species("hydrogen", "[HH]")


def _step(reactants: list[Species], products: list[Species], reaction_family: str) -> ElementaryStep:
    return ElementaryStep(reactants=reactants, products=products, reaction_family=reaction_family)


PATHWAYS = {
    "A_Core_Maillard_Ribose_Gly": [
        _step([RIBOSE, GLYCINE], [RIBOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([RIBOSE_GLY_SCHIFF], [RIBOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([RIBOSE_GLY_AMADORI], [DEOXYOSONE_3, GLYCINE], "Enolisation"),
        _step([DEOXYOSONE_3], [FURFURAL, WATER, WATER], "Sugar_Dehydration"),
    ],
    "A_Core_Maillard_Glucose_Gly": [
        _step([GLUCOSE, GLYCINE], [GLUCOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([GLUCOSE_GLY_SCHIFF], [GLUCOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([GLUCOSE_GLY_AMADORI], [GLUCOSE_DEOXYOSONE_3, GLYCINE], "Enolisation"),
        _step([GLUCOSE_DEOXYOSONE_3], [HMF, WATER, WATER], "Sugar_Dehydration"),
    ],
    "B_Strecker_Leu": [
        _step([PYRUVALDEHYDE, LEUCINE], [METHYLBUTANAL_3, AMINOACETONE, CO2], "Strecker_Degradation"),
    ],
    "C_S_Maillard_FFT": [
        _step([CYSTEINE, WATER], [ACETALDEHYDE, H2S, NH3, CO2], "Cysteine_Degradation"),
        _step([RIBOSE], [FURFURAL, WATER, WATER, WATER], "Sugar_Dehydration"),
        _step([FURFURAL, H2S, HYDROGEN], [FFT, WATER], "Thiol_Addition"),
    ],
    "D_Offflavour_Trapping_Gly": [
        _step([HEXANAL, GLYCINE], [HEXANAL_GLY_SCHIFF, WATER], "Lipid_Schiff_Base"),
    ],
    "D_Offflavour_Trapping_Lys": [
        _step([HEXANAL, LYSINE], [HEXANAL_LYS_SCHIFF, WATER], "Lipid_Schiff_Base"),
    ],
    "E_DHA_Competition": [
        _step([CYSTEINE], [DHA, H2S], "Beta_Elimination"),
        _step([DHA, LYSINE], [LAL], "DHA_Crosslinking"),
    ],
}


PATHWAY_METADATA = {
    "A_Core_Maillard_Ribose_Gly": {
        "target": FURFURAL,
        "type": "desirable",
        "consumes": ["D-ribose", "glycine"],
        "toxicity_flag": None,
    },
    "A_Core_Maillard_Glucose_Gly": {
        "target": HMF,
        "type": "desirable",
        "consumes": ["D-glucose", "glycine"],
        "toxicity_flag": "5-Hydroxymethylfurfural (HMF)",
    },
    "B_Strecker_Leu": {
        "target": METHYLBUTANAL_3,
        "type": "desirable",
        "consumes": ["pyruvaldehyde", "L-leucine"],
        "toxicity_flag": None,
    },
    "C_S_Maillard_FFT": {
        "target": FFT,
        "type": "desirable",
        "consumes": ["L-cysteine", "D-ribose"],
        "toxicity_flag": None,
    },
    "D_Offflavour_Trapping_Gly": {
        "target": HEXANAL_GLY_SCHIFF,
        "type": "masking",
        "consumes": ["hexanal", "glycine"],
        "toxicity_flag": None,
    },
    "D_Offflavour_Trapping_Lys": {
        "target": HEXANAL_LYS_SCHIFF,
        "type": "masking",
        "consumes": ["hexanal", "L-lysine"],
        "toxicity_flag": None,
    },
    "E_DHA_Competition": {
        "target": LAL,
        "type": "competing",
        "consumes": ["L-cysteine", "L-lysine"],
        "toxicity_flag": "Lysinoalanine (LAL)",
    },
}
