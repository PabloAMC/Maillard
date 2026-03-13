"""
Correction factors for amino acid reactivity in protein matrices
vs. free amino acid model systems.

These are empirical correction factors derived from literature comparing
model system predictions to protein matrix experiments. They should be
recalibrated as you generate your own experimental data.
"""

from dataclasses import dataclass
from enum import Enum

class ProteinType(Enum):
    FREE_AMINO_ACID = "free"
    PEA_CONCENTRATE = "pea_conc"    # ~60% protein, fibrous matrix
    PEA_ISOLATE = "pea_iso"         # ~85% protein
    SOY_CONCENTRATE = "soy_conc"
    SOY_ISOLATE = "soy_iso"
    MYCOPROTEIN = "myco"

@dataclass
class MatrixCorrection:
    protein_type: ProteinType
    lysine_accessibility: float     # fraction of total lysine reactive
    cysteine_accessibility: float   # usually lower due to disulfide burial
    volatile_retention: float       # fraction escaping matrix (rest is bound)
    source: str

MATRIX_CORRECTIONS = {
    ProteinType.FREE_AMINO_ACID: MatrixCorrection(
        protein_type=ProteinType.FREE_AMINO_ACID,
        lysine_accessibility=1.0,
        cysteine_accessibility=1.0,
        volatile_retention=1.0,
        source="model system, no correction"
    ),
    ProteinType.PEA_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.PEA_ISOLATE,
        # Maillard et al. 2017 JAFC: ~35-45% reactive lysine in pea isolate
        lysine_accessibility=0.40,
        # Cysteine more buried in legume storage proteins (legumin/vicilin)
        cysteine_accessibility=0.25,
        # Fat/protein binding reduces headspace by ~40-60%
        volatile_retention=0.50,
        source="Maillard 2017 JAFC + Keller 2020 Food Chem estimates"
    ),
    ProteinType.PEA_CONCENTRATE: MatrixCorrection(
        protein_type=ProteinType.PEA_CONCENTRATE,
        lysine_accessibility=0.30,
        cysteine_accessibility=0.20,
        volatile_retention=0.35,
        source="Estimated from pea isolate values with fiber correction"
    ),
    ProteinType.SOY_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.SOY_ISOLATE,
        lysine_accessibility=0.45,
        cysteine_accessibility=0.30,
        volatile_retention=0.55,
        source="Aaslyng 2009 + Raza 2021 estimates"
    ),
}

def apply_matrix_correction(
    predicted_concentrations: dict[str, float],
    reactive_amino_acids: dict[str, float],
    protein_type: ProteinType,
    denaturation_state: float = 0.5,  # 0=native, 1=fully denatured
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Scale predicted volatile concentrations and reactive AA concentrations
    by matrix accessibility factors.

    denaturation_state: extrusion/heating increases accessibility.
    1.0 (fully denatured) approaches the free AA model system.
    """
    corr = MATRIX_CORRECTIONS.get(
        protein_type, MATRIX_CORRECTIONS[ProteinType.FREE_AMINO_ACID]
    )
    lys_factor = corr.lysine_accessibility + denaturation_state * (
        1.0 - corr.lysine_accessibility
    )
    cys_factor = corr.cysteine_accessibility + denaturation_state * (
        1.0 - corr.cysteine_accessibility
    )

    corrected_aa = {}
    for aa, conc in reactive_amino_acids.items():
        if aa.lower() == "lysine":
            corrected_aa[aa] = conc * lys_factor
        elif aa.lower() == "cysteine":
            corrected_aa[aa] = conc * cys_factor
        else:
            corrected_aa[aa] = conc * (lys_factor + cys_factor) / 2.0

    corrected_volatiles = {
        compound: conc * corr.volatile_retention
        for compound, conc in predicted_concentrations.items()
    }
    return corrected_volatiles, corrected_aa
