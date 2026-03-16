"""
Correction factors for amino acid reactivity in protein matrices
vs. free amino acid model systems.

These are empirical correction factors derived from literature comparing
model system predictions to protein matrix experiments. They should be
recalibrated as you generate your own experimental data.

Current legume-matrix anchoring inside the repo:
- pea: Prigent 2024 TNBS/DTNB + Schneider 2023 OPA envelope
- soy: repo literature synthesis around soy glycinin/beta-conglycinin burial,
  sulfur limitation, extrusion-driven accessibility changes, and soy
  protein-polysaccharide conjugate trapping documented in
  docs/Maillard_Plant_based.md and
  docs/Elicit - Maillard Pathways in Plant-Based Cooking - Report.md
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
    lysine_accessibility_native: float
    lysine_accessibility_denatured: float
    cysteine_accessibility_native: float
    cysteine_accessibility_denatured: float
    volatile_retention: float       # fraction escaping matrix (rest is bound)
    volatile_retention_native: float
    volatile_retention_denatured: float
    source: str


@dataclass(frozen=True)
class EffectiveMatrixCorrection:
    protein_type: ProteinType
    lysine_accessibility: float
    cysteine_accessibility: float
    volatile_retention: float
    source: str


_LYSINE_IDENTIFIERS = {
    "lysine",
    "l-lysine",
    "nccccc(n)c(=o)o",
}

_CYSTEINE_IDENTIFIERS = {
    "cysteine",
    "l-cysteine",
    "nc(cs)c(=o)o",
}

MATRIX_CORRECTIONS = {
    ProteinType.FREE_AMINO_ACID: MatrixCorrection(
        protein_type=ProteinType.FREE_AMINO_ACID,
        lysine_accessibility=1.0,
        cysteine_accessibility=1.0,
        lysine_accessibility_native=1.0,
        lysine_accessibility_denatured=1.0,
        cysteine_accessibility_native=1.0,
        cysteine_accessibility_denatured=1.0,
        volatile_retention=1.0,
        volatile_retention_native=1.0,
        volatile_retention_denatured=1.0,
        source="model system, no correction"
    ),
    ProteinType.PEA_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.PEA_ISOLATE,
        # TNBS/OPA envelope in the repo suggests commercial pea isolates sit well
        # below free-AA accessibility even after wet-heat treatment.
        lysine_accessibility=0.39,
        cysteine_accessibility=0.04,
        lysine_accessibility_native=0.32,
        lysine_accessibility_denatured=0.46,
        cysteine_accessibility_native=0.01,
        cysteine_accessibility_denatured=0.07,
        volatile_retention=0.50,
        volatile_retention_native=0.42,
        volatile_retention_denatured=0.58,
        source="Prigent 2024 TNBS/DTNB + Schneider 2023 OPA envelope; PPI free SH near zero even after processing"
    ),
    ProteinType.PEA_CONCENTRATE: MatrixCorrection(
        protein_type=ProteinType.PEA_CONCENTRATE,
        lysine_accessibility=0.29,
        cysteine_accessibility=0.03,
        lysine_accessibility_native=0.22,
        lysine_accessibility_denatured=0.36,
        cysteine_accessibility_native=0.005,
        cysteine_accessibility_denatured=0.055,
        volatile_retention=0.35,
        volatile_retention_native=0.27,
        volatile_retention_denatured=0.43,
        source="Pea-isolate accessibility envelope with added fiber burial penalty"
    ),
    ProteinType.SOY_CONCENTRATE: MatrixCorrection(
        protein_type=ProteinType.SOY_CONCENTRATE,
        lysine_accessibility=0.35,
        cysteine_accessibility=0.16,
        lysine_accessibility_native=0.30,
        lysine_accessibility_denatured=0.40,
        cysteine_accessibility_native=0.11,
        cysteine_accessibility_denatured=0.21,
        volatile_retention=0.45,
        volatile_retention_native=0.37,
        volatile_retention_denatured=0.53,
        source=(
            "Soy concentrate envelope derived from soy-isolate globulin burial "
            "(glycinin/beta-conglycinin) with added fiber/concentrate penalty; "
            "repo anchors: Kutzli 2020 / Naik 2021 synthesis plus soy-polysaccharide "
            "conjugation trapping notes in docs/Maillard_Plant_based.md"
        )
    ),
    ProteinType.SOY_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.SOY_ISOLATE,
        # Soy is kept above pea because the repo literature synthesis treats soy
        # isolates as lysine-rich globulins whose extrusion state can expose more
        # reactive sites, while still remaining well below free-AA behavior and
        # sulfur-limited relative to true meat-like precursor systems.
        lysine_accessibility=0.45,
        cysteine_accessibility=0.20,
        lysine_accessibility_native=0.39,
        lysine_accessibility_denatured=0.51,
        cysteine_accessibility_native=0.14,
        cysteine_accessibility_denatured=0.26,
        volatile_retention=0.55,
        volatile_retention_native=0.47,
        volatile_retention_denatured=0.63,
        source=(
            "Soy isolate envelope anchored to repo literature on glycinin/beta-conglycinin "
            "burial, lysine-rich but sulfur-poor plant-protein composition, and extrusion- "
            "altered reactive-site accessibility; see docs/Maillard_Plant_based.md and "
            "docs/Elicit - Maillard Pathways in Plant-Based Cooking - Report.md "
            "(Kutzli 2020 / Naik 2021 synthesis)"
        )
    ),
    ProteinType.MYCOPROTEIN: MatrixCorrection(
        protein_type=ProteinType.MYCOPROTEIN,
        lysine_accessibility=0.50,
        cysteine_accessibility=0.24,
        lysine_accessibility_native=0.42,
        lysine_accessibility_denatured=0.58,
        cysteine_accessibility_native=0.18,
        cysteine_accessibility_denatured=0.30,
        volatile_retention=0.60,
        volatile_retention_native=0.54,
        volatile_retention_denatured=0.66,
        source="Conservative mycoprotein estimate pending dedicated calibration"
    ),
}


def clamp_denaturation_state(denaturation_state: float) -> float:
    return max(0.0, min(1.0, float(denaturation_state)))


def resolve_matrix_correction(
    protein_type: ProteinType,
    denaturation_state: float = 0.5,
) -> EffectiveMatrixCorrection:
    corr = MATRIX_CORRECTIONS.get(
        protein_type, MATRIX_CORRECTIONS[ProteinType.FREE_AMINO_ACID]
    )
    state = clamp_denaturation_state(denaturation_state)
    lys_factor = corr.lysine_accessibility_native + state * (
        corr.lysine_accessibility_denatured - corr.lysine_accessibility_native
    )
    cys_factor = corr.cysteine_accessibility_native + state * (
        corr.cysteine_accessibility_denatured - corr.cysteine_accessibility_native
    )
    volatile_retention = corr.volatile_retention_native + state * (
        corr.volatile_retention_denatured - corr.volatile_retention_native
    )
    return EffectiveMatrixCorrection(
        protein_type=protein_type,
        lysine_accessibility=lys_factor,
        cysteine_accessibility=cys_factor,
        volatile_retention=volatile_retention,
        source=corr.source,
    )


def _classify_reactive_aa_key(key: str) -> str | None:
    normalized = key.strip().lower()
    if normalized in _LYSINE_IDENTIFIERS:
        return "lysine"
    if normalized in _CYSTEINE_IDENTIFIERS:
        return "cysteine"
    return None

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
    corr = resolve_matrix_correction(
        protein_type,
        denaturation_state=denaturation_state,
    )

    corrected_aa = {}
    for aa, conc in reactive_amino_acids.items():
        aa_class = _classify_reactive_aa_key(aa)
        if aa_class == "lysine":
            corrected_aa[aa] = conc * corr.lysine_accessibility
        elif aa_class == "cysteine":
            corrected_aa[aa] = conc * corr.cysteine_accessibility
        else:
            corrected_aa[aa] = conc * (corr.lysine_accessibility + corr.cysteine_accessibility) / 2.0

    corrected_volatiles = {
        compound: conc * corr.volatile_retention
        for compound, conc in predicted_concentrations.items()
    }
    return corrected_volatiles, corrected_aa
