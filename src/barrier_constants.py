"""
src/barrier_constants.py — Centralised FAST-mode heuristic barrier constants.

These values are approximate activation energies (kcal/mol) for each
Maillard reaction family, sourced from published DFT/experimental data
and cross-checked against GFN2-xTB NEB estimates.

They are used by both `inverse_design.py` and `run_pipeline.py` for the
instant FAST-mode rankings.  Update this single file when new data
is available — both call-sites import from here.

Sources
-------
* Yaylayan & Huyghues-Despointes 1994 (Schiff base condensation)
* Martins & van Boekel 2003 (Amadori kinetics)
* Hofmann & Schieberle 2000 (Strecker degradation)
* Wedzicha 1984 (Cysteine thermolysis)
* Hodge 1953; Nursten 2005 (overall Maillard kinetics)
* Maillard_meat.md, Maillard_Plant_based.md (project literature reviews)
"""

from typing import Dict, Tuple

# Mapping: reaction-family substring (lowercased) → barrier in kcal/mol.
# The substring is matched with `in` against step.reaction_family.lower().
# Order matters: first match wins.  More specific patterns come first.
#
# Format: { pattern: (barrier_kcal, source_note) }
#
FAST_BARRIERS: Dict[str, Tuple[float, str]] = {
    # ── Sugar prerequisite ──────────────────────────────────────────
    "ring":         ( 5.0,  "Ring opening is near-barrierless (hemiacetal ⇌ open-chain)"),

    # ── Core Maillard cascade ───────────────────────────────────────
    "schiff":       (15.0,  "Yaylayan 1994: condensation ΔG‡ ≈ 12–20 kcal/mol; midpoint"),
    "amadori":      (23.0,  "Martins 2003: 1,2-proton shift ΔG‡ ≈ 20–28; midpoint"),
    "heyns":        (24.0,  "Ketose analogue of Amadori; slightly higher barrier"),

    # ── Sulfur pathways ─────────────────────────────────────────────
    "cysteine":     (30.0,  "Wedzicha 1984: thermolysis ΔG‡ ≈ 20–30; upper range"),
    "thiol_addition_trimolecular": (15.0, "Trimolecular H2S-mediated thiolation; collision-limited"),
    "thiohemiacetal": (18.0, "Formation of furfural-thiohemiacetal; fast nucleophilic addition"),
    "thiol_dehydration": (15.0, "Dehydration of thiohemiacetal to furfurylthiol; very fast"),
    "thiol":        (15.0,  "Furfural + H₂S thiol addition is fast; literature 10–18"),

    # ── Enolisation / dehydration ───────────────────────────────────
    "enolisation":  (28.0,  "Nursten 2005: enolisation is rate-limiting in advanced Maillard"),
    "dehydration":  (28.0,  "Coupled with enolisation; same approximate range"),

    # ── Strecker cascade ────────────────────────────────────────────
    "strecker":     (20.0,  "Hofmann 2000: decarboxylation ΔG‡ ≈ 18–25; lowered from 22"),
    "aminoketone":  (22.0,  "Pyrazine dimerisation follows Strecker; slightly higher"),

    # ── Retro-aldol ─────────────────────────────────────────────────
    "retro":        (32.0,  "Hodge 1953: C-C bond cleavage is high-barrier; softened from 35"),
    "thiazole":     (20.0,  "Lipid thiazole condensation; comparable to Strecker"),

    # ── DHA / β-elimination ─────────────────────────────────────────
    "beta":         (18.0,  "β-elimination of Ser/Cys; moderate barrier"),
    "dha":          (18.0,  "DHA crosslinking with Lys; same range as β-elim"),

    # ── Lipid trapping & Synergy ──────────────────────────────────────
    "lipid":        (14.0,  "Lipid aldehyde Schiff base trapping; fast condensation"),
    "synergy":      (18.0,  "Lipid-Maillard synergy (alkylthiazole/pyrazine) is highly favourable"),

    # ── Thermal / additive degradation ──────────────────────────────
    "thermal":      (25.0,  "Thiamine thermolysis; moderate barrier"),
    "additive":     (25.0,  "Generic additive degradation"),
    "glutathione":  (22.0,  "GSH peptide bond cleavage"),
}

# Default barrier when no family pattern matches
DEFAULT_BARRIER: float = 40.0

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})


def get_barrier(reaction_family: str) -> float:
    """Return the FAST-mode barrier for a reaction family string.

    Performs case-insensitive substring matching against the ordered
    ``FAST_BARRIERS`` dict.  Returns ``DEFAULT_BARRIER`` if no pattern
    matches.
    """
    if not reaction_family:
        return DEFAULT_BARRIER
    fm = reaction_family.lower()
    for pattern, (barrier, _note) in FAST_BARRIERS.items():
        if pattern in fm:
            return barrier
    return DEFAULT_BARRIER

