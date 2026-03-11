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

import yaml
import math
from pathlib import Path
from typing import Dict, Tuple, Optional

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
ARRHENIUS_FILE = ROOT / "data" / "lit" / "arrhenius_params.yml"

# Exact Mapping: normalized reaction family name → barrier in kcal/mol.
# Replaces fragile substring matching.
FAST_BARRIERS: Dict[str, Tuple[float, str]] = {
    # ── Sugar prerequisite ──────────────────────────────────────────
    "mutarotation":         ( 5.0,  "Ring opening is near-barrierless (hemiacetal ⇌ open-chain)"),
    "ring_opening":         ( 5.0,  "Ring opening is near-barrierless (hemiacetal ⇌ open-chain)"),

    # ── Core Maillard cascade ───────────────────────────────────────
    "schiff_condensation":  (15.0,  "Yaylayan 1994: condensation ΔG‡ ≈ 12–20 kcal/mol; midpoint"),
    "schiff_base_hydrolysis":(8.0,  "Schiff base reversion; fast"),
    "amadori_rearrangement":(23.0,  "Martins 2003: 1,2-proton shift ΔG‡ ≈ 20–28; midpoint"),
    "heyns_rearrangement":  (24.0,  "Ketose analogue of Amadori; slightly higher barrier"),

    # ── Sulfur pathways ─────────────────────────────────────────────
    "cysteine_thermolysis": (30.0,  "Wedzicha 1984: thermolysis ΔG‡ ≈ 20–30; upper range"),
    "thiol_addition_trimolecular": (15.0, "Trimolecular H2S-mediated thiolation; collision-limited"),
    "thiohemiacetal_formation": (18.0, "Formation of furfural-thiohemiacetal; fast nucleophilic addition"),
    "thiol_dehydration":    (15.0, "Dehydration of thiohemiacetal to furfurylthiol; very fast"),
    "thiol_addition":       (15.0,  "Furfural + H₂S thiol addition is fast; literature 10–18"),

    # ── Enolisation / dehydration ───────────────────────────────────
    "1,2-enolisation":      (28.0,  "Nursten 2005: enolisation is rate-limiting in advanced Maillard"),
    "2,3-enolisation":      (28.0,  "Nursten 2005: enolisation is rate-limiting in advanced Maillard"),
    "dehydration":          (28.0,  "Coupled with enolisation; same approximate range"),

    # ── Strecker cascade ────────────────────────────────────────────
    "strecker_degradation": (20.0,  "Hofmann 2000: decarboxylation ΔG‡ ≈ 18–25; lowered from 22"),
    "aminoketone_condensation": (22.0,  "Pyrazine dimerisation follows Strecker; slightly higher"),

    # ── Retro-aldol ─────────────────────────────────────────────────
    "retro_aldol":          (32.0,  "Hodge 1953: C-C bond cleavage is high-barrier; softened from 35"),
    "lipid_thiazole":       (20.0,  "Lipid thiazole condensation; comparable to Strecker"),

    # ── DHA / β-elimination ─────────────────────────────────────────
    "beta_elimination":     (18.0,  "β-elimination of Ser/Cys; moderate barrier"),
    "dha_crosslinking":     (18.0,  "DHA crosslinking with Lys; same range as β-elim"),

    # ── Lipid trapping & Synergy ──────────────────────────────────────
    "lipid_condensation":   (14.0,  "Lipid aldehyde Schiff base trapping; fast condensation"),
    "lipid_strecker_synergy": (18.0,  "Lipid-Maillard synergy (alkylthiazole/pyrazine) is highly favourable"),

    # ── Thermal / additive degradation ──────────────────────────────
    "thiamine_degradation": (25.0,  "Thiamine thermolysis; moderate barrier"),
    "additive_degradation": (25.0,  "Generic additive degradation"),
    "glutathione_cleavage": (22.0,  "GSH peptide bond cleavage"),
}

# Default barrier when no family pattern matches
DEFAULT_BARRIER: float = 40.0

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})


def get_barrier(reaction_family: str) -> float:
    """Return the FAST-mode barrier for a reaction family string.

    Performs exact matching against the normalized ``FAST_BARRIERS`` dict.
    Returns ``DEFAULT_BARRIER`` if no pattern matches.
    """
    if not reaction_family:
        return DEFAULT_BARRIER
        
    fm = reaction_family.lower().replace(" ", "_").replace("-", "_")
    
    # Normalize some common variants
    # Check exact match first
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0]
        
    if "enolisation" in fm:
        if "1,2" in reaction_family or "1_2" in fm: fm = "1,2-enolisation"
        elif "2,3" in reaction_family or "2_3" in fm: fm = "2,3-enolisation"
        else: fm = "1,2-enolisation" # Default
    elif "schiff" in fm:
        fm = "schiff_condensation" if "hydrolysis" not in fm and "reversion" not in fm else "schiff_base_hydrolysis"
    elif "retro" in fm: fm = "retro_aldol"
    elif "lipid" in fm and "synergy" in fm: fm = "lipid_strecker_synergy"
    elif "lipid" in fm: fm = "lipid_condensation"
    elif "synergy" in fm: fm = "lipid_strecker_synergy"
    elif "strecker" in fm: fm = "strecker_degradation"
    elif "amadori" in fm: fm = "amadori_rearrangement"
    elif "heyns" in fm: fm = "heyns_rearrangement"
    elif "cysteine" in fm or "thermolysis" in fm: fm = "cysteine_thermolysis"
    elif "thiol" in fm and "addition" in fm: fm = "thiol_addition"
    elif "pyrazine" in fm or "aminoketone" in fm: fm = "aminoketone_condensation"
    elif "thiazole" in fm: fm = "lipid_thiazole"
    elif "beta" in fm: fm = "beta_elimination"
    elif "ring" in fm or "mutarotation" in fm: fm = "ring_opening"
        
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0]
    return DEFAULT_BARRIER

# Global cache for arrhenius parameters
_ARRHENIUS_CACHE = None

def get_arrhenius_params(family: str) -> Optional[Tuple[float, float]]:
    """
    Retrieve literature-calibrated (A_value, Ea_kcal_mol) for a reaction family.
    Returns None if no data is found for the family.
    A_value is returned in 1/s or L/mol.s.
    Ea is returned in kcal/mol for consistency with our kinetics engine.
    """
    global _ARRHENIUS_CACHE
    if _ARRHENIUS_CACHE is None:
        if ARRHENIUS_FILE.exists():
            with open(ARRHENIUS_FILE, "r") as f:
                data = yaml.safe_load(f)
                _ARRHENIUS_CACHE = data.get("arrhenius_data", {})
        else:
            _ARRHENIUS_CACHE = {}
            
    if not family:
        return None
        
    fm = family.lower().replace(" ", "_").replace("-", "_")
    
    # Normalize to YAML keys
    yaml_key = None
    if "schiff" in fm: yaml_key = "schiff_condensation"
    elif "amadori" in fm: yaml_key = "amadori"
    elif "enolisation" in fm or ("dehydration" in fm and "thiol" not in fm): yaml_key = "dehydration"
    elif "strecker" in fm: yaml_key = "strecker"
    elif "pyrazine" in fm or "aminoketone" in fm: yaml_key = "pyrazine_condensation"
    elif "cysteine" in fm or "thermolysis" in fm: yaml_key = "cysteine_thermolysis"
    elif "thiol_addition" in fm: yaml_key = "thiol_addition"
    elif "retro" in fm: yaml_key = None # No retro_aldol data yet in yaml, only aggregates
    elif "beta" in fm or "dha" in fm: yaml_key = "beta_elimination_dha"
    elif "thiamine" in fm: yaml_key = "thiamine_degradation"
    elif "mutarotation" in fm or "ring" in fm: yaml_key = "mutarotation"
    
    if yaml_key and yaml_key in _ARRHENIUS_CACHE:
        entry = _ARRHENIUS_CACHE[yaml_key]
        A = float(entry.get("A_value", 0.0))
        Ea_kj = float(entry.get("Ea_kj_mol", 0.0))
        
        # Convert kJ/mol to kcal/mol
        Ea_kcal = Ea_kj / 4.184
        
        # If A_value is NaN or 0.0 (placeholder), return None so caller falls back to A=1e13
        if math.isnan(A) or A <= 0.0:
            return None
            
        return A, Ea_kcal
        
    return None


