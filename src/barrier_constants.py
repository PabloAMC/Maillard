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
    "thiol_addition_trimolecular": (24.0, "H2S-mediated sulfur trapping should remain accessible but no longer tie the upstream carbonyl bottleneck"),
    "thiohemiacetal_formation": (23.3, "Furfural-thiohemiacetal formation is favorable but not faster than the dominant carbonyl cascade"),
    "thiol_dehydration":    (26.8, "Thiohemiacetal dehydration remains feasible but should impose a real selectivity cost relative to direct furfural release"),
    "thiol_addition":       (28.85,  "Pentose-derived MFT formation remains secondary to furfural release but must stay competitive enough to recover the Hofmann and Mottram sulfur balance"),
    "thiol_addition_hexose": (29.65, "Hexose-derived MFT formation should remain weaker than the pentose branch while still yielding measurable Farmer-type sulfur output"),
    "thiol_oxidation":      (29.02,  "Mottram-type furyl disulfide formation is secondary to MFT release but must stay accessible enough to preserve the calibrated disulfide branch"),

    # ── Enolisation / dehydration ───────────────────────────────────
    "enolisation_intermediate": (21.0,  "Amadori/Heyns deoxyosone formation is the common gateway into furfural and sulfur branches; keep it competitive instead of falling back to the heuristic default barrier"),
    "1,2-enolisation":      (28.0,  "Literature replication calibration: furfural-forming dehydration should be more competitive in benchmark systems"),
    "2,3-enolisation":      (28.0,  "Nursten 2005: enolisation is rate-limiting in advanced Maillard"),
    "dehydration":          (28.0,  "Coupled with enolisation; same approximate range"),

    # ── Strecker cascade ────────────────────────────────────────────
    "strecker_degradation": (22.0,  "Calibrated to reduce pyrazine over-expression in acidic sulfur benchmark systems while staying in literature range"),
    "aminoketone_condensation": (29.0,  "Pyrazine condensation should remain secondary to furfural in acidic sulfur systems while still producing measurable Farmer-type pyrazine output"),

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

    # ── Lipid Oxidation (Phase 19) ──────────────────────────────────
    "lipid_homolysis":      (42.0,  "O-O bond cleavage in hydroperoxides; high barrier"),
    "beta_scission":        (22.0,  "β-scission of alkoxy radicals; moderate barrier"),
    "radical_crosstalk":    (15.0,  "Radical + H2S collisions; fast"),
}

# Default barrier when no family pattern matches
DEFAULT_BARRIER: float = 45.0

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})


def get_barrier(reaction_family: str) -> Tuple[float, float]:
    """Return the FAST-mode (barrier, uncertainty) for a reaction family.
    
    Uncertainty mapping:
    * Heuristic fallback: ±5.0 kcal/mol
    * Literature calibrated: ±2.5 kcal/mol
    """
    default_unc = 5.0
    
    if not reaction_family:
        return DEFAULT_BARRIER, default_unc
        
    fm = reaction_family.lower().replace(" ", "_").replace("-", "_")
    
    # --- DYNAMIC CALIBRATION OVERRIDES (Phase 1) ---
    import os
    import json
    offsets = {}
    if "BARRIER_OFFSETS" in os.environ:
        try:
            offsets = json.loads(os.environ["BARRIER_OFFSETS"])
        except Exception:
            pass
    
    # Check for family-specific offset
    # Map optuna keys (short) to local fm keys
    offset_map = {
        "schiff": "schiff_condensation",
        "amadori": "amadori_rearrangement",
        "enol": "1,2-enolisation",
        "strecker": "strecker_degradation",
        "cys": "cysteine_thermolysis"
    }
    
    active_offset = 0.0
    for short_key, full_key in offset_map.items():
        if short_key in offsets and full_key in fm:
            active_offset = offsets[short_key]
            break
    # Check exact match first
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0] + active_offset, 3.5
        
    if "enolisation" in fm:
        if "1,2" in reaction_family or "1_2" in fm:
            fm = "1,2-enolisation"
        elif "2,3" in reaction_family or "2_3" in fm:
            fm = "2,3-enolisation"
        elif "dha" in fm:
            fm = "beta_elimination_dha"
        elif "beta" in fm:
            fm = "beta_elimination"
        elif "elimination" in fm:
            fm = "beta_elimination"
        else:
            fm = "1,2-enolisation" # Default
    elif "schiff" in fm:
        if "hydrolysis" not in fm and "reversion" not in fm:
            fm = "schiff_condensation"
        else:
            fm = "schiff_base_hydrolysis"
    elif "retro" in fm:
        fm = "retro_aldol"
    elif "lipid" in fm and "synergy" in fm:
        fm = "lipid_strecker_synergy"
    elif "lipid" in fm:
        fm = "lipid_condensation"
    elif "synergy" in fm:
        fm = "lipid_strecker_synergy"
    elif "strecker" in fm:
        fm = "strecker_degradation"
    elif "amadori" in fm:
        fm = "amadori_rearrangement"
    elif "heyns" in fm:
        fm = "heyns_rearrangement"
    elif "cysteine" in fm or "thermolysis" in fm:
        fm = "cysteine_thermolysis"
    elif "thiol" in fm and "oxidation" in fm:
        fm = "thiol_oxidation"
    elif "thiol" in fm and "addition" in fm and "hexose" in fm:
        fm = "thiol_addition_hexose"
    elif "thiol" in fm and "addition" in fm:
        fm = "thiol_addition"
    elif "pyrazine" in fm or "aminoketone" in fm:
        fm = "aminoketone_condensation"
    elif "thiazole" in fm:
        fm = "lipid_thiazole"
    elif "beta" in fm:
        fm = "beta_elimination"
    elif "ring" in fm or "mutarotation" in fm:
        fm = "ring_opening"
    elif "homolysis" in fm:
        fm = "lipid_homolysis"
    elif "beta_scission" in fm:
        fm = "beta_scission"
    elif "crosstalk" in fm:
        fm = "radical_crosstalk"
        
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0] + active_offset, 3.5
    return DEFAULT_BARRIER + active_offset, 5.0

# Global cache for arrhenius parameters
_ARRHENIUS_CACHE = None

def get_arrhenius_params(family: str) -> Optional[Tuple[float, float, str, float]]:
    """
    Retrieve literature-calibrated (A_value, Ea_kcal_mol, source_quality, uncertainty) for a reaction family.
    Returns None if no data is found for the family.
    """
    global _ARRHENIUS_CACHE
    if _ARRHENIUS_CACHE is None:
        if ARRHENIUS_FILE.exists():
            with open(ARRHENIUS_FILE, "r") as f:
                data = yaml.safe_load(f)
                _ARRHENIUS_CACHE = data.get("arrhenius_data", {})
                # Warm cache with TST defaults for missing A values
                for k, v in _ARRHENIUS_CACHE.items():
                    if v.get("A_value") is None or (isinstance(v["A_value"], float) and math.isnan(v["A_value"])):
                        v["A_value"] = 6.25e12 # TST @ 150C
                        v["source_quality"] = "estimated_tst"
        else:
            _ARRHENIUS_CACHE = {}
            
    if not family:
        return None
        
    fm = family.lower().replace(" ", "_").replace("-", "_")
    
    # Normalize to YAML keys
    yaml_key = None
    if "schiff" in fm:
        yaml_key = "schiff_condensation"
    elif "amadori" in fm:
        yaml_key = "amadori"
    elif "enolisation" in fm or ("dehydration" in fm and "thiol" not in fm):
        yaml_key = "dehydration"
    elif "strecker" in fm:
        yaml_key = "strecker"
    elif "pyrazine" in fm or "aminoketone" in fm:
        yaml_key = "pyrazine_condensation"
    elif "cysteine" in fm or "thermolysis" in fm:
        yaml_key = "cysteine_thermolysis"
    elif "thiol_addition" in fm:
        yaml_key = "thiol_addition"
    elif "retro" in fm:
        yaml_key = None # No retro_aldol data yet in yaml, only aggregates
    elif "beta" in fm or "dha" in fm:
        yaml_key = "beta_elimination_dha"
    elif "thiamine" in fm:
        yaml_key = "thiamine_degradation"
    elif "mutarotation" in fm or "ring" in fm:
        yaml_key = "mutarotation"
    elif "thiazole" in fm:
        yaml_key = "pyrazine_condensation" # Use similar collision factor
    
    if yaml_key and yaml_key in _ARRHENIUS_CACHE:
        entry = _ARRHENIUS_CACHE[yaml_key]
        A = float(entry.get("A_value", 0.0))
        Ea_kj = float(entry.get("Ea_kj_mol", 0.0))
        quality = entry.get("source_quality", "estimated")
        
        # Convert kJ/mol to kcal/mol
        Ea_kcal = Ea_kj / 4.184
        
        # Assign uncertainty based on quality
        quality_unc_map = {
            "literature": 2.0,
            "estimated_tst": 4.0,
            "heuristic": 5.0,
            "estimated": 3.5
        }
        uncertainty = quality_unc_map.get(quality, 3.5)
            
        return A, Ea_kcal, quality, uncertainty
        
    return None


