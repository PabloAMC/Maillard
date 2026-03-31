"""
src/barrier_constants.py — Centralised FAST-mode heuristic barrier constants.

These values are approximate activation energies (kcal/mol) for each
Maillard reaction family, sourced from published DFT/experimental data
and cross-checked against GFN2-xTB NEB estimates.

They are used by both `pipeline.py` and `run_pipeline.py` for the
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

import json
import yaml
import math
from pathlib import Path
from typing import Dict, Tuple, Optional

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
ARRHENIUS_FILE = ROOT / "data" / "lit" / "arrhenius_params.yml"
REFINEMENT_PATCH_FILE = ROOT / "data" / "lit" / "refinement_surrogate_patches.json"

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

    # ── MACE-OFF24 Derived Barriers (Families 11-16) ───────────────
    # Retained for reference; the routing table below prefers higher tiers.
    "hexanal_radical_quench_mlp_derived": (7.89, "MACE-OFF24 screening estimate"),
    "mft_protein_noncovalent_mlp_derived": (25.69, "MACE-OFF24 screening estimate"),
    "quinone_cys_michael_mlp_derived": (20.33, "MACE-OFF24 screening estimate"),
    "quinone_lys_schiff_mlp_derived": (7.41, "MACE-OFF24 screening estimate"),
    "aa_ring_open_dicarbonyl_mlp_derived": (18.10, "MACE-OFF24 screening estimate"),
    "pe_schiff_base_mlp_derived": (15.20, "MACE-OFF24 screening estimate"),
    "pe_amadori_mlp_derived": (8.94, "MACE-OFF24 screening estimate"),
    "melanoidin_radical_trapping_mlp_derived": (26.29, "MACE-OFF24 screening estimate"),
    "lysinoalanine_crosslink_mlp_derived": (31.75, "MACE-OFF24 screening estimate"),
    "furosine_amadori_hydrolysis_mlp_derived": (5.33, "MACE-OFF24 screening estimate"),

    # ── xTB-derived Barriers (Step C4/C5 promotion — xtb_derived_gfn2 tier) ─
    # These 4 reactions have xTB path searches in data/geometries/xtb_inputs/;
    # the _dft entries will be written by ingest_dft_c4_c5_results.py once
    # the wB97M-V//r2SCAN-3c jobs complete.  Until then these values represent
    # the best available GFN2-xTB guided estimates from the path search.
    # Uncertainty floor: ±10 kcal/mol (factor 3–8 in concentration).
    "hexanal_radical_quench_xtb_derived": (
        7.58,
        "GFN2-xTB NEB path search (xtb_derived_gfn2); ±10 kcal/mol; use for ranking only"
    ),
    "quinone_cys_michael_xtb_derived": (
        17.84,
        "GFN2-xTB NEB path search (xtb_derived_gfn2); ±10 kcal/mol; use for ranking only"
    ),
    "aa_ring_open_dicarbonyl_xtb_derived": (
        7.58,   # Ea=31.7 kJ/mol from HCW = 7.577 kcal/mol — literature anchor
        "HCW literature Ea=31.7 kJ/mol (literature_derived_transfer); ±5 kcal/mol; bounded calibration"
    ),
    "lysinoalanine_crosslink_xtb_derived": (
        28.90,
        "GFN2-xTB NEB path search (xtb_derived_gfn2); ±10 kcal/mol; use for ranking only"
    ),

    # ── Literature-Derived Barriers (Step C4/C5 — literature_derived_transfer) ─
    # pe_schiff_base and pe_amadori have well-established literature Ea values
    # from SLR 15 (phospholipid-amine Maillard).  These are used directly and
    # are better than the MACE estimates.  The DFT refinement will validate them.
    # Uncertainty: ±5 kcal/mol on transfer from model compound literature.
    "pe_schiff_base_lit_derived": (
        22.21,  # 92.9 kJ/mol ÷ 4.184
        "SLR 15 literature Ea=92.9 kJ/mol (literature_derived_transfer); ±5 kcal/mol; bounded calibration"
    ),
    "pe_amadori_lit_derived": (
        19.81,  # 82.9 kJ/mol ÷ 4.184
        "SLR 15 literature Ea=82.9 kJ/mol (literature_derived_transfer); ±5 kcal/mol; bounded calibration"
    ),

    # ── DFT-Anchor Placeholders (will be filled by ingest_dft_c4_c5_results.py) ─
    # These keys are pre-registered so the routing table works before the DFT
    # jobs complete.  Values are None-signalled via the special sentinel 99.0.
    "hexanal_radical_quench_dft": (
        99.0,  # sentinel: not yet available
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
    "quinone_cys_michael_dft": (
        99.0,
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
    "aa_ring_open_dicarbonyl_dft": (
        99.0,
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
    "pe_schiff_base_dft": (
        99.0,
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
    "pe_amadori_dft": (
        99.0,
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
    "lysinoalanine_crosslink_dft": (
        99.0,
        "selective_dft_anchor pending: run scripts/run_dft_c4_c5.py"
    ),
}

# ── Evidence-tier metadata for the 6 C4/C5 DFT targets ─────────────────────
# Provides machine-readable current_tier, target_tier, uncertainty_kj, and
# promotion_ceiling per reaction.  Updated by ingest_dft_c4_c5_results.py
# when wB97M-V//r2SCAN-3c jobs complete.
DFT_ANCHOR_METADATA: Dict[str, Dict[str, object]] = {
    "hexanal_radical_quench": {
        "current_tier": "xtb_derived_gfn2",
        "target_tier": "selective_dft_anchor",
        "slr_family": "11",
        "active_key": "hexanal_radical_quench_xtb_derived",
        "dft_key": "hexanal_radical_quench_dft",
        "uncertainty_kj": 42.0,
        "promotion_ceiling": "ranking_only",
        "honest_label": "GFN2-xTB estimate, ±factor 3–8 in concentration, ranking only",
    },
    "quinone_cys_michael": {
        "current_tier": "xtb_derived_gfn2",
        "target_tier": "selective_dft_anchor",
        "slr_family": "13",
        "active_key": "quinone_cys_michael_xtb_derived",
        "dft_key": "quinone_cys_michael_dft",
        "uncertainty_kj": 42.0,
        "promotion_ceiling": "bounded_calibration",
        "honest_label": "GFN2-xTB estimate, ±factor 3–8 in concentration, ranking only",
    },
    "aa_ring_open_dicarbonyl": {
        "current_tier": "literature_derived_transfer",
        "target_tier": "selective_dft_anchor",
        "slr_family": "14",
        "active_key": "aa_ring_open_dicarbonyl_xtb_derived",
        "dft_key": "aa_ring_open_dicarbonyl_dft",
        "uncertainty_kj": 20.0,
        "promotion_ceiling": "bounded_calibration",
        "honest_label": "HCW literature Ea=31.7 kJ/mol, ±factor 2–5; bounded calibration",
    },
    "pe_schiff_base": {
        "current_tier": "literature_derived_transfer",
        "target_tier": "selective_dft_anchor",
        "slr_family": "15",
        "active_key": "pe_schiff_base_lit_derived",
        "dft_key": "pe_schiff_base_dft",
        "uncertainty_kj": 20.9,
        "promotion_ceiling": "bounded_calibration",
        "honest_label": "SLR 15 literature Ea=92.9 kJ/mol, ±factor 2–5; bounded calibration",
    },
    "pe_amadori": {
        "current_tier": "literature_derived_transfer",
        "target_tier": "selective_dft_anchor",
        "slr_family": "15",
        "active_key": "pe_amadori_lit_derived",
        "dft_key": "pe_amadori_dft",
        "uncertainty_kj": 20.9,
        "promotion_ceiling": "bounded_calibration",
        "honest_label": "SLR 15 literature Ea=82.9 kJ/mol, ±factor 2–5; bounded calibration",
    },
    "lysinoalanine_crosslink": {
        "current_tier": "xtb_derived_gfn2",
        "target_tier": "selective_dft_anchor",
        "slr_family": "12",
        "active_key": "lysinoalanine_crosslink_xtb_derived",
        "dft_key": "lysinoalanine_crosslink_dft",
        "uncertainty_kj": 42.0,
        "promotion_ceiling": "ranking_only",
        "honest_label": "GFN2-xTB estimate, ±factor 3–8 in concentration, ranking only",
    },
}


def evidence_tier(reaction_key: str) -> str:
    """Return the current evidence tier for a C4/C5 reaction key, or 'unknown'."""
    meta = DFT_ANCHOR_METADATA.get(reaction_key)
    if meta is None:
        return "unknown"
    # If the DFT placeholder is no longer sentinel (99.0), it has been filled
    dft_key = str(meta.get("dft_key", ""))
    if dft_key and dft_key in FAST_BARRIERS:
        val = FAST_BARRIERS[dft_key][0]
        if val < 99.0:
            return "selective_dft_anchor"
    return str(meta.get("current_tier", "unknown"))


# Default barrier when no family pattern matches
DEFAULT_BARRIER: float = 45.0
DEFAULT_REFERENCE_PREEXPONENTIAL: float = 1e11
ARRHENIUS_R_KCAL: float = 0.001987

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})


def _normalize_family_key(reaction_family: Optional[str]) -> str:
    if not reaction_family:
        return ""
    return reaction_family.lower().replace(" ", "_").replace("-", "_")


def _canonical_fast_family(reaction_family: Optional[str]) -> str:
    fm = _normalize_family_key(reaction_family)
    if not fm:
        return ""
    if fm in FAST_BARRIERS:
        return fm

    if "enolisation" in fm:
        if "1,2" in str(reaction_family) or "1_2" in fm:
            return "1,2-enolisation"
        if "2,3" in str(reaction_family) or "2_3" in fm:
            return "2,3-enolisation"
        if "dha" in fm:
            return "beta_elimination"
        if "beta" in fm or "elimination" in fm:
            return "beta_elimination"
        return "1,2-enolisation"
    if "schiff" in fm:
        if "hydrolysis" not in fm and "reversion" not in fm:
            return "schiff_condensation"
        return "schiff_base_hydrolysis"
    if "retro" in fm:
        return "retro_aldol"
    if "lipid" in fm and "synergy" in fm:
        return "lipid_strecker_synergy"
    if "lipid" in fm:
        return "lipid_condensation"
    if "synergy" in fm:
        return "lipid_strecker_synergy"
    if "strecker" in fm:
        return "strecker_degradation"
    if "amadori" in fm:
        return "amadori_rearrangement"
    if "heyns" in fm:
        return "heyns_rearrangement"
    if "cysteine" in fm or "thermolysis" in fm:
        return "cysteine_thermolysis"
    if "thiol" in fm and "oxidation" in fm:
        return "thiol_oxidation"
    if "thiol" in fm and "addition" in fm and "hexose" in fm:
        return "thiol_addition_hexose"
    if "thiol" in fm and "addition" in fm:
        return "thiol_addition"
    if "pyrazine" in fm or "aminoketone" in fm:
        return "aminoketone_condensation"
    if "thiazole" in fm:
        return "lipid_thiazole"
    if "beta" in fm:
        return "beta_elimination"
    if "ring" in fm or "mutarotation" in fm:
        return "ring_opening"
    if "homolysis" in fm:
        return "lipid_homolysis"
    if "beta_scission" in fm:
        return "beta_scission"
    if "crosstalk" in fm:
        return "radical_crosstalk"
    return fm


def _arrhenius_yaml_key(family: Optional[str]) -> Optional[str]:
    canonical_family = _canonical_fast_family(family)
    if not canonical_family:
        return None
    yaml_key_map = {
        "schiff_condensation": "schiff_condensation",
        "schiff_base_hydrolysis": None,
        "amadori_rearrangement": "amadori",
        "heyns_rearrangement": "amadori",
        "1,2-enolisation": "enolisation",
        "2,3-enolisation": "enolisation",
        "enolisation_intermediate": "enolisation",
        "dehydration": "dehydration",
        "strecker_degradation": "strecker",
        "aminoketone_condensation": "pyrazine_condensation",
        "cysteine_thermolysis": "cysteine_thermolysis",
        "thiol_addition": "thiol_addition",
        "thiol_addition_hexose": "thiol_addition",
        "retro_aldol": "retro_aldol",
        "beta_elimination": "beta_elimination_dha",
        "thiamine_degradation": "thiamine_degradation",
        "ring_opening": "mutarotation",
        "mutarotation": "mutarotation",
        "lipid_thiazole": "pyrazine_condensation",
        
        # ── Families 11-16: tiered routing (prefer _dft > _lit/_xtb > _mlp) ──
        # The ingest_dft_c4_c5_results.py script writes the _dft YAML entries.
        # Until then, _lit_derived and _xtb_derived are used.
        "hexanal_radical_quench": "hexanal_radical_quench_xtb_derived",
        "mft_protein_noncovalent": "mft_protein_noncovalent_mlp_derived",
        "quinone_cys_michael": "quinone_cys_michael_xtb_derived",
        "quinone_lys_schiff": "quinone_lys_schiff_mlp_derived",
        "aa_ring_open_dicarbonyl": "aa_ring_open_dicarbonyl_xtb_derived",
        "pe_schiff_base": "pe_schiff_base_lit_derived",
        "pe_amadori": "pe_amadori_lit_derived",
        "melanoidin_radical_trapping": "melanoidin_radical_trapping_mlp_derived",
        "lysinoalanine_crosslink": "lysinoalanine_crosslink_xtb_derived",
        "furosine_amadori_hydrolysis": "furosine_amadori_hydrolysis_mlp_derived",
    }
    # If the family key is already a fully-qualified FAST_BARRIERS key, map it to itself
    if canonical_family in yaml_key_map.values():
        return canonical_family

    return yaml_key_map.get(canonical_family)


_REFINEMENT_PATCH_CACHE: Optional[Dict[str, float]] = None
_REFINEMENT_PATCH_MTIME: Optional[float] = None


def _load_refinement_surrogate_offsets() -> Dict[str, float]:
    global _REFINEMENT_PATCH_CACHE, _REFINEMENT_PATCH_MTIME
    if not REFINEMENT_PATCH_FILE.exists():
        _REFINEMENT_PATCH_CACHE = {}
        _REFINEMENT_PATCH_MTIME = None
        return {}

    mtime = REFINEMENT_PATCH_FILE.stat().st_mtime
    if _REFINEMENT_PATCH_CACHE is not None and _REFINEMENT_PATCH_MTIME == mtime:
        return dict(_REFINEMENT_PATCH_CACHE)

    try:
        with open(REFINEMENT_PATCH_FILE, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        _REFINEMENT_PATCH_CACHE = {}
        _REFINEMENT_PATCH_MTIME = mtime
        return {}

    offsets = payload.get("accepted_offsets", {}) if isinstance(payload, dict) else {}
    if not isinstance(offsets, dict):
        offsets = {}
    normalized = {str(key): float(value) for key, value in offsets.items()}
    _REFINEMENT_PATCH_CACHE = normalized
    _REFINEMENT_PATCH_MTIME = mtime
    return dict(normalized)


def get_barrier(reaction_family: str) -> Tuple[float, float]:
    """Return the FAST-mode (barrier, uncertainty) for a reaction family.
    
    Uncertainty mapping:
    * Heuristic fallback: ±5.0 kcal/mol
    * Literature calibrated: ±2.5 kcal/mol
    """
    default_unc = 5.0
    
    if not reaction_family:
        return DEFAULT_BARRIER, default_unc
        
    fm = _normalize_family_key(reaction_family)
    
    # --- DYNAMIC CALIBRATION OVERRIDES (Phase 1) ---
    import os
    offsets = _load_refinement_surrogate_offsets()
    if "BARRIER_OFFSETS" in os.environ:
        try:
            runtime_offsets = json.loads(os.environ["BARRIER_OFFSETS"])
            if isinstance(runtime_offsets, dict):
                offsets.update({str(key): float(value) for key, value in runtime_offsets.items()})
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
    if fm in offsets:
        active_offset = float(offsets[fm])
    for short_key, full_key in offset_map.items():
        if short_key in offsets and full_key in fm:
            active_offset = offsets[short_key]
            break
    # Check exact match first
    if fm in FAST_BARRIERS:
        return FAST_BARRIERS[fm][0] + active_offset, 3.5
    fm = _canonical_fast_family(reaction_family)

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
        
    yaml_key = _arrhenius_yaml_key(family)

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
            "estimated": 3.5,
            "mlp_screen_mace": 15.0, # MACE-OFF24 derived uncertainty floor
            "xtb_derived_gfn2": 10.0,
            "selective_dft_anchor": 3.0
        }
        # explicit YAML uncertainty trumps mapped quality uncertainty
        uncertainty = float(entry.get("uncertainty_kj", quality_unc_map.get(quality, 3.5)))
            
        return A, Ea_kcal, quality, uncertainty
        
    return None


def get_reference_pre_exponential(family: Optional[str] = None) -> float:
    params = get_arrhenius_params(family or "")
    if params is None:
        return DEFAULT_REFERENCE_PREEXPONENTIAL
    return float(params[0])


def arrhenius_rate_constant(
    barrier_kcal: float,
    temperature_kelvin: float,
    family: Optional[str] = None,
    multiplier: float = 1.0,
) -> float:
    if barrier_kcal >= 99.0:
        return 0.0

    pre_exponential = get_reference_pre_exponential(family)
    exponent = -float(barrier_kcal) / (ARRHENIUS_R_KCAL * float(temperature_kelvin))
    return pre_exponential * math.exp(exponent) * max(float(multiplier), 0.0)


def effective_barrier_from_rate_constant(
    rate_constant: float,
    temperature_kelvin: float,
    family: Optional[str] = None,
    *,
    fallback_barrier: float = 99.0,
) -> float:
    if rate_constant <= 0.0:
        return fallback_barrier

    pre_exponential = get_reference_pre_exponential(family)
    if pre_exponential <= 0.0:
        return fallback_barrier

    ratio = rate_constant / pre_exponential
    if ratio <= 0.0:
        return fallback_barrier

    return max(0.0, -ARRHENIUS_R_KCAL * float(temperature_kelvin) * math.log(ratio))


