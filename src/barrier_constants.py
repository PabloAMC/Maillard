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
* data/Gemini_Deep_Research/maillard_meat.md, data/Gemini_Deep_Research/maillard_plant_based.md (project literature reviews)
"""

import json
import yaml
import math
from pathlib import Path
from typing import Any, Dict, Tuple, Optional, Sequence

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
ARRHENIUS_FILE = ROOT / "data" / "lit" / "arrhenius_params.yml"
REFINEMENT_PATCH_FILE = ROOT / "data" / "lit" / "refinement_surrogate_patches.json"
COMPUTATIONAL_PRIORS_FILE = ROOT / "data" / "lit" / "computational_priors.json"

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
DEFAULT_REFERENCE_PREEXPONENTIAL: float = 1e11
ARRHENIUS_R_KCAL: float = 0.001987

# Heme catalyst barrier reduction (kcal/mol)
HEME_CATALYST_REDUCTION: float = 5.0
HEME_CATALYST_FAMILIES = frozenset({"Strecker_Degradation", "Aminoketone_Condensation", "Lipid_Strecker_Synergy"})

DONOR_REACTIVITY_MULTIPLIERS: Dict[str, Dict[str, float]] = {
    "schiff_condensation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.94,
    },
    "amadori_rearrangement": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.94,
    },
    "heyns_rearrangement": {
        "phosphorylated": 1.04,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.98,
    },
    "enolisation_intermediate": {
        "phosphorylated": 1.12,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.92,
    },
    "1,2-enolisation": {
        "phosphorylated": 1.12,
        "pentose": 1.00,
        "fructose": 1.00,
        "glucose": 0.92,
    },
    "2,3-enolisation": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.95,
    },
    "dehydration": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.94,
    },
    "strecker_degradation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.88,
    },
    "aminoketone_condensation": {
        "phosphorylated": 1.10,
        "pentose": 1.00,
        "fructose": 1.02,
        "glucose": 0.88,
    },
    "thiol_addition": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.90,
    },
    "thiohemiacetal_formation": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.92,
    },
    "thiol_dehydration": {
        "phosphorylated": 1.08,
        "pentose": 1.00,
        "fructose": 0.98,
        "glucose": 0.92,
    },
}

_DONOR_PRIORITY = {
    "phosphorylated": 4,
    "pentose": 3,
    "fructose": 2,
    "glucose": 1,
}

# Reaction keys for which a literature- or family-anchored barrier is
# carried as part of the computational-gap kinetic priors.  Originally
# these were intended as DFT anchors; after retirement of unreliable
# xTB/DFT-derived values (2026-04-21) the keys are kept for runtime
# wiring but the underlying tiers are literature/family/surrogate.
# `asparagine_sugar_explicit_water_cluster` is excluded here because it
# is wired through the safety pipeline (`src/safety.py`) using the
# Knol 2009 lumped Ea, not through the DFT anchor metadata.
_DFT_ANCHOR_REACTION_KEYS = (
    "hexanal_radical_quench",
    "quinone_cys_michael",
    "aa_ring_open_dicarbonyl",
    "pe_schiff_base",
    "pe_amadori",
    "lysinoalanine_crosslink",
)


def _fallback_dft_anchor_metadata() -> Dict[str, Dict[str, Any]]:
    return {
        "hexanal_radical_quench": {
            "slr_family": "11",
            "current_tier": "no_literature_anchor",
            "target_tier": "wet_lab_anchor",
            "active_key": "hexanal_radical_quench_no_anchor",
            "active_barrier_kcal": None,
            "dft_key": "hexanal_radical_quench_dft",
            "uncertainty_kj": None,
            "promotion_ceiling": "wet_lab_required",
            "honest_label": "No reliable kinetic anchor; xTB-derived value retired 2026-04-21. Runtime falls back to literature_runtime barrier_kj_mol=31.72 default.",
        },
        "quinone_cys_michael": {
            "slr_family": "13",
            "current_tier": "literature_family_surrogate",
            "target_tier": "literature_derived_transfer",
            "active_key": "quinone_cys_michael_thiol_addition_family",
            "active_barrier_kcal": 6.93,
            "dft_key": "quinone_cys_michael_dft",
            "uncertainty_kj": 15.0,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "Family surrogate (Michael thiol addition, ~29 kJ/mol literature); xTB value retired 2026-04-21",
        },
        "aa_ring_open_dicarbonyl": {
            "slr_family": "14",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "aa_ring_open_dicarbonyl_hcw_family14_surrogate",
            "active_barrier_kcal": 7.58,
            "dft_key": "aa_ring_open_dicarbonyl_dft",
            "uncertainty_kj": 20.0,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "HCW Family 14 surrogate for an upstream ascorbic-acid dicarbonyl source, DFT validation pending, ±factor 2-5; bounded calibration",
        },
        "pe_schiff_base": {
            "slr_family": "15",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "pe_schiff_base_lit_derived",
            "active_barrier_kcal": 22.21,
            "dft_key": "pe_schiff_base_dft",
            "uncertainty_kj": 20.9,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "SLR 15 literature Ea=92.9 kJ/mol, ±factor 2-5; bounded calibration",
        },
        "pe_amadori": {
            "slr_family": "15",
            "current_tier": "literature_derived_transfer",
            "target_tier": "selective_dft_anchor",
            "active_key": "pe_amadori_lit_derived",
            "active_barrier_kcal": 19.81,
            "dft_key": "pe_amadori_dft",
            "uncertainty_kj": 20.9,
            "promotion_ceiling": "bounded_calibration",
            "honest_label": "SLR 15 literature Ea=82.9 kJ/mol, ±factor 2-5; bounded calibration",
        },
        "lysinoalanine_crosslink": {
            "slr_family": "12",
            "current_tier": "family_rule_surrogate",
            "target_tier": "selective_dft_anchor",
            "active_key": "lysinoalanine_crosslink_dha_family_surrogate",
            "active_barrier_kcal": 16.0,
            "dft_key": "lysinoalanine_crosslink_dft",
            "uncertainty_kj": 42.0,
            "promotion_ceiling": "ranking_only",
            "honest_label": "DHA-plus-lysine family surrogate, assumes prior dehydroalanine formation, ±factor 5-15; ranking only",
        },
        "asparagine_sugar_explicit_water_cluster": {
            "slr_family": "13_safety",
            "current_tier": "literature_derived_transfer",
            "target_tier": "literature_derived_transfer",
            "active_key": "asparagine_sugar_explicit_water_cluster_knol2009",
            "active_barrier_kcal": 30.83,
            "dft_key": "asparagine_sugar_explicit_water_cluster_dft",
            "uncertainty_kj": 10.0,
            "promotion_ceiling": "literature_anchor",
            "honest_label": "Knol et al. 2009 Ea=129 kJ/mol (DOI 10.1016/j.foodchem.2009.11.049); direct lumped Asn+sugar->acrylamide literature anchor",
        },
    }


def _load_dft_anchor_metadata() -> Dict[str, Dict[str, Any]]:
    fallback = _fallback_dft_anchor_metadata()
    try:
        with open(COMPUTATIONAL_PRIORS_FILE, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return fallback

    entries = {
        str(row.get("reaction_key", "")).strip(): row
        for row in payload.get("dft_kinetic_priors", {}).get("entries", [])
        if str(row.get("reaction_key", "")).strip()
    }
    metadata: Dict[str, Dict[str, Any]] = {}
    for reaction_key in _DFT_ANCHOR_REACTION_KEYS:
        row = entries.get(reaction_key)
        if row is None:
            metadata[reaction_key] = dict(fallback[reaction_key])
            continue
        barrier_kcal = row.get("barrier_kcal_mol")
        active_barrier_kcal = (
            float(barrier_kcal)
            if barrier_kcal is not None
            else fallback[reaction_key]["active_barrier_kcal"]
        )
        uncertainty_raw = row.get("uncertainty_kj", fallback[reaction_key]["uncertainty_kj"])
        uncertainty_kj = float(uncertainty_raw) if uncertainty_raw is not None else None
        metadata[reaction_key] = {
            "slr_family": str(row.get("slr_family", fallback[reaction_key]["slr_family"])),
            "current_tier": str(row.get("current_tier", fallback[reaction_key]["current_tier"])),
            "target_tier": str(row.get("target_tier", fallback[reaction_key]["target_tier"])),
            "active_key": str(row.get("active_arrhenius_key", fallback[reaction_key]["active_key"])),
            "active_barrier_kcal": active_barrier_kcal,
            "dft_key": f"{reaction_key}_dft",
            "uncertainty_kj": uncertainty_kj,
            "promotion_ceiling": str(row.get("promotion_ceiling", fallback[reaction_key]["promotion_ceiling"])),
            "honest_label": str(row.get("honest_label", fallback[reaction_key]["honest_label"])),
        }
    return metadata


DFT_ANCHOR_METADATA = _load_dft_anchor_metadata()


def _normalize_donor_context_token(value: str) -> str:
    return " ".join(str(value).strip().lower().replace("_", " ").replace("-", " ").split())


def infer_carbohydrate_donor_identity(reactant_labels: Optional[Sequence[str]] = None) -> Optional[str]:
    if not reactant_labels:
        return None

    detected: Dict[str, int] = {}
    for raw_label in reactant_labels:
        token = _normalize_donor_context_token(str(raw_label))
        if not token:
            continue
        if any(item in token for item in ["ribose 5 phosphate", "glucose 6 phosphate", "fructose 6 phosphate", "r5p"]):
            detected["phosphorylated"] = max(detected.get("phosphorylated", 0), _DONOR_PRIORITY["phosphorylated"])
            continue
        if any(item in token for item in ["ribose", "xylose", "arabinose"]):
            detected["pentose"] = max(detected.get("pentose", 0), _DONOR_PRIORITY["pentose"])
            continue
        if "fructose" in token:
            detected["fructose"] = max(detected.get("fructose", 0), _DONOR_PRIORITY["fructose"])
            continue
        if "glucose" in token:
            detected["glucose"] = max(detected.get("glucose", 0), _DONOR_PRIORITY["glucose"])

    if not detected:
        return None
    return max(detected.items(), key=lambda item: item[1])[0]


def get_donor_reactivity_multiplier(
    reaction_family: Optional[str],
    *,
    reactant_labels: Optional[Sequence[str]] = None,
    donor_identity: Optional[str] = None,
) -> float:
    donor = str(donor_identity or infer_carbohydrate_donor_identity(reactant_labels) or "").strip().lower()
    if not donor:
        return 1.0

    family_key = _canonical_fast_family(reaction_family)
    if not family_key:
        return 1.0

    family_multipliers = DONOR_REACTIVITY_MULTIPLIERS.get(family_key, {})
    return float(family_multipliers.get(donor, 1.0))


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
    }
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
            "estimated": 3.5
        }
        uncertainty = quality_unc_map.get(quality, 3.5)
            
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


