from __future__ import annotations

import json
from typing import Any, Dict, List, Mapping

from src.barrier_constants import FAST_BARRIERS, get_arrhenius_params
from src.family_deviation_audit import build_family_deviation_audit_artifact
from src.family_ingestion_plan import load_family_ingestion_plan
from src import data_paths


REFINEMENT_PATCH_FILE = data_paths.REFINEMENT_SURROGATE_PATCHES

FAMILY_BARRIER_LANES: Dict[str, List[str]] = {
    "amino_acid_sugar_core": [
        "schiff_condensation",
        "amadori_rearrangement",
        "enolisation_intermediate",
        "dehydration",
        "strecker_degradation",
        "aminoketone_condensation",
        "beta_elimination",
    ],
    "lipid_oxidation_and_carbonylic_crosstalk": [
        "lipid_homolysis",
        "beta_scission",
        "lipid_condensation",
        "lipid_strecker_synergy",
        "radical_crosstalk",
    ],
    "thiamine_fragmentation_support": [
        "thiamine_degradation",
        "thiol_addition",
        "thiol_oxidation",
    ],
    "nucleotide_and_ribose_support": [
        "mutarotation",
        "ring_opening",
        "dehydration",
    ],
    "glutathione_and_peptide_support": [
        "glutathione_cleavage",
        "cysteine_thermolysis",
        "thiol_addition",
    ],
    "alternative_protein_matrix_scope": [],
    "carbonyl_donor_hierarchy": [
        "mutarotation",
        "schiff_condensation",
        "amadori_rearrangement",
        "heyns_rearrangement",
        "enolisation_intermediate",
        "dehydration",
    ],
    "off_note_and_maillard_suppression": [
        "lipid_condensation",
        "lipid_homolysis",
        "beta_scission",
        "radical_crosstalk",
    ],
    "carbohydrate_pyrolysis_and_caramelization": [
        "mutarotation",
        "dehydration",
        "retro_aldol",
    ],
    "fermentation_pretreatment": [],
}


def _load_refinement_offsets() -> Dict[str, float]:
    if not REFINEMENT_PATCH_FILE.exists():
        return {}
    with open(REFINEMENT_PATCH_FILE, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    accepted = payload.get("accepted_offsets", {}) if isinstance(payload, dict) else {}
    if not isinstance(accepted, dict):
        return {}
    return {str(key): float(value) for key, value in accepted.items()}


def _barrier_lane_stage(*, explicit_count: int, arrhenius_count: int, quantitative_points: int) -> str:
    if explicit_count == 0:
        return "no_explicit_barrier_lane"
    if quantitative_points > 0:
        return "benchmark_exercised_barrier_lane"
    if arrhenius_count > 0:
        return "literature_anchored_barrier_lane"
    return "fast_heuristic_barrier_lane"


def _last_step_resolution(family_id: str) -> str:
    if family_id in {
        "lipid_oxidation_and_carbonylic_crosstalk",
        "off_note_and_maillard_suppression",
    }:
        return "hexanal_nonanal_internal_calibration_route_now_visible"
    return "no_new_barrier_constant_landed_in_last_step"


def _prediction_readout(row: Mapping[str, Any]) -> str:
    quantitative_points = int(row.get("quantitative_point_count", 0))
    max_ratio = row.get("max_ratio")
    if quantitative_points <= 0:
        return "No compound-level quantitative readout yet"
    if max_ratio is None:
        return "Quantitative points exist but no max-ratio summary is available"
    return f"{quantitative_points} quantitative points; current family max ratio {float(max_ratio):.3f}x"


def build_family_barrier_progress_artifact() -> Dict[str, Any]:
    plan = load_family_ingestion_plan()
    deviation_payload = build_family_deviation_audit_artifact()
    deviation_by_family = {
        str(row.get("chemistry_family", "unknown")): row
        for row in deviation_payload.get("family_deviation_table", [])
    }
    refinement_offsets = _load_refinement_offsets()

    rows = []
    for family in plan.get("families", []):
        family_id = str(family.get("family_id", "unknown"))
        mapped_lanes = list(FAMILY_BARRIER_LANES.get(family_id, []))
        explicit_fast_lanes = [lane for lane in mapped_lanes if lane in FAST_BARRIERS]
        arrhenius_lanes = [lane for lane in mapped_lanes if get_arrhenius_params(lane) is not None]
        patched_lanes = [lane for lane in mapped_lanes if lane in refinement_offsets]
        deviation_row = deviation_by_family.get(family_id, {})
        quantitative_points = int(deviation_row.get("quantitative_point_count", 0))

        row = {
            "slr_family": str(family.get("slr_family", "99")).zfill(2),
            "chemistry_family": family_id,
            "display_name": str(family.get("display_name", family_id)),
            "strategic_posture": str(family.get("strategic_posture", "unknown")),
            "runtime_concept": str(family.get("runtime_concept", "unknown")),
            "target_count": len(list(family.get("target_compounds_or_state_variables", []))),
            "mapped_barrier_lanes": mapped_lanes,
            "explicit_fast_barrier_count": len(explicit_fast_lanes),
            "arrhenius_anchor_count": len(arrhenius_lanes),
            "refinement_patch_count": len(patched_lanes),
            "benchmark_count": int(deviation_row.get("benchmark_count", 0)),
            "quantitative_point_count": quantitative_points,
            "max_ratio": deviation_row.get("max_ratio"),
            "mean_abs_log10_error": deviation_row.get("mean_abs_log10_error"),
            "barrier_lane_stage": _barrier_lane_stage(
                explicit_count=len(explicit_fast_lanes),
                arrhenius_count=len(arrhenius_lanes),
                quantitative_points=quantitative_points,
            ),
            "last_step_resolution": _last_step_resolution(family_id),
        }
        row["prediction_readout"] = _prediction_readout(row)
        rows.append(row)

    rows.sort(key=lambda item: (item.get("slr_family", "99"), item.get("chemistry_family", "unknown")))
    return {
        "summary": {
            "family_count": len(rows),
            "families_with_explicit_fast_barriers": sum(1 for row in rows if int(row.get("explicit_fast_barrier_count", 0)) > 0),
            "families_with_arrhenius_anchors": sum(1 for row in rows if int(row.get("arrhenius_anchor_count", 0)) > 0),
            "families_with_quantitative_prediction_readout": sum(1 for row in rows if int(row.get("quantitative_point_count", 0)) > 0),
            "families_without_explicit_barrier_lane": [
                row.get("chemistry_family", "unknown") for row in rows if int(row.get("explicit_fast_barrier_count", 0)) == 0
            ],
            "policy": "family_barrier_progress_must_distinguish_explicit_barrier_lanes_from_quantitative_prediction_readout_and_from_matrix_or_upstream_state_only_families",
        },
        "families": rows,
    }


def render_family_barrier_progress_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Barrier Progress",
        "",
        "This artifact tracks which chemistry families have explicit FAST barrier lanes, which are literature-anchored, and which already show benchmark-linked prediction readout.",
        "",
        "| SLR | Family | Barrier Stage | FAST Lanes | Arrhenius Anchors | Refinement Patches | Quantitative Points | Max Ratio | Last Step Resolution |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in payload.get("families", []):
        max_ratio = row.get("max_ratio")
        lines.append(
            f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {row.get('barrier_lane_stage', 'unknown')} | {int(row.get('explicit_fast_barrier_count', 0))} | {int(row.get('arrhenius_anchor_count', 0))} | {int(row.get('refinement_patch_count', 0))} | {int(row.get('quantitative_point_count', 0))} | {('-' if max_ratio is None else f'{float(max_ratio):.3f}')} | {row.get('last_step_resolution', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Prediction Readout",
        "",
        "| Family | Runtime Concept | Targets Tracked | Prediction Readout |",
        "| --- | --- | ---: | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('display_name', 'unknown')} | {row.get('runtime_concept', 'unknown')} | {int(row.get('target_count', 0))} | {row.get('prediction_readout', 'unknown')} |"
        )

    lines.append("")
    lines.append(f"Policy: {payload.get('summary', {}).get('policy', 'unknown')}")
    return "\n".join(lines) + "\n"