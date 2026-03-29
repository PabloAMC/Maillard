from __future__ import annotations

from typing import Any, Dict, List, Mapping

from src.barrier_constants import DFT_ANCHOR_METADATA, FAST_BARRIERS
from src.family_barrier_progress import build_family_barrier_progress_artifact
from src.family_ingestion_plan import load_family_ingestion_plan
from src.mlp_adoption_contract import load_mlp_candidates
from src.reaction_benchmark import load_reaction_benchmark_entries


def _plan_lookup() -> Dict[str, Mapping[str, Any]]:
    payload = load_family_ingestion_plan()
    return {
        str(row.get("family_id", "")).strip(): row
        for row in payload.get("families", [])
        if str(row.get("family_id", "")).strip()
    }


def _reaction_lookup() -> Dict[str, Any]:
    return {
        str(entry.reaction_family).strip(): entry
        for entry in load_reaction_benchmark_entries()
    }


def _literature_status(row: Mapping[str, Any]) -> str:
    explicit_fast = int(row.get("explicit_fast_barrier_count", 0))
    arrhenius = int(row.get("arrhenius_anchor_count", 0))
    if explicit_fast == 0:
        return "not_a_primary_barrier_family"
    if arrhenius >= explicit_fast:
        return "literature_arrhenius_available"
    if arrhenius > 0:
        return "partially_literature_anchored"
    return "encoded_without_arrhenius_anchor"


def _benchmark_visibility(benchmark_visible_lanes: List[str], runtime_concept: str) -> str:
    if benchmark_visible_lanes:
        return "benchmark_visible_reaction_family"
    if runtime_concept in {
        "matrix_family_scope_extension",
        "upstream_pretreatment_node",
        "thiamine_fragmentation_support",
        "sulfur_peptide_support",
        "nucleotide_and_ribose_support",
    }:
        return "support_or_state_lane_not_yet_reaction_benchmarked"
    return "not_yet_benchmark_visible"


def _xtb_status(literature_status: str, benchmark_visibility: str) -> str:
    if benchmark_visibility == "support_or_state_lane_not_yet_reaction_benchmarked":
        return "not_primary_use_case_for_xTB"
    if benchmark_visibility == "not_yet_benchmark_visible":
        return "defer_until_benchmark_visible"
    if literature_status == "literature_arrhenius_available":
        return "reserve_for_reopened_benchmark_residuals"
    return "screen_missing_or_unanchored_barriers_first"


def _dft_status(literature_status: str, benchmark_visibility: str, quantitative_points: int) -> str:
    if benchmark_visibility == "support_or_state_lane_not_yet_reaction_benchmarked":
        return "not_primary_use_case_for_family_wide_dft"
    if benchmark_visibility == "not_yet_benchmark_visible":
        return "defer_until_benchmark_visible"
    if literature_status == "literature_arrhenius_available" and quantitative_points > 0:
        return "selective_only_if_residuals_reopen"
    return "selective_dft_for_benchmark_visible_missing_or_unanchored_steps"


def _candidate_map_by_role_and_motif() -> Dict[str, Dict[str, List[str]]]:
    mapping: Dict[str, Dict[str, List[str]]] = {}
    for candidate in load_mlp_candidates():
        role_map = mapping.setdefault(str(candidate.proposed_role), {})
        for motif in candidate.target_motif_families:
            role_map.setdefault(str(motif), []).append(str(candidate.candidate_id))
    return mapping


def _mlp_roles(motif_families: List[str]) -> Dict[str, List[str]]:
    role_map = _candidate_map_by_role_and_motif()
    payload: Dict[str, List[str]] = {}
    for role, motif_map in role_map.items():
        matched = sorted(
            {
                candidate_id
                for motif in motif_families
                for candidate_id in motif_map.get(motif, [])
            }
        )
        if matched:
            payload[role] = matched
    return payload


def _next_compute_action(*, literature_status: str, benchmark_visibility: str, family_id: str) -> str:
    if benchmark_visibility == "support_or_state_lane_not_yet_reaction_benchmarked":
        return "stay_literature_and_calibration_first_no_family_wide_qm"
    if benchmark_visibility == "not_yet_benchmark_visible":
        return "defer_compute_escalation_until_the_family_is_benchmark_visible"
    if family_id == "lipid_oxidation_and_carbonylic_crosstalk":
        return "close_literature_and_observable_gaps_then_xTB_screen_lipid_crosstalk_steps_before_selective_DFT"
    if literature_status == "literature_arrhenius_available":
        return "keep_the_family_literature_first_and_use_xTB_or_DFT_only_for_reopened_benchmark_residuals"
    return "fill_missing_literature_anchors_then_xTB_screen_then_run_selective_DFT_only_for_decisive_steps"


def build_c4_c5_dft_status() -> List[Dict[str, Any]]:
    """Return per-reaction DFT status for the 6 C4/C5 targets.

    Reads from ``DFT_ANCHOR_METADATA`` (barrier_constants.py) and resolves
    whether the *_dft sentinel has been filled by ingest_dft_c4_c5_results.py.
    This is the canonical machine-readable surface for reaction-level DFT tier.
    """
    rows: List[Dict[str, Any]] = []
    for reaction_key, meta in DFT_ANCHOR_METADATA.items():
        dft_key = str(meta.get("dft_key", ""))
        dft_filled = False
        dft_barrier_kcal: Any = None
        if dft_key and dft_key in FAST_BARRIERS:
            val = FAST_BARRIERS[dft_key][0]
            if val < 99.0:
                dft_filled = True
                dft_barrier_kcal = val

        # Resolve the active barrier from the routing table
        active_key = str(meta.get("active_key", ""))
        active_barrier_kcal: Any = None
        if active_key and active_key in FAST_BARRIERS:
            active_barrier_kcal = FAST_BARRIERS[active_key][0]

        rows.append({
            "reaction_key":       reaction_key,
            "slr_family":         str(meta.get("slr_family", "??")).zfill(2),
            "current_tier":       str(meta.get("current_tier", "unknown")),
            "target_tier":        str(meta.get("target_tier", "selective_dft_anchor")),
            "active_key":         active_key,
            "active_barrier_kcal": active_barrier_kcal,
            "dft_key":            dft_key,
            "dft_filled":         dft_filled,
            "dft_barrier_kcal":   dft_barrier_kcal,
            "uncertainty_kj":     float(meta.get("uncertainty_kj", 42.0)),
            "promotion_ceiling":  str(meta.get("promotion_ceiling", "ranking_only")),
            "honest_label":       str(meta.get("honest_label", "")),
        })

    rows.sort(key=lambda r: (r["slr_family"], r["reaction_key"]))
    return rows


def build_dft_coverage_map_artifact() -> Dict[str, Any]:
    barrier_progress = build_family_barrier_progress_artifact()
    barrier_lookup = {
        str(row.get("chemistry_family", "")).strip(): row
        for row in barrier_progress.get("families", [])
        if str(row.get("chemistry_family", "")).strip()
    }
    plan_lookup = _plan_lookup()
    reaction_lookup = _reaction_lookup()

    rows: List[Dict[str, Any]] = []
    for family_id, plan_row in plan_lookup.items():
        barrier_row = barrier_lookup.get(family_id, {})
        mapped_lanes = list(barrier_row.get("mapped_barrier_lanes", []))
        benchmark_entries = [reaction_lookup[lane] for lane in mapped_lanes if lane in reaction_lookup]
        benchmark_visible_lanes = [str(entry.reaction_family) for entry in benchmark_entries]
        motif_families = sorted({str(entry.motif_family) for entry in benchmark_entries})
        runtime_concept = str(plan_row.get("runtime_concept", "unknown"))
        literature_status = _literature_status(barrier_row)
        benchmark_visibility = _benchmark_visibility(benchmark_visible_lanes, runtime_concept)
        quantitative_points = int(barrier_row.get("quantitative_point_count", 0))
        rows.append(
            {
                "slr_family": str(plan_row.get("slr_family", "99")).zfill(2),
                "family_id": family_id,
                "display_name": str(plan_row.get("display_name", family_id)),
                "strategic_posture": str(plan_row.get("strategic_posture", "unknown")),
                "runtime_concept": runtime_concept,
                "literature_status": literature_status,
                "benchmark_visibility": benchmark_visibility,
                "benchmark_visible_lanes": benchmark_visible_lanes,
                "motif_families": motif_families,
                "explicit_fast_barrier_count": int(barrier_row.get("explicit_fast_barrier_count", 0)),
                "arrhenius_anchor_count": int(barrier_row.get("arrhenius_anchor_count", 0)),
                "quantitative_point_count": quantitative_points,
                "xtb_status": _xtb_status(literature_status, benchmark_visibility),
                "dft_status": _dft_status(literature_status, benchmark_visibility, quantitative_points),
                "mlp_roles": _mlp_roles(motif_families),
                "next_compute_action": _next_compute_action(
                    literature_status=literature_status,
                    benchmark_visibility=benchmark_visibility,
                    family_id=family_id,
                ),
            }
        )

    rows.sort(key=lambda row: (row["slr_family"], row["family_id"]))

    c4_c5_status = build_c4_c5_dft_status()
    dft_filled_count = sum(1 for r in c4_c5_status if r["dft_filled"])

    return {
        "summary": {
            "family_count": len(rows),
            "benchmark_visible_family_count": sum(
                1 for row in rows if row.get("benchmark_visibility") == "benchmark_visible_reaction_family"
            ),
            "literature_first_family_count": sum(
                1 for row in rows if row.get("literature_status") == "literature_arrhenius_available"
            ),
            "c4_c5_reaction_count": len(c4_c5_status),
            "c4_c5_dft_filled_count": dft_filled_count,
            "c4_c5_dft_pending_count": len(c4_c5_status) - dft_filled_count,
            "policy": "dft_and_mlp_coverage_must_follow_literature_first_then_xTB_then_selective_DFT_with_MLPs_limited_to_bounded_offline_roles",
        },
        "families": rows,
        "c4_c5_reactions": c4_c5_status,
    }


def render_dft_coverage_map_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# DFT Coverage Map",
        "",
        "| SLR | Family | Literature Status | Benchmark Visibility | xTB Status | DFT Status | MLP Roles | Next Compute Action |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        mlp_roles = row.get("mlp_roles", {})
        role_summary = "; ".join(
            f"{role}: {', '.join(candidates)}"
            for role, candidates in sorted(mlp_roles.items())
        ) or "none"
        lines.append(
            f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {row.get('literature_status', 'unknown')} | {row.get('benchmark_visibility', 'unknown')} | {row.get('xtb_status', 'unknown')} | {row.get('dft_status', 'unknown')} | {role_summary} | {row.get('next_compute_action', 'unknown')} |"
        )

    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Families tracked: {int(summary.get('family_count', 0))}",
        f"Benchmark-visible families: {int(summary.get('benchmark_visible_family_count', 0))}",
        f"Literature-first families: {int(summary.get('literature_first_family_count', 0))}",
        f"Policy: {summary.get('policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"