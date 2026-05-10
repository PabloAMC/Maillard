from __future__ import annotations

from typing import Any, Dict, Mapping

from src.benchmark_validation import assess_matrix_benchmark_evidence, load_benchmark
from src.primary_benchmark_campaign import build_matrix_primary_benchmark_campaign


EXTERNAL_BENCHMARKS = {
    "pea_iso": "pea_isolate_40C_PratapSingh2021",
    "soy_iso": "soy_isolate_40C_PratapSingh2021",
}


def _repo_benchmark_path(benchmark_id: str) -> str:
    return f"data/benchmarks/{benchmark_id}.json"


def build_pea_soy_external_evidence_artifact() -> Dict[str, Any]:
    campaign = build_matrix_primary_benchmark_campaign()
    campaign_by_matrix = {row["matrix"]: row for row in campaign.get("arms", [])}

    lanes = []
    for protein_type, benchmark_id in EXTERNAL_BENCHMARKS.items():
        bench = load_benchmark(_repo_benchmark_path(benchmark_id))
        evidence = assess_matrix_benchmark_evidence(bench)
        campaign_row = campaign_by_matrix.get(protein_type, {})
        lanes.append(
            {
                "protein_type": protein_type,
                "external_benchmark_id": benchmark_id,
                "external_target_profile": evidence.target_profile,
                "external_data_status": evidence.external_data_status,
                "external_reference_signal_origin": evidence.reference_signal_origin,
                "promotion_ready_today": bool(evidence.promotable),
                "current_external_anchor_compounds": sorted((bench.get("measured_volatiles") or {}).keys()),
                "mixed_meaty_positive_external_present": evidence.target_profile in {"mixed", "meaty_positive"},
                "required_external_package": {
                    "benchmark_id": campaign_row.get("benchmark_id"),
                    "temperature_c": campaign_row.get("protocol_temperature_c"),
                    "ph": campaign_row.get("protocol_ph"),
                    "time_points_min": list(campaign_row.get("protocol_time_points_min", [])),
                    "desirable_targets": list(campaign_row.get("required_desirable_targets", [])),
                    "adverse_targets": list(campaign_row.get("required_adverse_targets", [])),
                },
                "package_effect_if_landed_today": {
                    "would_close_requirements": list(campaign_row.get("would_close_requirements", [])),
                    "remaining_after_package": list(campaign_row.get("remaining_after_protocol", [])),
                },
                "promotion_blocker": str(campaign_row.get("promotion_blocker", "unknown")),
                "next_best_action": "land_external_quantitative_mixed_matrix_benchmark_package",
            }
        )

    return {
        "summary": {
            "lane_count": len(lanes),
            "external_mixed_meaty_positive_ready": sum(1 for row in lanes if row.get("mixed_meaty_positive_external_present")),
            "promotion_ready_today": sum(1 for row in lanes if row.get("promotion_ready_today")),
            "policy": "ambient_external_matrix_only_benchmarks_strengthen_external_anchor_coverage_but_do_not_replace_the_missing_mixed_meaty_positive_external_package",
        },
        "lanes": lanes,
    }


def render_pea_soy_external_evidence_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Pea Soy External Evidence",
        "",
        "| Protein | External Benchmark | Target Profile | External Data Status | Promotion Ready Today | Mixed Meaty-Positive External Present | Next Best Action |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("lanes", []):
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {row.get('external_benchmark_id', 'unknown')} | {row.get('external_target_profile', 'unknown')} | {row.get('external_data_status', 'unknown')} | {'yes' if row.get('promotion_ready_today', False) else 'no'} | {'yes' if row.get('mixed_meaty_positive_external_present', False) else 'no'} | {row.get('next_best_action', 'unknown')} |"
        )

    lines.extend(["", "## Required External Package", "", "| Protein | Benchmark Candidate | Temp C | pH | Time Points | Desirable Targets | Adverse Targets |", "| --- | --- | ---: | ---: | --- | --- | --- |"])
    for row in payload.get("lanes", []):
        package = row.get("required_external_package", {})
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {package.get('benchmark_id', 'unknown')} | {'' if package.get('temperature_c') is None else f'{float(package.get('temperature_c')):.1f}'} | {'' if package.get('ph') is None else f'{float(package.get('ph')):.1f}'} | {', '.join(str(item) for item in package.get('time_points_min', [])) or 'none'} | {', '.join(str(item) for item in package.get('desirable_targets', [])) or 'none'} | {', '.join(str(item) for item in package.get('adverse_targets', [])) or 'none'} |"
        )

    lines.extend(["", "## Promotion Delta If Package Lands", "", "| Protein | Would Close | Still Remaining | Current Promotion Blocker |", "| --- | --- | --- | --- |"])
    for row in payload.get("lanes", []):
        effect = row.get("package_effect_if_landed_today", {})
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {', '.join(str(item) for item in effect.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in effect.get('remaining_after_package', [])) or 'none'} | {row.get('promotion_blocker', 'unknown')} |"
        )

    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Lanes tracked: {int(summary.get('lane_count', 0))}",
        f"External mixed/meaty-positive ready lanes: {int(summary.get('external_mixed_meaty_positive_ready', 0))}",
        f"Promotion-ready lanes today: {int(summary.get('promotion_ready_today', 0))}",
        f"Policy: {summary.get('policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"