from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping


def _lane_toggle_magnitude(adjustment: Mapping[str, Any]) -> float:
    return (
        abs(float(adjustment.get("target_score_delta", 0.0) or 0.0))
        + abs(float(adjustment.get("maillard_closure_delta", 0.0) or 0.0))
        + abs(float(adjustment.get("off_flavour_risk_delta", 0.0) or 0.0))
    )


def build_family_lane_sensitivity_payload(flavor_axis_summary: Mapping[str, Any]) -> Dict[str, Any]:
    family_lane_summary = flavor_axis_summary.get("family_lane_summary", {}) or {}
    lane_adjustments = (flavor_axis_summary.get("family_lane_adjustments", {}) or {}).get("per_lane", {}) or {}
    rows: List[Dict[str, Any]] = []
    for slr_family, lane in family_lane_summary.items():
        adjustment = lane_adjustments.get(slr_family, {}) or {}
        rows.append(
            {
                "slr_family": str(slr_family),
                "display_name": str(lane.get("display_name", slr_family)),
                "active": bool(lane.get("active", False)),
                "strategic_posture": str(lane.get("strategic_posture", "unknown")),
                "target_score_delta": float(adjustment.get("target_score_delta", 0.0) or 0.0),
                "maillard_closure_delta": float(adjustment.get("maillard_closure_delta", 0.0) or 0.0),
                "off_flavour_risk_delta": float(adjustment.get("off_flavour_risk_delta", 0.0) or 0.0),
                "toggle_magnitude": float(_lane_toggle_magnitude(adjustment)),
            }
        )
    rows.sort(key=lambda row: (-float(row.get("toggle_magnitude", 0.0)), row.get("slr_family", "99")))
    return {
        "summary": {
            "family_lane_count": len(rows),
            "active_family_lane_count": sum(1 for row in rows if row.get("active")),
            "sensitivity_policy": "family_lane_sensitivity_tracks_runtime_toggle_impact_not_barrier_offsets",
        },
        "family_lanes": rows,
    }


def build_multi_run_family_lane_sensitivity_payload(rows: Iterable[Mapping[str, Any]]) -> Dict[str, Any]:
    payloads = [build_family_lane_sensitivity_payload(row) for row in rows]
    aggregate: Dict[str, Dict[str, Any]] = {}
    for payload in payloads:
        for row in payload.get("family_lanes", []):
            bucket = aggregate.setdefault(
                str(row.get("slr_family", "99")),
                {
                    "slr_family": str(row.get("slr_family", "99")),
                    "display_name": str(row.get("display_name", "unknown")),
                    "run_count": 0,
                    "active_run_count": 0,
                    "mean_toggle_magnitude": 0.0,
                },
            )
            bucket["run_count"] += 1
            bucket["active_run_count"] += int(bool(row.get("active", False)))
            bucket["mean_toggle_magnitude"] += float(row.get("toggle_magnitude", 0.0))
    result_rows: List[Dict[str, Any]] = []
    for row in aggregate.values():
        row["mean_toggle_magnitude"] = float(row["mean_toggle_magnitude"]) / max(int(row["run_count"]), 1)
        result_rows.append(row)
    result_rows.sort(key=lambda row: (-float(row.get("mean_toggle_magnitude", 0.0)), row.get("slr_family", "99")))
    return {
        "summary": {
            "run_count": len(payloads),
            "family_lane_count": len(result_rows),
            "sensitivity_policy": "family_lane_sensitivity_tracks_runtime_toggle_impact_not_barrier_offsets",
        },
        "family_lanes": result_rows,
    }


def render_family_lane_sensitivity_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Lane Sensitivity",
        "",
        "| SLR | Family Lane | Active | Target Δ | Closure Δ | Off-flavour Δ | Toggle Magnitude |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in payload.get("family_lanes", []):
        lines.append(
            f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {bool(row.get('active', False))} | {float(row.get('target_score_delta', 0.0)):+.2f} | {float(row.get('maillard_closure_delta', 0.0)):+.2f} | {float(row.get('off_flavour_risk_delta', 0.0)):+.2f} | {float(row.get('toggle_magnitude', row.get('mean_toggle_magnitude', 0.0))):.2f} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Family lanes summarized: {int(summary.get('family_lane_count', 0))}",
        f"Sensitivity policy: {summary.get('sensitivity_policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"