from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from src.xtb_path_quality import assess_xtb_path_quality


ROOT = Path(__file__).resolve().parents[1]
MANAGED_READY_STATUSES = {"completed", "completed_cached"}

PROXY_LANES: tuple[Dict[str, str], ...] = (
    {
        "target_id": "pe_schiff_base",
        "family": "15",
        "display_name": "PE Schiff base",
        "geometry_dir": "data/geometries/xtb_inputs/pe_schiff_base",
        "execution_artifact": "results/computational_gap_refinement/pe_schiff_base_xtb_execution.json",
        "promotion_rule": "Family 15 pair gate: both PE steps need xTB path plus clean quality gate.",
    },
    {
        "target_id": "pe_amadori",
        "family": "15",
        "display_name": "PE Amadori",
        "geometry_dir": "data/geometries/xtb_inputs/pe_amadori",
        "execution_artifact": "results/computational_gap_refinement/pe_amadori_xtb_execution.json",
        "promotion_rule": "Family 15 pair gate: both PE steps need xTB path plus clean quality gate.",
    },
)


def build_computational_gap_proxy_readiness_artifact() -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    for lane in PROXY_LANES:
        runner_dir = ROOT / lane["geometry_dir"]
        quality = assess_xtb_path_quality(runner_dir)
        execution_artifact = ROOT / lane["execution_artifact"]
        execution_row: Dict[str, Any] = {}
        if execution_artifact.exists():
            execution_payload = json.loads(execution_artifact.read_text(encoding="utf-8"))
            execution_row = dict((execution_payload.get("jobs") or [{}])[0])
        rows.append(
            {
                **lane,
                "frame_count": int(quality["frame_count"]),
                "has_path_bundle": bool(quality["has_path_bundle"]),
                "has_ts_guess": bool(quality["has_ts_guess"]),
                "quality_gate_passed": bool(quality["quality_gate_passed"]),
                "failure_markers": list(quality["failure_markers"]),
                "synthesized_outputs": [path.relative_to(ROOT).as_posix() for path in quality["synthesized_outputs"]],
                "managed_execution_present": execution_artifact.exists(),
                "managed_execution_status": str(execution_row.get("status", "not_run")),
                "managed_quality_gate_passed": execution_row.get("quality_gate_passed", bool(quality["quality_gate_passed"])),
            }
        )

    family_15_pair_ready = all(
        row["managed_execution_present"]
        and row["managed_execution_status"] in MANAGED_READY_STATUSES
        and bool(row.get("managed_quality_gate_passed"))
        and row["has_path_bundle"]
        and row["has_ts_guess"]
        and row["quality_gate_passed"]
        for row in rows
    )

    for row in rows:
        if not row["managed_execution_present"]:
            row["promotion_status"] = "proxy_only_unmanaged_execution"
            row["promotion_blocker"] = "managed_xtb_execution_missing"
        elif row["managed_execution_status"] not in MANAGED_READY_STATUSES:
            row["promotion_status"] = "proxy_only_execution_status"
            row["promotion_blocker"] = row["managed_execution_status"]
        elif not row["has_ts_guess"] or not row["has_path_bundle"]:
            row["promotion_status"] = "proxy_only_missing_path"
            row["promotion_blocker"] = "missing_xtb_path_or_ts"
        elif not row["quality_gate_passed"]:
            row["promotion_status"] = "proxy_only_quality_warning"
            row["promotion_blocker"] = "xtb_quality_warnings"
        elif not family_15_pair_ready:
            row["promotion_status"] = "proxy_only_pair_gate"
            row["promotion_blocker"] = "family_15_pair_gate"
        else:
            row["promotion_status"] = "candidate_for_formal_queue"
            row["promotion_blocker"] = ""

    return {
        "summary": {
            "proxy_lane_count": len(rows),
            "family_15_pair_ready": family_15_pair_ready,
            "formal_candidate_count": sum(1 for row in rows if row["promotion_status"] == "candidate_for_formal_queue"),
            "quality_warning_count": sum(1 for row in rows if row["promotion_blocker"] == "xtb_quality_warnings"),
        },
        "lanes": rows,
    }


def render_computational_gap_proxy_readiness_markdown(payload: Dict[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Computational Gap Proxy Readiness",
        "",
        "This artifact distinguishes proxy-useful xTB lanes from lanes that are clean enough to be promoted into the formal DFT queue.",
        "",
        f"- Proxy lanes tracked: {int(summary.get('proxy_lane_count', 0))}",
        f"- Family 15 pair gate ready: {bool(summary.get('family_15_pair_ready', False))}",
        f"- Formal queue candidates: {int(summary.get('formal_candidate_count', 0))}",
        f"- Lanes with xTB quality warnings: {int(summary.get('quality_warning_count', 0))}",
        "",
        "| Lane | Family | Frames | Path | TS | Quality Gate | Promotion Status | Blocker |",
        "| --- | --- | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("lanes", []):
        lines.append(
            f"| {row.get('display_name', row.get('target_id', 'unknown'))} | {row.get('family', 'unknown')} | {int(row.get('frame_count', 0))} | {bool(row.get('has_path_bundle', False))} | {bool(row.get('has_ts_guess', False))} | {bool(row.get('quality_gate_passed', False))} | {row.get('promotion_status', '')} | {row.get('promotion_blocker', '') or '-'} |"
        )
    return "\n".join(lines) + "\n"