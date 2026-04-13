from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional
from pathlib import Path

from src.benchmark_validation import DEFAULT_TARGET_TAG, build_matrix_observable_closure_audit
from src.refinement_campaign import build_cheap_screening_artifact, build_selective_dft_plan
from src.refinement_watchlist import build_refinement_watchlist


def _watchlist_lookup(payload: Mapping[str, Any]) -> Dict[str, Mapping[str, Any]]:
    return {
        str(row.get("reaction_family", "")).strip().lower(): row
        for row in payload.get("candidates", [])
        if str(row.get("reaction_family", "")).strip()
    }


def build_offline_refinement_governance_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    closure_audit = build_matrix_observable_closure_audit(benchmark_files, target_tag=target_tag)
    cheap_payload = build_cheap_screening_artifact(benchmark_files, target_tag=target_tag)
    dft_payload = build_selective_dft_plan(benchmark_files, target_tag=target_tag)
    watchlist_payload = build_refinement_watchlist(benchmark_files, target_tag=target_tag)

    cheap_lookup = _watchlist_lookup(cheap_payload)
    dft_lookup = _watchlist_lookup(dft_payload)
    watchlist_lookup = _watchlist_lookup(watchlist_payload)
    benchmark_lookup = {
        str(item.get("benchmark_id", "")): item
        for item in closure_audit.get("benchmarks", [])
        if str(item.get("benchmark_id", "")).strip()
    }

    benchmark_rows: List[Dict[str, Any]] = []
    for row in closure_audit.get("mechanistic_refinement_watchlist", []):
        benchmark_id = str(row.get("benchmark_id", "unknown"))
        benchmark_meta = benchmark_lookup.get(benchmark_id, {})
        family_names = [str(family).strip() for family in row.get("candidate_reaction_families", []) if str(family).strip()]
        if not family_names:
            family_names = [
                str(candidate.get("reaction_family", "")).strip()
                for candidate in cheap_payload.get("candidates", [])
                if benchmark_id in set(candidate.get("linked_mechanistic_benchmarks", []))
                and str(candidate.get("reaction_family", "")).strip()
            ]
        candidate_families: List[Dict[str, Any]] = []
        seen: set[str] = set()
        for family_name in family_names:
            if not family_name:
                continue
            family_key = family_name.lower()
            if family_key in seen:
                continue
            seen.add(family_key)
            cheap_row = cheap_lookup.get(family_key, {})
            dft_row = dft_lookup.get(family_key, {})
            watchlist_row = watchlist_lookup.get(family_key, {})
            candidate_families.append(
                {
                    "reaction_family": family_name,
                    "watchlist_decision": str(watchlist_row.get("decision", "not_ranked")),
                    "cheap_screen_decision": str(cheap_row.get("screen_decision", "not_screened")),
                    "cheap_screen_improvement": float(cheap_row.get("total_improvement", 0.0) or 0.0),
                    "no_escalation_reason": str(cheap_row.get("no_escalation_reason", "not_screened")),
                    "dft_decision": str(dft_row.get("decision", "not_planned")),
                    "priority_tier": str(dft_row.get("priority_tier", "not_planned")),
                }
            )
        benchmark_rows.append(
            {
                "benchmark_id": benchmark_id,
                "protein_type": str(benchmark_meta.get("protein_type", row.get("protein_type", "unknown"))),
                "expected_decision_change": str(row.get("expected_decision_change", "unknown")),
                "target_compounds": list(row.get("target_compounds", [])),
                "observable_blockers": list(row.get("observable_blockers", []))
                or [
                    str(compound_row.get("compound", "unknown"))
                    for compound_row in benchmark_meta.get("compounds", [])
                    if str(compound_row.get("closure_action", "")) == "mechanistic_blocker"
                ],
                "candidate_reaction_families": candidate_families,
                "offline_compute_gate": "hold_observable_first",
            }
        )

    approved_jobs = list(dft_payload.get("offline_jobs", []))
    approved_job_families = sorted(
        {
            str(job.get("reaction_id", "")).split("family::", 1)[-1]
            for job in approved_jobs
            if str(job.get("reaction_id", "")).startswith("family::")
        }
    )
    advance_count = int(cheap_payload.get("summary", {}).get("advance", 0))
    run_now_count = int(dft_payload.get("summary", {}).get("run_now", 0))
    governing_status = "hold_observable_first"
    if run_now_count > 0:
        governing_status = "selective_refinement_allowed"

    blockers = [
        "observable closure still names the benchmark-visible compounds that must move before mechanistic escalation",
        "cheap-first screening produced no benchmark-visible improvement for the active mechanistic-priority families",
        "no selective DFT jobs are approved until a cheap-first perturbation improves benchmark-visible diagnostics",
    ]
    if run_now_count > 0:
        blockers[-1] = "approved selective DFT jobs remain sparse and benchmark-scoped"

    return {
        "summary": {
            "governing_status": governing_status,
            "mechanistic_priority_benchmark_count": len(benchmark_rows),
            "benchmarks_with_named_targets": sum(1 for row in benchmark_rows if row["target_compounds"]),
            "cheap_first_advance_count": advance_count,
            "approved_offline_job_count": run_now_count,
            "approved_offline_job_families": approved_job_families,
            "policy": "continue_only_when_observable_closure_names_the_targets_and_cheap_first_screening_improves_benchmark_visible_diagnostics",
        },
        "blockers": blockers,
        "mechanistic_priority_benchmarks": benchmark_rows,
        "approved_offline_jobs": approved_jobs,
    }


def render_offline_refinement_governance_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Offline Refinement Governance",
        "",
        f"Governing status: {summary.get('governing_status', 'unknown')}",
        f"Mechanistic-priority benchmarks: {int(summary.get('mechanistic_priority_benchmark_count', 0))}",
        f"Benchmarks with named target compounds: {int(summary.get('benchmarks_with_named_targets', 0))}",
        f"Cheap-first advances: {int(summary.get('cheap_first_advance_count', 0))}",
        f"Approved offline jobs: {int(summary.get('approved_offline_job_count', 0))}",
        f"Policy: {summary.get('policy', 'unknown')}",
        "",
        "## Benchmark Gates",
        "",
        "| Benchmark | Protein Type | Target Compounds | Expected Decision Change | Observable Blockers | Offline Gate |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("mechanistic_priority_benchmarks", []):
        lines.append(
            f"| {row.get('benchmark_id', 'unknown')} | {row.get('protein_type', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('target_compounds', [])) or 'none'} | "
            f"{row.get('expected_decision_change', 'unknown')} | "
            f"{'; '.join(str(item) for item in row.get('observable_blockers', [])) or 'none'} | "
            f"{row.get('offline_compute_gate', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Family Gate",
        "",
        "| Benchmark | Reaction Family | Watchlist | Cheap Screen | Improvement | DFT Decision | Reason |",
        "| --- | --- | --- | --- | ---: | --- | --- |",
    ])
    for row in payload.get("mechanistic_priority_benchmarks", []):
        benchmark_id = str(row.get("benchmark_id", "unknown"))
        families = row.get("candidate_reaction_families", [])
        if not families:
            lines.append(f"| {benchmark_id} | none | none | none | 0.00 | none | no candidate families were named |")
            continue
        for family in families:
            lines.append(
                f"| {benchmark_id} | {family.get('reaction_family', 'unknown')} | {family.get('watchlist_decision', 'unknown')} | "
                f"{family.get('cheap_screen_decision', 'unknown')} | {float(family.get('cheap_screen_improvement', 0.0)):.2f} | "
                f"{family.get('dft_decision', 'unknown')} | {family.get('no_escalation_reason', 'unknown')} |"
            )

    lines.extend([
        "",
        "## Blockers",
        "",
    ])
    for blocker in payload.get("blockers", []):
        lines.append(f"- {blocker}")
    if not payload.get("approved_offline_jobs"):
        lines.extend([
            "",
            "No selective DFT jobs are currently approved.",
        ])
    return "\n".join(lines) + "\n"