from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional

from src.benchmark_types import (
    BenchmarkIndexEntry,
    BenchmarkSummary,
    BenchmarkTargetSnapshot,
    MatrixBenchmarkAssertion,
    MatrixBenchmarkBranchDelta,
    MatrixBenchmarkDelta,
    MatrixBenchmarkEvidence,
    MatrixPromotionFamilyStatus,
    ThermodynamicGatingAudit,
)


def render_matrix_branch_deltas_markdown(
    rows: Iterable[MatrixBenchmarkBranchDelta],
    *,
    base_ref: str,
) -> str:
    delta_rows = list(rows)
    lines = [
        f"# Matrix Benchmark Branch Comparison vs {base_ref}",
        "",
        "| Benchmark | Compound | Change | Current Path | Base Path | Current Ref Origin | Base Ref Origin | Current Source Origin | Base Source Origin | Current Data Status | Base Data Status | Current Predicted ppb | Base Predicted ppb | Δ Predicted ppb | Δ Ratio |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in delta_rows:
        current_predicted = f"{row.current_predicted_ppb:.3f}" if row.current_predicted_ppb is not None else "n/a"
        baseline_predicted = f"{row.baseline_predicted_ppb:.3f}" if row.baseline_predicted_ppb is not None else "n/a"
        predicted_delta = f"{row.predicted_delta_ppb:.3f}" if row.predicted_delta_ppb is not None else "n/a"
        ratio_delta = f"{row.ratio_delta:.3f}" if row.ratio_delta is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.compound} | {row.change_type} | {row.current_execution_path} | {row.baseline_execution_path} | {row.current_reference_signal_origin} | {row.baseline_reference_signal_origin} | {row.current_source_origin} | {row.baseline_source_origin} | {row.current_external_data_status} | {row.baseline_external_data_status} | {current_predicted} | {baseline_predicted} | {predicted_delta} | {ratio_delta} |"
        )
    lines.extend([
        "",
        f"Changed rows: {len(delta_rows)}",
        f"Added rows: {sum(1 for row in delta_rows if row.change_type == 'added')}",
        f"Removed rows: {sum(1 for row in delta_rows if row.change_type == 'removed')}",
        f"Modified rows: {sum(1 for row in delta_rows if row.change_type == 'modified')}",
        f"Metadata-only changes: {sum(1 for row in delta_rows if row.change_type == 'metadata_changed')}",
    ])
    if delta_rows and not any(row.baseline_present for row in delta_rows):
        lines.extend([
            f"Base ref {base_ref} exposes no matrix benchmark rows under the current evaluator.",
            "All reported rows are additions from the current branch/worktree.",
        ])
    return "\n".join(lines) + "\n"


def render_family_lane_validation_markdown(payload: Dict[str, Any]) -> str:
    lines = [
        "# Family Lane Validation",
        "",
        "| SLR | Family | Posture | Benchmarks | Strict Ready | Supported | Payload Roles | Execution Paths |",
        "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in payload.get("families", []):
        execution_paths = (
            ", ".join(f"{key}={value}" for key, value in row.get("execution_paths", {}).items())
            or "none"
        )
        lines.append(
            f"| {row.get('slr_family', '') or 'n/a'} | {row.get('chemistry_family', 'unknown')}"
            f" | {row.get('strategic_posture', 'unknown')} | {int(row.get('benchmark_count', 0))}"
            f" | {int(row.get('strict_ready_count', 0))} | {int(row.get('supported_count', 0))}"
            f" | {', '.join(str(item) for item in row.get('payload_roles', [])) or 'none'}"
            f" | {execution_paths} |"
        )
    lines.extend([
        "",
        "## Lane Summary",
        "",
        "| Execution Path | Benchmarks | Strict Ready | Supported | Chemistry Families | Payload Roles |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ])
    for row in payload.get("lanes", []):
        lines.append(
            f"| {row.get('execution_path', 'unknown')} | {int(row.get('benchmark_count', 0))}"
            f" | {int(row.get('strict_ready_count', 0))} | {int(row.get('supported_count', 0))}"
            f" | {', '.join(str(item) for item in row.get('chemistry_families', [])) or 'none'}"
            f" | {', '.join(str(item) for item in row.get('payload_roles', [])) or 'none'} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Benchmarks summarized: {int(summary.get('benchmark_count', 0))}",
        f"Chemistry families summarized: {int(summary.get('family_count', 0))}",
        f"Execution lanes summarized: {int(summary.get('lane_count', 0))}",
    ])
    return "\n".join(lines) + "\n"


def render_matrix_benchmark_deltas_markdown(rows: Iterable[MatrixBenchmarkDelta]) -> str:
    deltas = list(rows)
    lines = [
        "# Matrix Benchmark Deltas",
        "",
        "| Benchmark | Protein | Path | Process State | Reference Origin | Ranking Contract | Compound | Role | Reference ppb | Predicted ppb | Abs Δ ppb | Δ % | Ratio | Calibration | Evidence | Fallback |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for row in deltas:
        pct = f"{100.0 * row.pct_delta:.1f}%" if row.pct_delta is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.protein_type} | {row.execution_path} | {row.process_state or 'n/a'} | {row.reference_signal_origin} | {row.ranking_contract_status} | {row.compound} | {row.role} | {row.reference_ppb:.3f} | {row.predicted_ppb:.3f} | {row.abs_delta_ppb:.3f} | {pct} | {row.ratio:.3f} | {row.calibration_source} | {row.calibration_evidence_strength} | {row.calibration_fallback_mode} |"
        )
    lines.extend([
        "",
        f"Delta rows: {len(deltas)}",
        f"Benchmarks covered: {len({row.benchmark_id for row in deltas})}",
    ])
    return "\n".join(lines) + "\n"


def render_thermodynamic_gating_audit_markdown(rows: Iterable[ThermodynamicGatingAudit]) -> str:
    audits = list(rows)
    lines = [
        "# Thermodynamic Gating Audit",
        "",
        "| Benchmark | Path | Applicable | Baseline Status | Gated Status | Δ MAE ppb | Δ Max Ratio | Material | Recommended Policy | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in audits:
        delta_mae = f"{row.delta_mae_ppb:.2f}" if row.delta_mae_ppb is not None else "n/a"
        delta_ratio = f"{row.delta_max_ratio:.3f}" if row.delta_max_ratio is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.execution_path} | {'yes' if row.applicable else 'no'} | {row.baseline_overall_status} | {row.gated_overall_status} | {delta_mae} | {delta_ratio} | {'yes' if row.material_improvement else 'no'} | {row.recommended_policy} | {row.notes} |"
        )
    lines.extend([
        "",
        f"Audited benchmarks: {len(audits)}",
        f"Material improvements: {sum(1 for row in audits if row.material_improvement)}",
    ])
    return "\n".join(lines) + "\n"


def render_benchmark_targets_markdown(
    snapshots: Iterable[BenchmarkTargetSnapshot],
    *,
    excluded_benchmark_ids: Optional[Iterable[str]] = None,
) -> str:
    rows = list(snapshots)
    excluded = list(excluded_benchmark_ids or [])
    lines = [
        "# Benchmark Targets",
        "",
        "| Benchmark | Target | Type | Roles | Panel Class | Evidence State | Proxy ppb | Observable ppb | Obs/Proxy | Matrix | Headspace | Class | Span | Depth | Headspace Class | Kaw 25C | Henry Name |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for row in rows:
        kaw = f"{row.henry_kaw_25c:.3e}" if row.henry_kaw_25c is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.target_name} | {row.target_type} | {', '.join(row.roles)} | {row.target_class} | {row.evidence_state} | {row.proxy_ppb:.3f} | {row.predicted_ppb:.3f} | {row.observable_ratio:.3f} | {row.matrix_factor:.3f} | {row.headspace_factor:.3f} | {row.volatile_class} | {row.span:.3f} | {row.depth} | {row.headspace_class} | {kaw} | {row.henry_source_name or 'n/a'} |"
        )

    lines.extend([
        "",
        f"Target rows: {len(rows)}",
        f"Low-headspace rows: {sum(1 for row in rows if row.headspace_class == 'low_headspace')}",
    ])
    if excluded:
        lines.extend([
            f"Excluded matrix-only benchmarks: {', '.join(sorted(excluded))}",
            "These benchmarks remain executable through summary/index artefacts, but they are deliberately omitted from target snapshots because they do not run through the free-precursor FAST target-ranking path.",
        ])
    return "\n".join(lines) + "\n"


def render_benchmark_index_markdown(entries: Iterable[BenchmarkIndexEntry]) -> str:
    rows = list(entries)
    lines = [
        "# Benchmark Index",
        "",
        "| Benchmark | Tier | Family | Protein | Process State | Execution Path | Engine | Cantera Role | Thermo Policy | Ranking Contract | Supported | Status | Strict Ready | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for entry in rows:
        notes = entry.reason or "indexed"
        lines.append(
            f"| {entry.benchmark_id} | {entry.tier} | {entry.family} | {entry.protein_type} | {entry.process_state or 'n/a'} | {entry.execution_path} | {entry.benchmark_engine} | {entry.cantera_role} | {entry.thermodynamic_gating_policy} | {entry.ranking_contract_status} | {'yes' if entry.supported else 'no'} | {entry.status} | {'yes' if entry.strict_ready else 'no'} | {notes} |"
        )
    return "\n".join(lines) + "\n"


def render_benchmark_summary_markdown(summaries: Iterable[BenchmarkSummary]) -> str:
    rows = list(summaries)
    lines = [
        "# Benchmark Summary",
        "",
        "| Benchmark | Tier | Family | Chemistry Families | Payload Roles | Protein | Process State | Execution Path | Engine | Cantera Role | Thermo Policy | Ranking Contract | Status | Strict Ready | Coverage | Pearson R | Max Ratio | Mean |log10 ratio| | MAE ppb | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for summary in rows:
        if not summary.supported:
            notes = summary.reason or "Unsupported"
            pearson = "n/a"
            max_ratio = "n/a"
            mean_log_error = "n/a"
            mae = "n/a"
            coverage = "0.0%"
            strict_ready = "no"
        else:
            notes = ", ".join(summary.blocking_issues) or "validated"
            pearson = f"{summary.pearson_r:.3f}" if summary.pearson_r is not None else "n/a"
            max_ratio = f"{summary.max_ratio:.3f}" if summary.max_ratio is not None else "n/a"
            mean_log_error = f"{summary.mean_abs_log10_error:.3f}" if summary.mean_abs_log10_error is not None else "n/a"
            mae = f"{summary.mae_ppb:.2f}" if summary.mae_ppb is not None else "n/a"
            coverage = f"{summary.coverage:.1%}"
            strict_ready = "yes" if summary.strict_ready else "no"

        lines.append(
            f"| {summary.benchmark_id} | {summary.tier} | {summary.family} | {', '.join(summary.chemistry_families) or 'none'} | {', '.join(summary.payload_roles) or 'none'} | {summary.protein_type} | {summary.process_state or 'n/a'} | {summary.execution_path} | {summary.benchmark_engine} | {summary.cantera_role} | {summary.thermodynamic_gating_policy} | {summary.ranking_contract_status} | {summary.overall_status} | {strict_ready} | {coverage} | {pearson} | {max_ratio} | {mean_log_error} | {mae} | {notes} |"
        )

    supported_count = sum(1 for summary in rows if summary.supported)
    pass_count = sum(1 for summary in rows if summary.overall_status in {"pass", "pass-no-ranking", "partial-pass"})
    strict_ready_count = sum(1 for summary in rows if summary.strict_ready)
    lines.extend([
        "",
        f"Supported benchmarks: {supported_count}/{len(rows)}",
        f"Benchmarks without blocking coverage/ranking gaps: {pass_count}/{len(rows)}",
        f"Strict-ready benchmarks: {strict_ready_count}/{len(rows)}",
    ])
    return "\n".join(lines) + "\n"


def render_matrix_benchmark_assertions_markdown(rows: Iterable[MatrixBenchmarkAssertion]) -> str:
    lines = [
        "# Matrix Benchmark Assertion Report",
        "",
        "| Benchmark | Protein | Path | State | Profile | Contract | Coverage | Min | Top-K | Hits | Top-K Status | Adverse | Max Ratio | Tol | Ratio Status | Status | Blocked |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- | --- |",
    ]
    for row in rows:
        max_ratio = f"{row.max_ratio:.3f}" if row.max_ratio is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.protein_type} | {row.execution_path} | {row.process_state or 'n/a'} | {row.target_profile} | {row.ranking_contract_status} | {row.coverage:.1%} | {row.min_coverage:.1%} | {row.top_k} | {row.top_k_hits} | {row.top_k_status} | {row.adverse_order_status} | {max_ratio} | {row.ratio_tolerance:.2f} | {row.ratio_status} | {row.overall_status} | {'yes' if row.strict_gate_blocked else 'no'} |"
        )
    return "\n".join(lines) + "\n"


def render_matrix_benchmark_evidence_markdown(rows: Iterable[MatrixBenchmarkEvidence]) -> str:
    lines = [
        "# Matrix Benchmark Evidence Audit Report",
        "",
        "| Benchmark | Protein | Path | State | Origin | Source | Target Profile | Data Status | Promotable | Blocker |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row.benchmark_id} | {row.protein_type} | {row.execution_path} | {row.process_state or 'n/a'} | {row.source_origin} | {row.source_reference} | {row.target_profile} | {row.external_data_status} | {'yes' if row.promotable else 'no'} | {row.promotion_blocker or 'none'} |"
        )
    return "\n".join(lines) + "\n"


def render_matrix_target_status_markdown(payload: Mapping[str, Any]) -> str:
    benchmark_rows = list(payload.get("benchmarks", []))
    lines = [
        "# Matrix Target Status",
        "",
        "| Benchmark | Protein | Path | Process State | Target Profile | Ref Signal | Evidence Origin | Quant Closed | Internal | Directional | Open | External Decision Ready | Evidence/Calibration Priority | Mechanistic Priority | Claim Posture | Blocker Class | Next Action | Best Computational Action |",
        "| --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in benchmark_rows:
        counts = row.get("support_counts", {})
        lines.append(
            f"| {row['benchmark_id']} | {row['protein_type']} | {row['execution_path']} | {row.get('process_state') or 'n/a'} | {row.get('target_profile', 'unknown')} | {row.get('reference_signal_origin', 'unknown')} | {row.get('source_origin', 'unknown')} | {int(counts.get('quantitative_closed', 0))} | {int(counts.get('internal_candidate', 0))} | {int(counts.get('directional_support', 0))} | {int(counts.get('open_gap', 0))} | {'yes' if row.get('promotion_ready', False) else 'no'} | {'yes' if row.get('evidence_or_calibration_priority_ready', False) else 'no'} | {'yes' if row.get('mechanistic_priority_ready', False) else 'no'} | {row.get('promotion_claim_posture', 'unknown')} | {row.get('blocker_class', 'unknown')} | {row.get('next_best_action', 'unknown')} | {row.get('best_computational_action', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Compound-Level Support",
        "",
        "| Benchmark | Compound | Role | Panel Class | Evidence State | Calibration Evidence | Support Status |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in benchmark_rows:
        for compound in row.get("compounds", []):
            lines.append(
                f"| {row['benchmark_id']} | {compound['compound']} | {compound['role']} | {compound.get('target_class', 'unknown')} | {compound.get('evidence_state', 'still_missing')} | {compound.get('calibration_evidence_strength', 'heuristic')} | {compound.get('support_status', 'open_gap')} |"
            )

    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Benchmarks covered: {int(summary.get('total_benchmarks', 0))}",
        f"Quantitative-support-ready benchmarks: {int(summary.get('quantitative_support_ready', 0))}",
        f"External-decision-ready benchmarks: {int(summary.get('promotion_ready', 0))}",
        f"Evidence/calibration-priority benchmarks: {int(summary.get('evidence_or_calibration_priority_ready', 0))}",
        f"Mechanistic-priority benchmarks: {int(summary.get('mechanistic_priority_ready', 0))}",
    ])
    return "\n".join(lines) + "\n"


def render_matrix_promotion_family_status_markdown(rows: Iterable[MatrixPromotionFamilyStatus]) -> str:
    lines = [
        "# Matrix Promotion Family Readiness",
        "",
        "| Protein | Off-Flavour Anchors | Meaty Candidates | External Meaty Anchors | Candidate Set Ready | External Unlocked | Blocker |",
        "| --- | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row.protein_type} | {row.off_flavour_anchor_count} | {row.meaty_candidate_count} | {row.external_meaty_anchor_count} | {'yes' if row.candidate_set_ready else 'no'} | {'yes' if row.external_assessment_unlocked else 'no'} | {row.blocker or 'none'} |"
        )
    return "\n".join(lines) + "\n"


def render_matrix_promotion_contract_markdown(payload: Mapping[str, Any]) -> str:
    rule = payload.get("promotion_rule", {})
    selected = payload.get("selected_promotion_target") or {}
    lines = [
        "# Matrix Promotion Contract",
        "",
        f"Contract id: {rule.get('contract_id', 'unknown')}",
        f"Minimum quantitative closed targets: {int(rule.get('minimum_quantitative_closed_targets', 0))}",
        f"Requires measured_volatiles: {'yes' if rule.get('requires_measured_volatiles', False) else 'no'}",
        f"Requires external quantitative origin: {'yes' if rule.get('requires_external_quantitative_origin', False) else 'no'}",
        f"Requires mixed or meaty-positive target profile: {'yes' if rule.get('requires_mixed_or_meaty_positive_target_profile', False) else 'no'}",
        f"Requires passing ranking contract: {'yes' if rule.get('requires_passing_ranking_contract', False) else 'no'}",
        f"Disallow internal-candidate support: {'yes' if rule.get('disallow_internal_candidate_support', False) else 'no'}",
        f"Disallow directional support: {'yes' if rule.get('disallow_directional_support', False) else 'no'}",
        "",
        "## Rule Notes",
        "",
    ]
    for note in rule.get("notes", []):
        lines.append(f"- {note}")

    if selected:
        lines.extend([
            "",
            "## Selected Promotion Target",
            "",
            f"- benchmark: {selected.get('benchmark_id', 'unknown')}",
            f"- protein: {selected.get('protein_type', 'unknown')}",
            f"- process_state: {selected.get('process_state', 'n/a')}",
            f"- target_profile: {selected.get('target_profile', 'unknown')}",
            f"- selection_policy: {selected.get('selection_policy', 'unknown')}",
        ])
        for rationale in selected.get("rationale", []):
            lines.append(f"- rationale: {rationale}")

    lines.extend([
        "",
        "## Benchmark Assessment",
        "",
        "| Benchmark | Protein | Process State | Target Profile | Promotion Ready | Blocker | Passed Requirements |",
        "| --- | --- | --- | --- | --- | --- | ---: |",
    ])
    for row in payload.get("benchmarks", []):
        passed = sum(1 for requirement in row.get("requirements", []) if requirement.get("passed"))
        lines.append(
            f"| {row.get('benchmark_id', 'unknown')} | {row.get('protein_type', 'unknown')} | {row.get('process_state', 'n/a')} | {row.get('target_profile', 'unknown')} | {'yes' if row.get('promotion_ready', False) else 'no'} | {row.get('promotion_blocker', 'unknown')} | {passed}/{len(row.get('requirements', []))} |"
        )

    lines.extend([
        "",
        "## Requirement Details",
        "",
        "| Benchmark | Requirement | Passed | Detail |",
        "| --- | --- | --- | --- |",
    ])
    for row in payload.get("benchmarks", []):
        for requirement in row.get("requirements", []):
            lines.append(
                f"| {row.get('benchmark_id', 'unknown')} | {requirement.get('label', 'unknown')} | {'yes' if requirement.get('passed', False) else 'no'} | {requirement.get('detail', 'unknown')} |"
            )
    return "\n".join(lines) + "\n"


def render_matrix_observable_closure_audit_markdown(payload: Mapping[str, Any]) -> str:
    selected = payload.get("selected_promotion_target") or {}
    lines = [
        "# Matrix Observable Closure Audit",
        "",
    ]
    if selected:
        lines.extend([
            f"Selected promotion target: {selected.get('benchmark_id', 'unknown')}",
            f"Selection policy: {selected.get('selection_policy', 'unknown')}",
            "",
        ])
    lines.extend([
        "| Benchmark | Protein | Process State | Target Profile | Promotion Blocker | Action Counts |",
        "| --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload.get("benchmarks", []):
        counts = row.get("closure_action_counts", {})
        rendered_counts = "; ".join(f"{key}={value}" for key, value in sorted(counts.items())) or "none"
        lines.append(
            f"| {row.get('benchmark_id', 'unknown')} | {row.get('protein_type', 'unknown')} | {row.get('process_state', 'n/a')} | {row.get('target_profile', 'unknown')} | {row.get('promotion_blocker', 'unknown')} | {rendered_counts} |"
        )

    lines.extend([
        "",
        "## Compound Closure Actions",
        "",
        "| Benchmark | Compound | Role | Panel Class | Evidence State | Calibration Evidence | Support Status | Closure Action |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload.get("benchmarks", []):
        for compound in row.get("compounds", []):
            lines.append(
                f"| {row.get('benchmark_id', 'unknown')} | {compound.get('compound', 'unknown')} | {compound.get('role', 'unknown')} | {compound.get('target_class', 'unknown')} | {compound.get('evidence_state', 'unknown')} | {compound.get('calibration_evidence_strength', 'unknown')} | {compound.get('support_status', 'unknown')} | {compound.get('closure_action', 'unknown')} |"
            )

    summary = payload.get("summary", {})
    rendered_summary = "; ".join(
        f"{key}={value}" for key, value in sorted((summary.get("closure_action_counts") or {}).items())
    ) or "none"
    lines.extend([
        "",
        f"Benchmarks audited: {int(summary.get('benchmarks_audited', 0))}",
        f"Mechanistic watchlist count: {int(summary.get('mechanistic_watchlist_count', 0))}",
        f"Closure action counts: {rendered_summary}",
    ])
    watchlist = payload.get("mechanistic_refinement_watchlist", [])
    if watchlist:
        lines.extend([
            "",
            "## Mechanistic Refinement Watchlist",
            "",
            "| Benchmark | Protein | Process State | Target Compounds | Expected Decision Change | Scope |",
            "| --- | --- | --- | --- | --- | --- |",
        ])
        for row in watchlist:
            lines.append(
                f"| {row.get('benchmark_id', 'unknown')} | {row.get('protein_type', 'unknown')} | {row.get('process_state', 'n/a')} | {', '.join(str(item) for item in row.get('target_compounds', [])) or 'none'} | {row.get('expected_decision_change', 'unknown')} | {row.get('allowed_scope', 'unknown')} |"
            )
    return "\n".join(lines) + "\n"


def render_matrix_experiment_intake_schema_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Experiment Intake Schema",
        "",
        f"Contract id: {payload.get('contract_id', 'unknown')}",
        f"Schema version: {payload.get('schema_version', 'unknown')}",
        f"Description: {payload.get('description', 'unknown')}",
        "",
        "## Required Fields",
        "",
        f"Top-level: {', '.join(str(item) for item in payload.get('required_top_level_fields', [])) or 'none'}",
        f"Conditions: {', '.join(str(item) for item in payload.get('required_conditions', [])) or 'none'}",
        f"Formulation: {', '.join(str(item) for item in payload.get('required_formulation_fields', [])) or 'none'}",
        f"Provenance: {', '.join(str(item) for item in payload.get('required_provenance_fields', [])) or 'none'}",
        "",
        "## Allowed Values",
        "",
        f"Source kinds: {', '.join(str(item) for item in payload.get('allowed_source_kinds', [])) or 'none'}",
        f"Protein types: {', '.join(str(item) for item in payload.get('allowed_protein_types', [])) or 'none'}",
        "",
        "## Policies",
        "",
    ]
    for note in payload.get("policies", []):
        lines.append(f"- {note}")
    return "\n".join(lines) + "\n"


def render_matrix_experiment_support_delta_markdown(payload: Mapping[str, Any]) -> str:
    experiment = payload.get("experiment", {})
    promotion = payload.get("promotion_assessment", {})
    support_delta = payload.get("support_delta", {})
    aligned = payload.get("aligned_benchmark", {})
    lines = [
        "# Matrix Experiment Support Delta",
        "",
        f"Experiment id: {experiment.get('experiment_id', 'unknown')}",
        f"Source kind: {experiment.get('source_kind', 'unknown')}",
        f"Protein: {experiment.get('protein_type', 'unknown')}",
        f"Process state: {experiment.get('process_state', 'unknown')}",
        f"Source reference: {experiment.get('source_reference', 'unknown')}",
        f"Aligned benchmark: {aligned.get('benchmark_id', 'none')}",
        "",
        "## Promotion Delta",
        "",
        f"Promotion ready before: {'yes' if promotion.get('promotion_ready_before', False) else 'no'}",
        f"Promotion ready after: {'yes' if promotion.get('promotion_ready_after', False) else 'no'}",
        f"Promotion claim allowed: {'yes' if promotion.get('promotion_claim_allowed', False) else 'no'}",
        f"Blocker before: {promotion.get('promotion_blocker_before', 'unknown')}",
        f"Blocker after: {promotion.get('promotion_blocker_after', 'unknown')}",
        f"Promotion claim policy: {promotion.get('promotion_claim_policy', 'unknown')}",
        f"Readiness change: {promotion.get('readiness_change', 'unknown')}",
        f"Landing recommendation: {promotion.get('landing_recommendation', 'unknown')}",
        "",
        "## Support Counts",
        "",
        f"Baseline: {support_delta.get('baseline_support_counts', {})}",
        f"Current: {support_delta.get('current_support_counts', {})}",
        f"Delta: {support_delta.get('delta_counts', {})}",
        "",
        "## Compound Support Delta",
        "",
        "| Compound | Role | Evidence State | Support Before | Support After | Delta | Ratio |",
        "| --- | --- | --- | --- | --- | --- | ---: |",
    ]
    for row in payload.get("compounds", []):
        lines.append(
            f"| {row.get('compound', 'unknown')} | {row.get('role', 'unknown')} | {row.get('evidence_state', 'unknown')} | {row.get('baseline_support_status', 'none')} | {row.get('support_status', 'unknown')} | {row.get('support_delta', 'unknown')} | {float(row.get('ratio', 0.0)):.3f} |"
        )
    return "\n".join(lines) + "\n"
