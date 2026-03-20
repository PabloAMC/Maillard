from __future__ import annotations

from typing import Any, Iterable, Mapping

import math
from typing import Any, Dict, Iterable, List, Optional, Tuple, TYPE_CHECKING

from src.pipeline import FormulationResult
from src.usability_reports import ValidatedEnvelopeReport, DomainWarning
from src.benchmark_types import (
    MatrixBenchmarkBranchDelta, 
    MatrixBenchmarkDelta, 
    ThermodynamicGatingAudit, 
    BenchmarkTargetSnapshot, 
    BenchmarkIndexEntry, 
    BenchmarkSummary,
    MatrixBenchmarkAssertion,
    MatrixPromotionFamilyStatus
)


def render_validated_envelope_markdown(report: 'ValidatedEnvelopeReport') -> str:
    lines = [
        "# Validated Envelope",
        "",
        f"Target tag: {report.target_tag}",
        f"Supported benchmarks: {report.supported_benchmarks}/{report.total_benchmarks}",
        f"Strict-ready benchmarks: {', '.join(report.strict_ready_benchmarks) if report.strict_ready_benchmarks else 'none'}",
        f"Matrix-only executable benchmarks: {', '.join(report.matrix_only_benchmarks) if report.matrix_only_benchmarks else 'none'}",
        "",
        "## Warnings",
    ]
    for warning in report.warnings:
        lines.append(f"- {warning}")
    lines.extend([
        "",
        "## Next Priorities",
    ])
    for item in report.next_priorities:
        lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def render_formulation_explainability_markdown(payload: Dict[str, object]) -> str:
    lines = [
        "# Formulation Explainability",
        "",
        f"Formulation: {payload['formulation_name']}",
        f"Protein type: {payload['protein_type']}",
        f"Target tag: {payload['target_tag']}",
        f"Minimize tag: {payload['minimize_tag']}",
        "",
        "## Matrix State",
    ]
    matrix = payload["matrix_explainability"]
    lines.extend([
        f"- Effective denaturation state: {matrix['effective_denaturation_state']:.3f}",
        f"- Lysine accessibility: {matrix['lysine_accessibility']:.3f}",
        f"- Cysteine accessibility: {matrix['cysteine_accessibility']:.3f}",
        f"- Bulk volatile retention: {matrix['bulk_volatile_retention']:.3f}",
        f"- Denaturation source: {matrix['denaturation_source']}",
        "",
        "## Scores",
    ])
    for key, value in payload["scores"].items():
        lines.append(f"- {key}: {value:.3f}")
    
    lines.append(render_projection_rows_markdown(payload["top_projection_rows"], heading="## Projection Rows", variant="detailed").rstrip())
    
    flavor_axis = payload.get("flavor_axis_summary", {})
    if flavor_axis:
        lines.append(render_flavor_axis_markdown(flavor_axis, heading="## Flavor Axis", variant="compact").rstrip())
    
    return "\n".join(lines) + "\n"


def render_projection_rows_markdown(
    rows: List[Dict[str, object]],
    *,
    heading: str,
    variant: str = "detailed",
) -> str:
    lines = [heading, ""]
    if variant == "compact":
        lines.extend([
            "| Compound | Panel Class | Evidence State | Reachability | Process State | Retention | Calibration Source | Observable Assumption | Evidence | Browning | Observable ppb |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | ---: |",
        ])
        for row in rows:
            lines.append(
                f"| {row['compound']} | {row.get('target_class', 'unknown')} | {row.get('evidence_state', 'still_missing')} | {row.get('reachability_status', 'merely_plausible')} | {row['process_state']} | {row['retention_runtime_mode']} | {row['calibration_source']} | {row.get('observable_assumption_summary', 'unknown')} | {row['calibration_evidence_strength']} | {float(row.get('browning_index', 0.0)):.2f} | {float(row['observable_ppb']):.2f} |"
            )
        lines.append("")
        return "\n".join(lines)

    lines.extend([
        "| Compound | Proxy ppb | Observable ppb | Obs/Proxy | Matrix | Dynamic | Headspace | Class | Panel Class | Evidence State | Reachability | Support Origin | Process | Retention | Calibration | Observable Assumption | Evidence | Fallback |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in rows:
        lines.append(
            f"| {row['compound']} | {float(row['proxy_ppb']):.3f} | {float(row['observable_ppb']):.3f} | {float(row['observable_ratio']):.3f} | {float(row['matrix_factor']):.3f} | {float(row.get('dynamic_retention_factor', 1.0)):.3f} | {float(row['headspace_factor']):.3f} | {row['volatile_class']} | {row.get('target_class', 'unknown')} | {row.get('evidence_state', 'still_missing')} | {row.get('reachability_status', 'merely_plausible')} | {row.get('support_origin', 'standard_matrix_support')} | {row['process_state']} | {row.get('retention_runtime_mode', 'static_class_profile')} | {row['calibration_source']} | {row.get('observable_assumption_summary', 'unknown')} | {row['calibration_evidence_strength']} | {row['calibration_fallback_mode']} |"
        )
    lines.append("")
    return "\n".join(lines)


def render_flavor_axis_markdown(
    flavor_axis: Dict[str, object],
    *,
    heading: str,
    variant: str = "detailed",
) -> str:
    lines = [heading]
    if variant == "compact":
        lines.extend([
            f"- Strecker balance score: {float(flavor_axis.get('strecker_balance_score', 0.0)):.3f}",
            f"- Pyrazine burden: {float(flavor_axis.get('pyrazine_burden', 0.0)):.3f}",
            f"- Thiamine pathway active: {flavor_axis.get('thiamine_pathway_active', False)}",
            f"- Thiamine source: {flavor_axis.get('thiamine_availability_source', 'unknown')}",
            "",
        ])
        return "\n".join(lines)

    lines.extend([
        f"- **strecker_balance_score:** {float(flavor_axis.get('strecker_balance_score', 0.0)):.2f}",
        f"- **strecker_gap_penalty:** {float(flavor_axis.get('strecker_gap_penalty', 0.0)):.2f}",
        f"- **pyrazine_signal_ppb:** {float(flavor_axis.get('pyrazine_signal_ppb', 0.0)):.2f}",
        f"- **pyrazine_propensity:** {float(flavor_axis.get('pyrazine_propensity', 0.0)):.2f}",
        f"- **pyrazine_burden:** {float(flavor_axis.get('pyrazine_burden', 0.0)):.2f}",
        f"- **pyrazine_penalty:** {float(flavor_axis.get('pyrazine_penalty', 0.0)):.2f}",
        f"- **furanone_support_score:** {float(flavor_axis.get('furanone_support_score', 1.0)):.2f}",
        f"- **furanone_penalty:** {float(flavor_axis.get('furanone_penalty', 0.0)):.2f}",
        f"- **thiamine_pathway_active:** {flavor_axis.get('thiamine_pathway_active', False)}",
        f"- **thiamine_availability_source:** {flavor_axis.get('thiamine_availability_source', 'unknown')}",
        f"- **thiamine_availability_explicit:** {flavor_axis.get('thiamine_availability_explicit', False)}",
        f"- **thiamine_provenance_mode:** {flavor_axis.get('thiamine_provenance_mode', 'inactive')}",
    ])
    expected_furanones = flavor_axis.get("furanone_expected", [])
    if expected_furanones:
        lines.append(f"- **furanone_expected:** {', '.join(str(item) for item in expected_furanones)}")
    observed_furanones = flavor_axis.get("furanone_observed", [])
    if observed_furanones:
        lines.append(f"- **furanone_observed:** {', '.join(str(item) for item in observed_furanones)}")
    missing_furanones = flavor_axis.get("furanone_missing", [])
    if missing_furanones:
        lines.append(f"- **furanone_missing:** {', '.join(str(item) for item in missing_furanones)}")
    lines.append(f"- **lincoln_crosstalk_prior:** {flavor_axis.get('lincoln_crosstalk_prior', {}).get('summary', '')}")
    lines.append("")
    return "\n".join(lines)


def render_domain_warnings_markdown(warnings: List['DomainWarning']) -> str:
    if not warnings:
        return "> [!NOTE]\n> Run is within the validated scientific envelope.\n"
        
    lines = ["## Scientific Domain Warnings", ""]
    for w in warnings:
        icon = "⚠️" if w.level == "WARNING" else "🚨"
        lines.append(f"- **{icon} {w.level} [{w.category}]**: {w.message}")
    return "\n".join(lines) + "\n"


def render_domain_warnings_cli(warnings: List['DomainWarning']):
    if not warnings:
        print("\n  ✅ Scientific Envelope: Trusted (matches PRIMARY benchmarks)")
        return

    print("\n  ⚠️  SCIENTIFIC DOMAIN WARNINGS:")
    for w in warnings:
        icon = "!" if w.level == "WARNING" else "!!"
        print(f"    [{icon}] {w.category.ljust(15)}: {w.message}")
    print("  " + "─" * 60)


def render_decision_summary_cli(result: 'FormulationResult', warnings: List['DomainWarning']):
    """
    Renders a premium, consolidated scientist-facing summary of the simulation run.
    """
    print("\n" + "═" * 80)
    print(" " * 25 + "📊 MAILLARD DECISION SUMMARY")
    print("═" * 80)

    confidence = getattr(result, "confidence_metadata", {}) or {}

    # ── Prediction Posture Banner ─────────────────────────────────────────────
    prediction_mode = confidence.get("prediction_mode", "")
    if prediction_mode == "benchmark_supported_quantitative":
        banner_icon, banner_label, banner_level = "✅", "QUANTITATIVE MODE", "Results are calibrated against PRIMARY benchmarks. Suitable for formulation triage."
    elif prediction_mode == "ranking_supported":
        banner_icon, banner_label, banner_level = "⚠️", "DIRECTIONAL MODE", "Results are directionally supported. Verify absolute concentrations experimentally."
    elif prediction_mode == "directional_only":
        banner_icon, banner_label, banner_level = "🟡", "DIRECTIONAL ONLY", "Confidence is low. Use for hypothesis generation, not decision-grade prioritization."
    else:
        banner_icon, banner_label, banner_level = "🔬", "EXPLORATORY MODE", "Results are speculative. No PRIMARY benchmark covers this formulation space."
    print(f"\n  {banner_icon} {banner_label}")
    print(f"      {banner_level}")
    print("  " + "─" * 76)

    if confidence:
        tier = str(confidence.get("tier", "unknown")).upper()
        score = float(confidence.get("score", 0.0))
        neighborhood = confidence.get("benchmark_neighborhood", "unknown")
        process_regime = confidence.get("process_regime", "unknown")
        posture = confidence.get("recommended_posture", "")
        decision_mode = confidence.get("decision_mode", "directional_hypothesis")
        print(f"\n  [0] DECISION CONFIDENCE: {tier} ({score:.0f}/100)")
        print(f"      Benchmark Basis  : {neighborhood}")
        print(f"      Decision Mode    : {decision_mode}")
        if process_regime not in {"unknown", "free_aqueous", "matrix_hydrated"}:
            print(f"      Process Regime   : {process_regime}")
        print(f"      Prediction Mode : {prediction_mode}")
        if process_regime in {"extrusion_like", "extrusion_heavy"}:
            panel = confidence.get("extrusion_observable_panel", {}) or {}
            if panel:
                panel_summary = []
                for category, label in (("meaty_positive", "meaty"), ("off_notes", "off-notes"), ("severity_markers", "severity")):
                    row = panel.get(category, {})
                    panel_summary.append(f"{label} {int(row.get('present_count', 0))}/{int(row.get('required_count', 0))}")
                readiness = "ready" if panel.get("minimum_panel_ready", False) else "incomplete"
                print(f"      Extrusion Panel  : {' | '.join(panel_summary)} ({readiness})")
        # P1: show accessibility profile when matrix formulation
        matrix_expl = getattr(result, "matrix_explainability", None)
        if matrix_expl and isinstance(matrix_expl, dict) and matrix_expl.get("accessibility_profile"):
            acc_profile = matrix_expl["accessibility_profile"]
            acc_warning = matrix_expl.get("accessibility_warning", False)
            acc_icon = "⚠️" if acc_warning else "✅"
            print(f"      Accessibility    : {acc_profile} {acc_icon}")
        if posture:
            print(f"      Recommended Use  : {posture}")
        for factor in confidence.get("dominant_factors", [])[:2]:
            print(f"      Why              : {factor}")
        calibration = confidence.get("calibration_diagnostics", {})
        if calibration:
            print(f"      Calibration      : {calibration.get('summary', '')}")
            axes = calibration.get("extrapolation_axes", [])
            if axes:
                print(f"      Extrapolation    : {', '.join(str(axis) for axis in axes)}")


    # 1. Scientific Envelope Section
    status_icon = "✅" if not warnings else "⚠️"
    status_label = "TRUSTED" if not warnings else "LIMITED"
    print(f"\n  [1] SCIENTIFIC ENVELOPE: {status_label} {status_icon}")
    if warnings:
        for w in warnings:
            dot = "!" if w.level == "WARNING" else "!!"
            print(f"    [{dot}] {w.category.ljust(15)}: {w.message}")
    
    # 2. Top-Line Metrics
    print("\n  [2] TOP-LINE METRICS:")
    print(f"      Target Score     : {float(result.target_score):.2f}")
    print(f"      Off-Flavour Risk : {float(result.off_flavour_risk):.2f}")
    print(f"      Safety Score     : {float(result.safety_score):.2f}")
    
    print("\n" + "═" * 80 + "\n")


def render_deep_explainability_cli(result: 'FormulationResult'):
    """
    Prints deep scientific explainability panels for the CLI.
    """
    # Panel A: Bottleneck Analysis
    print(f"  [A] BOTTLENECK ANALYSIS:")
    b_name = getattr(result, 'bottleneck_precursor', 'none').upper()
    sev = getattr(result, 'bottleneck_severity', 0.0)
    
    status = "OPTIMAL" if sev < 0.2 else ("CONSTRAINED" if sev < 0.6 else "SEVERE BOTTLENECK")
    print(f"      Limiting Precursor : {b_name} ({status})")
    
    # Average Matrix Retention Loss
    target_metadata = []
    for canon, meta in result.projection_metadata.items():
        if result.predicted_ppb.get(canon, 0) > 0:
            target_metadata.append(meta)
            
    avg_loss = 0.0
    if target_metadata:
        avg_loss = 100.0 * (1.0 - sum(float(m.get("matrix_factor", 1.0)) for m in target_metadata) / len(target_metadata))

    print(f"      Avg. Matrix Loss   : {avg_loss:.1f}% potential flavor trapped")

    # Panel B: Precursor Attribution
    print(f"\n  [B] PRECURSOR ATTRIBUTION (Yield Drivers):")
    attrib = getattr(result, 'precursor_contributions', {})
    if attrib:
        total_attrib = sum(attrib.values())
        sorted_attrib = sorted(attrib.items(), key=lambda x: x[1], reverse=True)[:3]
        for name, val in sorted_attrib:
            share = 100.0 * val / total_attrib if total_attrib > 0 else 0.0
            print(f"      {name.ljust(18)} : {share:.1f}% of observable yield")
    else:
        print(f"      No attribution data available.")

    # Panel C: Physical Suppression
    suppressed = getattr(result, 'suppressed_compounds', [])
    if suppressed:
        print(f"\n  [C] PHYSICAL SUPPRESSION (Lost Potential):")
        for item in suppressed[:2]:
            cause = item['primary_cause'].upper()
            reduction = item['reduction_factor'] * 100.0
            print(f"      {item['name'].ljust(18)} : -{reduction:.1f}% yield due to {cause}")

    # Panel D: Performance Drivers
    print(f"\n  [D] KEY AROMA DRIVERS:")
    top_targets = sorted(
        [(m.get("compound", "unknown"), m.get("observable_ppb", 0.0)) for m in result.projection_metadata.values()],
        key=lambda x: x[1],
        reverse=True
    )[:2]
    
    if top_targets:
        targets_str = ", ".join([f"{t[0]} ({t[1]:.1f} ppb)" for t in top_targets])
        print(f"      Top Yield Drivers  : {targets_str}")

    from src.projection_utils import generate_intervention_hint
    print(f"\n  [💡] INTERVENTION HINT:")
    hint = generate_intervention_hint(result)
    print(f"      {hint}")

    confidence = getattr(result, "confidence_metadata", {}) or {}
    compound_rows = confidence.get("compound_confidence", [])
    if compound_rows:
        print(f"\n  [E] COMPOUND CONFIDENCE:")
        for row in compound_rows[:3]:
            mode_label = str(row.get('prediction_mode', 'unknown')).replace('_', ' ').upper()
            cal_src = str(row.get('calibration_source', 'unknown'))
            cal_evid = str(row.get('calibration_evidence_strength', 'heuristic'))
            print(
                f"      {row['compound'].ljust(22)}: {str(row['tier']).upper()} ({float(row['score']):.0f}/100)  "
                f"obs {float(row['observable_ppb']):.1f} ppb  [{mode_label}]"
            )
            print(
                f"      {'':22}  class: {str(row.get('target_class', 'unknown'))} | state: {str(row.get('evidence_state', 'still_missing'))}"
            )
            print(
                f"      {'':22}  cal: {cal_src} | evidence: {cal_evid}"
            )

    # [G] Accessibility warning — fires when matrix trapping dominates a top compound
    accessibility_dominated = [
        row for row in compound_rows
        if row.get("accessibility_dominated") and float(row.get("observable_ppb", 0.0)) > 0.0
    ]
    if accessibility_dominated:
        print(f"\n  [G] ACCESSIBILITY WARNING:")
        for row in accessibility_dominated[:2]:
            mf = float(row.get("matrix_factor", 1.0))
            hf = float(row.get("headspace_factor", 1.0))
            print(
                f"      {row['compound'].ljust(22)}: matrix_factor={mf:.2f}, headspace_factor={hf:.2f}"
            )
            print(
                f"      {'':22}  Accessibility (not chemistry) limits this compound's yield."
            )
            print(
                f"      {'':22}  → Increasing denaturation state or raising pH may improve release."
            )

    ranking_drivers = confidence.get("sensitivity_summary", {})
    sensitivity = ranking_drivers if isinstance(ranking_drivers, dict) else {}
    ranking_drivers = sensitivity.get("ranking_drivers", [])
    if ranking_drivers:
        print(f"\n  [F] SENSITIVITY SUMMARY:")
        for item in ranking_drivers[:2]:
            print(
                f"      {item['input'].ljust(18)} : Δdecision {float(item['decision_delta']):+.2f}, "
                f"Δsafety {float(item['safety_delta']):+.2f} ({item['perturbation']})"
            )


def render_matrix_branch_deltas_markdown(
    rows: Iterable['MatrixBenchmarkBranchDelta'],
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


def render_matrix_benchmark_deltas_markdown(rows: Iterable['MatrixBenchmarkDelta']) -> str:
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


def render_thermodynamic_gating_audit_markdown(rows: Iterable['ThermodynamicGatingAudit']) -> str:
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
    snapshots: Iterable['BenchmarkTargetSnapshot'],
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


def render_benchmark_index_markdown(entries: Iterable['BenchmarkIndexEntry']) -> str:
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


def render_benchmark_summary_markdown(summaries: Iterable['BenchmarkSummary']) -> str:
    rows = list(summaries)
    lines = [
        "# Benchmark Summary",
        "",
        "| Benchmark | Tier | Family | Protein | Process State | Execution Path | Engine | Cantera Role | Thermo Policy | Ranking Contract | Status | Strict Ready | Coverage | Pearson R | Max Ratio | Mean |log10 ratio| | MAE ppb | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
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
            f"| {summary.benchmark_id} | {summary.tier} | {summary.family} | {summary.protein_type} | {summary.process_state or 'n/a'} | {summary.execution_path} | {summary.benchmark_engine} | {summary.cantera_role} | {summary.thermodynamic_gating_policy} | {summary.ranking_contract_status} | {summary.overall_status} | {strict_ready} | {coverage} | {pearson} | {max_ratio} | {mean_log_error} | {mae} | {notes} |"
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


def render_provenance_markdown(provenance: Dict[str, Any]) -> str:
    repository = provenance.get("repository", {})
    scientific_surface = provenance.get("scientific_surface", {})
    lines = ["## 5. Provenance\n"]
    lines.append(f"- **artifact_kind:** {provenance.get('artifact_kind', 'unknown')}\n")
    lines.append(f"- **generated_at:** {provenance.get('generated_at', 'unknown')}\n")
    lines.append(f"- **generator:** {provenance.get('generator', {}).get('command', '')}\n")
    lines.append(
        f"- **repository:** {repository.get('name', 'unknown')} | branch {repository.get('branch', 'unknown')} | commit {repository.get('short_commit', 'unknown')} | dirty {repository.get('dirty', False)}\n"
    )
    lines.append(f"- **input_fingerprint_sha256:** {provenance.get('input_fingerprint_sha256', '')}\n")
    if scientific_surface:
        lines.append("- **scientific_surface:**\n")
        for key, value in scientific_surface.items():
            lines.append(f"  - {key}: {value}\n")
    campaign = provenance.get("campaign")
    if campaign:
        lines.append(f"- **campaign_name:** {campaign.get('name', 'unknown')}\n")
    lines.append("\n")
    return "".join(lines)


def render_matrix_benchmark_assertions_markdown(rows: Iterable['MatrixBenchmarkAssertion']) -> str:
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


def render_matrix_benchmark_evidence_markdown(rows: Iterable['MatrixBenchmarkEvidence']) -> str:
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
        "| Benchmark | Protein | Path | Process State | Target Profile | Ref Signal | Quant Closed | Internal | Directional | Open | External Decision Ready | Mechanistic Priority | Next Action |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for row in benchmark_rows:
        counts = row.get("support_counts", {})
        lines.append(
            f"| {row['benchmark_id']} | {row['protein_type']} | {row['execution_path']} | {row.get('process_state') or 'n/a'} | {row.get('target_profile', 'unknown')} | {row.get('reference_signal_origin', 'unknown')} | {int(counts.get('quantitative_closed', 0))} | {int(counts.get('internal_candidate', 0))} | {int(counts.get('directional_support', 0))} | {int(counts.get('open_gap', 0))} | {'yes' if row.get('promotion_ready', False) else 'no'} | {'yes' if row.get('mechanistic_priority_ready', False) else 'no'} | {row.get('next_best_action', 'unknown')} |"
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
        f"Mechanistic-priority benchmarks: {int(summary.get('mechanistic_priority_ready', 0))}",
    ])
    return "\n".join(lines) + "\n"


def render_matrix_promotion_family_status_markdown(rows: Iterable['MatrixPromotionFamilyStatus']) -> str:
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
