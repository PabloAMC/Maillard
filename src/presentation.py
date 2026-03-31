from __future__ import annotations

from typing import Any, Iterable, Mapping

import math
from typing import Any, Dict, Iterable, List, Optional, Tuple, TYPE_CHECKING

from src.pipeline import FormulationResult
from src.usability_reports import ValidatedEnvelopeReport, DomainWarning
from src.benchmark_markdown import render_matrix_branch_deltas_markdown
from src.benchmark_markdown import (
    render_benchmark_index_markdown,
    render_benchmark_summary_markdown,
    render_benchmark_targets_markdown,
    render_matrix_benchmark_assertions_markdown,
    render_matrix_benchmark_deltas_markdown,
    render_matrix_benchmark_evidence_markdown,
    render_matrix_experiment_intake_schema_markdown,
    render_matrix_experiment_support_delta_markdown,
    render_matrix_observable_closure_audit_markdown,
    render_matrix_promotion_contract_markdown,
    render_matrix_promotion_family_status_markdown,
    render_matrix_target_status_markdown,
    render_thermodynamic_gating_audit_markdown,
)
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
    def _format_regimes(row: Mapping[str, object]) -> str:
        regimes = row.get("modeling_regimes", [])
        if not isinstance(regimes, list) or not regimes:
            return "unknown"
        return ", ".join(str(item) for item in regimes)

    if variant == "compact":
        lines.extend([
            "| Compound | Panel Class | Role | Kind | Evidence State | Reachability | Process State | Retention | Calibration Source | Observable Assumption | Evidence | Browning | Observable ppb |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | ---: |",
        ])
        for row in rows:
            lines.append(
                f"| {row['compound']} | {row.get('target_class', 'unknown')} | {row.get('panel_role', 'unknown')} | {row.get('observable_kind', 'unknown')} | {row.get('evidence_state', 'still_missing')} | {row.get('reachability_status', 'merely_plausible')} | {row['process_state']} | {row['retention_runtime_mode']} | {row['calibration_source']} | {row.get('observable_assumption_summary', 'unknown')} | {row['calibration_evidence_strength']} | {float(row.get('browning_index', 0.0)):.2f} | {float(row['observable_ppb']):.2f} |"
            )
        lines.append("")
        return "\n".join(lines)

    lines.extend([
        "| Compound | Proxy ppb | Observable ppb | Obs/Proxy | Matrix | Dynamic | Headspace | Class | Panel Class | Panel Role | Observable Kind | Regimes | Evidence State | Reachability | Support Origin | Process | Retention | Calibration | Observable Assumption | Evidence | Fallback |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in rows:
        lines.append(
            f"| {row['compound']} | {float(row['proxy_ppb']):.3f} | {float(row['observable_ppb']):.3f} | {float(row['observable_ratio']):.3f} | {float(row['matrix_factor']):.3f} | {float(row.get('dynamic_retention_factor', 1.0)):.3f} | {float(row['headspace_factor']):.3f} | {row['volatile_class']} | {row.get('target_class', 'unknown')} | {row.get('panel_role', 'unknown')} | {row.get('observable_kind', 'unknown')} | {_format_regimes(row)} | {row.get('evidence_state', 'still_missing')} | {row.get('reachability_status', 'merely_plausible')} | {row.get('support_origin', 'standard_matrix_support')} | {row['process_state']} | {row.get('retention_runtime_mode', 'static_class_profile')} | {row['calibration_source']} | {row.get('observable_assumption_summary', 'unknown')} | {row['calibration_evidence_strength']} | {row['calibration_fallback_mode']} |"
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
    family_lane_summary = flavor_axis.get("family_lane_summary", {})
    active_family_lanes = [
        family_lane_summary[slr_family]
        for slr_family in flavor_axis.get("active_family_lanes", [])
        if slr_family in family_lane_summary
    ]
    if variant == "compact":
        glutathione_lane = family_lane_summary.get("05", {})
        lines.extend([
            f"- Strecker balance score: {float(flavor_axis.get('strecker_balance_score', 0.0)):.3f}",
            f"- Pyrazine burden: {float(flavor_axis.get('pyrazine_burden', 0.0)):.3f}",
            f"- Thiamine pathway active: {flavor_axis.get('thiamine_pathway_active', False)}",
            f"- Thiamine source: {flavor_axis.get('thiamine_availability_source', 'unknown')}",
        ])
        if glutathione_lane:
            lines.append(f"- Glutathione support active: {glutathione_lane.get('glutathione_active', False)}")
        if active_family_lanes:
            lines.append(f"- Active family lanes: {', '.join(str(row.get('display_name', row.get('slr_family', 'unknown'))) for row in active_family_lanes)}")
        lines.append("")
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
    glutathione_lane = family_lane_summary.get("05", {})
    if glutathione_lane:
        lines.append(f"- **glutathione_support_active:** {glutathione_lane.get('glutathione_active', False)}")
        lines.append(f"- **peptide_support_active:** {glutathione_lane.get('peptide_support_active', False)}")
        lines.append(f"- **sulfur_peptide_support_score:** {float(glutathione_lane.get('sulfur_peptide_support_score', 0.0)):.2f}")
    upstream_contract = flavor_axis.get("family_upstream_contract", {}) or {}
    if upstream_contract:
        lines.append(f"- **effective_runtime_ph:** {float(upstream_contract.get('effective_pH', 0.0)):.2f}" if upstream_contract.get("effective_pH") is not None else "- **effective_runtime_ph:** unknown")
        lines.append(f"- **dominant_donor_class:** {upstream_contract.get('dominant_donor_class', 'none')}")
        lines.append(f"- **donor_limited:** {upstream_contract.get('donor_limited', False)}")
        donor_pool_factors = upstream_contract.get("donor_pool_factors", {}) or {}
        if donor_pool_factors:
            lines.append(
                "- **donor_pool_factors:** "
                + ", ".join(f"{name}={float(value):.2f}" for name, value in donor_pool_factors.items())
            )
        if upstream_contract.get("pretreatment_interventions"):
            lines.append(
                "- **pretreatment_interventions:** "
                + ", ".join(str(item) for item in upstream_contract.get("pretreatment_interventions", []))
            )
        if upstream_contract.get("added_precursors"):
            lines.append(
                "- **upstream_added_precursors:** "
                + ", ".join(str(item) for item in upstream_contract.get("added_precursors", []))
            )
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
    if active_family_lanes:
        lines.append("- **active_family_lanes:** " + ", ".join(str(row.get("display_name", row.get("slr_family", "unknown"))) for row in active_family_lanes))
        for row in active_family_lanes:
            lines.append(
                f"- **family_lane_{row.get('slr_family', 'unknown')}:** {row.get('display_name', 'unknown')} | posture={row.get('strategic_posture', 'unknown')} | {row.get('summary', '')}"
            )
    lane_adjustments = flavor_axis.get("family_lane_adjustments", {})
    if lane_adjustments:
        lines.append(f"- **family_target_score_delta:** {float(lane_adjustments.get('target_score_delta', 0.0)):.2f}")
        lines.append(f"- **family_maillard_closure_delta:** {float(lane_adjustments.get('maillard_closure_delta', 0.0)):.2f}")
        lines.append(f"- **family_off_flavour_risk_delta:** {float(lane_adjustments.get('off_flavour_risk_delta', 0.0)):.2f}")
    family_prior_bundle = flavor_axis.get("family_prior_bundle", {})
    if isinstance(family_prior_bundle, dict) and family_prior_bundle:
        lines.append(
            "- **family_prior_bundle:** "
            + ", ".join(f"{family}={len(rows)}" for family, rows in sorted(family_prior_bundle.items()))
        )
    lipid_lane = family_lane_summary.get("02", {})
    if lipid_lane:
        lines.append(
            f"- **lipid_benchmark_ready_targets:** {', '.join(str(item) for item in lipid_lane.get('benchmark_ready_targets', [])) or 'none'}"
        )
        lines.append(
            f"- **lipid_crosstalk_priors:** {', '.join(str(item) for item in lipid_lane.get('competition_prior_ids', [])) or 'none'}"
        )
        lines.append(
            f"- **lipid_maillard_closure_pressure:** {float(lipid_lane.get('maillard_closure_pressure', 0.0)):.2f}"
        )
    for marker in flavor_axis.get("family_state_markers", []):
        lines.append(
            f"- **state_marker_{marker.get('marker_id', 'unknown')}:** {marker.get('display_name', 'unknown')} | role={marker.get('panel_role', 'unknown')} | kind={marker.get('observable_kind', 'unknown')} | influence={marker.get('influence_mode', 'unknown')} | summary={marker.get('state_value_summary', marker.get('summary', ''))}"
        )
    lines.append("")
    return "\n".join(lines)


def render_family_role_explanation_markdown(explanation: Dict[str, Any], *, heading: str = "## Family Role Explanation") -> str:
    """
    Renders the output of build_family_role_explanation() as structured markdown.

    Answers three questions a scientist should be able to ask:
      - Which families DROVE the result (score-impacting + benchmark-anchored)?
      - Which families were only MODIFIERS (active but score-neutral)?
      - Which families are MISSING or TRANSFERRED (not active this run)?
    """
    summary = explanation.get("summary", {})
    lines = [
        heading,
        "",
        f"Active family lanes: **{summary.get('active_lane_count', 0)}** of "
        f"**{summary.get('total_canonical_family_count', 0)}** canonical families",
        f"Drivers: {summary.get('driver_count', 0)} | "
        f"Modifiers: {summary.get('modifier_count', 0)} | "
        f"Missing / transferred: {summary.get('missing_or_transferred_count', 0)}",
        "",
    ]

    drivers = explanation.get("drivers", [])
    if drivers:
        lines += [
            "### Drivers — families that shaped this result",
            "",
            "| SLR | Family | Posture | Δ target | Δ off-flavour | Δ closure | Payloads | Description |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
        for row in drivers:
            lines.append(
                f"| {row['slr_family']} | {row['display_name']} | {row['strategic_posture']} "
                f"| {row['target_score_delta']:+.4f} | {row['off_flavour_risk_delta']:+.4f} "
                f"| {row['maillard_closure_delta']:+.4f} | {row['primary_payload_count']} "
                f"| {row.get('summary', '')} |"
            )
        lines.append("")

    modifiers = explanation.get("modifiers", [])
    if modifiers:
        lines += [
            "### Modifiers — active but score-neutral lanes",
            "",
            "| SLR | Family | Posture | Description |",
            "| --- | --- | --- | --- |",
        ]
        for row in modifiers:
            lines.append(
                f"| {row['slr_family']} | {row['display_name']} | {row['strategic_posture']} "
                f"| {row.get('summary', '')} |"
            )
        lines.append("")

    missing = explanation.get("missing_or_transferred", [])
    if missing:
        lines += [
            "### Missing / transferred — not active in this evaluation",
            "",
            "| SLR | Family | Payloads | Reason |",
            "| --- | --- | ---: | --- |",
        ]
        for row in missing:
            lines.append(
                f"| {row['slr_family']} | {row['display_name']} "
                f"| {row['primary_payload_count']} | {row.get('reason', '')} |"
            )
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


