from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.benchmark_validation import BenchmarkSummary, summarize_benchmarks
from src.inverse_design import FormulationResult, InverseDesigner
from src.smirks_engine import ReactionConditions
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT


@dataclass(frozen=True)
class ValidatedEnvelopeReport:
    target_tag: str
    strict_ready_benchmarks: List[str]
    supported_benchmarks: int
    total_benchmarks: int
    matrix_only_benchmarks: List[str]
    warnings: List[str]
    next_priorities: List[str]


def build_validated_envelope_report(target_tag: str = "meaty") -> ValidatedEnvelopeReport:
    summaries = summarize_benchmarks(target_tag=target_tag)
    strict_ready = [summary.benchmark_id for summary in summaries if summary.strict_ready]
    matrix_only = [summary.benchmark_id for summary in summaries if summary.execution_path == "matrix_only"]
    warnings = []
    if matrix_only:
        warnings.append(
            "Matrix benchmarks are executable intake/headspace checks, but remain outside the strict gate and target snapshots."
        )
    warnings.append(
        "Benchmark-facing concentrations still use the FAST observable projection; Cantera remains diagnostic-reference-only."
    )
    warnings.append(
        "Peptide-bound and intact-protein accessibility remain outside the validated precursor envelope."
    )
    warnings.append(
        "The validated plant-matrix scope is currently limited to pea/soy matrix-only systems and not yet a broad release gate."
    )
    next_priorities = [
        "Expose matrix-state and projection explainability in user-facing artifacts.",
        "Promote domain-of-applicability warnings into routine CLI/report outputs.",
        "Replace bulk matrix retention with broader compound-aware observability across plant-matrix systems.",
    ]
    return ValidatedEnvelopeReport(
        target_tag=target_tag,
        strict_ready_benchmarks=strict_ready,
        supported_benchmarks=sum(1 for summary in summaries if summary.supported),
        total_benchmarks=len(summaries),
        matrix_only_benchmarks=matrix_only,
        warnings=warnings,
        next_priorities=next_priorities,
    )


def render_validated_envelope_markdown(report: ValidatedEnvelopeReport) -> str:
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


def build_formulation_explainability_payload(
    formulation: Dict[str, object],
    result: FormulationResult,
    *,
    target_tag: str,
    minimize_tag: str,
) -> Dict[str, object]:
    top_rows = []
    seen = set()
    for name, observable in sorted(result.predicted_ppb.items(), key=lambda item: item[1], reverse=True):
        if not isinstance(name, str) or name in seen:
            continue
        if name in result.projection_metadata:
            metadata = result.projection_metadata[name]
        else:
            continue
        seen.add(name)
        top_rows.append({
            "compound": name,
            "proxy_ppb": metadata.get("proxy_ppb", observable),
            "observable_ppb": metadata.get("observable_ppb", observable),
            "observable_ratio": metadata.get("proxy_to_observable_ratio", 1.0),
            "matrix_factor": metadata.get("matrix_factor", 1.0),
            "headspace_factor": metadata.get("headspace_factor", 1.0),
            "volatile_class": metadata.get("volatile_class", "other"),
        })
        if len(top_rows) >= 8:
            break

    return {
        "formulation_name": formulation.get("name", "unknown"),
        "target_tag": target_tag,
        "minimize_tag": minimize_tag,
        "protein_type": formulation.get("protein_type", "free"),
        "matrix_explainability": result.matrix_explainability,
        "scores": {
            "target_score": result.target_score,
            "off_flavour_risk": result.off_flavour_risk,
            "safety_score": result.safety_score,
            "texture_risk": result.texture_risk,
            "lysine_budget": result.lysine_budget,
            "trapping_efficiency": result.trapping_efficiency,
        },
        "top_projection_rows": top_rows,
        "detected_targets": result.detected_targets,
        "detected_minimize": result.detected_minimize,
    }


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
    lines.extend([
        "",
        "## Projection Rows",
        "",
        "| Compound | Proxy ppb | Observable ppb | Obs/Proxy | Matrix | Headspace | Class |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload["top_projection_rows"]:
        lines.append(
            f"| {row['compound']} | {row['proxy_ppb']:.3f} | {row['observable_ppb']:.3f} | {row['observable_ratio']:.3f} | {row['matrix_factor']:.3f} | {row['headspace_factor']:.3f} | {row['volatile_class']} |"
        )
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class DomainWarning:
    level: str  # "WARNING" or "CAUTION"
    category: str
    message: str


class DomainOfValidityChecker:
    """
    Checks if a formulation and conditions are within the "trusted envelope" 
    of validated benchmarks and physical assumptions.
    """
    
    TRUSTED_PRECURSORS = {
        "cysteine", "l-cysteine", "ribose", "d-ribose", 
        "glucose", "d-glucose", "asparagine", "l-asparagine"
    }
    
    MAX_TRUSTED_TEMP_C = 180.0
    MIN_TRUSTED_PH = 3.0
    MAX_TRUSTED_PH = 9.0

    def __init__(self, target_tag: str = "meaty"):
        self.target_tag = target_tag

    def check(
        self, 
        precursor_names: List[str], 
        protein_type: str, 
        temp_c: float, 
        ph: float,
        aw: float = 1.0
    ) -> List[DomainWarning]:
        warnings = []
        
        # 1. Matrix Check
        if protein_type != "free":
            warnings.append(DomainWarning(
                level="CAUTION",
                category="MATRIX",
                message=f"Matrix '{protein_type}' uses speculative accessibility scaling; PRIMARY benchmarks are free-precursor only."
            ))

        # 2. Precursor Support Check
        unsupported_precursors = []
        for name in precursor_names:
            norm_name = name.strip().lower()
            if norm_name not in self.TRUSTED_PRECURSORS:
                unsupported_precursors.append(name)
        
        if unsupported_precursors:
            warnings.append(DomainWarning(
                level="WARNING",
                category="PRECURSORS",
                message=f"Sparse benchmark analogies: {', '.join(unsupported_precursors)} lack PRIMARY quantitative validation."
            ))

        # 3. Process Condition Check
        if temp_c > self.MAX_TRUSTED_TEMP_C:
            warnings.append(DomainWarning(
                level="CAUTION",
                category="SEVERITY",
                message=f"Temperature {temp_c}°C exceeds the validated Arrhenius envelope (> {self.MAX_TRUSTED_TEMP_C}°C)."
            ))
            
        if ph < self.MIN_TRUSTED_PH or ph > self.MAX_TRUSTED_PH:
            warnings.append(DomainWarning(
                level="WARNING",
                category="PHYSICOCHEMICAL",
                message=f"pH {ph} is outside the trusted kinetic range ({self.MIN_TRUSTED_PH}-{self.MAX_TRUSTED_PH})."
            ))

        return warnings


def render_domain_warnings_markdown(warnings: List[DomainWarning]) -> str:
    if not warnings:
        return "> [!NOTE]\n> Run is within the validated scientific envelope.\n"
        
    lines = ["## Scientific Domain Warnings", ""]
    for w in warnings:
        icon = "⚠️" if w.level == "WARNING" else "🚨"
        lines.append(f"- **{icon} {w.level} [{w.category}]**: {w.message}")
    return "\n".join(lines) + "\n"


def render_domain_warnings_cli(warnings: List[DomainWarning]):
    if not warnings:
        print("\n  ✅ Scientific Envelope: Trusted (matches PRIMARY benchmarks)")
        return

    print("\n  ⚠️  SCIENTIFIC DOMAIN WARNINGS:")
    for w in warnings:
        icon = "!" if w.level == "WARNING" else "!!"
        print(f"    [{icon}] {w.category.ljust(15)}: {w.message}")
    print("  " + "─" * 60)


def render_decision_summary_cli(result: "FormulationResult", warnings: List[DomainWarning]):
    """
    Renders a premium, consolidated scientist-facing summary of the simulation run.
    """
    print("\n" + "═" * 80)
    print(" " * 25 + "📊 MAILLARD DECISION SUMMARY")
    print("═" * 80)

    # 1. Scientific Envelope Section
    status_icon = "✅" if not warnings else "⚠️"
    status_label = "TRUSTED" if not warnings else "LIMITED"
    print(f"\n  [1] SCIENTIFIC ENVELOPE: {status_label} {status_icon}")
    if warnings:
        for w in warnings:
            dot = "!" if w.level == "WARNING" else "!!"
            print(f"      {dot} [{w.category}] {w.message}")
    else:
        print("      Run matches PRIMARY quantitative benchmarks (Cys/Glc/Rib/Asn).")

    # 2. Flavor Profile Section (Desirable & Penalties)
    print("\n  [2] SENSORY PROJECTION:")
    
    # Sort targets by ppb descending
    targets = []
    for canon, ppb in result.predicted_ppb.items():
        # Heuristic to find the label if we have projection_metadata
        meta = result.projection_metadata.get(canon, {})
        # This is a bit tricky as FormulationResult is defined elsewhere.
        # But we can assume common keys.
        targets.append((canon, ppb))
    
    # Actually, the radar data from FormulationResult is probably better
    # radar: Dict[str, tuple] (score, count)
    radar_desirable = sorted(
        [(tag, score) for tag, (score, cnt) in result.radar.items() if score > 0 and tag != "beany"],
        key=lambda x: x[1],
        reverse=True
    )
    
    top_flavor = radar_desirable[0][0].upper() if radar_desirable else "NEUTRAL"
    print(f"      Dominant Profile : {top_flavor}")
    
    if radar_desirable[:3]:
        notes = ", ".join([f"{t[0]}({t[1]:.1f})" for t in radar_desirable[:3]])
        print(f"      Key Targets      : {notes}")
    
    beany_score = result.radar.get("beany", (0.0, 0))[0]
    if beany_score > 0:
        print(f"      Off-Flavor Risk  : BEANY ({beany_score:.1f})")

    # 3. Safety & Compliance
    safety_icon = "✅" if result.safety_score == 0 else "🚨"
    print(f"\n  [3] SAFETY & COMPLIANCE: {safety_icon}")
    if result.flagged_toxics:
        print(f"      WARNING: Detected {', '.join(result.flagged_toxics)}")
    else:
        print(f"      No regulated toxic markers detected at process thresholds.")

    # 4. Matrix & Physicality
    print(f"\n  [4] MATRIX & PHYSICALITY:")
    p_type = result.matrix_explainability.get("protein_type", "unknown")
    denat = result.matrix_explainability.get("effective_denaturation_state", 0.0)
    print(f"      Protein Matrix   : {p_type} (Denaturation: {denat:.1%})")
    
    # Lysine budget if available
    if result.lysine_budget > 0:
        print(f"      Lysine Budget    : {result.lysine_budget:.1f}% consumed by DHA pathways")

    print("\n" + "═" * 80 + "\n")


def render_deep_explainability_cli(result: "FormulationResult"):
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
    # We estimate this from projection_metadata of the target profile
    target_metadata = []
    for canon, meta in result.projection_metadata.items():
        # Check if it's a target (using ppb as proxy for activity)
        if result.predicted_ppb.get(canon, 0) > 0:
            target_metadata.append(meta)
            
    avg_loss = 0.0
    if target_metadata:
        avg_loss = 100.0 * (1.0 - sum(m.get("matrix_factor", 1.0) for m in target_metadata) / len(target_metadata))

    print(f"      Avg. Matrix Loss   : {avg_loss:.1f}% potential flavor trapped")

    # Panel B: Performance Drivers
    print(f"\n  [B] PERFORMANCE DRIVERS:")
    # Top 2 targets
    top_targets = sorted(
        [(m.get("compound", "unknown"), m.get("observable_ppb", 0.0)) for m in result.projection_metadata.values()],
        key=lambda x: x[1],
        reverse=True
    )[:2]
    
    if top_targets:
        targets_str = ", ".join([f"{t[0]} ({t[1]:.1f} ppb)" for t in top_targets])
        print(f"      Top Yield Drivers  : {targets_str}")

    # Off-flavor drivers if any
    beany_drivers = []
    for name, ppb in result.predicted_ppb.items():
        if "hexanal" in name.lower() or "nonanal" in name.lower() or "decadienal" in name.lower():
             beany_drivers.append((name, ppb))
    
    if beany_drivers:
        beany_drivers.sort(key=lambda x: x[1], reverse=True)
        off_str = ", ".join([f"{d[0]} ({d[1]:.1f} ppb)" for d in beany_drivers[:2]])
        print(f"      Off-Flavor Drivers : {off_str}")

    print(f"\n  [💡] INTERVENTION HINT:")
    hint = generate_intervention_hint(result)
    print(f"      {hint}")


def generate_intervention_hint(result: "FormulationResult") -> str:
    """
    Generates an actionable intervention strategy based on simulation results.
    """
    # 1. Yield Bottleneck
    sev = getattr(result, 'bottleneck_severity', 0.0)
    b_prec = getattr(result, 'bottleneck_precursor', 'none')
    if sev > 0.5:
        return f"Yield limited by {b_prec}. Increase its molar ratio to boost flavor."
        
    # 2. Matrix Loss
    target_metadata = [m for m in result.projection_metadata.values() if m.get("observable_ppb", 0) > 0]
    avg_matrix = sum(m.get("matrix_factor", 1.0) for m in target_metadata) / len(target_metadata) if target_metadata else 1.0
    if avg_matrix < 0.5:
        return "High matrix retention. Increase heating time or pH to recover trapped volatiles."
        
    # 3. Off-flavors
    if result.off_flavour_risk > 15.0:
        return "High beany risk. Check lipid precursor purity or add oxidation inhibitors."
        
    # 4. Success
    if result.target_score > 30.0:
        return "Formulation balanced. Consider sensory panel validation."
        
    return "Check precursors and process conditions for missing flavor potential."