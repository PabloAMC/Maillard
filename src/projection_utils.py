from __future__ import annotations

import datetime
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, TYPE_CHECKING
import subprocess

from src.pipeline import FormulationResult
from src.projection_metadata import ProjectionMetadataRow, normalize_projection_metadata_row


def _support_origin_from_projection_meta(meta: Mapping[str, object]) -> str:
    process_state = str(meta.get("process_state", "unknown"))
    calibration_state = str(meta.get("calibration_process_state", "unknown"))
    fallback = str(meta.get("calibration_fallback_mode", "class_level"))
    source = str(meta.get("calibration_source", "class_fallback")).lower()
    strength = str(meta.get("calibration_evidence_strength", "heuristic")).lower()
    extrusion_state = process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"}

    if extrusion_state:
        if fallback == "nearest_process_state" or calibration_state not in {"unknown", process_state}:
            return "lower_regime_transfer"
        if strength == "heuristic" or source == "class_fallback":
            return "extrusion_extrapolation"
        return "extrusion_specific_support"

    return "standard_matrix_support"


def _reachability_payload(meta: Mapping[str, object]) -> Dict[str, object]:
    evidence_state = str(meta.get("evidence_state", "still_missing")).lower()
    source = str(meta.get("calibration_source", "class_fallback")).lower()
    strength = str(meta.get("calibration_evidence_strength", "heuristic")).lower()
    fallback = str(meta.get("calibration_fallback_mode", "class_level")).lower()
    notes_blob = " ".join(
        [
            source,
            str(meta.get("calibration_notes", "")).lower(),
            str(meta.get("retention_runtime_mode", "")).lower(),
        ]
    )

    direct_anchor = evidence_state in {"externally_benchmarked", "internally_benchmarked", "conditional_calibration"} or (
        strength == "literature_anchored" and fallback == "compound_specific"
    )
    transferred_support = evidence_state == "transferred_prior" or (
        strength in {"conditional_literature_anchored", "process_state_mismatch", "directional_transferred", "class_anchored"}
        or fallback in {"nearest_process_state", "compound_specific_process_state", "class_level"}
        or "transfer" in source
        or "carryover" in source
        or "ratio" in source
    )
    computational_refinement = any(token in notes_blob for token in ["dft", "xtb", "qm", "semiempirical", "computational", "refinement"])

    if direct_anchor:
        return {
            "chemically_reachable": True,
            "reachability_status": "chemically_reachable",
            "reachability_basis": "direct_anchor",
        }
    if transferred_support or computational_refinement:
        return {
            "chemically_reachable": True,
            "reachability_status": "conditionally_reachable",
            "reachability_basis": "transferred_or_refined_support",
        }
    return {
        "chemically_reachable": False,
        "reachability_status": "merely_plausible",
        "reachability_basis": "mechanistic_surrogate_only",
    }


def _observable_assumption_summary(meta: Mapping[str, object]) -> str:
    retention_mode = str(meta.get("retention_runtime_mode", "static_class_profile"))
    fallback_mode = str(meta.get("calibration_fallback_mode", "class_level"))
    support_origin = _support_origin_from_projection_meta(meta)
    accessibility_warning = bool(meta.get("accessibility_warning", False))
    parts = [retention_mode, fallback_mode, support_origin]
    if accessibility_warning:
        parts.append("accessibility_dominated")
    return " | ".join(parts)


def _panel_modeling_regimes(meta: Mapping[str, object]) -> List[str]:
    values = meta.get("modeling_regimes", [])
    if not isinstance(values, list):
        return []
    return [str(item) for item in values if str(item)]


def _resolve_projection_metadata_row(
    target: Mapping[str, object],
    metadata: Mapping[str, ProjectionMetadataRow],
) -> Optional[ProjectionMetadataRow]:
    target_name = str(target.get("name", "")).strip()
    observable_fallback = float(target.get("concentration", 0.0) or 0.0)

    target_projection = target.get("projection")
    if isinstance(target_projection, Mapping) and target_projection:
        return normalize_projection_metadata_row(
            target_projection,
            compound_fallback=target_name or "unknown",
            observable_ppb_fallback=observable_fallback,
        )

    direct_match = metadata.get(target_name)
    if direct_match:
        return normalize_projection_metadata_row(
            direct_match,
            compound_fallback=target_name or "unknown",
            observable_ppb_fallback=observable_fallback,
        )

    lowered_target_name = target_name.lower()
    for raw_row in metadata.values():
        compound_name = str(raw_row.get("compound", "")).strip()
        if compound_name.lower() == lowered_target_name:
            return normalize_projection_metadata_row(
                raw_row,
                compound_fallback=target_name or compound_name or "unknown",
                observable_ppb_fallback=observable_fallback,
            )
    return None


def build_projection_rows(
    result: 'FormulationResult',
) -> List[Dict[str, object]]:
    """
    Builds a list of projection metadata rows for a given recommendation result.
    This is the core data-building logic for reporting.
    """
    rows: List[Dict[str, object]] = []
    metadata = result.projection_metadata
    
    # Sort by concentration descending
    targets = sorted(
        result.targets,
        key=lambda t: float(t.get("concentration", 0.0)),
        reverse=True
    )

    for target in targets:
        name = str(target.get("name", "unknown"))
        meta = _resolve_projection_metadata_row(target, metadata)
        if not meta:
            continue
            
        rows.append({
            "compound": str(meta.get("compound", name)),
            "proxy_ppb": meta.get("proxy_ppb"),
            "observable_ppb": meta.get("observable_ppb"),
            "observable_ratio": meta.get("observable_ratio", meta.get("proxy_to_observable_ratio")),
            "matrix_factor": meta.get("matrix_factor"),
            "dynamic_retention_factor": meta.get("dynamic_retention_factor"),
            "headspace_factor": meta.get("headspace_factor"),
            "volatile_class": meta.get("volatile_class"),
            "process_state": meta.get("process_state"),
            "retention_runtime_mode": meta.get("retention_runtime_mode"),
            "calibration_source": meta.get("calibration_source"),
            "calibration_evidence_strength": meta.get("calibration_evidence_strength"),
            "calibration_fallback_mode": meta.get("calibration_fallback_mode"),
            "evidence_state": meta.get("evidence_state"),
            "target_class": meta.get("target_class"),
            "panel_role": meta.get("panel_role"),
            "observable_kind": meta.get("observable_kind"),
            "modeling_regimes": _panel_modeling_regimes(meta),
            "chemistry_family": meta.get("chemistry_family"),
            "supporting_families": list(meta.get("supporting_families", [])),
            "observable_panel_tags": list(meta.get("observable_panel_tags", [])),
            "browning_index": meta.get("browning_index", target.get("browning_index", 0.0)),
            "support_origin": _support_origin_from_projection_meta(meta),
            "observable_assumption_summary": _observable_assumption_summary(meta),
            **_reachability_payload(meta),
        })

    if rows:
        return rows

    for meta in sorted(
        metadata.values(),
        key=lambda row: float(row.get("observable_ppb", 0.0) or 0.0),
        reverse=True,
    ):
        normalized = normalize_projection_metadata_row(
            meta,
            compound_fallback=str(meta.get("compound", "unknown")),
            observable_ppb_fallback=float(meta.get("observable_ppb", 0.0) or 0.0),
        )
        rows.append({
            "compound": str(normalized.get("compound", "unknown")),
            "proxy_ppb": normalized.get("proxy_ppb"),
            "observable_ppb": normalized.get("observable_ppb"),
            "observable_ratio": normalized.get("observable_ratio", normalized.get("proxy_to_observable_ratio")),
            "matrix_factor": normalized.get("matrix_factor"),
            "dynamic_retention_factor": normalized.get("dynamic_retention_factor"),
            "headspace_factor": normalized.get("headspace_factor"),
            "volatile_class": normalized.get("volatile_class"),
            "process_state": normalized.get("process_state"),
            "retention_runtime_mode": normalized.get("retention_runtime_mode"),
            "calibration_source": normalized.get("calibration_source"),
            "calibration_evidence_strength": normalized.get("calibration_evidence_strength"),
            "calibration_fallback_mode": normalized.get("calibration_fallback_mode"),
            "evidence_state": normalized.get("evidence_state"),
            "target_class": normalized.get("target_class"),
            "panel_role": normalized.get("panel_role"),
            "observable_kind": normalized.get("observable_kind"),
            "modeling_regimes": _panel_modeling_regimes(normalized),
            "chemistry_family": normalized.get("chemistry_family"),
            "supporting_families": list(normalized.get("supporting_families", [])),
            "observable_panel_tags": list(normalized.get("observable_panel_tags", [])),
            "browning_index": normalized.get("browning_index", 0.0),
            "support_origin": _support_origin_from_projection_meta(normalized),
            "observable_assumption_summary": _observable_assumption_summary(normalized),
            **_reachability_payload(normalized),
        })
    return rows


def build_artifact_provenance(
    artifact_kind: str,
    output_dir: Path,
    inputs: Any,
    campaign_metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Builds robust provenance metadata for an artifact, including git state and scientific surface.
    """
    from src.reporting import _safe_git_output, _build_scientific_surface, SCHEMA_VERSION
    from src.artifact_io import repo_root
    import json
    import hashlib
    import sys
    import platform
    import shlex

    root = repo_root()
    status_text = _safe_git_output(root, ["status", "--porcelain"]) or ""
    changed_paths = []
    for line in status_text.splitlines():
        path_text = line[3:].strip()
        if path_text:
            changed_paths.append(path_text)

    serialized_inputs = json.dumps(inputs, sort_keys=True, default=str)
    provenance: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "artifact_kind": artifact_kind,
        "generated_at": datetime.datetime.now().isoformat(),
        "generator": {
            "entrypoint": Path(sys.argv[0]).name,
            "argv": sys.argv[1:],
            "command": shlex.join(sys.argv),
        },
        "runtime": {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "repository": {
            "name": root.name,
            "root": str(root),
            "branch": _safe_git_output(root, ["rev-parse", "--abbrev-ref", "HEAD"]) or "unknown",
            "commit": _safe_git_output(root, ["rev-parse", "HEAD"]) or "unknown",
            "short_commit": _safe_git_output(root, ["rev-parse", "--short", "HEAD"]) or "unknown",
            "dirty": bool(changed_paths),
            "changed_file_count": len(changed_paths),
            "changed_files_sample": changed_paths[:10],
        },
        "input_fingerprint_sha256": hashlib.sha256(serialized_inputs.encode("utf-8")).hexdigest(),
        "output_directory": str(output_dir.relative_to(root) if output_dir.is_relative_to(root) else output_dir),
        "scientific_surface": _build_scientific_surface(root),
    }
    if campaign_metadata:
        provenance["campaign"] = campaign_metadata
    return provenance


def generate_intervention_hint(result: 'FormulationResult') -> str:
    """
    Generates an actionable intervention strategy based on simulation results.
    """
    # 1. Yield Bottleneck
    sev = getattr(result, 'bottleneck_severity', 0.0)
    b_prec = getattr(result, 'bottleneck_precursor', 'none')
    if sev > 0.6:
        return f"Simulation shows severe {b_prec} depletion. Increase its ratio to boost aromatics."
        
    # 2. Physical Suppression
    suppressed = getattr(result, 'suppressed_compounds', [])
    if suppressed:
        top = suppressed[0]
        if top['reduction_factor'] > 0.8:
            if top['primary_cause'] == 'headspace':
                return f"Major headspace loss for {top['name']}. Optimize processing (time/temp) to release volatiles."
            else:
                return f"{top['name']} is heavily trapped in the matrix. Try increasing pH or denaturation state."

    # 3. Precursor Attribution / Balance
    attrib = getattr(result, 'precursor_contributions', {})
    if attrib:
        sorted_attrib = sorted(attrib.items(), key=lambda x: x[1], reverse=True)
        if len(sorted_attrib) > 1:
            top_name, top_val = sorted_attrib[0]
            next_name, next_val = sorted_attrib[1]
            if top_val > 5.0 * next_val:
                return f"Yield is overwhelmingly driven by {top_name}. Add {next_name} to diversify the profile."

    # 4. Off-flavors
    if result.off_flavour_risk > 15.0:
        return "High lipid-beany risk. Reduce lipid precursors or add specific antioxidants."
        
    # 5. Success
    if result.target_score > 30.0:
        return "High-performance formulation. Proceed to sensory validation."
        
    return "Balanced profile but low intensity. Consider increasing temperature or total precursor load."
