from __future__ import annotations

import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, TYPE_CHECKING
import subprocess

from src.pipeline import FormulationResult
from src.projection_metadata import ProjectionMetadataRow


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
        name = target.get("name")
        meta: Optional['ProjectionMetadataRow'] = metadata.get(name)
        if not meta:
            continue
            
        rows.append({
            "compound": name,
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
            "browning_index": target.get("browning_index", 0.0),
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
    from src.reporting import _repo_root, _safe_git_output, _build_scientific_surface, SCHEMA_VERSION
    import json
    import hashlib
    import sys
    import platform
    import shlex

    root = _repo_root()
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
