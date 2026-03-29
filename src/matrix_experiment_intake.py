from __future__ import annotations

import json
import datetime as dt
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import yaml

from src.benchmark_validation import (
    DEFAULT_TARGET_TAG,
    _matrix_compound_support_status,
    _projection_metadata_for_match,
    assess_matrix_benchmark_evidence,
    build_matrix_target_status_artifact,
    evaluate_benchmark_payload,
    get_benchmark_files,
    get_matrix_ranking_contract,
    get_benchmark_metadata,
    load_benchmark,
    summarize_evaluation_for_benchmark,
)
from src.matrix_targets import get_compound_panel_entry


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEMA_PATH = ROOT / "data" / "protocols" / "matrix_experiment_intake_schema.json"
DEFAULT_EXAMPLE_PATH = ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml"

_SOURCE_KINDS = {"external_literature", "internal_experiment", "synthetic_diagnostic"}
_PROTEIN_TYPES = {"pea_iso", "soy_iso", "myco", "free"}
_SUPPORT_RANK = {
    "open_gap": 0,
    "directional_support": 1,
    "internal_candidate": 2,
    "quantitative_closed": 3,
}


def load_matrix_experiment_intake_schema(path: Optional[Path | str] = None) -> Dict[str, Any]:
    schema_path = Path(path) if path is not None else DEFAULT_SCHEMA_PATH
    with open(schema_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_matrix_experiment_intake(path: Path | str) -> Dict[str, Any]:
    payload_path = Path(path)
    text = payload_path.read_text(encoding="utf-8")
    if payload_path.suffix.lower() in {".yaml", ".yml"}:
        payload = yaml.safe_load(text)
    else:
        payload = json.loads(text)
    if not isinstance(payload, dict):
        raise ValueError("Matrix experiment intake payload must be a mapping")
    return payload


def normalize_matrix_experiment_intake(payload: Mapping[str, Any]) -> Dict[str, Any]:
    schema = load_matrix_experiment_intake_schema()
    normalized = dict(payload)

    for field in schema.get("required_top_level_fields", []):
        if field not in normalized:
            raise ValueError(f"Missing required top-level field: {field}")

    source_kind = str(normalized.get("source_kind", "")).strip()
    if source_kind not in _SOURCE_KINDS:
        raise ValueError(f"Unsupported source_kind: {source_kind}")
    protein_type = str(normalized.get("protein_type", "")).strip()
    if protein_type not in _PROTEIN_TYPES:
        raise ValueError(f"Unsupported protein_type: {protein_type}")

    conditions = dict(normalized.get("conditions") or {})
    for field in schema.get("required_conditions", []):
        if field not in conditions:
            raise ValueError(f"Missing required condition field: {field}")
        conditions[field] = float(conditions[field])
    normalized["conditions"] = conditions

    formulation = dict(normalized.get("formulation") or {})
    precursors = dict(formulation.get("precursors") or {})
    if not precursors:
        raise ValueError("Formulation must include at least one precursor entry")
    for name, row in list(precursors.items()):
        entry = dict(row or {})
        if "concentration_mM" not in entry:
            raise ValueError(f"Precursor {name} is missing concentration_mM")
        entry["concentration_mM"] = float(entry["concentration_mM"])
        precursors[name] = entry
    formulation["precursors"] = precursors
    normalized["formulation"] = formulation

    measured = dict(normalized.get("measured_volatiles") or {})
    if not measured:
        raise ValueError("Payload must include at least one measured volatile")
    for name, row in list(measured.items()):
        entry = dict(row or {})
        if "conc_ppb" not in entry:
            raise ValueError(f"Measured volatile {name} is missing conc_ppb")
        entry["conc_ppb"] = float(entry["conc_ppb"])
        if "uncertainty_pct" in entry and entry["uncertainty_pct"] is not None:
            entry["uncertainty_pct"] = float(entry["uncertainty_pct"])
        measured[name] = entry
    normalized["measured_volatiles"] = measured

    provenance = dict(normalized.get("provenance") or {})
    for field in schema.get("required_provenance_fields", []):
        if field not in provenance:
            raise ValueError(f"Missing provenance field: {field}")
    for key, value in list(provenance.items()):
        if isinstance(value, (dt.date, dt.datetime)):
            provenance[key] = value.isoformat()
    normalized["provenance"] = provenance

    normalized["experiment_id"] = str(normalized.get("experiment_id", "matrix_experiment_payload")).strip()
    normalized["source_kind"] = source_kind
    normalized["protein_type"] = protein_type
    normalized["process_state"] = str(normalized.get("process_state", "unknown")).strip()
    normalized["matrix_format"] = str(normalized.get("matrix_format", "unspecified")).strip()
    normalized["benchmark_alignment"] = dict(normalized.get("benchmark_alignment") or {})
    normalized["comparison_contract"] = dict(normalized.get("comparison_contract") or {})
    normalized["analytical_context"] = dict(normalized.get("analytical_context") or {})
    normalized["denaturation_state"] = float(normalized.get("denaturation_state", 0.5))
    return normalized


def _infer_role_from_panel(compound: str) -> str:
    panel_entry = get_compound_panel_entry(compound) or {}
    target_class = str(panel_entry.get("target_class", "unknown"))
    if target_class in {"adverse_lipid_markers", "severity_markers", "safety_markers"}:
        return "adverse_marker"
    return "desirable_marker"


def _default_comparison_contract(payload: Mapping[str, Any]) -> Dict[str, Any]:
    measured = payload.get("measured_volatiles", {})
    sorted_targets = sorted(
        measured.items(),
        key=lambda item: float((item[1] or {}).get("conc_ppb", 0.0)),
        reverse=True,
    )
    observable_targets = []
    adverse_markers = []
    for index, (compound, _row) in enumerate(sorted_targets, start=1):
        role = _infer_role_from_panel(str(compound))
        if role == "adverse_marker":
            adverse_markers.append(str(compound))
        observable_targets.append(
            {
                "name": str(compound),
                "role": role,
                "expected_rank": index,
                "direction": "lower" if role == "adverse_marker" else "higher",
            }
        )
    return {
        "observable_targets": observable_targets,
        "adverse_markers": adverse_markers,
        "citation_provenance": [str(payload.get("provenance", {}).get("source_reference", "unspecified"))],
        "notes": "Auto-derived from the measured panel because no explicit comparison_contract was supplied.",
    }


def build_matrix_experiment_benchmark_payload(payload: Mapping[str, Any]) -> Dict[str, Any]:
    normalized = normalize_matrix_experiment_intake(payload)
    comparison_contract = normalized.get("comparison_contract") or _default_comparison_contract(normalized)
    precursors = normalized.get("formulation", {}).get("precursors", {})
    matrix_only = len(precursors) == 1
    execution_path = "matrix_only" if matrix_only else "matrix_precursor_augmented"
    family = "matrix_headspace" if matrix_only else "matrix_precursor_augmented"
    tier = "PRIMARY" if normalized.get("source_kind") == "external_literature" else "SECONDARY"
    provenance = normalized.get("provenance", {})

    return {
        "benchmark_id": normalized.get("experiment_id"),
        "source_doi": provenance.get("source_doi") if normalized.get("source_kind") == "external_literature" else None,
        "source_metadata": {
            "origin": provenance.get("origin", normalized.get("source_kind")),
            "measurement_date": provenance.get("measurement_date"),
            "generator": "matrix_experiment_intake",
            "notes": provenance.get("notes") or normalized.get("analytical_context", {}).get("notes"),
        },
        "precursors": precursors,
        "conditions": {
            "temp_C": normalized["conditions"]["temp_C"],
            "ph": normalized["conditions"]["ph"],
            "water_activity": normalized["conditions"]["water_activity"],
            "time_min": normalized["conditions"]["time_min"],
        },
        "metadata": {
            "tier": tier,
            "family": family,
            "execution_path": execution_path,
            "notes": "Generated from matrix experiment intake payload.",
        },
        "process_metadata": {
            "state": normalized.get("process_state"),
            "preheating": normalized.get("benchmark_alignment", {}).get("preheating", "unspecified"),
            "extrusion_history": normalized.get("benchmark_alignment", {}).get("extrusion_history", "unspecified"),
            "matrix_format": normalized.get("matrix_format", "unspecified"),
        },
        "matrix_ranking_contract": comparison_contract,
        "measured_volatiles": normalized.get("measured_volatiles", {}),
        "protein_type": normalized.get("protein_type"),
        "denaturation_state": normalized.get("denaturation_state", 0.5),
    }


def materialize_matrix_experiment_benchmark(
    payload_or_path: Mapping[str, Any] | Path | str,
    *,
    output_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    if isinstance(payload_or_path, (Path, str)):
        payload = load_matrix_experiment_intake(payload_or_path)
    else:
        payload = dict(payload_or_path)

    benchmark = build_matrix_experiment_benchmark_payload(payload)
    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(benchmark, indent=2), encoding="utf-8")
    return benchmark


def select_aligned_matrix_benchmark(
    payload: Mapping[str, Any],
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> Optional[Path]:
    normalized = normalize_matrix_experiment_intake(payload)
    target_benchmark_id = str(normalized.get("benchmark_alignment", {}).get("target_benchmark_id", "")).strip()
    bench_files = [Path(item) for item in (benchmark_files or get_benchmark_files())]
    if target_benchmark_id:
        for bench_file in bench_files:
            bench = load_benchmark(bench_file)
            if str(bench.get("benchmark_id", "")) == target_benchmark_id:
                return Path(bench_file)

    measured_names = {str(name).strip().lower() for name in normalized.get("measured_volatiles", {})}

    def _origin_priority(bench: Mapping[str, Any]) -> int:
        origin = str((bench.get("source_metadata") or {}).get("origin", "")).strip().lower()
        return {
            "external_literature": 4,
            "internal_reproducibility_candidate": 3,
            "internal_experiment": 2,
            "synthetic_diagnostic": 1,
        }.get(origin, 0)

    best_path: Optional[Path] = None
    best_score: tuple[int, int, int, int, int, str] = (-1, -1, -1, -1, -1, "")
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue
        same_protein = int(str(bench.get("protein_type", "")) == normalized.get("protein_type"))
        same_state = int(str((bench.get("process_metadata") or {}).get("state", "")) == normalized.get("process_state"))
        same_path = int(str(metadata.execution_path) == ("matrix_only" if len(normalized.get("formulation", {}).get("precursors", {})) == 1 else "matrix_precursor_augmented"))
        contract = get_matrix_ranking_contract(bench)
        overlap = len(measured_names.intersection({str(item.get("name", "")).strip().lower() for item in contract.get("observable_targets", [])}))
        score = (same_protein, same_state, same_path, overlap, _origin_priority(bench), str(bench.get("benchmark_id", "")))
        if score > best_score:
            best_score = score
            best_path = Path(bench_file)
    return best_path


def _build_support_rows_for_benchmark(
    bench: Mapping[str, Any],
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    evaluation = evaluate_benchmark_payload(dict(bench), benchmark_id=str(bench.get("benchmark_id", "matrix_experiment_payload")), target_tag=target_tag)
    summary = summarize_evaluation_for_benchmark(
        evaluation,
        dict(bench),
        protein_type=str(bench.get("protein_type", "free")),
    )
    evidence = assess_matrix_benchmark_evidence(dict(bench))
    contract = get_matrix_ranking_contract(dict(bench))
    adverse_markers = {str(item).strip().lower() for item in contract.get("adverse_markers", [])}
    compounds = []
    counts = {
        "quantitative_closed": 0,
        "internal_candidate": 0,
        "directional_support": 0,
        "open_gap": 0,
    }
    for item in contract.get("observable_targets", []):
        compound_name = str(item.get("name", "")).strip()
        if not compound_name:
            continue
        matched = next((comparison for comparison in evaluation.comparisons if comparison.compound == compound_name), None)
        if matched is None:
            continue
        meta = _projection_metadata_for_match(evaluation, matched)
        role = str(item.get("role", "adverse_marker" if compound_name.lower() in adverse_markers else "desirable_marker"))
        support_status = _matrix_compound_support_status(
            evidence_state=str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
            calibration_evidence_strength=str(meta.get("calibration_evidence_strength", "heuristic")),
            reference_signal_origin=summary.reference_signal_origin,
            source_origin=evidence.source_origin,
        )
        counts[support_status] += 1
        compounds.append(
            {
                "compound": compound_name,
                "role": role,
                "target_class": str(meta.get("target_class", item.get("target_class", "unknown"))),
                "evidence_state": str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
                "calibration_evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                "support_status": support_status,
                "measured_ppb": float(matched.measured_ppb),
                "predicted_ppb": float(matched.predicted_ppb),
                "ratio": float(matched.ratio),
            }
        )

    promotion_ready = (
        evidence.target_profile in {"mixed", "meaty_positive"}
        and summary.ranking_contract_status == "pass"
        and counts["quantitative_closed"] >= 2
        and counts["internal_candidate"] == 0
        and counts["directional_support"] == 0
    )
    if evidence.target_profile not in {"mixed", "meaty_positive"}:
        blocker = "benchmark lacks meaty-positive targets"
    elif summary.ranking_contract_status != "pass":
        blocker = "ranking contract not yet passing"
    elif counts["quantitative_closed"] < 2:
        blocker = "insufficient externally measured target closure"
    elif counts["internal_candidate"] > 0 or counts["directional_support"] > 0:
        blocker = "depends on internal or transferred support"
    else:
        blocker = "none"

    return {
        "summary": summary,
        "evidence": evidence,
        "support_counts": counts,
        "compounds": compounds,
        "promotion_ready": promotion_ready,
        "promotion_blocker": blocker,
    }


def build_matrix_experiment_support_delta_artifact(
    payload_or_path: Mapping[str, Any] | Path | str,
    *,
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    payload = load_matrix_experiment_intake(payload_or_path) if isinstance(payload_or_path, (Path, str)) else dict(payload_or_path)
    normalized = normalize_matrix_experiment_intake(payload)
    intake_benchmark = build_matrix_experiment_benchmark_payload(normalized)
    current = _build_support_rows_for_benchmark(intake_benchmark, target_tag=target_tag)

    aligned_path = select_aligned_matrix_benchmark(normalized, benchmark_files=benchmark_files)
    baseline_row: Optional[Dict[str, Any]] = None
    if aligned_path is not None:
        baseline_payload = build_matrix_target_status_artifact([aligned_path], target_tag=target_tag)
        if baseline_payload.get("benchmarks"):
            baseline_row = baseline_payload["benchmarks"][0]

    baseline_compounds = {
        str(row.get("compound", "")).strip().lower(): row
        for row in (baseline_row or {}).get("compounds", [])
    }

    delta_rows = []
    delta_counts = {"strengthened": 0, "unchanged": 0, "weakened": 0, "new_compound": 0}
    for row in current.get("compounds", []):
        compound_key = str(row.get("compound", "")).strip().lower()
        baseline_compound = baseline_compounds.get(compound_key)
        baseline_status = str((baseline_compound or {}).get("support_status", "open_gap")) if baseline_compound is not None else None
        current_status = str(row.get("support_status", "open_gap"))
        if baseline_status is None:
            delta_kind = "new_compound"
        else:
            delta_value = _SUPPORT_RANK.get(current_status, 0) - _SUPPORT_RANK.get(baseline_status, 0)
            if delta_value > 0:
                delta_kind = "strengthened"
            elif delta_value < 0:
                delta_kind = "weakened"
            else:
                delta_kind = "unchanged"
        delta_counts[delta_kind] += 1
        delta_rows.append(
            {
                **row,
                "baseline_support_status": baseline_status,
                "support_delta": delta_kind,
            }
        )

    promotion_before = bool((baseline_row or {}).get("promotion_ready", False))
    promotion_after = bool(current.get("promotion_ready", False))
    blocker_before = str((baseline_row or {}).get("promotion_blocker", "none")) if baseline_row is not None else "none"
    blocker_after = str(current.get("promotion_blocker", "none"))
    source_kind = str(normalized.get("source_kind", "internal_experiment"))
    promotion_claim_allowed = source_kind == "external_literature"
    if promotion_claim_allowed:
        promotion_claim_policy = "external_literature_can_upgrade_promotion_readiness_if_the_benchmark_contract_is_satisfied"
    elif source_kind == "internal_experiment":
        promotion_claim_policy = "internal_measurements_strengthen_calibration_and_comparator_evidence_but_do_not_unlock_external_decision_ready_claims"
    else:
        promotion_claim_policy = "diagnostic_inputs_do_not_unlock_external_decision_ready_claims"

    if not promotion_before and promotion_after:
        readiness_change = "promoted_to_external_decision_ready"
    elif delta_counts["strengthened"] > 0:
        readiness_change = "evidence_strengthened_not_yet_promotable"
    elif blocker_before != blocker_after:
        readiness_change = "promotion_blocker_shifted"
    else:
        readiness_change = "no_material_change"

    if source_kind == "external_literature" and promotion_after:
        landing_recommendation = "land_in_benchmark_payload_and_calibration_review"
    elif source_kind == "external_literature":
        landing_recommendation = "land_in_benchmark_candidate_or_blocker_registry"
    elif source_kind == "internal_experiment":
        landing_recommendation = "land_in_internal_candidate_or_calibration_payload_not_promotion_claim"
    else:
        landing_recommendation = "retain_as_diagnostic_only"

    return {
        "schema_version": "1.0",
        "description": "Support-delta comparison artifact for a matrix experiment intake payload against the current model and the closest aligned matrix benchmark.",
        "experiment": {
            "experiment_id": normalized.get("experiment_id"),
            "source_kind": source_kind,
            "protein_type": normalized.get("protein_type"),
            "process_state": normalized.get("process_state"),
            "source_reference": normalized.get("provenance", {}).get("source_reference"),
        },
        "aligned_benchmark": {
            "benchmark_id": (baseline_row or {}).get("benchmark_id"),
            "path": str(aligned_path) if aligned_path is not None else None,
        },
        "promotion_assessment": {
            "promotion_ready_before": promotion_before,
            "promotion_ready_after": promotion_after,
            "promotion_claim_allowed": promotion_claim_allowed,
            "promotion_blocker_before": blocker_before,
            "promotion_blocker_after": blocker_after,
            "promotion_claim_policy": promotion_claim_policy,
            "readiness_change": readiness_change,
            "landing_recommendation": landing_recommendation,
        },
        "support_delta": {
            "baseline_support_counts": (baseline_row or {}).get("support_counts", {}),
            "current_support_counts": current.get("support_counts", {}),
            "delta_counts": delta_counts,
        },
        "compounds": delta_rows,
    }


def calibrate_from_intake(
    payload_or_path: Mapping[str, Any] | Path | str,
    *,
    dry_run: bool = False
) -> Optional[Dict[str, Any]]:
    from src.matrix_calibration_optimizer import calibrate_matrix_constants
    
    if isinstance(payload_or_path, (Path, str)):
        payload = load_matrix_experiment_intake(payload_or_path)
    else:
        payload = dict(payload_or_path)
        
    benchmark = build_matrix_experiment_benchmark_payload(payload)
    protein_type = str(benchmark.get("protein_type", "free"))
    
    if dry_run:
        import logging
        logging.getLogger(__name__).info(f"Dry run calibration for {protein_type} using experiment {benchmark.get('benchmark_id')}")
        return {"dry_run": True, "benchmark_id": benchmark.get("benchmark_id")}
        
    return calibrate_matrix_constants([benchmark], protein_type)