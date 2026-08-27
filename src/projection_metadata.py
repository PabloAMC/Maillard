from __future__ import annotations

from typing import Any, Dict, Mapping, TypeAlias, TypedDict


class ProjectionMetadataRow(TypedDict, total=False):
    compound: str
    proxy_ppb: float
    observable_ppb: float
    observable_ratio: float
    proxy_to_observable_ratio: float
    matrix_factor: float
    base_matrix_factor: float
    class_matrix_factor: float
    dynamic_retention_factor: float
    retention_runtime_mode: str
    retention_reference_sources: list[str]
    reversible_release_factor: float
    temporal_attenuation_factor: float
    extrusion_moisture_factor: float
    extrusion_structure_factor: float
    headspace_factor: float
    calibration_factor: float
    melanoidin_trapping_factor: float
    browning_index: float
    browning_narrative: str
    volatile_class: str
    process_state: str
    accessibility_profile: str
    accessibility_warning: bool
    accessibility_dominant_source: str
    calibration_source: str
    calibration_process_state: str
    calibration_evidence_strength: str
    calibration_fallback_mode: str
    evidence_state: str
    target_class: str
    chemistry_family: str | None
    supporting_families: list[str]
    observable_panel_tags: list[str]
    panel_role: str
    observable_kind: str
    modeling_regimes: list[str]
    decision_panel_source: str
    calibration_observable_factor: float | None
    calibration_notes: str
    limiting_precursor_molar: float
    projection_load_factor: float
    projection_temperature_factor: float
    projection_time_factor: float
    projection_severity: float
    projection_kinetic_drive: float
    projection_conversion_extent: float
    volatile_yield_fraction: float
    total_volatile_budget_molar: float
    projection_strategy_name: str
    projection_precursor_unit: str
    projection_ppb_basis: str
    projection_limiting_pool_to_molar_factor: float
    projection_baseline_volatile_yield_fraction: float
    projection_conversion_ceiling_fraction: float
    projection_apparent_activation_energy_kj_mol: float
    projection_reference_temperature_kelvin: float
    projection_reference_conversion_time_min: float
    projection_ppb_conversion_factor: float
    projection_strategy_notes: str


ProjectionMetadataMap: TypeAlias = Dict[str, ProjectionMetadataRow]


def normalize_projection_metadata_row(
    row: Mapping[str, Any],
    *,
    compound_fallback: str,
    observable_ppb_fallback: float = 0.0,
) -> ProjectionMetadataRow:
    observable_ppb = float(row.get("observable_ppb", observable_ppb_fallback))
    proxy_ppb = float(row.get("proxy_ppb", observable_ppb))
    proxy_ratio = row.get("proxy_to_observable_ratio")
    if proxy_ratio is None:
        proxy_ratio = 1.0 if proxy_ppb <= 0.0 else observable_ppb / proxy_ppb

    normalized: ProjectionMetadataRow = {
        "compound": str(row.get("compound", compound_fallback)),
        "proxy_ppb": proxy_ppb,
        "observable_ppb": observable_ppb,
        "proxy_to_observable_ratio": float(proxy_ratio),
        "matrix_factor": float(row.get("matrix_factor", 1.0)),
        "base_matrix_factor": float(row.get("base_matrix_factor", row.get("matrix_factor", 1.0))),
        "class_matrix_factor": float(row.get("class_matrix_factor", 1.0)),
        "dynamic_retention_factor": float(row.get("dynamic_retention_factor", 1.0)),
        "retention_runtime_mode": str(row.get("retention_runtime_mode", "static_class_profile")),
        "retention_reference_sources": list(row.get("retention_reference_sources", [])),
        "reversible_release_factor": float(row.get("reversible_release_factor", 1.0)),
        "temporal_attenuation_factor": float(row.get("temporal_attenuation_factor", 1.0)),
        "extrusion_moisture_factor": float(row.get("extrusion_moisture_factor", 1.0)),
        "extrusion_structure_factor": float(row.get("extrusion_structure_factor", 1.0)),
        "headspace_factor": float(row.get("headspace_factor", 1.0)),
        "calibration_factor": float(row.get("calibration_factor", row.get("calibration_observable_factor") or 1.0)),
        "melanoidin_trapping_factor": float(row.get("melanoidin_trapping_factor", 1.0)),
        "browning_index": float(row.get("browning_index", row.get("projection_severity", 0.0))),
        "browning_narrative": str(row.get("browning_narrative", "")),
        "volatile_class": str(row.get("volatile_class", "other")),
        "process_state": str(row.get("process_state", "unknown")),
        "accessibility_profile": str(row.get("accessibility_profile", "unknown")),
        "accessibility_warning": bool(row.get("accessibility_warning", False)),
        "accessibility_dominant_source": str(row.get("accessibility_dominant_source", "unknown")),
        "calibration_source": str(row.get("calibration_source", "class_fallback")),
        "calibration_process_state": str(row.get("calibration_process_state", row.get("process_state", "unknown"))),
        "calibration_evidence_strength": str(row.get("calibration_evidence_strength", "heuristic")),
        "calibration_fallback_mode": str(row.get("calibration_fallback_mode", "class_level")),
        "evidence_state": str(row.get("evidence_state", "still_missing")),
        "target_class": str(row.get("target_class", "unknown")),
        "chemistry_family": row.get("chemistry_family"),
        "supporting_families": list(row.get("supporting_families", [])),
        "observable_panel_tags": list(row.get("observable_panel_tags", [])),
        "panel_role": str(row.get("panel_role", "unknown")),
        "observable_kind": str(row.get("observable_kind", "unknown")),
        "modeling_regimes": list(row.get("modeling_regimes", [])),
        "decision_panel_source": str(row.get("decision_panel_source", "")),
        "calibration_observable_factor": row.get("calibration_observable_factor"),
        "calibration_notes": str(row.get("calibration_notes", "")),
    }

    for key, value in row.items():
        if key not in normalized:
            normalized[key] = value
    return normalized


def make_projection_metadata_row(
    *,
    compound: str,
    proxy_ppb: float,
    observable_ppb: float,
    extras: Mapping[str, Any] | None = None,
) -> ProjectionMetadataRow:
    payload: Dict[str, Any] = {
        "compound": compound,
        "proxy_ppb": float(proxy_ppb),
        "observable_ppb": float(observable_ppb),
    }
    if extras:
        payload.update(dict(extras))
    return normalize_projection_metadata_row(
        payload,
        compound_fallback=compound,
        observable_ppb_fallback=float(observable_ppb),
    )