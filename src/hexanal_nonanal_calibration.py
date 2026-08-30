from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from src.benchmark_validation import load_benchmark
from src.artifact_io import repo_root


@dataclass(frozen=True)
class CalibrationLane:
    lane_id: str
    protein_type: str
    primary_benchmark: str
    reference_benchmark: str


HEXANAL_NONANAL_RATIO_BOUNDS = (0.5, 2.0)

CALIBRATION_LANES = (
    CalibrationLane(
        lane_id="pea_protocol_pilot_vs_internal2026",
        protein_type="pea_iso",
        primary_benchmark="pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026",
        reference_benchmark="pea_isolate_ribose_cysteine_100C_45min_Internal2026",
    ),
    CalibrationLane(
        lane_id="soy_protocol_pilot_vs_internal2026",
        protein_type="soy_iso",
        primary_benchmark="soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026",
        reference_benchmark="soy_isolate_ribose_cysteine_100C_45min_Internal2026",
    ),
)


def _benchmark_path(benchmark_id: str) -> Path:
    return repo_root() / "data" / "benchmarks" / f"{benchmark_id}.json"


def _volatile_entry(benchmark: Mapping[str, Any], compound: str) -> Dict[str, Any]:
    measured = benchmark.get("measured_volatiles") or {}
    reference = benchmark.get("reference_volatiles") or {}
    entry = measured.get(compound)
    if isinstance(entry, Mapping):
        return dict(entry)
    entry = reference.get(compound)
    if isinstance(entry, Mapping):
        return dict(entry)
    return {}


def _conc_ppb(benchmark: Mapping[str, Any], compound: str) -> Optional[float]:
    entry = _volatile_entry(benchmark, compound)
    value = entry.get("conc_ppb")
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _uncertainty_pct(benchmark: Mapping[str, Any], compound: str) -> Optional[float]:
    entry = _volatile_entry(benchmark, compound)
    value = entry.get("uncertainty_pct")
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _ratio(primary_value: Optional[float], reference_value: Optional[float]) -> Optional[float]:
    if primary_value is None or reference_value is None or reference_value <= 0.0:
        return None
    return primary_value / reference_value


def _in_bounds(value: Optional[float], bounds: tuple[float, float]) -> bool:
    if value is None:
        return False
    lower, upper = bounds
    return lower <= value <= upper


def _format_measurement(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{float(value):.6g}"


def build_hexanal_nonanal_calibration_artifact() -> Dict[str, Any]:
    lane_rows = []
    passed_count = 0
    prediction_change_rows = []
    for lane in CALIBRATION_LANES:
        primary = load_benchmark(_benchmark_path(lane.primary_benchmark))
        reference = load_benchmark(_benchmark_path(lane.reference_benchmark))

        hexanal_primary = _conc_ppb(primary, "Hexanal")
        hexanal_reference = _conc_ppb(reference, "Hexanal")
        nonanal_primary = _conc_ppb(primary, "Nonanal")
        nonanal_reference = _conc_ppb(reference, "Nonanal")

        hexanal_ratio = _ratio(hexanal_primary, hexanal_reference)
        nonanal_ratio = _ratio(nonanal_primary, nonanal_reference)
        hexanal_passed = _in_bounds(hexanal_ratio, HEXANAL_NONANAL_RATIO_BOUNDS)
        nonanal_passed = _in_bounds(nonanal_ratio, HEXANAL_NONANAL_RATIO_BOUNDS)
        passed = hexanal_passed and nonanal_passed
        if passed:
            passed_count += 1

        for compound_name, primary_value, reference_value, ratio_value, in_bounds in (
            ("Hexanal", hexanal_primary, hexanal_reference, hexanal_ratio, hexanal_passed),
            ("Nonanal", nonanal_primary, nonanal_reference, nonanal_ratio, nonanal_passed),
        ):
            prediction_change_rows.append(
                {
                    "lane_id": lane.lane_id,
                    "protein_type": lane.protein_type,
                    "compound": compound_name,
                    "internal2026_ppb": reference_value,
                    "protocolpilot2026_ppb": primary_value,
                    "ratio": ratio_value,
                    "closure_state": "calibration_closed" if in_bounds else "calibration_hazard",
                    "next_decision_gate": (
                        "retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence"
                        if in_bounds
                        else "revisit_headspace_retention_and_lipid_oxidation_calibration"
                    ),
                }
            )

        lane_rows.append(
            {
                "lane_id": lane.lane_id,
                "protein_type": lane.protein_type,
                "primary_benchmark_id": lane.primary_benchmark,
                "reference_benchmark_id": lane.reference_benchmark,
                "closure_action": "calibration_closed" if passed else "calibration_hazard",
                "passed": passed,
                "acceptance_bounds": {
                    "hexanal_ratio": list(HEXANAL_NONANAL_RATIO_BOUNDS),
                    "nonanal_ratio": list(HEXANAL_NONANAL_RATIO_BOUNDS),
                },
                "compounds": {
                    "Hexanal": {
                        "primary_ppb": hexanal_primary,
                        "reference_ppb": hexanal_reference,
                        "primary_uncertainty_pct": _uncertainty_pct(primary, "Hexanal"),
                        "reference_uncertainty_pct": _uncertainty_pct(reference, "Hexanal"),
                        "ratio": hexanal_ratio,
                        "in_bounds": hexanal_passed,
                    },
                    "Nonanal": {
                        "primary_ppb": nonanal_primary,
                        "reference_ppb": nonanal_reference,
                        "primary_uncertainty_pct": _uncertainty_pct(primary, "Nonanal"),
                        "reference_uncertainty_pct": _uncertainty_pct(reference, "Nonanal"),
                        "ratio": nonanal_ratio,
                        "in_bounds": nonanal_passed,
                    },
                },
                "scientific_basis": [
                    "ProtocolPilot2026 is a SYNTHETIC diagnostic payload whose volatile values were frozen from the model's own run path (audit 2026-08-26); it is not measured evidence of any kind.",
                    "Internal2026 is the frozen internal reference lane for matrix-calibration drift checks; it is likewise model output, not a measurement.",
                    "Closure here therefore means only that two snapshots of the model agree with each other (a reproducibility/drift check); it carries zero calibration evidence and does not unlock any promotion.",
                ],
                "next_best_action": (
                    "retain_internal_calibration_route_and_seek_external_quantitative_mixed_matrix_evidence"
                    if passed
                    else "revisit_headspace_retention_and_lipid_oxidation_calibration"
                ),
            }
        )

    return {
        "summary": {
            "lane_count": len(lane_rows),
            "closed_lane_count": passed_count,
            "hazard_lane_count": len(lane_rows) - passed_count,
            "closed_marker_count": sum(1 for row in prediction_change_rows if row.get("closure_state") == "calibration_closed"),
            "marker_count": len(prediction_change_rows),
            "default_ratio_bounds": list(HEXANAL_NONANAL_RATIO_BOUNDS),
            "policy": "hexanal_nonanal_protocol_pilot_to_internal2026_ratio_check_closes_internal_calibration_routes_but_not_external_promotion_claims",
            "comparator_nature": "synthetic_model_output_vs_synthetic_model_output",
            "evidence_disclaimer": (
                "Both lanes compare frozen model output against frozen model output "
                "(audit 2026-08-26); 'closed' is a reproducibility statement, not "
                "agreement with any measurement."
            ),
        },
        "lanes": lane_rows,
        "prediction_change_cascade": prediction_change_rows,
    }


def render_hexanal_nonanal_calibration_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Hexanal Nonanal Calibration Closure",
        "",
        "| Lane | Protein | Primary Benchmark | Reference Benchmark | Closure Action | Hexanal Ratio | Hexanal In Bounds | Nonanal Ratio | Nonanal In Bounds | Next Best Action |",
        "| --- | --- | --- | --- | --- | ---: | --- | ---: | --- | --- |",
    ]
    for row in payload.get("lanes", []):
        hexanal = row.get("compounds", {}).get("Hexanal", {})
        nonanal = row.get("compounds", {}).get("Nonanal", {})
        hexanal_ratio = hexanal.get("ratio")
        nonanal_ratio = nonanal.get("ratio")
        lines.append(
            f"| {row.get('lane_id', 'unknown')} | {row.get('protein_type', 'unknown')} | {row.get('primary_benchmark_id', 'unknown')} | {row.get('reference_benchmark_id', 'unknown')} | {row.get('closure_action', 'unknown')} | {'' if hexanal_ratio is None else f'{float(hexanal_ratio):.3f}'} | {'yes' if hexanal.get('in_bounds', False) else 'no'} | {'' if nonanal_ratio is None else f'{float(nonanal_ratio):.3f}'} | {'yes' if nonanal.get('in_bounds', False) else 'no'} | {row.get('next_best_action', 'unknown')} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Lanes checked: {int(summary.get('lane_count', 0))}",
            f"Closed lanes: {int(summary.get('closed_lane_count', 0))}",
            f"Hazard lanes: {int(summary.get('hazard_lane_count', 0))}",
            f"Closed markers: {int(summary.get('closed_marker_count', 0))} / {int(summary.get('marker_count', 0))}",
            f"Default ratio bounds: {summary.get('default_ratio_bounds', [])}",
            f"Policy: {summary.get('policy', 'unknown')}",
        ]
    )

    lines.extend(
        [
            "",
            "## Prediction Validation Chain",
            "",
            "| Lane | Protein | Compound | Internal2026 ppb | ProtocolPilot2026 ppb | Ratio | Closure State | Next Decision Gate |",
            "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for row in payload.get("prediction_change_cascade", []):
        ratio_value = row.get("ratio")
        lines.append(
            f"| {row.get('lane_id', 'unknown')} | {row.get('protein_type', 'unknown')} | {row.get('compound', 'unknown')} | {_format_measurement(row.get('internal2026_ppb'))} | {_format_measurement(row.get('protocolpilot2026_ppb'))} | {'' if ratio_value is None else f'{float(ratio_value):.3f}'} | {row.get('closure_state', 'unknown')} | {row.get('next_decision_gate', 'unknown')} |"
        )
    return "\n".join(lines) + "\n"