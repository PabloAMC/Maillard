from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src.family_validation_overview import build_family_validation_overview_artifact
from src.literature_family_registry import build_family_payload_coverage_artifact
from src.literature_family_registry import iter_matrix_decision_panel_entries
from src.family_ingestion_plan import load_family_ingestion_plan
from src.family_promotion_state import build_family_promotion_state_artifact
from src.matrix_prior_registry import summarize_matrix_prior_bundle


ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"
BENCHMARK_DIR = ROOT / "data" / "benchmarks"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
PROCESS_GAP_REGISTRY_PATH = DATA_LIT_DIR / "process_gap_registry.json"
PROCESS_STATE_CALIBRATIONS_PATH = DATA_LIT_DIR / "process_state_calibrations.json"
COMPUTATIONAL_PRIORS_PATH = DATA_LIT_DIR / "computational_priors.json"
SAFETY_REFERENCE_PAYLOADS_PATH = DATA_LIT_DIR / "safety_reference_payloads.json"


READY_STATUS_DEFAULT_TEMPLATE_KIND = {
    "ready_for_reference_encoding": "safety_payload",
    "ready_for_calibration_encoding": "process_state_calibration",
    "ready_for_intake_encoding": "benchmark_payload",
    "ready_for_directional_prior_encoding": "computational_prior",
}

ARTIFACT_TYPE_TO_TEMPLATE_KIND = {
    "benchmark": "benchmark_payload",
    "process_state_calibration": "process_state_calibration",
    "computational_prior": "computational_prior",
    "directional_prior": "computational_prior",
    "safety_reference": "safety_payload",
    "flavor_reference_payload": "flavor_reference_payload",
    "retention_payload": "retention_payload",
    "structural_gap_entry": "structural_gap_entry",
}

PROMOTION_QUEUE_SLR_FAMILIES = {"03", "04", "05", "06", "07", "10", "11", "12", "13", "14", "15", "16"}

PROMOTION_POSTURE_WEIGHTS = {
    "immediate_expansion_lane": 5.0,
    "upstream_pretreatment_lane": 4.0,
    "high_value_support_lane": 3.5,
    "matrix_scope_lane": 2.0,
}


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _to_repo_relative(path: Path, root: Path = ROOT) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _normalize_name(name: str) -> str:
    lowered = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    return "_".join(part for part in lowered.split("_") if part)


def _load_benchmark_payloads() -> Dict[str, Dict[str, Any]]:
    payloads: Dict[str, Dict[str, Any]] = {}
    for bench_path in sorted(BENCHMARK_DIR.glob("*.json")):
        payload = _load_json(bench_path)
        benchmark_id = str(payload.get("benchmark_id", bench_path.stem))
        payloads[benchmark_id] = payload
    return payloads


def _iter_prior_entries(payload: Mapping[str, Any]) -> Iterable[Dict[str, Any]]:
    for section_name, section_payload in payload.items():
        if not isinstance(section_payload, list):
            continue
        for entry in section_payload:
            if isinstance(entry, Mapping):
                row = dict(entry)
                row.setdefault("_section_name", str(section_name))
                yield row


def _iter_payload_entries(payload: Any, *, section_name: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    if isinstance(payload, list):
        for entry in payload:
            if isinstance(entry, Mapping):
                row = dict(entry)
                if section_name is not None:
                    row.setdefault("_section_name", section_name)
                yield row
        return
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            if isinstance(value, list):
                next_section = str(key) if section_name is None else section_name
                yield from _iter_payload_entries(value, section_name=next_section)
            elif isinstance(value, Mapping):
                next_section = str(key) if section_name is None else section_name
                yield from _iter_payload_entries(value, section_name=next_section)


def _canonical_template_kind_for_artifact(artifact_type: str) -> str:
    normalized = str(artifact_type).strip()
    return ARTIFACT_TYPE_TO_TEMPLATE_KIND.get(normalized, normalized)


def _infer_target_payload_types(entry: Mapping[str, Any]) -> List[str]:
    payload_types: List[str] = []
    for artifact in entry.get("runtime_artifacts", []) or []:
        template_kind = _canonical_template_kind_for_artifact(str(artifact.get("artifact_type", "")))
        if template_kind and template_kind not in payload_types:
            payload_types.append(template_kind)
    if payload_types:
        return payload_types
    fallback = READY_STATUS_DEFAULT_TEMPLATE_KIND.get(str(entry.get("status", "")), "reference_payload")
    return [fallback]


def _resolve_primary_template_kind(entry: Mapping[str, Any]) -> str:
    payload_types = _infer_target_payload_types(entry)
    return payload_types[0] if payload_types else "reference_payload"


def _artifact_exists(artifact: Mapping[str, Any], *, root: Path = ROOT) -> bool:
    path = root / str(artifact.get("path", ""))
    if not path.exists():
        return False
    artifact_id = str(artifact.get("artifact_id", "")).strip()
    artifact_type = str(artifact.get("artifact_type", "")).strip()
    if not artifact_id or artifact_type == "benchmark":
        return True
    payload = _load_json(path)
    if artifact_type == "intake_registry_entry":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("eligible_references", []))
    if artifact_type == "slr_incorporation_ledger":
        return any(str(entry.get("paper_id", "")) == artifact_id for entry in payload.get("entries", []))
    if artifact_type == "process_state_calibration":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("entries", []))
    if artifact_type == "safety_reference":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("entries", []))

    sections = artifact.get("sections", [])
    if isinstance(sections, list) and sections:
        for section_name in sections:
            for entry in _iter_payload_entries(payload.get(str(section_name), []), section_name=str(section_name)):
                if str(entry.get("id", entry.get("protein_type", ""))) == artifact_id:
                    return True
        return False

    if isinstance(payload, Mapping):
        for entry in _iter_payload_entries(payload):
            if str(entry.get("id", entry.get("protein_type", ""))) == artifact_id:
                return True
    return False


def _ready_reference_rows(root: Path = ROOT) -> List[Dict[str, Any]]:
    intake = _load_json(root / _to_repo_relative(INTAKE_REGISTRY_PATH, ROOT))
    rows: List[Dict[str, Any]] = []
    for entry in intake.get("eligible_references", []):
        status = str(entry.get("status", ""))
        if status not in READY_STATUS_DEFAULT_TEMPLATE_KIND:
            continue
        runtime_artifacts = []
        for artifact in entry.get("runtime_artifacts", []) or []:
            artifact_row = dict(artifact)
            artifact_row["template_kind"] = _canonical_template_kind_for_artifact(str(artifact_row.get("artifact_type", "")))
            artifact_row["exists"] = _artifact_exists(artifact_row, root=root)
            runtime_artifacts.append(artifact_row)
        all_exist = bool(runtime_artifacts) and all(bool(item.get("exists", False)) for item in runtime_artifacts)
        encoding_status = "encoded_runtime_artifact" if all_exist else "template_required"
        target_payload_types = _infer_target_payload_types(entry)
        rows.append(
            {
                "id": str(entry.get("id", "unknown")),
                "citation": str(entry.get("citation", "unknown")),
                "kind": str(entry.get("kind", "unknown")),
                "chemistry_family": str(entry.get("chemistry_family", "")),
                "slr_family_source": str(entry.get("slr_family_source", "")),
                "source_payload_role": str(entry.get("payload_role", "unknown")),
                "target_payload_types": target_payload_types,
                "matrix_family": str(entry.get("matrix_family", "unknown")),
                "status": status,
                "template_kind": _resolve_primary_template_kind(entry),
                "requires_primary_data": bool(entry.get("requires_primary_data", False)),
                "target_modules": [str(item) for item in entry.get("target_modules", []) or []],
                "encoding_status": encoding_status,
                "runtime_artifacts": runtime_artifacts,
            }
        )
    rows.sort(key=lambda row: (row["template_kind"], row["id"]))
    return rows


def _resolve_protein_type(matrix_family: str) -> str:
    normalized = _normalize_name(matrix_family)
    if normalized.startswith("pea"):
        return "pea_iso"
    if normalized.startswith("soy"):
        return "soy_iso"
    if "myco" in normalized:
        return "myco"
    return "free"


def _suggest_reference_section(entry: Mapping[str, Any]) -> str:
    chemistry_family = str(entry.get("chemistry_family", "")).strip()
    observable_tags = {str(item).strip() for item in entry.get("observable_panel_tags", []) or []}
    if chemistry_family in {"amino_acid_sugar_core", "thiamine_fragmentation_support", "glutathione_and_peptide_support"}:
        return "sulfur_reference_anchors"
    if chemistry_family == "carbonyl_donor_hierarchy":
        return "strecker_reference_anchors"
    if chemistry_family == "carbohydrate_pyrolysis_and_caramelization":
        return "furanone_reference_anchors" if "furanone" in observable_tags else "carbonyl_reference_anchors"
    if chemistry_family in {"lipid_oxidation_and_carbonylic_crosstalk", "off_note_and_maillard_suppression"}:
        return "carbonyl_reference_anchors"
    return "reference_section_to_define"


def _suggest_retention_section(entry: Mapping[str, Any]) -> str:
    observable_tags = {str(item).strip() for item in entry.get("observable_panel_tags", []) or []}
    if "off_note" in observable_tags or "adverse_marker" in observable_tags:
        return "aldehydes"
    return "sulfur_compounds"


def _build_template_for_entry(entry: Mapping[str, Any]) -> Dict[str, Any]:
    entry_id = str(entry.get("id", "unknown"))
    template_kind = _resolve_primary_template_kind(entry)
    key_values = dict(entry.get("key_values", {}))
    protein_type = _resolve_protein_type(str(entry.get("matrix_family", "unknown")))
    base_template = {
        "entry_id": entry_id,
        "template_kind": template_kind,
        "source_payload_role": str(entry.get("payload_role", "unknown")),
        "target_payload_types": _infer_target_payload_types(entry),
        "chemistry_family": str(entry.get("chemistry_family", "")),
        "slr_family_source": str(entry.get("slr_family_source", "")),
        "matrix_family": str(entry.get("matrix_family", "unknown")),
        "observable_panel_tags": [str(item) for item in entry.get("observable_panel_tags", []) or []],
        "process_state_scope": [str(item) for item in entry.get("process_state_scope", []) or []],
    }

    if template_kind == "benchmark_payload":
        tracked = dict(key_values.get("tracked_uht_markers_ug_per_l", {}))
        measured_volatiles = {
            compound: {"conc_ppb": float(value), "uncertainty_pct": 10.0}
            for compound, value in tracked.items()
        }
        observable_targets = [
            {
                "name": compound,
                "role": "adverse_marker",
                "expected_rank": index + 1,
                "direction": "higher",
            }
            for index, (compound, _value) in enumerate(sorted(tracked.items(), key=lambda item: float(item[1]), reverse=True)[:3])
        ]
        output_name = "pea_isolate_uht_140C_Trikusuma2019.json"
        return {
            **base_template,
            "suggested_output_path": f"data/benchmarks/{output_name}",
            "benchmark_id": "pea_isolate_uht_140C_Trikusuma2019",
            "protein_type": protein_type,
            "conditions": {
                "temp_C": float(key_values.get("uht_temp_C", 140.0)),
                "ph": 7.1,
                "water_activity": 0.98,
                "time_min": float(key_values.get("uht_hold_seconds", 6.0)) / 60.0,
            },
            "matrix_ranking_contract": {
                "observable_targets": observable_targets,
                "adverse_markers": [item["name"] for item in observable_targets],
                "calibration_mode": "compound_specific_headspace",
            },
            "measured_volatiles": measured_volatiles,
        }

    if template_kind == "process_state_calibration":
        return {
            **base_template,
            "suggested_output_path": "data/lit/process_state_calibrations.json",
            "protein_type": protein_type,
            "kind": str(entry.get("kind", "calibration_reference")),
            "numeric_anchors": key_values,
            "source_citation": str(entry.get("citation", "unknown")),
        }

    if template_kind == "computational_prior":
        return {
            **base_template,
            "suggested_output_path": "data/lit/computational_priors.json",
            "target_section": "strecker_crosstalk_priors",
            "protein_type": protein_type,
            "required_context": ["lipids", "polyphenols"],
            "key_values": key_values,
        }

    if template_kind == "flavor_reference_payload":
        return {
            **base_template,
            "suggested_output_path": "data/lit/flavor_reference_payloads.json",
            "target_section": _suggest_reference_section(entry),
            "source_citation": str(entry.get("citation", "unknown")),
            "key_values": key_values,
        }

    if template_kind == "retention_payload":
        return {
            **base_template,
            "suggested_output_path": "data/lit/retention_reference_payloads.json",
            "target_section": _suggest_retention_section(entry),
            "target_matrix_family": str(entry.get("matrix_family", "unknown")),
            "source_citation": str(entry.get("citation", "unknown")),
            "key_values": key_values,
        }

    return {
        **base_template,
        "suggested_output_path": "data/lit/safety_reference_payloads.json",
        "target_section": "entries",
        "analyte": "acrylamide",
        "numeric_anchors": key_values,
        "source_citation": str(entry.get("citation", "unknown")),
    }


def build_runtime_templates(root: Path = ROOT) -> List[Dict[str, Any]]:
    intake = _load_json(root / _to_repo_relative(INTAKE_REGISTRY_PATH, ROOT))
    templates: List[Dict[str, Any]] = []
    ready_rows = {row["id"]: row for row in _ready_reference_rows(root)}
    for entry in intake.get("eligible_references", []):
        entry_id = str(entry.get("id", "unknown"))
        if entry_id not in ready_rows:
            continue
        template = _build_template_for_entry(entry)
        template["encoding_status"] = ready_rows[entry_id]["encoding_status"]
        templates.append(template)
    templates.sort(key=lambda row: (row["template_kind"], row["entry_id"]))
    return templates


def build_payload_queue_review(root: Path = ROOT) -> Dict[str, Any]:
    ready_rows = _ready_reference_rows(root)
    queue_by_payload_type: Dict[str, int] = {}
    queue_by_chemistry_family: Dict[str, Dict[str, int]] = {}
    for row in ready_rows:
        family = str(row.get("chemistry_family", "unknown")) or "unknown"
        family_counts = queue_by_chemistry_family.setdefault(family, {})
        for payload_type in row.get("target_payload_types", []) or []:
            queue_by_payload_type[payload_type] = queue_by_payload_type.get(payload_type, 0) + 1
            family_counts[payload_type] = family_counts.get(payload_type, 0) + 1

    family_rows = []
    for chemistry_family in sorted(queue_by_chemistry_family):
        counts = queue_by_chemistry_family[chemistry_family]
        family_rows.append(
            {
                "chemistry_family": chemistry_family,
                "payload_type_counts": dict(sorted(counts.items())),
                "ready_reference_count": sum(int(value) for value in counts.values()),
            }
        )
    return {
        "queue_by_payload_type": dict(sorted(queue_by_payload_type.items())),
        "queue_by_chemistry_family": family_rows,
    }


def _promotion_benchmarkability_weight(status: str) -> float:
    normalized = str(status).strip()
    if "benchmark_ready" in normalized:
        return 3.0
    if "calibration_ready" in normalized:
        return 2.0
    if "bounded_quantitative_support" in normalized:
        return 1.5
    if "mostly_directional" in normalized:
        return 0.5
    return 1.0


def _direct_matrix_decision_bonus(plan_row: Mapping[str, Any]) -> float:
    targets = {str(item).strip().lower() for item in plan_row.get("target_compounds_or_state_variables", []) or []}
    bonus = 0.0
    if targets.intersection({"ribose", "xylose", "glucose", "fructose", "ribose-5-phosphate"}):
        bonus += 1.0
    if targets.intersection({"2-methyl-3-furanthiol", "2-furfurylthiol", "bis(2-methyl-3-furyl) disulfide"}):
        bonus += 1.0
    if targets.intersection({"fermentation_pretreatment_state", "free_amino_acid_enrichment", "nucleotide_enrichment"}):
        bonus += 0.25
    return bonus


def _minimum_runtime_landing_kind(
    plan_row: Mapping[str, Any],
    coverage_row: Mapping[str, Any],
    validation_row: Mapping[str, Any],
) -> str:
    preferred = [str(item) for item in plan_row.get("preferred_payload_types", []) or []]
    if "benchmark_payload" in preferred and int(validation_row.get("benchmark_count", 0)) == 0:
        return "benchmark_payload"
    if "process_state_calibration" in preferred and int(coverage_row.get("total_primary_payload_count", 0)) == 0:
        return "process_state_calibration"
    if "benchmark_payload" in preferred and int(validation_row.get("quantitative_point_count", 0)) == 0:
        return "benchmark_payload"
    return "bounded_prior"


def _runtime_landing_contract_note(landing_kind: str) -> str:
    if landing_kind == "benchmark_payload":
        return "must land as a benchmark-linked payload that changes the matrix decision surface or trust posture"
    if landing_kind == "process_state_calibration":
        return "must land as a runtime-consumed process-state calibration, not narrative-only markdown"
    return "must land as a bounded prior with explicit uncertainty_posture and process_state_applicability"


def _family_promotion_blocker(
    plan_row: Mapping[str, Any],
    coverage_row: Mapping[str, Any],
    validation_row: Mapping[str, Any],
    promotion_row: Mapping[str, Any],
    landing_kind: str,
) -> str:
    promotion_state = str(promotion_row.get("promotion_state", "unknown"))
    if promotion_state == "benchmark_linked_support":
        return "already promoted to benchmark-linked support; standalone primary benchmark payload is still missing"
    if promotion_state == "calibration_grade":
        return "already promoted to calibration-grade support; direct benchmark-linked closure is still missing"
    if promotion_state == "near_quantitative":
        return "already near-quantitative"
    blocked_by = [str(item) for item in plan_row.get("blocked_by", []) or []]
    if blocked_by:
        return "; ".join(blocked_by)
    if landing_kind == "benchmark_payload" and int(validation_row.get("benchmark_count", 0)) == 0:
        return "no benchmark-linked closure yet for a family that directly changes formulation choices"
    if landing_kind == "process_state_calibration" and int(coverage_row.get("total_primary_payload_count", 0)) == 0:
        return "no process-state calibration encoded yet"
    if landing_kind == "bounded_prior" and int(coverage_row.get("total_primary_payload_count", 0)) == 0:
        return "no bounded prior encoded with explicit uncertainty posture"
    if int(validation_row.get("quantitative_point_count", 0)) == 0:
        return "runtime support exists but no benchmark-visible decision change is yet proven"
    return "none"


def build_s11_c_family_promotion_queue(root: Path = ROOT) -> Dict[str, Any]:
    plan = load_family_ingestion_plan()
    coverage_payload = build_family_payload_coverage_artifact()
    promotion_payload = build_family_promotion_state_artifact()
    coverage_rows = {
        str(row.get("family_id", "unknown")): row
        for row in coverage_payload.get("families", [])
    }
    promotion_rows = {
        str(row.get("chemistry_family", "unknown")): row
        for row in promotion_payload.get("families", [])
    }
    validation_payload = build_family_validation_overview_artifact()
    validation_rows = {
        str(row.get("chemistry_family", "unknown")): row
        for row in validation_payload.get("families", [])
    }

    candidates: List[Dict[str, Any]] = []
    for family in plan.get("families", []):
        slr_family = str(family.get("slr_family", "")).zfill(2)
        if slr_family not in PROMOTION_QUEUE_SLR_FAMILIES:
            continue

        family_id = str(family.get("family_id", "unknown"))
        coverage_row = coverage_rows.get(family_id, {})
        validation_row = validation_rows.get(family_id, {})
        promotion_row = promotion_rows.get(family_id, {})
        decision_panel_entry_count = len(list(iter_matrix_decision_panel_entries(family=family_id)))
        impact_score = (
            float(PROMOTION_POSTURE_WEIGHTS.get(str(family.get("strategic_posture", "unknown")), 1.5))
            + _promotion_benchmarkability_weight(str(family.get("benchmarkability_status", "unknown")))
            + _direct_matrix_decision_bonus(family)
            + (1.25 * float(min(decision_panel_entry_count, 4)))
            + (0.10 * float(min(int(coverage_row.get("total_supporting_payload_count", 0)), 20)))
        )
        landing_kind = _minimum_runtime_landing_kind(family, coverage_row, validation_row)
        blocker = _family_promotion_blocker(family, coverage_row, validation_row, promotion_row, landing_kind)
        candidates.append(
            {
                "slr_family": slr_family,
                "family_id": family_id,
                "display_name": str(family.get("display_name", family_id)),
            "current_promotion_state": str(promotion_row.get("promotion_state", "unknown")),
                "strategic_posture": str(family.get("strategic_posture", "unknown")),
                "benchmarkability_status": str(family.get("benchmarkability_status", "unknown")),
                "decision_panel_entry_count": decision_panel_entry_count,
                "primary_payload_count": int(coverage_row.get("total_primary_payload_count", 0)),
                "supporting_payload_count": int(coverage_row.get("total_supporting_payload_count", 0)),
                "benchmark_count": int(validation_row.get("benchmark_count", 0)),
                "quantitative_point_count": int(validation_row.get("quantitative_point_count", 0)),
                "minimum_runtime_landing": landing_kind,
                "runtime_landing_contract": _runtime_landing_contract_note(landing_kind),
                "promotion_blocker": blocker,
                "reject_narrative_only": True,
                "impact_score": round(impact_score, 2),
            }
        )

    candidates.sort(key=lambda row: (-float(row.get("impact_score", 0.0)), row.get("slr_family", "99")))
    for index, row in enumerate(candidates):
        if row.get("current_promotion_state") in {"benchmark_linked_support", "calibration_grade", "near_quantitative"}:
            row["promotion_recommendation"] = "already_promoted"
        elif index == 0:
            row["promotion_recommendation"] = "select_now"
        elif index == 1:
            row["promotion_recommendation"] = "second_if_first_blocked"
        else:
            row["promotion_recommendation"] = "defer"

    selected_family = candidates[0] if candidates else None
    fallback_family = candidates[1] if len(candidates) > 1 else None
    return {
        "selection_policy": "rank only families 03/04/05/06/07/10 by direct matrix decision impact, require runtime landing, and reject narrative-only promotion",
        "selected_family": selected_family,
        "fallback_family": fallback_family,
        "candidates": candidates,
    }


def build_matrix_prior_review() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for protein_type in ["pea_iso", "soy_iso", "myco"]:
        bundle = summarize_matrix_prior_bundle(protein_type)
        rows.append(
            {
                "protein_type": protein_type,
                "has_accessibility_window": "accessibility_window" in bundle,
                "has_denaturation_heuristic": "denaturation_heuristic" in bundle,
                "has_volatile_class_profile": "volatile_class_profile" in bundle,
                "has_matrix_correction": "matrix_correction" in bundle,
                "provenance_tiers": sorted(
                    {str(item.get("provenance_tier", "unknown")) for item in bundle.values()}
                ),
                "uncertainty_postures": sorted(
                    {str(item.get("uncertainty_posture", "unknown")) for item in bundle.values()}
                ),
                "process_state_applicability": sorted(
                    {
                        str(state)
                        for item in bundle.values()
                        for state in item.get("process_state_applicability", [])
                    }
                ),
            }
        )
    return rows


def build_literature_gap_review(root: Path = ROOT) -> Dict[str, List[Dict[str, Any]]]:
    intake = _load_json(root / _to_repo_relative(INTAKE_REGISTRY_PATH, ROOT))
    process_gaps = _load_json(root / _to_repo_relative(PROCESS_GAP_REGISTRY_PATH, ROOT))
    intake_rows = []
    for entry in intake.get("structural_gaps", []):
        intake_rows.append(
            {
                "gap_id": str(entry.get("id", "unknown")),
                "priority": str(entry.get("priority", "unknown")),
                "requires_primary_data": bool(entry.get("requires_primary_data", True)),
                "why": str(entry.get("why", "")),
            }
        )
    process_rows = []
    for entry in process_gaps.get("entries", []):
        process_rows.append(
            {
                "gap_id": str(entry.get("gap_id", "unknown")),
                "gap_type": str(entry.get("gap_type", "unknown")),
                "wet_lab_requirement": str(entry.get("wet_lab_requirement", "unknown")),
                "blocks_modules": [str(item) for item in entry.get("blocks_modules", []) or []],
                "computational_fallback": str(entry.get("computational_fallback", "")),
            }
        )
    return {
        "intake_structural_gap_review": intake_rows,
        "process_gap_review": process_rows,
    }


def build_literature_learning_loop_payload(root: Path = ROOT) -> Dict[str, Any]:
    ready_rows = _ready_reference_rows(root)
    templates = build_runtime_templates(root)
    prior_review = build_matrix_prior_review()
    gap_review = build_literature_gap_review(root)
    family_payload_coverage = build_family_payload_coverage_artifact()
    payload_queue_review = build_payload_queue_review(root)
    s11_c_family_promotion_queue = build_s11_c_family_promotion_queue(root)
    encoded_count = sum(1 for row in ready_rows if row["encoding_status"] == "encoded_runtime_artifact")
    return {
        "schema_version": "1.0",
        "source": _to_repo_relative(INTAKE_REGISTRY_PATH, root),
        "ready_reference_rows": ready_rows,
        "runtime_templates": templates,
        "matrix_prior_review": prior_review,
        "family_payload_coverage": family_payload_coverage,
        "payload_queue_review": payload_queue_review,
        "s11_c_family_promotion_queue": s11_c_family_promotion_queue,
        **gap_review,
        "summary": {
            "ready_reference_count": len(ready_rows),
            "encoded_runtime_reference_count": encoded_count,
            "template_queue_count": len(templates),
            "matrix_family_count": len({row["matrix_family"] for row in ready_rows}),
            "matrix_prior_families": [row["protein_type"] for row in prior_review],
            "families_with_primary_payload_support": int(family_payload_coverage.get("summary", {}).get("families_with_primary_payload_support", 0)),
            "payload_type_queue": dict(payload_queue_review.get("queue_by_payload_type", {})),
            "selected_s11_c_family": (s11_c_family_promotion_queue.get("selected_family") or {}).get("family_id"),
            "intake_structural_gap_count": len(gap_review["intake_structural_gap_review"]),
            "process_gap_count": len(gap_review["process_gap_review"]),
        },
    }


def render_literature_learning_loop_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Literature Learning Loop",
        "",
        "## Ready References",
        "",
        "| ID | Kind | Chemistry Family | Source Payload Role | Target Payload Types | Matrix Family | Template Kind | Encoding Status | Runtime Artifacts |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("ready_reference_rows", []):
        artifact_labels = ", ".join(
            f"{item.get('artifact_type', 'artifact')}:{item.get('artifact_id', 'unknown')}"
            for item in row.get("runtime_artifacts", [])
        ) or "none"
        lines.append(
            f"| {row['id']} | {row['kind']} | {row['chemistry_family']} | {row['source_payload_role']} | {', '.join(row.get('target_payload_types', [])) or 'none'} | {row['matrix_family']} | {row['template_kind']} | {row['encoding_status']} | {artifact_labels} |"
        )

    lines.extend([
        "",
        "## Payload Queue Review",
        "",
        "| Chemistry Family | Ready References | Payload Type Counts |",
        "| --- | ---: | --- |",
    ])
    for row in payload.get("payload_queue_review", {}).get("queue_by_chemistry_family", []):
        counts = ", ".join(f"{key}:{value}" for key, value in row.get("payload_type_counts", {}).items()) or "none"
        lines.append(
            f"| {row['chemistry_family']} | {row['ready_reference_count']} | {counts} |"
        )

    lines.extend([
        "",
        "## Family Payload Coverage",
        "",
        "| SLR | Family | Primary Payloads | Supporting Payloads | Runtime Support |",
        "| --- | --- | ---: | ---: | --- |",
    ])
    for row in payload.get("family_payload_coverage", {}).get("families", []):
        lines.append(
            f"| {row['slr_family']} | {row['family_id']} | {row['total_primary_payload_count']} | {row['total_supporting_payload_count']} | {row['has_runtime_support']} |"
        )

    queue = payload.get("s11_c_family_promotion_queue", {}) or {}
    selected_family = queue.get("selected_family") or {}
    fallback_family = queue.get("fallback_family") or {}
    if queue.get("candidates"):
        lines.extend([
            "",
            "## S11.C Family Promotion Queue",
            "",
            f"Selected family: {selected_family.get('family_id', 'unknown')}",
            f"Fallback family: {fallback_family.get('family_id', 'none')}",
            f"Selection policy: {queue.get('selection_policy', 'unknown')}",
            "",
            "| SLR | Family | Impact | Recommendation | Landing | Primary Payloads | Benchmarks | Quantitative Points | Blocker |",
            "| --- | --- | ---: | --- | --- | ---: | ---: | ---: | --- |",
        ])
        for row in queue.get("candidates", []):
            lines.append(
                f"| {row.get('slr_family', '99')} | {row.get('family_id', 'unknown')} | {float(row.get('impact_score', 0.0)):.2f} | {row.get('promotion_recommendation', 'defer')} | {row.get('minimum_runtime_landing', 'unknown')} | {int(row.get('primary_payload_count', 0))} | {int(row.get('benchmark_count', 0))} | {int(row.get('quantitative_point_count', 0))} | {row.get('promotion_blocker', 'unknown')} |"
            )

    lines.extend([
        "",
        "## Matrix Prior Review",
        "",
        "| Protein Type | Accessibility | Denaturation | Matrix Correction | Uncertainty | Process States |",
        "| --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload.get("matrix_prior_review", []):
        lines.append(
            f"| {row['protein_type']} | {row['has_accessibility_window']} | {row['has_denaturation_heuristic']} | {row['has_matrix_correction']} | {', '.join(row['uncertainty_postures'])} | {', '.join(row['process_state_applicability'])} |"
        )

    lines.extend([
        "",
        "## Structural Gaps",
        "",
    ])
    for row in payload.get("intake_structural_gap_review", []):
        lines.append(f"- intake gap {row['gap_id']}: primary_data={row['requires_primary_data']} priority={row['priority']}")
    for row in payload.get("process_gap_review", []):
        lines.append(f"- process gap {row['gap_id']}: wet_lab_requirement={row['wet_lab_requirement']}")

    return "\n".join(lines) + "\n"