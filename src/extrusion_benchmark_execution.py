from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import yaml

from src import data_paths
from src.extrusion import build_extrusion_process_profile, compute_extrusion_headspace_adjustment
from src.matrix_experiment_intake import build_matrix_experiment_support_delta_artifact


ROOT = data_paths.REPO_ROOT


def _under_root(root: Path, default: Path) -> Path:
    """``default`` when ``root`` is the repository checkout; the same repo-relative file
    re-rooted under ``root`` otherwise (tests and scripts pass explicit roots)."""
    return default if root == data_paths.REPO_ROOT else root / data_paths.rel(default)


def _load_yaml(path: Path | str) -> Dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("Workbook payload must be a mapping")
    return payload


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _load_benchmark(benchmark_id: str, root: Path = ROOT) -> Dict[str, Any]:
    path = _under_root(root, data_paths.benchmark_path(benchmark_id))
    payload = _load_json(path)
    if not payload:
        raise ValueError(f"Benchmark id not found: {benchmark_id}")
    return payload


def _contains_placeholder(value: Any) -> bool:
    if isinstance(value, str):
        return value.startswith("REPLACE_WITH_")
    if isinstance(value, Mapping):
        return any(_contains_placeholder(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_placeholder(item) for item in value)
    return False


def _normalize_float_map(rows: Mapping[str, Any]) -> Dict[str, float]:
    values: Dict[str, float] = {}
    for key, value in rows.items():
        values[str(key)] = float(value)
    return values


def _matrix_precursors(protein_type: str) -> Dict[str, Dict[str, float]]:
    protein_label = "Soy Protein Isolate" if protein_type == "soy_iso" else "Pea Protein Isolate"
    return {
        protein_label: {"concentration_mM": 1000.0},
        "D-Ribose": {"concentration_mM": 1.0},
        "L-Cysteine": {"concentration_mM": 1.0},
    }


def _matrix_contract_for_extrusion(measured_volatiles: Mapping[str, float]) -> Dict[str, Any]:
    desirable = [
        "2-methyl-3-furanthiol",
        "2-furfurylthiol",
        "bis(2-methyl-3-furyl) disulfide",
        "2,5-dimethylpyrazine",
    ]
    adverse = ["Hexanal", "2-pentylfuran", "furfural"]
    observable_targets = []
    rank = 1
    for name in desirable:
        if name in measured_volatiles:
            observable_targets.append({"name": name, "role": "desirable_marker", "expected_rank": rank, "direction": "higher"})
            rank += 1
    for name in adverse:
        if name in measured_volatiles:
            observable_targets.append({"name": name, "role": "adverse_marker", "expected_rank": rank, "direction": "lower"})
            rank += 1
    return {
        "observable_targets": observable_targets,
        "adverse_markers": [name for name in ["furfural", "Hexanal", "Nonanal", "2-pentylfuran"] if name in measured_volatiles],
        "citation_provenance": ["Extrusion external closure workbook landed through the S17.1 workbook processor."],
        "notes": "Acrylamide and direct damage markers are tracked in the extrusion closure summary, not in the matrix support-delta comparator.",
    }


def build_extrusion_external_closure_execution_artifact(
    workbook_or_path: Mapping[str, Any] | Path | str,
    *,
    root: Path = ROOT,
    allow_placeholders: bool = False,
) -> Dict[str, Any]:
    workbook = _load_yaml(workbook_or_path) if isinstance(workbook_or_path, (Path, str)) else dict(workbook_or_path)
    experiments = list(workbook.get("experiments", []))
    if not experiments:
        raise ValueError("Extrusion closure workbook must include at least one experiment")

    experiment_rows = []
    for experiment in experiments:
        if not allow_placeholders and _contains_placeholder(experiment):
            raise ValueError(f"Workbook experiment {experiment.get('experiment_id', 'unknown')} still contains placeholder values")

        measured_volatiles = _normalize_float_map(experiment.get("measured_volatiles_ppb", {}))
        if not measured_volatiles:
            raise ValueError("Extrusion closure workbook experiment must include measured_volatiles_ppb")
        measured_damage = _normalize_float_map(experiment.get("measured_damage_markers", {}))
        damage_complete = all(name in measured_damage for name in ["reactive_lysine_fraction", "furosine_mg_per_kg", "lysinoalanine_mg_per_kg"])
        source_kind = str(experiment.get("source_kind", "external_literature"))
        target_process = dict(experiment.get("target_process_settings", {}))
        protein_type = str(experiment.get("protein_type", workbook.get("selected_matrix", "soy_iso")))
        intake_payload = {
            "experiment_id": str(experiment.get("experiment_id", "extrusion_external_closure")),
            "source_kind": source_kind,
            "protein_type": protein_type,
            "process_state": str(experiment.get("process_state", "extrusion_structured")),
            "conditions": {
                "temp_C": float(target_process.get("effective_temperature_celsius", target_process.get("barrel_temperature_celsius", 145.0))),
                "ph": 5.8,
                "water_activity": 0.75,
                "time_min": max(0.1, float(target_process.get("mean_residence_time_seconds", 30.0)) / 60.0),
            },
            "formulation": {"precursors": _matrix_precursors(protein_type)},
            "measured_volatiles": {
                name: {"conc_ppb": value, "uncertainty_pct": 10.0}
                for name, value in measured_volatiles.items()
                if str(name) != "acrylamide"
            },
            "provenance": dict(experiment.get("provenance", {})),
            "matrix_format": "extruded protein matrix with ribose and cysteine benchmark precursors",
            "benchmark_alignment": dict(experiment.get("benchmark_alignment", {})),
            "comparison_contract": _matrix_contract_for_extrusion(measured_volatiles),
            "analytical_context": {
                "headspace_method": "HS-SPME-GC-MS",
                "quantification_mode": "internal_standard_calibrated",
                "replicates": 3,
                "non_detect_policy": "report_lod_and_do_not_backfill",
                "internal_standards": ["thiol_specific_internal_standard", "hexanal-d12"],
                "notes": "Derived from the extrusion external closure workbook.",
            },
            "denaturation_state": 0.85,
        }
        support_delta = build_matrix_experiment_support_delta_artifact(intake_payload)
        experiment_rows.append(
            {
                "experiment_id": intake_payload["experiment_id"],
                "arm_id": str(experiment.get("arm_id", "unknown")),
                "source_kind": source_kind,
                "intake_payload": intake_payload,
                "support_delta": support_delta,
                "direct_damage_markers": measured_damage,
                "direct_damage_complete": damage_complete,
                "acrylamide_ppb": measured_volatiles.get("acrylamide"),
            }
        )

    all_damage_complete = all(row["direct_damage_complete"] for row in experiment_rows)
    promotion_after = [bool(row["support_delta"]["promotion_assessment"]["promotion_ready_after"]) for row in experiment_rows]
    source_kinds = sorted({str(row.get("source_kind", "unknown")) for row in experiment_rows})
    if source_kinds == ["synthetic_diagnostic"]:
        status = "diagnostic_example_only"
    else:
        status = "ready_for_external_landing_review" if all_damage_complete else "incomplete_direct_damage_bundle"
    return {
        "summary": {
            "package_id": str(workbook.get("package_id", "extrusion_external_closure_workbook")),
            "selected_matrix": str(workbook.get("selected_matrix", "soy_iso")),
            "experiment_count": len(experiment_rows),
            "source_kinds": source_kinds,
            "all_direct_damage_markers_complete": all_damage_complete,
            "all_support_delta_payloads_promotion_ready": all(promotion_after),
            "status": status,
        },
        "experiments": experiment_rows,
    }


def render_extrusion_external_closure_execution_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    lines = [
        "# Extrusion External Closure Execution",
        "",
        f"Package id: {summary.get('package_id', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Experiments: {int(summary.get('experiment_count', 0))}",
        f"All direct damage markers complete: {bool(summary.get('all_direct_damage_markers_complete', False))}",
        f"All support-delta payloads promotion ready: {bool(summary.get('all_support_delta_payloads_promotion_ready', False))}",
        f"Status: {summary.get('status', 'unknown')}",
        "",
        "| Experiment | Arm | Damage Complete | Promotion Ready After | Blocker After |",
        "| --- | --- | --- | --- | --- |",
    ]
    for experiment in payload.get("experiments", []):
        assessment = experiment.get("support_delta", {}).get("promotion_assessment", {})
        lines.append(
            f"| {experiment.get('experiment_id', 'unknown')} | {experiment.get('arm_id', 'unknown')} | {bool(experiment.get('direct_damage_complete', False))} | {bool(assessment.get('promotion_ready_after', False))} | {assessment.get('promotion_blocker_after', 'unknown')} |"
        )
    return "\n".join(lines) + "\n"


def _reference_signal(benchmark_id: str, root: Path = ROOT) -> Dict[str, float]:
    benchmark = _load_benchmark(benchmark_id, root)
    signal_map = benchmark.get("measured_volatiles") or benchmark.get("reference_volatiles") or {}
    return {
        str(name): float(dict(row or {}).get("conc_ppb", 0.0) or 0.0)
        for name, row in signal_map.items()
    }


def _compute_follow_on_metrics(experiment: Mapping[str, Any], reference_signal: Mapping[str, float], process_profile: Mapping[str, Any]) -> Dict[str, float]:
    feed = _normalize_float_map(experiment.get("feed_reference_assays", {}))
    post = _normalize_float_map(experiment.get("post_extrusion_process_state_assays", {}))
    measured = _normalize_float_map(experiment.get("measured_volatiles_ppb", {}))
    mft_name = "2-methyl-3-furanthiol"
    pyrazine_name = "2,5-dimethylpyrazine"
    disulfide_name = "bis(2-methyl-3-furyl) disulfide"

    free_sh_retention = post["ellman_free_sh_umol_per_g"] / max(feed["pre_extrusion_free_sh_umol_per_g"], 1e-12)
    mft_reference = max(reference_signal.get(mft_name, 0.0), 1e-12)
    pyrazine_reference = max(reference_signal.get(pyrazine_name, 0.0), 1e-12)
    mft_recovery = measured[mft_name] / mft_reference
    pyrazine_recovery = measured[pyrazine_name] / pyrazine_reference
    return {
        "free_sh_retention_fraction": free_sh_retention,
        "disulfide_pressure_proxy": 1.0 - free_sh_retention,
        "sulfur_to_pyrazine_retention_ratio": mft_recovery / max(pyrazine_recovery, 1e-12),
        "furyl_disulfide_to_mft_ratio": measured[disulfide_name] / max(measured[mft_name], 1e-12),
        "furosine_mg_per_kg": post["furosine_mg_per_kg"],
        "mft_recovery_fraction": mft_recovery,
        "pyrazine_recovery_fraction": pyrazine_recovery,
        "predicted_mft_headspace_factor": float(compute_extrusion_headspace_adjustment(mft_name, process_profile).get("combined_headspace_factor", 1.0)),
        "predicted_pyrazine_headspace_factor": float(compute_extrusion_headspace_adjustment(pyrazine_name, process_profile).get("combined_headspace_factor", 1.0)),
    }


def build_extrusion_disulfide_follow_on_execution_artifact(
    workbook_or_path: Mapping[str, Any] | Path | str,
    *,
    root: Path = ROOT,
    allow_placeholders: bool = False,
) -> Dict[str, Any]:
    workbook = _load_yaml(workbook_or_path) if isinstance(workbook_or_path, (Path, str)) else dict(workbook_or_path)
    experiments = list(workbook.get("experiments", []))
    if len(experiments) < 2:
        raise ValueError("5.8 follow-on workbook must include at least two experiments so furosine gradient can be evaluated")

    rows = []
    for experiment in experiments:
        if not allow_placeholders and _contains_placeholder(experiment):
            raise ValueError(f"Workbook experiment {experiment.get('experiment_id', 'unknown')} still contains placeholder values")
        reference_benchmark_id = str(experiment.get("reference_benchmark_id", "")).strip()
        if not reference_benchmark_id:
            raise ValueError("5.8 follow-on workbook experiment must define reference_benchmark_id")
        reference_signal = _reference_signal(reference_benchmark_id, root)
        target_process = dict(experiment.get("target_process_settings", {}))
        process_profile = build_extrusion_process_profile(
            base_temperature_celsius=float(target_process.get("effective_temperature_celsius", 145.0) or 145.0),
            water_activity=0.75,
            protein_type=str(workbook.get("selected_matrix", "soy_iso")),
            sme_kj_per_kg=float(target_process.get("sme_kj_per_kg", 0.0) or 0.0),
            moisture_regime="hme",
        )
        metrics = _compute_follow_on_metrics(experiment, reference_signal, process_profile)
        rows.append(
            {
                "experiment_id": str(experiment.get("experiment_id", "unknown")),
                "arm_id": str(experiment.get("arm_id", "unknown")),
                "reference_benchmark_id": reference_benchmark_id,
                "metrics": metrics,
            }
        )

    rows.sort(key=lambda row: float(row["metrics"].get("furosine_mg_per_kg", 0.0)))
    furosine_gradient = rows[-1]["metrics"]["furosine_mg_per_kg"] - rows[0]["metrics"]["furosine_mg_per_kg"]
    source_kinds = sorted({str(experiment.get("source_kind", workbook.get("source_kind", "unknown"))) for experiment in experiments})
    return {
        "summary": {
            "package_id": str(workbook.get("package_id", "extrusion_disulfide_follow_on_workbook")),
            "selected_matrix": str(workbook.get("selected_matrix", "soy_iso")),
            "experiment_count": len(rows),
            "source_kinds": source_kinds,
            "furosine_damage_gradient": furosine_gradient,
            "status": "diagnostic_example_only" if source_kinds == ["synthetic_diagnostic"] else "ready_for_5_8_model_fitting_review",
        },
        "experiments": rows,
    }


def render_extrusion_disulfide_follow_on_execution_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    lines = [
        "# Extrusion 5.8 Follow-On Execution",
        "",
        f"Package id: {summary.get('package_id', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Experiments: {int(summary.get('experiment_count', 0))}",
        f"Furosine damage gradient: {float(summary.get('furosine_damage_gradient', 0.0)):.3f}",
        f"Status: {summary.get('status', 'unknown')}",
        "",
        "| Experiment | Arm | Free SH Retention | Disulfide Pressure | Sulfur/Pyrazine Ratio | Furyl Disulfide/MFT | MFT Predicted | Pyrazine Predicted |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for experiment in payload.get("experiments", []):
        metrics = dict(experiment.get("metrics", {}))
        lines.append(
            f"| {experiment.get('experiment_id', 'unknown')} | {experiment.get('arm_id', 'unknown')} | {float(metrics.get('free_sh_retention_fraction', 0.0)):.3f} | {float(metrics.get('disulfide_pressure_proxy', 0.0)):.3f} | {float(metrics.get('sulfur_to_pyrazine_retention_ratio', 0.0)):.3f} | {float(metrics.get('furyl_disulfide_to_mft_ratio', 0.0)):.3f} | {float(metrics.get('predicted_mft_headspace_factor', 0.0)):.3f} | {float(metrics.get('predicted_pyrazine_headspace_factor', 0.0)):.3f} |"
        )
    return "\n".join(lines) + "\n"


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(dict(payload), indent=2), encoding="utf-8")


def _write_yaml(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(yaml.safe_dump(dict(payload), sort_keys=False, allow_unicode=False), encoding="utf-8")


def export_extrusion_external_closure_execution(
    workbook_path: Path | str,
    output_dir: Path | str,
    *,
    root: Path = ROOT,
) -> Dict[str, Any]:
    workbook = _load_yaml(workbook_path)
    artifact = build_extrusion_external_closure_execution_artifact(workbook, root=root)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    stem = Path(workbook_path).stem
    _write_json(destination / f"{stem}_execution.json", artifact)
    (destination / f"{stem}_execution.md").write_text(render_extrusion_external_closure_execution_markdown(artifact), encoding="utf-8")
    for row in artifact.get("experiments", []):
        experiment_id = str(row.get("experiment_id", "unknown"))
        _write_yaml(destination / f"{experiment_id}_intake.yaml", row.get("intake_payload", {}))
        _write_json(destination / f"{experiment_id}_support_delta.json", row.get("support_delta", {}))
    return artifact


def export_extrusion_disulfide_follow_on_execution(
    workbook_path: Path | str,
    output_dir: Path | str,
    *,
    root: Path = ROOT,
) -> Dict[str, Any]:
    workbook = _load_yaml(workbook_path)
    artifact = build_extrusion_disulfide_follow_on_execution_artifact(workbook, root=root)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    stem = Path(workbook_path).stem
    _write_json(destination / f"{stem}_execution.json", artifact)
    (destination / f"{stem}_execution.md").write_text(render_extrusion_disulfide_follow_on_execution_markdown(artifact), encoding="utf-8")
    return artifact


def build_example_extrusion_external_closure_workbook(root: Path = ROOT) -> Dict[str, Any]:
    from src.extrusion_benchmark_landing import build_extrusion_external_closure_workbook

    workbook = build_extrusion_external_closure_workbook(root)
    workbook["package_id"] = "spi_extrusion_external_closure_diagnostic_example_2026"
    workbook["status"] = "synthetic_diagnostic_example"
    for index, experiment in enumerate(workbook.get("experiments", [])):
        experiment["source_kind"] = "synthetic_diagnostic"
        experiment["provenance"] = {
            "origin": "synthetic_diagnostic",
            "source_reference": f"synthetic_diagnostic_spi_extrusion_arm_{index + 1}",
            "source_doi": "synthetic_diagnostic",
            "measurement_date": "2026-04-08",
            "notes": "Diagnostic example generated from repo templates. Not real wet-lab data.",
        }
        experiment["measured_damage_markers"] = {
            "reactive_lysine_fraction": 0.73 - 0.05 * index,
            "furosine_mg_per_kg": 24.0 + 4.0 * index,
            "lysinoalanine_mg_per_kg": 88.0 + 9.0 * index,
        }
        experiment["measured_volatiles_ppb"] = {
            "2-methyl-3-furanthiol": 0.0038 - 0.0004 * index,
            "2-furfurylthiol": 0.0054 - 0.0005 * index,
            "bis(2-methyl-3-furyl) disulfide": 0.00052 + 0.00008 * index,
            "2,5-dimethylpyrazine": 0.00031 - 0.00002 * index,
            "Hexanal": 0.00018 + 0.00003 * index,
            "2-pentylfuran": 0.00011 + 0.00002 * index,
            "furfural": 0.21 + 0.03 * index,
            "acrylamide": 0.00007 + 0.00001 * index,
        }
    return workbook


def build_example_extrusion_disulfide_follow_on_workbook(root: Path = ROOT) -> Dict[str, Any]:
    from src.extrusion_benchmark_landing import build_extrusion_disulfide_follow_on_workbook

    workbook = build_extrusion_disulfide_follow_on_workbook(root)
    workbook["package_id"] = "spi_extrusion_5_8_diagnostic_example_2026"
    workbook["status"] = "synthetic_diagnostic_example"
    workbook["source_kind"] = "synthetic_diagnostic"
    profiles = [
        {
            "feed_reference_assays": {
                "pre_extrusion_free_sh_umol_per_g": 12.0,
                "pre_extrusion_free_amino_groups_umol_per_g": 18.0,
            },
            "post_extrusion_process_state_assays": {
                "ellman_free_sh_umol_per_g": 8.4,
                "opa_free_amino_groups_umol_per_g": 14.4,
                "furosine_mg_per_kg": 24.0,
                "post_extrusion_ph": 5.7,
            },
            "measured_volatiles_ppb": {
                "2-methyl-3-furanthiol": 0.0038,
                "2-furfurylthiol": 0.0053,
                "bis(2-methyl-3-furyl) disulfide": 0.00052,
                "2,5-dimethylpyrazine": 0.00031,
                "Hexanal": 0.00018,
                "2-pentylfuran": 0.00011,
            },
        },
        {
            "feed_reference_assays": {
                "pre_extrusion_free_sh_umol_per_g": 12.0,
                "pre_extrusion_free_amino_groups_umol_per_g": 18.0,
            },
            "post_extrusion_process_state_assays": {
                "ellman_free_sh_umol_per_g": 6.6,
                "opa_free_amino_groups_umol_per_g": 13.5,
                "furosine_mg_per_kg": 31.0,
                "post_extrusion_ph": 5.6,
            },
            "measured_volatiles_ppb": {
                "2-methyl-3-furanthiol": 0.0031,
                "2-furfurylthiol": 0.0048,
                "bis(2-methyl-3-furyl) disulfide": 0.00061,
                "2,5-dimethylpyrazine": 0.00027,
                "Hexanal": 0.00021,
                "2-pentylfuran": 0.00013,
            },
        },
    ]
    for experiment, profile in zip(workbook.get("experiments", []), profiles):
        experiment["source_kind"] = "synthetic_diagnostic"
        experiment.update(profile)
    return workbook


def export_extrusion_diagnostic_examples_bundle(output_dir: Path | str, *, root: Path = ROOT) -> Dict[str, Any]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)

    closure_workbook = build_example_extrusion_external_closure_workbook(root)
    closure_yaml = destination / "extrusion_external_closure_diagnostic_example.yaml"
    _write_yaml(closure_yaml, closure_workbook)
    closure_execution = export_extrusion_external_closure_execution(closure_yaml, destination, root=root)

    follow_on_workbook = build_example_extrusion_disulfide_follow_on_workbook(root)
    follow_on_yaml = destination / "extrusion_disulfide_follow_on_diagnostic_example.yaml"
    _write_yaml(follow_on_yaml, follow_on_workbook)
    follow_on_execution = export_extrusion_disulfide_follow_on_execution(follow_on_yaml, destination, root=root)

    bundle = {
        "summary": {
            "status": "diagnostic_examples_generated",
            "output_dir": str(destination),
        },
        "artifacts": {
            "closure_workbook": str(closure_yaml),
            "follow_on_workbook": str(follow_on_yaml),
            "closure_execution_status": closure_execution.get("summary", {}).get("status", "unknown"),
            "follow_on_execution_status": follow_on_execution.get("summary", {}).get("status", "unknown"),
        },
    }
    _write_json(destination / "extrusion_diagnostic_examples_bundle.json", bundle)
    return bundle