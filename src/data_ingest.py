from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

import pandas as pd
import yaml

from src import data_paths
from src.benchmark_validation import get_benchmark_files
from src.matrix_experiment_intake import (
    build_matrix_experiment_benchmark_payload,
    build_matrix_experiment_support_delta_artifact,
    load_matrix_experiment_intake,
    normalize_matrix_experiment_intake,
)
from src.presentation import render_matrix_experiment_support_delta_markdown
from src.uncertainty_propagation import propagate_benchmarks


DEFAULT_OUTPUT_DIR = data_paths.VALIDATION_DIR
#: Preview runs (no --confirm) must not write into the tracked validation
#: directory. Before 2026-08-27 a plain preview dropped four files into
#: results/validation/ uninvited, which looked like validated artifacts.
DEFAULT_PREVIEW_OUTPUT_DIR = data_paths.INGEST_PREVIEWS_DIR


@dataclass(frozen=True)
class IngestArtifacts:
    preview_json: Path
    preview_md: Path
    support_delta_json: Path
    support_delta_md: Path
    intake_yaml: Optional[Path]


_COLUMN_ALIASES: Dict[str, Sequence[str]] = {
    "compound": (
        "compound",
        "compound name",
        "analyte",
        "analyte name",
        "volatile",
        "volatile name",
        "marker",
        "name",
    ),
    "conc_ppb": (
        "conc ppb",
        "concentration ppb",
        "measured ppb",
        "observed ppb",
        "ppb",
        "value ppb",
        "headspace ppb",
    ),
    "uncertainty_pct": (
        "uncertainty pct",
        "uncertainty %",
        "rsd %",
        "cv %",
        "error %",
        "relative uncertainty %",
    ),
    "temp_C": ("temp c", "temperature c", "temperature", "temp", "temperature deg c"),
    "ph": ("ph",),
    "water_activity": ("water activity", "aw", "a_w"),
    "time_min": ("time min", "time minutes", "minutes", "duration min", "time"),
    "source_reference": ("source reference", "reference", "citation", "source"),
    "source_doi": ("source doi", "doi"),
    "measurement_date": ("measurement date", "date", "run date"),
    "notes": ("notes", "comment", "comments"),
    "target_benchmark_id": ("target benchmark id", "benchmark id", "aligned benchmark id"),
}


def _normalize_header(value: str) -> str:
    token = "".join(ch.lower() if ch.isalnum() else " " for ch in str(value))
    return " ".join(token.split())


def _alias_score(header: str, alias: str) -> float:
    normalized_header = _normalize_header(header)
    normalized_alias = _normalize_header(alias)
    if not normalized_header or not normalized_alias:
        return 0.0
    if normalized_header == normalized_alias:
        return 1.0
    if normalized_alias in normalized_header or normalized_header in normalized_alias:
        return 0.94
    header_tokens = set(normalized_header.split())
    alias_tokens = set(normalized_alias.split())
    if header_tokens and alias_tokens:
        overlap = len(header_tokens.intersection(alias_tokens)) / len(alias_tokens)
        if overlap >= 0.75:
            return 0.9
        if overlap >= 0.5:
            return 0.78
    return 0.0


def _best_column(columns: Sequence[str], key: str) -> Optional[str]:
    best: tuple[float, Optional[str]] = (0.0, None)
    for column in columns:
        for alias in _COLUMN_ALIASES.get(key, ()):
            score = _alias_score(column, alias)
            if score > best[0]:
                best = (score, str(column))
    return best[1] if best[0] >= 0.75 else None


def _read_tabular_rows(path: Path, *, sheet_name: Optional[str] = None) -> List[Dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        frame = pd.read_csv(path)
    elif suffix in {".xlsx", ".xls"}:
        frame = pd.read_excel(path, sheet_name=sheet_name or 0)
    else:
        raise ValueError(f"Unsupported tabular ingest format: {path.suffix}")
    frame = frame.where(pd.notnull(frame), None)
    return frame.to_dict(orient="records")


def _load_source(path: Path, *, sheet_name: Optional[str] = None) -> Any:
    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        text = path.read_text(encoding="utf-8")
        return yaml.safe_load(text)
    if suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    return _read_tabular_rows(path, sheet_name=sheet_name)


def _parse_precursor_args(values: Sequence[str]) -> Dict[str, Dict[str, float]]:
    precursors: Dict[str, Dict[str, float]] = {}
    for item in values:
        name, sep, raw_value = str(item).partition("=")
        if not sep:
            raise ValueError(f"Invalid --precursor value '{item}'. Use NAME=CONCENTRATION_MILLIMOLAR.")
        precursors[name.strip()] = {"concentration_mM": float(raw_value)}
    return precursors


def _resolve_scalar(
    cli_value: Any,
    file_value: Any,
    *,
    required: bool = False,
    default: Any = None,
    field_name: str,
) -> Any:
    value = cli_value if cli_value not in {None, ""} else file_value
    if value in {None, ""}:
        value = default
    if required and value in {None, ""}:
        raise ValueError(f"Missing required ingest field: {field_name}")
    return value


def _metadata_from_rows(rows: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    if not rows:
        return {}
    columns = [str(column) for column in rows[0].keys()]
    mapped: Dict[str, str] = {}
    for key in _COLUMN_ALIASES:
        column = _best_column(columns, key)
        if column is not None:
            mapped[key] = column
    first = rows[0]
    metadata: Dict[str, Any] = {}
    for key, column in mapped.items():
        metadata[key] = first.get(column)
    metadata["_mapped_columns"] = mapped
    return metadata


def _build_payload_from_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    file_path: Path,
    protein_type: Optional[str],
    process_state: Optional[str],
    experiment_id: Optional[str],
    source_kind: str,
    evidence_class: str,
    temp_c: Optional[float],
    ph: Optional[float],
    water_activity: Optional[float],
    time_min: Optional[float],
    matrix_format: Optional[str],
    source_reference: Optional[str],
    source_doi: Optional[str],
    measurement_date: Optional[str],
    notes: Optional[str],
    target_benchmark_id: Optional[str],
    precursors: Mapping[str, Mapping[str, float]],
) -> Dict[str, Any]:
    if not rows:
        raise ValueError("Ingest file does not contain any rows")

    metadata = _metadata_from_rows(rows)
    mapped = metadata.get("_mapped_columns", {})
    compound_column = mapped.get("compound")
    conc_column = mapped.get("conc_ppb")
    if compound_column is None or conc_column is None:
        raise ValueError("Could not map the compound and ppb columns from the ingest file")

    uncertainty_column = mapped.get("uncertainty_pct")
    measured_volatiles: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        compound = str(row.get(compound_column, "") or "").strip()
        if not compound:
            continue
        raw_value = row.get(conc_column)
        if raw_value in {None, ""}:
            continue
        measured_entry: Dict[str, Any] = {"conc_ppb": float(raw_value)}
        if uncertainty_column is not None and row.get(uncertainty_column) not in {None, ""}:
            measured_entry["uncertainty_pct"] = float(row.get(uncertainty_column))
        measured_volatiles[compound] = measured_entry
    if not measured_volatiles:
        raise ValueError("No measured volatiles could be mapped from the ingest rows")

    resolved_precursors = dict(precursors)
    if not resolved_precursors:
        raise ValueError("At least one --precursor NAME=CONCENTRATION_MILLIMOLAR entry is required for tabular ingest")

    aligned_benchmark_id = _resolve_scalar(
        target_benchmark_id,
        metadata.get("target_benchmark_id"),
        required=False,
        default=None,
        field_name="target_benchmark_id",
    )
    payload: Dict[str, Any] = {
        "experiment_id": _resolve_scalar(
            experiment_id,
            None,
            required=False,
            default=file_path.stem,
            field_name="experiment_id",
        ),
        "source_kind": source_kind,
        "evidence_class": evidence_class,
        "protein_type": _resolve_scalar(
            protein_type,
            None,
            required=True,
            field_name="protein_type",
        ),
        "process_state": _resolve_scalar(
            process_state,
            None,
            required=True,
            field_name="process_state",
        ),
        "matrix_format": _resolve_scalar(
            matrix_format,
            None,
            required=False,
            default="instrument_import",
            field_name="matrix_format",
        ),
        "conditions": {
            "temp_C": float(_resolve_scalar(temp_c, metadata.get("temp_C"), required=True, field_name="temp_C")),
            "ph": float(_resolve_scalar(ph, metadata.get("ph"), required=True, field_name="ph")),
            "water_activity": float(
                _resolve_scalar(water_activity, metadata.get("water_activity"), required=True, field_name="water_activity")
            ),
            "time_min": float(_resolve_scalar(time_min, metadata.get("time_min"), required=True, field_name="time_min")),
        },
        "formulation": {"precursors": resolved_precursors},
        "measured_volatiles": measured_volatiles,
        "provenance": {
            "origin": source_kind,
            "source_reference": _resolve_scalar(
                source_reference,
                metadata.get("source_reference"),
                required=True,
                field_name="source_reference",
            ),
            "source_doi": _resolve_scalar(
                source_doi,
                metadata.get("source_doi"),
                required=False,
                default=None,
                field_name="source_doi",
            ),
            "measurement_date": _resolve_scalar(
                measurement_date,
                metadata.get("measurement_date"),
                required=False,
                default="unspecified",
                field_name="measurement_date",
            ),
            "notes": _resolve_scalar(
                notes,
                metadata.get("notes"),
                required=False,
                default="Imported from scientist-facing ingest CLI.",
                field_name="notes",
            ),
        },
    }
    if aligned_benchmark_id:
        payload["benchmark_alignment"] = {"target_benchmark_id": aligned_benchmark_id}
    return payload


def build_intake_payload(
    file_path: Path | str,
    *,
    protein_type: Optional[str] = None,
    process_state: Optional[str] = None,
    experiment_id: Optional[str] = None,
    source_kind: str = "internal_experiment",
    evidence_class: str = "calibration_candidate",
    temp_c: Optional[float] = None,
    ph: Optional[float] = None,
    water_activity: Optional[float] = None,
    time_min: Optional[float] = None,
    matrix_format: Optional[str] = None,
    source_reference: Optional[str] = None,
    source_doi: Optional[str] = None,
    measurement_date: Optional[str] = None,
    notes: Optional[str] = None,
    target_benchmark_id: Optional[str] = None,
    precursor_entries: Sequence[str] = (),
    sheet_name: Optional[str] = None,
) -> Dict[str, Any]:
    path = Path(file_path)
    payload_or_rows = _load_source(path, sheet_name=sheet_name)
    if isinstance(payload_or_rows, dict) and payload_or_rows.get("measured_volatiles"):
        payload = dict(payload_or_rows)
        if protein_type:
            payload["protein_type"] = protein_type
        if process_state:
            payload["process_state"] = process_state
        if experiment_id:
            payload["experiment_id"] = experiment_id
        if source_kind:
            payload["source_kind"] = source_kind
        if evidence_class:
            payload["evidence_class"] = evidence_class
        if matrix_format:
            payload["matrix_format"] = matrix_format
        payload.setdefault("provenance", {})
        if source_reference:
            payload["provenance"]["source_reference"] = source_reference
        if source_doi is not None:
            payload["provenance"]["source_doi"] = source_doi
        if measurement_date:
            payload["provenance"]["measurement_date"] = measurement_date
        if notes:
            payload["provenance"]["notes"] = notes
        if target_benchmark_id:
            payload.setdefault("benchmark_alignment", {})["target_benchmark_id"] = target_benchmark_id
        if precursor_entries:
            payload.setdefault("formulation", {})["precursors"] = _parse_precursor_args(precursor_entries)
        return normalize_matrix_experiment_intake(payload)

    if isinstance(payload_or_rows, list):
        payload = _build_payload_from_rows(
            payload_or_rows,
            file_path=path,
            protein_type=protein_type,
            process_state=process_state,
            experiment_id=experiment_id,
            source_kind=source_kind,
            evidence_class=evidence_class,
            temp_c=temp_c,
            ph=ph,
            water_activity=water_activity,
            time_min=time_min,
            matrix_format=matrix_format,
            source_reference=source_reference,
            source_doi=source_doi,
            measurement_date=measurement_date,
            notes=notes,
            target_benchmark_id=target_benchmark_id,
            precursors=_parse_precursor_args(precursor_entries),
        )
        return normalize_matrix_experiment_intake(payload)

    raise ValueError("Unsupported ingest payload structure")


def _coverage_tuple(payload: Mapping[str, Any]) -> tuple[Optional[int], Optional[int], Optional[float]]:
    summary = payload.get("summary", {})
    hits = summary.get("ci_coverage_hits")
    matched = summary.get("matched_compound_count")
    rate = summary.get("ci_coverage_rate")
    return (
        int(hits) if hits is not None else None,
        int(matched) if matched is not None else None,
        float(rate) if rate is not None else None,
    )


def _uncertainty_lookup(payload: Mapping[str, Any]) -> Dict[str, Dict[str, Mapping[str, Any]]]:
    lookup: Dict[str, Dict[str, Mapping[str, Any]]] = {}
    for benchmark in payload.get("benchmarks", []) or []:
        benchmark_id = str(benchmark.get("benchmark_id", ""))
        if not benchmark_id:
            continue
        lookup[benchmark_id] = {
            str(row.get("compound", "")): row
            for row in benchmark.get("compounds", []) or []
            if str(row.get("compound", "")).strip()
        }
    return lookup


def _ci_width_log10(row: Mapping[str, Any]) -> Optional[float]:
    p5 = float(row.get("predicted_p5", 0.0) or 0.0)
    p95 = float(row.get("predicted_p95", 0.0) or 0.0)
    if p5 <= 0.0 or p95 <= 0.0:
        return None
    return math.log10(p95 / p5)


def _median_shift_dex(baseline_row: Mapping[str, Any], candidate_row: Mapping[str, Any]) -> Optional[float]:
    baseline = float(baseline_row.get("predicted_p50", 0.0) or 0.0)
    candidate = float(candidate_row.get("predicted_p50", 0.0) or 0.0)
    if baseline <= 0.0 or candidate <= 0.0:
        return None
    return math.log10(candidate / baseline)


def build_ingest_preview(
    payload: Mapping[str, Any],
    *,
    preview_samples: int = 80,
    preview_seed: int = 0,
) -> Dict[str, Any]:
    normalized = normalize_matrix_experiment_intake(payload)
    support_delta = build_matrix_experiment_support_delta_artifact(normalized)
    candidate_benchmark = build_matrix_experiment_benchmark_payload(normalized)
    benchmark_id = str(candidate_benchmark.get("benchmark_id", normalized.get("experiment_id", "ingest_payload")))
    benchmark_files = list(get_benchmark_files())

    baseline_uncertainty = propagate_benchmarks(
        benchmark_files=benchmark_files,
        n_samples=preview_samples,
        seed=preview_seed,
    )

    coverage_preview: Dict[str, Any]
    if str(candidate_benchmark.get("metadata", {}).get("execution_path", "")) == "matrix_only":
        coverage_preview = {
            "participates_in_trust_loop": False,
            "reason": "matrix_only benchmarks are outside the current trust-loop headline surface",
            "before": None,
            "after": None,
        }
        candidate_uncertainty = baseline_uncertainty
    else:
        with TemporaryDirectory(prefix="maillard_ingest_preview_") as tmp_dir:
            candidate_path = Path(tmp_dir) / f"{benchmark_id}.json"
            candidate_path.write_text(json.dumps(candidate_benchmark, indent=2, sort_keys=True), encoding="utf-8")
            candidate_uncertainty = propagate_benchmarks(
                benchmark_files=[*benchmark_files, candidate_path],
                n_samples=preview_samples,
                seed=preview_seed,
            )
        before_hits, before_total, before_rate = _coverage_tuple(baseline_uncertainty)
        after_hits, after_total, after_rate = _coverage_tuple(candidate_uncertainty)
        coverage_preview = {
            "participates_in_trust_loop": True,
            "before": {
                "inside_ci_count": before_hits,
                "matched_compound_count": before_total,
                "coverage_rate": before_rate,
            },
            "after": {
                "inside_ci_count": after_hits,
                "matched_compound_count": after_total,
                "coverage_rate": after_rate,
            },
            "delta": {
                "inside_ci_count": (after_hits - before_hits) if None not in {before_hits, after_hits} else None,
                "matched_compound_count": (after_total - before_total) if None not in {before_total, after_total} else None,
                "coverage_rate": (after_rate - before_rate) if None not in {before_rate, after_rate} else None,
            },
        }

    baseline_lookup = _uncertainty_lookup(baseline_uncertainty)
    candidate_lookup = _uncertainty_lookup(candidate_uncertainty)
    aligned_benchmark_id = str(support_delta.get("aligned_benchmark", {}).get("benchmark_id", "") or "")
    aligned_rows = baseline_lookup.get(aligned_benchmark_id, {})
    candidate_rows = candidate_lookup.get(benchmark_id, {})

    comparison_rows: List[Dict[str, Any]] = []
    for compound, candidate_row in candidate_rows.items():
        baseline_row = aligned_rows.get(compound)
        baseline_width = _ci_width_log10(baseline_row) if baseline_row is not None else None
        candidate_width = _ci_width_log10(candidate_row)
        width_delta = None
        if baseline_width is not None and candidate_width is not None:
            width_delta = baseline_width - candidate_width
        median_shift = _median_shift_dex(baseline_row, candidate_row) if baseline_row is not None else None
        comparison_rows.append(
            {
                "compound": compound,
                "aligned_benchmark_id": aligned_benchmark_id or None,
                "baseline_ci_width_log10": baseline_width,
                "candidate_ci_width_log10": candidate_width,
                "ci_width_delta_log10": width_delta,
                "candidate_predicted_p50": float(candidate_row.get("predicted_p50", 0.0) or 0.0),
                "baseline_predicted_p50": float(baseline_row.get("predicted_p50", 0.0) or 0.0) if baseline_row is not None else None,
                "median_shift_dex": median_shift,
            }
        )

    top_tightening = [
        row for row in sorted(
            comparison_rows,
            key=lambda row: float(row.get("ci_width_delta_log10", float("-inf")) or float("-inf")),
            reverse=True,
        )
        if row.get("ci_width_delta_log10") is not None and float(row["ci_width_delta_log10"]) > 0.0
    ][:5]
    top_median_shifts = [
        row for row in sorted(
            comparison_rows,
            key=lambda row: abs(float(row.get("median_shift_dex", 0.0) or 0.0)),
            reverse=True,
        )
        if row.get("median_shift_dex") is not None and abs(float(row["median_shift_dex"])) >= 0.3
    ][:5]

    return {
        "experiment_id": normalized.get("experiment_id"),
        "protein_type": normalized.get("protein_type"),
        "process_state": normalized.get("process_state"),
        "benchmark_preview": {
            "candidate_benchmark_id": benchmark_id,
            "benchmarks_added": [benchmark_id],
            "aligned_benchmark_id": aligned_benchmark_id or None,
        },
        "coverage_preview": coverage_preview,
        "promotion_preview": support_delta.get("promotion_assessment", {}),
        "support_delta": support_delta.get("support_delta", {}),
        "compound_deltas": {
            "top_strengthened": [
                row for row in support_delta.get("compounds", []) if row.get("support_delta") == "strengthened"
            ][:5],
            "top_weakened": [
                row for row in support_delta.get("compounds", []) if row.get("support_delta") == "weakened"
            ][:5],
            "top_envelope_tightening_vs_aligned": top_tightening,
            "top_median_shift_vs_aligned": top_median_shifts,
        },
        "audit": {
            "preview_samples": preview_samples,
            "preview_seed": preview_seed,
            "baseline_uncertainty_summary": baseline_uncertainty.get("summary", {}),
            "candidate_uncertainty_summary": candidate_uncertainty.get("summary", {}),
        },
    }


def render_ingest_preview_markdown(preview: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Ingest Preview",
        "",
        f"- Experiment: {preview.get('experiment_id', 'unknown')}",
        f"- Protein type: {preview.get('protein_type', 'unknown')}",
        f"- Process state: {preview.get('process_state', 'unknown')}",
    ]
    benchmark_preview = preview.get("benchmark_preview", {})
    lines.append(f"- Candidate benchmark: {benchmark_preview.get('candidate_benchmark_id', 'unknown')}")
    if benchmark_preview.get("aligned_benchmark_id"):
        lines.append(f"- Closest aligned benchmark: {benchmark_preview.get('aligned_benchmark_id')}")
    coverage_preview = preview.get("coverage_preview", {})
    if coverage_preview.get("participates_in_trust_loop"):
        before = coverage_preview.get("before", {})
        after = coverage_preview.get("after", {})
        lines.append(
            f"- Trust-loop coverage: {before.get('inside_ci_count', 0)}/{before.get('matched_compound_count', 0)} -> {after.get('inside_ci_count', 0)}/{after.get('matched_compound_count', 0)}"
        )
    else:
        lines.append(f"- Trust-loop coverage: n/a ({coverage_preview.get('reason', 'not applicable')})")

    promotion = preview.get("promotion_preview", {})
    lines.extend(
        [
            "",
            "## Promotion Preview",
            "",
            f"- Ready before: {promotion.get('promotion_ready_before', False)}",
            f"- Ready after: {promotion.get('promotion_ready_after', False)}",
            f"- Readiness change: {promotion.get('readiness_change', 'unknown')}",
            f"- Landing recommendation: {promotion.get('landing_recommendation', 'unknown')}",
        ]
    )

    compound_deltas = preview.get("compound_deltas", {})
    sections = [
        ("Top Strengthened Compounds", compound_deltas.get("top_strengthened", []), "support_delta"),
        ("Top Weakened Compounds", compound_deltas.get("top_weakened", []), "support_delta"),
        ("Top Envelope Tightening Vs Aligned Benchmark", compound_deltas.get("top_envelope_tightening_vs_aligned", []), "ci_width_delta_log10"),
        ("Top Median Shift Vs Aligned Benchmark", compound_deltas.get("top_median_shift_vs_aligned", []), "median_shift_dex"),
    ]
    for heading, rows, key in sections:
        lines.extend(["", f"## {heading}", ""])
        if not rows:
            lines.append("- None")
            continue
        for row in rows:
            if key == "support_delta":
                lines.append(f"- {row.get('compound', 'unknown')}: {row.get('support_delta', 'unknown')}")
            else:
                lines.append(f"- {row.get('compound', 'unknown')}: {float(row.get(key, 0.0) or 0.0):+.2f}")
    lines.append("")
    return "\n".join(lines)


def write_ingest_artifacts(
    payload: Mapping[str, Any],
    preview: Mapping[str, Any],
    *,
    output_dir: Path | str = DEFAULT_OUTPUT_DIR,
    persist_intake: bool,
) -> IngestArtifacts:
    normalized = normalize_matrix_experiment_intake(payload)
    experiment_id = str(normalized.get("experiment_id", "matrix_ingest"))
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)

    preview_json = destination / f"{experiment_id}_ingest_preview.json"
    preview_md = destination / f"{experiment_id}_ingest_preview.md"
    preview_json.write_text(json.dumps(preview, indent=2, sort_keys=True), encoding="utf-8")
    preview_md.write_text(render_ingest_preview_markdown(preview), encoding="utf-8")

    support_delta = build_matrix_experiment_support_delta_artifact(normalized)
    support_delta_json = destination / f"{experiment_id}_support_delta.json"
    support_delta_md = destination / f"{experiment_id}_support_delta.md"
    support_delta_json.write_text(json.dumps(support_delta, indent=2, sort_keys=True), encoding="utf-8")
    support_delta_md.write_text(render_matrix_experiment_support_delta_markdown(support_delta), encoding="utf-8")

    intake_yaml: Optional[Path] = None
    if persist_intake:
        intake_yaml = destination / f"{experiment_id}.yaml"
        intake_yaml.write_text(yaml.safe_dump(dict(normalized), sort_keys=False), encoding="utf-8")

    return IngestArtifacts(
        preview_json=preview_json,
        preview_md=preview_md,
        support_delta_json=support_delta_json,
        support_delta_md=support_delta_md,
        intake_yaml=intake_yaml,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Normalize a scientist-facing experiment file into the matrix intake contract.")
    parser.add_argument("--file", required=True, help="Input CSV, Excel, YAML, or JSON file.")
    parser.add_argument("--protein-type", choices=["pea_iso", "soy_iso", "wheat_gluten", "myco", "free"], help="Protein matrix type.")
    parser.add_argument("--process-state", help="Matrix process state, e.g. extrusion_structured.")
    parser.add_argument("--experiment-id", help="Stable ID for the canonical intake payload.")
    parser.add_argument("--source-kind", default="internal_experiment", choices=["external_literature", "internal_experiment", "synthetic_diagnostic"])
    parser.add_argument("--evidence-class", default="calibration_candidate", choices=["calibration_candidate", "external_validation_only", "diagnostic_only"])
    parser.add_argument("--temp-c", type=float)
    parser.add_argument("--ph", type=float)
    parser.add_argument("--water-activity", type=float)
    parser.add_argument("--time-min", type=float)
    parser.add_argument("--matrix-format")
    parser.add_argument("--source-reference")
    parser.add_argument("--source-doi")
    parser.add_argument("--measurement-date")
    parser.add_argument("--notes")
    parser.add_argument("--target-benchmark-id")
    parser.add_argument("--precursor", action="append", default=[], help="Repeatable NAME=CONCENTRATION_MILLIMOLAR pair.")
    parser.add_argument("--sheet-name", help="Excel sheet name to ingest. Defaults to the first sheet.")
    parser.add_argument("--preview-samples", type=int, default=80)
    parser.add_argument("--preview-seed", type=int, default=0)
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Where to write the artifacts. Defaults to "
            f"'{data_paths.rel(DEFAULT_PREVIEW_OUTPUT_DIR)}' in preview mode and "
            f"'{data_paths.rel(DEFAULT_OUTPUT_DIR)}' with --confirm."
        ),
    )
    parser.add_argument(
        "--confirm",
        action="store_true",
        help=(
            "Persist the canonical intake YAML (in addition to the preview and "
            "support-delta artifacts) into the validation directory. This writes ONE "
            "YAML file; it does not rebuild the benchmark panel or regenerate "
            "validation artifacts — run the generators for that."
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    payload = build_intake_payload(
        args.file,
        protein_type=args.protein_type,
        process_state=args.process_state,
        experiment_id=args.experiment_id,
        source_kind=args.source_kind,
        evidence_class=args.evidence_class,
        temp_c=args.temp_c,
        ph=args.ph,
        water_activity=args.water_activity,
        time_min=args.time_min,
        matrix_format=args.matrix_format,
        source_reference=args.source_reference,
        source_doi=args.source_doi,
        measurement_date=args.measurement_date,
        notes=args.notes,
        target_benchmark_id=args.target_benchmark_id,
        precursor_entries=args.precursor,
        sheet_name=args.sheet_name,
    )
    preview = build_ingest_preview(
        payload,
        preview_samples=args.preview_samples,
        preview_seed=args.preview_seed,
    )
    if args.output_dir:
        output_dir: Path | str = args.output_dir
    else:
        output_dir = DEFAULT_OUTPUT_DIR if args.confirm else DEFAULT_PREVIEW_OUTPUT_DIR

    artifacts = write_ingest_artifacts(
        payload,
        preview,
        output_dir=output_dir,
        persist_intake=bool(args.confirm),
    )

    print(f"Experiment: {preview.get('experiment_id')}")
    print(f"Candidate benchmark: {preview.get('benchmark_preview', {}).get('candidate_benchmark_id')}")
    promotion = preview.get("promotion_preview", {})
    print(
        "Promotion readiness: "
        f"{promotion.get('promotion_ready_before', False)} -> {promotion.get('promotion_ready_after', False)} "
        f"({promotion.get('readiness_change', 'unknown')})"
    )
    coverage = preview.get("coverage_preview", {})
    if coverage.get("participates_in_trust_loop"):
        before = coverage.get("before", {})
        after = coverage.get("after", {})
        print(
            "Trust-loop coverage: "
            f"{before.get('inside_ci_count', 0)}/{before.get('matched_compound_count', 0)} -> "
            f"{after.get('inside_ci_count', 0)}/{after.get('matched_compound_count', 0)}"
        )
    else:
        print(f"Trust-loop coverage: n/a ({coverage.get('reason', 'not applicable')})")
    print(f"Preview JSON: {artifacts.preview_json}")
    print(f"Preview Markdown: {artifacts.preview_md}")
    print(f"Support delta JSON: {artifacts.support_delta_json}")
    print(f"Support delta Markdown: {artifacts.support_delta_md}")
    if artifacts.intake_yaml is not None:
        print(f"Canonical intake YAML: {artifacts.intake_yaml}")
        print(
            "Note: --confirm wrote that one YAML file. It did NOT rebuild the benchmark "
            "panel or regenerate validation artifacts — run the generators in "
            "scripts/generators/ (or ./scripts/docker_maillard.sh summary) for that."
        )
    else:
        print(
            "Preview only. No canonical intake YAML was written and nothing under "
            f"{data_paths.rel(DEFAULT_OUTPUT_DIR)}/ was touched. Re-run with --confirm to persist the intake YAML."
        )
    return 0
