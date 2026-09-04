from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping

from src import data_paths
from src import data_access
from src.literature_intake_registry import build_literature_backlog_artifact


_PRIORITY_WEIGHTS = {
    "critical": 6.0,
    "high": 4.0,
    "medium": 2.5,
    "low": 1.0,
}

_WORKSTREAMS = [
    {
        "id": "primary_matrix_benchmark_package",
        "title": "Primary PPI/SPI benchmark package",
        "phase_order": 1,
        "intake_gap_ids": [
            "ppi_meaty_positive_matrix_benchmark",
            "spi_meaty_positive_matrix_benchmark",
            "ppi_spi_time_series",
            "meaty_off_flavour_safety_tradeoff_panel",
        ],
        "process_gap_ids": [
            "intact_spi_ppi_quantified_mft_fft",
            "pipi_spi_pyrazine_absolute_quantitation",
            "aqueous_pipi_spi_acrylamide_kinetics",
            "ellman_opa_dsc_same_experiment",
        ],
        "repo_ready_inputs": [
            "minimum_primary_experiment contract already encoded",
            "PPI/SPI primary benchmark protocol already documented",
            "pea and soy benchmark case-study notes already documented",
            "matrix experiment intake contract already available",
        ],
        "deliverables": [
            "shareable wet-lab request for benchmark-grade PPI and SPI runs",
            "ingestion-ready contract for same-run meaty, off-note, and safety markers",
            "time-series data package that can promote matrix lanes without more literature hunting",
        ],
        "why_now": "No curated runtime-ready or benchmark-ready citations remain, so the next leverage is the existing wet-lab contract rather than more reading.",
    },
    {
        "id": "extrusion_benchmark_translation",
        "title": "Minimum viable extrusion benchmark",
        "phase_order": 2,
        "intake_gap_ids": [
            "meaty_off_flavour_safety_tradeoff_panel",
            "mft_fft_matrix_retention",
        ],
        "process_gap_ids": [
            "intact_spi_ppi_quantified_mft_fft",
            "aqueous_pipi_spi_acrylamide_kinetics",
            "ellman_opa_dsc_same_experiment",
        ],
        "repo_ready_inputs": [
            "extrusion process model already exposes SME, moisture regime, and sequential zone profiles",
            "furosine and LAL pre-extrusion baselines already encoded in the extrusion layer",
            "soy high-temperature safety instrumentation is already traceable in safety references",
        ],
        "deliverables": [
            "one-protein, two-SME benchmark request with MFT, hexanal, and furosine in the same run",
            "results-facing protocol artifact that can be handed to a lab team",
            "bridge from aqueous benchmark closure into extrusion confidence claims",
        ],
        "why_now": "Extrusion is the dominant scientist-facing process gap, but it should land immediately after the primary matrix package so the first extrudate campaign inherits the same analyte contract.",
    },
    {
        "id": "retention_follow_on_closure",
        "title": "Retention and release follow-on",
        "phase_order": 3,
        "intake_gap_ids": [
            "mft_fft_matrix_retention",
        ],
        "process_gap_ids": [
            "ellman_opa_dsc_same_experiment",
        ],
        "repo_ready_inputs": [
            "bounded computational closure candidate already declared in the intake registry",
            "matrix retention priors already exist for directional reporting",
        ],
        "deliverables": [
            "bounded retention update once first matrix or extrusion runs exist",
            "calibrated follow-on for die-loss or release modeling instead of a free-standing solver",
        ],
        "why_now": "This is the right third step because retention needs the benchmark panel above to avoid fitting against proxies only.",
    },
]

_DEFERRED_MODELING_REVIEW = [
    {
        "item_id": "5.7",
        "title": "Bidirectional Lipid-Maillard Crosstalk",
        "decision": "defer",
        "reason": "Improves tradeoff realism, but the first unlock still depends more on same-run benchmark data than on a broader coupled solver.",
        "revisit_after": "after primary matrix and extrusion benchmark panels exist",
    },
    {
        "item_id": "5.8",
        "title": "Disulfide Bond Evolution / MFT Retention",
        "decision": "next_after_s17_baseline",
        "reason": "This is the only deferred modeling item directly adjacent to extrusion, but it should follow the first benchmark protocol because the required SH and damage observables can be captured in that run.",
        "revisit_after": "after the first SPI extrusion benchmark produces Ellman, furosine, and volatile data",
    },
    {
        "item_id": "5.10",
        "title": "Sunflower Chlorogenic Acid Off-Note",
        "decision": "defer",
        "reason": "Important for sunflower-specific formulation work, but it does not unlock the current PPI/SPI benchmark or extrusion validation path.",
        "revisit_after": "once the core SPI or PPI extrusion benchmark is closed",
    },
    {
        "item_id": "5.11",
        "title": "Transport / Diffusion Model for Volatile Release",
        "decision": "defer",
        "reason": "A diffusion model is a second-order release refinement; the repo still needs direct benchmark release data before replacing scalar retention corrections.",
        "revisit_after": "after die-exit or cooling release measurements exist",
    },
]


def _priority_score(priority: str) -> float:
    return float(_PRIORITY_WEIGHTS.get(str(priority).strip().lower(), 1.0))


def _under_root(root: Path, path: Path) -> Path:
    """``path`` (a ``data_paths`` constant) re-rooted under a scratch ``root``."""
    if Path(root).resolve() == data_paths.REPO_ROOT:
        return path
    return root / data_paths.rel(path)


def _build_gap_maps(root: Path) -> tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    intake = data_access.load_json(_under_root(root, data_paths.BENCHMARK_INTAKE_REGISTRY))
    process = data_access.load_json(_under_root(root, data_paths.PROCESS_GAP_REGISTRY))
    intake_gaps = {
        str(entry.get("id", "unknown")): dict(entry)
        for entry in intake.get("structural_gaps", [])
        if isinstance(entry, Mapping)
    }
    process_gaps = {
        str(entry.get("gap_id", "unknown")): dict(entry)
        for entry in process.get("entries", [])
        if isinstance(entry, Mapping)
    }
    return intake_gaps, process_gaps


def _rank_workstreams(
    intake_gaps: Mapping[str, Mapping[str, Any]],
    process_gaps: Mapping[str, Mapping[str, Any]],
) -> List[Dict[str, Any]]:
    ranked: List[Dict[str, Any]] = []
    for definition in _WORKSTREAMS:
        intake_rows = [dict(intake_gaps[gap_id]) for gap_id in definition["intake_gap_ids"] if gap_id in intake_gaps]
        process_rows = [dict(process_gaps[gap_id]) for gap_id in definition["process_gap_ids"] if gap_id in process_gaps]
        intake_score = sum(_priority_score(str(row.get("priority", "low"))) for row in intake_rows)
        process_score = sum(1.0 + 0.25 * len(row.get("blocks_modules", []) or []) for row in process_rows)
        ready_input_bonus = 0.75 * len(definition.get("repo_ready_inputs", []))
        phase_bonus = max(0.0, 4.0 - float(definition.get("phase_order", 4)))
        unlock_score = round(intake_score + process_score + ready_input_bonus + phase_bonus, 2)
        near_miss_entries = []
        for row in intake_rows:
            for item in row.get("near_miss_candidates", []) or []:
                if not isinstance(item, Mapping):
                    continue
                entry_id = str(item.get("entry_id", "unknown"))
                if entry_id not in near_miss_entries:
                    near_miss_entries.append(entry_id)
        ranked.append(
            {
                "workstream_id": definition["id"],
                "title": definition["title"],
                "phase_order": int(definition.get("phase_order", 99)),
                "unlock_score": unlock_score,
                "intake_gap_ids": [row.get("id", row.get("gap_id", "unknown")) for row in intake_rows],
                "process_gap_ids": [row.get("gap_id", "unknown") for row in process_rows],
                "intake_gap_count": len(intake_rows),
                "process_gap_count": len(process_rows),
                "repo_ready_inputs": list(definition.get("repo_ready_inputs", [])),
                "deliverables": list(definition.get("deliverables", [])),
                "why_now": str(definition.get("why_now", "")),
                "near_miss_entries": near_miss_entries,
            }
        )

    ranked.sort(key=lambda row: (-float(row.get("unlock_score", 0.0)), int(row.get("phase_order", 99)), row.get("workstream_id", "unknown")))
    for index, row in enumerate(ranked):
        row["recommendation"] = "run_now" if index == 0 else "next" if index == 1 else "later"
    return ranked


def build_structural_unlock_triage(root: Path = data_paths.REPO_ROOT) -> Dict[str, Any]:
    backlog = build_literature_backlog_artifact(root)
    intake_gaps, process_gaps = _build_gap_maps(root)
    ranked = _rank_workstreams(intake_gaps, process_gaps)
    top = ranked[0] if ranked else None
    nxt = ranked[1] if len(ranked) > 1 else None
    return {
        "schema_version": "1.0",
        "source_paths": [
            _under_root(root, data_paths.BENCHMARK_INTAKE_REGISTRY).relative_to(root).as_posix(),
            _under_root(root, data_paths.PROCESS_GAP_REGISTRY).relative_to(root).as_posix(),
        ],
        "summary": {
            "ready_runtime_backlog_count": int(backlog.get("summary", {}).get("ready_runtime_count", 0)),
            "ready_benchmark_backlog_count": int(backlog.get("summary", {}).get("ready_benchmark_count", 0)),
            "citation_backlog_exhausted": (
                int(backlog.get("summary", {}).get("ready_runtime_count", 0)) == 0
                and int(backlog.get("summary", {}).get("ready_benchmark_count", 0)) == 0
            ),
            "recommended_now": None if top is None else top.get("workstream_id"),
            "recommended_next": None if nxt is None else nxt.get("workstream_id"),
        },
        "minimum_primary_experiment": dict(backlog.get("minimum_primary_experiment", {})),
        "ranked_workstreams": ranked,
        "deferred_modeling_review": list(_DEFERRED_MODELING_REVIEW),
    }


def render_structural_unlock_triage_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    minimum_primary_experiment = dict(payload.get("minimum_primary_experiment", {}))
    lines = [
        "# Structural Unlock Triage",
        "",
        f"Citation backlog exhausted: {summary.get('citation_backlog_exhausted', False)}",
        f"Ready runtime backlog: {summary.get('ready_runtime_backlog_count', 0)}",
        f"Ready benchmark backlog: {summary.get('ready_benchmark_backlog_count', 0)}",
        f"Recommended now: {summary.get('recommended_now', 'none')}",
        f"Recommended next: {summary.get('recommended_next', 'none')}",
        "",
        "## Ranked Workstreams",
        "",
        "| Workstream | Score | Recommendation | Intake Gaps | Process Gaps | Near Miss Citations |",
        "| --- | ---: | --- | ---: | ---: | --- |",
    ]
    for row in payload.get("ranked_workstreams", []):
        near_misses = ", ".join(row.get("near_miss_entries", [])) or "none"
        lines.append(
            f"| {row.get('workstream_id', 'unknown')} | {float(row.get('unlock_score', 0.0)):.2f} | {row.get('recommendation', 'later')} | {int(row.get('intake_gap_count', 0))} | {int(row.get('process_gap_count', 0))} | {near_misses} |"
        )

    if minimum_primary_experiment:
        lines.extend([
            "",
            "## Existing Wet-Lab Contract",
            "",
            f"Matrices: {', '.join(minimum_primary_experiment.get('matrices', [])) or 'none'}",
            f"Precursors: {', '.join(f'{key}={value}' for key, value in minimum_primary_experiment.get('exogenous_precursors', {}).items()) or 'none'}",
            f"Analytical panel: {', '.join(minimum_primary_experiment.get('analytical_panel', [])) or 'none'}",
            f"Instrumentation: {minimum_primary_experiment.get('instrumentation', 'unknown')}",
        ])

    for row in payload.get("ranked_workstreams", []):
        lines.extend([
            "",
            f"## {row.get('title', row.get('workstream_id', 'unknown'))}",
            "",
            f"Why now: {row.get('why_now', 'unknown')}",
            f"Repo-ready inputs: {', '.join(row.get('repo_ready_inputs', [])) or 'none'}",
            f"Deliverables: {', '.join(row.get('deliverables', [])) or 'none'}",
            f"Intake gaps: {', '.join(row.get('intake_gap_ids', [])) or 'none'}",
            f"Process gaps: {', '.join(row.get('process_gap_ids', [])) or 'none'}",
        ])

    lines.extend([
        "",
        "## Deferred Modeling Review",
        "",
        "| Item | Decision | Revisit After | Reason |",
        "| --- | --- | --- | --- |",
    ])
    for row in payload.get("deferred_modeling_review", []):
        lines.append(
            f"| {row.get('item_id', 'unknown')} {row.get('title', '')} | {row.get('decision', 'defer')} | {row.get('revisit_after', 'unknown')} | {row.get('reason', '')} |"
        )

    return "\n".join(lines) + "\n"