"""S22.2 — Convert top-N VoI candidates into concrete lab requests.

For each requested candidate we emit:

* `data/protocols/requested_<request_id>.yaml` — pre-filled intake YAML in
  the same schema used by `data/protocols/example_matrix_experiment_intake.yaml`,
  but with a `status: pending_lab` field and empty `measured_volatiles` ready
  for lab fill-in.
* `results/validation/experiment_requests/<request_id>.md` — human protocol
  combining the VoI rationale with the matching DOE_TEMPLATES recipe from
  `src.doe_generator`.

The function is intentionally small. It reuses, but does not mutate, the
existing benchmark/registry artefacts, and never auto-promotes anything —
the request is *requested*, not *fulfilled*.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

import yaml

from src.benchmark_validation import (
    get_benchmark_files,
    load_benchmark,
)
from src.doe_generator import DOE_TEMPLATES
from src.experiment_value import ExperimentCandidate, infer_matrix_family

ROOT = Path(__file__).resolve().parents[1]
PROTOCOLS_DIR = ROOT / "data" / "protocols"
REQUESTS_DIR = ROOT / "results" / "validation" / "experiment_requests"


def _display_path(path: Path) -> str:
    """Repo-relative path when possible, otherwise the path itself."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


@dataclass(frozen=True)
class ExperimentRequest:
    request_id: str
    candidate_rank: int
    benchmark_id: str
    compound: str
    voi_score: float
    suggested_doe_template: str
    intake_yaml_path: Path
    protocol_md_path: Path
    budget_label: Optional[str] = None
    goal: Optional[str] = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _slug(text: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return s or "request"


def _benchmark_index() -> Dict[str, Dict[str, Any]]:
    """Map benchmark_id -> raw benchmark dict for fast lookup."""
    index: Dict[str, Dict[str, Any]] = {}
    for path in get_benchmark_files():
        try:
            bench = load_benchmark(path)
        except Exception:
            continue
        bid = str(bench.get("benchmark_id") or path.stem)
        index[bid] = bench
    return index


def _matches_filter(value: Optional[str], needle: Optional[str]) -> bool:
    if needle is None:
        return True
    if value is None:
        return False
    return needle.lower() in str(value).lower()


# ---------------------------------------------------------------------------
# Analytical-context defaults (mirror intake schema's `analytical_context`)
# ---------------------------------------------------------------------------
#
# These keep the CRO loop tight: the protocol Markdown asks for exactly the
# fields that `data/protocols/matrix_experiment_intake_schema.json` already
# accepts, and the pre-filled intake YAML carries the same defaults so a
# returned payload lands cleanly via `scripts/ingest_results.py` without
# extra translation. We never claim a measurement here; status stays
# `pending_lab` until the lab fills `measured_volatiles`.

# Map DoE template -> default analytical_context block. Each entry is a
# best-default starting point; the lab is expected to confirm or override
# in the returned YAML, not us.
_TEMPLATE_ANALYTICAL_CONTEXT: Dict[str, Dict[str, Any]] = {
    "blocking_benchmark_gap": {
        "headspace_method": "HS-SPME-GC-MS",
        "quantification_mode": "internal_standard_calibrated",
        "replicates": 3,
        "non_detect_policy": "report_lod_and_do_not_backfill",
    },
    "missing_absolute_anchor": {
        "headspace_method": "HS-SPME-GC-MS",
        "quantification_mode": "internal_standard_calibrated",
        "replicates": 3,
        "non_detect_policy": "report_lod_and_do_not_backfill",
    },
    "missing_positive_flavor_anchor": {
        "headspace_method": "SPME-GC-MS_polar_fiber",
        "quantification_mode": "internal_standard_calibrated",
        "replicates": 3,
        "non_detect_policy": "report_lod_and_do_not_backfill",
    },
    "missing_kinetic_dataset": {
        "headspace_method": "LC-MS_MS_timecourse_quench",
        "quantification_mode": "internal_standard_calibrated",
        "replicates": 3,
        "non_detect_policy": "report_lod_and_do_not_backfill",
    },
    "missing_process_state_bundle": {
        "headspace_method": "simultaneous_DSC_OPA_Ellman",
        "quantification_mode": "absolute_calibrated_reagent_assay",
        "replicates": 3,
        "non_detect_policy": "report_lod_and_do_not_backfill",
    },
}

_DEFAULT_ANALYTICAL_CONTEXT: Dict[str, Any] = {
    "headspace_method": "HS-SPME-GC-MS",
    "quantification_mode": "internal_standard_calibrated",
    "replicates": 3,
    "non_detect_policy": "report_lod_and_do_not_backfill",
}

# Suggest matched isotopically-labeled internal standards based on compound
# name substrings. Conservative on purpose — when nothing matches we ask the
# CRO to nominate one rather than guessing wrong.
_INTERNAL_STANDARD_HINTS: Sequence[tuple[str, str]] = (
    ("furanthiol", "13C-2-methyl-3-furanthiol"),
    ("furfurylthiol", "13C-2-furfurylthiol"),
    ("methional", "d3-methional"),
    ("hexanal", "hexanal-d12"),
    ("nonanal", "nonanal-d17"),
    ("pyrazine", "13C-2,5-dimethylpyrazine"),
    ("hemf", "13C-HEMF"),
    ("dmhf", "13C-DMHF"),
    ("furosine", "13C6-furosine"),
    ("acrylamide", "13C3-acrylamide"),
)


def _suggest_internal_standards(compound: str, template_key: str) -> List[str]:
    """Return suggested isotopically-labeled internal standards for `compound`.

    Falls back to a generic placeholder so the CRO row is never silently
    empty — silence here would let a returned payload skip the calibration
    contract entirely.
    """
    name = (compound or "").lower()
    matched: List[str] = []
    for needle, label in _INTERNAL_STANDARD_HINTS:
        if needle in name and label not in matched:
            matched.append(label)
    if not matched:
        matched = [f"compound_specific_internal_standard_for_{_slug(compound)}"]
    if template_key == "blocking_benchmark_gap" and "hexanal-d12" not in matched:
        # Multi-factorial template covers both meaty + lipid bands; ensure
        # the lipid anchor is present.
        matched.append("hexanal-d12")
    return matched


def _default_analytical_context(
    *, template_key: str, compound: str
) -> Dict[str, Any]:
    base = dict(_TEMPLATE_ANALYTICAL_CONTEXT.get(template_key, _DEFAULT_ANALYTICAL_CONTEXT))
    base["internal_standards"] = _suggest_internal_standards(compound, template_key)
    base["notes"] = (
        "Defaults emitted by next-experiment; lab is expected to confirm or "
        "override before returning the intake YAML."
    )
    return base


def _build_intake_payload(
    *,
    request_id: str,
    candidate: ExperimentCandidate,
    bench: Optional[Mapping[str, Any]],
    goal: Optional[str],
    budget_label: Optional[str],
) -> Dict[str, Any]:
    conditions = dict((bench or {}).get("conditions", {}) or {})
    formulation = dict((bench or {}).get("precursors", {}) or {})
    payload: Dict[str, Any] = {
        "experiment_id": request_id,
        "status": "pending_lab",
        "source_kind": "model_requested_experiment",
        "protein_type": (bench or {}).get("protein_type", "free"),
        "process_state": (bench or {}).get("process_state")
        or (bench or {}).get("metadata", {}).get("process_state")
        or "aqueous_pre_extrusion_model",
        "matrix_format": (bench or {}).get("metadata", {}).get("matrix_format"),
        "conditions": conditions,
        "formulation": {"precursors": formulation},
        "measured_volatiles": {
            candidate.compound: {
                "conc_ppb": None,
                "uncertainty_pct": None,
                "_target_role": "primary_request_target",
            }
        },
        "request_metadata": {
            "originating_voi_rank": candidate.rank,
            "voi_score": candidate.voi_score,
            "voi_rationale": candidate.rationale,
            "envelope_miss_log10": candidate.envelope_miss_log10,
            "ci_width_log10": candidate.ci_width_log10,
            "decision_relevance": candidate.decision_relevance,
            "p5_predicted_ppb": candidate.predicted_p5,
            "p50_predicted_ppb": candidate.predicted_p50,
            "p95_predicted_ppb": candidate.predicted_p95,
            "previous_measured_ppb": candidate.measured_ppb,
            "odour_threshold_ug_per_kg": candidate.odour_threshold_ug_per_kg,
            "suggested_doe_template": candidate.suggested_doe_template,
            "goal": goal,
            "budget_label": budget_label,
        },
        "benchmark_alignment": {
            "target_benchmark_id": candidate.benchmark_id,
            "notes": "Re-measure to constrain the model envelope at this benchmark.",
        },
        "analytical_context": _default_analytical_context(
            template_key=candidate.suggested_doe_template,
            compound=candidate.compound,
        ),
    }
    return payload


def _render_protocol_markdown(
    *,
    request_id: str,
    candidate: ExperimentCandidate,
    template_key: str,
    intake_path: Path,
    bench: Optional[Mapping[str, Any]],
    goal: Optional[str],
    budget_label: Optional[str],
) -> str:
    template = DOE_TEMPLATES.get(template_key, {})
    factors = template.get("factors", []) or []
    instructions = template.get("instructions", "")
    method = template.get("method", "")
    bench_ctx_lines: List[str] = []
    if bench is not None:
        cond = bench.get("conditions", {}) or {}
        bench_ctx_lines = [
            f"- Benchmark target: `{candidate.benchmark_id}`",
            f"- Protein type: `{bench.get('protein_type', 'free')}`",
            f"- Conditions: {json.dumps(cond)}",
        ]
    bench_ctx = "\n".join(bench_ctx_lines) if bench_ctx_lines else "- (No upstream benchmark context found.)"

    factor_lines = "\n".join(f"  - {f}" for f in factors) if factors else "  - (use template defaults)"
    width = (
        f"{candidate.ci_width_log10:.2f} dex"
        if candidate.ci_width_log10 is not None
        else "n/a"
    )
    miss = (
        f"{candidate.envelope_miss_log10:.2f} dex"
        if not candidate.inside_ci
        else "inside 90% CI"
    )
    odt = (
        f"{candidate.odour_threshold_ug_per_kg:g} µg/kg"
        if candidate.odour_threshold_ug_per_kg is not None
        else "unknown"
    )
    goal_line = f"- Goal: {goal}" if goal else ""
    budget_line = f"- Budget label: {budget_label}" if budget_label else ""

    analytical = _default_analytical_context(
        template_key=template_key, compound=candidate.compound
    )
    is_lines = "\n".join(f"  - `{label}`" for label in analytical["internal_standards"])
    analytical_block = (
        f"- `headspace_method`: `{analytical['headspace_method']}`\n"
        f"- `quantification_mode`: `{analytical['quantification_mode']}`\n"
        f"- `replicates`: `{analytical['replicates']}` (minimum)\n"
        f"- `non_detect_policy`: `{analytical['non_detect_policy']}`\n"
        f"- `internal_standards`:\n{is_lines}"
    )

    # 2026-08-27 (Wave I). The dynamic-range line used to quote ONLY the model's own
    # predicted midpoint. Where the benchmark already carries a published measurement,
    # that is the number a lab must be able to see: on the rank-1 card the model midpoint
    # was 0.017 ppb against a published 13 ppb -- a ~760x window -- and a CRO calibrating
    # to the model's figure could have set a range that cannot detect the real value.
    # Lead with the measurement, and state the model's figure as the disagreement it is.
    _measured = float(getattr(candidate, "measured_ppb", 0.0) or 0.0)
    _predicted = float(candidate.predicted_p50 or 0.0)
    if _measured > 0.0 and _predicted > 0.0:
        _fold = max(_measured, _predicted) / max(min(_measured, _predicted), 1e-12)
        dynamic_range_line = (
            "- [ ] Confirm target compound identity and calibrate over a range spanning "
            f"BOTH the published value (≈ {_measured:g} ppb) and the model's 90% CI "
            f"midpoint (≈ {_predicted:g} ppb). They differ by {_fold:.3g}× — that "
            "disagreement is *why* this experiment is ranked, so a calibration range "
            "covering only one of them cannot resolve it.\n"
        )
    elif _measured > 0.0:
        dynamic_range_line = (
            "- [ ] Confirm target compound identity and calibrate around the published "
            f"value (≈ {_measured:g} ppb); the model predicts effectively zero here.\n"
        )
    else:
        dynamic_range_line = (
            "- [ ] Confirm target compound identity and expected dynamic range. No "
            "published value exists for this row, so the model's 90% CI midpoint "
            f"(≈ {_predicted:g} ppb) is the only guide — a starting point, not a target.\n"
        )

    cro_checklist = (
        dynamic_range_line
        + "- [ ] Procure or confirm availability of the suggested isotopically-labeled "
        "internal standards above; substitute and **note the swap** in `analytical_context.notes`.\n"
        f"- [ ] Run ≥ {analytical['replicates']} biological replicates plus a process blank and a "
        "matrix blank.\n"
        "- [ ] Measure and report LoD and LoQ; do **not** backfill non-detects "
        "(`non_detect_policy: report_lod_and_do_not_backfill`).\n"
        "- [ ] Quantify against the internal standards (NOT semi-quant external "
        "calibration); set `quantification_mode: internal_standard_calibrated`.\n"
        "- [ ] Record the measured `conc_ppb` and `uncertainty_pct` (1σ relative) "
        f"under `measured_volatiles.{candidate.compound}` in the pre-filled intake YAML.\n"
        "- [ ] Fill `provenance.source_doi` (or internal lab batch ID), "
        "`provenance.measurement_date`, and `provenance.notes` (instrument, method file).\n"
        "- [ ] Set `status: measured` and return the YAML; do not edit `request_metadata` "
        "(it carries the upstream VoI rationale)."
    )

    return f"""# Experiment Request `{request_id}`

VoI rank **#{candidate.rank}**, score **{candidate.voi_score:.2f}**.
Generated from the current Monte-Carlo benchmark envelope.

## Why this experiment

- Compound: **{candidate.compound}**
- Envelope miss: {miss}
- 90% CI width: {width}
- ODT: {odt}; decision relevance {candidate.decision_relevance:.2f}
- Rationale: {candidate.rationale}
{goal_line}
{budget_line}

## Benchmark context

{bench_ctx}

## Suggested protocol (`{template_key}`)

- Method: **{method}**
- Factors:
{factor_lines}

### Instructions

{instructions}

## Analytical context (mirror in intake YAML `analytical_context`)

{analytical_block}

## CRO send-to-lab checklist

{cro_checklist}

## Data return

- Pre-filled intake YAML: `{_display_path(intake_path)}`
- Fill in `measured_volatiles.{candidate.compound}.conc_ppb` (and any
  recommended_same_run_extensions) with the lab readout, set `status:
  measured`, and re-run the propagation to update the envelope.
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_requests(
    candidates: Sequence[ExperimentCandidate],
    *,
    top_n: int = 5,
    protein_type: Optional[str] = None,
    matrix_filter: Optional[Sequence[str]] = None,
    goal: Optional[str] = None,
    budget_label: Optional[str] = None,
    protocols_dir: Path = PROTOCOLS_DIR,
    requests_dir: Path = REQUESTS_DIR,
    benchmark_index: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> List[ExperimentRequest]:
    """Materialise the top-N candidates into intake YAML + protocol markdown."""
    if benchmark_index is None:
        benchmark_index = _benchmark_index()
    protocols_dir.mkdir(parents=True, exist_ok=True)
    requests_dir.mkdir(parents=True, exist_ok=True)

    wanted_matrices = (
        {str(m).strip().lower() for m in matrix_filter if str(m).strip()}
        if matrix_filter
        else set()
    )

    requests: List[ExperimentRequest] = []
    for candidate in candidates:
        if len(requests) >= top_n:
            break
        bench = benchmark_index.get(candidate.benchmark_id)
        if not _matches_filter(
            (bench or {}).get("protein_type") if bench else None,
            protein_type,
        ):
            continue
        if wanted_matrices:
            inferred = infer_matrix_family(candidate.benchmark_id, bench)
            if inferred.lower() not in wanted_matrices:
                continue
        request_id = _slug(
            f"requested_{candidate.benchmark_id}_{candidate.compound}_rank{candidate.rank}"
        )
        intake_path = protocols_dir / f"{request_id}.yaml"
        protocol_path = requests_dir / f"{request_id}.md"

        intake_payload = _build_intake_payload(
            request_id=request_id,
            candidate=candidate,
            bench=bench,
            goal=goal,
            budget_label=budget_label,
        )
        intake_path.write_text(
            yaml.safe_dump(intake_payload, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        protocol_md = _render_protocol_markdown(
            request_id=request_id,
            candidate=candidate,
            template_key=candidate.suggested_doe_template,
            intake_path=intake_path,
            bench=bench,
            goal=goal,
            budget_label=budget_label,
        )
        protocol_path.write_text(protocol_md, encoding="utf-8")

        requests.append(
            ExperimentRequest(
                request_id=request_id,
                candidate_rank=candidate.rank,
                benchmark_id=candidate.benchmark_id,
                compound=candidate.compound,
                voi_score=candidate.voi_score,
                suggested_doe_template=candidate.suggested_doe_template,
                intake_yaml_path=intake_path,
                protocol_md_path=protocol_path,
                budget_label=budget_label,
                goal=goal,
            )
        )
    return requests


def render_index_markdown(requests: Sequence[ExperimentRequest]) -> str:
    if not requests:
        return "# Experiment Requests\n\n_No requests generated for the current filter._\n"
    lines = [
        "# Experiment Requests",
        "",
        f"Generated {len(requests)} request(s).",
        "",
        "| Rank | Request ID | Benchmark | Compound | VoI | DoE template | Intake YAML | Protocol |",
        "| ---: | --- | --- | --- | ---: | --- | --- | --- |",
    ]
    for r in requests:
        lines.append(
            f"| {r.candidate_rank} | `{r.request_id}` | `{r.benchmark_id}` | "
            f"{r.compound} | {r.voi_score:.2f} | `{r.suggested_doe_template}` | "
            f"`{_display_path(r.intake_yaml_path)}` | "
            f"`{_display_path(r.protocol_md_path)}` |"
        )
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_index(
    requests: Sequence[ExperimentRequest],
    *,
    requests_dir: Path = REQUESTS_DIR,
) -> Path:
    requests_dir.mkdir(parents=True, exist_ok=True)
    path = requests_dir / "index.md"
    path.write_text(render_index_markdown(requests), encoding="utf-8")
    return path
