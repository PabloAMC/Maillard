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
from src.experiment_value import ExperimentCandidate

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
