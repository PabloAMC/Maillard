from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chemistry_benchmark_validator import (
    build_adoption_decisions_from_assessment,
    build_mlp_assessment_artifact,
    render_mlp_assessment_markdown,
)
from src.mlp_adoption_contract import build_adoption_note_payload, render_adoption_note_markdown
from src.mlp_external_benchmarks import build_external_mlp_landscape_payload, render_external_mlp_landscape_markdown
from src.reaction_benchmark import build_reaction_benchmark_artifact, render_reaction_benchmark_markdown
from src.results_db import ResultsDB


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build the MLP reaction benchmark, assessment, external-landscape "
            "and adoption notes into results/validation/, and record the "
            "adoption decisions in the results DB."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "results" / "validation"),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    benchmark_payload = build_reaction_benchmark_artifact()
    assessment_payload = build_mlp_assessment_artifact()
    external_payload = build_external_mlp_landscape_payload()
    decisions = build_adoption_decisions_from_assessment(assessment_payload)
    adoption_payload = build_adoption_note_payload(
        decisions,
        benchmark_set_id=str(assessment_payload["summary"]["benchmark_set_id"]),
    )

    files = {
        output_dir / "mlp_reaction_benchmark.json": benchmark_payload,
        output_dir / "mlp_assessment.json": assessment_payload,
        output_dir / "mlp_external_mlp_landscape.json": external_payload,
        output_dir / "mlp_adoption_notes.json": adoption_payload,
    }
    for path, payload in files.items():
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    (output_dir / "mlp_reaction_benchmark.md").write_text(
        render_reaction_benchmark_markdown(benchmark_payload),
        encoding="utf-8",
    )
    (output_dir / "mlp_assessment.md").write_text(
        render_mlp_assessment_markdown(assessment_payload),
        encoding="utf-8",
    )
    (output_dir / "mlp_external_mlp_landscape.md").write_text(
        render_external_mlp_landscape_markdown(external_payload),
        encoding="utf-8",
    )
    (output_dir / "mlp_adoption_notes.md").write_text(
        render_adoption_note_markdown(adoption_payload),
        encoding="utf-8",
    )

    db = ResultsDB()
    for decision in decisions:
        db.add_ml_adoption_decision({
            "candidate_id": decision.candidate_id,
            "model_family": decision.model_family,
            "model_name": decision.model_name,
            "proposed_role": decision.proposed_role,
            "decision": decision.decision,
            "benchmark_set_id": decision.benchmark_set_id,
            "coverage_ratio": decision.coverage_ratio,
            "rank_correlation": decision.rank_correlation,
            "mean_abs_error_kcal": decision.mean_abs_error_kcal,
            "max_abs_error_kcal": decision.max_abs_error_kcal,
            "stop_reasons": decision.stop_reasons,
            "rationale": decision.rationale,
            "fallback_comparator": decision.fallback_comparator,
            "benchmark_visible_gap": decision.benchmark_visible_gap,
            "approved_for_default": decision.approved_for_default,
        })


if __name__ == "__main__":
    main()