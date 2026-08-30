from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.geometry_benchmark import build_geometry_benchmark_artifact, render_geometry_benchmark_markdown
from src.geometry_benchmark_validator import build_geometry_assessment_artifact, render_geometry_assessment_markdown


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build the MLP geometry benchmark and its assessment; writes "
            "results/validation/mlp_geometry_benchmark.{json,md} and "
            "mlp_geometry_assessment.{json,md}."
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

    benchmark_payload = build_geometry_benchmark_artifact()
    assessment_payload = build_geometry_assessment_artifact()

    files = {
        output_dir / "mlp_geometry_benchmark.json": benchmark_payload,
        output_dir / "mlp_geometry_assessment.json": assessment_payload,
    }
    for path, payload in files.items():
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    (output_dir / "mlp_geometry_benchmark.md").write_text(
        render_geometry_benchmark_markdown(benchmark_payload),
        encoding="utf-8",
    )
    (output_dir / "mlp_geometry_assessment.md").write_text(
        render_geometry_assessment_markdown(assessment_payload),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()