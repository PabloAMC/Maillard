#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_refinement import (  # noqa: E402
    build_computational_gap_dft_ingestion_artifact,
    build_computational_gap_dft_job_manifest,
    render_computational_gap_dft_ingestion_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="results/computational_gap_refinement/computational_gap_dft_job_manifest.json")
    parser.add_argument("--execution", default="results/computational_gap_refinement/computational_gap_dft_execution.json")
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)
    if manifest_path.exists():
        manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        manifest_payload = build_computational_gap_dft_job_manifest()

    payload = build_computational_gap_dft_ingestion_artifact(
        manifest_payload=manifest_payload,
        execution_path=Path(args.execution),
    )

    markdown_path = output_dir / "computational_gap_dft_ingestion_report.md"
    json_path = output_dir / "computational_gap_dft_ingestion_report.json"
    markdown_path.write_text(render_computational_gap_dft_ingestion_markdown(payload), encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown_path.read_text(encoding="utf-8"))
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())