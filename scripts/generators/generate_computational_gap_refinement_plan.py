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
    build_computational_gap_dft_job_manifest,
    build_computational_gap_refinement_plan_artifact,
    build_computational_gap_xtb_job_manifest,
    render_computational_gap_refinement_plan_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--manifest-dir", default="results/computational_gap_refinement")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    manifest_dir = Path(args.manifest_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    plan_payload = build_computational_gap_refinement_plan_artifact()
    xtb_manifest = build_computational_gap_xtb_job_manifest(plan_payload)
    dft_manifest = build_computational_gap_dft_job_manifest(plan_payload)

    plan_md = output_dir / "computational_gap_refinement_plan.md"
    plan_json = output_dir / "computational_gap_refinement_plan.json"
    xtb_json = manifest_dir / "computational_gap_xtb_job_manifest.json"
    dft_json = manifest_dir / "computational_gap_dft_job_manifest.json"

    plan_md.write_text(render_computational_gap_refinement_plan_markdown(plan_payload), encoding="utf-8")
    plan_json.write_text(json.dumps(plan_payload, indent=2), encoding="utf-8")
    xtb_json.write_text(json.dumps(xtb_manifest, indent=2), encoding="utf-8")
    dft_json.write_text(json.dumps(dft_manifest, indent=2), encoding="utf-8")

    print(plan_md.read_text(encoding="utf-8"))
    print(f"Wrote {plan_md}")
    print(f"Wrote {plan_json}")
    print(f"Wrote {xtb_json}")
    print(f"Wrote {dft_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())