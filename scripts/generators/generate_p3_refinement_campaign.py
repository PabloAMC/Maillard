#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.refinement_campaign import (
    build_cheap_screening_artifact,
    build_global_sensitivity_artifact,
    build_refinement_impact_artifact,
    build_selective_dft_plan,
    render_cheap_screening_markdown,
    render_global_sensitivity_markdown,
    render_refinement_impact_markdown,
    render_selective_dft_plan_markdown,
)
from src.p3_refinement_governance import (
    build_p3_refinement_governance_artifact,
    render_p3_refinement_governance_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    parser.add_argument("--patch-path", default="data/lit/refinement_surrogate_patches.json")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    patch_path = Path(args.patch_path)
    patch_path.parent.mkdir(parents=True, exist_ok=True)

    global_payload = build_global_sensitivity_artifact(target_tag=args.target_tag)
    cheap_payload = build_cheap_screening_artifact(target_tag=args.target_tag)
    dft_payload = build_selective_dft_plan(target_tag=args.target_tag)
    impact_payload = build_refinement_impact_artifact(target_tag=args.target_tag)
    governance_payload = build_p3_refinement_governance_artifact(target_tag=args.target_tag)

    files = {
        output_dir / "p3_global_sensitivity.md": render_global_sensitivity_markdown(global_payload),
        output_dir / "cheap_refinement_screening.md": render_cheap_screening_markdown(cheap_payload),
        output_dir / "selective_dft_plan.md": render_selective_dft_plan_markdown(dft_payload),
        output_dir / "refinement_impact.md": render_refinement_impact_markdown(impact_payload),
        output_dir / "p3_refinement_governance.md": render_p3_refinement_governance_markdown(governance_payload),
    }
    json_payloads = {
        output_dir / "p3_global_sensitivity.json": global_payload,
        output_dir / "cheap_refinement_screening.json": cheap_payload,
        output_dir / "selective_dft_plan.json": dft_payload,
        output_dir / "p3_offline_dft_jobs.json": dft_payload.get("offline_jobs", []),
        output_dir / "refinement_impact.json": impact_payload,
        output_dir / "p3_refinement_governance.json": governance_payload,
    }

    for path, content in files.items():
        path.write_text(content, encoding="utf-8")
    for path, payload in json_payloads.items():
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    patch_path.write_text(json.dumps(impact_payload.get("patch", {}), indent=2), encoding="utf-8")

    print(files[output_dir / "p3_global_sensitivity.md"])
    print(files[output_dir / "cheap_refinement_screening.md"])
    print(files[output_dir / "selective_dft_plan.md"])
    print(files[output_dir / "refinement_impact.md"])
    print(files[output_dir / "p3_refinement_governance.md"])
    print(f"Wrote {output_dir / 'p3_global_sensitivity.md'}")
    print(f"Wrote {output_dir / 'p3_global_sensitivity.json'}")
    print(f"Wrote {output_dir / 'cheap_refinement_screening.md'}")
    print(f"Wrote {output_dir / 'cheap_refinement_screening.json'}")
    print(f"Wrote {output_dir / 'selective_dft_plan.md'}")
    print(f"Wrote {output_dir / 'selective_dft_plan.json'}")
    print(f"Wrote {output_dir / 'p3_offline_dft_jobs.json'}")
    print(f"Wrote {output_dir / 'refinement_impact.md'}")
    print(f"Wrote {output_dir / 'refinement_impact.json'}")
    print(f"Wrote {output_dir / 'p3_refinement_governance.md'}")
    print(f"Wrote {output_dir / 'p3_refinement_governance.json'}")
    print(f"Wrote {patch_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())