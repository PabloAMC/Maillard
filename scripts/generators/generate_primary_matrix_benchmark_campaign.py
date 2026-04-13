#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.primary_benchmark_campaign import (  # noqa: E402
    export_matrix_primary_benchmark_campaign,
    export_primary_matrix_external_package,
    export_primary_matrix_external_package_intake_template,
    render_matrix_primary_benchmark_campaign_markdown,
    render_primary_matrix_external_package_markdown,
    render_primary_matrix_external_package_intake_template_markdown,
)


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    campaign = export_matrix_primary_benchmark_campaign(str(output_dir), root=ROOT)
    package = export_primary_matrix_external_package(str(output_dir), root=ROOT)
    intake_template = export_primary_matrix_external_package_intake_template(str(output_dir), root=ROOT)
    print(render_matrix_primary_benchmark_campaign_markdown(campaign))
    print(render_primary_matrix_external_package_markdown(package))
    print(render_primary_matrix_external_package_intake_template_markdown(intake_template))
    print(f"Wrote {output_dir / 'matrix_primary_benchmark_campaign.md'}")
    print(f"Wrote {output_dir / 'matrix_primary_benchmark_campaign.json'}")
    print(f"Wrote {output_dir / 'primary_matrix_external_package.md'}")
    print(f"Wrote {output_dir / 'primary_matrix_external_package.json'}")
    print(f"Wrote {output_dir / 'primary_matrix_external_package_intake_template.md'}")
    print(f"Wrote {output_dir / 'primary_matrix_external_package_intake_template.yaml'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())