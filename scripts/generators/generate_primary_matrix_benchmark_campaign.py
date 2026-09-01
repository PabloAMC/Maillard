#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.primary_benchmark_campaign import (  # noqa: E402
    export_matrix_primary_benchmark_campaign,
    export_primary_matrix_external_package,
    export_primary_matrix_external_package_intake_template,
    render_matrix_primary_benchmark_campaign_markdown,
    render_primary_matrix_external_package_markdown,
    render_primary_matrix_external_package_intake_template_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Export the primary matrix benchmark campaign, its external package "
            "and the package intake template into results/validation/."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(data_paths.VALIDATION_DIR),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
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