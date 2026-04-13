#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extrusion_benchmark_execution import (  # noqa: E402
    export_extrusion_disulfide_follow_on_execution,
    render_extrusion_disulfide_follow_on_execution_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True)
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    artifact = export_extrusion_disulfide_follow_on_execution(args.workbook, args.output_dir, root=ROOT)
    print(render_extrusion_disulfide_follow_on_execution_markdown(artifact))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())