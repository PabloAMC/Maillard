#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extrusion_benchmark_execution import export_extrusion_diagnostic_examples_bundle  # noqa: E402


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    payload = export_extrusion_diagnostic_examples_bundle(output_dir, root=ROOT)
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())