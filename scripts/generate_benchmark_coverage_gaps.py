#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.generate_benchmark_coverage_gaps import _build_rows, _render_markdown, main


if __name__ == "__main__":
    raise SystemExit(main())