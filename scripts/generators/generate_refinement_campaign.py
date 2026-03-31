#!/usr/bin/env python3
"""Legacy wrapper for selective refinement campaign generation."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators.generate_selective_refinement_campaign import main


if __name__ == "__main__":
    raise SystemExit(main())