#!/usr/bin/env python
"""
FREEZE THE FIT / HOLD-OUT WAVE GENERATORS (2026-09-03).

``scripts/generators/generate_kinetic_core_b*_{fit,holdout,reports,scorers}.py`` and
``probe_amine_fate_b2_4.py`` are the record of how every frozen parameter of the kinetic
core was produced. They are not re-run; no test executes them; a change to any of them is
a NEW WAVE, not an edit. This script writes their SHA-256 into
``results/validation/wave_generators_manifest.json`` and
``tests/scientific/test_wave_generators_frozen.py`` fails when a file no longer matches --
so a wave generator cannot be edited without regenerating the manifest in the same commit
and saying why in ``scripts/generators/WAVES.md``.

Usage:
    python scripts/generators/build_wave_manifest.py            # rewrite the manifest
    python scripts/generators/build_wave_manifest.py --check    # exit 1 if anything moved
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

MANIFEST = data_paths.VALIDATION_DIR / "wave_generators_manifest.json"
GENERATORS = ROOT / "scripts" / "generators"
PATTERNS = ("generate_kinetic_core_b*.py", "probe_amine_fate_b2_4.py")


def wave_files():
    seen = set()
    for pattern in PATTERNS:
        for path in sorted(GENERATORS.glob(pattern)):
            if path.name.endswith(("_laplace.py", "_profile.py")):
                continue  # derived-artifact generators, not waves
            if path not in seen:
                seen.add(path)
                yield path


def current() -> dict:
    return {
        data_paths.rel(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in wave_files()
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    now = current()
    if args.check:
        if not MANIFEST.exists():
            print(f"{data_paths.rel(MANIFEST)} is missing"); return 1
        recorded = json.loads(MANIFEST.read_text(encoding="utf-8"))["files"]
        moved = sorted(set(now) ^ set(recorded)) + sorted(k for k in now if k in recorded and now[k] != recorded[k])
        if moved:
            print("wave generators changed without a new manifest:\n  " + "\n  ".join(moved)); return 1
        print(f"wave manifest: {len(now)} generators frozen, none changed"); return 0
    MANIFEST.write_text(
        json.dumps(
            {
                "artifact": "wave_generators_manifest",
                "rule": (
                    "these files are the record of how the kinetic core's frozen parameters were "
                    "produced; they are not re-run and a change to any of them is a new wave, "
                    "recorded in scripts/generators/WAVES.md with the manifest rebuilt"
                ),
                "files": now,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {data_paths.rel(MANIFEST)}: {len(now)} wave generators frozen")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
