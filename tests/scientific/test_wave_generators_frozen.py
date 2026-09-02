"""The fit / hold-out wave generators are frozen (2026-09-03, owner decision, step 3).

They are the record of how every frozen kinetic-core parameter was produced; nothing
re-runs them and no test executes them. A change to any of them is a NEW WAVE: rebuild
the manifest (scripts/generators/build_wave_manifest.py) and record the wave in
scripts/generators/WAVES.md in the same commit. This test turns a silent edit into a
red build.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "results" / "validation" / "wave_generators_manifest.json"


def test_wave_generators_match_the_manifest():
    proc = subprocess.run(
        [sys.executable, "scripts/generators/build_wave_manifest.py", "--check"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_the_manifest_covers_every_wave_generator():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = set(manifest["files"])
    on_disk = {
        str(p.relative_to(ROOT))
        for p in (ROOT / "scripts" / "generators").glob("generate_kinetic_core_b*.py")
        if not p.name.endswith(("_laplace.py", "_profile.py"))
    } | {"scripts/generators/probe_amine_fate_b2_4.py"}
    assert files == on_disk, sorted(files ^ on_disk)
    assert "new wave" in manifest["rule"]
