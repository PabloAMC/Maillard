"""Smoke tests for the command-line entry points that survive retirement step B5.

2026-09-03: `scripts/run_pipeline.py`, `run_campaign.py`, `run_cantera_kinetics.py`,
`optimize_formulation.py`, `explain_formulation.py` and `ingest_results.py` were deleted with
the legacy lane. The front door is `scripts/maillard.py`; the two core artifacts have their
own generators. Each must at least parse `--help` and exit 0 from a fresh interpreter.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

ENTRY_POINTS = (
    "scripts/maillard.py",
    "scripts/generators/generate_core_panel_scores.py",
    "scripts/generators/generate_core_prediction_uncertainty.py",
    "scripts/generators/generate_model_card.py",
    "scripts/generators/build_data_readme.py",
)


@pytest.mark.parametrize("script", ENTRY_POINTS)
def test_entry_point_help_exits_zero(script):
    proc = subprocess.run(
        [sys.executable, script, "--help"], cwd=str(ROOT), capture_output=True, text=True, timeout=300
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.strip()


def test_maillard_verbs_are_the_five_core_verbs():
    proc = subprocess.run(
        [sys.executable, "scripts/maillard.py", "--help"], cwd=str(ROOT), capture_output=True, text=True, timeout=300
    )
    for verb in ("compare", "predict", "explain", "rank", "score"):
        assert verb in proc.stdout
    assert "--lane" not in proc.stdout
