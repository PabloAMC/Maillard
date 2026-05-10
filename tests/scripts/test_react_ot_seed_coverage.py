import json
from pathlib import Path

import pytest

pytest.importorskip("ase")

from scripts import report_react_ot_seed_coverage as coverage


REACTANT = """3
reactant
O 0.0 0.0 0.0
H 0.0 0.0 0.96
H 0.93 0.0 -0.24
"""

# Same atom ordering, displaced enough that the Kabsch-aligned RMSD
# vs the reactant clearly exceeds the 0.30 Å gate.
SEED_OK = """3
ts
O 0.0 0.0 0.0
H 0.0 0.0 2.40
H 2.20 0.0 -0.24
"""

# Same atom ordering, near-identical to reactant -> warn (rmsd < 0.3)
SEED_TOO_CLOSE = """3
ts
O 0.0 0.0 0.0
H 0.0 0.0 0.97
H 0.94 0.0 -0.24
"""

XTB_TS = """3
xtb_ts
O 0.0 0.0 0.0
H 0.0 0.0 1.35
H 1.25 0.0 -0.24
"""


def _stage_target(
    target: str,
    seed_dir: Path,
    xtb_inputs: Path,
    *,
    seed_xyz: str | None,
    reactant_xyz: str | None = REACTANT,
    xtb_ts_xyz: str | None = XTB_TS,
) -> None:
    target_dir = xtb_inputs / target
    target_dir.mkdir(parents=True, exist_ok=True)
    if reactant_xyz is not None:
        (target_dir / "reactant.xyz").write_text(reactant_xyz)
    if xtb_ts_xyz is not None:
        (target_dir / "xtbpath_ts.xyz").write_text(xtb_ts_xyz)
    if seed_xyz is not None:
        seed_dir.mkdir(parents=True, exist_ok=True)
        (seed_dir / f"{target}_react_ot_seed.xyz").write_text(seed_xyz)


def test_evaluate_target_ok(tmp_path):
    seed_dir = tmp_path / "results"
    xtb_inputs = tmp_path / "xtb_inputs"
    target = "lysinoalanine_crosslink"
    _stage_target(target, seed_dir, xtb_inputs, seed_xyz=SEED_OK)

    # Bypass REPO_ROOT-based path stripping in the evaluator:
    # _evaluate_target stores paths relative to REPO_ROOT for display.
    coverage.REPO_ROOT = tmp_path  # type: ignore[attr-defined]

    record = coverage._evaluate_target(target, seed_dir, xtb_inputs)
    assert record["status"] == "ok"
    assert record["atom_count_match"]
    assert record["rmsd_vs_reactant_passes_min_gate"]
    assert record["rmsd_vs_reactant_angstrom"] >= coverage.MIN_REACTANT_TS_RMSD_ANGSTROM
    assert isinstance(record["rmsd_vs_xtb_ts_angstrom"], float)


def test_evaluate_target_warn_when_too_close_to_reactant(tmp_path):
    seed_dir = tmp_path / "results"
    xtb_inputs = tmp_path / "xtb_inputs"
    target = "pe_schiff_base"
    _stage_target(target, seed_dir, xtb_inputs, seed_xyz=SEED_TOO_CLOSE)
    coverage.REPO_ROOT = tmp_path  # type: ignore[attr-defined]

    record = coverage._evaluate_target(target, seed_dir, xtb_inputs)
    assert record["status"] == "warn"
    assert not record["rmsd_vs_reactant_passes_min_gate"]


def test_evaluate_target_missing_seed(tmp_path):
    seed_dir = tmp_path / "results"
    xtb_inputs = tmp_path / "xtb_inputs"
    target = "aa_ring_open_dicarbonyl"
    _stage_target(target, seed_dir, xtb_inputs, seed_xyz=None)
    coverage.REPO_ROOT = tmp_path  # type: ignore[attr-defined]

    record = coverage._evaluate_target(target, seed_dir, xtb_inputs)
    assert record["status"] == "missing_seed"
    assert record["seed_present"] is False


def test_build_report_summary_and_markdown(tmp_path):
    seed_dir = tmp_path / "results"
    xtb_inputs = tmp_path / "xtb_inputs"
    _stage_target("lysinoalanine_crosslink", seed_dir, xtb_inputs, seed_xyz=SEED_OK)
    _stage_target("pe_schiff_base", seed_dir, xtb_inputs, seed_xyz=SEED_TOO_CLOSE)
    _stage_target("aa_ring_open_dicarbonyl", seed_dir, xtb_inputs, seed_xyz=None)
    _stage_target("asparagine_sugar_explicit_water_cluster", seed_dir, xtb_inputs, seed_xyz=SEED_OK)
    coverage.REPO_ROOT = tmp_path  # type: ignore[attr-defined]

    report = coverage.build_report(seed_dir, xtb_inputs, list(coverage.ELIGIBLE_TARGETS))
    summary = report["summary"]
    assert summary["n_targets"] == 4
    assert summary["n_ok"] == 2
    assert summary["n_warn"] == 1
    assert summary["n_missing_seed"] == 1
    assert report["trust_posture"]["is_runtime_authority"] is False

    md = coverage.render_markdown(report)
    assert "React-OT seed coverage report" in md
    assert "lysinoalanine_crosslink" in md
    assert "missing seeds" in md


def test_orchestrator_help_runs():
    """Smoke: importing the orchestrator and parsing --help must succeed."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "scripts/orchestrate_react_ot_colab.py", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "React-OT Colab loop" in result.stdout or "orchestrator" in result.stdout
