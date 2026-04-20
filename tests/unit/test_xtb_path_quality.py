from pathlib import Path

from src.xtb_path_quality import assess_xtb_path_quality


def test_assess_xtb_path_quality_materializes_outputs_and_detects_scf_warnings(tmp_path: Path):
    runner_dir = tmp_path / "pe_schiff_base"
    runner_dir.mkdir()
    frame = "3\nframe\nH 0 0 0\nH 0 0 1\nH 0 1 0\n"
    (runner_dir / "xtbpath_1.xyz").write_text(frame, encoding="utf-8")
    (runner_dir / "xtbpath_2.xyz").write_text(frame, encoding="utf-8")
    (runner_dir / "xtbopt.log").write_text(
        "SCF not converged, aborting\nSetup of Coulomb evaluator failed\n",
        encoding="utf-8",
    )

    payload = assess_xtb_path_quality(runner_dir, materialize_missing=True)

    assert payload["frame_count"] == 2
    assert payload["has_path_bundle"] is True
    assert payload["has_ts_guess"] is True
    assert payload["quality_gate_passed"] is False
    assert "SCF not converged, aborting" in payload["failure_markers"]
    assert len(payload["synthesized_outputs"]) == 2