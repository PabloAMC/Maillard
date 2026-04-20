from pathlib import Path

import src.computational_gap_proxy_readiness as proxy


def test_proxy_readiness_artifact_enforces_family_15_pair_gate(monkeypatch, tmp_path: Path):
    root = tmp_path
    pe_schiff = root / "data/geometries/xtb_inputs/pe_schiff_base"
    pe_amadori = root / "data/geometries/xtb_inputs/pe_amadori"
    execution_dir = root / "results/computational_gap_refinement"
    pe_schiff.mkdir(parents=True)
    pe_amadori.mkdir(parents=True)
    execution_dir.mkdir(parents=True)

    frame = "3\nframe\nH 0 0 0\nH 0 0 1\nH 0 1 0\n"
    for directory in [pe_schiff, pe_amadori]:
        (directory / "xtbpath_1.xyz").write_text(frame, encoding="utf-8")
        (directory / "xtbpath.xyz").write_text(frame, encoding="utf-8")
        (directory / "xtbpath_ts.xyz").write_text(frame, encoding="utf-8")

    execution_payload = {
        "jobs": [
            {
                "status": "completed",
                "quality_gate_passed": True,
            }
        ]
    }
    (execution_dir / "pe_schiff_base_xtb_execution.json").write_text(__import__("json").dumps(execution_payload), encoding="utf-8")
    (execution_dir / "pe_amadori_xtb_execution.json").write_text(__import__("json").dumps(execution_payload), encoding="utf-8")

    (pe_schiff / "xtbopt.log").write_text("SCF not converged, aborting\n", encoding="utf-8")
    (pe_amadori / "xtbopt.log").write_text("normal termination of xtb\n", encoding="utf-8")

    monkeypatch.setattr(proxy, "ROOT", root)
    payload = proxy.build_computational_gap_proxy_readiness_artifact()
    by_id = {row["target_id"]: row for row in payload["lanes"]}

    assert payload["summary"]["family_15_pair_ready"] is False
    assert by_id["pe_schiff_base"]["promotion_status"] == "proxy_only_quality_warning"
    assert by_id["pe_amadori"]["promotion_status"] == "proxy_only_pair_gate"


def test_proxy_readiness_accepts_completed_cached_as_managed_execution(monkeypatch, tmp_path: Path):
    root = tmp_path
    pe_schiff = root / "data/geometries/xtb_inputs/pe_schiff_base"
    pe_amadori = root / "data/geometries/xtb_inputs/pe_amadori"
    execution_dir = root / "results/computational_gap_refinement"
    pe_schiff.mkdir(parents=True)
    pe_amadori.mkdir(parents=True)
    execution_dir.mkdir(parents=True)

    frame = "3\nframe\nH 0 0 0\nH 0 0 1\nH 0 1 0\n"
    for directory in [pe_schiff, pe_amadori]:
        (directory / "xtbpath_1.xyz").write_text(frame, encoding="utf-8")
        (directory / "xtbpath.xyz").write_text(frame, encoding="utf-8")
        (directory / "xtbpath_ts.xyz").write_text(frame, encoding="utf-8")
        (directory / "xtbopt.log").write_text("normal termination of xtb\n", encoding="utf-8")

    execution_payload = {
        "jobs": [
            {
                "status": "completed_cached",
                "quality_gate_passed": True,
            }
        ]
    }
    (execution_dir / "pe_schiff_base_xtb_execution.json").write_text(__import__("json").dumps(execution_payload), encoding="utf-8")
    (execution_dir / "pe_amadori_xtb_execution.json").write_text(__import__("json").dumps(execution_payload), encoding="utf-8")

    monkeypatch.setattr(proxy, "ROOT", root)
    payload = proxy.build_computational_gap_proxy_readiness_artifact()

    assert payload["summary"]["family_15_pair_ready"] is True
    assert payload["summary"]["formal_candidate_count"] == 2