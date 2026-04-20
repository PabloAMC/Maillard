import json
import sys
from pathlib import Path

from scripts import run_computational_gap_xtb as runner


def test_runner_reuses_cached_outputs_in_execute_mode(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "ROOT", tmp_path)

    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    runner_dir = tmp_path / "xtb_case"
    runner_dir.mkdir()

    reactant_path = runner_dir / "reactant.xyz"
    product_path = runner_dir / "product.xyz"
    run_script_path = runner_dir / "run_xtb.sh"
    xtbpath_path = runner_dir / "xtbpath.xyz"
    ts_path = runner_dir / "xtbpath_ts.xyz"

    reactant_path.write_text("2\nreactant\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n", encoding="utf-8")
    product_path.write_text("2\nproduct\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n", encoding="utf-8")
    run_script_path.write_text("#!/bin/bash\nexit 1\n", encoding="utf-8")
    xtbpath_path.write_text("2\npath\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n", encoding="utf-8")
    ts_path.write_text("2\nts\nH 0.0 0.0 0.0\nH 0.0 0.0 0.75\n", encoding="utf-8")

    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "lysinoalanine_crosslink",
                        "reaction_key": "lysinoalanine_crosslink",
                        "status": "ready",
                        "runner_script": str(run_script_path.relative_to(runner.ROOT)),
                        "required_inputs": [
                            str(reactant_path.relative_to(runner.ROOT)),
                            str(product_path.relative_to(runner.ROOT)),
                        ],
                        "expected_outputs": [
                            str(xtbpath_path.relative_to(runner.ROOT)),
                            str(ts_path.relative_to(runner.ROOT)),
                        ],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_xtb.py",
            "--manifest",
            str(manifest_path),
            "--output",
            str(output_path),
            "--execute",
        ],
    )

    assert runner.main() == 0

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["summary"]["completed_count"] == 1
    assert payload["summary"]["completed_cached_count"] == 1
    assert payload["summary"]["failed_count"] == 0
    assert payload["jobs"][0]["status"] == "completed_cached"
    assert payload["jobs"][0]["cached_outputs_used"] is True
    assert payload["jobs"][0]["quality_gate_passed"] is True


def test_runner_supports_explicit_proxy_target_and_target_specific_output(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.chdir(tmp_path)

    proxy_dir = tmp_path / "data/geometries/xtb_inputs/pe_amadori"
    proxy_dir.mkdir(parents=True)
    (proxy_dir / "reactant.xyz").write_text("2\nreactant\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n", encoding="utf-8")
    (proxy_dir / "product.xyz").write_text("2\nproduct\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n", encoding="utf-8")
    (proxy_dir / "run_xtb.sh").write_text("#!/bin/bash\nexit 1\n", encoding="utf-8")
    (proxy_dir / "xtbpath.xyz").write_text("2\npath\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n", encoding="utf-8")
    (proxy_dir / "xtbpath_ts.xyz").write_text("2\nts\nH 0.0 0.0 0.0\nH 0.0 0.0 0.75\n", encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_xtb.py",
            "--manifest",
            str(tmp_path / "missing_manifest.json"),
            "--target",
            "pe_amadori",
            "--execute",
        ],
    )

    assert runner.main() == 0

    output_path = tmp_path / "results/computational_gap_refinement/pe_amadori_xtb_execution.json"
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["summary"]["job_count"] == 1
    assert payload["jobs"][0]["target_id"] == "pe_amadori"
    assert payload["jobs"][0]["status"] == "completed_cached"