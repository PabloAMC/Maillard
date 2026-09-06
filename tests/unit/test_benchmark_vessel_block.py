"""
Programme step R1 (2026-09-06): every panel bundle carries a DECIDED vessel block, the
arithmetic behind the scorecard's oxygen line is what the docstring says, and the three
systems the review compared come out where the review put them.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import vessel
from src.kinetic_core.panel import load_bundle


def _panel_paths():
    yield from sorted(data_paths.BENCHMARKS_DIR.glob("*.json"))
    yield from sorted(data_paths.EXTERNAL_VALIDATION_DIR.glob("*.json"))
    yield from sorted(data_paths.MAILLARD_PATH_HOLDOUT_DIR.glob("*.json"))


PANEL = list(_panel_paths())


def test_the_panel_is_the_size_the_completion_script_declares():
    assert len(PANEL) == 39


@pytest.mark.parametrize("path", PANEL, ids=lambda p: p.stem)
def test_every_panel_bundle_has_a_decided_vessel_block(path: Path):
    block = (load_bundle(path).get("conditions") or {}).get("vessel")
    assert isinstance(block, dict), f"{path.stem}: no vessel block"
    for field in vessel.REQUIRED_FIELDS:
        assert field in block, f"{path.stem}: vessel block lacks {field}"
    assert block["atmosphere"] in vessel.ATMOSPHERES, path.stem
    assert block["water_source"] in vessel.WATER_SOURCES, path.stem
    assert block["provenance_class"] in vessel.PROVENANCE_CLASSES, path.stem
    assert len(block["provenance_note"]) >= 120, path.stem
    if block["provenance_class"] == "primary_source_pdf":
        assert "data/articles/" in block["provenance_note"], path.stem
    if block["provenance_class"] == "repo_verbatim_methods_quote":
        note = block["provenance_note"].upper()
        assert "NOT ON DISK" in note and "SECOND-HAND" in note, path.stem
    if block["atmosphere"] == "not_applicable":
        assert block["fill_mL"] is None and block["vessel_mL"] is None, path.stem
    if block["fill_mL"] is not None and block["vessel_mL"] is not None:
        assert block["vessel_mL"] > block["fill_mL"], path.stem


def test_the_completion_script_is_current():
    import subprocess
    import sys
    done = subprocess.run(
        [sys.executable, "scripts/generators/complete_benchmark_vessel_fields.py", "--check"],
        cwd=data_paths.REPO_ROOT, capture_output=True, text=True,
    )
    assert done.returncode == 0, done.stderr + done.stdout


def test_o2_per_ml_of_headspace_is_the_ideal_gas_number():
    # 0.20946 * 101325 Pa * 1e-6 m3 / (8.314 * 293.15) = 8.706e-6 mol
    assert vessel.O2_MMOL_PER_ML_HEADSPACE == pytest.approx(8.706e-3, rel=1e-3)


def _record(stem: str):
    for path in PANEL:
        if path.stem == stem:
            return vessel.oxygen_record(load_bundle(path))
    raise AssertionError(stem)


def test_the_three_systems_the_review_compared():
    """Yiltirak ~2.0, Bolton ~2.0, Hofmann Table 1 ~0.26 mol O2 per mol thiol."""
    y = _record("mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026")
    assert y["status"] == "computed"
    assert y["headspace_mL"] == 17.0
    assert y["o2_mmol"] == pytest.approx(0.149, abs=0.002)
    assert y["thiol_mmol"] == pytest.approx(0.075)
    assert y["o2_to_thiol"] == pytest.approx(2.0, abs=0.05)
    h = _record("hofmann1998_ribose_cysteine_145C_20min_pH5")
    assert h["headspace_mL"] == 100.0
    assert h["thiol_mmol"] == pytest.approx(3.3)
    assert h["o2_to_thiol"] == pytest.approx(0.27, abs=0.02)
    b = _record("thiamine_cys_glucose_120C_Bolton1994")
    assert b["status"] == "computed"
    assert 1.5 < b["o2_to_thiol"] < 2.6


def test_fed_hydrogen_sulfide_counts_as_a_thiol_and_thiamine_does_not():
    r = _record("hofmann1998_norfuraneol_h2s_145C_20min_pH5")
    assert r["thiol_mmol"] == pytest.approx(1.0, rel=0.05)  # 1 mmol H2S in 50 mL
    b = load_bundle(next(p for p in PANEL if p.stem == "thiamine_cys_glucose_120C_Bolton1994"))
    only_cys = sum(
        v["concentration_mM"] for k, v in b["precursors"].items() if "cystein" in k.lower()
    ) * 33.3 / 1000.0
    assert vessel.thiol_mmol(b, 33.3) == pytest.approx(only_cys)


def test_statuses_cover_every_non_computed_case():
    seen = {vessel.oxygen_record(load_bundle(p))["status"] for p in PANEL}
    assert {"computed", "ambiguous", "open", "continuous", "not_applicable"} <= seen
    assert "missing" not in seen


def test_format_cell():
    assert vessel.format_o2_to_thiol({"status": "open"}) == "open"
    assert vessel.format_o2_to_thiol({"status": "computed", "o2_to_thiol": 1.984, "o2_mmol": 0.15}) == "1.98"
    assert vessel.format_o2_to_thiol({"status": "computed", "o2_to_thiol": None, "o2_mmol": 0.12}) == "0.120 mmol O2, no thiol"
