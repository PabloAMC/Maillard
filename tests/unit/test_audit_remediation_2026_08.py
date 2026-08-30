"""Regression tests for the 2026-08-26 audit remediation (tasks/audit_remediation.md)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_experiment_intake import calibrate_from_intake
from src.thermo import QuasiHarmonicCorrector


def test_calibrate_from_intake_refuses_holdout_evidence():
    payload = {
        "experiment_id": "holdout_leak_attempt",
        "source_kind": "external_literature",
        "evidence_class": "external_validation_only",
        "protein_type": "pea_iso",
        "process_state": "ambient_slurry",
        "conditions": {"temp_C": 40.0, "ph": 6.0, "water_activity": 0.95, "time_min": 10.0},
        "formulation": {"precursors": {"Pea Protein Isolate": {"concentration_mM": 1000.0}}},
        "measured_volatiles": {"Hexanal": {"conc_ppb": 1260.0}},
        "provenance": {
            "origin": "external_literature",
            "source_reference": "hold-out bundle",
            "source_doi": "10.0000/holdout",
            "measurement_date": "2020-01-01",
        },
    }
    with pytest.raises(ValueError, match="external_validation_only"):
        calibrate_from_intake(payload, dry_run=True)


# test_evaluator_reads_water_activity_key was removed 2026-08-26 together with the
# dead-code island (src/benchmark_registry.py, benchmark_evaluator.py,
# benchmark_reporting.py, benchmark_assertions.py, benchmark_markdown.py). The test
# asserted on the source text of benchmark_evaluator.py, which no longer exists; the
# water-activity key contract it guarded lived only in that duplicate lane.


def test_qrrho_rotor_entropy_uses_reduced_moment():
    calc = QuasiHarmonicCorrector(temp_k=298.15)
    # Grimme 2012 reference magnitudes: ~11.9 J/mol/K at 100 cm^-1, ~18.6 at 20 cm^-1.
    assert calc._rotor_entropy(100.0) == pytest.approx(11.9, abs=0.5)
    assert calc._rotor_entropy(20.0) == pytest.approx(18.6, abs=0.7)
    # The broken implementation returned a frequency-independent ~45.9 J/mol/K.
    assert calc._rotor_entropy(100.0) != pytest.approx(calc._rotor_entropy(20.0), abs=1.0)


def test_qrrho_caps_low_mode_entropy_below_harmonic():
    # The entire point of QRRHO: the free-rotor branch must CAP the harmonic
    # divergence for low modes, never exceed it.
    calc = QuasiHarmonicCorrector(temp_k=298.15)
    for freq in (10.0, 20.0, 50.0):
        assert calc._rotor_entropy(freq) < calc._harmonic_entropy(freq)


def test_marker_injection_units_are_mM():
    # The injected lipid markers must land in the same unit basis (mM) as the
    # precursor concentrations they share initial_concentrations with.
    # "benchmark_evaluator" dropped 2026-08-26: the dead-code island it belonged to
    # was deleted, leaving benchmark_validation as the only live copy.
    for module in ("benchmark_validation",):
        source = (ROOT / "src" / f"{module}.py").read_text()
        assert "conc_mM = molar_conc * 1000.0" in source, module
        # 2-Pentylfuran is the C9 pentyl form, not the C8 butyl form.
        assert '"2-Pentylfuran": "CCCCCC1=CC=CO1"' in source, module
        assert '"CCCCC1=CC=CO1"' not in source, module
