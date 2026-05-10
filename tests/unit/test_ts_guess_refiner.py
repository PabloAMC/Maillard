"""Unit tests for P0 — TS-guess quality gate & refinement orchestrator.

All xTB calls are mocked; these tests validate the decision logic:
  - Happy path: good guess passes through unchanged
  - Gate triggers when fmax > threshold
  - Candidate rejected when n_imag != 1
  - Candidate rejected when energy is absurd
  - Cascade: relaxed_scan tried before cineb
  - Original accepted when bad-but-not-invalid (n_imag=1, high fmax)
  - Original rejected when clearly a minimum (n_imag=0)
  - Policy "off" skips everything
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from src.ts_guess_refiner import RefinementOutcome, refine_if_needed

_DUMMY_XYZ = "3\ntest\nH 0 0 0\nH 0 0 1\nH 0 1 0\n"
_DUMMY_REACTANT = "3\nreactant\nH 0 0 0\nH 0 0 2\nH 0 2 0\n"
_DUMMY_PRODUCT = "3\nproduct\nH 0 0 0\nH 0 0 0.9\nH 0 0.9 0\n"


# Reactant energy reference: -10.5 Eh.  TS probes should be within
# [-5, +200] kcal/mol of the reactant to pass the energy sanity check.
# 0.05 Eh ≈ 31 kcal/mol (normal barrier).
_REACTANT_EH = -10.5
_TS_EH = -10.45  # ~31 kcal/mol above reactant


def _make_probe(fmax=0.3, energy_eh=_TS_EH, n_imag=1, lowest_freq_cm=-300.0, highest_imag_cm=-300.0):
    return {
        "fmax_ev_ang": fmax,
        "energy_eh": energy_eh,
        "n_imag": n_imag,
        "lowest_freq_cm": lowest_freq_cm,
        "highest_imag_cm": highest_imag_cm,
    }


def _make_scan_result(fmax=0.2, energy_ev=-270.0):
    return {
        "xyz": _DUMMY_XYZ,
        "energy_ev": energy_ev,
        "fmax_ev_ang": fmax,
        "scan_distance": 1.5,
        "scan_points": [],
        "forming_bond": [0, 1],
    }


class TestHappyPath:
    """Good TS guess should pass through unchanged."""

    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_good_guess_passes(self, mock_probe):
        mock_probe.side_effect = [
            _make_probe(fmax=0.3, n_imag=1, energy_eh=_TS_EH),  # TS probe
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant probe
        ]
        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT,
            threshold_fmax=1.0,
        )
        assert result.accepted is True
        assert result.final_method == "original"
        assert result.xyz == _DUMMY_XYZ

    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_policy_off_skips(self, mock_probe):
        mock_probe.side_effect = [
            _make_probe(fmax=0.3, n_imag=1, energy_eh=_TS_EH),
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),
        ]
        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT,
            policy="off",
        )
        assert result.accepted is True
        assert result.final_method == "original"
        assert "off" in result.reason


class TestGateTriggering:
    """Gate triggers when fmax > threshold or n_imag != 1."""

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_high_fmax_triggers_refinement(self, mock_probe, mock_scan, mock_neb):
        # TS probe: bad fmax; relaxed scan produces good candidate
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # initial TS
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            _make_probe(fmax=0.5, n_imag=1, energy_eh=_TS_EH),  # scan candidate
        ]
        mock_scan.return_value = _make_scan_result(fmax=0.5)

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.final_method == "relaxed_scan"
        assert result.accepted is True
        mock_scan.assert_called_once()


class TestCandidateRejection:
    """Candidates with n_imag != 1 or absurd energy should be rejected."""

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_scan_candidate_with_zero_imag_accepted_if_fmax_improves(self, mock_probe, mock_scan, mock_neb):
        """Scan max with n_imag=0 accepted when fmax improves (scan max = TS by construction)."""
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=0, energy_eh=_TS_EH),  # initial TS — minimum
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            # Scan probe: n_imag=0, energy may mismatch reactant, but fmax improved
            _make_probe(fmax=0.2, n_imag=0, energy_eh=_REACTANT_EH - 3.0),
        ]
        mock_scan.return_value = _make_scan_result()

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        # Scan accepted: fmax improved 2.5→0.2 (scan max = TS region)
        assert result.accepted is True
        assert result.final_method == "relaxed_scan"

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_scan_candidate_with_zero_imag_and_no_fmax_improvement_rejected(self, mock_probe, mock_scan, mock_neb):
        """Scan with n_imag=0 and no fmax improvement stays on original."""
        mock_probe.side_effect = [
            _make_probe(fmax=0.15, n_imag=1, energy_eh=_TS_EH),  # initial TS — already low fmax
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            _make_probe(fmax=0.3, n_imag=0, energy_eh=_TS_EH),  # scan — n_imag=0 AND worse fmax
        ]
        mock_scan.return_value = _make_scan_result(fmax=0.3)

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
            policy="always",
        )
        # Scan not better → original kept (n_imag=1, sane)
        assert result.accepted is True
        assert result.final_method == "original"

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_absurd_energy_rejected(self, mock_probe, mock_scan, mock_neb):
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # initial TS (sane)
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            # Scan candidate: energy way too high (>200 kcal over reactant)
            _make_probe(fmax=0.2, n_imag=1, energy_eh=_REACTANT_EH + 1.0),  # +627 kcal/mol!
        ]
        mock_scan.return_value = _make_scan_result()

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.final_method == "original"  # scan candidate rejected

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan", return_value=None)
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_original_minimum_rejected_when_no_alternatives(self, mock_probe, mock_scan, mock_neb):
        """If original is a minimum (0 imag) and no alternatives work → reject."""
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=0, energy_eh=_TS_EH),  # initial TS is a minimum
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
        ]

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.accepted is False
        assert "invalid" in result.reason.lower() or "rejected" in result.reason.lower() or "failed" in result.reason.lower()


class TestCascadeOrder:
    """Relaxed scan tried before CI-NEB."""

    @patch("src.ts_guess_refiner.run_xtb_cineb")
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_scan_wins_when_better_than_neb(self, mock_probe, mock_scan, mock_neb):
        # Scan succeeds with better fmax → scan wins even though NEB also tried
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # initial TS
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            _make_probe(fmax=0.3, n_imag=1, energy_eh=_TS_EH),  # scan candidate
            _make_probe(fmax=0.5, n_imag=1, energy_eh=_TS_EH),  # NEB candidate (worse)
        ]
        mock_scan.return_value = _make_scan_result(fmax=0.3)
        mock_neb.return_value = {
            "xyz": _DUMMY_XYZ,
            "energy_ev": -270.0,
            "fmax_ev_ang": 0.5,
        }

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.final_method == "relaxed_scan"
        # NEB is always tried now
        mock_neb.assert_called_once()

    @patch("src.ts_guess_refiner.run_xtb_cineb")
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan", return_value=None)
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_neb_wins_when_scan_fails(self, mock_probe, mock_scan, mock_neb):
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # initial TS
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            _make_probe(fmax=0.4, n_imag=1, energy_eh=_TS_EH),  # NEB candidate
        ]
        mock_neb.return_value = {
            "xyz": _DUMMY_XYZ,
            "energy_ev": -270.0,
            "fmax_ev_ang": 0.4,
        }

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.final_method == "cineb"
        mock_neb.assert_called_once()


class TestBadGuessPolicy:
    """'Better bad guess than no guess' — unless clearly invalid."""

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan", return_value=None)
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_bad_but_valid_guess_accepted(self, mock_probe, mock_scan, mock_neb):
        """High fmax but n_imag=1 and sane energy → accept original."""
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # bad fmax but has 1 imag + sane energy
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
        ]

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
        )
        assert result.accepted is True
        assert result.final_method == "original"


class TestAlwaysPolicy:
    """policy='always' should refine even good guesses."""

    @patch("src.ts_guess_refiner.run_xtb_cineb", return_value=None)
    @patch("src.ts_guess_refiner.run_xtb_relaxed_scan")
    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_always_refines_good_guess(self, mock_probe, mock_scan, mock_neb):
        mock_probe.side_effect = [
            _make_probe(fmax=0.3, n_imag=1, energy_eh=_TS_EH),  # initial: already good
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
            _make_probe(fmax=0.15, n_imag=1, energy_eh=_TS_EH),  # scan candidate (even better)
        ]
        mock_scan.return_value = _make_scan_result(fmax=0.15)

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT, _DUMMY_PRODUCT,
            threshold_fmax=1.0,
            policy="always",
        )
        assert result.final_method == "relaxed_scan"


class TestNoProductXYZ:
    """Without product_xyz, scan and NEB are skipped."""

    @patch("src.ts_guess_refiner.probe_ts_guess_xtb")
    def test_no_product_no_refinement(self, mock_probe):
        mock_probe.side_effect = [
            _make_probe(fmax=2.5, n_imag=1, energy_eh=_TS_EH),  # initial TS: bad fmax, sane energy
            _make_probe(fmax=0.1, n_imag=0, energy_eh=_REACTANT_EH),  # reactant
        ]

        result = refine_if_needed(
            _DUMMY_XYZ, _DUMMY_REACTANT,
            product_xyz=None,
            threshold_fmax=1.0,
        )
        # No product → no scan/NEB, but n_imag=1 → accept original
        assert result.accepted is True
        assert result.final_method == "original"


class TestDetectFormingBond:
    """detect_forming_bond should prefer heavy-heavy pairs over X-H pairs."""

    def test_prefers_heavy_pair_over_proton_transfer(self):
        from src.xtb_backend import detect_forming_bond

        # Reactant: C at origin, N far away (3.2 Å), H even farther (4.8 Å)
        reactant = (
            "3\nreactant\n"
            "C  0.0  0.0  0.0\n"
            "N  3.2  0.0  0.0\n"
            "H  4.8  0.0  0.0\n"
        )
        # Product: C-N bond formed (1.47 Å), H migrated to C (1.1 Å)
        # H has larger Δd (3.7) vs C-N Δd (1.73), but C-N is heavy-heavy.
        product = (
            "3\nproduct\n"
            "C  0.0  0.0  0.0\n"
            "N  1.47 0.0  0.0\n"
            "H  1.1  0.0  0.0\n"
        )
        pair = detect_forming_bond(reactant, product)
        assert pair == (0, 1), f"Expected heavy pair (0,1)=(C,N), got {pair}"

    def test_falls_back_to_h_pair_when_no_heavy(self):
        from src.xtb_backend import detect_forming_bond

        # Only H-X bonds form (pure proton transfer)
        reactant = "2\nreactant\nO  0.0  0.0  0.0\nH  3.0  0.0  0.0\n"
        product = "2\nproduct\nO  0.0  0.0  0.0\nH  0.96 0.0  0.0\n"
        pair = detect_forming_bond(reactant, product)
        assert pair == (0, 1)
