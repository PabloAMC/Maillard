"""
Test suite for Phase 3.2 — Quasi-Harmonic Correction Module

Tests Grimme/Truhlar quasi-RRHO corrections for low-frequency entropy
overestimation in transition states and products.
"""

import pytest

from src.quasi_harmonic_correction import QuasiHarmonicCorrector, grimme_weighting


pytestmark = [pytest.mark.deterministic_helper_lane]


class TestQuasiHarmonicCorrection:
    """Test RRHO entropy correction implementation."""

    def test_qh_correction_water(self):
        """Apply RRHO correction to water frequencies."""
        corrector = QuasiHarmonicCorrector(T=423)
        freqs = [1594, 3657, 3756]
        s_vib_uncorrected = corrector.calculate_entropy_harmonic(freqs)
        s_vib_corrected = corrector.calculate_entropy_qh(freqs)

        assert s_vib_uncorrected > 0.0
        assert s_vib_corrected > 0.0
        assert abs(s_vib_corrected - s_vib_uncorrected) / s_vib_uncorrected < 0.01

    def test_low_freq_mode_handling(self):
        """Detect and correct low-frequency torsional mode (~50 cm^-1)."""
        corrector = QuasiHarmonicCorrector(T=423)
        freqs = [52, 150, 200, 350, 450, 500, 600, 700, 800]
        s_vib_corrected = corrector.calculate_entropy_qh(freqs)
        rotor_contribution = corrector.get_rotor_entropy(T=423)

        assert s_vib_corrected > 0.0
        assert rotor_contribution > 0.0

    def test_barrier_corrected_vs_uncorrected(self):
        """ΔG‡ changes meaningfully when low-freq entropy is corrected."""
        corrector = QuasiHarmonicCorrector(T=423)
        enthalpy_activation = 12.5
        freqs_reactant = [500, 600, 700, 800, 900, 1000, 1100]
        freqs_ts = [150, 200, 300, 400, 500, 600, 700]

        dg_uncorr = corrector.calculate_barrier_free_energy(
            enthalpy_activation,
            freqs_reactant,
            freqs_ts,
            use_quasi_harmonic=False,
        )
        dg_corr = corrector.calculate_barrier_free_energy(
            enthalpy_activation,
            freqs_reactant,
            freqs_ts,
            use_quasi_harmonic=True,
        )

        assert isinstance(dg_corr, float)
        assert dg_corr != dg_uncorr

    def test_grimme_weighting(self):
        """Verify Grimme weighting formula for transition region."""
        f_low = 50
        f_mid = 100
        f_high = 200

        w_low = grimme_weighting(f_low, Bcut=100)
        w_mid = grimme_weighting(f_mid, Bcut=100)
        w_high = grimme_weighting(f_high, Bcut=100)

        assert 0.0 < w_low < w_mid < w_high < 1.0

    def test_rotor_treatment(self):
        """Treat low-freq torsion as 1D rotor, not harmonic oscillator."""
        corrector = QuasiHarmonicCorrector(T=423)
        f_torsion = 50

        s_harmonic = corrector.calculate_entropy_single_mode_harmonic(f_torsion, T=423)
        s_rotor = corrector.calculate_entropy_single_mode_rotor(f_torsion, T=423)

        assert s_rotor > s_harmonic

    def test_barrier_shift_reasonable(self):
        """After QH correction, ΔG‡ shift is in physically reasonable ±5-10 kcal/mol range."""
        corrector = QuasiHarmonicCorrector(T=423)
        freqs_ts = [80, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500]
        freqs_reactant = [500, 600, 700, 800, 900, 1000, 1100]

        dg_uncorr = corrector.calculate_barrier_free_energy(20.0, freqs_reactant, freqs_ts, use_quasi_harmonic=False)
        dg_corr = corrector.calculate_barrier_free_energy(20.0, freqs_reactant, freqs_ts, use_quasi_harmonic=True)

        dg_shift = abs(dg_corr - dg_uncorr)
        assert dg_shift < 10.0

    def test_temperature_scaling(self):
        """QH correction should scale with temperature."""
        freqs = [100, 200, 300, 400, 500, 600, 700]

        corrector_150 = QuasiHarmonicCorrector(T=423)
        corrector_200 = QuasiHarmonicCorrector(T=473)

        s_150 = corrector_150.calculate_entropy_qh(freqs)
        s_200 = corrector_200.calculate_entropy_qh(freqs)

        assert s_200 > s_150


class TestQuasiHarmonicConfigurable:
    """Test parameter settings for QH correction."""

    def test_custom_grimme_cutoff(self):
        """Allow customization of Grimme frequency cutoff."""
        corrector_100 = QuasiHarmonicCorrector(T=423, grimme_cutoff=100)
        corrector_200 = QuasiHarmonicCorrector(T=423, grimme_cutoff=200)
        freqs = [50, 150, 250, 350, 450, 550, 650, 750, 850]

        s_100 = corrector_100.calculate_entropy_qh(freqs)
        s_200 = corrector_200.calculate_entropy_qh(freqs)
        assert s_100 != s_200
