from __future__ import annotations

import math
from typing import Iterable, List


PLANCK = 6.62607015e-34
BOLTZMANN = 1.380649e-23
SPEED_OF_LIGHT_CM_S = 29979245800.0
AVOGADRO = 6.02214076e23
JOULE_TO_KCAL = 1.0 / 4184.0


def grimme_weighting(freq_cm1: float, Bcut: float = 100.0) -> float:
    frequency = float(freq_cm1)
    cutoff = float(Bcut)
    if frequency <= 0.0 or cutoff <= 0.0:
        return 0.0
    return 1.0 / (1.0 + (cutoff / frequency) ** 4)


class QuasiHarmonicCorrector:
    def __init__(self, T: float = 298.15, grimme_cutoff: float = 100.0, rotor_moment_of_inertia: float = 1.0e-44):
        self.temperature_k = float(T)
        self.grimme_cutoff = float(grimme_cutoff)
        self.rotor_moment_of_inertia = float(rotor_moment_of_inertia)

    def calculate_entropy_single_mode_harmonic(self, freq_cm1: float, T: float | None = None) -> float:
        temperature = self.temperature_k if T is None else float(T)
        frequency = float(freq_cm1)
        if frequency <= 0.0 or temperature <= 0.0:
            return 0.0
        x = (PLANCK * frequency * SPEED_OF_LIGHT_CM_S) / (BOLTZMANN * temperature)
        if x > 700.0:
            return 0.0
        exp_x = math.exp(x)
        exp_neg_x = math.exp(-x)
        entropy_j_mol_k = AVOGADRO * BOLTZMANN * (x / (exp_x - 1.0) - math.log(1.0 - exp_neg_x))
        return entropy_j_mol_k * JOULE_TO_KCAL

    def calculate_entropy_single_mode_rotor(self, freq_cm1: float, T: float | None = None) -> float:
        temperature = self.temperature_k if T is None else float(T)
        frequency = float(freq_cm1)
        if frequency <= 0.0 or temperature <= 0.0:
            return 0.0
        prefactor = (8.0 * math.pi**3 * self.rotor_moment_of_inertia * BOLTZMANN * temperature) / (PLANCK**2)
        entropy_j_mol_k = AVOGADRO * BOLTZMANN * (0.5 + 0.5 * math.log(prefactor))
        return entropy_j_mol_k * JOULE_TO_KCAL

    def get_rotor_entropy(self, T: float | None = None) -> float:
        return self.calculate_entropy_single_mode_rotor(self.grimme_cutoff / 2.0, T=T)

    def calculate_entropy_harmonic(self, freqs_cm1: Iterable[float]) -> float:
        return sum(self.calculate_entropy_single_mode_harmonic(freq) for freq in freqs_cm1 if float(freq) > 0.0)

    def calculate_entropy_qh(self, freqs_cm1: Iterable[float]) -> float:
        corrected = 0.0
        for freq in freqs_cm1:
            frequency = float(freq)
            if frequency <= 0.0:
                continue
            s_harm = self.calculate_entropy_single_mode_harmonic(frequency)
            s_rot = self.calculate_entropy_single_mode_rotor(frequency)
            weight = grimme_weighting(frequency, Bcut=self.grimme_cutoff)
            corrected += (weight * s_harm) + ((1.0 - weight) * s_rot)
        return corrected

    def calculate_barrier_free_energy(
        self,
        enthalpy_activation_kcal_mol: float,
        reactant_freqs_cm1: Iterable[float],
        ts_freqs_cm1: Iterable[float],
        *,
        use_quasi_harmonic: bool,
    ) -> float:
        entropy_fn = self.calculate_entropy_qh if use_quasi_harmonic else self.calculate_entropy_harmonic
        delta_s = entropy_fn(ts_freqs_cm1) - entropy_fn(reactant_freqs_cm1)
        return float(enthalpy_activation_kcal_mol) - (self.temperature_k * delta_s)


__all__: List[str] = ["QuasiHarmonicCorrector", "grimme_weighting"]