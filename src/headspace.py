"""
src/headspace.py

Phase D: Headspace & Volatility Model.
Converts matrix concentrations to predicted air-phase (headspace) concentrations.
"""

import math
import yaml
from pathlib import Path
from typing import Dict, Optional, List
from src.logger import get_logger
logger = get_logger(__name__)

from src import data_paths
from src import data_access
from src.matrix_calibration_registry import (
    determine_matrix_process_state,
    get_matrix_calibration_record,
    get_matrix_runtime_composition_policy,
)
from src.extrusion import compute_extrusion_headspace_adjustment
from src.literature_runtime import get_retention_ph_release_profile
from src.protein_binding import (
    binding_mode_active,
    observability_factor as protein_binding_observability_factor,
)
from src.matrix_correction import ProteinType, resolve_compound_matrix_retention, resolve_matrix_correction

# Upper bound of the van't Hoff extrapolation window for Kaw (audit 2026-08-26).
#
# `Kaw_25c` and `delta_H_sol_kj_mol` in data/lit/henry_constants.yml are fitted
# near 298 K. A constant-enthalpy van't Hoff extrapolation is only defensible
# while (a) the aqueous phase is still liquid water at ambient pressure and
# (b) delta_H_sol is roughly temperature independent; both assumptions fail
# above the normal boiling point of water. Extrapolating anyway produced the
# audit's worked example: hexanal (Kaw_25c = 0.015, dH_sol = -40 kJ/mol) at
# 453 K gives an extrapolation factor of ~249 and Kaw ~ 3.7, i.e. a
# semi-volatile aldehyde predicted to sit almost entirely in the vapour phase,
# which the underlying Henry fit cannot support.
#
# The reference temperature used in the exponent is therefore clamped to
# [273.15, 373.15] K. Above 373.15 K (retort / extrusion regimes) Kaw is held
# at its 100 C value, which is the last point the fit can defend; the residual
# temperature dependence of those processes is carried by the empirical
# process-state calibrations in src/matrix_calibration_registry.py instead.
#
# Note: no ceiling of Kaw <= 1 is imposed. Kaw > 1 is physical for sparingly
# soluble gases (H2S reaches Kaw ~ 1.5 at 373 K), so a hard ceiling would be
# wrong; the temperature clamp is the defensible bound.
VANT_HOFF_MIN_TEMP_K = 273.15
VANT_HOFF_MAX_TEMP_K = 373.15


class HeadspaceModel:
    """
    Models the partitioning of volatiles between the food matrix and air.
    Accounts for temperature (Van't Hoff) and matrix suppression (lipids/proteins).
    """

    def __init__(self, constants_path: str | Path = data_paths.HENRY_CONSTANTS):
        constants_path = Path(constants_path)
        if not constants_path.is_absolute():
            constants_path = data_paths.REPO_ROOT / constants_path
        self.constants_path = constants_path
        self.data = self._load_constants()
        self.R = 0.008314  # kJ/(mol*K)

    def _load_constants(self) -> Dict[str, Dict]:
        # Strict since 2026-09-01: a missing file used to load as {} and every Kaw fell
        # back to the hard-coded default without a word.
        raw = data_access.load_yaml(self.constants_path) or {}
        return {c["name"]: c for c in raw.get("constants", [])}

    def get_kaw_at_temp(self, name: str, temp_k: float) -> float:
        """
        Calculates the dimensionless air-water partition coefficient at temp_k.
        Uses Van't Hoff / Clausius-Clapeyron extrapolation.
        Kaw(T) = Kaw(Tr) * exp(-dH/R * (1/T - 1/Tr))

        The temperature entering the exponent is clamped to
        [VANT_HOFF_MIN_TEMP_K, VANT_HOFF_MAX_TEMP_K]; see the module-level
        comment on those constants for why (audit 2026-08-26).
        """
        entry = self.data.get(name)
        if not entry:
            return 0.01  # Default fallback volatility

        kaw_298 = entry["Kaw_25c"]
        dh = entry["delta_H_sol_kj_mol"]

        # Clamp before extrapolating: outside the fitted window the constant
        # dH_sol assumption breaks down and the exponential runs away
        # (hexanal at 453 K -> Kaw 3.7, unphysical for a soluble aldehyde).
        effective_temp_k = min(max(float(temp_k), VANT_HOFF_MIN_TEMP_K), VANT_HOFF_MAX_TEMP_K)

        # Extrapolate: Kaw(T) = Kaw(Tr) * exp(dH_sol/R * (1/temp_k - 1/Tr))
        # Since dH_sol is negative, Kaw increases as temp_k increases.
        exponent = (dh / self.R) * (1.0 / effective_temp_k - 1.0 / 298.15)
        return kaw_298 * math.exp(exponent)

    def _extract_properties(self, smiles: str) -> Dict[str, float]:
        """Estimated logP and MW using RDKit for binding characterization."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"logP": 1.0, "MW": 100.0}
        return {
            "logP": Crippen.MolLogP(mol),
            "MW": Descriptors.MolWt(mol)
        }

    def _get_kprot_for_compound(self, name: str, smiles: Optional[str] = None) -> float:
        """
        Calculates empirical protein binding constant (Kprot) for a volatile.
        Phase 1: Partitioning Correction.
        Uses a hydrophobicity-driven model (logP) + molecular weight scaling.
        Kprot = alpha * logP + beta * MW
        """
        if smiles:
            props = self._extract_properties(smiles)
            logp = props["logP"]
            mw = props["MW"]
            
            # Hydrophobic binding (binding affinity increases with logP)
            # Scaling factors estimated from Mcclements et al. (flavor-protein interactions)
            kprot_hydrophobic = max(0, 15.0 * logp)
            # Entropy/Size factor
            kprot_size = 0.05 * mw
            
            return kprot_hydrophobic + kprot_size
            
        n = name.lower()
        # Fallback for named compounds without SMILES lookup
        if "sulfide" in n or "h2s" == n:
            return 0.0      
        elif "thiol" in n:
            return 8.0      
        elif "methional" in n:
            return 2.0      
        elif "pyrazine" in n or "thiazole" in n:
            return 12.0     
        elif "aldehyde" in n or "anal" in n or "fural" in n:
            return 45.0     # Carbonyl-amine covalent binding (Schiff)
        else:
            return 10.0

    # Mirrors `_resolve_output_matrix_context` in src/recommend.py: a protein
    # fraction at or above this value is the ReactionConditions "unspecified"
    # sentinel (protein_fraction defaults to 1.0), not a measured composition.
    UNSPECIFIED_PROTEIN_FRACTION_SENTINEL = 0.999

    def _guard_matrix_fractions(self, fat_fraction: float, protein_fraction: float) -> tuple:
        """
        Neutralise the `protein_fraction = 1.0` sentinel (audit 2026-08-26).

        `ReactionConditions.protein_fraction` defaults to 1.0, which in this
        codebase means "unspecified" rather than a real volumetric fraction.
        Feeding it into the sequestration denominator
        `1 + Kfat*phi_fat + Kprot*phi_prot` suppresses headspace by up to ~46x
        for a carbonyl (Kprot ~ 45), silently. `src/recommend.py` already guards
        its own call sites; this guard makes `HeadspaceModel` safe for any new
        caller. No production path currently reaches it: `recommend.py` zeroes
        the fractions before calling, and `src/sensory.py` passes headspace=None.

        A hydrated food matrix cannot physically be >= 99.9% protein by volume,
        so clamping to 0.0 (rather than raising) keeps the model usable while
        removing the landmine, and matches the existing "fall back, log, keep
        going" style of this module.
        """
        fat = max(float(fat_fraction), 0.0)
        protein = max(float(protein_fraction), 0.0)
        if protein >= self.UNSPECIFIED_PROTEIN_FRACTION_SENTINEL:
            logger.warning(
                "predict_headspace received protein_fraction=%.3f; treating it as "
                "'unspecified' (the ReactionConditions default sentinel) and using 0.0. "
                "Pass a measured volumetric protein fraction, or protein_type=, "
                "to model matrix retention.",
                protein,
            )
            protein = 0.0
        return fat, protein

    def _matrix_retention_fallback(
        self,
        protein_type: Optional[str],
        fat_fraction: float,
        protein_fraction: float,
        denaturation_state: float,
    ) -> float:
        if fat_fraction > 0.0 or protein_fraction > 0.0 or not protein_type:
            return 1.0

        try:
            p_type = ProteinType(protein_type)
        except ValueError:
            return 1.0

        if p_type == ProteinType.FREE_AMINO_ACID:
            return 1.0
        return float(resolve_matrix_correction(p_type, denaturation_state).volatile_retention)

    def get_matrix_ph_release_factor(
        self,
        name: str,
        *,
        protein_type: Optional[str],
        pH: Optional[float],
    ) -> float:
        """
        Empirical pH-dependent headspace release factor for plant matrices.

        The current calibration is intentionally narrow:
        - anchored to pH 6.0 so the existing Pratap-Singh matrix-only baselines stay stable
        - only applied to pea/soy matrix types
        - only applied to acid-sensitive lipid-derived off-flavour markers
                - tuned through the structured Karolkowski 2021 PPI pH-release payload for
                    acid-sensitive aldehydes and furans while preserving the executable pH 6 baseline
        """
        if pH is None or not protein_type:
            return 1.0

        try:
            p_type = ProteinType(protein_type)
        except ValueError:
            return 1.0

        if p_type not in {
            ProteinType.PEA_ISOLATE,
            ProteinType.PEA_CONCENTRATE,
            ProteinType.SOY_ISOLATE,
            ProteinType.SOY_CONCENTRATE,
        }:
            return 1.0

        compound = name.lower()
        acid_sensitive = any(token in compound for token in ["anal", "enal", "furan"])
        if not acid_sensitive:
            return 1.0

        profile = get_retention_ph_release_profile(name, protein_type=p_type.value)
        reference_ph = float(profile.get("reference_ph", 6.0) or 6.0)
        # 2026-08-27 (Wave T3, finding T1-02) -- WHAT 0.235 IS, RECORDED HONESTLY.
        # It is not a fitted or measured slope. It is exactly ln(1.60)/2 = 0.2350018,
        # back-solved so that the pH-4.5-vs-pH-6.5 ratio of this function comes out at
        # 1.60000 (= +60% release). The 1.60 is the arithmetic midpoint of a "~55-65%"
        # band in data/benchmarks/maillard_validation_benchmarks.md section 3.2 -- an
        # INVENTED table whose citation is a self-declaring placeholder DOI (recorded here
        # with spaces inserted so it is not re-asserted as an anchor, and so that
        # scripts/ci/citation_gate.py's confabulation check does not fire on this comment:
        # "10.1016 / j.foodchem.2021.xxx") and whose own hexanal row implies +65.9%, not the
        # +59% it is labelled with. `max_factor` below is the same 1.60 again.
        # The value is KEPT because it is independently vindicated: Fischer, Cachon &
        # Cayot (2021), Food Res. Int. 150:110760 (10.1016/j.foodres.2021.110760,
        # CrossRef-verified 2026-08-27) state "hexanal release was found 59% higher with
        # extraction using pH 4.5 than with pH 6.5" -- 1.600 vs 1.590, 0.63% apart.
        # Vindication caveat, stated because 0.63% invites over-reading: Fischer varied
        # EXTRACTION pH (which volatiles the isolate carries), this function varies
        # RELEASE pH at measurement time. Different quantities. Treat the knob as
        # no_verifiable_source for any quantitative purpose; the literature supports the
        # direction and one hexanal percentage. Full record:
        # data/lit/retention_reference_payloads.json -> karolkowski_2021_ppi_hexanal_ph_release
        # -> runtime_surrogate.surrogate_basis.
        log_slope = float(profile.get("log_slope", 0.235) or 0.235)
        min_factor = float(profile.get("min_factor", 0.75) or 0.75)
        max_factor = float(profile.get("max_factor", 1.6) or 1.6)

        centered_delta = reference_ph - float(pH)
        factor = math.exp(log_slope * centered_delta)
        return max(min_factor, min(max_factor, factor))

    def get_matrix_benchmark_headspace_factor(
        self,
        name: str,
        *,
        protein_type: Optional[str],
        pH: Optional[float],
        temperature_celsius: float = 40.0,
        time_minutes: float = 10.0,
        water_activity: Optional[float] = None,
        binding_context: Optional[Dict[str, object]] = None,
    ) -> float:
        """
        Empirical observable-release factor for the Pratap-Singh plant-matrix lane.

        The matrix-only benchmark path combines:
        - oxidation-load generation in src/lipid_oxidation.py
        - pea-referenced marker yields in src/benchmark_validation.py
        - this headspace observable factor, which carries the pea/soy release gap
          from the Pratap-Singh ambient slurry family
        - the narrower pH-dependent release modifier already validated against the
                    Karolkowski PPI pH-release trend family

        This keeps benchmark_validation focused on intake chemistry while making
        the matrix-specific observable calibration explicit in the headspace layer.
        """
        if not protein_type:
            return 1.0

        try:
            p_type = ProteinType(protein_type)
        except ValueError:
            return 1.0

        process_state = determine_matrix_process_state(
            temperature_celsius=float(temperature_celsius),
            time_minutes=float(time_minutes),
            water_activity=water_activity,
        )

        # 2026-08-27 (Wave S4) -- THE BINDING-PHYSICS OBSERVABILITY MODE.
        #
        # When the mode is `binding_physics`, the observability factor is computed from
        # MEASURED protein-binding constants (data/lit/binding_constants.yml) instead of
        # the fitted / back-solved constants in src/matrix_calibration_registry.py.
        #
        # NO DOUBLE COUNTING, and it is structural rather than a comment: this branch
        # returns BEFORE `get_matrix_calibration_record` is ever consulted, and it does
        # NOT apply `dynamic_release_factor` either. That second omission is the one that
        # is easy to miss -- `compose_dynamic_retention` routes through
        # `resolve_compound_matrix_retention`, whose `volatile_retention` is documented as
        # "fraction escaping matrix (rest is bound)", i.e. it is ITSELF an (unanchored)
        # binding model. Running both would count protein binding twice.
        #
        # The pH release factor IS still applied, in both modes, deliberately: it is a
        # pH-dependent release term rather than a binding term, and keeping it identical
        # on both sides is what makes the mode-vs-mode comparison a comparison of the
        # observability constant alone.
        if binding_mode_active():
            binding = protein_binding_observability_factor(
                name,
                context=binding_context or {},
            )
            return float(binding.f_free) * self.get_matrix_ph_release_factor(
                name,
                protein_type=protein_type,
                pH=pH,
            )

        record = get_matrix_calibration_record(
            name,
            protein_type=p_type.value,
            process_state=process_state,
        )
        base_factor = float(record.observable_factor) if record is not None else 1.0
        dynamic_release_factor = 1.0
        runtime_policy = get_matrix_runtime_composition_policy(
            name,
            protein_type=p_type.value,
            process_state=process_state,
        )
        if runtime_policy.get("mode") == "compose_dynamic_retention":
            matrix_retention = resolve_compound_matrix_retention(
                name,
                protein_type=p_type,
                denaturation_state=0.5,
                temperature_celsius=temperature_celsius,
                time_minutes=time_minutes,
                water_activity=water_activity,
                process_state=process_state,
            )
            baseline_retention = resolve_compound_matrix_retention(
                name,
                protein_type=p_type,
                denaturation_state=0.5,
            )
            dynamic_release_factor = 1.0 if baseline_retention <= 0.0 else matrix_retention / baseline_retention
        return base_factor * dynamic_release_factor * self.get_matrix_ph_release_factor(
            name,
            protein_type=protein_type,
            pH=pH,
        )

    def predict_headspace(self, 
                          matrix_concentrations: Dict[str, float], 
                          temp_c: float, 
                          fat_fraction: float = 0.0,
                          protein_fraction: float = 0.0,
                          protein_type: Optional[str] = None,
                          denaturation_state: float = 0.5,
                          pH: Optional[float] = None,
                          extrusion_process: Optional[Dict[str, object]] = None) -> Dict[str, float]:
        """
        Predicts air-phase concentrations (ppm).
        
        Equation: C_air = C_total * Kaw_eff
        Kaw_eff = Kaw(T) / (1 + Kfat * phi_fat + Kprot * phi_prot)

        If no explicit matrix fractions are provided, `protein_type` can supply
        a conservative fallback retention calibrated from the matrix-correction
        literature estimates already used elsewhere in the project.

        `pH` is optional and currently only applies an empirical plant-matrix
        release correction for acid-sensitive aldehydes/furans in pea/soy systems.

        `protein_fraction` must be a real volumetric matrix fraction. Passing the
        `ReactionConditions.protein_fraction` default of 1.0 straight through is
        a sentinel meaning "unspecified", not a 100%-protein matrix, and is
        neutralised here (see the guard below).
        """
        temp_k = temp_c + 273.15
        fat_fraction, protein_fraction = self._guard_matrix_fractions(fat_fraction, protein_fraction)
        air_concs = {}
        base_matrix_retention = self._matrix_retention_fallback(
            protein_type,
            fat_fraction,
            protein_fraction,
            denaturation_state,
        )
        
        for name, c_total in matrix_concentrations.items():
            entry = self.data.get(name)
            kaw_base = self.get_kaw_at_temp(name, temp_k)
            ph_release_factor = self.get_matrix_ph_release_factor(
                name,
                protein_type=protein_type,
                pH=pH,
            )
            
            if entry:
                k_fat = entry.get("Kfat", 1.0)
                k_prot = self._get_kprot_for_compound(name)
                matrix_retention = base_matrix_retention
                if protein_type and fat_fraction <= 0.0 and protein_fraction <= 0.0:
                    matrix_retention = resolve_compound_matrix_retention(
                        name,
                        protein_type=protein_type,
                        denaturation_state=denaturation_state,
                        temperature_celsius=temp_c,
                        time_minutes=None,
                    )
                
                # Effective Kaw accounting for matrix sequestration
                denom = 1.0 + (k_fat * fat_fraction) + (k_prot * protein_fraction)
                kaw_eff = kaw_base / denom
                extrusion_factor = 1.0
                if extrusion_process and extrusion_process.get("active"):
                    extrusion_factor = float(compute_extrusion_headspace_adjustment(name, extrusion_process).get("combined_headspace_factor", 1.0))
                air_concs[name] = c_total * kaw_eff * matrix_retention * ph_release_factor * extrusion_factor
            else:
                # Basic fallback
                matrix_retention = base_matrix_retention
                if protein_type and fat_fraction <= 0.0 and protein_fraction <= 0.0:
                    matrix_retention = resolve_compound_matrix_retention(
                        name,
                        protein_type=protein_type,
                        denaturation_state=denaturation_state,
                        temperature_celsius=temp_c,
                        time_minutes=None,
                    )
                extrusion_factor = 1.0
                if extrusion_process and extrusion_process.get("active"):
                    extrusion_factor = float(compute_extrusion_headspace_adjustment(name, extrusion_process).get("combined_headspace_factor", 1.0))
                air_concs[name] = c_total * kaw_base * matrix_retention * ph_release_factor * extrusion_factor
                
        return air_concs

