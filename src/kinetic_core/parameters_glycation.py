"""
src/kinetic_core/parameters_glycation.py -- THE GLYCATION ARM (Build Wave B20, 2026-09-09).

Five steps on the trunk lane (network.GLYCATION_REACTIONS) that make the protein's bound lysine a
reactant: glucose + bound lysine -> bound fructosyl-lysine (k_glyc, second order); the bound Amadori
compound -> CML (k_flp_cml), -> CEL via methylglyoxal lumped (k_flp_cel), -> 3-deoxyglucosone + the
lysine back (k_flp_decay, the dominant loss); CML -> the melanoidin pools (k_cml_loss). The bound
lysine pool comes from the spec's protein loading and the matrix table's amine density
(data/species/protein_matrices.yml) times the declared available fraction, the centre of the
table's availability band; with no loading the pool is zero and the arm is inert.

WHAT IS FITTED (wave B20, pre-registered in results/validation/kinetic_core_b20_prereg.md): the
five log10 rate constants at the trunk's 100 C reference, against the ten rate constants Nguyen
2016 printed for sodium caseinate 30 g/L + glucose 150 mmol/L in 0.1 M phosphate pH 6.8 at 120 and
130 C (nguyen2016_extraction.md Table 1, system M1), each row weighted by its printed interval. The
fitted values below are FROZEN LITERALS asserted equal to the B20 fit report by
tests/unit/test_kinetic_core_b20.py (the B7 / B18 pattern).

WHAT IS DECLARED. The barriers: Nguyen's two temperatures are ten degrees apart with overlapping
intervals and give Q10 values from 0.3 to 2.9, which the authors themselves call apparent, so no
barrier is fitted. Each step takes a MEASURED barrier from the nearest measured step: the
glycation the trunk's Amadori formation (Martins 2005, 96.8 kJ/mol); the Amadori decay the trunk's
Amadori -> 3-deoxyglucosone (97.1); CML formation Berk 2021's fructosyl-lysine -> CML (113 kJ/mol,
three temperatures in sesame); CEL formation Berk 2021's methylglyoxal + bound lysine -> CEL (92);
the CML loss a flat barrier (Nguyen's pair falls with temperature, which no barrier reproduces;
flagged, the Kocadagli glyoxal-sink precedent). The glyoxal -> CML route is NOT written: it fits to
zero in Nguyen 2016, Berk 2021 and Hamzalioglu 2026. The available fraction of the lysine sites is
the matrix table's band centre, and the band travels on every answer as a declared interval.
"""
from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .parameters import AW_OF_MEASUREMENT, KineticParameter

#: The trunk's reference temperature for the stored constants (parameters.T_REF_K is 373.15 K).
GLYCATION_FIT_PH = 6.8
#: Declared barriers, kJ/mol, each a MEASURED value from the nearest measured step (module docstring).
EA_GLYC_KJ_MOL = 96.8          # Martins 2005 k_schiff (glucose + glycine -> Amadori), the trunk's own
EA_FLP_DECAY_KJ_MOL = 97.1     # Martins 2005 k_ama_tdg (Amadori -> 3-deoxyglucosone), the trunk's own
EA_FLP_CML_KJ_MOL = 113.0      # Berk 2021 k8, fructosyl-lysine -> CML, 180-220 C, R2 0.966
EA_FLP_CEL_KJ_MOL = 92.0       # Berk 2021 k15, methylglyoxal + bound lysine -> CEL, R2 0.999
EA_CML_LOSS_KJ_MOL = 0.0       # declared flat: Nguyen's k11 pair (0.29 -> 0.077 /min from 120 to 130 C) has no barrier
#: The molar masses the concentration report uses (g/mol): lysine for the bound residue, free
#: fructosyl-lysine, CML, CEL.
MOLAR_MASS_G_PER_MOL: Mapping[str, float] = {"LYSP": 146.19, "FLP": 308.33, "CML": 204.22, "CEL": 218.25}

# ---------------------------------------------------------------------------
# THE FROZEN B20 VALUES (the B20 fit report's optimum, 2026-09-09: cost 17.2 on 10 rows, reduced
# chi-square 3.4, every coordinate identified; asserted against the report by the unit test).
# Prior centres: Nguyen 2016 M1 brought from 120 / 130 C to 100 C with the declared barriers
# (k3 1.5e-4 L/(mmol min); k7 8.8e-3, k9 2.3e-3, k8 5.2e-2, k11 0.29 / 0.077 per minute).
# ---------------------------------------------------------------------------
FROZEN_B20: Mapping[str, float] = {
    "log10_k_glyc_100C": -4.707341524577862,
    "log10_k_flp_cml_100C": -3.317402852045118,
    "log10_k_flp_cel_100C": -3.652055191735323,
    "log10_k_flp_decay_100C": -1.856048900182166,
    "log10_k_cml_loss_100C": -0.9827379909354617,
}

_NGUYEN = ("Nguyen, van der Fels-Klerx & van Boekel 2016, Food Chem. 192:125 (doi 10.1016/j.foodchem.2015.06.110), "
           "Table 1 system M1: sodium caseinate 30 g/L (about 16 mmol/L lysine residues) + glucose 150 mmol/L, 0.1 M "
           "phosphate pH 6.8, 120 and 130 C, multiresponse fit with 95 % HPD intervals")
_NGUYEN_DOSSIER = "nguyen2016_extraction.md sec. 4 (the rate table); results/validation/kinetic_core_b20_prereg.md"
_BERK = ("Berk, Hamzalioglu & Gokmen 2021 (sesame seed, 180-220 C): k8 fructosyl-lysine -> CML Ea 113 kJ/mol; k15 "
         "methylglyoxal + bound lysine -> CEL Ea 92 kJ/mol; berk2021_extraction.md sec. 4")

GLYCATION_AVAILABILITY_CAVEAT = (
    "GLYCATION (B20): the bound-lysine pool is the matrix's amine density times the centre of its declared "
    "availability band; the band's corners are carried as the interval on every CML, CEL and fructosyl-lysine "
    "answer. The rates are one laboratory's (casein in water at 120-130 C) and the barriers are declared "
    "from the nearest measured step; the CML loss has no barrier."
)
GLYCATION_NO_PROTEIN_REASON = (
    "GLYCATION TARGETS need a protein: CML, CEL and fructosyl-lysine are made on protein-bound lysine, "
    "and this spec states no protein loading (protein_g_per_l with a matrix on file, or protein_sites "
    "with an amine density). Free lysine is not this arm's substrate (it resolves to the acrylamide lane "
    "as an amine). Refused rather than answered with a structural zero."
)


def _param(key: str, transformation: str, log10_k: float, ea: float, order: int, unit: str, note: str,
           flags: Tuple[str, ...] = ()) -> KineticParameter:
    return KineticParameter(
        key=key, transformation=transformation, k_ref=10.0 ** float(log10_k), ea_kj_mol=float(ea),
        unit=unit, order=order, evidence_class="derived_from_fit_data",
        source_anchor=_NGUYEN, dossier_anchor=_NGUYEN_DOSSIER,
        conditions="water, 0.1 M phosphate pH 6.8, 120-130 C, casein-bound lysine; the barrier is declared (see flags)",
        ph_of_measurement=GLYCATION_FIT_PH, temperature_range_c=(120.0, 130.0),
        rate_transfer="licensed_at_measurement_ph_only", aw_of_measurement=AW_OF_MEASUREMENT,
        flags=("b20_glycation", "fitted_wave_b20", "barrier_declared_not_fitted") + flags, note=note,
    )


def with_fitted_glycation(log10_k_glyc: float, log10_k_flp_cml: float, log10_k_flp_cel: float,
                          log10_k_flp_decay: float, log10_k_cml_loss: float) -> Dict[str, KineticParameter]:
    """The glycation block at arbitrary log10 constants (the fit generator's hook and the report reader's)."""
    return {
        "k_glyc": _param("k_glyc", "glucose + bound lysine -> bound fructosyl-lysine (Schiff base + Amadori, lumped)",
                         log10_k_glyc, EA_GLYC_KJ_MOL, 2, "L/(mmol*min)",
                         "B20 fit rows: Nguyen 2016 k3 at 120 and 130 C. Barrier: the trunk's Amadori formation (Martins 2005).",
                         ("barrier_from_martins_k_schiff",)),
        "k_flp_cml": _param("k_flp_cml", "bound fructosyl-lysine -> CML + C4 fragments (oxidative cleavage)",
                            log10_k_flp_cml, EA_FLP_CML_KJ_MOL, 1, "1/min",
                            "B20 fit rows: Nguyen 2016 k7 at 120 and 130 C. Barrier: " + _BERK + " (k8).",
                            ("barrier_from_berk2021",)),
        "k_flp_cel": _param("k_flp_cel", "bound fructosyl-lysine -> CEL + C3 fragments (via methylglyoxal, lumped)",
                            log10_k_flp_cel, EA_FLP_CEL_KJ_MOL, 1, "1/min",
                            "B20 fit rows: Nguyen 2016 k9 at 120 and 130 C (the 120 C interval spans zero). Barrier: " + _BERK + " (k15).",
                            ("barrier_from_berk2021",)),
        "k_flp_decay": _param("k_flp_decay", "bound fructosyl-lysine -> 3-deoxyglucosone + bound lysine (Amadori decay)",
                              log10_k_flp_decay, EA_FLP_DECAY_KJ_MOL, 1, "1/min",
                              "B20 fit rows: Nguyen 2016 k8 ('AP -> MRPs') at 120 and 130 C. Barrier: the trunk's Amadori -> 3-DG (Martins 2005).",
                              ("barrier_from_martins_k_ama_tdg",)),
        "k_cml_loss": _param("k_cml_loss", "CML -> melanoidin pools (loss)",
                             log10_k_cml_loss, EA_CML_LOSS_KJ_MOL, 1, "1/min",
                             "B20 fit rows: Nguyen 2016 k11 at 120 and 130 C, which FALL with temperature; the barrier is declared "
                             "flat (the Kocadagli glyoxal-sink precedent), flagged.",
                             ("ea_declared_zero", "rows_fall_with_temperature")),
    }


GLYCATION_PARAMETERS: Mapping[str, KineticParameter] = with_fitted_glycation(
    FROZEN_B20["log10_k_glyc_100C"], FROZEN_B20["log10_k_flp_cml_100C"], FROZEN_B20["log10_k_flp_cel_100C"],
    FROZEN_B20["log10_k_flp_decay_100C"], FROZEN_B20["log10_k_cml_loss_100C"],
)
GLYCATION_KEYS: Tuple[str, ...] = tuple(GLYCATION_PARAMETERS)
#: The order of the fitted coordinates in the report block ``glycation``.
GLYCATION_COORDINATES: Tuple[str, ...] = tuple(FROZEN_B20)


def available_fraction(band: Tuple[float, float]) -> float:
    """The declared available fraction of the matrix's lysine sites: the centre of its stated band."""
    lo, hi = float(band[0]), float(band[1])
    return 0.5 * (lo + hi)


GLYCATION_WISHLIST: Mapping[str, str] = {
    "k_glyc": "glucose + bound lysine in a PLANT isolate (pea, soy) at two temperatures in water: Nguyen 2016 is casein",
    "k_flp_cml": "a barrier in water: Berk 2021's 113 kJ/mol is a dry seed at 180-220 C; Nguyen's two points give none",
    "k_cml_loss": "any CML loss measured at more than one temperature; the pair on disk falls with temperature",
    "availability": "the fraction of an isolate's lysine that glycates (furosine after a defined cook), per isolate",
}
