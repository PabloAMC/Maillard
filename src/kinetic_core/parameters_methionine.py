"""
src/kinetic_core/parameters_methionine.py -- THE METHIONINE CHAIN (Build Wave B22, 2026-09-09).

Four steps on the trunk lane (network.METHIONINE_REACTIONS): the Strecker deamination of methionine
on glyoxal and on methylglyoxal (the B18 step with methionine as the amino acid: the aminoketone is
the dicarbonyl's, glycine's AKG / AKM by construction; the aldehyde is methional), the retro-Michael
release of methanethiol from methional, and the disulfide on an apparent second-order constant.

WHAT IS FITTED (wave B22, pre-registered in results/validation/kinetic_core_b22_prereg.md): the
identity ratio of methionine to glycine on the two Strecker steps (one log10 factor applied to both
B18 constants, B18's barriers and pH term kept); log10 k_mtal_msh at 100 C and its barrier (free,
20-150 kJ/mol); log10 k_msh_dmds at 100 C with its barrier declared from Pan 2025's apparent value
for the disulfide (81 kJ/mol). Rows: Pan 2025 Table 2's nine zero-order constants (methional,
methanethiol, dimethyl disulfide at 100 / 120 / 140 C), each modelled as the mean formation rate
over 30-600 s in Pan's pot. The fitted values below are FROZEN LITERALS asserted equal to the B22
fit report by tests/unit/test_kinetic_core_b22.py.

WHAT IS DECLARED. (i) The unit of Pan's constants, micromoles per litre per second, is inferred
(the paper prints none; its two endpoint levels fix it) -- the wave stands or falls with it.
(ii) Sucrose (44 mmol/L in Pan's pot) is omitted: the trunk has no sucrose. (iii) The identity
ratio is conditional on the trunk's dicarbonyl supply in Pan's pot. (iv) The disulfide constant is
apparent and pot-specific. (v) Methionine is charged, declared, as glycine at the same molarity for
the Amadori chemistry. Dimethyl trisulfide is not written: it needs hydrogen sulfide, which the
sugar path does not carry.
"""
from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .parameters import AW_OF_MEASUREMENT, KineticParameter
from .parameters_pyrazine import FROZEN_B18

METHIONINE_FIT_PH = 6.2
#: Pan 2025's apparent barrier for the disulfide (three temperatures, unit-free), declared for k_msh_dmds.
EA_MSH_DMDS_KJ_MOL = 81.0
#: The band the release step's barrier is fitted inside.
EA_MTAL_MSH_BAND_KJ_MOL: Tuple[float, float] = (20.0, 150.0)
MOLAR_MASS_G_PER_MOL: Mapping[str, float] = {"MET": 149.21, "MTAL": 104.17, "MSH": 48.11, "DMDS": 94.20}

# ---------------------------------------------------------------------------
# THE FROZEN B22 VALUES. Replaced by the fit; asserted against the report.
# Prior centres: identity ratio 1 (log10 0); k_mtal_msh from Pan's methanethiol / methional rate ratio at
# 120 C (0.14) over a methional lifetime of minutes; k_msh_dmds from the disulfide / methanethiol ratio.
# ---------------------------------------------------------------------------
#: THE B22 FIT REPORT'S OPTIMUM, kept as the record and asserted against the report by the unit test. It is NOT
#: installed: the wave did not ship (kinetic_core_b22_prereg.md sec. 6). The identity ratio sits on its
#: ceiling (a hundred times glycine's constants) and Pan's methional rates are still 3.6 to 5.6 decades
#: below the printed ones while Deng 2022's pot comes out 1.5 to 2.9 decades too high: methional does not
#: track free dicarbonyl times methionine with glycine's constants.
FROZEN_B22: Mapping[str, float] = {
    "log10_identity_ratio_met_over_gly": 1.9999999999999998,
    "log10_k_mtal_msh_100C": -0.2814742548055645,
    "ea_mtal_msh_kj_mol": 20.000000254221657,
    "log10_k_msh_dmds_100C": 1.9999999999999998,
}
#: The wave's verdict. While False, the operative set carries the four constants at ZERO (the steps exist
#: and carry no flux) and a request for methional, methanethiol from methional or dimethyl disulfide is
#: refused by name with the verdict; the B17 precedent.
METHIONINE_SHIPPED = False
#: The inert values the operative set carries while the wave is not shipped (a ratio of zero, zero rates).
INERT_B22: Mapping[str, float] = {
    "log10_identity_ratio_met_over_gly": -300.0,
    "log10_k_mtal_msh_100C": -300.0,
    "ea_mtal_msh_kj_mol": 64.0,
    "log10_k_msh_dmds_100C": -300.0,
}
METHIONINE_COORDINATES: Tuple[str, ...] = tuple(FROZEN_B22)

_PAN = ("Pan et al. 2025 (methionine 0.268 mmol/L + fructose 111 + glucose 83 + sucrose 44 mmol/L, 50 mmol/L citrate pH 6.2, "
        "100 / 120 / 140 C, 30-600 s; Table 2 zero-order constants, unit inferred as umol L-1 s-1 from the printed endpoints); "
        "pan2025_extraction.md sec. 4")
METHIONINE_CAVEAT = (
    "METHIONINE CHAIN (B22): methional comes from methionine's Strecker step on the trunk's glyoxal and methylglyoxal at "
    "an identity ratio to glycine FITTED on one laboratory's zero-order rates in a fruit-sugar pot (Pan 2025), whose unit "
    "the paper does not print (inferred from its endpoints); the ratio is conditional on the trunk's dicarbonyl supply "
    "in that pot. Methanethiol is the release from methional; the disulfide runs on an APPARENT constant with no oxidant "
    "tracked. Dimethyl trisulfide is not made here: it needs hydrogen sulfide, which the sugar path does not carry."
)
METHIONINE_NOT_SHIPPED_REASON = (
    "the methionine chain (wave B22) did not ship: fitted on Pan 2025's rates, the identity ratio to glycine ran to its "
    "ceiling (a hundredfold) and methional was still 3.6 to 5.6 decades below the printed rates in that pot while it came "
    "out 1.5 to 2.9 decades too high in Deng 2022's methionine + glucose pot; methional does not form as free dicarbonyl "
    "times methionine with glycine's Strecker constants (kinetic_core_b22_prereg.md sec. 6). The steps stay in the "
    "network at zero and the targets are refused rather than answered with a structure the data refuted."
)
METHIONINE_NO_DMTS_REASON = (
    "dimethyl trisulfide needs hydrogen sulfide (Chin & Lindsay 1994 detect none without it), which the sugar path does "
    "not carry; the sulfur lane makes hydrogen sulfide from cysteine but has no methionine. Refused rather than invented."
)


def _p(key, transformation, k_ref, ea, order, unit, note, flags):
    return KineticParameter(
        key=key, transformation=transformation, k_ref=float(k_ref), ea_kj_mol=float(ea), unit=unit, order=order,
        evidence_class="derived_from_fit_data", source_anchor=_PAN,
        dossier_anchor="pan2025_extraction.md sec. 4; deng2022_extraction.md; results/validation/kinetic_core_b22_prereg.md",
        conditions="water, citrate pH 6.2, 100-140 C, methionine 0.27 mmol/L in a 240 mmol/L hexose pool",
        ph_of_measurement=METHIONINE_FIT_PH, temperature_range_c=(100.0, 140.0), rate_transfer="licensed_at_measurement_ph_only",
        aw_of_measurement=AW_OF_MEASUREMENT, flags=("b22_methionine", "fitted_wave_b22") + tuple(flags), note=note,
    )


def with_fitted_methionine(log10_ratio: float, log10_k_mtal_msh: float, ea_mtal_msh: float, log10_k_msh_dmds: float) -> Dict[str, KineticParameter]:
    """The methionine block at arbitrary values (the fit generator's hook and the report reader's)."""
    def _pow(v: float) -> float:
        # the inert record (log10 -300) must be an exact zero, not a denormal: "no flux" means none
        return 0.0 if float(v) < -100.0 else 10.0 ** float(v)

    r = _pow(log10_ratio)
    k_go = r * 10.0 ** FROZEN_B18["log10_k_go_ak_100C"]
    k_mgo = r * 10.0 ** FROZEN_B18["log10_k_mgo_ak_100C"]
    return {
        "k_go_met": _p("k_go_met", "glyoxal + methionine -> aminoacetaldehyde + methional + CO2 (Strecker, net)", k_go,
                       FROZEN_B18["ea_go_ak_kj_mol"], 2, "L/(mmol*min)",
                       "B22: glycine's B18 constant times the fitted identity ratio; B18's barrier and pH term.",
                       ("identity_ratio_on_b18", "conditional_on_trunk_dicarbonyl_supply", "unit_of_source_inferred")),
        "k_mgo_met": _p("k_mgo_met", "methylglyoxal + methionine -> aminoacetone + methional + CO2 (Strecker, net)", k_mgo,
                        FROZEN_B18["ea_mgo_ak_kj_mol"], 2, "L/(mmol*min)",
                        "B22: glycine's B18 constant times the same identity ratio.",
                        ("identity_ratio_on_b18", "conditional_on_trunk_dicarbonyl_supply", "unit_of_source_inferred")),
        "k_mtal_msh": _p("k_mtal_msh", "methional -> methanethiol + acrolein (retro-Michael)", _pow(log10_k_mtal_msh),
                         float(ea_mtal_msh), 1, "1/min", "B22 fit rows: Pan 2025's three methanethiol rates; barrier fitted within 20-150 kJ/mol.",
                         ("barrier_fitted", "unit_of_source_inferred")),
        "k_msh_dmds": _p("k_msh_dmds", "2 methanethiol -> dimethyl disulfide (APPARENT; no oxidant tracked)", _pow(log10_k_msh_dmds),
                         EA_MSH_DMDS_KJ_MOL, 2, "L/(mmol*min)",
                         "B22 fit rows: Pan 2025's three disulfide rates; barrier declared from Pan's apparent 81 kJ/mol. Pot-specific.",
                         ("apparent_constant_internal_oxidant", "barrier_declared_from_pan2025", "no_transfer_claimed")),
    }


METHIONINE_PARAMETERS: Mapping[str, KineticParameter] = with_fitted_methionine(
    *[(FROZEN_B22 if METHIONINE_SHIPPED else INERT_B22)[k] for k in METHIONINE_COORDINATES])
METHIONINE_KEYS: Tuple[str, ...] = tuple(METHIONINE_PARAMETERS)
#: The two Strecker steps share B18's pH term (trunk_conditions applies it to PYRAZINE_PH_STEPS + these).
METHIONINE_PH_STEPS: Tuple[str, ...] = ("k_go_met", "k_mgo_met")

METHIONINE_WISHLIST: Mapping[str, str] = {
    "k_go_met": "methionine on FED glyoxal or methylglyoxal at two temperatures in water (Zhou 2024's design with methionine), to free the ratio from the trunk's supply",
    "k_mtal_msh": "methional heated alone in buffer at two temperatures with methanethiol measured (Yao 2025's pot with numbers)",
    "k_msh_dmds": "methanethiol loss with the oxidant stated (Chin & Lindsay 1994 at cooking temperature)",
    "unit": "Pan 2025's Table 2 unit, from the authors",
}
