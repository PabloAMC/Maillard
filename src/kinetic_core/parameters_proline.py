"""
src/kinetic_core/parameters_proline.py -- 2-ACETYL-1-PYRROLINE FROM PROLINE (Build Wave B24, 2026-09-09).

Two steps on the trunk lane (network.PROLINE_REACTIONS): the Strecker of proline on methylglyoxal,
which gives 1-pyrroline (the ring nitrogen stays; hydroxyacetone and CO2 leave to the fragment pool),
and the acylation of 1-pyrroline by methylglyoxal to 2-acetyl-1-pyrroline (Hofmann & Schieberle 1998b:
hydrated methylglyoxal adds at C-2, the aldehyde carbon leaves as CO2, air oxidises the pyrrolidine).

WHAT IS FITTED (wave B24, pre-registered in results/validation/kinetic_core_b24_prereg.md): the two
log10 constants at 100 C, on Hofmann & Schieberle 1998b's fed-pyrroline yields (Table 7, experiments
1 and 2) and proline + methylglyoxal yields (Table 9, three ratios), all at 100 C, pH 7, 0.5 M
phosphate, 30 min. The values below are FROZEN LITERALS asserted equal to the B24 fit report by
tests/unit/test_kinetic_core_b24.py.

WHAT IS DECLARED. No barrier is measured for either step: the acylation carries Chan & Reineccius
1994's whole-cascade apparent barrier for 2-acetyl-1-pyrroline (60.2 kJ/mol, glucose + amino acids
in water, 75-115 C; flagged as apparent), the proline Strecker carries B18's methylglyoxal Strecker
barrier (114.9 kJ/mol, measured on glycine). No loss of 1-pyrroline or of the product is written
(Hofmann's excess-pyrroline experiment shows one; reported). Proline is charged, declared, as glycine
for the Amadori chemistry. The product's odour threshold is not on disk, so no odour-activity value
is computed for it.
"""
from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .parameters import AW_OF_MEASUREMENT, KineticParameter
from .parameters_pyrazine import FROZEN_B18

PROLINE_FIT_PH = 7.0
EA_PYRL_AP_KJ_MOL = 60.2              # Chan & Reineccius 1994 (RSC chapter), 2-acetyl-1-pyrroline, whole cascade, apparent
EA_MGO_PRO_KJ_MOL = float(FROZEN_B18["ea_mgo_ak_kj_mol"])   # B18's methylglyoxal Strecker barrier, measured on glycine
MOLAR_MASS_G_PER_MOL: Mapping[str, float] = {"PRO": 115.13, "PYRL": 69.11, "AP": 111.14}

#: The B24 fit report's optimum (asserted against the report by the unit test). Prior centres: Hofmann's
#: bilinear lower bounds (acylation ~9e-4 L/(mmol min) at 100 C; the proline Strecker from Table 9's
#: yield per proline over 30 min at 4 mmol/L methylglyoxal, through the acylation).
#: THE B24 FIT REPORT'S OPTIMUM, kept as the record and asserted against the report by the unit test. NOT
#: installed: the wave did not ship (kinetic_core_b24_prereg.md sec. 6). The fed-pyrroline acylation rows fit
#: (within 0.3 dex) but the whole chain from proline does not: the model's yield rises linearly with the
#: methylglyoxal charge where Hofmann's saturates, because the 1-pyrroline sinks (the tetrahydropyridine
#: branch through hydroxyacetone, the pyrroline's own loss) are not written.
FROZEN_B24: Mapping[str, float] = {
    "log10_k_pyrl_ap_100C": -2.561776899512073,
    "log10_k_mgo_pro_100C": -6.620784383439292,
}
#: The wave's verdict. While False the two constants are carried at ZERO and the targets are refused by name.
PROLINE_SHIPPED = False
INERT_B24: Mapping[str, float] = {"log10_k_pyrl_ap_100C": -300.0, "log10_k_mgo_pro_100C": -300.0}
PROLINE_COORDINATES: Tuple[str, ...] = tuple(FROZEN_B24)

# ===========================================================================
# WAVE B24b (2026-09-10): THE BRANCH THAT REFUSED B24
# ===========================================================================
# B24's own outcome names this wave and what it needs. Two more species (hydroxyacetone and the
# tetrahydropyridine) and three coordinates: the branch constant, the pH gate on it, and the
# 1-pyrroline loss B24 had none of.
#
# THE ONE THING THAT MAKES THIS FITTABLE, and it is a statement in a source rather than an
# inference: Schieberle & Hofmann 2005 say in words that hydroxyacetone gives ONLY the
# tetrahydropyridine and methylglyoxal gives ONLY 2-acetyl-1-pyrroline. So the two products are not
# competing rates on one substrate. They are competing claims on the METHYLGLYOXAL -- the amino acid
# turns it into hydroxyacetone and commits it to one branch, the pyrroline acylates it and commits
# it to the other. That is why an excess of proline drives the tetrahydropyridine and an excess of
# methylglyoxal drives the pyrroline product, and why the fit target is a RATIO.
EA_HA_ATHP_KJ_MOL = EA_PYRL_AP_KJ_MOL   # declared equal to the acylation's: no barrier is measured
                                        # for this step anywhere, and both are condensations of
                                        # 1-pyrroline with a small carbonyl in the same pot.
#: Schieberle & Hofmann 2005 Table 2, the tetrahydropyridine against pH, as RATIOS to the pH-7 rung.
#: The pH-3 rung is CENSORED ("<0.1 ug") and enters as a one-sided bound, never as a point.
ATHP_PH_LADDER_UG: Mapping[str, float] = {"5.0": 0.9, "7.0": 10.8, "9.0": 38.4}
ATHP_PH_LADDER_CENSORED: Mapping[str, float] = {"3.0": 0.1}
ATHP_PH_SOURCE = ("Schieberle & Hofmann 2005 (ACS Symp. Ser.) Table 2, re-read from 200-dpi rasters "
                  "because the OCR layer dropped the pH 9.0 row; identical charge, buffer, time and "
                  "temperature to Hofmann & Schieberle 1998b Table 4, which it confirms to the digit. "
                  "schieberle2005_extraction.md")
#: Hofmann & Schieberle 1998b Table 9: the AP : ATHP molar ratio across a hundredfold methylglyoxal
#: ladder against 400 mmol/L proline. THE SHAPE THIS WAVE IS FITTED ON.
AP_ATHP_SWITCH: Mapping[float, float] = {4.0: 0.16, 40.0: 0.51, 400.0: 12.8}
AP_ATHP_SWITCH_SOURCE = ("Hofmann & Schieberle 1998b Table 9, proline 400 mmol/L + methylglyoxal "
                         "4 / 40 / 400 mmol/L, 0.5 M phosphate pH 7, 100 C, 30 min; a ratio inside "
                         "one analysis, so the response factor and the extraction cancel")
#: Experiment 3: 1-pyrroline in fivefold excess over methylglyoxal. What sizes the pyrroline loss.
PYRL_EXCESS_ROW: Mapping[str, float] = {"pyrl_mmol_l": 10.0, "mgo_mmol_l": 2.0, "ap_molpct_of_mgo": 0.33}
_HOFMANN = ("Hofmann & Schieberle 1998b, J. Agric. Food Chem. 46:2270 (doi 10.1021/jf970990g): Table 7 (1-pyrroline 2 mmol/L + "
            "methylglyoxal 10 or 2 mmol/L, 0.5 M phosphate pH 7, 100 C, 30 min: 28.7 and 5.3 mol %), Table 9 (proline 400 mmol/L + "
            "methylglyoxal 4 / 40 / 400 mmol/L: 0.0058 / 0.0125 / 0.0179 mol % of proline); hofmann1998b_extraction.md sec. 4")
PROLINE_CAVEAT = (
    "2-ACETYL-1-PYRROLINE (B24): two constants fitted on one laboratory's 30-minute yields at 100 C and pH 7 (fed "
    "1-pyrroline, and proline + methylglyoxal at three ratios); NO barrier is measured for either step (the acylation "
    "carries a whole-cascade apparent 60 kJ/mol, the proline Strecker glycine's 115), no loss of 1-pyrroline or of the "
    "product is written, proline is charged as glycine for the Amadori chemistry, and the product's odour threshold "
    "is not on disk. A temperature other than 100 C is an extrapolation on declared barriers."
)


PROLINE_NOT_SHIPPED_REASON = (
    "2-acetyl-1-pyrroline (wave B24) did not ship: the acylation of fed 1-pyrroline by methylglyoxal fits Hofmann & "
    "Schieberle 1998b within 0.3 dex, but the whole chain from proline does not (the yield rises linearly with the "
    "methylglyoxal charge in the model and saturates in the source, 1.3 dex low at 4 mmol/L and 1.2 dex high at 400), "
    "because the 1-pyrroline sinks are not written: the tetrahydropyridine branch through hydroxyacetone and the "
    "pyrroline's own loss (the source's excess-pyrroline pot is 2.1 dex below the model). The steps stay in the "
    "network at zero and the target is refused rather than answered with a chain the source refutes "
    "(kinetic_core_b24_prereg.md sec. 6)."
)


def _p(key, transformation, log10_k, ea, note, flags):
    return KineticParameter(
        key=key, transformation=transformation, k_ref=(0.0 if float(log10_k) < -100.0 else 10.0 ** float(log10_k)), ea_kj_mol=float(ea), unit="L/(mmol*min)", order=2,
        evidence_class="derived_from_fit_data", source_anchor=_HOFMANN,
        dossier_anchor="hofmann1998b_extraction.md sec. 4; chan1994b_extraction.md sec. 4; results/validation/kinetic_core_b24_prereg.md",
        conditions="water, 0.5 M phosphate pH 7, 100 C, 30 min; fed 1-pyrroline and fed methylglyoxal",
        ph_of_measurement=PROLINE_FIT_PH, temperature_range_c=(100.0, 100.0), rate_transfer="licensed_at_measurement_ph_only",
        aw_of_measurement=AW_OF_MEASUREMENT, flags=("b24_proline", "fitted_wave_b24", "barrier_declared_not_fitted") + tuple(flags), note=note,
    )


def with_fitted_proline_b24b(log10_k_ha_athp: float, log10_k_pyrl_loss: float) -> Dict[str, KineticParameter]:
    """B24b's two new constants. Zero until a B24b report supplies them."""
    return {
        "k_ha_athp": _p("k_ha_athp", "1-pyrroline + hydroxyacetone -> 2-acetyltetrahydropyridine", log10_k_ha_athp,
                        EA_HA_ATHP_KJ_MOL,
                        "B24b fit rows: Hofmann & Schieberle 1998b Table 9's three AP : ATHP ratios, and "
                        "Schieberle & Hofmann 2005 Table 2's pH ladder as ratios to pH 7. Barrier DECLARED "
                        "equal to the acylation's; none is measured for this step anywhere.",
                        ("barrier_declared_equal_to_acylation", "exclusive_branch_stated_by_source")),
        "k_pyrl_loss": _p("k_pyrl_loss", "1-pyrroline -> melanoidin pools (its own loss, first order)", log10_k_pyrl_loss,
                          EA_HA_ATHP_KJ_MOL,
                          "B24b fit row: Hofmann & Schieberle 1998b Table 7 experiment 3, 1-pyrroline in "
                          "fivefold excess giving 0.33 mol % of the methylglyoxal where B24 made 41. The "
                          "source measures a DISAPPEARANCE and names no product, so the sink is accounting "
                          "and not mechanism.",
                          ("first_order_is_the_crudest_form", "no_product_named_by_source")),
    }


def with_fitted_proline(log10_k_pyrl_ap: float, log10_k_mgo_pro: float) -> Dict[str, KineticParameter]:
    """The proline block at arbitrary values (the fit generator's hook and the report reader's)."""
    return {
        "k_pyrl_ap": _p("k_pyrl_ap", "1-pyrroline + methylglyoxal -> 2-acetyl-1-pyrroline + CO2 (acylation; net)", log10_k_pyrl_ap,
                        EA_PYRL_AP_KJ_MOL, "B24 fit rows: Hofmann 1998b Table 7 experiments 1 and 2. Barrier: Chan & Reineccius 1994's whole-cascade 60.2 kJ/mol (apparent).",
                        ("barrier_from_chan1994b_whole_cascade",)),
        "k_mgo_pro": _p("k_mgo_pro", "methylglyoxal + proline -> 1-pyrroline + hydroxyacetone + CO2 (Strecker of proline; net)", log10_k_mgo_pro,
                        EA_MGO_PRO_KJ_MOL, "B24 fit rows: Hofmann 1998b Table 9 (three methylglyoxal : proline ratios). Barrier: B18's methylglyoxal Strecker (glycine).",
                        ("barrier_from_b18_mgo_strecker",)),
    }


#: B24b's optimum, once its report exists. Inert until then, exactly as B24's two are.
FROZEN_B24B: Mapping[str, float] = {}
B24B_SHIPPED = False
INERT_B24B: Mapping[str, float] = {"log10_k_ha_athp_100C": -300.0, "log10_k_pyrl_loss_100C": -300.0}
PROLINE_B24B_COORDINATES: Tuple[str, ...] = tuple(INERT_B24B)

PROLINE_PARAMETERS: Mapping[str, KineticParameter] = {
    **with_fitted_proline(*[(FROZEN_B24 if PROLINE_SHIPPED else INERT_B24)[k] for k in PROLINE_COORDINATES]),
    **with_fitted_proline_b24b(
        *[(FROZEN_B24B if B24B_SHIPPED else INERT_B24B)[k] for k in PROLINE_B24B_COORDINATES]),
}
PROLINE_KEYS: Tuple[str, ...] = tuple(PROLINE_PARAMETERS)
#: The proline Strecker step takes B18's pH term (an amine-dependent step).
#: B24b adds k_ha_athp: it condenses 1-pyrroline (an amine) with a carbonyl, the same shape B18's
#: term was fitted for. The transfer is DECLARED and CHECKED against a ladder it was not fitted on
#: (Schieberle & Hofmann 2005 Table 2); the check is reported on the B24b ship rule.
PROLINE_PH_STEPS: Tuple[str, ...] = ("k_mgo_pro", "k_ha_athp")

PROLINE_WISHLIST: Mapping[str, str] = {
    "k_pyrl_ap": "1-pyrroline + methylglyoxal at a second temperature (Hofmann's pot at 80 or 120 C) for a measured barrier",
    "k_mgo_pro": "proline + a fed dicarbonyl with 1-pyrroline measured, at two temperatures",
    "loss": "1-pyrroline's own loss (Hofmann's excess-pyrroline suppression) and 2-acetyl-1-pyrroline's stability in water",
    "threshold": "the product's odour threshold in water from a paper on disk (Buttery 1983)",
}
