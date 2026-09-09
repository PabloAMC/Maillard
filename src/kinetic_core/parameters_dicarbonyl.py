"""
src/kinetic_core/parameters_dicarbonyl.py -- THE DICARBONYL TRIO (Build Wave B13, 2026-09-07).

Glucosone, glyoxal and diacetyl on the trunk lane, with the constants Kocadagli & Gokmen
2016 (J. Agric. Food Chem. 64:6446, doi 10.1021/acs.jafc.6b01862) fitted to their amine-free
glucose glass at 160 / 180 / 200 C, re-referenced from the paper's T_b = 180 C to the core's
100 C through the same helper Build Wave B7 used for the furanic channel
(`parameters_furanic._kocadagli`). Dossier: `kocadagli2016jafc_extraction.md` sec. 4 (Table 2,
glucose system, k_b in min^-1 x 10^3 with 95 % HPD).

WHY THESE THREE. Glyoxal is A CML precursor by structure and methylglyoxal (already on the
trunk) the CEL precursor: the panel's two AGE rows are refused today for want of the species.
CORRECTION (2026-09-09, Nguyen 2016, Berk 2021, Hamzalioglu 2026 read): in every aqueous
multiresponse fit on disk the glyoxal -> CML route fits to ZERO, or a thousand to ten thousand
times below the Amadori -> CML route (casein + glucose or lactose at 120-140 C, sesame at
180-220 C, whole milk at 110-140 C); CML comes from the Amadori compound's oxidative cleavage.
A CML row that charges glyoxal as the precursor will under-predict; the Amadori route is the
one to write when the CML row is taken up (tasks/data_restructure_plan.md, glycation log). Diacetyl
is the buttery odorant and the 3-mercapto-2-butanone precursor Yiltirak 2026 quantifies.
Glucosone is the oxidative entry that makes glyoxal at all.

WHAT IS NOT CLAIMED. No fit row of this repository measures any of the three; the constants
are one laboratory's, one amine-free matrix's, 160-200 C, extrapolated downward. Two sinks
carry the authors' own decisions (glyoxal barrier fixed to zero; diacetyl sink rate zero) and
are flagged. The wishlist lists the measurement that would replace each. The sulfur network
does NOT carry these steps (B9's topology is frozen); a sulfur wave that wants diacetyl for
3-mercapto-2-butanone adopts them then.
"""
from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .parameters import KineticParameter
from .parameters_furanic import _kocadagli

#: Kocadagli 2016 JAFC Table 2, GLUCOSE system: (k_b x 1000 at 180 C, HPD), (Ea, HPD).
DICARBONYL_PARAMETERS: Mapping[str, KineticParameter] = {
    "k_glc_g": _kocadagli(
        "k_glc_g", "glucose -> glucosone", 0.069, 125.9, 4.9, 9,
        flags=("b13_dicarbonyl", "oxidative_entry_in_air"),
        note="k_b 0.069 +/- 0.005 x 1e-3 /min at 180 C; Ea 125.9 +/- 4.9. The smallest "
             "entry on the trunk by two decades; at 100-145 C it is 1e-6 to 1e-5 /min.",
    ),
    "k_g_go": _kocadagli(
        "k_g_go", "glucosone -> glyoxal + C4", 737.0, 93.8, 6.4, 10,
        flags=("b13_dicarbonyl",),
        note="k_b 737 +/- 58.9 x 1e-3 /min at 180 C; Ea 93.8 +/- 6.4. Glucosone is transient: "
             "this step is four decades faster than its formation.",
    ),
    "k_odg_da": _kocadagli(
        "k_odg_da", "1-deoxyglucosone -> diacetyl + C2", 12.2, 150.8, 8.8, 12,
        flags=("b13_dicarbonyl",),
        note="k_b 12.2 +/- 1.12 x 1e-3 /min at 180 C; Ea 150.8 +/- 8.8 -- the steepest barrier on "
             "the trunk, so diacetyl is a high-temperature product. 2026-09-09: 466x below a "
             "second laboratory at 160 C -- SECOND_LABORATORY_2016.",
    ),
    "k_go_sink": _kocadagli(
        "k_go_sink", "glyoxal -> unassigned (P3)", 32.6, 0.0, None, 15,
        flags=("b13_dicarbonyl", "ea_fixed_to_zero_by_authors"),
        note="k_b 32.6 +/- 8.83 x 1e-3 /min at 180 C; the authors FIXED the barrier to zero "
             "during estimation, so the sink runs at its 180 C rate at every temperature. "
             "Declared, flagged; the wishlist asks for a glyoxal loss rate at two temperatures. "
             "2026-09-09: a second laboratory now supplies three temperatures and the ZERO "
             "BARRIER IS REFUTED while the rate agrees to 1.87x -- SECOND_LABORATORY_2016.",
    ),
    "k_da_sink": _kocadagli(
        "k_da_sink", "diacetyl -> unassigned (P5)", 0.0, 0.0, None, 17,
        flags=("b13_dicarbonyl", "rate_zero_in_source"),
        note="Kocadagli step 17: 0 +/- 0 (blank Ea). Diacetyl accumulates in the source's "
             "glass; carried at zero as a PREDICTION the data may reject. 2026-09-09: A "
             "SECOND LABORATORY HAS REJECTED IT -- 130e-3 /min in roasted hazelnut. The "
             "value stays until a wave installs one; see SECOND_LABORATORY_2016.",
    ),
}

DICARBONYL_KEYS: Tuple[str, ...] = tuple(DICARBONYL_PARAMETERS)

# ===========================================================================
# THE SECOND LABORATORY (2026-09-09, from the reading audit)
# ===========================================================================
# Until this date every constant on the trunk came from ONE laboratory, one amine-free
# glass, 160-200 C. Goncuoglu Tas & Gokmen 2016 (whole Tombul hazelnuts, 5 g, 150 / 160 /
# 170 C, 15-120 min, multiresponse fit over 26 steps) is a second laboratory, a second
# matrix and a real food. Its Table 1 is NOT installed here and NOTHING below changes a
# shipped value -- a refit is a wave, and this is the record that says which way the wave
# would push. Rates are as printed, per minute; where the shipped constant lives at another
# temperature it is transported by its own barrier and that is marked.
#
# WHAT AGREES. Three constants inside a factor of two, across two laboratories, two
# matrices and a 20 C gap: `k_tdg_ddg` 1.5x, `k_ddg_hmf` 1.13x, and `k_go_sink`'s RATE
# 1.87x. That is the first cross-laboratory agreement the trunk has ever had and it is the
# more important half of this record.
#
# WHAT DOES NOT. Two shipped decisions are refuted outright and two constants disagree by
# orders of magnitude. Both refuted decisions were flagged as decisions when they shipped,
# which is why they are checkable now.
SECOND_LABORATORY_2016: Mapping[str, Mapping[str, object]] = {
    "k_da_sink": {
        "shipped": "0 /min, Ea blank: the source measured no diacetyl loss in its glass",
        "measured_elsewhere": "54 / 130 / 106 x 1e-3 /min at 150 / 160 / 170 C",
        "verdict": "THE PREDICTION IS REJECTED. The constant ships at zero as 'a prediction the "
                   "data may reject' and a second laboratory has now rejected it: diacetyl is "
                   "consumed in a real matrix. The rate is also NON-MONOTONE in temperature over "
                   "three points, so this measurement supplies a size and not a barrier.",
        "anchor": "goncuoglu2016_extraction.md sec. 3 Table 1 step k24",
    },
    "k_go_sink": {
        "shipped": "32.6e-3 /min with the BARRIER FIXED TO ZERO by its authors, so it runs at its "
                   "180 C rate at every temperature",
        "measured_elsewhere": "18 / 61 / 290 x 1e-3 /min at 150 / 160 / 170 C",
        "verdict": "THE RATE AGREES (1.87x at 160 C) AND THE ZERO BARRIER IS REFUTED: a sixteenfold "
                   "rise over 20 C is not a zero barrier. The window is too narrow to put a credible "
                   "barrier in its place -- a three-point refit gives about 216 kJ/mol, which is not "
                   "believable -- so the wishlist entry stands unchanged and is now evidenced rather "
                   "than merely prudent.",
        "anchor": "goncuoglu2016_extraction.md sec. 3 Table 1 step k25",
    },
    "k_odg_da": {
        "shipped": "12.2e-3 /min at 180 C with Ea 150.8, which transports to 1.92e-3 /min at 160 C",
        "measured_elsewhere": "371 / 895 / 1073 x 1e-3 /min at 150 / 160 / 170 C",
        "verdict": "466x APART AT 160 C, and the two credible intervals do not come within two "
                   "decades of each other. Diacetyl is made far faster in a roasting nut than in an "
                   "amine-free glass. Nothing here says which is right for a plant-protein cook; it "
                   "says the constant is matrix-dependent and the model carries one matrix.",
        "anchor": "goncuoglu2016_extraction.md sec. 3 Table 1 step k17",
    },
    "k_hmf_self": {
        "shipped": "8.97e-7 /min, Ea zero by declaration (0.9 % lost in 7 days at 5 C), in "
                   "parameters_furanic.py",
        "measured_elsewhere": "12 / 21 / 103 x 1e-3 /min at 150 / 160 / 170 C",
        "verdict": "23 000x APART. The furanic channel already prints 'EXPECT HMF TO BE "
                   "OVER-PREDICTED' and names the empty 50-150 C window as the reason; this puts a "
                   "measured size on that warning. Half-life at 160 C: about 33 minutes measured "
                   "against about 1.5 years shipped.",
        "anchor": "goncuoglu2016_extraction.md sec. 3 Table 1 step k26",
    },
}
#: What this laboratory does NOT settle, said plainly so nobody reads the table as a refit.
SECOND_LABORATORY_2016_LIMITS = (
    "A roasting hazelnut is a dry, lipid-rich, whole-tissue matrix and the shipped constants come "
    "from an aqueous amine-free glass; the authors say themselves that Arrhenius fails across their "
    "own three temperatures. The dicarbonyl ORDER also inverts between the two: glyoxal about equal "
    "to methylglyoxal above 3-deoxyglucosone here, the reverse in aqueous Leitzen 2021. So these "
    "numbers size a disagreement; they do not replace anything, and a wave that installs any of them "
    "has to say which matrix it claims to be modelling."
)

#: What would replace each declared decision (read by the wishlist through the flags).
DICARBONYL_WISHLIST: Mapping[str, str] = {
    "k_go_sink": "glyoxal loss from a glucose/glycine pot at two temperatures (a barrier for the sink)",
    "k_da_sink": "diacetyl loss from a glucose/glycine pot heated alone (the source measured none)",
    "k_glc_g": "glucosone in a glucose/glycine solution at 100-145 C (the entry is extrapolated from a 160-200 C glass)",
}

# ===========================================================================
# BUILD WAVE B21 (2026-09-09): THE AQUEOUS GLUCOSONE ROUTE TO GLYOXAL
# ===========================================================================
# In water the glucosone comes from the Amadori compound, not from the sugar: Hamzalioglu 2026
# (whole milk, lactulosyl-lysine, 110-140 C, multiresponse fit) is the only aqueous entry on disk,
# with first-order constants (basis-free) and a measured barrier. One new step, r_ama_g
# (AMA -> G + Gly), and an AQUEOUS value for the glucosone -> glyoxal constant k_g_go that replaces
# the glass value in the operative set (the glass entry above stays as the record and as the
# "before"). Pre-registered in results/validation/kinetic_core_b21_prereg.md; the fitted values
# below are FROZEN LITERALS asserted equal to the B21 fit report by tests/unit/test_kinetic_core_b21.py.
EA_AMA_G_KJ_MOL = 75.9            # Hamzalioglu 2026 Table 2, lactulosyl-lysine -> glucosone, +/- 21.1 (measured)
EA_G_GO_AQUEOUS_KJ_MOL = 4.2      # Hamzalioglu 2026 Table 2, glucosone -> glyoxal, +/- 15.7 (consistent with zero; the glass says 93.8)
#: The declared transfer from lactulosyl-lysine in milk to fructosyl-glycine in water, in decades.
AQUEOUS_TRANSFER_BAND_DECADES = 0.5
_HAMZALIOGLU = ("Hamzalioglu, Kocadagli & Gokmen 2026, J. Agric. Food Chem. (whole milk, lactose + casein-bound lysine, 110 / 120 / 130 / "
                "140 C, 0.5-5 min, multiresponse fit): Table 1 steps 4 (LacLys -> glucosone) and 10 (glucosone -> glyoxal), Table 2 barriers; "
                "hamzalioglu2026_extraction.md sec. 4")
#: The B21 fit report's optimum (2026-09-09: cost 3.5 on six rows, both coordinates identified); asserted against the report by the unit test.
FROZEN_B21: Mapping[str, float] = {
    "log10_k_ama_g_100C": -2.0339128364819725,
    "log10_k_g_go_aqueous_100C": -0.5104286909031046,
}
AQUEOUS_GLYOXAL_COORDINATES: Tuple[str, ...] = tuple(FROZEN_B21)
AQUEOUS_GLYOXAL_CAVEAT = (
    "GLYOXAL SUPPLY (B21): the glucosone that makes glyoxal comes from the Amadori compound at a rate fitted on "
    "one laboratory's milk constants (lactulosyl-lysine, 110-140 C) and carried to fructosyl-glycine in water as "
    "a declared transfer (+/- 0.5 dex on every glyoxal, glucosone and pyrazine answer); the glucosone -> glyoxal "
    "barrier is that laboratory's, consistent with zero; the dry-glass glyoxal sink (B13) is unchanged. The route "
    "drains the Amadori compound at the size of the trunk's own Amadori -> 3-deoxyglucosone step, and Martins "
    "2005's Amadori series wants a fifth less of that drain (median error 0.035 -> 0.093 dex): a joint refit is "
    "the next pre-registration on this route."
)


def with_aqueous_glyoxal(log10_k_ama_g: float, log10_k_g_go: float) -> Dict[str, KineticParameter]:
    """The aqueous glyoxal-supply block at arbitrary values (the fit generator's hook and the report reader's)."""
    common = dict(evidence_class="derived_from_fit_data", source_anchor=_HAMZALIOGLU,
                  dossier_anchor="hamzalioglu2026_extraction.md sec. 4; results/validation/kinetic_core_b21_prereg.md",
                  conditions="water (milk serum), 110-140 C; first-order constants, no water basis needed; lactulosyl-lysine -> "
                             "fructosyl-glycine declared (+/- 0.5 dex)",
                  ph_of_measurement=6.7, temperature_range_c=(110.0, 140.0), rate_transfer="licensed_at_measurement_ph_only", unit="1/min", order=1)
    return {
        "k_ama_g": KineticParameter(
            key="k_ama_g", transformation="Amadori (DFG) -> glucosone + Gly (oxidative cleavage; returns the amine)",
            k_ref=10.0 ** float(log10_k_ama_g), ea_kj_mol=EA_AMA_G_KJ_MOL,
            flags=("b21_aqueous_glyoxal", "fitted_wave_b21", "barrier_measured_hamzalioglu2026", "transfer_declared_laclys_to_fructosylglycine"),
            note="B21 fit rows: Hamzalioglu 2026's four LacLys -> glucosone constants (the 130 C one wide). The only aqueous glucosone entry on disk.",
            **common),
        "k_g_go": KineticParameter(
            key="k_g_go", transformation="glucosone -> glyoxal + C4 residue (AQUEOUS value; the glass value is the record)",
            k_ref=10.0 ** float(log10_k_g_go), ea_kj_mol=EA_G_GO_AQUEOUS_KJ_MOL,
            flags=("b21_aqueous_glyoxal", "fitted_wave_b21", "barrier_consistent_with_zero", "replaces_b13_glass_value_in_the_operative_set"),
            note="B21 fit rows: Hamzalioglu 2026's two determinate glucosone -> glyoxal constants (120, 130 C); the 110 and 140 C ones "
                 "are indeterminate and reported. The glass barrier (93.8) is the other reading and is reported beside it.",
            **common),
    }


AQUEOUS_GLYOXAL_PARAMETERS: Mapping[str, KineticParameter] = with_aqueous_glyoxal(
    FROZEN_B21["log10_k_ama_g_100C"], FROZEN_B21["log10_k_g_go_aqueous_100C"])
AQUEOUS_GLYOXAL_KEYS: Tuple[str, ...] = tuple(AQUEOUS_GLYOXAL_PARAMETERS)
#: The glass value of k_g_go, kept for the ship rule's "before" and for the record.
GLASS_K_G_GO: KineticParameter = DICARBONYL_PARAMETERS["k_g_go"]

__all__ = ["DICARBONYL_KEYS", "DICARBONYL_PARAMETERS", "DICARBONYL_WISHLIST", "AQUEOUS_GLYOXAL_PARAMETERS", "AQUEOUS_GLYOXAL_KEYS",
           "AQUEOUS_GLYOXAL_COORDINATES", "AQUEOUS_GLYOXAL_CAVEAT", "FROZEN_B21", "with_aqueous_glyoxal", "GLASS_K_G_GO"]
