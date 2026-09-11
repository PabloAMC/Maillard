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

import math

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

# ---------------------------------------------------------------------------
# ENV-B13 (2026-09-10): the disputed sinks, at a drawn value
# ---------------------------------------------------------------------------
# The Monte-Carlo envelope had NO prior row for any of these constants, so every published interval
# asserted them with certainty -- including two whose own authors flagged them as decisions and
# which a second laboratory has since refuted. This is the hook the envelope moves them through.
# It changes nothing on the default path: `DISPUTED_SINK_KEYS` is the whole surface, and with no
# override the module's own literals are used exactly as before.
#
# THE BANDS ARE NOT A REFIT. Each spans the two laboratories at a common temperature, and the
# centre stays where it shipped, because which laboratory is right for a plant-protein cook is not
# settled by either of them -- an aqueous amine-free glass against a dry, lipid-rich whole nut,
# disagreeing on the dicarbonyl ORDER as well as the rates.
_R = 8.314462618          # J/(mol K)
_T_REF_K = 373.15         # the module's k_ref reference, 100 C
#: Where a DRAWN barrier pivots. Only k_go_sink needs one: its two laboratories agree on the RATE
#: at 160 C and disagree only about the barrier, so that is the temperature the draw must preserve.
_BARRIER_DRAW_ANCHOR_T_K: Mapping[str, float] = {"k_go_sink": 433.15}

DISPUTED_SINK_KEYS: Tuple[str, ...] = ("k_da_sink", "k_go_sink", "k_odg_da", "k_hmf_self")

#: (log10 k at 100 C lo, hi) and (Ea lo, hi) for each. `None` means "not in dispute, do not sample".
#: Every endpoint below is traceable to a printed number; see SECOND_LABORATORY_2016 for the sources.
DISPUTED_SINK_BANDS: Mapping[str, Mapping[str, object]] = {
    "k_da_sink": {
        "log10_k_100C": (-12.0, -1.94),
        "ea_kj_mol": None,
        "basis": "shipped at EXACTLY ZERO with a blank barrier, which -12 stands in for; the upper "
                 "end is Goncuoglu Tas 2016's 130e-3 /min at 160 C carried to 100 C at this "
                 "module's own steepest trunk barrier. Flat between them says 'somewhere between "
                 "nothing and what the second laboratory measured', which is the state of knowledge.",
    },
    "k_go_sink": {
        "log10_k_100C": None,
        "ea_kj_mol": (0.0, 150.8),
        "basis": "the RATE is left alone: the two laboratories agree to 1.87x at 160 C, which is "
                 "the best cross-laboratory agreement on this trunk. The BARRIER is the disputed "
                 "part -- its authors FIXED it to zero and a sixteenfold rise over 20 C refutes "
                 "that. The upper end is the steepest barrier this module itself carries "
                 "(k_odg_da, 150.8). The 20 C window is too narrow to fit a credible barrier and "
                 "this band does not pretend to be one; it is the honest width of not knowing.",
    },
    "k_odg_da": {
        "log10_k_100C": (-4.63, -1.96),
        "ea_kj_mol": None,
        "basis": "the 466x disagreement at 160 C, expressed at 100 C on the shipped barrier. The "
                 "two credible intervals do not come within two decades of each other.",
    },
    "k_hmf_self": {
        "log10_k_100C": (-6.05, -0.96),
        "ea_kj_mol": None,
        "basis": "from the shipped 8.97e-7 /min (one temperature, barrier zero by declaration) up "
                 "to Gokmen 2012's 0.111 /min at 180 C, the widest of the three readings. The "
                 "barrier stays at zero: one temperature each, so no Arrhenius is licensed.",
    },
}
#: ENV-B34 (2026-09-11): the 3-deoxyglucosone limb and the amine-free sugar entries, banded on the
#: SOURCE'S OWN PRINTED 95 % HPD (Kocadagli & Gokmen 2016 Table 2, glucose system, reparameterised
#: Arrhenius; kocadagli2016jafc_extraction.md sec. 4). No centre moves. The rate band is the relative
#: HPD on k_b applied to the shipped 100 C value; the barrier band is the printed Ea +/- HPD. Drawn
#: independently, the same simplification ENV-B13 made. Wave B34 found k_tdg_ddg 32x too slow in
#: water on the one aqueous pot that measures its product; until then it sat in AGREEING_SINK_KEYS
#: because two laboratories agree on it to 1.5x -- both of them measuring it DRY, and this paper's
#: own NaCl column printing a barrier of 117.7 +/- 11.1 for the same step against 36.9 +/- 6.3 in
#: the glucose column. Agreement is not accuracy when both laboratories share a matrix.
HPD_SINK_BANDS: Mapping[str, Mapping[str, object]] = {
    "k_glc_tdg": {"k_rel_hpd": 2.44 / 4.19, "ea_hpd_kj_mol": 52.7,
                  "basis": "Table 2 step 3, k_b 4.19 +/- 2.44 (x1e-3), Ea 107.2 +/- 52.7"},
    "k_tdg_ddg": {"k_rel_hpd": 3.39 / 30.5, "ea_hpd_kj_mol": 6.3,
                  "basis": "Table 2 step 4, k_b 30.5 +/- 3.39 (x1e-3), Ea 36.9 +/- 6.3; the NaCl column prints Ea 117.7 +/- 11.1 for the same step"},
    "k_fru_int": {"k_rel_hpd": 22.8 / 330.0, "ea_hpd_kj_mol": 6.6,
                  "basis": "Table 2 step 6, k_b 330 +/- 22.8 (x1e-3), Ea 100.4 +/- 6.6"},
    "k_fru_odg": {"k_rel_hpd": 0.40 / 2.11, "ea_hpd_kj_mol": 21.8,
                  "basis": "Table 2 step 8, k_b 2.11 +/- 0.40 (x1e-3), Ea 99.3 +/- 21.8"},
}
HPD_SINK_KEYS: Tuple[str, ...] = tuple(HPD_SINK_BANDS)

#: The three that are NOT sampled, listed so a reader finds a row rather than a silence.
AGREEING_SINK_KEYS: Tuple[str, ...] = ("k_glc_g", "k_g_go", "k_ddg_hmf")
AGREEING_SINK_REASON = (
    "NOT SAMPLED, and not by oversight. k_ddg_hmf (1.13x) and k_go_sink's RATE (1.87x) are "
    "cross-laboratory agreements inside a factor of two, and k_ddg_hmf is a timescale bracket with a "
    "declared zero barrier and no HPD to draw from. k_glc_g and k_g_go have one determination each "
    "and no second laboratory to disagree with them. A band invented for a constant nobody disputes "
    "would be a fabricated interval. ENV-B34 (2026-09-11) REMOVED k_tdg_ddg from this list: its 1.5x "
    "cross-laboratory agreement was two DRY matrices agreeing with each other, wave B34 measured it 32x "
    "too slow in water, and its own paper prints a threefold barrier disagreement between its glucose "
    "and NaCl columns. It carries the printed HPD now (HPD_SINK_BANDS)."
)


def with_disputed_sinks(overrides: Mapping[str, Mapping[str, float]]) -> Dict[str, KineticParameter]:
    """
    Rebuild the disputed sinks at drawn values. ENV-B13's only hook into this module.

    ``overrides`` maps a key in ``DISPUTED_SINK_KEYS`` to ``{"log10_k_100C": x}`` and/or
    ``{"ea_kj_mol": y}``. Anything absent keeps the shipped literal.
    """
    from dataclasses import replace

    out: Dict[str, KineticParameter] = {}
    for key, block in overrides.items():
        if key not in DISPUTED_SINK_KEYS and key not in HPD_SINK_KEYS:
            raise KeyError(f"{key} is not a banded sink; ENV-B13 moves {DISPUTED_SINK_KEYS} and ENV-B34 {HPD_SINK_KEYS}")
        base = DICARBONYL_PARAMETERS.get(key)
        if base is None:
            from .parameters_furanic import FURANIC_PARAMETERS

            base = FURANIC_PARAMETERS[key]
        fields: Dict[str, object] = {}
        if block.get("log10_k_100C") is not None:
            fields["k_ref"] = 10.0 ** float(block["log10_k_100C"])
        if block.get("ea_kj_mol") is not None:
            ea = float(block["ea_kj_mol"])
            fields["ea_kj_mol"] = ea
            if key in _BARRIER_DRAW_ANCHOR_T_K and block.get("log10_k_100C") is None:
                # A DRAWN BARRIER MUST PIVOT ABOUT WHERE THE EVIDENCE IS, NOT ABOUT 100 C.
                # k_go_sink's two laboratories agree to 1.87x at 160 C and disagree only about the
                # barrier. Holding k_ref at 100 C while the barrier swings would throw the rate at
                # 160 C across decades and destroy the one agreement this trunk has. So the drawn
                # barrier is applied with the rate held at the anchor temperature, and k_ref at
                # 100 C is recomputed from it.
                t_anchor = _BARRIER_DRAW_ANCHOR_T_K[key]
                k_anchor = base.k_ref * math.exp(
                    -(base.ea_kj_mol * 1000.0 / _R) * (1.0 / t_anchor - 1.0 / _T_REF_K))
                fields["k_ref"] = k_anchor / math.exp(
                    -(ea * 1000.0 / _R) * (1.0 / t_anchor - 1.0 / _T_REF_K))
        out[key] = replace(base, **fields) if fields else base
    return out

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
        "measured_elsewhere": "12 / 21 / 103 x 1e-3 /min at 150 / 160 / 170 C (hazelnut); and about "
                              "0.111 /min at 180 C from a THIRD reading, Gokmen 2012",
        "verdict": "23 000x APART on the hazelnut lumped sink, and 1.2e5x apart on Gokmen 2012 "
                   "('approximately 67 % of the HMF was lost within 10 min' with equimolar "
                   "asparagine at 180 C, giving 0.111 /min and a half-life near 6 minutes, mine). "
                   "The furanic channel already prints 'EXPECT HMF TO BE OVER-PREDICTED' and names "
                   "the empty 50-150 C window as the reason; these put two measured sizes on that "
                   "warning. Half-life at 160 C: about 33 minutes measured against about 1.5 years "
                   "shipped. THREE CAUTIONS, none of them small. (i) Gokmen's is the AMINE sink, "
                   "HMF plus asparagine, not self-degradation, so it is not the same quantity this "
                   "constant names. (ii) It is one temperature, so no barrier follows; the 90 to "
                   "180 C series in the same paper was run for ACRYLAMIDE and its 138.78 kJ/mol "
                   "must not be imported as the HMF-loss barrier. (iii) Gokmen and Goncuoglu are "
                   "co-authors, so the 1.08x agreement between their two numbers is within-group "
                   "convergence and not independent replication.",
        "anchor": "goncuoglu2016_extraction.md sec. 3 Table 1 step k26; gokmen2012_extraction.md",
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

# ---------------------------------------------------------------------------
# B39 (2026-09-11): THE FED 3-DEOXYGLUCOSONE TRIANGLE
# ---------------------------------------------------------------------------
# Mittelmaier et al. 2011 (mittelmaier2010_extraction.md) feed pure 3-DG, pure 3,4-DGE and pure
# 3-DGal at 120 C, pH 5, in water, and follow all three: the dehydration is reversible and the enone
# hydrates to either C4 epimer. The trunk carried 3-DG -> 3,4-DGE one way, 3,4-DGE -> HMF as the only
# exit, and no epimer. Five coordinates, fitted on that paper's six printed maxima and shares and on
# Zhang 2021's within-study 3,4-DDG/3-DG ratios (90-110 C): the existing k_tdg_ddg and k_ddg_hmf, and
# three new steps. ONE barrier is declared for the three new steps, equal to the forward step's
# 36.9 kJ/mol, because the source measures at one temperature and inventing three would be three
# fabricated numbers. Pre-registration: results/validation/kinetic_core_b39_prereg.md.
#
# INERT UNTIL THE FIT SHIPS: with SHIPPED_B39 False the three new steps run at k = 0 and the two
# existing constants keep their shipped literals, so every prediction is bit-for-bit what it was.
EA_FED_3DEOXY_KJ_MOL = 36.9   # the forward step's own barrier (Kocadagli 2016 Table 2 step 4), declared for all three new steps
FED_3DEOXY_COORDINATES: Tuple[str, ...] = (
    "log10_k_tdg_ddg_100C", "log10_k_ddg_tdg_100C", "log10_k_ddg_dgal_100C", "log10_k_dgal_ddg_100C", "log10_k_ddg_hmf_100C",
)
_FED_3DEOXY_KEY_OF = {
    "log10_k_tdg_ddg_100C": "k_tdg_ddg", "log10_k_ddg_tdg_100C": "k_ddg_tdg", "log10_k_ddg_dgal_100C": "k_ddg_dgal",
    "log10_k_dgal_ddg_100C": "k_dgal_ddg", "log10_k_ddg_hmf_100C": "k_ddg_hmf",
}
_MITTELMAIER = ("Mittelmaier, Funfrocken, Fenn, Berlich & Pischetsrieder 2011, Anal. Bioanal. Chem. 399:1689 (fed 3-DG / 3,4-DGE / 3-DGal, "
                "~200 uM, PD model pH 5, 120 C, 0-120 min); Zhang, Sun, Pu, Zhang, Sun & Zhao 2021, Food Sci. Nutr. 9:290 (0.3 M glucose in "
                "water, 90-110 C, within-study 3,4-DDG/3-DG ratio); mittelmaier2010_extraction.md, zhang2020_extraction.md")
#: True since wave B41 (2026-09-11). B39 fitted the triangle and did not ship (the fed peak came
#: 3x early: the 3-DG exits had no pH term); B40 added the term to both exits and did not ship (the
#: Leitzen methylglyoxal row rejected the fragmentation exit's term); B41 shipped the formic-acid
#: exit's term with this refit: kinetic_core_b41_ship_rule.json SHIP.
SHIPPED_B39: bool = True
#: The B41 fit report's optimum (cost 0.81 on twelve rows, chi2_red 0.12, fed peak at 23 min; all
#: five coordinates PINNED). Asserted equal to the report by tests/unit/test_kinetic_core_b39.py.
#: The inert "before" was {None, -30, -30, -30, None}: the furanic literals and k = 0 on the new steps.
FROZEN_B39: Mapping[str, float] = {
    "log10_k_tdg_ddg_100C": -1.9325002090039471,
    "log10_k_ddg_tdg_100C": -1.9647477052223088,
    "log10_k_ddg_dgal_100C": -1.8115924725920176,
    "log10_k_dgal_ddg_100C": -1.697993368507553,
    "log10_k_ddg_hmf_100C": -1.4503764398702186,
}
SHIPPING_FIT_REPORT = "results/validation/kinetic_core_b41_fit_report.json"


def with_fed_3deoxy(log10: Mapping[str, float]) -> Dict[str, KineticParameter]:
    """
    The fed-triangle block at arbitrary values (the fit generator's hook and the envelope's).
    ``log10`` maps a coordinate in FED_3DEOXY_COORDINATES to log10 k at 100 C; a None keeps the
    furanic module's literal for the two existing constants.
    """
    from dataclasses import replace
    from .parameters_furanic import FURANIC_PARAMETERS

    out: Dict[str, KineticParameter] = {}
    fitted = SHIPPED_B39
    for coord in FED_3DEOXY_COORDINATES:
        key = _FED_3DEOXY_KEY_OF[coord]
        value = log10.get(coord)
        if key in ("k_tdg_ddg", "k_ddg_hmf"):
            base = FURANIC_PARAMETERS[key]
            if value is None:
                out[key] = base
            else:
                flags = tuple(base.flags) + (("fitted_wave_b39",) if fitted or value is not None else ())
                out[key] = replace(base, k_ref=10.0 ** float(value), flags=flags,
                                   evidence_class="derived_from_fit_data" if fitted else base.evidence_class,
                                   source_anchor=(_MITTELMAIER if fitted else base.source_anchor))
            continue
        k_ref = 0.0 if value is None or value <= -29.0 else 10.0 ** float(value)
        transformation = {
            "k_ddg_tdg": "3,4-dideoxyglucosone-3-ene + H2O -> 3-deoxyglucosone (reverse hydration)",
            "k_ddg_dgal": "3,4-dideoxyglucosone-3-ene + H2O -> 3-deoxygalactosone (epimer hydration)",
            "k_dgal_ddg": "3-deoxygalactosone -> 3,4-dideoxyglucosone-3-ene (epimer dehydration)",
        }[key]
        out[key] = KineticParameter(
            key=key, transformation=transformation, k_ref=k_ref, ea_kj_mol=EA_FED_3DEOXY_KJ_MOL,
            evidence_class="derived_from_fit_data" if fitted else "structural_constant",
            source_anchor=_MITTELMAIER if fitted else "B39 structural step, k = 0 until the fit ships",
            dossier_anchor="mittelmaier2010_extraction.md sec. 2-4; results/validation/kinetic_core_b39_prereg.md",
            conditions="water, pH 5, 120 C, fed 200 uM; one declared barrier for the three new steps",
            ph_of_measurement=5.0, temperature_range_c=(120.0, 120.0), rate_transfer="licensed_at_measurement_ph_only",
            unit="1/min", order=1,
            flags=("b39_fed_3deoxy", "barrier_declared_equal_to_forward") + (("fitted_wave_b39",) if fitted else ("inert_until_b39_ships",)),
            note="Mittelmaier 2011 proves the step by feeding its product; the rate is this wave's fit, the barrier the forward step's.",
        )
    return out


FED_3DEOXY_PARAMETERS: Mapping[str, KineticParameter] = with_fed_3deoxy(FROZEN_B39)
FED_3DEOXY_KEYS: Tuple[str, ...] = tuple(FED_3DEOXY_PARAMETERS)

__all__ = ["FED_3DEOXY_COORDINATES", "FED_3DEOXY_PARAMETERS", "FED_3DEOXY_KEYS", "FROZEN_B39", "SHIPPED_B39", "with_fed_3deoxy", "EA_FED_3DEOXY_KJ_MOL", "DICARBONYL_KEYS", "DICARBONYL_PARAMETERS", "DICARBONYL_WISHLIST", "AQUEOUS_GLYOXAL_PARAMETERS", "AQUEOUS_GLYOXAL_KEYS",
           "AQUEOUS_GLYOXAL_COORDINATES", "AQUEOUS_GLYOXAL_CAVEAT", "FROZEN_B21", "with_aqueous_glyoxal", "GLASS_K_G_GO"]
