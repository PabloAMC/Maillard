"""
src/kinetic_core/parameters_furanic.py

EVERY CONSTANT OF THE FURANIC CHANNEL (Build Wave B7, 2026-08-29), WITH ITS
SOURCE, ITS AUDIT VERDICT AND ITS PROHIBITION.
=============================================================================

Module 6 of the kinetic-core rebuild: the HMF node and the DMHF/HEMF node.
The species live in ``species.py`` (trunk block) and ``species_sulfur.py``
(the two sulfur-coupled sinks); the topology lives in ``network.py`` and
``sulfur.py``; the analysis layer lives in ``furanic.py``.

WHAT IS DIFFERENT ABOUT THIS MODULE
-----------------------------------
Almost nothing in it is fitted. The HMF node is INGESTED WHOLE from one
declared-FIT system (Kocadagli & Gokmen 2016 JAFC, the amine-free glucose
melt) and its sink is INGESTED WHOLE from one declared-FIT dedicated study
(Hamzalioglu & Gokmen 2018). The DMHF node has exactly ONE fitted number in
the entire channel -- ``k_dpo_af`` -- calibrated on Blank, Fay, Lakner &
Schlosser 1997's isotope-dilution cells, and everything else about DMHF is a
declared assumption with its consequence written down.

THE FOUR THINGS THIS FILE REFUSES TO DO, EACH WITH THE SOURCE THAT FORBIDS IT
-----------------------------------------------------------------------------
1. **No activation energy on any FURANONE edge, from any source.** All five
   papers of the K5b cluster are SINGLE-TEMPERATURE (Blank 90 C, Wang & Ho
   120 C, Shu & Ho 160 C, Apriyantono reflux, Poisson a confounded roast ramp).
   K5b sec. 7.1: *"A barrier fitted to a single-temperature dataset is a rate
   constant wearing an Arrhenius costume."* Every furanone edge here therefore
   carries ``Ea = 0`` BY DECLARATION and the temperature dependence of DMHF is
   inherited ENTIRELY from the measured deoxyosone-formation steps upstream of
   it. That is a modelling assumption, it is banded, and the band is priced by
   re-integration rather than nominated.
2. **No fit of any kind to Shu & Ho 1988's 6.0 % GC area.** K5b sec. 8.6 names
   this failure mode by name: it is the ``thiol_addition_pentodiulose`` /
   ``cys_ribose_140C_Hofmann1998`` class, with the ADDED defect that 6.0 % is
   not even a concentration. Edge C ships with rate EXACTLY ZERO.
3. **No Hamzalioglu Ea extrapolated above 50 C.** The data span 5-50 C on
   three points with R^2 0.87 and a slope resting on the single point at which
   the authors themselves declare the pseudo-first-order assumption
   compromised. Cooking is 120-200 C. The constant is CLAMPED at 323.15 K.
4. **No branch fraction between HMF's two limbs.** K5a MUST-NOT #1. The share
   is whatever the dynamic fructose and 3-deoxyglucosone pools make it, and
   there is no scalar in this file that could encode it.

UNITS: mmol/L, minutes, Kelvin, kJ/mol -- B1's units, unchanged. Rate
constants are k at T_REF_K = 373.15 K (100 C), because that is the reference
temperature ``KineticParameter.k_at`` reparameterises to; every source below
publishes at a different reference and the conversion is done ONCE, here, in
``_at_100C``, and is unit-tested.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Optional, Tuple

from .parameters import R_KJ, T_REF_K, KineticParameter

# ---------------------------------------------------------------------------
# 0. Unit arithmetic -- the two places a factor of 1000 hides in this channel
# ---------------------------------------------------------------------------

#: A second-order constant printed in ``M^-1 day^-1`` -> ``L/(mmol*min)``.
#: 1 M^-1 = 1e-3 L/mmol and 1 day = 1440 min, so the factor is 1e-3/1440.
#: Hamzalioglu 2018 is the only second-order source in the channel and it
#: prints M^-1 day^-1, so this conversion is on the critical path for the only
#: measured HMF sink in existence.
M_INV_DAY_TO_L_PER_MMOL_MIN: float = 1.0e-3 / 1440.0


def m_inv_day_to_core_units(k: float) -> float:
    """M^-1 day^-1 -> L/(mmol*min). Unit-tested; never inlined at a call site."""
    return float(k) * M_INV_DAY_TO_L_PER_MMOL_MIN


#: ``ug per mmol of sugar`` and ``mg per mol of substrate`` ARE THE SAME UNIT.
#: 1 ug/mmol = 1e-6 g / 1e-3 mol = 1e-3 g/mol = 1 mg/mol, exactly.
#:
#: THIS IS THE CHANNEL'S SHARPEST NUMERICAL TRAP and it is recorded here rather
#: than in prose because it looks like a units difference and is not: Blank
#: 1997 prints ``ug/mmol of SUGAR`` and Wang & Ho 2008 prints ``mg/mol of
#: METHYLGLYOXAL``. The units are identical; THE DENOMINATORS ARE DIFFERENT
#: MOLECULES. Blank's 2.6 and Wang & Ho's 13 are on the same numeric scale and
#: are NOT comparable, and an ingestion that "converted" between them would
#: silently divide by the wrong pool.
UG_PER_MMOL_EQUALS_MG_PER_MOL: bool = True
MU_M_HAZARD: str = (
    "ug/mmol (Blank 1997, per mmol SUGAR) and mg/mol (Wang & Ho 2008, per mol "
    "METHYLGLYOXAL) are the SAME unit to machine precision -- 1 ug/mmol = "
    "1 mg/mol exactly. They differ only in what the denominator counts. Any "
    "code that converts between them is wrong twice: once because there is "
    "nothing to convert, and once because it has changed the reference pool "
    "without saying so."
)


def _at_100C(k_at_reference: float, ea_kj_mol: float, reference_k: float) -> float:
    """
    Re-reference a published ``k`` from its own reference temperature to 100 C.

    Sources in this channel publish at four different reference temperatures
    (Kocadagli reparameterises to T_b = 453.15 K; Hamzalioglu publishes an
    Arrhenius A, i.e. an infinite-temperature intercept; Blank and Wang & Ho
    publish a single temperature and no constant at all). Doing the algebra in
    one place, once, is the only defence against the sign error that this
    exponent invites.
    """
    return float(k_at_reference) * math.exp(
        -(float(ea_kj_mol) / R_KJ)
        * (1.0 / T_REF_K - 1.0 / float(reference_k))
    )


def _from_arrhenius_A(a: float, ea_kj_mol: float) -> float:
    """k at 100 C from an Arrhenius pre-exponential A and Ea (same unit as A)."""
    return float(a) * math.exp(-(float(ea_kj_mol) / R_KJ) / T_REF_K)


# ---------------------------------------------------------------------------
# 1. THE HMF NODE -- Kocadagli & Gokmen 2016, JAFC 10.1021/acs.jafc.6b01862
# ---------------------------------------------------------------------------
# ** THE FILE-SWAP TRAP, recorded because Amendment 12 makes it binding. **
# `data/articles/Kocada2016.pdf`   (SHORTER stem) = JAFC 10.1021/acs.jafc.6b01862
#                                  = glucose +/- NaCl caramelization = THIS source.
# `data/articles/Kocadagli2016.pdf` (LONGER stem) = Food Chem 10.1016/j.foodchem
#                                  .2016.05.150 = glucose/wheat flour = NOT this
#                                  source, and its text layer is cipher-garbled
#                                  so no grep can see its content.
#
# DECLARED ROLE (Amendment 12): the GLUCOSE system is the module's primary FIT.
# THE GLUCOSE-NaCl SYSTEM IS A HOLD-OUT AND NOT ONE NUMBER FROM IT IS IN THIS
# FILE. The firewall test in tests/unit/test_kinetic_core_b7.py greps this
# module for every NaCl-arm literal.
#
# REPARAMETERISATION, verbatim from the paper: "k(T) = k_b * exp[-(Ea/R)(1/T -
# 1/T_b)]", T_b = 180 C = 453.15 K, k_b in min^-1 x 10^3. Every k_b below is
# therefore divided by 1000 and re-referenced to 100 C by ``_at_100C``.

_KOCADAGLI_TB_K = 453.15
_KOCADAGLI_SOURCE = (
    "Kocadagli & Gokmen 2016, J. Agric. Food Chem. "
    "10.1021/acs.jafc.6b01862, Table 2 (GLUCOSE system), reparameterised "
    "Arrhenius at T_b = 180 C"
)
_KOCADAGLI_CONDITIONS = (
    "freeze-dried amorphous D-glucose melt, AMINE-FREE, 160/180/200 C, closed "
    "vials; pH not measured and not measurable in a melt; a_w not reported. "
    "The only amine-free system in the HMF corpus, which is why it is the FIT "
    "panel: its sugar chemistry cannot be confounded by the amino-acid lane."
)
_KOCADAGLI_DOSSIER = (
    "data/lit/extraction_dossiers/kocadagli2016jafc_extraction.md sec. 4 "
    "(Table 2) and sec. 3 (Table 1); consolidated in k5a_hmf_synthesis.md "
    "sec. 7.1 rows 10-12"
)


def _kocadagli(
    key: str,
    transformation: str,
    k_b_x1000: float,
    ea: float,
    ea_ci: Optional[float],
    step: int,
    *,
    flags: Tuple[str, ...] = (),
    note: str = "",
    evidence_class: str = "measured_rate",
) -> KineticParameter:
    return KineticParameter(
        key=key,
        transformation=transformation,
        k_ref=_at_100C(k_b_x1000 * 1.0e-3, ea, _KOCADAGLI_TB_K),
        ea_kj_mol=ea,
        unit="1/min",
        order=1,
        evidence_class=evidence_class,
        source_anchor=f"{_KOCADAGLI_SOURCE}, step {step}",
        dossier_anchor=_KOCADAGLI_DOSSIER,
        conditions=_KOCADAGLI_CONDITIONS,
        ph_of_measurement=None,
        temperature_range_c=(160.0, 200.0),
        rate_transfer="not_licensed",
        ea_ci95_kj_mol=ea_ci,
        flags=("amine_free_glucose_melt", "extrapolated_below_160C") + flags,
        note=note,
    )


#: Hamzalioglu & Gokmen 2018, Food Chem. 240:354-360,
#: doi 10.1016/j.foodchem.2017.07.131. THE ONLY MEASURED HMF SINK IN THE CORPUS
#: THAT SURVIVES AUDIT (K5a sec. 4.1 refuses the other seven, one of them
#: because the authors tested and rejected it themselves).
_HAMZALIOGLU_SOURCE = (
    "Hamzalioglu & Gokmen 2018, Food Chem. 240:354-360, "
    "doi 10.1016/j.foodchem.2017.07.131, Table 1 HIGH-MOISTURE arm, "
    "second-order constant derived as k = k'/[AA] with the stated "
    "[AA] = 30 umol / 1.5 mL = 20 mM"
)
_HAMZALIOGLU_CONDITIONS = (
    "aqueous 0.05 % benzoic acid, pH 3.5, HMF + amino acid, 5 / 25 / 50 C, "
    "7-day storage, pseudo-first-order at 20 mM amino acid"
)

#: THE PRE-EXPONENTIAL USED HERE IS THE DOSSIER'S REFIT, NOT THE PUBLISHED ONE.
#: Amendment 12 makes this binding, and K5a sec. 6.5 is the reason: all SIX of
#: this paper's activation energies reproduce from its own Table 1 to four
#: decimal places, and only TWO of its six pre-exponentials do. Two diagnosable
#: errors -- a sign flip on every negative intercept, and a SWAP of the
#: Coffee-Cys / Coffee-Lys pair. This is the THIRD audited case in the corpus
#: of a correct Ea bolted to a defective prefactor.
#:
#: HMF-Cys is one of the two rows whose published A does reproduce (23980.59 vs
#: refit 24115.1, 0.56 % apart). The refit value is used anyway, because the
#: instruction is to use the refit set and because using the published value
#: for the rows that happen to pass and the refit for the rows that fail would
#: be a selection nobody could audit.
_HAMZALIOGLU_A_PSEUDO_FIRST_ORDER_PER_DAY: float = 24115.1
_HAMZALIOGLU_EA_KJ_MOL: float = 29.675
_HAMZALIOGLU_AA_LOADING_M: float = 0.020

#: The second-order Arrhenius A, in core units. A_2nd = A_pseudo / [AA].
_HMF_CYS_A_CORE_UNITS: float = m_inv_day_to_core_units(
    _HAMZALIOGLU_A_PSEUDO_FIRST_ORDER_PER_DAY / _HAMZALIOGLU_AA_LOADING_M
)

#: THE HARD CEILING. Above this the constant is HELD, not extrapolated.
#: K5a sec. 7.3, "any Hamzalioglu Ea extrapolated above 50 C" -- PROHIBITED
#: DERIVATION, in the register of k3 sec. C.1.
HMF_SINK_NO_EXTRAPOLATION_ABOVE_K: float = 323.15

#: Applied by ``sulfur.sulfur_rate_constants_at``. One entry; a mapping rather
#: than a special case so that the next capped constant cannot be added by
#: editing an ``if``.
TEMPERATURE_CAP_K: Mapping[str, float] = {
    "k_hmf_cys": HMF_SINK_NO_EXTRAPOLATION_ABOVE_K,
}

#: HMF self-degradation. Hamzalioglu's model-free control: HMF ALONE in
#: aqueous acid at pH 3.5 and 5 C loses 0.9 % in 7 days. That is ONE
#: temperature, so it licenses NO activation energy, and the constant is
#: carried with Ea = 0 BY DECLARATION.
#:
#: THE CONSEQUENCE IS STATED RATHER THAN HIDDEN, because it is the module's
#: largest single weakness: with Ea = 0 this sink is 9e-7 /min at every
#: temperature, i.e. utterly negligible at cooking temperature. Combined with
#: K5a declared gap G2 ("no HMF-sink constant exists above 50 C that survives
#: audit -- the 50-150 C window is empty"), THE MODEL HAS NO EFFECTIVE HMF SINK
#: AT COOKING TEMPERATURE and must therefore be expected to OVER-PREDICT HMF.
#: That expectation is pre-registered, not discovered.
_HMF_SELF_DEGRADATION_FRACTION_7D_5C: float = 0.009
_HMF_SELF_DEGRADATION_K_PER_MIN: float = (
    -math.log(1.0 - _HMF_SELF_DEGRADATION_FRACTION_7D_5C) / (7.0 * 24.0 * 60.0)
)


# ---------------------------------------------------------------------------
# 2. THE FURANONE (DMHF) EDGES -- and the single fitted number in the channel
# ---------------------------------------------------------------------------
# Edge A  (intact skeleton, via acetylformoin)
#   pentose:  DPO + Gly -> AF          k_dpo_af   <-- THE ONLY FITTED CONSTANT
#   hexose:   ODG       -> AF          k_odg_af   <-- DERIVED, see below
#             AF        -> DMHF        k_af_dmhf  <-- DECLARED, not rate-limiting
# Edge B  (methylglyoxal, C3 + C3, acetylformoin-free)
#             2 MGO     -> DMHF        k_mgo_dmhf <-- DIGITISED PRIOR ONLY
# Edge C  (cysteine / H2S sink)
#             DMHF + H2S -> DMHFS      k_dmhf_h2s <-- EXACTLY ZERO
#
# WHY Ea = 0 ON ALL OF THEM, AND WHAT THAT ACTUALLY MEANS
# -------------------------------------------------------
# It does NOT mean "DMHF formation has no barrier". It means the barrier is
# carried by the step UPSTREAM: the 1-deoxyosone supply, whose rate constants
# are MEASURED (Martins' k_ama_odg, Ea 107 +/- 7; Kocadagli's amine-free
# Fru -> 1-DG, Ea 99.3 +/- 21.8). Setting the partition's own Ea to zero says
# "the acetylformoin branch takes a temperature-INDEPENDENT share of the
# deoxyosone flux", which is a statable, falsifiable assumption -- and it is
# the only one available, because no furanone Ea exists anywhere.
#
# THE BAND. The assumption is priced by re-integrating at partition Ea =
# +/- FURANONE_PARTITION_EA_BAND_KJ_MOL and reporting the span, exactly as B6
# prices its Q10. The band's width is itself a declared choice, not a
# measurement, and it is round on purpose.
FURANONE_PARTITION_EA_KJ_MOL: float = 0.0
FURANONE_PARTITION_EA_BAND_KJ_MOL: float = 50.0

FURANONE_EA_ASSUMPTION: Mapping[str, Any] = {
    "value_kj_mol": FURANONE_PARTITION_EA_KJ_MOL,
    "band_kj_mol": FURANONE_PARTITION_EA_BAND_KJ_MOL,
    "class": "declared_assumption",
    "basis": (
        "NO activation energy for any furanone family exists in the accessible "
        "literature (research_round3_channels.md sec. 0d, re-affirmed by "
        "k5b_dmhf_synthesis.md sec. 7.1 after reading five more papers). All "
        "five K5b papers are single-temperature. The partition's Ea is "
        "therefore set to zero by declaration and the DMHF prediction inherits "
        "its whole temperature dependence from the MEASURED deoxyosone-"
        "formation steps upstream."
    ),
    "warning": (
        "DMHF/furanone: the formation edge carries NO activation energy of its "
        "own -- none exists in the literature, on any edge, from any of the "
        "five papers in the cluster, all of which are single-temperature. Its "
        "temperature dependence is INHERITED from the measured 1-deoxyosone "
        "step. The reported interval spans a +/-50 kJ/mol partition barrier "
        "and is priced by re-integration, not nominated."
    ),
}


def furanone_band_factor(temperature_k: float, ea_kj_mol: float) -> float:
    """The multiplicative factor a partition barrier ``ea`` puts on the rate."""
    return math.exp(
        -(float(ea_kj_mol) / R_KJ)
        * (1.0 / float(temperature_k) - 1.0 / T_REF_K)
    )


#: EDGE C SHIPS AT EXACTLY ZERO. Not "small", not "estimated": zero, with the
#: species and the balanced stoichiometry present so that the channel is
#: structurally complete and so that the day a magnitude arrives there is a
#: place to put it.
EDGE_C_ZERO_BY_DECLARATION: Mapping[str, str] = {
    "rate": "0.0 L/(mmol*min), exactly, at every temperature and every pH",
    "why": (
        "Shu & Ho 1988 is the only fed-precursor DMHF + cysteine experiment in "
        "the corpus and it reports NO residual DMHF, no conversion, no mass "
        "balance and no molar yield of any product -- only GC AREA PERCENT of "
        "57 volatiles at three initial pH values in an UNBUFFERED bomb. "
        "k5b_dmhf_synthesis.md sec. 8.6 names, by name, the failure mode a "
        "constant fitted to its 6.0 % would repeat: "
        "``thiol_addition_pentodiulose`` was fitted to "
        "``cys_ribose_140C_Hofmann1998``'s 342/200 ppb, Wave S2b proved that "
        "target was a repo-internal derivation rather than a measurement, and "
        "Wave S2c had to revert it. 6.0 % is not even a concentration."
    ),
    "what_would_change_it": (
        "Haleva-Toledo et al. 1999, JAFC 47:4140-4145 -- the only identified "
        "source that quantifies DMHF inhibition by cysteine against pH IN A "
        "BUFFER. It is item 1 of the K5b ranked fetch list and it also carries "
        "5-HMF and 5-methylfurfural in the same buffer, so one acquisition "
        "would serve both nodes of this module."
    ),
}

#: Wang & Ho 2008 Figure 1, MG-alone arm, pH 5, 120 C, 60 min: 13 mg DMHF per
#: mol of methylglyoxal from a 1.4 M MG charge. DIGITISED FROM A BAR CHART by
#: the K5b wave (the plot region has no text layer at all), external-standard
#: HPLC with no recovery correction, pH hold unstated. K5b sec. 8.5:
#: "three transmission defects deep ... usable as a LOOSE (>=3x) hold-out band,
#: not as a contract." It is carried as a LEVEL PRIOR so that Edge B is present
#: and non-zero, and it is flagged on every artefact.
_WANG_HO_MG_CHARGE_MMOL_PER_L: float = 1400.0
_WANG_HO_DMHF_MG_PER_MOL_MG: float = 13.0
_WANG_HO_TIME_MIN: float = 60.0
_WANG_HO_DMHF_MW: float = 128.13


def _edge_b_k_at_120C() -> float:
    """
    ``2 MGO -> DMHF`` in L/(mmol*min), from Wang & Ho's single measured yield.

    Small-conversion algebra, stated so it can be checked: with
    d[DMHF]/dt = k [MGO]^2 and a conversion of order 1e-4, [MGO] is constant to
    four figures over the hour, so [DMHF] = k [MGO]0^2 t.
    """
    mol_dmhf_per_mol_mg = (
        _WANG_HO_DMHF_MG_PER_MOL_MG * 1.0e-3 / _WANG_HO_DMHF_MW
    )
    dmhf_mmol_per_l = mol_dmhf_per_mol_mg * _WANG_HO_MG_CHARGE_MMOL_PER_L
    return dmhf_mmol_per_l / (
        _WANG_HO_MG_CHARGE_MMOL_PER_L ** 2 * _WANG_HO_TIME_MIN
    )


EDGE_B_K_AT_120C_L_PER_MMOL_MIN: float = _edge_b_k_at_120C()


#: AF -> DMHF. A DECLARED, UNCONSTRAINED assumption in the same register as
#: ``thiol_addition_pentodiulose`` and ``deoxyosone_reduction``: acetylformoin
#: is asserted NOT to be rate-limiting, so the constant is set to a round
#: multiple of the parent deoxyosone's own measured sink. NOTHING IS FITTED TO
#: IT and the fit report carries a three-decade sensitivity sweep showing that
#: the DMHF prediction is insensitive to it, which is the only defence a number
#: with no source can have.
AF_NOT_RATE_LIMITING_MULTIPLE: float = 10.0

#: The measured competing sink whose temperature dependence Edge A inherits:
#: Martins' 1-deoxyglucosone -> acetic acid, k = 1.45 /min at 100 C,
#: Ea = 76 +/- 4 kJ/mol, independently corroborated by Knol 2010 Table 2
#: (75 +/- 10) to within 1 kJ/mol. Naming it here rather than reaching into
#: MARTINS_M4 at call time keeps the inheritance auditable.
_ODG_REFERENCE_K_100C: float = 1.45
_ODG_REFERENCE_EA: float = 76.0


# ---------------------------------------------------------------------------
# 3. THE FROZEN FIT -- one number, and the report that produced it
# ---------------------------------------------------------------------------
# ``k_dpo_af`` is the ONLY fitted constant in the furanic channel. It is
# calibrated on Blank, Fay, Lakner & Schlosser 1997 (JAFC 45:2642-2648,
# article ID JF960997I, verified from the primary document), Table 1's
# PENTOSE + GLYCINE cells -- a declared FIT set (Amendment 12: "blank1997's 39
# SIDA cells = FIT").
#
# The value below is FROZEN. ``scripts/generators/generate_kinetic_core_b7_fit.py``
# recomputes it from the declared FIT rows and writes
# ``results/validation/kinetic_core_b7_fit_report.json``; a unit test asserts
# the report and this literal agree to 1e-9. Carrying the literal here rather
# than reading the report at import time is deliberate: the furanic block hangs
# on the TRUNK, so every B1/B2/B3 generator would otherwise need a B7 report to
# exist before it could run.
FROZEN_K_DPO_AF_L_PER_MMOL_MIN: float = 4.028960042285363e-06

#: Blank's amine loading, and the basis of the pentose -> hexose transfer.
#: Blank 1997 sec. 2: 1 M pentose + 1 M amino acid. The HEXOSE edge is given
#: the PENTOSE edge's PSEUDO-FIRST-ORDER constant at that loading:
#:      k_odg_af [1/min] = k_dpo_af [L/(mmol*min)] * 1000 [mmol/L]
#:
#: THIS IS A DECLARED ASSUMPTION AND IT IS THE WEAKEST LINK IN THE DMHF NODE.
#: K5b declared gap, verbatim: "There is NO absolute hexose DMHF yield in any
#: of the five papers. The hexose intact-C6 STRUCTURE is settled twice over
#: (Wang & Ho CAMOLA; Poisson coffee CAMOLA with an internal positive control)
#: and its MAGNITUDE is measured NOWHERE." The structure is therefore honoured
#: exactly -- the hexose edge is first order in the deoxyosone alone, because
#: the C6 skeleton stays intact and needs no Strecker carbon -- and the level
#: is transferred with this one sentence and flagged on every hexose answer.
BLANK_AMINE_LOADING_MMOL_PER_L: float = 1000.0
BLANK_SUGAR_LOADING_MMOL_PER_L: float = 1000.0
BLANK_TEMPERATURE_C: float = 90.0
BLANK_TIME_MIN: float = 60.0
BLANK_PH: float = 6.0
BLANK_STATED_MAX_SD_FRACTION: float = 0.10


# ---------------------------------------------------------------------------
# 4. The parameter registry
# ---------------------------------------------------------------------------


def _declared(
    key: str,
    transformation: str,
    k_ref: float,
    ea: float,
    order: int,
    *,
    source_anchor: str,
    dossier_anchor: str,
    conditions: str,
    ph: Optional[float],
    t_range: Tuple[float, float],
    evidence_class: str,
    flags: Tuple[str, ...],
    note: str,
) -> KineticParameter:
    return KineticParameter(
        key=key,
        transformation=transformation,
        k_ref=k_ref,
        ea_kj_mol=ea,
        unit="1/min" if order == 1 else "L/(mmol*min)",
        order=order,
        evidence_class=evidence_class,
        source_anchor=source_anchor,
        dossier_anchor=dossier_anchor,
        conditions=conditions,
        ph_of_measurement=ph,
        temperature_range_c=t_range,
        rate_transfer="not_licensed",
        flags=flags,
        note=note,
    )


MEASURED_FURANIC: Mapping[str, KineticParameter] = {
    # ---- the 3-DG limb ---------------------------------------------------
    "k_glc_tdg": _kocadagli(
        "k_glc_tdg", "Glc -> 3-deoxyglucosone (amine-free caramelization)",
        4.19, 107.2, 52.7, 3,
        flags=("hpd_58_percent_of_estimate", "amine_independent_route"),
        note=(
            "THE AMINE-FREE ENTRY TO THE 3-DG LIMB, and the reason the module "
            "can answer a glucose-only pot at all: B1's trunk makes "
            "3-deoxyglucosone ONLY from the Amadori compound, so before B7 a "
            "sugar-only system had no 3-DG and therefore no 3-DG limb. Its "
            "95 % HPD is 58 % of the estimate on the k_b side, which is large "
            "but does not reach the k3 sec. C.6 refusal threshold (interval "
            "spanning zero); the flag travels."
        ),
    ),
    "k_tdg_ddg": _kocadagli(
        "k_tdg_ddg", "3-deoxyglucosone -> 3,4-dideoxyglucosone",
        30.5, 36.9, 6.3, 4,
        flags=("plateau_caveat_2.1x_over_40C", "product_semi_quantitated"),
        note=(
            "THE RATE-DETERMINING STEP OF THE 3-DG LIMB, established "
            "independently in TWO matrices: 'almost 5 times higher ... "
            "3,4-dideoxyglucosone formation from 3-deoxyglucosone was the rate "
            "determining step' (Goncuoglu Tas 2017, hazelnut) and 'dehydration "
            "of 3-deoxyglucosone to 3,4-dideoxyglucosone is a rate-limiting "
            "step' (Kocadagli Food Chem 2016, wheat flour). K5a C3. "
            "Its Ea moves k only 2.1x over 40 C, which is the Ma-2022 plateau "
            "caveat's neighbourhood; K5a sec. 6.3 grades it USE-Q rather than "
            "refusing it, and the flag travels."
        ),
    ),
    "k_ddg_hmf": _kocadagli(
        "k_ddg_hmf", "3,4-dideoxyglucosone -> HMF",
        119.0, 0.0, None, 5,
        evidence_class="bounded_from_a_timescale_bracket",
        flags=("ea_author_fixed_to_zero", "reactant_semi_quantitated",
               "no_usable_ea_exists_for_this_edge"),
        note=(
            "Ea = 0 IS THE AUTHORS' OWN VALUE, not this module's choice. Their "
            "footnote, verbatim: 'Zero activation energy (Ea) indicates that "
            "the reaction rate constant (k) of the elementary step does not "
            "follow Arrhenius equation and the Ea was fixed to zero during "
            "parameter estimation.' Their k for this step runs 160 -> 110 -> "
            "137 (x1e-3 /min) across 160 -> 180 -> 200 C: non-monotonic, and a "
            "3-point refit gives Ea = -7.0 kJ/mol with R^2 = 0.189. "
            "K5a sec. 7.3 refuses every published Ea for this edge in every "
            "paper of the cluster; there is no usable one anywhere. "
            "evidence_class is 'bounded_from_a_timescale_bracket' rather than "
            "'measured_rate' precisely because a constant with no temperature "
            "dependence is a bracket at the measured temperatures and nothing "
            "more."
        ),
    ),
    "k_tdg_mgo": _kocadagli(
        "k_tdg_mgo", "3-deoxyglucosone -> methylglyoxal + C3 residue",
        304.0, 84.8, 6.9, 11,
        flags=("amine_independent_route", "mgo_parent_is_matrix_dependent"),
        note=(
            "The amine-free methylglyoxal source, and the feed for DMHF's "
            "Edge B in a sugar-only pot. K5a sec. 5 row 3 records that THE "
            "PARENT OF MGO SWITCHES BETWEEN MATRICES IN THE SAME LAB -- 3-DG "
            "in the melt and in juice, 1-DG in flour, hazelnut and nuts/seeds "
            "-- and the authors reconcile it themselves: 'the main source of "
            "methylglyoxal may quantitatively depend on the amount of "
            "precursor alpha-dicarbonyl compound formed'. B1's trunk already "
            "carries the Amadori -> MGO route; this adds the melt's, so both "
            "parents exist and the split is set by the pools rather than by a "
            "constant."
        ),
    ),
    # ---- the fructose limb: A PAIR, never separable ----------------------
    "k_fru_int": _kocadagli(
        "k_fru_int", "Fru -> undetermined cyclic intermediate (Int)",
        330.0, 100.4, 6.6, 6,
        flags=("intermediate_pool_unmeasured", "ratio_only_with_k_int_hmf"),
        note=(
            "The single best-behaved HMF-lane Ea in the whole Gokmen corpus: "
            "the published global fit gives 100.4 +/- 6.6 and an independent "
            "3-point refit of the paper's own Table 1 gives 100.5 with "
            "R^2 = 1.000. It is ALSO the only one: K5a sec. 6.2 shows that in "
            "all three REAL-MATRIX systems this same limb collapses (wheat "
            "flour Ea = -97.6, R^2 = 0.322; hazelnut non-monotonic; five "
            "nuts/seeds with k = 0 in 6 of 10 cells). "
            "RATIO-ONLY. [Int] is never measured, so only the PRODUCT "
            "k_int_hmf*[Int] is identified by the data and this constant is "
            "meaningful only paired with its partner. K5a C2."
        ),
    ),
    "k_int_hmf": _kocadagli(
        "k_int_hmf", "Int -> HMF (the fructose-limb terminal step)",
        1.84, 151.4, 34.3, 7,
        flags=("intermediate_pool_unmeasured", "ratio_only_with_k_fru_int"),
        note=(
            "The fructose limb's RATE-DETERMINING STEP -- the mirror image of "
            "the 3-DG limb, whose RDS is the FIRST step (K5a C3/C4). Gursul "
            "Aktag, verbatim: 'formation of FFC from fructose was found to be "
            "the fast step and the rate determining step was the HMF formation "
            "from FFC.' "
            "Ea 151.4 +/- 34.3 published, 145.0 by refit, R^2 = 1.000. Label "
            "it 'apparent Ea for a lumped fructose dehydration over an "
            "UNMEASURED intermediate pool in an amine-free glucose melt' and "
            "NEVER 'the HMF barrier'."
        ),
    ),
    "k_fru_odg": _kocadagli(
        "k_fru_odg", "Fru -> 1-deoxyglucosone (2,3-enolisation, amine-free)",
        2.11, 99.3, 21.8, 8,
        flags=("amine_independent_route", "odg_parent_is_matrix_dependent"),
        note=(
            "THE AMINE-FREE ENTRY TO THE FURANONE LIMB. K5a sec. 5 row 4: the "
            "parent of 1-DG has FOUR different answers in four matrices from "
            "ONE lab -- fructose only in the melt, the Amadori compound only "
            "in flour, both in hazelnut, and not detected at all in acidic "
            "juice. B1's trunk carries the Amadori parent; this adds the "
            "melt's. Without it a glucose/alanine pot has no 1-deoxyosone at "
            "all (the trunk's only amine is glycine) and DMHF would be a "
            "structural zero for a bookkeeping reason rather than a chemical "
            "one."
        ),
    ),
    # ---- the HMF sinks ---------------------------------------------------
    "k_hmf_self": _declared(
        "k_hmf_self", "HMF -> unassigned products (self-degradation)",
        _HMF_SELF_DEGRADATION_K_PER_MIN, 0.0, 1,
        source_anchor=(
            "Hamzalioglu & Gokmen 2018, sec. 3.1, model-free control: HMF "
            "alone, no amino acid, 0.9 % lost in 7 days at 5 C and pH 3.5"
        ),
        dossier_anchor=(
            "data/lit/extraction_dossiers/hamzalioglu2018_extraction.md "
            "sec. 7 / S1; k5a_hmf_synthesis.md C7"
        ),
        conditions="aqueous 0.05 % benzoic acid, pH 3.5, 5 C, 7 days",
        ph=3.5,
        t_range=(5.0, 5.0),
        evidence_class="measured_rate",
        flags=("single_temperature_no_ea_licensed", "ea_zero_by_declaration",
               "negligible_at_cooking_temperature"),
        note=(
            "ONE TEMPERATURE, so NO activation energy is licensed and Ea is "
            "zero by declaration. THE CONSEQUENCE IS THE MODULE'S LARGEST "
            "WEAKNESS AND IT IS PRE-REGISTERED, NOT DISCOVERED: at 9e-7 /min "
            "this sink removes 3e-5 of the HMF over a 30-minute bake, so the "
            "model effectively has NO HMF sink at cooking temperature. K5a "
            "declared gap G2 says why nothing better exists -- Hamzalioglu "
            "covers 5-50 C, the hazelnut first-order sink covers 150-160 C but "
            "is a lumped decay with no amine dependence, and THE 50-150 C "
            "WINDOW IS EMPTY. Expect HMF to be OVER-predicted. "
            "K5a C7 also warns that this sink and the amine sink are NOT "
            "additive: 'presence of amino compound in reaction medium limited "
            "the decomposition of HMF ... thus dimerization of HMF was "
            "inhibited'. Summing them, as this network does, therefore "
            "over-counts -- in the direction that makes the model's HMF "
            "over-prediction WORSE, not better, so the simplification cannot "
            "be flattering it."
        ),
    ),
    "k_hmf_cys": _declared(
        "k_hmf_cys", "HMF + cysteine -> adduct pool (the measured HMF sink)",
        _from_arrhenius_A(_HMF_CYS_A_CORE_UNITS, _HAMZALIOGLU_EA_KJ_MOL),
        _HAMZALIOGLU_EA_KJ_MOL, 2,
        source_anchor=_HAMZALIOGLU_SOURCE,
        dossier_anchor=(
            "data/lit/extraction_dossiers/hamzalioglu2018_extraction.md "
            "sec. 4 and sec. 6 (the refit); k5a_hmf_synthesis.md sec. 4.2 and "
            "7.1 rows 1, 4, 6"
        ),
        conditions=_HAMZALIOGLU_CONDITIONS,
        ph=3.5,
        t_range=(5.0, 50.0),
        evidence_class="measured_rate",
        flags=("refit_prefactor_not_published_one", "clamped_above_50C",
               "stoichiometry_not_1_to_1", "competitive_not_first_order",
               "selectivity_collapses_in_a_low_moisture_matrix"),
        note=(
            "THE CORPUS'S FIRST GENUINE SECOND-ORDER HMF-SINK CONSTANT: "
            "3.95 / 5.15 / 23.3 M^-1 day^-1 at 5 / 25 / 50 C, derived from the "
            "published pseudo-first-order k' and the stated 20 mM loading. "
            "Four caveats travel with it and each changes a prediction: "
            "(i) AMINE IDENTITY OUTRANKS TEMPERATURE: cysteine at 5 C removes "
            "more HMF than lysine at 50 C, so a 45 C temperature gap loses to "
            "the choice of amine. This module carries the CYSTEINE row only, "
            "because cysteine is the only amine the sulfur lane tracks. "
            "(ii) That water-phase selectivity COLLAPSES to near-parity in a "
            "low-moisture roasted-coffee matrix -- a declared HOLD-OUT whose "
            "numbers are deliberately NOT written here, and which this module "
            "has no moisture term with which to predict. "
            "(iii) The stoichiometry is NOT 1:1 (four adducts confirmed). "
            "(iv) At 160-200 C the edge is OUT-COMPETED for the amine pool by "
            "GO, 1-DG and MGO (Kocadagli's Table 1: 616 / 373 / 115 against a "
            "REJECTED HMF + AA step), so it is written bimolecular and "
            "competitive rather than as a fixed first-order decay -- K5a "
            "MUST-NOT #4."
        ),
    ),
    # ---- the furanone edges ---------------------------------------------
    "k_odg_af": _declared(
        "k_odg_af", "1-deoxyglucosone -> acetylformoin (HEXOSE, intact C6)",
        FROZEN_K_DPO_AF_L_PER_MMOL_MIN * BLANK_AMINE_LOADING_MMOL_PER_L,
        FURANONE_PARTITION_EA_KJ_MOL, 1,
        source_anchor=(
            "DERIVED, not measured. The pentose edge's PSEUDO-FIRST-ORDER "
            "constant at Blank 1997's 1 M amine loading, transferred to the "
            "hexose series. Structure from Wang & Ho 2008 sec. 4.1 (CAMOLA: "
            "[13C6]:[12C6] = 1:1 with [13C1]-[13C5] ABSENT, i.e. 100 % intact "
            "C6) and Poisson 2019 T9/T10 (intact share 87.3-100 % at all nine "
            "roast times, with 2,3-butanedione's intact share collapsing "
            "25.4 -> 0.4 % in the SAME runs as an internal positive control)"
        ),
        dossier_anchor=(
            "k5b_dmhf_synthesis.md sec. 6a rows A19-A22 and sec. 10 "
            "('There is NO absolute hexose DMHF yield in any of the five "
            "papers')"
        ),
        conditions=(
            "transferred; no hexose DMHF magnitude is measured anywhere in the "
            "K5b cluster"
        ),
        ph=None,
        t_range=(90.0, 90.0),
        evidence_class="bounded_from_a_timescale_bracket",
        flags=("hexose_magnitude_unmeasured", "declared_transfer_from_pentose",
               "ea_zero_by_declaration", "widest_band_in_the_module"),
        note=(
            "FIRST ORDER IN THE DEOXYOSONE ALONE, and that is the structure "
            "rather than a convenience: the hexose route keeps the sugar's own "
            "C6 skeleton and needs no Strecker carbon, which is exactly what "
            "the two CAMOLA experiments establish. The LEVEL is the assumption."
        ),
    ),
    "k_af_dmhf": _declared(
        "k_af_dmhf", "acetylformoin -> DMHF (reduction)",
        AF_NOT_RATE_LIMITING_MULTIPLE * _ODG_REFERENCE_K_100C,
        _ODG_REFERENCE_EA, 1,
        source_anchor=(
            "DECLARED, UNCONSTRAINED. Route from Hofmann & Schieberle 2001 via "
            "Wang & Ho 2008 (ref 17) and Poisson 2019 (ref 19); NO rate, "
            "conversion or time course is published for it anywhere"
        ),
        dossier_anchor=(
            "k5b_dmhf_synthesis.md sec. 6a row A26 ('rate constant / Ea / time "
            "course: NONE anywhere') and R6"
        ),
        conditions="none -- this constant has no measurement behind it",
        ph=None,
        t_range=(90.0, 90.0),
        evidence_class="bounded_from_a_timescale_bracket",
        flags=("unconstrained_declared_assumption", "not_rate_limiting",
               "sensitivity_swept_in_the_fit_report"),
        note=(
            "IN THE SAME REGISTER AS ``thiol_addition_pentodiulose`` AND "
            "``deoxyosone_reduction``: a number with no source. The assumption "
            "it encodes is 'acetylformoin does not accumulate', set as a round "
            "10x the parent deoxyosone's own measured sink so the value is "
            "legible. NOTHING IS FITTED TO IT; the B7 fit report sweeps it "
            "over three decades and reports the DMHF sensitivity, which is the "
            "only defence available to a constant of this class."
        ),
    ),
    "k_mgo_dmhf": _declared(
        "k_mgo_dmhf", "2 methylglyoxal -> DMHF (Edge B, C3+C3)",
        EDGE_B_K_AT_120C_L_PER_MMOL_MIN, FURANONE_PARTITION_EA_KJ_MOL, 2,
        source_anchor=(
            "Wang & Ho 2008, JAFC 56:7405-7409, Figure 1, MG-alone arm: "
            "13 mg DMHF per mol MG from 1.4 M MG in 0.5 M phosphate at pH 5, "
            "120 C, 60 min"
        ),
        dossier_anchor=(
            "data/lit/extraction_dossiers/wang2008_extraction.md sec. 4.2; "
            "k5b_dmhf_synthesis.md sec. 6b rows B1, B6-B11 and sec. 8.5"
        ),
        conditions="1.4 M methylglyoxal, 0.5 M phosphate, pH 5 set-point "
                   "(hold not stated), 120 C, 60 min",
        ph=5.0,
        t_range=(120.0, 120.0),
        evidence_class="bounded_from_a_timescale_bracket",
        flags=("digitised_from_a_bar_chart", "level_prior_only",
               "single_temperature_no_ea_licensed", "ea_zero_by_declaration",
               "acetylformoin_free_by_construction"),
        note=(
            "STRUCTURALLY SEPARATE FROM EDGE A, and that separation is a hard "
            "measured result rather than a modelling choice: feeding "
            "[13C6]glucose + [12C3]MG gives a 4:1 mixture of [13C6]- and "
            "[12C6]-DMHF and NO [12C6]acetylformoin, so the two routes share "
            "no intermediate; and only [13C6] and [12C6] isotopomers are ever "
            "seen, never [13C3], so no DMHF is assembled from one "
            "hexose-derived C3 plus one MG-derived C3. This resolves an "
            "ambiguity Poisson 2019 explicitly leaves open (its pathway b vs "
            "c). The network honours it by construction: acetylformoin appears "
            "in NO reaction that has methylglyoxal as a reactant. "
            "THE LEVEL IS A PRIOR AND NOTHING ELSE. K5b sec. 8.5: 'three "
            "transmission defects deep -- usable as a loose (>=3x) hold-out "
            "band, not as a contract.'"
        ),
    ),
    "k_dmhf_h2s": _declared(
        "k_dmhf_h2s", "DMHF + H2S -> 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone",
        0.0, 0.0, 2,
        source_anchor=(
            "Shu & Ho 1988, JAFC 36:801-803, Table I -- STRUCTURE ONLY. The "
            "paper reports a GC area percent and no consumption of any kind"
        ),
        dossier_anchor=(
            "data/lit/extraction_dossiers/shu1988_extraction.md; "
            "k5b_dmhf_synthesis.md sec. 3.4, 6c rows C1-C12 and sec. 8.6"
        ),
        conditions="0.1 M DMHF + 0.1 M cysteine-HCl, UNBUFFERED, initial pH "
                   "2.2 / 5.1 / 7.1, 160 C, 30 min, Parr bomb",
        ph=None,
        t_range=(160.0, 160.0),
        evidence_class="structural_constant",
        flags=("rate_is_exactly_zero_by_declaration", "prohibited_to_fit",
               "structure_only"),
        note=(
            "SEE ``EDGE_C_ZERO_BY_DECLARATION``. The edge exists so that the "
            "channel is structurally complete -- it is the exact C6 "
            "counterpart of the repo's existing C5 "
            "``furanone_reductive_sulfhydrylation`` -- and it moves nothing. "
            "Its pH behaviour is real and measured as a SHAPE (markers at "
            "6.0 % and 5.8 % at pH 2.2 and 5.1, not detected at 7.1) but the "
            "authors' own reading is that the pH-7.1 zeros are SECONDARY "
            "CONSUMPTION rather than absence of coupling, so even the shape "
            "tests net survival and not coupling rate."
        ),
    ),
}


#: The pentose Edge A constant, built from the frozen fit. Second order because
#: the pentose route needs the Strecker aldehyde's carbon: C5 + C1 -> C6.
#: Blank 1997 Table 5 measures the apparent reaction order in amino acid at
#: n = 0.35-0.71 -- SUB-LINEAR -- so writing it first order in the amine is a
#: KNOWN, DECLARED structural mis-specification, recorded as hold-out H8 and
#: expected to fail. Recording it is the point.
def _k_dpo_af(value: float) -> KineticParameter:
    return _declared(
        "k_dpo_af",
        "1-deoxypentosone + glycine -> acetylformoin (PENTOSE, Strecker C1)",
        value, FURANONE_PARTITION_EA_KJ_MOL, 2,
        source_anchor=(
            "FITTED -- the only fitted constant in the furanic channel. "
            "Calibrated on Blank, Fay, Lakner & Schlosser 1997, JAFC "
            "45:2642-2648 (article ID JF960997I), Table 1 PENTOSE + GLYCINE "
            "HDMF cells: arabinose 5.1, xylose 2.6, ribose 3.6 ug per mmol of "
            "sugar, isotope-dilution assay against synthesised [13C2]HDMF, "
            "max SD <= 10 %, 90 C / 1 h / 0.2 M phosphate / pH 6 / 1:1"
        ),
        dossier_anchor=(
            "data/lit/extraction_dossiers/blank1997_extraction.md sec. 2; "
            "k5b_dmhf_synthesis.md sec. 6a rows A1-A5 and its FIT proposal "
            "table row F1; "
            "FIT_HOLDOUT_DECLARATION.md Amendment 12"
        ),
        conditions=(
            "1 M pentose + 1 M glycine, 0.2 M phosphate, pH 6.0 set-point "
            "(hold not stated, no final pH reported), 90 C, 1 h"
        ),
        ph=BLANK_PH,
        t_range=(BLANK_TEMPERATURE_C, BLANK_TEMPERATURE_C),
        evidence_class="derived_from_fit_data",
        flags=("the_only_fitted_constant_in_the_channel",
               "single_temperature_no_ea_licensed", "ea_zero_by_declaration",
               "first_order_in_amine_but_measured_n_is_0.35_to_0.71",
               "one_generic_pentose_cannot_resolve_the_sugar_ordering"),
        note=(
            "TWO DECLARED LIMITATIONS, both of which will show as fit "
            "residuals rather than being hidden. "
            "(i) The sulfur lane carries ONE generic aldopentose, so the "
            "1.96x arabinose > ribose > xylose spread Blank measures cannot be "
            "reproduced and is reported as a residual -- the same treatment "
            "B2 gives Hofmann's 1.38x ribose/xylose gap. "
            "(ii) The step is first order in the amine and Blank's Table 5 "
            "measures n = 0.35-0.71. That is hold-out H8, and it is declared "
            "in advance as a KNOWN structural defect rather than discovered "
            "later."
        ),
    )


FITTED_FURANIC: Mapping[str, KineticParameter] = {
    "k_dpo_af": _k_dpo_af(FROZEN_K_DPO_AF_L_PER_MMOL_MIN),
}


def with_fitted_furanic(k_dpo_af: float) -> Dict[str, KineticParameter]:
    """The furanic fitted block at an arbitrary value, for the fit generator."""
    fitted = _k_dpo_af(float(k_dpo_af))
    hexose = MEASURED_FURANIC["k_odg_af"]
    return {
        "k_dpo_af": fitted,
        # The hexose edge is DERIVED from the pentose one and must move with it
        # or the declared transfer rule silently stops holding during a fit.
        "k_odg_af": KineticParameter(
            **{
                **{
                    field: getattr(hexose, field)
                    for field in hexose.__dataclass_fields__
                },
                "k_ref": float(k_dpo_af) * BLANK_AMINE_LOADING_MMOL_PER_L,
            }
        ),
    }


FURANIC_PARAMETERS: Mapping[str, KineticParameter] = {
    **MEASURED_FURANIC,
    **FITTED_FURANIC,
}


# ---------------------------------------------------------------------------
# 5. THE PROHIBITION REGISTER -- read before adding anything to this file
# ---------------------------------------------------------------------------

PROHIBITED_DERIVATIONS: Mapping[str, str] = {
    "any_furanone_activation_energy": (
        "All five papers of the DMHF cluster are SINGLE-TEMPERATURE. A barrier "
        "fitted to any of them would absorb every unmodelled factor in the "
        "system and would not transfer. K5b sec. 7.1, and "
        "research_round3_channels.md sec. 0d before it."
    ),
    "thiol_addition_dmhf_fitted_to_shu_6_percent": (
        "K5b sec. 8.6, by name. Shu & Ho's 6.0 % is a GC AREA PERCENT with no "
        "internal standard, in an unbuffered bomb with initial-only pH labels, "
        "and it is not a concentration. This is the "
        "``thiol_addition_pentodiulose`` failure class with an extra defect."
    ),
    "hamzalioglu_ea_above_50C": (
        "Data span 5-50 C; cooking is 120-200 C. Three points, R^2 0.80-0.91, "
        "slope resting on the single point at which the authors declare the "
        "pseudo-first-order assumption compromised. The constant is CLAMPED, "
        "not extrapolated. K5a sec. 7.3, register of k3 sec. C.1."
    ),
    "gursul_aktag_2020_table_2_activation_energies": (
        "ALL 43 REFUSED. Not one reproduces from the paper's own Table 1; six "
        "are mathematically underivable (k fixed at 0.0 at 27 C, so ln k is "
        "-infinity); the R^2 column is impossible for a 2-point fit; seven "
        "sign flips. K5a sec. 7.3."
    ),
    "goncuoglu_tas_2017_activation_energies": (
        "Range 0-1174 kJ/mol with six zeros, values live in a Table S4 that is "
        "not on disk, and the authors disclaim Arrhenius twice. K5a sec. 7.3."
    ),
    "lee_2024_arrhenius_constants": (
        "Its ka/kd are matrix-to-vapour TRANSPORT and oven-wall deposition, "
        "not chemistry. research_round3_channels.md sec. A.3.4. Only the "
        "scheme sentence is quotable."
    ),
    "lee_and_nagy_1990_36x_fructose_advantage": (
        "Cited, never read, and contradicted 6-9x by the measurement of the "
        "paper that cites it under identical stated conditions (Agcam 2022: "
        "4-6x). K5a sec. 7.3."
    ),
    "perez_locas_yaylayan_90_10_fructose_glucose_split": (
        "A branch fraction quoted second-hand in a paper that measured neither "
        "3-DG nor the fructofuranosyl cation. Attributing it to Agcam 2022 "
        "would be a citation error, and hard-coding it would violate K5a "
        "MUST-NOT #1."
    ),
    "any_fixed_branch_fraction_between_the_two_hmf_limbs": (
        "Six papers, six matrices, verdicts spanning 'fructose limb dominant' "
        "to '3-DG limb dominant by infinity', and every paper that names WHY "
        "its limb wins names a SUPPLY reason rather than a rate constant. K5a "
        "sec. 3.1 Rule 1 and MUST-NOT #1. There is no scalar in this module "
        "that could encode one, and a unit test asserts it."
    ),
    "comparing_a_constant_on_Int_with_one_on_3DG": (
        "[Int] is unmeasured, so its constants are identified only up to the "
        "pool scale and are not commensurable in magnitude with constants on a "
        "measured species. K5a C2 / MUST-NOT #3."
    ),
    "any_dft_barrier": (
        "Owner policy, repository-wide. No computational evidence class "
        "exists; ``parameters.assert_no_dft`` enforces it."
    ),
}


def assert_no_dft_furanic(
    registry: Mapping[str, KineticParameter] = FURANIC_PARAMETERS,
) -> None:
    """Owner policy: no DFT-derived barrier may enter this registry."""
    from .parameters import assert_no_dft

    assert_no_dft(registry)


assert_no_dft_furanic()


def furanic_registry_metadata() -> Dict[str, Any]:
    """Everything a report needs to say where each furanic constant came from."""
    return {
        "module": "kinetic_core.parameters_furanic",
        "wave": "B7",
        "n_parameters": len(FURANIC_PARAMETERS),
        "n_fitted": len(FITTED_FURANIC),
        "fitted_keys": sorted(FITTED_FURANIC),
        "parameters": {
            key: parameter.as_metadata()
            for key, parameter in sorted(FURANIC_PARAMETERS.items())
        },
        "prohibited_derivations": dict(PROHIBITED_DERIVATIONS),
        "edge_c": dict(EDGE_C_ZERO_BY_DECLARATION),
        "furanone_ea_assumption": dict(FURANONE_EA_ASSUMPTION),
        "temperature_cap_K": dict(TEMPERATURE_CAP_K),
        "mu_m_hazard": MU_M_HAZARD,
    }


__all__ = [
    "AF_NOT_RATE_LIMITING_MULTIPLE",
    "BLANK_AMINE_LOADING_MMOL_PER_L",
    "BLANK_PH",
    "BLANK_STATED_MAX_SD_FRACTION",
    "BLANK_SUGAR_LOADING_MMOL_PER_L",
    "BLANK_TEMPERATURE_C",
    "BLANK_TIME_MIN",
    "EDGE_B_K_AT_120C_L_PER_MMOL_MIN",
    "EDGE_C_ZERO_BY_DECLARATION",
    "FITTED_FURANIC",
    "FROZEN_K_DPO_AF_L_PER_MMOL_MIN",
    "FURANIC_PARAMETERS",
    "FURANONE_EA_ASSUMPTION",
    "FURANONE_PARTITION_EA_BAND_KJ_MOL",
    "FURANONE_PARTITION_EA_KJ_MOL",
    "HMF_SINK_NO_EXTRAPOLATION_ABOVE_K",
    "MEASURED_FURANIC",
    "MU_M_HAZARD",
    "M_INV_DAY_TO_L_PER_MMOL_MIN",
    "PROHIBITED_DERIVATIONS",
    "TEMPERATURE_CAP_K",
    "assert_no_dft_furanic",
    "furanic_registry_metadata",
    "furanone_band_factor",
    "furanic_registry_metadata",
    "m_inv_day_to_core_units",
    "with_fitted_furanic",
]
