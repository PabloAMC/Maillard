"""
src/kinetic_core/parameters.py

THE SINGLE PARAMETER REGISTRY OF THE MASS-ACTION KINETIC CORE.
==============================================================

Every rate constant the network can use is declared here, once, with:

  * its VALUE and UNIT in the module's native units (mmol/L, min, K, kJ/mol);
  * its SOURCE ANCHOR down to the table and the dossier line;
  * the CONDITIONS it was measured under, including ``ph_of_measurement``;
  * an ``evidence_class`` and a ``rate_transfer`` licence flag, both of which
    are carried through to runtime metadata by ``as_metadata()``.

THREE STANDING POLICIES, ENFORCED HERE
--------------------------------------
1. **NO DFT.** No barrier in this file is computational. ``evidence_class`` is
   restricted to the four measured/derived classes below and there is no
   pathway by which a computed barrier can enter; ``assert_no_dft()`` is called
   at import and is also asserted by the unit tests. (Owner policy, pinned:
   never run or trust DFT barriers; measured kinetics only.)

2. **pH IS FIXED, NOT MODELLED.** There is NO pH term in this module, because
   the corpus contains no licensed pH dependency for the trunk. What it
   contains is the opposite: five independent sign-crossings
   (k3_final_parameter_inventory.md sec. B.2, "NO family-level pH term can
   pass"), and Ajandouz 2008's explicit refusal to transfer BROWNING
   activation energies off its own pH (sec. C.13, quoted in A.1.1). Every
   parameter therefore carries ``ph_of_measurement`` and the network is valid
   only at that pH. ``check_ph_homogeneity()`` refuses to assemble a network
   whose operative parameters were measured at different pH values.

3. **NO WATER-ACTIVITY TERM.** Nothing in the Module 4 fit corpus varies a_w:
   the entire Martins/Knol trunk is dilute aqueous. An a_w dependence would be
   invented, so there is none. ``AW_OF_MEASUREMENT`` records the single value
   the whole registry sits at.

WHAT IS *NOT* IN THIS FILE
--------------------------
The two Module 4 HOLD-OUT numbers -- Martins' epsilon = 0.64 L/(mmol*cm) and
Knol's epsilon = 282 L/(mol*cm) -- are deliberately absent. They live only in
the hold-out scorer, are read only at scoring time, and cannot be reached from
the runtime network. See ``docs/reference/FIT_HOLDOUT_DECLARATION.md`` D.6.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Mapping, Optional, Tuple

# Gas constant, kJ/(mol K) -- matches the kJ/mol unit Ea is carried in.
R_KJ = 8.314462618e-3

#: Reference temperature of the reparameterised Arrhenius form, 100 C.
#: This is Martins' own T_av (thesis eqs 5.3-5.6) and the midpoint of the
#: 80/100/120 C fit corpus. Over an 80-120 C window A and Ea are so correlated
#: that a fitted A is individually meaningless; (k_ref, Ea) is near-orthogonal.
T_REF_K = 373.15

#: The one water activity the whole registry sits at. Dilute aqueous solution.
AW_OF_MEASUREMENT = 1.0

#: Allowed evidence classes. There is no computational class, by policy.
EVIDENCE_CLASSES = (
    "measured_rate",          # a rate constant printed in a paper
    "measured_activation_energy",  # an Ea printed in a paper
    "derived_from_fit_data",  # estimated here by least squares on FIT rows
    "structural_constant",    # a stoichiometry or atom count, not a rate
)


@dataclass(frozen=True)
class KineticParameter:
    """One rate constant, with everything needed to defend it."""

    key: str
    transformation: str
    #: k at T_REF_K, in `unit`
    k_ref: Optional[float]
    ea_kj_mol: Optional[float]
    unit: str
    #: 1 = first order, 2 = second order (bimolecular)
    order: int
    evidence_class: str
    source_anchor: str
    dossier_anchor: str
    conditions: str
    ph_of_measurement: Optional[float]
    temperature_range_c: Tuple[float, float]
    #: "licensed" | "not_licensed" | "licensed_at_measurement_ph_only"
    rate_transfer: str
    k_ref_ci95: Optional[float] = None
    ea_ci95_kj_mol: Optional[float] = None
    aw_of_measurement: float = AW_OF_MEASUREMENT
    flags: Tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.evidence_class not in EVIDENCE_CLASSES:
            raise ValueError(
                f"{self.key}: evidence_class {self.evidence_class!r} not in "
                f"{EVIDENCE_CLASSES}. Computational (DFT-derived) barriers are "
                f"refused by owner policy."
            )
        if self.order not in (1, 2):
            raise ValueError(f"{self.key}: order must be 1 or 2, got {self.order}")
        expected_unit = "1/min" if self.order == 1 else "L/(mmol*min)"
        if self.unit != expected_unit:
            raise ValueError(
                f"{self.key}: order {self.order} requires unit {expected_unit!r}, "
                f"got {self.unit!r}"
            )

    def k_at(self, temperature_k: float) -> float:
        """k(T) = k_ref * exp(-(Ea/R)(1/T - 1/T_ref)), same unit as k_ref."""
        import math

        if self.k_ref is None or self.ea_kj_mol is None:
            raise ValueError(f"{self.key}: parameter is not numerically populated")
        return float(self.k_ref) * math.exp(
            -(float(self.ea_kj_mol) / R_KJ)
            * (1.0 / float(temperature_k) - 1.0 / T_REF_K)
        )

    def as_metadata(self) -> Dict[str, Any]:
        """Runtime metadata: every provenance flag travels with the number."""
        payload = asdict(self)
        payload["flags"] = list(self.flags)
        payload["temperature_range_c"] = list(self.temperature_range_c)
        payload["reference_temperature_K"] = T_REF_K
        payload["dft_derived"] = False
        return payload


# ---------------------------------------------------------------------------
# THE MEASURED BACKBONE -- Martins & van Boekel 2005, Table 2 (model M4)
# ---------------------------------------------------------------------------
# Martins, S.I.F.S. & van Boekel, M.A.J.S. (2005) "A kinetic model for the
# glucose/glycine Maillard reaction pathways", Food Chemistry 90(1-2):257-269,
# doi 10.1016/j.foodchem.2004.04.006 = Chapter 5 of Martins, S.I.F.S. (2003)
# PhD thesis, Wageningen University, ISBN 90-5808-823-5, Table 5.2 (thesis
# p. 122). Reparameterised Arrhenius, X = k at T_av = 100 C.
#
# SYSTEM: glucose 0.2 mol/L + glycine 0.2 mol/L in 0.1 mol/L phosphate buffer,
# pH 6.8 (initial, uncontrolled during heating), 10 mL in screw-capped Schott
# tubes, oil bath, 80/90/100/110/120 C.
#
# DECLARED ROLE: FIT. "Martins 2005 T2, all 10 steps + HPDs" is a FIT row of
# docs/reference/FIT_HOLDOUT_DECLARATION.md D.6 (Module 4). The BROWNING
# READOUT of step 9 -- its epsilon and the melanoidin response it converts --
# is the HOLD-OUT and is not in this file.
#
# TRANSCRIPTION PROVENANCE, stated because two transcriptions exist:
#   * steps 1 and 9 are taken from the K3 inventory, which re-read them from a
#     300 dpi render after finding that the PDF text layer STRIPS THE MINUS
#     SIGN FROM EVERY EXPONENT in this table (inventory lines 116-120). Those
#     two rows are the only ones the inventory prints with an uncertainty.
#   * steps 2-8 and 10 are taken from the repository's existing transcription
#     in scripts/generators/generate_trunk_rate_calibration.py (MARTINS_M4),
#     whose Ea are the printed ROUNDED values. Their k-side 95% HPDs are NOT
#     transcribed anywhere in the repository; the flag
#     'k_hpd_not_transcribed' says so on every affected row.
#   * the two transcriptions disagree slightly on step 1: the inventory gives
#     Ea 96.8 +/- 2.8 and X 1.6e-5, the repo table 97.0 +/- 3.0 and 1.61e-5.
#     The inventory's value is used, being the re-read one; the difference is
#     0.2 kJ/mol and is recorded rather than smoothed away.

_MARTINS_SOURCE = (
    "Martins & van Boekel 2005, Food Chem 90(1-2):257-269 Table 2 "
    "(= Martins 2003 thesis Table 5.2, p. 122), model M4"
)
_MARTINS_CONDITIONS = (
    "glucose 200 mmol/L + glycine 200 mmol/L, 0.1 mol/L phosphate, pH 6.8 "
    "initial (uncontrolled during heating), aqueous, screw-capped glass tubes, "
    "oil bath, 80-120 C"
)


def _martins(
    key: str,
    transformation: str,
    k_ref: float,
    ea: float,
    ea_ci: Optional[float],
    order: int,
    step: int,
    *,
    k_ci: Optional[float] = None,
    dossier: str = "",
    flags: Tuple[str, ...] = (),
    note: str = "",
) -> KineticParameter:
    return KineticParameter(
        key=key,
        transformation=transformation,
        k_ref=k_ref,
        ea_kj_mol=ea,
        unit="1/min" if order == 1 else "L/(mmol*min)",
        order=order,
        evidence_class="measured_rate",
        source_anchor=f"{_MARTINS_SOURCE}, step {step}",
        dossier_anchor=dossier or (
            "data/lit/extraction_dossiers/k3_final_parameter_inventory.md "
            "sec. A.1 line 123 ('Martins Table 2, all 10 steps')"
        ),
        conditions=_MARTINS_CONDITIONS,
        ph_of_measurement=6.8,
        temperature_range_c=(80.0, 120.0),
        rate_transfer="licensed_at_measurement_ph_only",
        k_ref_ci95=k_ci,
        ea_ci95_kj_mol=ea_ci,
        flags=flags,
        note=note,
    )


MARTINS_M4: Mapping[str, KineticParameter] = {
    "k_schiff": _martins(
        "k_schiff", "Glc + Gly -> Schiff base (condensation)",
        1.6e-5, 96.8, 2.8, 2, 1,
        dossier="k3_final_parameter_inventory.md sec. A.1 lines 116-117",
        flags=("composite_step_split_structurally", "minus_sign_stripped_in_pdf_text_layer"),
        note=(
            "Martins' step 1 is the ONE-STEP lumped condensation Glc + Gly -> DFG. "
            "The network writes it as two steps (condensation then rearrangement) "
            "because the Schiff base must exist as a state for any later module to "
            "attach to it -- but the SOURCE REFUSES THE SPLIT (see "
            "SCHIFF_AMADORI_SPLIT below), so this k is the composite's "
            "rate-determining constant and the composite is what is parameterised. "
            "The repo's shipped `schiff_condensation` Ea = 97.0 is CONFIRMED by this "
            "row; its shipped A = 1.5e11 L/(mol*s) is NOT (it over-predicts its own "
            "cited source by 14.8x; inventory line 118)."
        ),
    ),
    "k_glc_fru": _martins(
        "k_glc_fru", "Glc -> Fru (Lobry de Bruyn-Alberda van Ekenstein)",
        1.64e-3, 123.0, 5.0, 1, 2, flags=("k_hpd_not_transcribed", "ea_printed_rounded")),
    "k_fru_glc": _martins(
        "k_fru_glc", "Fru -> Glc (reverse isomerisation)",
        9.15e-3, 93.0, 3.0, 1, 3, flags=("k_hpd_not_transcribed", "ea_printed_rounded")),
    "k_ama_tdg": _martins(
        "k_ama_tdg", "Amadori (DFG) -> 3-deoxyglucosone + Gly",
        1.11e-2, 97.0, 2.0, 1, 4, flags=("k_hpd_not_transcribed", "ea_printed_rounded")),
    "k_tdg_fa": _martins(
        "k_tdg_fa", "3-deoxyglucosone -> formic acid (+ unmeasured C5 residue)",
        3.45e-2, 30.0, 9.0, 1, 5,
        flags=("k_hpd_not_transcribed", "ea_printed_rounded", "ea_conflicts_with_knol_2010"),
        note=(
            "Ea 30 +/- 9 kJ/mol is the lowest barrier in the Martins set and it "
            "CONFLICTS with Knol 2010 T2's formic-acid formation Ea of 84 +/- 14 "
            "kJ/mol (a FIT Module 4 row per the declaration's orchestrator decision "
            "2). The intervals do not overlap. Martins' value is operative here "
            "because it was measured on this system in this temperature window; the "
            "Knol value is carried in CROSS_LAB_COMPARATORS and the conflict is "
            "reported, not averaged."
        )),
    "k_ama_mgo": _martins(
        "k_ama_mgo", "Amadori (DFG) -> methylglyoxal + Gly (+ unmeasured C3 residue)",
        7.08e-3, 125.0, 5.0, 1, 6, flags=("k_hpd_not_transcribed", "ea_printed_rounded")),
    "k_ama_odg": _martins(
        "k_ama_odg", "Amadori (DFG) -> 1-deoxyglucosone + Gly",
        1.57e-2, 107.0, 7.0, 1, 7, flags=("k_hpd_not_transcribed", "ea_printed_rounded")),
    "k_odg_aa": _martins(
        "k_odg_aa", "1-deoxyglucosone -> acetic acid (+ unmeasured C4 residue)",
        1.45, 76.0, 4.0, 1, 8,
        flags=("k_hpd_not_transcribed", "ea_printed_rounded", "ea_agrees_with_knol_2010"),
        note="Knol 2010 T2's acetic-acid Ea, 75 +/- 10 kJ/mol, agrees within 1 kJ/mol.",
    ),
    "k_tdg_mel": _martins(
        "k_tdg_mel", "3-deoxyglucosone + Gly -> melanoidin (THE MASS SINK)",
        8.12e-4, 95.2, 2.3, 2, 9, k_ci=1.7e-5,
        dossier="k3_final_parameter_inventory.md sec. A.1 lines 119-120",
        flags=("mass_sink", "browning_readout_is_holdout"),
        note=(
            "THE FIRST BROWNING RATE PARAMETER THIS LANE HAS EVER HAD (inventory "
            "line 119, 'Z2 #15'). The RATE is a declared FIT row; the EPSILON that "
            "converts its product to an absorbance, and the browning response "
            "itself, are the declared HOLD-OUT and are not in this file. "
            "CIRCULARITY WARNING, stated once and carried into the hold-out report: "
            "Martins fitted this k to the very melanoidin response that is held out, "
            "so a hold-out evaluation that USES this k is a reproducibility check, "
            "not an out-of-sample test. The genuinely out-of-sample variant "
            "re-estimates the sink from the reactant side only."
        ),
    ),
    "k_fru_acids": _martins(
        "k_fru_acids", "Fru -> formic acid + acetic acid (+ unmeasured C3 residue)",
        4.41e-5, 237.0, 36.0, 1, 10,
        flags=("k_hpd_not_transcribed", "ea_printed_rounded", "ea_poorly_determined"),
        note="Ea 237 +/- 36 kJ/mol is the loosest in the set; the 80-120 C window "
             "is narrow for a barrier this large."),
}


# ---------------------------------------------------------------------------
# THE SCHIFF/AMADORI SPLIT -- a structural constant, not a measurement
# ---------------------------------------------------------------------------
# The build spec asks for Glc <-> Schiff <-> Amadori with the Amadori step 44.9x
# faster. The SOURCE REFUSES THE SPLIT: Martins tested the E1 (Schiff)
# intermediate three ways and then removing it entirely "fitted the data equally
# well"; "E1 is not a rate-determining step" (k3_final_parameter_inventory.md
# line 129, 'STRUCTURAL -- REFUSE THE SPLIT ... Carry ONE composite').
#
# THE RESOLUTION IMPLEMENTED HERE, stated plainly because it is a deviation:
# the Schiff base is carried as a STATE (later modules need something to attach
# to) but not as an independent PARAMETER. The rearrangement is pinned at
# 44.9x the condensation's pseudo-first-order rate at the experiment's own
# 200 mmol/L amine loading, which is the ratio the repository's own S3
# calibration recovered from these same FIT data
# (results/validation/trunk_rate_calibration_refit.md, "ratio 44.9"). Because
# the rearrangement is irreversible and has no competing sink, EVERY molecule
# that condenses eventually reaches the Amadori pool, so the pair carries the
# same TOTAL flux as Martins' one-step step 1.
#
# IT IS NOT FREE, AND THE SIZE OF THAT IS MEASURED RATHER THAN ASSERTED. At the
# pinned ratio the Schiff pool holds up about 2.5 mmol/L of amine at 100 C.
# Raising the ratio 10x (44.9 -> 449) moves the Amadori pool by -1.1% and the
# melanoidin sink by +7.5% at 100 C / 120 min, because releasing that amine
# faster feeds the BIMOLECULAR melanoidin step. Under 10% for a 10x change is
# small enough that no conclusion in this wave rests on the ratio -- it cannot
# produce the 1.45x browning bias the hold-out reports -- but it is not zero,
# and the unit test pins it rather than pretending otherwise.
SCHIFF_AMADORI_SPLIT: Dict[str, Any] = {
    "ratio_amadori_over_schiff_pseudo_first_order": 44.9,
    "amine_loading_mmol_L_for_the_ratio": 200.0,
    "evidence_class": "derived_from_fit_data",
    "source_anchor": (
        "results/validation/trunk_rate_calibration_refit.md (repository audit wave "
        "S3, profile likelihood over the Schiff/Amadori ordering), fitted to the "
        "SAME Martins FIT rows this module uses"
    ),
    "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 129",
    "source_verdict_on_the_split": (
        "REFUSED. Martins 2005 T1: removing the Schiff intermediate entirely "
        "'fitted the data equally well'; 'E1 is not a rate-determining step'."
    ),
    "consequence": (
        "the split is STRUCTURAL ONLY. It is not identifiable from any data in "
        "the corpus and no result in this module may be read as evidence for it."
    ),
    "reverse_steps": (
        "NOT IMPLEMENTED. The build spec writes the condensation as reversible "
        "(Glc <-> Schiff <-> Amadori). No reverse rate constant for either step "
        "exists in the corpus, so writing one would be an invented functional "
        "form. Glc <-> Fru IS reversible here because BOTH directions are "
        "measured (Martins steps 2 and 3). Declared structural gap."
    ),
}


# ---------------------------------------------------------------------------
# THE FITTED EXTENSION -- placeholders, populated by the fit
# ---------------------------------------------------------------------------
# These four steps have NO measured rate constant anywhere in the corpus. They
# exist because the network must give every intermediate a CONSUMPTION term:
# Martins' scheme lets methylglyoxal, formic acid and acetic acid accumulate
# forever, which is the exact defect ("nothing is ever consumed") this rebuild
# was commissioned to fix.
#
# Their values are estimated by least squares on the DECLARED FIT ROWS ONLY
# (the nine non-browning Martins responses) and written back as
# evidence_class='derived_from_fit_data'. Until that happens their k_ref and Ea
# are None and the network REFUSES TO INTEGRATE with them -- there is no
# silent default.

_FITTED_TEMPLATE: Dict[str, Tuple[str, int, str, str]] = {
    # key: (transformation, order, why it exists, where the carbon goes)
    "k_glc_frag": (
        "Glc -> unassigned fragment carbon (amine-independent sugar loss)", 1,
        "The amine-independent lane. Its EXISTENCE is measured -- Martins 2005 "
        "reports that removing the amine INVERTS the formic:acetic ratio "
        "(25% AA / 5% FA with glycine -> 1.2% FA / 0.65% AA without; inventory "
        "B2.9) -- but no rate constant for it exists. Bounded to allow ~zero so "
        "the data may reject it.",
        "FRAG_C",
    ),
    "k_mgo_mel": (
        "methylglyoxal -> melanoidin carbon", 1,
        "Methylglyoxal is the most reactive alpha-dicarbonyl in the trunk and "
        "Martins' scheme gives it NO sink at all. Routed into the melanoidin "
        "pool because alpha-dicarbonyl incorporation into the polymer is the "
        "only fate the corpus documents for it.",
        "MEL_C",
    ),
    "k_fa_frag": (
        "formic acid -> unassigned fragment carbon (decomposition)", 1,
        "Formic acid decomposition. Bounded to allow ~zero.",
        "FRAG_C",
    ),
    "k_aa_frag": (
        "acetic acid -> unassigned fragment carbon (decomposition)", 1,
        "Acetic acid decomposition. Bounded to allow ~zero.",
        "FRAG_C",
    ),
}

#: Search bounds for the fitted extension: (log10 k_ref at 100 C) and Ea kJ/mol.
#: Deliberately wide, and the fit is started from random points inside them, so
#: that any agreement with a literature value is a RESULT, not an initialisation.
FITTED_BOUNDS_LOG10K: Mapping[str, Tuple[float, float]] = {
    "k_glc_frag": (-8.0, -1.0),
    "k_mgo_mel": (-6.0, 1.0),
    "k_fa_frag": (-8.0, 0.0),
    "k_aa_frag": (-8.0, 0.0),
}
FITTED_EA_BOUNDS: Tuple[float, float] = (20.0, 260.0)

FITTED_KEYS: Tuple[str, ...] = tuple(_FITTED_TEMPLATE)


def fitted_placeholders() -> Dict[str, KineticParameter]:
    """The four fitted steps, unpopulated. Integration refuses them as-is."""
    out: Dict[str, KineticParameter] = {}
    for key, (transformation, order, why, sink) in _FITTED_TEMPLATE.items():
        out[key] = KineticParameter(
            key=key,
            transformation=transformation,
            k_ref=None,
            ea_kj_mol=None,
            unit="1/min" if order == 1 else "L/(mmol*min)",
            order=order,
            evidence_class="derived_from_fit_data",
            source_anchor=(
                "ESTIMATED IN THIS MODULE by least squares on the declared FIT rows "
                "of data/lit/timeseries/martins2005_glucose_glycine_80_100_120C_pH68.yml "
                "(the nine non-browning responses). No literature value exists."
            ),
            dossier_anchor=(
                "docs/reference/FIT_HOLDOUT_DECLARATION.md D.6 Module 4 "
                "('Martins 2005 T2, all 10 steps + HPDs' = FIT)"
            ),
            conditions=_MARTINS_CONDITIONS,
            ph_of_measurement=6.8,
            temperature_range_c=(80.0, 120.0),
            rate_transfer="not_licensed",
            flags=("fitted_here", "no_literature_value", f"carbon_sink:{sink}"),
            note=why,
        )
    return out


def with_fitted_values(
    fitted: Mapping[str, Tuple[float, float]]
) -> Dict[str, KineticParameter]:
    """Return the fitted steps populated with ``{key: (k_ref, Ea)}``."""
    from dataclasses import replace

    out = fitted_placeholders()
    for key, (k_ref, ea) in fitted.items():
        if key not in out:
            raise KeyError(f"{key!r} is not a fitted step; it is measured or unknown")
        out[key] = replace(out[key], k_ref=float(k_ref), ea_kj_mol=float(ea))
    return out


# ---------------------------------------------------------------------------
# CROSS-LAB COMPARATORS AND PRIORS -- carried, NOT operative
# ---------------------------------------------------------------------------
# These are measured numbers that the declaration assigns to Module 4 but that
# do NOT drive the ODE, either because a same-system value already does or
# because their transfer is not licensed. They are reported by the fit report
# so the conflicts are visible.
CROSS_LAB_COMPARATORS: Tuple[Dict[str, Any], ...] = (
    {
        "quantity": "Ea, sugar isomerisation",
        "value_kj_mol": 61.0,
        "ci95_kj_mol": 8.0,
        "source_anchor": "Knol et al. 2010, Table 2",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 125",
        "conditions": "aqueous asparagine/glucose, 120-200 C",
        "ph_of_measurement": None,
        "declared_role": "FIT (Module 4), per FIT_HOLDOUT_DECLARATION.md sec.5 decision 2",
        "operative": False,
        "why_not_operative": (
            "Martins measures the same transformation on the same system this "
            "module integrates, in the module's own 80-120 C window, and measures "
            "BOTH directions (123 +/- 5 forward, 93 +/- 3 reverse). Knol's 61 +/- 8 "
            "is a NET lumped isomerisation over 120-200 C. CONFLICT: 61 +/- 8 does "
            "not overlap either Martins direction. Reported, not averaged."
        ),
    },
    {
        "quantity": "Ea, acetic acid formation",
        "value_kj_mol": 75.0,
        "ci95_kj_mol": 10.0,
        "source_anchor": "Knol et al. 2010, Table 2",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 126",
        "conditions": "aqueous asparagine/glucose, 120-200 C",
        "ph_of_measurement": None,
        "declared_role": "FIT (Module 4)",
        "operative": False,
        "why_not_operative": (
            "AGREES with the operative Martins step 8 value (76 +/- 4) to within "
            "1 kJ/mol, across two labs, two amines and two temperature windows. "
            "This is the strongest cross-lab corroboration in Module 4 and it "
            "costs nothing, so it is reported as such rather than merged."
        ),
    },
    {
        "quantity": "Ea, formic acid formation",
        "value_kj_mol": 84.0,
        "ci95_kj_mol": 14.0,
        "source_anchor": "Knol et al. 2010, Table 2",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 127",
        "conditions": "aqueous asparagine/glucose, 120-200 C",
        "ph_of_measurement": None,
        "declared_role": "FIT (Module 4)",
        "operative": False,
        "why_not_operative": (
            "DIRECT CONFLICT with the operative Martins step 5 value (30 +/- 9). "
            "The intervals are 54 kJ/mol apart and do not overlap. Martins is "
            "operative because it is the same system in the same window; the "
            "conflict is a standing owner question, not something this module "
            "resolves."
        ),
    },
    {
        "quantity": "Ea, condensation (Glc + Asn -> Schiff base)",
        "value_kj_mol": 57.6,
        "ci95_kj_mol": 8.0,
        "source_anchor": "Knol et al. 2005, Table 1",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 124",
        "conditions": "aqueous asparagine/glucose, 120-200 C",
        "ph_of_measurement": None,
        "declared_role": "FIT (Module 3 trunk); listed here because it collides with Module 4",
        "operative": False,
        "why_not_operative": (
            "CONFLICT flagged by the inventory itself (Z2 #17): 57.6 +/- 8.0 versus "
            "Martins' 96.8 +/- 2.8 for the same condensation, intervals "
            "non-overlapping. Different amine (asparagine vs glycine) and a "
            "different window. Owner question before either is called 'the' "
            "condensation barrier."
        ),
    },
    {
        "quantity": "Ea, browning A420 (glucose + casein / + BSA, alkaline)",
        "value_kj_mol": None,
        "value_set_kj_mol": [164.0, 120.0, 130.0, 126.0, 92.0, 95.0],
        "source_anchor": "Ajandouz et al. 2008, Table 3 p. 1250",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1.1 lines 147-148, 168-170",
        "conditions": "glucose 0.2 M + casein or BSA 5 mg/mL, pH 8.0 or 9.7, 60-100 C",
        "ph_of_measurement": 8.0,
        "declared_role": "Module 9, FIT as priors on Ea only",
        "operative": False,
        "rate_transfer": "not_licensed",
        "why_not_operative": (
            "THE ALKALINE-pH WALL. Ajandouz licenses Ea transfer to pH 5-7 only for "
            "glucose-loss and amino-loss, and EXPLICITLY NOT FOR BROWNING, whose Ea "
            "fall 15-29% between pH 8.0 and 9.7 (inventory sec. C.13). Carried as an "
            "unvalidated prior with its measurement pH attached, exactly as "
            "FIT_HOLDOUT_DECLARATION.md sec.5 decision 3 requires. It is a "
            "third-lab CONTEXT for Martins' melanoidin Ea 95.2, not a referee for it."
        ),
    },
    {
        "quantity": "Ea, beta-elimination (cysteine)",
        "value_kj_mol": 123.0,
        "ci95_kj_mol": None,
        "source_anchor": "Zheng & Ho 1994, Tables I/III/V",
        "dossier_anchor": "k3_final_parameter_inventory.md sec. A.1 line 130",
        "conditions": "aqueous, pH 9",
        "ph_of_measurement": 9.0,
        "declared_role": "FIT (Module 4), 'kinetic reference ... not a benchmark source'",
        "operative": False,
        "rate_transfer": "not_licensed",
        "why_not_operative": (
            "Cysteine is not in the glucose/glycine trunk; this constant belongs to "
            "the sulfur module that plugs into this core. It is registered here so "
            "that module inherits it WITH its pH 9 label rather than re-deriving it. "
            "It also condemns the repo's shipped beta_elimination_dha Ea of 79.5 "
            "(43.5 kJ/mol below the measured aqueous barrier, with a silent NaN "
            "prefactor). "
            "GAP REPORTED: the declaration assigns 'Zheng 1994 Tables I/III/V "
            "(36 k + 8 Ea)' to FIT, but this single value is the ONLY Zheng number "
            "transcribed anywhere in the repository's dossiers. The other 43 are "
            "not available to be used."
        ),
    },
)


# ---------------------------------------------------------------------------
# Policy checks
# ---------------------------------------------------------------------------


def assert_no_dft(registry: Mapping[str, KineticParameter] = MARTINS_M4) -> None:
    """
    Owner policy, pinned: no DFT-derived barrier may enter this module.

    The dataclass already refuses an unknown evidence_class; this is the
    belt-and-braces check the tests call, and it also catches a barrier that
    was smuggled in under a measured label with a computational anchor.
    """
    banned = ("dft", "b3lyp", "wb97", "m06", "cbs-qb3", "g4", "computed", "calculated")
    for key, param in registry.items():
        haystack = " ".join(
            [param.source_anchor, param.dossier_anchor, param.note, param.evidence_class]
        ).lower()
        for token in banned:
            if token in haystack:
                raise ValueError(
                    f"{key}: parameter provenance mentions {token!r}. DFT-derived "
                    f"barriers are refused by owner policy."
                )


def check_ph_homogeneity(
    registry: Mapping[str, KineticParameter]
) -> float:
    """
    Refuse to assemble a network from parameters measured at different pH.

    There is no pH term in this module (see the module docstring, policy 2), so
    mixing pH values would silently assert a transferability the corpus
    explicitly denies. Returns the single pH the network is valid at.
    """
    values = {
        p.ph_of_measurement for p in registry.values() if p.ph_of_measurement is not None
    }
    if not values:
        raise ValueError("no parameter declares a pH of measurement")
    if len(values) > 1:
        raise ValueError(
            f"parameters span several pH values {sorted(values)} and this module "
            f"has NO pH term. Mixing them would invent a transfer the corpus "
            f"refuses (k3_final_parameter_inventory.md sec. B.2: 'NO family-level "
            f"pH term can pass')."
        )
    return float(next(iter(values)))


#: The pH the whole operative network is valid at, checked at import.
NETWORK_PH: float = check_ph_homogeneity(MARTINS_M4)

assert_no_dft(MARTINS_M4)


def registry_metadata(
    parameters: Mapping[str, KineticParameter]
) -> Dict[str, Any]:
    """Full runtime metadata block for a parameter set."""
    return {
        "reference_temperature_K": T_REF_K,
        "network_ph": NETWORK_PH,
        "aw_of_measurement": AW_OF_MEASUREMENT,
        "ph_term_present": False,
        "aw_term_present": False,
        "ph_term_absent_because": (
            "no licensed pH dependency exists for the trunk. "
            "k3_final_parameter_inventory.md sec. B.2 records five independent "
            "sign-crossings and concludes 'NO family-level pH term can pass'; "
            "Ajandouz 2008 explicitly refuses to transfer browning Ea off its own "
            "pH (sec. C.13). The network is pH-FIXED at the value every operative "
            "parameter was measured at."
        ),
        "aw_term_absent_because": (
            "nothing in the Module 4 fit corpus varies water activity; the whole "
            "Martins/Knol trunk is dilute aqueous. An a_w term would be invented."
        ),
        "dft_free": True,
        "parameters": {k: p.as_metadata() for k, p in sorted(parameters.items())},
        "schiff_amadori_split": dict(SCHIFF_AMADORI_SPLIT),
        "cross_lab_comparators": [dict(row) for row in CROSS_LAB_COMPARATORS],
    }
