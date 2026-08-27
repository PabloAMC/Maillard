from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from src.barrier_constants import arrhenius_rate_constant, get_donor_reactivity_multiplier
from src.extrusion import build_extrusion_process_profile, normalize_moisture_regime


# ── Amine-nucleophile matcher exemption (2026-08-27, Wave H) ──────────────────
# `_ionization_correction` and `_water_activity_correction` classify a reaction
# family by SUBSTRING.  Both were written for the amine-nucleophile steps of the
# cascade — Schiff condensation, the Amadori/Heyns rearrangement, Strecker — where
# the rate really does track the fraction of free (deprotonated) amine and where
# Labuza's aw optimum was measured.
#
# Wave G1 introduced `Enolisation_2_3_Amadori`: the 2,3-enolisation /
# beta-elimination that RELEASES the amino acid from the Amadori compound and
# opens the accepted 1-deoxyosone -> norfuraneol -> MFT route (van den Ouweland &
# Peer 1975).  It has no amine nucleophile — the amine is a leaving group — but
# its NAME contains "amadori", so it silently collected both corrections.  At the
# Hofmann1998 conditions (pH 5.0, aw 0.98) that is 1e-3 * 0.62 = 6.2e-4, i.e.
# +6.06 kcal/mol of effective barrier, which pushed the real MFT route ~1600x
# below the demoted one-step shortcut that Wave G1 had just replaced it with: the
# flagship chemistry fix was inert, and the two knobs Wave H was asked to refit
# (`thiol_addition_norfuraneol`, `furanone_cyclisation`) had exactly zero
# derivative on every prediction.
#
# The exemption is deliberately narrow: it fires only on families that are
# themselves enolisations / eliminations / deaminations, i.e. steps in which the
# nitrogen leaves.  Of the ~35 families the engine emits, `Enolisation_2_3_Amadori`
# is the only one that both matches an amine-nucleophile keyword and is such a
# step, so this changes the classification of exactly one family.  Measured
# prediction impact at the shipped barriers: none — the norfuraneol route becomes
# the selected MFT path with the identical span and flux the shortcut had, so the
# panel is byte-identical (verified 2026-08-27).  What it restores is
# identifiability of the new route's barriers.
_AMINE_LEAVING_GROUP_FAMILY_MARKERS = ("enolisation", "elimination", "deamination")


def _releases_rather_than_attacks_with_the_amine(normalized_family: str) -> bool:
    """True for families in which the amine is a leaving group, not a nucleophile."""
    return any(marker in normalized_family for marker in _AMINE_LEAVING_GROUP_FAMILY_MARKERS)


# ── Wave S1b (2026-08-27): the pH / water-activity routing repair ─────────────
#
# WHY THIS BLOCK EXISTS.  Wave S2 built a 52-claim directional-accuracy panel
# (`docs/validation/directional_claims_panel.yml`) — the first measurement in this
# repo of the ORDINAL question the model is actually used for ("which way does X
# move when I turn this knob").  The model scored 16/19 on sugar identity,
# temperature, time and precursor loading, and 2/10 on pH and water activity —
# worse than a coin.  The failure was not a modelling failure.  It was three
# ROUTING defects, each verified by inspection before it was touched:
#
#   1. `ReactionConditions.get_ph_multiplier` — which encodes the documented
#      enolisation route selection and is unit-tested for the correct signs —
#      was NEVER CALLED on the prediction path.  `benchmark_validation`
#      line ~662 goes through `get_rate_constant()`, which applied
#      `_ionization_correction` and `_water_activity_correction` and not
#      `get_ph_multiplier`.  Confirmed by grep: the only callers were
#      `kinetics.py`, `pathway_ranker.py` and `cantera_export.py`, none of which
#      is reachable from `evaluate_benchmark_payload`.  Written, tested, never run.
#   2. The pyrazine branch of `_ionization_correction` keyed on the SUBSTRING
#      "pyrazine".  MEASURED: of the 29 reaction families this engine actually
#      emits across `data/benchmarks/`, NOT ONE contains that substring.  The
#      pyrazine-forming step here is `Aminoketone_Condensation` (two aminoacetones
#      -> 2,5-dimethylpyrazine + 2 H2O + H2, `reaction_templates.py:1000`).  The
#      branch returned 1.0 at every pH and the model moved dimethylpyrazine the
#      WRONG WAY with pH against two independent direct measurements.
#   3. `_water_activity_correction` had the right Labuza shape and reached 3 of
#      those 29 families, none of them on the furan/HMF track, and its "furfural"
#      branch was DEAD for the same substring reason (no emitted family contains
#      "furfural" either).
#
# THE ROOT CAUSE IS THE SUBSTRING MATCHING ITSELF, and it is a repeat offender in
# this repo: `Furanone_Strecker_Reduction` (Wave P ledger) and the Wave I
# offset-key fix are the same bug wearing different hats.  So the sets below are
# EXPLICIT FAMILY NAMES.  Every name in them was verified to be emitted by the
# engine; every name NOT in them is justified in the comment above its set.
#
# WHICH pH TERM APPLIES WHERE — the two functions are DIFFERENT PHYSICS, and the
# reason this had to be reasoned out rather than just wired is that naive wiring
# would have applied both to the same family:
#
#   `_ionization_correction` is a REAGENT-AVAILABILITY term: the Henderson-
#   Hasselbalch fraction of the nucleophile that is present as the free base and
#   can therefore attack at all.  It is monotone increasing in pH and it belongs
#   to the families in which a nitrogen lone pair does the attacking.
#
#   `get_ph_multiplier` is a ROUTE-SELECTION term: which way the Amadori compound
#   enolises.  1,2-enolisation is acid-catalysed and opens the 3-deoxyosone ->
#   furfural / HMF track; 2,3-enolisation is base-catalysed and opens the
#   1-deoxyosone -> reductone -> pyrazine / furanone track.  That is a branch-point
#   partition, not a reagent count, and it is the pH physics the amine term cannot
#   express.
#
# So the two sets are DISJOINT BY CONSTRUCTION and the disjointness is asserted at
# import time below — each family gets the pH physics exactly once.
#
# `get_ph_multiplier` is deliberately NOT wired to the families its own substrings
# would have swept up ("thiol", "thio", "cysteine", "furan", "condensation").
# Those matches are substring accidents of the kind this block exists to remove:
# the acid preference of the furan track is a property of the ENOLISATION BRANCH
# POINT, and applying a 4.9x acid boost again at every downstream thiol addition,
# cyclodehydration and oxidation step would compound one physical effect five or
# six times along a single route.  Worse, "condensation" would have claimed
# `Aminoketone_Condensation` — the pyrazine step — for the acid-peaked Schiff
# Gaussian, i.e. defect 2 all over again in the opposite direction.

#: Families whose rate tracks the free (deprotonated) alpha-amino nucleophile.
#: pKa 8.0.  This set reproduced EXACTLY what the old "amadori"/"strecker"/"schiff"
#: substrings reached among emitted families — it was a de-substring-ing, not a
#: widening.  `Enolisation_2_3_Amadori` is excluded for the Wave H reason above
#: (the amine LEAVES; it does not attack).
#:
#: WAVE T4 (2026-08-27) — ONE DELIBERATE WIDENING, and it is labelled as such
#: because S1b's "not a widening" claim above no longer holds verbatim.
#: `heyns_rearrangement` is added.  It is the KETOSE analogue of the Amadori
#: rearrangement (`barrier_constants.py:146`, "Ketose analogue of Amadori"),
#: emitted by `reaction_templates.py:60` for fructose + an amino acid, carrying
#: the same alpha-amino nitrogen through the same Schiff-base -> 1,2-shift
#: rearrangement.  The old `"amadori"` substring did not match the string
#: `heyns_rearrangement`, so this gap PREDATES Wave S1b and S1b did not create
#: it — what S1b created was the false claim at `get_ph_multiplier`'s docstring
#: that no Heyns family is emitted.  The gap's measured size: at pH 5 the
#: rearrangement of a ketose was running at 1208x the corrected rate of the
#: physically identical aldose step, because it escaped BOTH this term and the
#: Labuza term below.  Two chemically identical-in-kind steps receiving
#: corrections that differ by three orders of magnitude is not a defensible
#: state, whichever of the two is right.
#:
#: HONEST CAVEAT, so the next reader does not mistake this for a derivation:
#: including `amadori_rearrangement` here is itself inherited from the legacy
#: substring, not re-derived from mechanism — for a rearrangement the amine is
#: ALREADY bound and it is the 1,2-proton shift that is rate-limiting, so the
#: free-base fraction is a proxy for the upstream condensation equilibrium
#: rather than for a nucleophilic attack.  This wave makes Heyns CONSISTENT with
#: Amadori; it does not settle whether that shared treatment is correct.  Doing
#: so means re-deriving the aldose case, which is a calibration decision with
#: owner sign-off, and is filed as [P] in `tasks/audit_remediation.md`.
_ALPHA_AMINO_NUCLEOPHILE_FAMILIES = frozenset({
    "schiff_base_formation",
    "lipid_schiff_base",
    "amadori_rearrangement",
    "heyns_rearrangement",   # Wave T4: the ketose twin of the line above
    "strecker_degradation",
    "lipid_strecker_synergy",
})

#: The pyrazine branch, given the family it was always meant to name.  Two
#: alpha-aminoketones self-condense; the attacking species is the free amine of
#: one of them, and an alpha-carbonyl is strongly electron-withdrawing, so that
#: amine is far more acidic than an alpha-amino acid's.  The pKa 6.5 is the value
#: the dead branch already carried and is NOT retuned here.
#: `Lipid_Thiazole_Condensation` is EXCLUDED: its nitrogen donor is ammonia
#: (pKa 9.25), a third pKa this function does not carry — left uncorrected rather
#: than corrected with the wrong constant.
_AMINOKETONE_NUCLEOPHILE_FAMILIES = frozenset({
    "aminoketone_condensation",
})

#: The enolisation branch point, and ONLY the branch point.  These are the three
#: families that ARE the 1,2- / 2,3-enolisation partition:
#:   `Enolisation_1_2`         3-deoxyosone -> furfural/HMF + 2 H2O  (acid arm)
#:   `Enolisation_2_3`         3-deoxyosone -> pyruvaldehyde + C2/C3 (base arm)
#:   `Enolisation_2_3_Amadori` Amadori -> 1-deoxyosone + amino acid  (base arm)
#: `Enolisation_Intermediate` (Amadori -> 3-deoxyosone) is EXCLUDED: it is the
#: step BEFORE the branch point, common to both arms, so a pH factor there would
#: scale both arms identically and express no selection.
#: `Fructofuranosyl_Dehydration` is EXCLUDED from the pH term: it is a direct
#: ketose dehydration that bypasses the Amadori branch point entirely, so the
#: 1,2-vs-2,3 partition does not describe it.
_ENOLISATION_ROUTE_PH_FAMILIES = frozenset({
    "enolisation_1_2",
    "enolisation_2_3",
    "enolisation_2_3_amadori",
})

assert not (_ENOLISATION_ROUTE_PH_FAMILIES & _ALPHA_AMINO_NUCLEOPHILE_FAMILIES), (
    "a family would receive the pH physics twice"
)
assert not (_ENOLISATION_ROUTE_PH_FAMILIES & _AMINOKETONE_NUCLEOPHILE_FAMILIES), (
    "a family would receive the pH physics twice"
)
assert not (_ALPHA_AMINO_NUCLEOPHILE_FAMILIES & _AMINOKETONE_NUCLEOPHILE_FAMILIES), (
    "a family would receive the amine-ionisation term twice"
)


# ── Water activity: which families the Labuza physics governs ─────────────────
#
# The criterion is MEASURED STOICHIOMETRY, not a keyword and not a judgement call.
# `<scratchpad>/s1b_water.py` counts net water produced per step for every family
# the engine emits over all of `data/benchmarks/`.  Every family is stoichio-
# metrically UNIFORM (one net water count across all its steps) except
# `Additive_Thermal_Degradation`, which is +2/0/-1/-2 and is therefore excluded:
# a single family-level factor cannot honestly represent it.
#
# The physics is mass action on the step's own water: a step that RELEASES water
# is pushed back by water (Le Chatelier — and for the imine/enol condensations
# here the reverse, hydrolysis, is genuinely fast in water); a step that CONSUMES
# water is rate-limited by water availability.  Neither arm is fitted.
#
# LIMITATION, stated rather than hidden: the inhibition factor is applied ONCE per
# step regardless of whether the step sheds one water or three, because the shipped
# `1.3 - aw` shape carries no stoichiometric exponent and adding one would be
# retuning the curve.  A three-water dehydration is therefore under-penalised
# relative to a one-water one.
#
# ALSO A LIMITATION: for families that lump a redox into the same step
# (`Deoxyosone_Reduction`, `Mercaptoketone_Formation`, `Thiol_Addition_H2`,
# `Furanone_Reductive_Opening`, `Furan_Ring_Aromatisation`, and
# `Aminoketone_Condensation`'s H2), the released water is partly a redox byproduct
# rather than purely a condensation equilibrium, so the reverse reaction the
# Le Chatelier argument invokes is not fully accessible.  They are included anyway,
# because excluding them would require exactly the kind of per-family judgement
# call this block replaced with a measurement — but the reader should know the
# argument is weaker for them.

#: Net water RELEASED (>= +1 per step, uniform): dehydrations and condensations,
#: inhibited by water.  Measured counts in parentheses.
_WATER_RELEASING_FAMILIES = frozenset({
    "aminoketone_condensation",              # +2  2 aminoacetone -> DMP + 2 H2O + H2
    "enolisation_1_2",                       # +2  3-deoxyosone -> furfural/HMF   <- THE FURAN TRACK
    "fructofuranosyl_dehydration",           # +3  ketose -> HMF direct
    "deoxyosone_reduction",                  # +1
    "furan_ring_aromatisation",              # +1
    "furanone_amino_acid_reduction",         # +1
    "furanone_cyclisation",                  # +1  1-deoxyosone -> norfuraneol + H2O
    "furanone_formation",                    # +1  (sibling of the above)
    "furanone_reductive_opening",            # +1
    "lipid_thiazole_condensation",           # +3
    "mercaptoketone_cyclodehydration",       # +2  aldol -> MFT + 2 H2O
    "mercaptoketone_formation",              # +1
    "safety_risk_acrylamide",                # +1
    "schiff_base_formation",                 # +1  the canonical condensation
    "lipid_schiff_base",                     # +1  same chemistry, SMIRKS tier
    "thiol_addition_h2",                     # +1
    "thiol_addition_hexose_legacy_shortcut",  # +3
    "thiol_addition_pentodiulose",           # +2
    "thiol_addition_norfuraneol",            # +1  (emitted only on the norfuraneol route)
    "thiol_dehydration",                     # +1  thiohemiacetal -> FFT + H2O
})

#: Net water CONSUMED (-1 per step, uniform): hydrolytic steps, rate-limited by
#: water availability.  Mass action, factor = aw, floored the same way the
#: dehydration arm is floored.
_WATER_CONSUMING_FAMILIES = frozenset({
    "cysteine_degradation",                  # -1  hydrolytic release of H2S
    "generalized_deamination",               # -1  hydrolytic deamination
})

#: Net ZERO water, and therefore NO mass-action term: `Amadori_Rearrangement`,
#: `Heyns_Rearrangement`, `Strecker_Degradation`, `Beta_Elimination`,
#: `Enolisation_2_3`, `Enolisation_2_3_Amadori`, `Enolisation_Intermediate`,
#: `Mercaptoketone_Aldol_Addition`, `Retro_Aldol_Fragmentation`,
#: `Thiohemiacetal_Formation`, `Thiol_Oxidation`.  The first three nevertheless
#: keep the EMPIRICAL Labuza peaked curve: that curve is not a mass-action term at
#: all, it is the measured overall-browning response (diffusion-limited below
#: aw ~0.4, diluted above ~0.8) and it was measured on exactly the
#: amine-condensation chemistry those families carry.  The curve is untouched.
#:
#: WAVE T4 (2026-08-27): `heyns_rearrangement` added, on the same MEASURED
#: criterion the block header states.  Counted off the actual emitted step, not
#: assumed — `SmirksEngine.enumerate(resolve_many(["D-Fructose", "Glycine"]))`
#: emits exactly one Heyns step:
#:     D-Fructose-Glycine-Schiff-base  ->  D-Fructose-Glycine-Heyns
#:     OCC(=NCC(=O)O)C(O)C(O)C(O)CO    ->  O=CC(NCC(=O)O)C(O)C(O)C(O)CO
#: One reactant, one product, no water on either side: net water 0, uniform over
#: its single step — byte-identical stoichiometry to `Amadori_Rearrangement`.  It
#: therefore takes NO mass-action term (it is in neither arm below) and takes the
#: empirical browning curve for the same reason Amadori does.
_LABUZA_EMPIRICAL_FAMILIES = frozenset({
    "amadori_rearrangement",
    "heyns_rearrangement",   # Wave T4: net water 0, measured; ketose twin of the above
    "strecker_degradation",
    "lipid_strecker_synergy",
})

assert not (_WATER_RELEASING_FAMILIES & _WATER_CONSUMING_FAMILIES), (
    "a family cannot both release and consume water"
)

# WAVE T4 (2026-08-27): the third leg of the same discipline, previously
# unasserted.  `_water_activity_correction` tests the empirical set FIRST and
# returns from whichever arm matches, so an overlap would not double-multiply —
# it would SILENTLY SHADOW, giving a water-shedding family the peaked browning
# curve instead of its mass-action term with no error anywhere.  That is the
# quieter half of the failure class this module keeps catching, so it is now an
# import-time error rather than a reading exercise.  Each family gets EXACTLY ONE
# water-activity treatment: empirical, releasing, consuming, or none.
assert not (_LABUZA_EMPIRICAL_FAMILIES & _WATER_RELEASING_FAMILIES), (
    "a family would take both the empirical browning curve and the dehydration "
    "mass-action term; the empirical arm would silently shadow the other"
)
assert not (_LABUZA_EMPIRICAL_FAMILIES & _WATER_CONSUMING_FAMILIES), (
    "a family would take both the empirical browning curve and the hydrolysis "
    "mass-action term; the empirical arm would silently shadow the other"
)


@dataclass
class ReactionConditions:
    """
    Environmental parameters governing the Maillard cascade.
    """
    def __init__(self,
                 pH: float = 6.0,
                 temperature_celsius: float = 120.0,
                 water_activity: float = 0.8,
                 fat_fraction: float = 0.0,
                 protein_fraction: float = 1.0,
                 dielectric_constant: float = 78.4, # Default: Water
                 solvent_name: str = "water",
                 matrix_fiber: float = 0.0, # Placeholder for blind spot
                 metal_catalyst: Optional[str] = None, # Placeholder for blind spot
                 protein_type: str = "free",
                 sme_kj_per_kg: float = 0.0,
                 moisture_regime: Optional[str] = None,
                 screw_speed_rpm: Optional[float] = None,
                 feed_rate_kg_per_h: Optional[float] = None,
                 die_exit_temperature_celsius: Optional[float] = None,
                 sterilization_temperature_celsius: Optional[float] = None,
                 sterilization_time_minutes: float = 0.0,
                 barrel_zone_temperatures: Optional[Sequence[float]] = None,
                 barrel_zone_time_fractions: Optional[Sequence[float]] = None,
                 temperature_profile: Optional[Sequence[Sequence[float]]] = None,
                 water_activity_profile: Optional[Sequence[Sequence[float]]] = None,
                 ):
        self.pH = pH
        self.temperature_celsius = temperature_celsius
        self.water_activity = water_activity
        self.fat_fraction = fat_fraction
        self.protein_fraction = protein_fraction
        self.dielectric_constant = dielectric_constant
        self.solvent_name = solvent_name
        self.matrix_fiber = matrix_fiber
        self.metal_catalyst = metal_catalyst
        self.protein_type = protein_type
        self.sme_kj_per_kg = float(sme_kj_per_kg)
        self.moisture_regime = moisture_regime
        self.screw_speed_rpm = None if screw_speed_rpm is None else float(screw_speed_rpm)
        self.feed_rate_kg_per_h = None if feed_rate_kg_per_h is None else float(feed_rate_kg_per_h)
        self.die_exit_temperature_celsius = None if die_exit_temperature_celsius is None else float(die_exit_temperature_celsius)
        self.sterilization_temperature_celsius = sterilization_temperature_celsius
        self.sterilization_time_minutes = float(sterilization_time_minutes)
        self.barrel_zone_temperatures = None if barrel_zone_temperatures is None else [float(value) for value in barrel_zone_temperatures]
        self.barrel_zone_time_fractions = None if barrel_zone_time_fractions is None else [float(value) for value in barrel_zone_time_fractions]
        self.temperature_profile = None if temperature_profile is None else [
            (float(pair[0]), float(pair[1])) for pair in temperature_profile
        ]
        self.water_activity_profile = None if water_activity_profile is None else [
            (float(pair[0]), float(pair[1])) for pair in water_activity_profile
        ]
        self.__post_init__()

    def __post_init__(self):
        """Set dielectric constant based on solvent name if provided."""
        presets = {
            "water": 78.4,
            "ethanol": 24.5,
            "methanol": 32.7,
            "lipid": 2.0,
            "benzene": 2.3,
            "dimethyl_sulfoxide": 46.7
        }
        if self.solvent_name.lower() in presets:
            self.dielectric_constant = presets[self.solvent_name.lower()]
        self.moisture_regime = normalize_moisture_regime(self.moisture_regime, self.water_activity)
        self._validate_profile(self.temperature_profile, "temperature_profile")
        self._validate_profile(self.water_activity_profile, "water_activity_profile")

    @staticmethod
    def _validate_profile(profile: Optional[Sequence[Sequence[float]]], profile_name: str) -> None:
        if profile is None:
            return
        if len(profile) < 2:
            raise ValueError(f"{profile_name} must include at least an initial point and a terminal point")
        last_time = None
        for pair in profile:
            time_point = float(pair[0])
            if last_time is not None and time_point <= last_time:
                raise ValueError(f"{profile_name} time points must be strictly increasing")
            last_time = time_point

    @property
    def extrusion_profile(self) -> dict[str, object]:
        return build_extrusion_process_profile(
            base_temperature_celsius=self.temperature_celsius,
            water_activity=self.water_activity,
            protein_type=self.protein_type,
            sme_kj_per_kg=self.sme_kj_per_kg,
            moisture_regime=self.moisture_regime,
            screw_speed_rpm=self.screw_speed_rpm,
            feed_rate_kg_per_h=self.feed_rate_kg_per_h,
            die_exit_temperature_celsius=self.die_exit_temperature_celsius,
            sterilization_temperature_celsius=self.sterilization_temperature_celsius,
            sterilization_time_minutes=self.sterilization_time_minutes,
            zone_temperatures=self.barrel_zone_temperatures,
            zone_time_fractions=self.barrel_zone_time_fractions,
        )

    @property
    def effective_temperature_celsius(self) -> float:
        return float(self.extrusion_profile.get("effective_temperature_celsius", self.temperature_celsius))

    @property
    def temperature_kelvin(self) -> float:
        return self.effective_temperature_celsius + 273.15

    @property
    def nominal_temperature_kelvin(self) -> float:
        return self.temperature_celsius + 273.15

    @property
    def has_dynamic_profile(self) -> bool:
        return bool(self.temperature_profile and len(self.temperature_profile) >= 2)

    def temperature_celsius_at_time(self, time_minutes: float) -> float:
        if not self.has_dynamic_profile:
            return float(self.temperature_celsius)
        profile = sorted(self.temperature_profile, key=lambda item: item[0])
        time_points = [point[0] for point in profile]
        temp_points = [point[1] for point in profile]
        return float(np.interp(float(time_minutes), time_points, temp_points))

    def temperature_celsius_at(self, time_minutes: float) -> float:
        return self.temperature_celsius_at_time(time_minutes)

    def water_activity_at(self, time_minutes: float) -> float:
        if not self.water_activity_profile or len(self.water_activity_profile) < 2:
            return float(self.water_activity)
        profile = sorted(self.water_activity_profile, key=lambda item: item[0])
        time_points = [point[0] for point in profile]
        aw_points = [point[1] for point in profile]
        return float(np.interp(float(time_minutes), time_points, aw_points))
        
    def _sigmoid(self, x: float, center: float, k: float) -> float:
        """Helper for sigmoid transitions."""
        import math
        try:
            return 1.0 / (1.0 + math.exp(-k * (x - center)))
        except OverflowError:
            return 0.0 if x < center else 1.0

    def _gaussian(self, x: float, center: float, sigma: float) -> float:
        """Helper for peaked responses (e.g. Schiff base)."""
        import math
        return math.exp(-0.5 * ((x - center) / sigma) ** 2)

    def get_ph_multiplier(self, reaction_family: str) -> float:
        """
        Calculates a kinetic multiplier based on pH using smooth sigmoid/Gaussian curves.
        
        Physics:
        - 1,2-enolisation (furans): Favored at acidic pH (pH < 6).
        - 2,3-enolisation (pyrazines): Favored at alkaline pH (pH > 6).
        - Schiff base: Optimal at pH 5.5 (amine nucleophilicity vs protonation balance).

        WAVE S1b (2026-08-27).  This function is now reachable from the prediction
        path, but ONLY through `_enolisation_route_ph_correction`, which admits the
        three enolisation branch-point families and nothing else.  The shapes and
        constants are untouched.  Its OTHER callers (`kinetics.py`,
        `pathway_ranker.py`, `cantera_export.py`) still call it directly with the
        raw substring matching below, and those paths do NOT feed
        `evaluate_benchmark_payload`.

        DEAD SUBSTRING KEYS, reported rather than silently deleted because the
        keys document intent:
          * "pyrazine"            - the pyrazine step is `Aminoketone_Condensation`
          * "nitrogen_heterocycle", "oxygen_heterocycle" - never used as names
          * "1,2" and "2,3"       - families spell these `1_2` / `2_3`; the
                                    comma forms match nothing (the underscore
                                    forms, also present, do the work).  NOT
                                    removable: `tests/unit/test_conditions.py`
                                    calls this function with the FAST_BARRIERS
                                    canonical spellings ("1,2-enolisation"),
                                    which the comma forms are what reach.
        "heyns" is LIVE AND LOAD-BEARING — DO NOT DELETE IT.

        WAVE T4 (2026-08-27) CORRECTION OF THIS DOCSTRING.  It used to list
        "heyns" among the dead keys with the gloss "no Heyns family is emitted".
        That was FALSE, and it was the most dangerous line in the file, because
        "heyns" is the ONLY key in branch 2 that matches `Heyns_Rearrangement`
        ("amadori" does not appear in the string "heyns_rearrangement"), so a
        cleanup agent acting on the sentence would have silently deleted the pH
        dependence of the ketose rearrangement from the three ungated lanes.
        The engine emits `Heyns_Rearrangement` for any ketose + amino acid
        (`reaction_templates.py:60`; `tests/unit/test_smirks_engine.py::
        TestHeynsRearrangement::test_heyns_fires` has asserted it fires since
        before Wave S1b).  The claim survived because the S1b census enumerated
        families over `data/benchmarks/`, and NO SHIPPED BENCHMARK USES A KETOSE
        — a true statement about the panel, written down as a statement about
        the engine.  The census enumerator now runs a fructose pool as well.

        The Schiff branch's "condensation" key is LIVE and DANGEROUS: it matches
        `Aminoketone_Condensation`, i.e. it would give the pyrazine step an
        acid-peaked Gaussian.  That is precisely why the gate exists.

        NO DOUBLE-APPLICATION, and it is worth stating because Wave T4 also added
        `heyns_rearrangement` to `_ALPHA_AMINO_NUCLEOPHILE_FAMILIES`.  Heyns is
        now treated EXACTLY as Amadori already was, in both lanes:
          * prediction lane (`get_rate_constant`) — `get_ph_multiplier` is reached
            only through `_enolisation_route_ph_correction`, which admits the
            three enolisation branch-point families and nothing else, so both
            Heyns and Amadori get the ionisation term ONCE and the route term
            never.  Asserted at import and in
            `tests/scientific/test_wave_s1b_ph_aw_routing_2026_08.py::
            test_no_family_can_receive_the_ph_physics_twice`.
          * ungated lanes (`kinetics.py`, `pathway_ranker.py`,
            `cantera_export.py`) — these call this function DIRECTLY and never
            call `_ionization_correction`, so there the substring is the only pH
            term either family receives.  Again once.
        """
        if not reaction_family:
            return 1.0
            
        fam = reaction_family.lower()
        
        # 1. Schiff Base: Gaussian peak at pH 5.5
        mult = 1.0
        if any(x in fam for x in ["schiff", "condensation"]):
            mult *= (1.0 + 2.0 * self._gaussian(self.pH, 5.5, 1.0))

        # 2. 2,3-enolisation / Pyrazine / Strecker / Amadori / Heyns (Alkaline favored)
        elif any(x in fam for x in ["2,3", "2_3", "pyrazine", "strecker", "amadori", "heyns", "nitrogen_heterocycle"]):
            base_mult = 0.2 + 8.0 * self._sigmoid(self.pH, 6.5, 1.5)
            mult *= max(0.01, base_mult)

        # 3. 1,2-enolisation / Furan / Thiol / Thio / Cysteine / [Generic Enolisation] (Acidic favored)
        elif any(x in fam for x in ["1,2", "1_2", "furan", "thiol", "thio", "cysteine", "enolisation", "oxygen_heterocycle"]):
            mult *= (1.0 + 4.0 * (1.0 - self._sigmoid(self.pH, 6.0, 2.0)))
            
        # 4. Metal Catalysis (Heme/Iron Enrichment)
        # Heme is a potent promoter of lipid oxidation and certain Maillard steps.
        if self.metal_catalyst and self.metal_catalyst.lower() == "heme":
            if "oxidation" in fam or "radical" in fam:
                mult *= 5.0  # Heme significantly accelerates radical generation
            elif "pyrazine" in fam:
                mult *= 1.5  # Iron catalysis of nitrogen cyclization
                
        return mult

    def get_rate_constant(
        self,
        pathway_type: str,
        ea_override_kcal: float = None,
        *,
        reactant_labels: Optional[Sequence[str]] = None,
        time_minutes: Optional[float] = None,
    ) -> float:
        """
        Arrhenius rate constant: k = A * exp(-Ea / RT)
        ea_override_kcal: barrier in kcal/mol from results_db.py
        
        Default Activation Energies (kJ/mol):
        - amadori_formation: 75.0
        - strecker_degradation: 85.0
        - thiol_formation: 95.0
        - pyrazine_formation: 110.0
        - furfural_formation: 70.0
        - acrylamide_formation: 130.0
        """
        if time_minutes is not None and self.has_dynamic_profile:
            T_K = self.temperature_celsius_at_time(float(time_minutes)) + 273.15
        else:
            T_K = self.temperature_kelvin
        
        # Ea in kJ/mol
        ACTIVATION_ENERGIES_KJ = {
            "amadori_formation": 75.0,
            "strecker_degradation": 85.0,
            "thiol_formation": 95.0,       
            "pyrazine_formation": 110.0,
            "furfural_formation": 70.0,
            "acrylamide_formation": 130.0,
        }
        
        if ea_override_kcal is not None:
            Ea_J = ea_override_kcal * 4184.0 # kcal/mol -> J/mol
        else:
            # Match family to default kJ/mol map
            fam = pathway_type.lower()
            ea_kj = 90.0 # Default fallback
            for key, val in ACTIVATION_ENERGIES_KJ.items():
                if key in fam or fam in key:
                    ea_kj = val
                    break
            Ea_J = ea_kj * 1000.0
        
        barrier_kcal = Ea_J / 4184.0
        
        # Apply pH ionization correction (free-nucleophile fraction)
        ph_factor = self._ionization_correction(pathway_type)

        # Apply the pH route-selection correction (1,2- vs 2,3-enolisation).
        #
        # WAVE S1b (2026-08-27): THIS IS THE WIRING FIX.  `get_ph_multiplier`
        # existed, was documented, and was unit-tested for the correct signs, and
        # no prediction had ever called it — the prediction path enters here.
        # It is connected at exactly ONE point, and only for the families in
        # `_ENOLISATION_ROUTE_PH_FAMILIES`, which is disjoint from the families
        # `_ionization_correction` claims (asserted at import).  See the block at
        # the top of this module for why the two are different physics.
        route_factor = self._enolisation_route_ph_correction(pathway_type)

        # Apply water activity correction (Labuza)
        aw_factor = self._water_activity_correction(pathway_type)
        donor_factor = get_donor_reactivity_multiplier(pathway_type, reactant_labels=reactant_labels)

        return arrhenius_rate_constant(
            barrier_kcal,
            T_K,
            family=pathway_type,
            multiplier=ph_factor * route_factor * aw_factor * donor_factor,
        )

    def _enolisation_route_ph_correction(self, pathway_type: str) -> float:
        """`get_ph_multiplier`, gated to the enolisation branch point.

        The gate is what makes the connection honest.  `get_ph_multiplier` matches
        by substring internally, and its substrings reach roughly half the emitted
        families — including `Aminoketone_Condensation`, which its "condensation"
        key would hand the ACID-peaked Schiff Gaussian, i.e. the exact inversion
        this wave exists to fix.  Only the three branch-point families are routed
        through it; for everything else this returns 1.0 and the family's pH
        physics is whatever `_ionization_correction` gives it.
        """
        if not pathway_type:
            return 1.0
        if pathway_type.lower() not in _ENOLISATION_ROUTE_PH_FAMILIES:
            return 1.0
        return self.get_ph_multiplier(pathway_type)

    def _ionization_correction(self, pathway_type: str) -> float:
        """
        Henderson-Hasselbalch based correction for reactive species ionization.
        Replaces arbitrary sigmoids.
        """
        fam = pathway_type.lower()
        if _releases_rather_than_attacks_with_the_amine(fam):
            return 1.0

        # WAVE S1b (2026-08-27): explicit family sets replace substring keys.
        # The pKa values and the Henderson-Hasselbalch form are UNCHANGED; what
        # changed is which families reach them.  See the module header.
        #
        # Amine reactions require the deprotonated form (R-NH2).
        if fam in _ALPHA_AMINO_NUCLEOPHILE_FAMILIES:
            pKa = 8.0  # Common alpha-amino pKa
            return 1.0 / (1.0 + 10**(pKa - self.pH))

        # Pyrazines often peak at slightly alkaline pH.  The engine's pyrazine
        # step is `Aminoketone_Condensation`; the old `"pyrazine" in fam` key
        # matched NO emitted family and this branch was dead at every pH.
        if fam in _AMINOKETONE_NUCLEOPHILE_FAMILIES:
            return 1.0 / (1.0 + 10**(6.5 - self.pH))

        return 1.0

    def _water_activity_correction(self, pathway_type: str) -> float:
        """
        Correction for water activity effect on Maillard reaction rate.
        Based on Labuza & Saltmarch (1981).
        Rate peaks around aw=0.65, drops at both extremes.

        WAVE S1b (2026-08-27) — ROUTING REPAIR, NOT A RESHAPING.  Both curves below
        are byte-identical to the ones this function has always carried.  What
        changed is reach: the peaked curve was hitting 3 of the 29 emitted families
        and the dehydration curve was keyed on the substring "furfural", which
        matches NO emitted family and so was dead.  Family membership is now
        decided by measured net water stoichiometry (see the module header), which
        is why the dehydration arm now covers the furan/HMF track it was written
        for, and why a water-CONSUMING arm exists at all.
        """
        aw = self.water_activity
        fam = pathway_type.lower()

        # NOTE (Wave S1b): the Wave H `_releases_rather_than_attacks_with_the_amine`
        # short-circuit that used to stand HERE has been removed, and removing it
        # was a bug fix, not a relaxation.  That guard matched the substring
        # "enolisation", so it returned 1.0 for `Enolisation_1_2` — the
        # 3-deoxyosone -> furfural/HMF dehydration, i.e. THE FURAN TRACK, the one
        # family this repair most needed to reach — before any set membership was
        # even consulted.  It was written to stop `Enolisation_2_3_Amadori` picking
        # up the AMINE Labuza curve from the substring "amadori"; the explicit sets
        # below do that strictly better, because `Enolisation_2_3_Amadori` is
        # net-zero in water and appears in NO set, so it still returns 1.0 and
        # `tests/unit/test_wave_h_2026_08.py` still pins it.  The helper is kept and
        # is still consulted by `_ionization_correction`.
        if fam in _LABUZA_EMPIRICAL_FAMILIES:
            # The empirical overall-browning response.  Not mass action: these
            # families are net-zero in water.
            if aw < 0.2:
                return 0.1
            elif aw < 0.65:
                # Linear ramp
                return 0.1 + (aw - 0.2) / (0.65 - 0.2) * 0.9
            else:
                # Product inhibition/dilution beyond peak
                return 1.0 - (aw - 0.65) / (1.0 - 0.65) * 0.4

        if fam in _WATER_RELEASING_FAMILIES:
            # Dehydration / condensation: water is a PRODUCT, so water pushes back.
            # This is the branch the "furfural" substring was reaching for.
            return max(0.1, 1.3 - aw)

        if fam in _WATER_CONSUMING_FAMILIES:
            # Hydrolysis: water is a REACTANT, so the rate tracks its activity.
            # Mass action, no free parameter; floored like the arm above.
            return max(0.1, aw)

        return 1.0
        
    def get_arrhenius_multiplier(self, activation_barrier_kcal: float) -> float:
        """
        Returns the relative rate multiplier e^(-Ea / RT).
        Uses R = 1.987 cal/(mol*K).
        """
        if activation_barrier_kcal <= 0:
            return 1.0
            
        R_kcal = 0.001987
        exponent = - (activation_barrier_kcal) / (R_kcal * self.temperature_kelvin)
        
        import math
        try:
            return math.exp(exponent)
        except OverflowError:
            return 0.0
            
    def get_water_activity_multiplier(self) -> float:
        """
        Calculates the kinetic multiplier for water activity (aw).
        
        Physics:
        Maillard reactions peak at aw ≈ 0.6-0.8 (typically around 0.7).
        - At low aw (glassy state), diffusion becomes limiting.
        - At high aw, water acts as a diluent and product inhibitor (Le Chatelier).
        
        Replaced the piecewise linear model with a continuous Gaussian curve:
        mult = exp(-0.5 * ((aw - 0.7) / 0.15)^2)
        """
        # Center peak around aw = 0.7, with a standard deviation of 0.15
        val = self._gaussian(self.water_activity, center=0.7, sigma=0.15)
        # Prevent it from dropping exactly to zero to avoid numerical singularities in ODEs
        return max(0.01, val)
