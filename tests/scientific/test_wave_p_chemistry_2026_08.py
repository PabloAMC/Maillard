"""Wave P (2026-08-27) — the six evidence-grounded chemistry changes, pinned.

Every assertion here pins a number the wave MEASURED, not a number it wanted. Where a
number got worse it is pinned worse, with the cause in the comment. The point of this
file is that a future change to any of the six lanes has to move a named pin and say why.

The six items:
  1. `thiol_addition_pentodiulose` refitted against Hofmann1998 (28.60 -> 26.35)
  2. the Hofmann & Schieberle C2 + C3 recombination lane to MFT
  3. norfuraneol -> 2,3-pentanedione -> 2-mercapto-3-pentanone
  4. nonanal cleaved from the OLEATE pool, not the linoleate pool
  5. fructose reaches HMF by its own ring-retained dehydration
  6. the DMHF `[HH]` pool gate removed in favour of the amino-acid reduction
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest
from rdkit import Chem

from src.barrier_constants import (
    DEFAULT_BARRIER,
    FAST_BARRIERS,
    _canonical_fast_family,
    get_barrier,
)
from src.benchmark_validation import evaluate_benchmark
from src.conditions import ReactionConditions
from src.lipid_oxidation import (
    MARKER_HYDROPEROXIDE_POOL,
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    predict_hexanal_generation,
)
from src.smirks_engine import SmirksEngine, Species

_COND = ReactionConditions(pH=5.5, temperature_celsius=150.0)

_RIBOSE = Species("D-Ribose", "O=CC(O)C(O)C(O)CO")
_GLUCOSE = Species("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO")
_FRUCTOSE = Species("D-Fructose", "OCC(O)C(O)C(O)C(=O)CO")
_CYSTEINE = Species("L-Cysteine", "NC(CS)C(=O)O")
_GLYCINE = Species("Glycine", "NCC(=O)O")

_HOFMANN = Path("data/benchmarks/cys_ribose_140C_Hofmann1998.json")

#: The six families Wave P introduced.
_NEW_FAMILIES = (
    "Mercaptoketone_Formation",
    "Mercaptoketone_Aldol_Addition",
    "Mercaptoketone_Cyclodehydration",
    "Furanone_Reductive_Opening",
    "Furanone_Amino_Acid_Reduction",
    "Fructofuranosyl_Dehydration",
)


def _enumerate(precursors, generations: int = 3):
    return SmirksEngine(conditions=_COND).enumerate(list(precursors), generations)


def _atoms(species_list):
    counts: collections.Counter = collections.Counter()
    for sp in species_list:
        mol = Chem.MolFromSmiles(sp.smiles)
        assert mol is not None, f"unparseable SMILES {sp.smiles!r}"
        for atom in Chem.AddHs(mol).GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts


# ── family plumbing ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("family", _NEW_FAMILIES)
def test_new_families_have_explicit_barrier_keys_and_do_not_fall_through(family):
    """The canonicalisation fallthrough is a KNOWN defect class in this table.

    Wave G1 fix 8 found eight emitted families silently inheriting
    ``DEFAULT_BARRIER = 45.0`` (a ~39,000-year half-life at 150 C), i.e. switched off
    without anyone noticing. `_canonical_fast_family` also has substring rules that
    would have caught two of these names in the wrong bucket. Both halves are pinned:
    the key resolves to ITSELF (not to a class analogue), and the value is not the
    default.
    """
    key = _canonical_fast_family(family)
    assert key == family.lower(), (
        f"{family} canonicalises to {key!r}, i.e. it is inheriting another family's "
        "barrier through a substring rule instead of carrying its own explicit key"
    )
    assert key in FAST_BARRIERS, f"{family} has no explicit FAST_BARRIERS entry"
    assert get_barrier(family)[0] != DEFAULT_BARRIER, (
        f"{family} falls through to DEFAULT_BARRIER"
    )
    rationale = FAST_BARRIERS[key][1]
    assert "ESTIMATED" in rationale, (
        f"{family}'s barrier must be labelled ESTIMATED — none of the Wave P values "
        "has a kinetic measurement behind it"
    )


def test_new_families_collect_no_accidental_ph_correction_from_their_name():
    """No Wave P family may pick up a pH correction it was not deliberately given.

    RE-PINNED AND NARROWED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR).
    NOT RELAXED: what changed is that the water-activity half of this guard was asserting
    the ABSENCE of a correction that has since been deliberately, explicitly ADDED, so
    keeping it would have been pinning a defect.

    ORIGINAL TEXT, kept verbatim because its reasoning is still the reason this test exists:
      "`src/conditions.py` classifies reactions for the pH/aw corrections by SUBSTRING. A
      family name containing "strecker", "amadori", "schiff", "pyrazine" or "furfural"
      silently picks up corrections that can be worth hundreds of x. The Wave P families are
      named so that they receive the SAME treatment as the sibling steps they extend or
      replace (`Furanone_Cyclisation`, `Furanone_Formation`), which receive none. This test
      is the guard on that: it fails the day someone renames one of them into a trigger word
      and changes a prediction by accident."

    WHAT WAVE S1b CHANGED. `src/conditions.py` no longer classifies by substring at all --
    both corrections now key on EXPLICIT FAMILY-NAME SETS, which removes the accident this
    test was written to catch at its source. The pH half is therefore STRENGTHENED into a
    real invariant (no Wave P family is in either amine set, so all of them must return
    exactly 1.0, by construction rather than by luck of spelling).

    The aw half is REPLACED by its correct successor below, because five of the six Wave P
    families are net water-RELEASING and are now deliberately in
    `_WATER_RELEASING_FAMILIES`: Mercaptoketone_Formation (+1), Mercaptoketone_Cyclodehydration
    (+2), Furanone_Reductive_Opening (+1), Furanone_Amino_Acid_Reduction (+1),
    Fructofuranosyl_Dehydration (+3). Only Mercaptoketone_Aldol_Addition is net-zero and
    still uncorrected. Membership is decided by MEASURED stoichiometry, and that measurement
    is what the successor test pins.
    """
    for family in _NEW_FAMILIES:
        fam = family.lower()
        assert _COND._ionization_correction(fam) == 1.0, (
            f"{family} picked up an amine-ionisation correction it was not given "
            f"deliberately -- check `_ALPHA_AMINO_NUCLEOPHILE_FAMILIES` and "
            f"`_AMINOKETONE_NUCLEOPHILE_FAMILIES` in src/conditions.py"
        )


def test_new_families_water_activity_membership_matches_their_measured_stoichiometry():
    """The aw correction reaches exactly the Wave P families that shed water, and no others.

    NEW 2026-08-27 (Wave S1b). This replaces the aw half of the guard above. Wave S2
    measured `_water_activity_correction` reaching 3 of the 29 families this engine emits,
    with its dehydration branch keyed on the substring "furfural" -- which matches NO emitted
    family, so the branch was dead. Membership is now decided by net water produced per step,
    counted directly off the enumerated steps.

    The numbers in the mapping below are that measured count. If a Wave P family's
    stoichiometry changes, this test fails and the family's membership must be re-derived --
    it must NOT be edited to match.
    """
    net_water_per_step = {
        "Mercaptoketone_Formation": +1,
        "Mercaptoketone_Aldol_Addition": 0,
        "Mercaptoketone_Cyclodehydration": +2,
        "Furanone_Reductive_Opening": +1,
        "Furanone_Amino_Acid_Reduction": +1,
        "Fructofuranosyl_Dehydration": +3,
    }
    assert set(net_water_per_step) == set(_NEW_FAMILIES)
    # aw 0.8 is the ReactionConditions default: the dehydration branch gives 1.3-0.8 = 0.5.
    assert _COND.water_activity == pytest.approx(0.8)
    for family, net_water in net_water_per_step.items():
        got = _COND._water_activity_correction(family.lower())
        if net_water > 0:
            assert got == pytest.approx(0.5), (
                f"{family} sheds {net_water} water per step, so it must carry the "
                f"dehydration inhibition; got {got}"
            )
        else:
            assert got == 1.0, (
                f"{family} is net-zero in water, so it must carry NO mass-action term; "
                f"got {got}"
            )


# ── item 1: the refit ──────────────────────────────────────────────────────────

def test_pentodiulose_barrier_is_the_wave_p_fit_and_carries_the_conversion_caveat():
    """28.60 kcal/mol, ESTIMATED again -- the Wave P fit is REVERTED and its record RETRACTED.

    Wave N shipped 28.60 (the un-fitted `thiol_addition` class value) and declined to
    refit, because refitting re-couples the only sulfur anchor to a fitted constant.
    Wave P did it with owner approval, AFTER the wave's chemistry additions so the fit
    sees the network that ships. The load-bearing part is not the number: it is that
    the rationale must carry the fit target's own `content_verification_note` verbatim,
    so nobody can read 26.35 as a measured barrier.

    RE-PINNED 2026-08-27 (Wave S2c). 26.35 -> 28.60. Wave N's stated worry -- that refitting
    "re-couples the only sulfur anchor to a fitted constant" -- turned out to understate the
    problem: there was no anchor. Wave S2b traced `cys_ribose_140C_Hofmann1998`'s MFT 342 /
    FFT 200 ppb to data/benchmarks/maillard_validation_benchmarks.md section 1.3, an
    abstract-reconstructed range table committed in c7efbbc, the SAME commit that created the
    benchmark JSON; both values are interior points of two invented, OVERLAPPING mol % bands
    (0.0300 mol % -> 342.5 -> 342 ppb; the FFT band's geometric mean 0.017321 mol % -> 197.8
    -> 200 ppb, on the file's own unattested 10 mM / MW 114.17 basis). So the fit was against
    the repo's own guess. The constant is reverted, results/validation/
    sulfur_barrier_refit_pentodiulose.{json,md} is RETRACTED, and this test goes back to
    asserting what it was written to assert: that this constant claims no provenance it does
    not have. THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.

    WHAT THIS TEST NOW OWNS: the rationale must keep the WHOLE history -- estimate, fit,
    revert -- because deleting the fit from the record would hide that it ever happened, and
    the Wave K caveat text must survive verbatim while being marked superseded (the problem
    was never the undocumented mol%->ppb conversion; it is that there is no measurement on the
    far end of it).
    """
    value, rationale = FAST_BARRIERS["thiol_addition_pentodiulose"]
    assert value == pytest.approx(28.60, abs=1e-9)
    assert value != pytest.approx(26.35, abs=1e-9), (
        "back at the Wave P fitted value. That fit's sole target was "
        "cys_ribose_140C_Hofmann1998, whose values are a repo-internal derivation -- "
        "refitting against it is circular. See tasks/audit_remediation.md '## Wave S2b'."
    )
    # The whole history, in order, still readable at the point of use.
    assert "ESTIMATED" in rationale
    assert "FITTED 2026-08-27 (Wave P item 1)" in rationale
    assert "REVERTED 2026-08-27 (Wave S2c)" in rationale
    assert "cys_ribose_140C_Hofmann1998" in rationale
    assert "ZERO ABSOLUTE LITERATURE ANCHORS" in rationale
    # The verbatim Wave K caveat, and the sentence that says what it means for THIS constant.
    # Both are KEPT and both are marked superseded rather than deleted.
    assert "mol%->ppb conversion" in rationale
    assert "NOT documented anywhere in this repo" in rationale
    assert "UNVERIFIED" in rationale
    assert "LOCALISED HERE" in rationale
    assert "SUPERSEDED" in rationale
    # The boundary honesty from the retired fit: the profile minimum was at the range floor,
    # so even the adopted value could not reach its (fabricated) target.
    assert "PROFILE MINIMUM SITS AT THE RANGE FLOOR" in rationale
    assert "sulfur_barrier_refit_pentodiulose" in rationale and "RETRACTED" in rationale


def test_hofmann1998_after_the_refit():
    """Both rows improved; both are still misses; the benchmark is still a fit target.

    Movement across Wave P, measured 2026-08-27:
        MFT  151.87 -> 242.38 ppb vs 342   (2.2519x under -> 1.4110x under)
        FFT  243.72 -> 217.98 ppb vs 200   (1.2186x over  -> 1.0900x over)
    FFT was NOT fitted. It co-moves because the two lanes draw on the same upstream
    sugar flux — which is exactly why ONE knob was fitted against these two rows and
    not two.
    """
    evaluation = evaluate_benchmark(_HOFMANN)
    predicted = {c.compound: c.predicted_ppb for c in evaluation.comparisons}
    # RE-PINNED 2026-08-27 (Wave S1). MFT 242.38 -> 283.59, FFT 217.99 -> 297.28.
    # CAUSE: the flux propagator became ADDITIVE over parallel channels
    # (`src/recommend.py::predict_from_steps`, `_route_channel_id`). BOTH compounds are
    # reached by two enumerated routes here, so both rose. NOTHING WAS RE-FITTED to
    # compensate: `thiol_addition_pentodiulose` is still the 26.35 this wave fitted, and
    # the wave that changed the propagator deliberately did not touch a barrier.
    # HALF OF THIS GOT WORSE and is pinned worse: MFT 1.4110x under -> 1.2060x under
    # (better), FFT 1.0900x over -> 1.4864x over (WORSE), benchmark MALE 0.0935 -> 0.1267
    # dex and max_ratio 1.4110 -> 1.4864, i.e. the contract now fails on BOTH of its
    # criteria again rather than on one. The two lanes share their upstream trunk, so a
    # barrier that pushed FFT down would push MFT down with it.
    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH ROUTING REPAIR). MFT 283.59 -> 154.85,
    # FFT 297.28 -> 267.50. NO BARRIER MOVED AND NOTHING WAS RE-FITTED.
    # CAUSE: `get_ph_multiplier` -- the enolisation route-selection term -- had never been
    # called on the prediction path and now is. At this benchmark's pH 5.0 it gives
    # `Enolisation_1_2` (the 3-deoxyosone -> furfural/FFT arm) a 4.5x acid boost and
    # `Enolisation_2_3_Amadori` (the 1-deoxyosone -> MFT arm) ~1.0, so the fixed volatile
    # budget moves from the MFT arm to the FFT/furfural arm.
    # BOTH HALVES ARE WORSE ON THE CONTRACT AND ARE PINNED WORSE: max_ratio 1.4864 ->
    # 2.2086, MALE 0.1267 -> 0.2352 dex. MFT went 1.2060x under -> 2.2086x under (WORSE);
    # FFT went 1.4864x over -> 1.3375x over (better). The untouched contract is
    # 1.45x / 0.09 dex and it now fails both criteria by more.
    # NOT A CONFLICT WITH A MEASUREMENT -- Wave S1b's first draft said it was, and Wave S2b
    # (same day) showed otherwise. The 342 / 200 ppb targets were derived INSIDE THIS
    # REPOSITORY from data/benchmarks/maillard_validation_benchmarks.md section 1.3, an
    # abstract-reconstructed range table committed in the SAME commit as the benchmark file;
    # both values are interior points of two INVENTED, OVERLAPPING bands (MFT 228-571 ppb,
    # FFT 114-342 ppb). The MFT > FFT ordering is midpoint selection, not measurement, and
    # the 1.45x / 0.09 dex contract is ~1.7x tighter than its own source band. The pH
    # mechanism and this degradation are both real; the yardstick is not. See the file's
    # `content_verification_note.wave_s2_followup` and '## Wave S2b' in
    # tasks/audit_remediation.md.
    # RE-PINNED 2026-08-27 (Wave S2c -- THE ANCHOR RETIREMENT). MFT 154.85 -> 78.09,
    # FFT 267.50 -> 293.67. CAUSE: `thiol_addition_pentodiulose` REVERTED 26.35 -> 28.60,
    # because the Wave P fit that produced 26.35 had exactly one target and that target is
    # not a measurement (see the module note above and tasks/audit_remediation.md
    # '## Wave S2b'). max_ratio 2.2086 -> 4.3797, MALE 0.2352 -> 0.4041 dex.
    # BOTH ROWS ARE WORSE AND BOTH ARE PINNED WORSE. Nothing was clawed back, and the
    # benchmark's 1.45x / 0.09 dex contract was RETIRED rather than widened -- it now
    # inherits the global free-precursor default 1.5x / 0.10 dex and fails that by more.
    # READ THE MISS FOR WHAT IT IS: an error against a yardstick this repository invented.
    # It is NOT evidence about the chemistry, in either direction, and it must never be
    # quoted as accuracy against literature.
    # RE-PINNED 2026-08-28 (Wave X -- THE NORFURANEOL ROUTE RETURNS AS A PARALLEL CHANNEL).
    # MFT 78.09 -> 119.08, FFT 293.67 -> 282.87. NO BARRIER MOVED: the Wave X fit against
    # Hofmann Table 4 was REJECTED by its own isotope gate and
    # `furanone_reductive_sulfhydrylation` ships at the un-fitted `thiol_addition` class
    # value 28.60 (results/validation/furanone_reductive_sulfhydrylation_refit_hofmann.json).
    # CAUSE: `norfuraneol + H2S + 2[H] -> MFT + 2 H2O` was re-added as a THIRD MFT channel
    # (src/reaction_templates.py::_norfuraneol_mft_parallel_route). Wave N had removed it on
    # Cerny & Davidek's spiking result; Hofmann & Schieberle 1998 Table 4 MEASURES it when
    # norfuraneol is FED (211.2 ug, 0.19 mol %, 14x the yield from ribose), and under the
    # Wave S1 additive propagator "real but minor in situ" is expressible. The isotope
    # constraint is now a REGRESSION TEST rather than the step's absence:
    # tests/scientific/test_wave_x_step_level_2026_08.py::
    # test_norfuraneol_route_stays_a_minor_share_of_ribose_mft_flux.
    # MFT is WORSE against this file's number and is pinned worse (4.38x under -> 2.87x
    # under is *closer* to 342, but 342 IS NOT A MEASUREMENT -- see the paragraph above --
    # so "closer" here means nothing about the chemistry); FFT co-moves DOWN because the two
    # lanes share upstream flux and MFT's share of a fixed budget grew.
    # THE NUMBER THAT MATTERS IS ON THE REAL ANCHOR, NOT HERE: on
    # hofmann1998_ribose_cysteine_145C_20min_pH5, whose 198 ppb IS a primary-source
    # measurement, MFT went 468.58 -> 713.74, i.e. 2.37x over -> 3.61x OVER. The wave made
    # the real panel row worse and says so in tasks/audit_remediation.md '## Wave X' (b).
    assert predicted["2-methyl-3-furanthiol"] == pytest.approx(119.08, rel=1e-3)
    assert predicted["2-furfurylthiol"] == pytest.approx(282.87, rel=1e-3)
    # Still under / still over — neither direction flipped.
    assert predicted["2-methyl-3-furanthiol"] < 342.0
    assert predicted["2-furfurylthiol"] > 200.0


# ── item 2: the C2 + C3 lane, and the finding it produced ──────────────────────

def test_c2_c3_recombination_lane_exists_and_is_exactly_balanced():
    """Hofmann & Schieberle 1998's highest-yielding MFT system, finally expressible."""
    steps = _enumerate([_RIBOSE, _CYSTEINE])
    families = [s.reaction_family for s in steps]
    for family in (
        "Mercaptoketone_Formation",
        "Mercaptoketone_Aldol_Addition",
        "Mercaptoketone_Cyclodehydration",
    ):
        assert family in families, f"{family} does not fire in ribose + cysteine"
    for step in steps:
        if step.reaction_family.startswith("Mercaptoketone"):
            assert _atoms(step.reactants) == _atoms(step.products), (
                f"{step.reaction_family} unbalanced"
            )


def test_hexose_system_also_reaches_the_c2_c3_lane():
    """Measured, and it corrected a guess.

    The lane needs glycolaldehyde, which `_retro_aldol_fragmentation` emits only from
    the PENTOSE deoxyosone — so the lane was expected to be pentose-only. It is not:
    in glucose + cysteine the identical species arrives by
    glyoxal -> Strecker -> 2-aminoethanal -> hydrolytic deamination. Pool identity is
    by canonical SMILES, not by label, which is what connects the longer route.
    """
    families = {s.reaction_family for s in _enumerate([_GLUCOSE, _CYSTEINE])}
    assert "Mercaptoketone_Cyclodehydration" in families


def test_the_second_mft_channel_now_contributes_after_the_wave_s1_propagator_fix():
    """WAVE P'S FINDING, AND ITS RESOLUTION. Renamed and inverted 2026-08-27 (Wave S1).

    WHAT WAVE P PINNED, verbatim from the version this replaces: the C2+C3 lane's
    "MEASURED CONTRIBUTION TO THE SHIPPED PREDICTION: EXACTLY ZERO", because
    `src/recommend.py` relaxed to the LOWEST-SPAN path per product and never summed
    parallel channels. Its own failure message said this failure would be the
    notification that the propagator had become additive. It has, and this is that
    re-pin.

    MEASURED NOW, on cys_ribose_140C_Hofmann1998 at unchanged constants:

        pentodiulose lane alone (C2+C3 disabled)   217.25 ppb
        both lanes                                 283.59 ppb

    NOTE WHAT THIS IS NOT. Wave P predicted 242.38 + 71.02 = 313.39 ppb "if the two MFT
    channels are genuinely independent". Two things falsify that arithmetic and both are
    measured, not argued. First, the two lanes DO share their rate-limiting step (the
    trunk `Amadori_Rearrangement`), so the naive sum was never the right combination —
    see tests/scientific/test_wave_s1_additive_flux_2026_08.py. Second, and more
    importantly, the volatile budget is FIXED: adding a channel to MFT changes MFT's
    SHARE of that budget, so the single-lane figures cannot simply be added, and the
    pentodiulose-alone number itself falls (242.38 -> 217.25) once the competing FFT
    channels also become additive. The 313.39 was never obtainable.
    """
    import src.smirks_engine as engine_module

    baseline = {
        c.compound: c.predicted_ppb for c in evaluate_benchmark(_HOFMANN).comparisons
    }
    original = engine_module._c2_c3_mft_recombination
    try:
        engine_module._c2_c3_mft_recombination = lambda pool: []
        without = {
            c.compound: c.predicted_ppb
            for c in evaluate_benchmark(_HOFMANN).comparisons
        }
    finally:
        engine_module._c2_c3_mft_recombination = original

    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH ROUTING REPAIR). 217.25 -> 128.40 and
    # 283.59 -> 154.85. NOT a propagator change: both lanes fell by the same mechanism,
    # the 4.5x acid boost that `Enolisation_1_2` now receives at this benchmark's pH 5.0,
    # which moves budget share from the MFT arm to the FFT/furfural arm. THE QUANTITY THIS
    # TEST EXISTS TO MEASURE IS UNCHANGED IN KIND: the C2+C3 lane still contributes, and
    # the both-lanes/pentodiulose-alone ratio is 1.2060x (was 1.3054x under Wave S1).
    # RE-PINNED 2026-08-27 (Wave S2c -- THE BARRIER REVERT). 128.40 -> 46.82 and
    # 154.85 -> 78.09. `thiol_addition_pentodiulose` 26.35 -> 28.60; both lanes run through
    # it, so both fell. THE PROPERTY THIS TEST OWNS IS UNCHANGED and in fact strengthened:
    # the C2+C3 lane still contributes, and its share of MFT ROSE (ratio 1.2060 -> 1.6678)
    # because raising the pentodiulose barrier makes the pentodiulose-only lane weaker
    # relative to the C2+C3 lane, which does not run that step.
    # RE-PINNED 2026-08-28 (Wave X -- A THIRD MFT CHANNEL). 46.82 -> 90.10 and 78.09 ->
    # 119.08. NO BARRIER MOVED (the Wave X fit was rejected by its isotope gate). CAUSE: the
    # norfuraneol -> MFT step returned as a parallel channel, so the "C2+C3 disabled" arm now
    # carries TWO lanes (pentodiulose + norfuraneol) rather than one, which is why the
    # without-figure nearly doubled. THE PROPERTY THIS TEST OWNS IS UNCHANGED: the C2+C3 lane
    # still contributes, and it must, or the propagator has gone back to winner-takes-all.
    assert without["2-methyl-3-furanthiol"] == pytest.approx(90.10, rel=1e-3)
    assert baseline["2-methyl-3-furanthiol"] == pytest.approx(119.08, rel=1e-3)
    assert baseline["2-methyl-3-furanthiol"] > without["2-methyl-3-furanthiol"], (
        "the C2+C3 lane must contribute to predicted MFT. If this fails the flux "
        "propagator has gone back to winner-takes-all selection."
    )
    # And it is a real contribution, not rounding: +20.6% on the flagship number.
    # RE-PINNED 2026-08-27 (Wave S1b): 1.3054 -> 1.2060. The C2+C3 lane's SHARE of MFT fell
    # because the pH routing repair moved budget away from the MFT arm as a whole; the lane
    # is still contributing, which is the property this test owns.
    # RE-PINNED 2026-08-27 (Wave S2c): 1.2060 -> 1.6678, i.e. +66.8% on the flagship number
    # rather than +20.6%. The lane got MORE important, not less, and for a mechanical reason:
    # the reverted barrier sits on the pentodiulose lane only.
    # RE-PINNED 2026-08-28 (Wave X): 1.6678 -> 1.3216, i.e. +32.2% on the flagship number
    # rather than +66.8%. The C2+C3 lane's SHARE fell, and it fell for an arithmetic reason
    # with no bearing on the lane itself: the denominator gained a whole extra channel (the
    # returned norfuraneol route). The lane's own contribution did not shrink.
    assert baseline["2-methyl-3-furanthiol"] / without["2-methyl-3-furanthiol"] == pytest.approx(
        1.3216, rel=1e-3
    )


# ── item 3: norfuraneol's real sulfur fate ─────────────────────────────────────

def test_norfuraneol_has_consumers_again_and_they_are_the_evidenced_ones():
    """Wave N retired norfuraneol from the MFT lane and left it a dead end.

    Cerny & Davidek 2003 (10.1021/jf026123f) assign it 2-mercapto-3-pentanone and NOT
    MFT; Whitfield & Mottram 1999 (10.1021/jf980980v) supply the 2,3-pentanedione
    intermediate. The isomer matters: 3-mercapto-2-pentanone has a second,
    ribose-derived origin, so this lane must NOT emit it.
    """
    steps = _enumerate([_RIBOSE, _CYSTEINE])
    consumers = {
        s.reaction_family
        for s in steps
        if any(r.label == "norfuraneol" for r in s.reactants)
    }
    # RE-PINNED 2026-08-28 (Wave X). Norfuraneol has a SECOND consumer again:
    # `Furanone_Reductive_Sulfhydrylation`, the Hofmann Figure 1 route to MFT, re-added as a
    # SLOW PARALLEL channel constrained by that paper's Table 4 (211.2 ug MFT from FED
    # norfuraneol, 14x the yield from ribose). This does NOT reinstate the claim Wave N
    # retired. Cerny & Davidek 2003's spiking experiment says the route is unimportant IN
    # SITU, not that it does not exist, and the in-situ share is now asserted to stay a
    # MINORITY in tests/scientific/test_wave_x_step_level_2026_08.py::
    # test_norfuraneol_route_stays_a_minor_share_of_ribose_mft_flux (measured 34.3%, ceiling
    # 0.50 from the source's word "mainly"). The barrier for the new step is NOT the retired
    # `thiol_addition_norfuraneol` 26.85; the fit that would have moved it was rejected by
    # that same gate. Full argument: tasks/audit_remediation.md '## Wave X' (a) and (b).
    # The assertion below still owns the property this test was written for -- that the
    # EVIDENCED norfuraneol fate is present -- and now also pins the exact size of the set,
    # so a third consumer cannot appear silently.
    assert consumers == {"Furanone_Reductive_Opening", "Furanone_Reductive_Sulfhydrylation"}

    products = {p.label for s in steps for p in s.products}
    assert "2,3-pentanedione" in products
    assert "2-mercapto-3-pentanone" in products
    assert "3-mercapto-2-pentanone" not in products, (
        "this lane must not emit the isomer Cerny & Davidek assign partly to ribose"
    )

    # And nothing on this lane may reach MFT: that was the Wave N correction.
    for step in steps:
        if step.reaction_family in {"Furanone_Reductive_Opening"}:
            assert not any(p.label == "2-methyl-3-furanthiol" for p in step.products)


# ── item 4: nonanal comes off the oleate pool ──────────────────────────────────

def test_nonanal_scales_with_oleate_and_hexanal_with_linoleate():
    """`oleic_acid_pct` was declared, populated, registered — and read by nothing."""
    assert MARKER_HYDROPEROXIDE_POOL["Nonanal"] == "oleate"
    assert MARKER_HYDROPEROXIDE_POOL["Hexanal"] == "linoleate"
    assert MARKER_HYDROPEROXIDE_POOL["2-Pentylfuran"] == "linoleate"

    out = predict_hexanal_generation(PEA_LIPID_PROFILE, temp_C=100.0, time_min=45.0)
    ratio = PEA_LIPID_PROFILE.oleic_acid_pct / PEA_LIPID_PROFILE.linoleic_acid_pct
    assert out["total_hydroperoxide_oleate"] == pytest.approx(
        out["total_hydroperoxide"] * ratio, rel=1e-12
    )
    # The branching ratio itself was NOT touched; only the pool it multiplies.
    assert out["nonanal"] == pytest.approx(out["total_hydroperoxide_oleate"] * 0.15, rel=1e-12)


def test_trikusuma_nonanal_recovery_broke_by_exactly_the_oleic_linoleic_ratio():
    """An honest regression, pinned as one, and NOT refitted away.

    `pea_iso/heated_matrix/nonanal` observable_factor was back-solved so this benchmark
    reproduced its own measured 24 ppb — against the pre-Wave-P linoleate-pool nonanal.
    After the substrate correction the row reads 10.56 ppb, i.e. 2.2727x under, which is
    EXACTLY 1 / (22.0 / 50.0). That exact arithmetic identity is the evidence that the
    factor was absorbing a substrate-assignment error and nothing else. Refitting it
    would put the error straight back into the same constant, so it was left alone; see
    the dated note on the record in src/matrix_calibration_registry.py.
    """
    evaluation = evaluate_benchmark(
        Path("data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json")
    )
    predicted = {c.compound: c.predicted_ppb for c in evaluation.comparisons}
    expected_ratio = PEA_LIPID_PROFILE.oleic_acid_pct / PEA_LIPID_PROFILE.linoleic_acid_pct
    assert predicted["nonanal"] == pytest.approx(24.0 * expected_ratio, rel=1e-4)
    # The two rows whose factors were NOT substrate-affected still recover exactly.
    assert predicted["hexanal"] == pytest.approx(782.0, rel=1e-4)
    assert predicted["2-pentylfuran"] == pytest.approx(163.0, rel=1e-4)


def test_soy_nonanal_uses_the_soy_oleic_fraction():
    out = predict_hexanal_generation(SOY_LIPID_PROFILE, temp_C=100.0, time_min=45.0)
    ratio = SOY_LIPID_PROFILE.oleic_acid_pct / SOY_LIPID_PROFILE.linoleic_acid_pct
    assert ratio == pytest.approx(23.0 / 53.0)
    assert out["total_hydroperoxide_oleate"] == pytest.approx(
        out["total_hydroperoxide"] * ratio, rel=1e-12
    )


# ── item 5: fructose reaches HMF by its own route ──────────────────────────────

def test_fructose_reaches_hmf_without_the_glucose_3_deoxyosone():
    """Perez Locas & Yaylayan 2008 exclude 3-deoxyglucosone for fructose.

    The distinction is TOPOLOGICAL: fructose keeps its own ring (Antal 1990). The Heyns
    product and the 3-deoxyosone are still formed — they feed the retro-aldol and
    Strecker lanes — but they no longer carry fructose's HMF.
    """
    steps = _enumerate([_FRUCTOSE, _GLYCINE])
    hmf_producers = {
        s.reaction_family for s in steps if any(p.label == "HMF" for p in s.products)
    }
    assert hmf_producers == {"Fructofuranosyl_Dehydration"}

    direct = next(s for s in steps if s.reaction_family == "Fructofuranosyl_Dehydration")
    assert [r.label for r in direct.reactants] == ["D-Fructose"]
    assert _atoms(direct.reactants) == _atoms(direct.products)

    # The glucose limb is untouched: it still reaches HMF through the 3-deoxyosone.
    glucose_hmf = {
        s.reaction_family
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "HMF" for p in s.products)
    }
    assert glucose_hmf == {"Enolisation_1_2"}


# ── item 6: the DMHF [HH] gate is gone ─────────────────────────────────────────

def test_furaneol_from_glucose_glycine_survives_disabling_the_pyrazine_lane():
    """Red-team H4, second half — the exact measurement Wave L1 published as a defect.

    Wave L1, monkey-patching `_aminoketone_condensation` to return no steps:
        [baseline]                    16 steps | DMHF steps 1 | H2 producers: [pyrazine]
        [aminoketone cond. DISABLED]  14 steps | DMHF steps 0 | H2 producers: []
    i.e. predicted furaneol from glucose was a downstream dependent of pyrazine
    chemistry, while Wang & Ho 2008's CAMOLA experiment establishes the route as real.

    After Wave P the token is not re-sourced, it is GONE: the amino acid is the
    reductant (Blank & Fay 1996; Kerler et al. 2010), the step balances exactly with
    no `[HH]`, and DMHF survives.
    """
    import src.smirks_engine as engine_module

    def dmhf_steps():
        return [
            s
            for s in _enumerate([_GLUCOSE, _GLYCINE])
            if any(p.label == "DMHF" for p in s.products)
        ]

    baseline = dmhf_steps()
    assert len(baseline) == 1

    original = engine_module._aminoketone_condensation
    try:
        engine_module._aminoketone_condensation = lambda pool: []
        without_pyrazine = dmhf_steps()
    finally:
        engine_module._aminoketone_condensation = original

    assert len(without_pyrazine) == 1, (
        "furaneol from glucose/glycine is contingent on pyrazine chemistry again"
    )


def test_the_dmhf_step_consumes_no_reducing_equivalent_token():
    steps = [
        s
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "DMHF" for p in s.products)
    ]
    assert steps, "no DMHF step at all"
    for step in steps:
        assert step.reaction_family == "Furanone_Amino_Acid_Reduction"
        assert not any(r.smiles == "[HH]" for r in step.reactants), (
            "the DMHF step is pool-gated on the reducing-equivalent token again"
        )
        assert _atoms(step.reactants) == _atoms(step.products)
        # The amino acid is the reductant: it leaves as its Strecker aldehyde + CO2 + NH3.
        labels = {p.label for p in step.products}
        assert {"DMHF", "formaldehyde", "CO2", "ammonia", "water"} <= labels


def test_hexose_dmhf_carbons_all_come_from_the_sugar():
    """Wang & Ho 2008's CAMOLA constraint: the C6 skeleton stays intact.

    A 1:1 [13C6]/[12C6]glucose mixture gave a 1:1 [13C6]/[12C6]DMHF mixture with no
    intermediate isotopologues. The amino-acid-coupled step preserves that: the amino
    acid's carbon leaves as its own aldehyde and as CO2, so all six DMHF carbons are
    still the hexose's.
    """
    steps = [
        s
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "DMHF" for p in s.products)
    ]
    step = steps[0]
    sugar_carbons = _atoms(
        [r for r in step.reactants if r.label == "hexose-1-deoxyosone"]
    )["C"]
    assert sugar_carbons == 6
    dmhf_carbons = _atoms([p for p in step.products if p.label == "DMHF"])["C"]
    assert dmhf_carbons == 6
