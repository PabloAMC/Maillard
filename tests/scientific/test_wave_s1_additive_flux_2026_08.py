"""Wave S1 (2026-08-27) — the two STRUCTURAL fixes in `src/recommend.py`, pinned.

  FIX 1  the flux propagator became ADDITIVE over parallel channels. It used to relax to
         the lowest-span path per product and keep THAT route's flux alone, so adding a
         second real route to a compound contributed exactly nothing (Wave P measured
         zero on the flagship MFT number).
  FIX 2  the compound-specific matrix calibration registry became REACHABLE on the
         `matrix_precursor_augmented` lane. Species injected there are labelled by
         canonical SMILES, so the name-keyed registry missed and the observability factor
         silently applied as 1.0 (Wave O finding (f)).

Every number here was MEASURED on the tree it pins. Where a number got worse it is pinned
worse, with the cause in the comment.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

import src.matrix_calibration_registry as registry
import src.recommend as recommend
from src.benchmark_validation import evaluate_benchmark
from src.pathway_extractor import ElementaryStep, Species
from src.recommend import Recommender, _route_channel_id

_HOFMANN = Path("data/benchmarks/cys_ribose_140C_Hofmann1998.json")
_SOY_SNAPSHOT = Path("data/benchmarks/soy_isolate_ribose_cysteine_100C_45min_Internal2026.json")


def _predicted(bench: Path) -> dict[str, float]:
    return {c.compound: float(c.predicted_ppb) for c in evaluate_benchmark(bench).comparisons}


# ── FIX 1: the propagator is additive, and it is additive over the RIGHT things ──────

#: A synthetic two-route network. `A` is the only precursor; `P` is reachable through the
#: independent intermediates `I1` and `I2`. Nothing here is chemistry — the point is the
#: propagator's arithmetic, isolated from any barrier table.
_A = Species("A", "CCCCO")
_I1 = Species("I1", "CCCCC")
_I2 = Species("I2", "CCCCCC")
_P = Species("P", "Cc1ccccc1")

_ROUTE_1 = [
    ElementaryStep(reactants=[_A], products=[_I1], reaction_family="Synthetic_Route_1_Entry"),
    ElementaryStep(reactants=[_I1], products=[_P], reaction_family="Synthetic_Route_1_Exit"),
]
_ROUTE_2 = [
    ElementaryStep(reactants=[_A], products=[_I2], reaction_family="Synthetic_Route_2_Entry"),
    ElementaryStep(reactants=[_I2], products=[_P], reaction_family="Synthetic_Route_2_Exit"),
]
_SYNTHETIC_BARRIERS = {
    "CCCCO->CCCCC": 20.0,
    "CCCCC->Cc1ccccc1": 25.0,
    "CCCCO->CCCCCC": 22.0,
    "CCCCCC->Cc1ccccc1": 24.0,
}
_SYNTHETIC_INITIAL = {"CCCCO": 1.0}


def _channel_flux(steps):
    result = Recommender().predict_from_steps(
        list(steps), dict(_SYNTHETIC_BARRIERS), dict(_SYNTHETIC_INITIAL)
    )
    return result["debug_channel_flux"]["Cc1ccccc1"]


def test_two_independent_routes_to_one_product_add():
    """THE PROPERTY WAVE S1 EXISTS TO CREATE.

    Before 2026-08-27 the relaxation kept the lowest-span route per product and discarded
    every other route's flux, so this assertion read `combined == max(alone_1, alone_2)`.
    Real kinetics sums parallel channels; the allocation layer then partitions a FIXED
    volatile budget over those sums, so the sum changes the SPLIT and never the total
    (pinned separately in `test_the_volatile_budget_still_caps_the_sum`).
    """
    alone_1 = _channel_flux(_ROUTE_1)
    alone_2 = _channel_flux(_ROUTE_2)
    combined = _channel_flux(_ROUTE_1 + _ROUTE_2)

    assert len(alone_1) == 1
    assert len(alone_2) == 1
    assert len(combined) == 2, "the two routes must be recognised as two channels"

    expected = sum(alone_1.values()) + sum(alone_2.values())
    assert sum(combined.values()) == pytest.approx(expected, rel=1e-12), (
        "parallel channels must ADD. If this fails the propagator has gone back to "
        "winner-takes-all selection."
    )
    # And it is genuinely a sum, not a coincidence of equal fluxes: the two routes carry
    # measurably different flux (3.01x apart here), so `max` and `sum` cannot agree.
    ratio = sum(alone_2.values()) / sum(alone_1.values())
    assert ratio == pytest.approx(3.014, rel=1e-3)


def test_the_same_route_is_not_counted_twice():
    """The double-counting guard: channel identity is the route's FULL ORDERED STEP-SET.

    A duplicated step is the same channel and contributes ONCE. This is what stops an
    enumeration that emits the same transformation twice from doubling a prediction.
    """
    once = _channel_flux(_ROUTE_1)
    twice = _channel_flux(_ROUTE_1 + list(_ROUTE_1))

    assert len(twice) == 1, "a duplicated route is one channel, not two"
    assert sum(twice.values()) == pytest.approx(sum(once.values()), rel=1e-12)


def test_channel_identity_is_the_step_set_not_the_rate_limiting_step():
    """WHY THE OBVIOUS DEDUPE RULE WAS REJECTED — and it was rejected on a MEASUREMENT.

    The natural guard is "two routes sharing a rate-limiting step are one channel, take
    the max". It was implemented and measured. Both MFT routes on
    `cys_ribose_140C_Hofmann1998` — the pentodiulose lane and the Hofmann & Schieberle
    C2+C3 lane — have the SAME highest-barrier step, `Amadori_Rearrangement` at 29.06
    kcal/mol, because it sits on the shared cysteine/ribose trunk that every route in the
    network passes through. Under that rule MFT kept its old 242.38 ppb exactly and the
    whole live panel moved 3 rows by less than 1.15x: winner-takes-all in all but name.

    It is also wrong on the physics. For X --(slow)--> Y then two fast branches Y --> P and
    one fast branch Y --> Q, the trunk fixes the total flux and the branches PARTITION it
    by conductance, so P's share is 2/3, not 1/2. This propagator's per-route weight is
    pool * exp(-span/RT) with exp(span/RT) = sum_i exp(barrier_i/RT), so a dominant trunk
    collapses every route's weight onto the same trunk value — and it is precisely SUMMING
    them that reproduces the 2/3. Taking the max returns 1/2.
    """
    assert _route_channel_id([]) == ""
    a = _route_channel_id([{"step_key": "x->y", "barrier": 30.0}, {"step_key": "y->p", "barrier": 10.0}])
    b = _route_channel_id([{"step_key": "x->y", "barrier": 30.0}, {"step_key": "y->q", "barrier": 10.0}])
    assert a != b, (
        "routes that share their rate-limiting step but end differently are DISTINCT "
        "channels under the shipped rule"
    )
    assert a == "x->y|y->p"


def test_the_hofmann_lane_rate_limiting_steps_are_measured_not_assumed():
    """The measurement behind the rule choice, pinned so nobody re-derives the wrong expectation.

    Wave P predicted MFT would land at 242.38 + 71.02 = 313.39 ppb "if the two MFT channels
    are genuinely independent". Their rate-limiting steps are NOT distinct — both are the
    Amadori rearrangement — so that arithmetic never applied. The shipped, budget-normalised
    answer is 283.59 ppb, below Wave P's estimate, because the second MFT channel competes
    for the same fixed volatile budget as the FFT channels that also grew.

    RENAMED 2026-08-27 (Wave S2c) from `test_both_hofmann_mft_routes_share_the_amadori_
    rate_limiting_step`, because after the `thiol_addition_pentodiulose` 26.35 -> 28.60 revert
    THEY NO LONGER SHARE IT: the MFT lane now bottlenecks on `Thiol_Addition_Pentodiulose`
    (effective 29.5354 kcal/mol at aw 0.98) while the FFT lane still bottlenecks on
    `Amadori_Rearrangement` (29.0603). The old name asserted a conclusion; the new one names
    what the test actually owns, which is that these steps are MEASURED off the shipped tree
    rather than assumed. Full numbers and consequences in the comment at the assertions.
    """
    result = None
    original = Recommender.predict_from_steps

    def spy(self, *args, **kwargs):
        nonlocal result
        result = original(self, *args, **kwargs)
        return result

    Recommender.predict_from_steps = spy
    try:
        evaluate_benchmark(_HOFMANN)
    finally:
        Recommender.predict_from_steps = original

    assert result is not None
    channels = result["debug_channel_flux"]["Cc1occc1S"]
    # RE-PINNED 2026-08-28 (Wave X): 2 -> 3. The norfuraneol -> MFT step returned as a THIRD
    # parallel channel (`Furanone_Reductive_Sulfhydrylation`), constrained by Hofmann &
    # Schieberle 1998 Table 4 and gated by the isotope share test in
    # tests/scientific/test_wave_x_step_level_2026_08.py. This is precisely the situation
    # this test exists to police -- more parallel routes to one product - and the additive
    # rule is unchanged: the channel id is still the route's FULL ordered step-set, so the
    # third route is distinct and its flux is SUMMED, not selected. The budget invariant in
    # test_the_volatile_budget_still_caps_the_sum is what stops that sum minting mass.
    assert len(channels) == 3, "MFT is reached by exactly three enumerated routes here"

    def slowest_family(canon: str) -> str:
        path = result["debug_paths"][canon]
        return max(path, key=lambda t: float(t["barrier"]))["family"]

    # RE-PINNED 2026-08-27 (Wave S2c -- THE BARRIER REVERT). THE MFT LANE NO LONGER SHARES
    # ITS RATE-LIMITING STEP WITH THE FFT LANE, and that is a real structural movement, not a
    # bookkeeping one. `thiol_addition_pentodiulose` went 26.35 -> 28.60 kcal/mol; at this
    # benchmark's aw 0.98 the water-activity term adds ~0.935 kcal/mol to that step, so its
    # EFFECTIVE barrier is 29.5354 -- above the shared Amadori trunk's 29.0603. Measured
    # effective barriers on the shipped tree:
    #     MFT: Schiff 21.6070 | Amadori 29.0603 | Enolisation_2_3_Amadori 28.0311 |
    #          Deoxyosone_Reduction 28.9354 | Thiol_Addition_Pentodiulose 29.5354  <- slowest
    #     FFT: Schiff 21.6070 | Amadori 29.0603  <- slowest | Enolisation_Intermediate 21.0000 |
    #          Enolisation_1_2 27.6964 | Thiohemiacetal_Formation 23.3000 |
    #          Thiol_Dehydration 27.7354
    # WHAT THIS DOES AND DOES NOT CHANGE. It does NOT reopen the propagator rule: the shipped
    # channel id is the route's FULL ordered step-set (see the test above), so two routes are
    # distinct whether or not they share a rate-limiting step, and the additive sum is
    # unaffected. It DOES retire this test's original argument -- "Wave P's 242.38 + 71.02
    # arithmetic never applied because both lanes bottleneck on the same trunk". That
    # argument is now only half true, and the surviving half is the decisive one: the
    # volatile budget is FIXED, so single-lane figures were never addable regardless of where
    # the bottleneck sits. The budget invariant is asserted independently in
    # test_the_volatile_budget_still_caps_the_sum.
    assert slowest_family("Cc1occc1S") == "Thiol_Addition_Pentodiulose"
    assert slowest_family("SCc1ccco1") == "Amadori_Rearrangement"


def test_the_volatile_budget_still_caps_the_sum():
    """MASS HONESTY. Adding channels must REDISTRIBUTE the allocation, never mint mass.

    The allocation layer converts activities to mole fractions of a fixed
    `total_volatile_budget_molar`, so the molar sum is the budget both before and after —
    measured here by running the SAME converged state through the additive proxy and
    through the pre-Wave-S1 winner-takes-all proxy (still reachable via
    `channel_flux_totals=None`) and comparing both against the budget.
    """
    seen: list[dict[str, float]] = []
    original = recommend._project_weighted_flux_to_ppb

    def spy(*args, **kwargs):
        strategy = kwargs.get("projection_strategy", recommend.DEFAULT_PROJECTION_STRATEGY)
        budget = kwargs["projection_budget"]
        old = original(*args, **{**kwargs, "channel_flux_totals": None})
        new = original(*args, **kwargs)

        def molar(alloc):
            return sum(
                v / (recommend._mw_from_smiles(k) * strategy.ppb_conversion_factor)
                for k, v in alloc.items()
            )

        seen.append(
            {
                "budget": float(budget.total_volatile_budget_molar),
                "molar_old": molar(old),
                "molar_new": molar(new),
            }
        )
        return new

    recommend._project_weighted_flux_to_ppb = spy
    try:
        evaluate_benchmark(_HOFMANN)
    finally:
        recommend._project_weighted_flux_to_ppb = original

    assert seen, "the projection layer must have run"
    for row in seen:
        assert row["molar_old"] == pytest.approx(row["budget"], rel=1e-12)
        assert row["molar_new"] == pytest.approx(row["budget"], rel=1e-12), (
            "the additive propagator allocated a different total molar mass than the "
            "budget. Adding channels must move the split, not the total."
        )


def test_hofmann_sulfur_pair_after_the_additive_propagator():
    """OLD -> NEW on the flagship benchmark, and HALF OF IT GOT WORSE.

    MFT  242.38 -> 283.59 ppb vs 342   |  1.4110x under -> 1.2060x under   (better)
    FFT  218.00 -> 297.28 ppb vs 200   |  1.0900x over  -> 1.4864x over    (WORSE)

    Both compounds gained a second channel, so both rose; FFT was already over. The
    benchmark's own contract (1.45x / 0.09 dex) is UNTOUCHED and now fails on max_ratio as
    well as on MALE. Nothing was tuned to claw this back — the two lanes share their
    upstream trunk, so any barrier that pushed FFT down would push MFT down with it.

    RE-PINNED 2026-08-27 (Wave S1b -- THE pH ROUTING REPAIR, a DIFFERENT change):
    MFT  283.59 -> 154.85 ppb vs 342   |  1.2060x under -> 2.2086x under   (WORSE)
    FFT  297.28 -> 267.50 ppb vs 200   |  1.4864x over  -> 1.3375x over    (better)

    That movement is NOT the propagator. `get_ph_multiplier` -- the enolisation
    route-selection term, which had never been called on the prediction path -- now applies,
    and at this benchmark's pH 5.0 it gives the 1,2-enolisation (FFT / furfural) arm a 4.5x
    acid boost that the 2,3-enolisation (MFT) arm does not get. Against a fixed volatile
    budget that moves share from MFT to FFT. NO BARRIER MOVED and nothing was refitted; the
    untouched contract now fails both criteria by more (max_ratio 2.2086, MALE 0.2352 dex).
    The ADDITIVE-PROPAGATOR property this test was written to own is unaffected and is
    asserted separately: both compounds are still reached by two channels each.

    RE-PINNED AGAIN 2026-08-27 (Wave S2c -- THE ANCHOR RETIREMENT, a THIRD kind of change):
    MFT  154.85 ->  78.09 ppb vs 342   |  2.2086x under -> 4.3797x under   (WORSE)
    FFT  267.50 -> 293.67 ppb vs 200   |  1.3375x over  -> 1.4684x over    (WORSE)

    This one is not a propagator change and not a routing change -- it is a PROVENANCE
    finding with a mechanical consequence. `thiol_addition_pentodiulose` is REVERTED
    26.35 -> 28.60 (the un-fitted Wave N class value) and the Wave P fit record is RETRACTED,
    because that refit's sole fit target was this benchmark and Wave S2b showed its
    342 / 200 ppb are a REPO-INTERNAL DERIVATION: interior points of two invented, overlapping
    mol % bands in data/benchmarks/maillard_validation_benchmarks.md section 1.3, committed in
    the same commit as the benchmark JSON. The benchmark's 1.45x / 0.09 dex contract is
    RETIRED (not widened) and its tier demoted PRIMARY -> REFERENCE; it now inherits the
    global free-precursor default 1.5x / 0.10 dex and fails that by more (4.3797 / 0.4041 dex).
    BOTH ROWS ARE WORSE AND BOTH ARE PINNED WORSE. Read them as errors against a yardstick
    this repository invented, not as evidence about the chemistry.
    THE SULFUR BRANCH NOW HAS ZERO ABSOLUTE LITERATURE ANCHORS.
    """
    predicted = _predicted(_HOFMANN)
    # RE-PINNED 2026-08-28 (Wave X -- A THIRD MFT CHANNEL, NO BARRIER MOVED).
    # MFT 78.0867 -> 119.0800, FFT 293.6735 -> 282.8675. CAUSE: the norfuraneol -> MFT step
    # returned as a parallel channel (`Furanone_Reductive_Sulfhydrylation`), constrained by
    # Hofmann & Schieberle 1998 Table 4 and held to a minority in-situ share by
    # tests/scientific/test_wave_x_step_level_2026_08.py. The Wave X barrier fit was REJECTED
    # by that gate, so `furanone_reductive_sulfhydrylation` ships at the un-fitted
    # `thiol_addition` class value 28.60 and no constant in the tree moved.
    # BOTH DIRECTIONS ARE PINNED HONESTLY: MFT rose because a channel was ADDED and the
    # propagator sums channels; FFT fell because the volatile budget is FIXED and MFT's
    # share of it grew. THE YARDSTICK IS STILL NOT A MEASUREMENT: 342 / 200 ppb are this
    # repository's own arithmetic (Wave S2b, confirmed from the primary source by Wave W).
    # For the movement against a REAL measurement, see
    # hofmann1998_ribose_cysteine_145C_20min_pH5: MFT 468.58 -> 713.74 ppb against 198,
    # i.e. 2.37x over -> 3.61x over. This wave made that row WORSE and reports it as such.
    assert predicted["2-methyl-3-furanthiol"] == pytest.approx(119.0800, rel=1e-4)
    assert predicted["2-furfurylthiol"] == pytest.approx(282.8675, rel=1e-4)


# ── FIX 2: the matrix calibration registry is reachable again ────────────────────────

def _patched_records(protein_type: str, process_state: str, compound: str, scale: float):
    patched = []
    for record in registry._MATRIX_CALIBRATION_RECORDS:
        if (
            record.protein_type == protein_type
            and record.process_state == process_state
            and record.compound == compound
        ):
            record = dataclasses.replace(
                record, observable_factor=record.observable_factor * scale
            )
        patched.append(record)
    return tuple(patched)


def test_the_augmented_lane_now_sees_the_matrix_calibration_registry():
    """WAVE O's OWN MEASUREMENT, INVERTED.

    Wave O proved the defect by changing the soy hexanal observability factor by 4.32x and
    finding the `soy_isolate_ribose_cysteine_100C_45min_Internal2026` snapshot BIT-IDENTICAL.
    The same perturbation must now move the prediction by exactly that factor.

    NOTE WHICH RECORD IS PERTURBED, and it is the second half of Wave O's finding. The
    snapshot's conditions (100 C / 45 min / aw 0.95) make
    `determine_matrix_process_state` return `intermediate_matrix`, whose fallback order is
    (`intermediate_matrix`, `ambient_slurry`) — so the record actually reached is the
    AMBIENT one, not the heated_matrix entry, and the declared
    `aqueous_pre_extrusion_model` state in the benchmark file is not what the runtime
    recomputes. That process-state mismatch is UNFIXED and still open; this test pins the
    lookup, not the state.
    """
    baseline = _predicted(_SOY_SNAPSHOT)["Hexanal"]

    original = registry._MATRIX_CALIBRATION_RECORDS
    try:
        registry._MATRIX_CALIBRATION_RECORDS = _patched_records(
            "soy_iso", "ambient_slurry", "hexanal", 4.32
        )
        perturbed = _predicted(_SOY_SNAPSHOT)["Hexanal"]
    finally:
        registry._MATRIX_CALIBRATION_RECORDS = original

    assert perturbed == pytest.approx(baseline * 4.32, rel=1e-9), (
        "the compound-specific matrix calibration registry is unreachable from the "
        "matrix_precursor_augmented lane again — the SMILES-vs-name lookup regressed."
    )
    # And the restore worked, so the rest of the suite sees the shipped constants.
    assert _predicted(_SOY_SNAPSHOT)["Hexanal"] == pytest.approx(baseline, rel=1e-12)


def test_smiles_labelled_species_resolve_to_their_registry_name():
    """The fix itself, at the boundary. `_registry_compound_name` keeps a real label and
    falls back to the compound-database name only when the label IS the SMILES."""
    target_lookup = {"CCCCCC=O": {"name": "Hexanal", "type": "competing"}}
    assert recommend._registry_compound_name("CCCCCC=O", "CCCCCC=O", target_lookup) == "Hexanal"
    assert recommend._registry_compound_name("CCCCCC=O", None, target_lookup) == "Hexanal"
    assert recommend._registry_compound_name("CCCCCC=O", "hexanal", target_lookup) == "hexanal"
    # Unknown species degrade to the pre-fix behaviour rather than raising.
    assert recommend._registry_compound_name("CCCCCC=O", "CCCCCC=O", {}) == "CCCCCC=O"


def test_registry_reachability_moves_exactly_the_four_internal_snapshots():
    """WHICH ROWS FIX 2 MOVED, measured 2026-08-27, ordered by size.

    pea  Hexanal x4.317250  (the Wave O refitted pea ambient hexanal factor)
    soy  Hexanal x9.540070  (the Wave O refitted soy ambient hexanal factor)
    soy  Nonanal x1.066667  (= 0.160 / 0.150, the soy-vs-pea ambient release ratio)
    pea  Nonanal unchanged  (the pea ambient nonanal factor is 1.0)

    Every other scored row in the panel is bit-identical: the two Pratap-Singh and the
    three Trikusuma rows run the `matrix_only` path, which passed compound NAMES to the
    registry and was never affected, and the eight external hold-out points likewise.
    """
    pea = _predicted(Path("data/benchmarks/pea_isolate_ribose_cysteine_100C_45min_Internal2026.json"))
    soy = _predicted(_SOY_SNAPSHOT)
    assert pea["Hexanal"] == pytest.approx(0.7425331, rel=1e-5)
    assert soy["Hexanal"] == pytest.approx(1.7006158, rel=1e-5)
    assert soy["Nonanal"] == pytest.approx(0.0425151, rel=1e-5)
    assert pea["Nonanal"] == pytest.approx(0.0389915, rel=1e-5)
