"""Wave S1b (2026-08-27) — the pH / water-activity ROUTING repair.

Wave S2 built this repository's first ORDINAL accuracy measurement
(`docs/validation/directional_claims_panel.yml`) and found the model at 16/19 on sugar
identity, temperature, time and precursor loading, and **2/10 on pH and water activity** —
worse than a coin, and systematically so. The cause was not modelling. It was three
ROUTING defects, each confirmed by inspection before it was touched:

 1. `ReactionConditions.get_ph_multiplier` — documented, unit-tested for the correct signs,
    and NEVER CALLED on the prediction path.
 2. `_ionization_correction`'s pyrazine branch keyed on the substring ``"pyrazine"``, which
    matches NONE of the families this engine emits, so it returned 1.0 at every pH.
 3. `_water_activity_correction` reached 3 of 29 emitted families, none on the furan/HMF
    track, and its dehydration branch keyed on the equally dead ``"furfural"``.

NO CORRECTION CURVE WAS RESHAPED AND NO CONSTANT WAS REFITTED. These tests pin the
ROUTING — which family gets which term, and how many times — not the shapes, which
`tests/unit/test_conditions.py` and `tests/unit/test_wave_h_2026_08.py` already own.

The substring matching is the root cause and it is a repeat offender in this repo
(`Furanone_Strecker_Reduction` in Wave P; the Wave I offset keys), so the strongest guard
here is `test_water_activity_membership_matches_measured_stoichiometry`, which re-derives
set membership from the enumerated steps rather than trusting the literal in the source.
"""

from __future__ import annotations

import collections
import json
from pathlib import Path

import pytest

import src.conditions as conditions_module
from src.benchmark_validation import (
    benchmark_to_conditions,
    benchmark_to_formulation,
    evaluate_benchmark_payload,
)
from src.chem_utils import canonicalize_smiles
from src.conditions import ReactionConditions
from src.precursor_resolver import resolve_many
from src.smirks_engine import SmirksEngine

ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = ROOT / "data" / "benchmarks"
WATER = canonicalize_smiles("O", fallback_to_original=True)


#: Precursor pools the SHIPPED BENCHMARK PANEL does not contain, enumerated here so
#: the census below measures the ENGINE rather than the panel.
#:
#: WAVE T4 (2026-08-27) — THIS LIST IS THE REPAIR. Wave S1b's census iterated
#: `data/benchmarks/*.json` and reported its result as a property of the engine.
#: NOT ONE SHIPPED BENCHMARK USES A KETOSE — `grep -io fructose` over the sugar
#: lists of `data/benchmarks/*.json` returns 0 — so the ketose branch of
#: `_amadori_pathway` (`src/reaction_templates.py:58-76`) was never enumerated,
#: and `Heyns_Rearrangement` was invisible to every guard in the repo. Wave S1b
#: consequently wrote down "no Heyns family is emitted" in `src/conditions.py`
#: and asserted it here, while `tests/unit/test_smirks_engine.py::
#: TestHeynsRearrangement::test_heyns_fires` had been asserting the opposite in
#: the same suite the whole time. A true statement about the panel, published as
#: a statement about the engine.
#:
#: ADD A POOL HERE whenever the engine gains a branch the panel cannot reach.
_ENGINE_ONLY_PRECURSOR_POOLS: tuple[tuple[str, ...], ...] = (
    ("D-Fructose", "Glycine"),  # the ketose branch -> Heyns_Rearrangement
)


def _emitted_family_net_water() -> dict[str, set[int]]:
    """Net water produced per step, per family, over every family the ENGINE can emit.

    This is the MEASUREMENT the water-activity family sets are derived from. It is
    recomputed here rather than imported, so the test cannot agree with the source by
    construction.

    SCOPE, and it is load-bearing (Wave T4): this enumerates the shipped benchmark
    panel PLUS `_ENGINE_ONLY_PRECURSOR_POOLS`. Restricting it to the panel, as Wave
    S1b did, silently scopes every claim built on it to the panel — see the note on
    that constant.
    """
    net: dict[str, set[int]] = collections.defaultdict(set)
    step_sets: list[list] = []
    for path in sorted(BENCHMARK_DIR.glob("*.json")):
        try:
            bench = json.loads(path.read_text(encoding="utf-8"))
            formulation = benchmark_to_formulation(bench)
            cond = benchmark_to_conditions(bench)
            names = (
                formulation["sugars"]
                + formulation["amino_acids"]
                + formulation.get("additives", [])
                + formulation.get("lipids", [])
            )
            step_sets.append(SmirksEngine(cond).enumerate(resolve_many(names), max_generations=4))
        except Exception:  # noqa: BLE001 — unsupported payloads are simply not evidence here
            continue
    for pool in _ENGINE_ONLY_PRECURSOR_POOLS:
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0, water_activity=0.6)
        step_sets.append(SmirksEngine(cond).enumerate(resolve_many(list(pool)), max_generations=4))
    for steps in step_sets:
        for step in steps:
            produced = sum(
                1
                for s in step.products
                if canonicalize_smiles(s.smiles, fallback_to_original=True) == WATER
            )
            consumed = sum(
                1
                for s in step.reactants
                if canonicalize_smiles(s.smiles, fallback_to_original=True) == WATER
            )
            net[(step.reaction_family or "unknown").lower()].add(produced - consumed)
    return dict(net)


# ── defect 1: the wiring, and the fact that it can only apply once ────────────


def test_get_ph_multiplier_is_reachable_from_the_prediction_path():
    """The whole point of Wave S1b: this function is no longer dead.

    Wave S2's finding, verbatim: "Written, tested, and never executed where it matters."
    `evaluate_benchmark_payload` reaches `conditions.get_rate_constant()`, which now routes
    the three enolisation branch-point families through `get_ph_multiplier`. If this ever
    goes back to 1.0 the pH physics has been disconnected again.
    """
    acidic = ReactionConditions(pH=4.0, water_activity=0.8)
    alkaline = ReactionConditions(pH=9.0, water_activity=0.8)

    # 1,2-enolisation opens 3-deoxyosone -> furfural / HMF and is ACID-catalysed.
    assert acidic._enolisation_route_ph_correction("Enolisation_1_2") > 4.0
    assert alkaline._enolisation_route_ph_correction("Enolisation_1_2") == pytest.approx(
        1.0, abs=0.05
    )
    # 2,3-enolisation opens 1-deoxyosone -> reductone -> pyrazine / furanone, BASE-catalysed.
    assert acidic._enolisation_route_ph_correction("Enolisation_2_3") < 0.5
    assert alkaline._enolisation_route_ph_correction("Enolisation_2_3") > 4.0

    # And the rate constant itself must carry it, which is what "on the prediction path"
    # means. Same barrier, same temperature, same aw: only pH differs.
    k_acid = acidic.get_rate_constant("Enolisation_1_2", ea_override_kcal=25.0)
    k_alkaline = alkaline.get_rate_constant("Enolisation_1_2", ea_override_kcal=25.0)
    assert k_acid > 4.0 * k_alkaline, (
        "get_rate_constant is not applying the enolisation route-selection term — this is "
        "exactly the Wave S2 defect, reintroduced"
    )


def test_the_gate_admits_only_the_enolisation_branch_point():
    """`get_ph_multiplier` matches by SUBSTRING internally, so the gate is load-bearing.

    Its substrings reach roughly half the emitted families. Two consequences the gate
    exists to prevent, both stated as assertions:

      * ``"condensation"`` matches `Aminoketone_Condensation` — the PYRAZINE step — and
        would hand it the ACID-peaked Schiff Gaussian, i.e. Wave S2's defect 2 again in the
        opposite direction.
      * ``"thiol"`` / ``"thio"`` / ``"furan"`` reach the whole downstream sulfur and furan
        track, so an ungated call would compound one physical effect (the acid preference of
        the branch point) five or six times along a single route.
    """
    cond = ReactionConditions(pH=4.0, water_activity=0.8)

    # The raw function DOES claim these, which is why it must not be called directly here.
    assert cond.get_ph_multiplier("Aminoketone_Condensation") > 1.0
    assert cond.get_ph_multiplier("Thiol_Oxidation") > 4.0
    assert cond.get_ph_multiplier("Furanone_Cyclisation") > 4.0

    # The gate refuses all of them.
    for family in (
        "Aminoketone_Condensation",
        "Thiol_Oxidation",
        "Thiol_Dehydration",
        "Thiohemiacetal_Formation",
        "Furanone_Cyclisation",
        "Cysteine_Degradation",
        "Schiff_Base_Formation",
        "Strecker_Degradation",
        "Amadori_Rearrangement",
        "Enolisation_Intermediate",
    ):
        assert cond._enolisation_route_ph_correction(family) == 1.0, (
            f"{family} is receiving the enolisation route-selection term. Only the three "
            f"branch-point families may."
        )


def test_no_family_can_receive_the_ph_physics_twice():
    """The two pH terms are different physics and must never both apply to one family.

    `_ionization_correction` is REAGENT AVAILABILITY (Henderson-Hasselbalch free-base
    fraction). `get_ph_multiplier` is ROUTE SELECTION (which way the Amadori compound
    enolises). Applying both to one family would count a single pH dependence twice.
    `src/conditions.py` asserts the disjointness at import; this asserts it behaviourally,
    over every family the engine can emit.
    """
    assert not (
        conditions_module._ENOLISATION_ROUTE_PH_FAMILIES
        & conditions_module._ALPHA_AMINO_NUCLEOPHILE_FAMILIES
    )
    assert not (
        conditions_module._ENOLISATION_ROUTE_PH_FAMILIES
        & conditions_module._AMINOKETONE_NUCLEOPHILE_FAMILIES
    )
    assert not (
        conditions_module._ALPHA_AMINO_NUCLEOPHILE_FAMILIES
        & conditions_module._AMINOKETONE_NUCLEOPHILE_FAMILIES
    )

    cond = ReactionConditions(pH=4.0, water_activity=0.8)
    for family in _emitted_family_net_water():
        ionisation = cond._ionization_correction(family)
        route = cond._enolisation_route_ph_correction(family)
        assert ionisation == 1.0 or route == 1.0, (
            f"{family} receives BOTH pH terms (ionisation={ionisation}, route={route}). "
            f"One physical dependence would be counted twice."
        )


# ── defect 2: the dead pyrazine key ───────────────────────────────────────────


def test_the_pyrazine_ionisation_branch_reaches_the_family_that_makes_pyrazines():
    """Wave S2: the ``"pyrazine"`` substring matched NO emitted family.

    The engine's pyrazine step is `Aminoketone_Condensation` (two aminoacetones ->
    2,5-dimethylpyrazine + 2 H2O + H2, `reaction_templates.py:1000`). The pKa 6.5 is the
    value the dead branch already carried and was NOT retuned.
    """
    emitted = set(_emitted_family_net_water())
    assert not any("pyrazine" in fam for fam in emitted), (
        "a family name now contains 'pyrazine' — re-read the substring history in "
        "src/conditions.py before assuming the old key would have worked"
    )
    assert "aminoketone_condensation" in emitted

    low = ReactionConditions(pH=4.0, water_activity=0.8)
    high = ReactionConditions(pH=9.0, water_activity=0.8)
    assert low._ionization_correction("Aminoketone_Condensation") < 0.05
    assert high._ionization_correction("Aminoketone_Condensation") > 0.9
    # pKa 6.5 exactly: half-ionised at pH 6.5.
    at_pka = ReactionConditions(pH=6.5, water_activity=0.8)
    assert at_pka._ionization_correction("Aminoketone_Condensation") == pytest.approx(0.5)


def test_dimethylpyrazine_now_rises_with_ph_as_two_independent_measurements_say():
    """The two directional-panel rows this fix was diagnosed against (PH-04, PH-06).

    Laemont & Barringer 2023 measure dimethylpyrazine isomers at 26.6 -> 37.4 -> 68.2 ppb
    over pH 4 -> 7 -> 9. Before Wave S1b the model gave 99.5 -> 16.8 -> 14.9, i.e. the wrong
    way, in every system tested. THE MAGNITUDES ARE NOT PINNED AND ARE NOT CLAIMED TO BE
    RIGHT — only the ORDINAL direction, which is the only thing the panel scores.
    """

    def dmp(ph: float) -> float:
        payload = {
            "benchmark_id": f"wave_s1b_dmp_ph{ph}",
            "precursors": {
                "Glycine": {"concentration_mM": 10.0},
                "D-Glucose": {"concentration_mM": 10.0},
            },
            "conditions": {
                "temp_C": 160.0,
                "ph": ph,
                "water_activity": 0.9,
                "time_min": 30.0,
            },
            "metadata": {"tier": "PROBE", "family": "wave_s1b_probe"},
            "measured_volatiles": {},
            "protein_type": "free",
        }
        evaluation = evaluate_benchmark_payload(payload, benchmark_id=payload["benchmark_id"])
        lowered = {k.lower(): v for k, v in evaluation.predicted_ppb.items()}
        return float(lowered["2,5-dimethylpyrazine"])

    acidic, neutral, alkaline = dmp(4.0), dmp(7.0), dmp(9.0)
    assert acidic < neutral < alkaline, (
        f"2,5-dimethylpyrazine vs pH is {acidic:.4g} / {neutral:.4g} / {alkaline:.4g} at "
        f"pH 4 / 7 / 9. It must RISE with pH — that is the direction two independent "
        f"measurements report, and getting it backwards was Wave S2's headline pH failure."
    )


# ── defect 3: the water-activity reach ────────────────────────────────────────


def test_water_activity_membership_matches_measured_stoichiometry():
    """THE STRONGEST GUARD IN THIS MODULE, and the reason it re-measures rather than imports.

    Membership in the water-activity sets is not a keyword and not a judgement call: a step
    that RELEASES water is pushed back by water, a step that CONSUMES water is rate-limited
    by it, and a net-zero step gets no mass-action term at all. This test recomputes the net
    water per step for every family the engine emits and checks the shipped sets against it.

    A family whose stoichiometry changes must have its membership RE-DERIVED. Do not edit
    the sets to make this pass.

    ONE DELIBERATE EXCLUSION, asserted rather than ignored: `Additive_Thermal_Degradation`
    is the only stoichiometrically NON-UNIFORM family (+2 / 0 / -1 / -2 across its steps),
    so a single family-level factor cannot honestly represent it and it is in neither set.
    """
    net = _emitted_family_net_water()
    releasing = conditions_module._WATER_RELEASING_FAMILIES
    consuming = conditions_module._WATER_CONSUMING_FAMILIES

    non_uniform = {fam for fam, vals in net.items() if len(vals) > 1}
    assert non_uniform == {"additive_thermal_degradation"}, (
        f"the set of stoichiometrically non-uniform families changed to {non_uniform}; "
        f"each new member needs its own decision, not a default"
    )

    mismatches = []
    for family, values in sorted(net.items()):
        in_releasing = family in releasing
        in_consuming = family in consuming
        if family in non_uniform:
            expected = not (in_releasing or in_consuming)
            want = "in neither set (non-uniform)"
        else:
            net_water = next(iter(values))
            if net_water > 0:
                expected, want = in_releasing and not in_consuming, "water-releasing"
            elif net_water < 0:
                expected, want = in_consuming and not in_releasing, "water-consuming"
            else:
                expected, want = not (in_releasing or in_consuming), "in neither set"
        if not expected:
            mismatches.append(
                f"{family}: net water {sorted(values)}, should be {want}, "
                f"releasing={in_releasing} consuming={in_consuming}"
            )
    assert not mismatches, "water-activity set membership contradicts the stoichiometry:\n" + "\n".join(
        mismatches
    )


def test_the_water_activity_correction_reaches_the_furan_track():
    """Wave S2's finding was that it "misses the furan track entirely". This is the fix.

    `Enolisation_1_2` IS the furan track: 3-deoxyosone -> furfural / HMF + 2 H2O. It is also
    the family a bug inside Wave S1b itself briefly excluded — the Wave H
    `_releases_rather_than_attacks_with_the_amine` substring guard matches "enolisation" and
    was returning 1.0 here before any set was consulted. This test is that bug's tombstone.
    """
    dry = ReactionConditions(pH=6.0, water_activity=0.30)
    wet = ReactionConditions(pH=6.0, water_activity=0.95)
    assert dry._water_activity_correction("Enolisation_1_2") == pytest.approx(1.0)
    assert wet._water_activity_correction("Enolisation_1_2") == pytest.approx(0.35)

    # Wave H's protected case survives, and survives for a BETTER reason than before: it is
    # net-zero in water and therefore in no set, rather than being rescued by a substring.
    assert wet._water_activity_correction("Enolisation_2_3_Amadori") == 1.0
    assert wet._ionization_correction("Enolisation_2_3_Amadori") == 1.0


def test_hydrolytic_families_are_promoted_by_water_not_inhibited_by_it():
    """The "hydrolyses the reverse" arm. Mass action, no free parameter.

    Two emitted families consume water (`Cysteine_Degradation`, the hydrolytic release of
    H2S; `Generalized_Deamination`). Their rate tracks water availability, so the factor is
    `aw` itself, floored the way the dehydration arm is floored.
    """
    dry = ReactionConditions(pH=6.0, water_activity=0.30)
    wet = ReactionConditions(pH=6.0, water_activity=0.95)
    for family in ("Cysteine_Degradation", "Generalized_Deamination"):
        assert dry._water_activity_correction(family) == pytest.approx(0.30)
        assert wet._water_activity_correction(family) == pytest.approx(0.95)
        assert wet._water_activity_correction(family) > dry._water_activity_correction(family)


# ── the dead-key census ───────────────────────────────────────────────────────


def test_the_substring_keys_that_matched_nothing_stay_documented_as_dead():
    """Seven keys were called dead in Wave S1b's audit. This pins the corrected census.

    They are reported rather than silently deleted because the keys document intent — but
    the test asserts they are still dead, so if a family is ever renamed into one of them
    the surprise happens here rather than in a prediction.

    WAVE T4 (2026-08-27) — ONE OF THE SEVEN WAS NEVER DEAD. ``"heyns"`` matches
    `Heyns_Rearrangement`, which `src/reaction_templates.py:60` emits for every
    ketose + amino acid pair. It is excluded from the loop below and asserted LIVE
    instead. The old census passed because `_emitted_family_net_water` enumerated
    only `data/benchmarks/`, and no shipped benchmark uses a ketose; the enumerator
    now also runs `_ENGINE_ONLY_PRECURSOR_POOLS`. This is not a relaxation — the
    assertion below is strictly stronger than "no hits", because it pins WHICH
    family the key reaches and that it is the sole matcher.
    """
    emitted = set(_emitted_family_net_water())
    for key in ("pyrazine", "furfural", "nitrogen_heterocycle", "oxygen_heterocycle", "1,2", "2,3"):
        hits = sorted(fam for fam in emitted if key in fam)
        assert not hits, (
            f"substring {key!r} now matches {hits}. It was DEAD when Wave S1b audited it, "
            f"and several corrections were written as though it were live — re-read "
            f"src/conditions.py before relying on either."
        )

    # THE KEY THAT WAS LIVE ALL ALONG. `get_ph_multiplier`'s branch 2 lists both
    # "amadori" and "heyns"; "amadori" does NOT appear in the string
    # "heyns_rearrangement", so "heyns" is the ONLY matcher the ketose
    # rearrangement has in the three ungated lanes (kinetics, pathway_ranker,
    # cantera_export). Deleting it on the strength of the old census would have
    # silently removed their pH dependence.
    assert sorted(fam for fam in emitted if "heyns" in fam) == ["heyns_rearrangement"]
    assert not any("heyns" in fam for fam in emitted if "amadori" in fam), (
        "if a family name ever carries both, the two branch-2 keys stop being "
        "independent and the 'sole matcher' argument above needs re-deriving"
    )

    # And the one that is live and dangerous, which is why the gate exists.
    assert "aminoketone_condensation" in {fam for fam in emitted if "condensation" in fam}


def test_the_ketose_rearrangement_is_corrected_like_its_aldose_twin():
    """Wave T4's headline. `Heyns_Rearrangement` used to escape BOTH pH and aw terms.

    The Heyns rearrangement is the ketose analogue of Amadori (`barrier_constants.py`
    says so in the barrier note itself), carries the same alpha-amino nitrogen through
    the same Schiff-base -> 1,2-shift, and is net-zero in water exactly as Amadori is
    (measured: one emitted step, one reactant, one product, no water either side). It
    was nevertheless in NONE of Wave S1b's family sets while Amadori was in two, so on
    any ketose formulation the rearrangement ran uncorrected: at pH 5 the ratio of the
    two families' total correction was 1208x.

    This pins EQUALITY rather than a magnitude. If the aldose treatment is ever
    re-derived (it is an open [P]: for a rearrangement the amine is already bound, so
    the free-base fraction is a proxy, not a mechanism), both must move together.
    """
    for ph in (5.0, 6.0, 7.0, 9.0):
        for aw in (0.3, 0.8, 0.95):
            cond = ReactionConditions(pH=ph, water_activity=aw)
            for term in ("_ionization_correction", "_water_activity_correction",
                         "_enolisation_route_ph_correction"):
                amadori = getattr(cond, term)("Amadori_Rearrangement")
                heyns = getattr(cond, term)("Heyns_Rearrangement")
                assert heyns == pytest.approx(amadori), (
                    f"{term} disagrees between the aldose and ketose rearrangements at "
                    f"pH {ph} / aw {aw}: Amadori {amadori:.6g} vs Heyns {heyns:.6g}. These "
                    f"are the same reaction on the other sugar; whichever value is right, "
                    f"they cannot differ."
                )

    # And it is a REAL correction, not two matching 1.0s.
    acidic = ReactionConditions(pH=5.0, water_activity=0.8)
    assert acidic._ionization_correction("Heyns_Rearrangement") < 0.01
    assert acidic._water_activity_correction("Heyns_Rearrangement") < 1.0


def test_each_family_takes_exactly_one_water_activity_treatment():
    """The aw counterpart of `test_no_family_can_receive_the_ph_physics_twice`.

    WAVE T4 (2026-08-27). `_water_activity_correction` returns from the FIRST arm that
    matches — empirical, releasing, consuming — so an overlap does not double-multiply,
    it SILENTLY SHADOWS: a water-shedding family placed in the empirical set would get
    the peaked browning curve instead of its mass-action term, with no error anywhere.
    That failure is quieter than the pH one and had no assertion until now.
    `src/conditions.py` asserts the disjointness at import; this asserts it over every
    family the engine can emit.
    """
    empirical = conditions_module._LABUZA_EMPIRICAL_FAMILIES
    releasing = conditions_module._WATER_RELEASING_FAMILIES
    consuming = conditions_module._WATER_CONSUMING_FAMILIES
    assert not (empirical & releasing)
    assert not (empirical & consuming)
    assert not (releasing & consuming)

    for family in _emitted_family_net_water():
        memberships = [
            name
            for name, group in (
                ("empirical", empirical), ("releasing", releasing), ("consuming", consuming)
            )
            if family in group
        ]
        assert len(memberships) <= 1, (
            f"{family} is in {memberships}. `_water_activity_correction` would return "
            f"from {memberships[0]} and the rest would be dead — one water-activity "
            f"treatment per family."
        )
