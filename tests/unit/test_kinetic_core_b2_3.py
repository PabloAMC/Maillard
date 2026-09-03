"""FROZEN-WAVE REGRESSION RECORD (labelled 2026-09-03, test audit). The wave generator this file
tests is frozen (scripts/generators/WAVES.md); these tests fail only if the frozen report, the
network or the parameter tables change. They are the contract of a finished wave, not live checks
of new behaviour.


tests/unit/test_kinetic_core_b2_3.py

BUILD WAVE B2.3 -- the three things this wave must not be able to lose.

B2.3 is a CONSERVATION FIX, not a modelling wave: no parameter, no rate, no fit
row. Its whole deliverable is that a class of defect becomes unreachable, and a
class of data gap becomes visible. So the tests here are not accuracy tests --
they are tests that the guard rails are load-bearing:

  1. THE CHARGE-CLOSURE VALIDATOR, with a PLANTED VIOLATION control. A validator
     that has never been seen to fail is not evidence of anything, so the
     control below reintroduces exactly the B2.2 defect on a copy of the
     network and asserts that the import-time check rejects it.
  2. BUNDLE-HASH INVARIANCE. Amendment 9 clause 2 lets B2.3 complete the
     CONDITION records of the 21 frozen bundles and nothing else. The hashes
     below were taken from the bundles BEFORE the completion script first ran
     and are compared against the bundles with `conditions.buffer` deleted, so
     any edit to a measured value, a compound list, an evidence_class, a role
     or one of the other three condition fields fails this file.
  3. THE BUFFER-BLOCK SCHEMA, including the rule that "not on disk" must
     produce `buffer_unknown` and never a plausible guess.

Plus the firewall, in the two forms B2.2 established: a literal grep and a
SYSTEMS walk.

They are unit tests. Anything that reads the frozen fit report skips itself if
it has not been generated, so the file runs in seconds either way.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from src.kinetic_core import acrylamide as A
from src.kinetic_core import network as N
from src.kinetic_core import ph_state as PH
from src.kinetic_core import sulfur as S
from src.kinetic_core.species_sulfur import SULFUR_INDEX, TERMINAL_POOLS

ROOT = Path(__file__).resolve().parents[2]
FIT_REPORT = ROOT / "results/validation/kinetic_core_b2_3_fit_report.json"
BUNDLE_DIR = ROOT / "data/benchmarks/external_validation"


def _frozen():
    if not FIT_REPORT.exists():  # pragma: no cover - environment dependent
        pytest.skip(f"{FIT_REPORT.name} not generated yet")
    return json.loads(FIT_REPORT.read_text())


# ===========================================================================
# 1. THE CHARGE-CLOSURE VALIDATOR
# ===========================================================================


def test_the_whole_network_passes_charge_closure_at_import():
    """
    It already ran at import -- `sulfur.py` and `acrylamide.py` both call it at
    module scope, so an unbalanced network cannot even be imported. Calling it
    again here is not redundant: it is what makes the FAILURE of the planted
    control below meaningful, by showing the same function returns cleanly on
    the real network.
    """
    PH.validate_charge_closure(S.FULL_REACTIONS, S.CENTRE_LEDGER)
    PH.validate_charge_closure(N.REACTIONS, N.TRUNK_CENTRE_LEDGER)
    PH.validate_charge_closure(
        A.FULL_ACRYLAMIDE_REACTIONS, A.ACRYLAMIDE_CENTRE_LEDGER)


def test_planted_violation_the_b2_2_defect_is_rejected():
    """
    THE CONTROL. Reintroduce the exact B2.2 defect -- a sink that deletes
    cysteine's carboxylate into an untitratable pool without declaring it --
    and assert the validator refuses it.

    `r_cys_thermal` is used because it is the first site B2.2's diagnosis sec.
    3 names, and the planted version is byte-for-byte B2.2's own stoichiometry.
    """
    from src.kinetic_core.network import Reaction

    planted = Reaction(
        "r_cys_thermal", {"Cys": 1},
        {"FRAG_C": 3, "FRAG_N": 1, "FRAG_S": 1},  # <- B2.2, verbatim
        "k_cys_thermal", "planted violation, test only")
    network = tuple(
        planted if r.key == "r_cys_thermal" else r for r in S.FULL_REACTIONS)

    with pytest.raises(ValueError) as excinfo:
        PH.validate_charge_closure(network, S.CENTRE_LEDGER)
    message = str(excinfo.value)
    assert "r_cys_thermal" in message
    assert "carboxyl" in message
    # the ledger already declares this step's AMINE, so the defect surfaces as
    # a declaration that no longer matches its stoichiometry -- which is the
    # right diagnosis and names the exact discrepancy
    assert "-1" in message and "0" in message


def test_planted_violation_an_undeclared_leak_names_the_consequence():
    """
    THE OTHER HALF OF THE CONTROL: a step that leaks and is in NO ledger entry
    at all. This is the path a NEW reaction would take, and its message has to
    explain WHY it matters -- a validator whose failure reads as a schema
    complaint gets suppressed; one that says "the pot drifts for a bookkeeping
    reason" gets fixed.
    """
    from src.kinetic_core.network import Reaction

    planted = Reaction(
        "r_dpo_nf", {"DPO": 1, "Cys": 1}, {"NF": 1, "FRAG_C": 3, "FRAG_N": 1,
                                           "FRAG_S": 1},
        "k_dpo_nf", "planted violation, test only")
    network = tuple(
        planted if r.key == "r_dpo_nf" else r for r in S.FULL_REACTIONS)

    with pytest.raises(ValueError) as excinfo:
        PH.validate_charge_closure(network, S.CENTRE_LEDGER)
    message = str(excinfo.value)
    assert "r_dpo_nf" in message
    assert "NOT" in message and "CENTRE_LEDGER" in message
    assert "bookkeeping" in message
    assert "11.4" in message and "4.9" in message, (
        "the message must carry the measured consequence, not just the rule")


def test_planted_violation_a_ledger_that_lies_is_rejected():
    """
    The second failure mode: the stoichiometry and its declaration drift apart.
    A ledger nobody re-checks is how a fix rots.
    """
    lying = dict(S.CENTRE_LEDGER)
    lying["r_tdg_fa"] = {"carboxyl": +7, "basis": "deliberately wrong"}
    with pytest.raises(ValueError, match="r_tdg_fa"):
        PH.validate_charge_closure(S.FULL_REACTIONS, lying)


def test_planted_violation_a_stale_ledger_entry_is_rejected():
    """The third: paperwork for a step that moves nothing."""
    stale = dict(S.CENTRE_LEDGER)
    stale["r_dpo_nf"] = {"carboxyl": +1, "basis": "dead paperwork"}
    with pytest.raises(ValueError, match="r_dpo_nf"):
        PH.validate_charge_closure(S.FULL_REACTIONS, stale)


def test_every_declared_movement_carries_a_written_basis():
    for ledger in (S.CENTRE_LEDGER, A.ACRYLAMIDE_CENTRE_LEDGER,
                   N.TRUNK_CENTRE_LEDGER):
        for key, entry in ledger.items():
            assert str(entry.get("basis", "")).strip(), key
            assert len(str(entry["basis"])) > 40, (
                f"{key}: a one-word basis is not a basis")


def test_no_step_deletes_a_carboxyl_without_declaring_it():
    """
    THE HEADLINE PROPERTY, stated positively. After B2.3 no reaction in either
    extended network destroys a modelled CARBOXYL centre except the two trunk
    steps that decompose an acid the trunk itself made.
    """
    allowed = {"r_fa_frag", "r_aa_frag", "a_cys_glc"}
    for network in (S.FULL_REACTIONS, A.FULL_ACRYLAMIDE_REACTIONS):
        for reaction in network:
            delta = PH.centre_delta(reaction.reactants, reaction.products)
            if delta["carboxyl"] < 0:
                assert reaction.key in allowed, (
                    f"{reaction.key} deletes a carboxylate; that is the B2.2 "
                    f"defect and it is not on the allowed list")


def test_the_centre_table_matches_the_charge_balance_it_claims_to_describe():
    """
    `TITRATABLE_CENTRES` must not drift away from `titratable_inventory`. A
    species has a centre in the table if and only if putting one mmol/L of it
    into the inventory produces a group -- otherwise the audit is policing a
    fiction.
    """
    drift = PH.PhDrift(acid_yield=1.0, arp_amine_pka=8.0)
    spec = PH.BufferSpec(kind="none", declared=True, source="test")
    for key, centres in PH.TITRATABLE_CENTRES.items():
        if key not in SULFUR_INDEX:
            continue  # B3-lane species, checked by its own network's validator
        inventory = PH.titratable_inventory({key: 1.0}, drift, spec)
        assert inventory, f"{key} claims centres {centres} but adds no group"
    # and the converse: a species the inventory CAN see must be in the table
    for key in ("Cys", "ARP", "TTCA", "FA", "AA", "ACID", "CBX"):
        assert key in PH.TITRATABLE_CENTRES, key


def test_the_untracked_list_is_pinned_so_nothing_joins_silently():
    """
    The declared-gap list is enumerated, not open. If a new titratable species
    enters the network, this test is what forces somebody to decide whether it
    belongs in the balance or in the gap list.
    """
    assert set(PH.UNTRACKED_TITRATABLE) == {
        "H2S", "Gly", "SB", "AMA", "MEL_N", "THI", "MFT", "FFT", "MESH",
    }
    for key, reason in PH.UNTRACKED_TITRATABLE.items():
        assert len(reason) > 60, f"{key}: a gap needs a stated reason"


def test_the_carried_carboxyl_pool_is_terminal_and_moves_only_the_ph():
    """
    `CBX` may never be a reactant. If it were, adding it would have changed a
    volatile prediction, and the claim that a conservation fix cannot move a
    volatile except through the pH would be false.
    """
    assert "CBX" in TERMINAL_POOLS
    for reaction in S.FULL_REACTIONS + A.FULL_ACRYLAMIDE_REACTIONS:
        assert "CBX" not in reaction.reactants, reaction.key
    producers = [r.key for r in S.FULL_REACTIONS if "CBX" in r.products]
    assert len(producers) >= 9, producers


def test_the_carried_carboxyl_enters_at_full_strength_and_the_lump_does_not():
    """
    THE DISTINCTION THAT JUSTIFIES A SECOND POOL. `ACID` is a lump whose acid
    yield is fitted, so it enters the charge balance scaled. `CBX` is a
    carboxyl carried out of a molecule that had one, so it enters unscaled.
    Halving the fitted yield must move the ACID contribution and leave CBX
    alone.
    """
    spec = PH.BufferSpec(kind="none", declared=True, source="test")
    full = PH.PhDrift(acid_yield=1.0, arp_amine_pka=8.0)
    half = PH.PhDrift(acid_yield=0.5, arp_amine_pka=8.0)

    def total(state, drift):
        return sum(c for _g, c in PH.titratable_inventory(state, drift, spec))

    assert total({"ACID": 10.0}, half) == pytest.approx(
        0.5 * total({"ACID": 10.0}, full))
    assert total({"CBX": 10.0}, half) == pytest.approx(
        total({"CBX": 10.0}, full))


def test_the_amine_fate_is_declared_once_and_says_what_would_falsify_it():
    basis = PH.AMINE_FATE_BASIS
    assert "3.22" in basis and "3.42" in basis, (
        "the declaration must name the FIT anchors it rests on")
    assert "9.25" in basis, "it must name the pKa that makes the arithmetic bite"
    # every ledger entry that destroys an amine must cite it, except the one
    # whose evidence is its own product's ring
    for key, entry in S.CENTRE_LEDGER.items():
        if int(entry.get("amine", 0)) >= 0:
            continue
        if key == "r_cys_actz":
            assert "thiazol" in entry["basis"].lower()
            continue
        assert basis in entry["basis"], (
            f"{key} destroys an amine without citing AMINE_FATE_BASIS")


def test_carbon_nitrogen_and_sulfur_still_balance_everywhere():
    """The three older invariants, unchanged by the fourth."""
    S.validate_sulfur_balance()
    N.validate_balance()
    A.validate_acrylamide_balance()


# ===========================================================================
# 2. BUNDLE-HASH INVARIANCE
# ===========================================================================
# SHA-256 of each bundle with `conditions.buffer` deleted, taken from the
# bundles BEFORE scripts/generators/complete_benchmark_buffer_fields.py first
# ran. Hashing EVERYTHING-EXCEPT rather than just the measured block is
# deliberate: it proves not only that the measured values are untouched but
# that nothing else moved either.

BUNDLE_BASELINES = {
    "external_validation_bi_2020_raw_pea_hexanal":
        "b8f9814c0c2c10701679de5bf0a64aa9b393df2ad45d368f4e0940ade6a58e5d",
    "external_validation_bi_2020_roasted_pea_hexanal":
        "7d5cf5174622face448ff0095c03b5813d9227481b05af6eb4e8343b05b41ffe",
    # re-pinned 2026-09-03: quantification_class + quantification_note added (measured block untouched)
    "external_validation_li_2026_spi_wg_hme_control":
        "7af20770d5e2ebb52817a51db8e25f29577cf09994b7113545115d4a1b33a2d5",
    # re-pinned 2026-09-03: quantification_class + quantification_note added (measured block untouched)
    "external_validation_liu_2023_ppi_offnote_baseline":
        "5e0b0a7ec0e56a5efff9c1d01225a5eac15c9129fe58fd08032bd1226cd15dcd",
    "mp_holdout_fructose_asparagine_180C_Lin2022":
        "85efa67b7c88adc2d8616d7ed601efe6f6f7763238b2aef237bc485c405d0505",
    "mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019":
        "af74d73ebf7b3a0e5ab06a6a5979c1bedc631927b156322605576f4d5006dd67",
    "mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019":
        "539151b2537cfad852a7c6e77821afbd73059294ce701f4752250da7806da7d5",
    "mp_holdout_glucose_asparagine_180C_10min_Chang2021":
        "b2a0a4f3eef22450e8c4f2d29643e16b052555dc74eb73f0cedc297b9c1cd531",
    "mp_holdout_glucose_asparagine_180C_30min_Chang2021":
        "b9adc2ffe8be995795fc816cc01823b53a5bb4d64876b64377ed45da9b6db532",
    "mp_holdout_glucose_asparagine_180C_30min_water_Chang2021":
        "531bb130e00c5aafe789a2d2e47867924a894df0cc086eaa91a8ed7eb4c8bc3d",
    "mp_holdout_glucose_asparagine_180C_Ye2024":
        "974846cfc0fa87780e208fbc882310c41efcd1caa1d15fc4ae527f648edd204b",
    "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021":
        "fd26c70a0d0020ebbe1fca8dfd63a8cd2668ac008c7aa35ea8e4b7bfc98b2f81",
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3":
        "0d18247402d397ef1c35a328d80886b111d6fe42ac12454e0847a82b1cead79e",
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7":
        "020e775444dfc3cf9b5c32c0d6c15c5645d41d171c34ddba987fea99e481ecf3",
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3":
        "a3d126ac3fab64c57f8a463659e3e89474a2405e152e6c72e33272b37e8b8927",
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7":
        "7c78672fbcba460489e0528dd6becc24a530bb67fde6cf73aaee5e6c6480a183",
    # re-pinned 2026-09-03 (B9): the bundle returned to the hold-out with its notes,
    # evidence_class and hold_out_history updated; measured block untouched.
    "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5":
        "008cfa372890a196f8b9183ae2ab84800b64238e182d905f3a134499300698b5",
    "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026":
        "03d56fd60b097c85c731c87181f5269a4ea9aa4e4bd38c3f5267cef5f66d8251",
    "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026":
        "3a51f4621291f6131f01aaa0e1a30b555dcb9ad92fa9aa4f8fa9df086d0996f9",
    "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026":
        "060ee66bbad61ba5cbdb8e683c325a5eed657094806920dc1223f792f381de82",
    "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026":
        "907d92f499f897fbdf596de4bd1eae530b45eaa4d03aaf1c8b842dbaca195ad6",
}


def _hash_without_buffer(payload):
    clone = json.loads(json.dumps(payload))
    conditions = clone.get("conditions")
    if isinstance(conditions, dict):
        conditions.pop("buffer", None)
    return hashlib.sha256(
        json.dumps(clone, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_there_are_exactly_twenty_one_frozen_bundles():
    assert len(list(BUNDLE_DIR.rglob("*.json"))) == 21  # xylose returned after B9 (2026-09-03)


def test_the_buffer_completion_changed_nothing_but_the_buffer():
    """
    THE PROOF AMENDMENT 9 ASKS FOR. Measured values, compound lists,
    evidence_class and roles are unchanged -- and so is everything else in
    every bundle except the one key the amendment authorises.
    """
    offences = []
    for path in sorted(BUNDLE_DIR.rglob("*.json")):
        payload = json.loads(path.read_text())
        key = payload.get("benchmark_id") or path.stem
        expected = BUNDLE_BASELINES.get(key)
        assert expected is not None, f"{key}: no baseline hash recorded"
        actual = _hash_without_buffer(payload)
        if actual != expected:
            offences.append(f"{key}: {expected[:12]} -> {actual[:12]}")
    assert not offences, (
        "the buffer completion touched something it may not touch:\n"
        + "\n".join(offences))


# ===========================================================================
# 3. THE BUFFER-BLOCK SCHEMA
# ===========================================================================

_ALLOWED_SCALES = {"bench_cooled", "at_temperature", "unstated"}
_ALLOWED_CLASSES = {
    "primary_source_pdf", "repo_verbatim_methods_quote", "unknown",
}


def _buffer_blocks():
    for path in sorted(BUNDLE_DIR.rglob("*.json")):
        payload = json.loads(path.read_text())
        key = payload.get("benchmark_id") or path.stem
        yield key, payload.get("conditions", {}).get("buffer")


def test_every_bundle_has_a_decided_buffer_block():
    for key, block in _buffer_blocks():
        assert block is not None, f"{key}: no buffer block"
        for field in ("species", "concentration_M", "initial_ph_scale",
                      "provenance_class", "provenance_note"):
            assert field in block, f"{key}: buffer block lacks {field}"
        assert block["initial_ph_scale"] in _ALLOWED_SCALES, key
        assert block["provenance_class"] in _ALLOWED_CLASSES, key
        assert len(block["provenance_note"]) > 120, (
            f"{key}: a provenance note without an anchor is not provenance")


def test_a_named_buffer_carries_a_molarity_and_an_unknown_never_does():
    for key, block in _buffer_blocks():
        species = block["species"]
        if species in ("none", "buffer_unknown"):
            assert block["concentration_M"] is None, key
        else:
            assert isinstance(block["concentration_M"], (int, float)), key
            assert block["concentration_M"] > 0.0, key


def test_an_unknown_buffer_must_say_the_source_is_not_available():
    """
    THE RULE AGAINST GUESSING. `buffer_unknown` is only legitimate where the
    evidence is absent, and the note has to say so in those words -- otherwise
    it is indistinguishable from laziness, and worse, a plausible guess could
    be relabelled as unknown to dodge scrutiny.
    """
    for key, block in _buffer_blocks():
        if block["species"] != "buffer_unknown":
            continue
        assert block["provenance_class"] == "unknown", key
        assert "NOT ON DISK" in block["provenance_note"].upper(), key


def test_a_positive_no_buffer_finding_is_distinct_from_unknown():
    """
    "the paper says water" and "we do not know" are different facts and the
    schema must keep them apart. At least one bundle must be each.
    """
    species = [block["species"] for _key, block in _buffer_blocks()]
    assert species.count("none") >= 1
    assert species.count("buffer_unknown") >= 1
    assert "none" != "buffer_unknown"


def test_a_second_hand_reading_says_so_in_its_own_note():
    for key, block in _buffer_blocks():
        if block["provenance_class"] != "repo_verbatim_methods_quote":
            continue
        note = block["provenance_note"].upper()
        assert "NOT ON DISK" in note and "SECOND-HAND" in note, key


def test_a_primary_reading_names_the_pdf_it_was_read_from():
    for key, block in _buffer_blocks():
        if block["provenance_class"] != "primary_source_pdf":
            continue
        assert "data/articles/" in block["provenance_note"], key


def test_the_ph_disagreements_are_reported_and_not_acted_on():
    """
    Five bundles' recorded `conditions.ph` disagrees with their own source.
    Amendment 9 licenses completing the CONDITION RECORD, not editing an
    executable condition -- so each disagreement must be written down and the
    pH must be left alone.
    """
    flagged = [key for key, block in _buffer_blocks()
               if "ph_disagreement" in block]
    assert len(flagged) >= 4, flagged
    for key, block in _buffer_blocks():
        if "ph_disagreement" not in block:
            continue
        assert len(block["ph_disagreement"]) > 40, key


# ===========================================================================
# 4. THE FIREWALL -- literal grep and SYSTEMS walk
# ===========================================================================
# DISCLOSURE, EXTENDED FOR THIS WAVE: completing the buffer fields required
# opening the frozen bundles, whose `content_verification.quoted_table` blocks
# print the Yiltirak measured values beside the methods quotation. Those values
# were ALREADY on B2.2's disclosed-as-seen list. What follows is the mechanical
# check that they were not used.

_HOLDOUT_LITERALS = (
    "5.907", "11.439", "60.400", "12.757", "73.157",
    "62.6", "8.195",
    "0.4 / 0.5 / 0.5",
    "6.88", "1.28", "3.29", "1.46", "2.4", "1.68", "1.71", "1.62",
    "525.62", "325.22", "50.07", "582.34", "436.63",
    "696.99", "813.65", "59.70", "553.0", "229.0",
)

_FIT_SIDE_FILES = (
    "src/kinetic_core/sulfur.py",
    "src/kinetic_core/parameters_sulfur.py",
    "src/kinetic_core/species_sulfur.py",
    "src/kinetic_core/ph_state.py",
    "src/kinetic_core/network.py",
    "scripts/generators/generate_kinetic_core_b2_3_fit.py",
)


def test_no_holdout_literal_appears_on_the_fit_side():
    import re

    offences = []
    for relative in _FIT_SIDE_FILES:
        text = (ROOT / relative).read_text()
        for literal in _HOLDOUT_LITERALS:
            # WAVE B2.4 -- ONE CHARACTER, AND THE REASON IS ON THE RECORD.
            # The literal list contains "2.4" (a Yiltirak hold-out fold). That
            # substring also occurs inside the WAVE NAME "B2.4", so as first
            # written this test made it impossible for any fit-side file to
            # mention the fourth wave of its own module -- which is a defect in
            # the guard, not in the prose. The added lookbehind exempts a match
            # glued to a letter and nothing else: a measured value is never
            # immediately preceded by a letter, so the guard loses no reach.
            pattern = re.compile(
                r"(?<![A-Za-z])(?<![0-9.eE])" + re.escape(literal)
                + r"(?![0-9eE])")
            for line_no, line in enumerate(text.splitlines(), 1):
                if not pattern.search(line):
                    continue
                upper = line.upper()
                if "HOLD-OUT" in upper or "IS NOT HERE" in upper:
                    continue
                offences.append(
                    f"{relative}:{line_no}: {literal!r} in {line.strip()[:90]}")
    assert not offences, (
        "hold-out literals leaked onto the fit side:\n" + "\n".join(offences))


def test_the_frozen_bundles_are_never_opened_by_the_fit_path():
    """
    B2.3 makes this test MORE necessary, not less: one script in this repo now
    legitimately opens all 21 bundles. It is not on the fit side and it must
    never become so.
    """
    io_tokens = ("open(", "read_text", "read_bytes", "json.load", "glob(",
                 "iterdir", "rglob")
    for relative in _FIT_SIDE_FILES + (
        "scripts/generators/generate_kinetic_core_b2_3_holdout.py",
    ):
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "external_validation" not in line:
                continue
            assert not any(token in line for token in io_tokens), (
                f"{relative}:{line_no} reads a frozen bundle: {line.strip()[:100]}"
            )


def test_the_fit_never_integrates_a_holdout_condition():
    """THE SYSTEMS WALK -- B2.2's, re-run against B2.3's own generator."""
    import scripts.generators.generate_kinetic_core_b2_3_fit as fit

    for name, spec in fit.SYSTEMS.items():
        t_c = float(spec["t_c"])
        ph = float(spec["ph"])
        if name.startswith("kang"):
            assert t_c in (100.0, 120.0), f"{name} at {t_c} C is the hold-out rung"
        if name.startswith("zhou") and not name.startswith("zhou_drift_"):
            assert ph == 7.0, f"{name} at pH {ph} is a held-out column"
        assert ph != 9.0, f"{name} sits on Sun 2019's held-out pH-9 column"
        assert ph != 6.5, f"{name} sits on Whitfield 2001's held-out pH"
        if name.startswith("hofmann_") and "pH5" not in name:
            raise AssertionError(f"unexpected Hofmann fit system {name}")
        if {"PENT", "Cys"} <= set(spec["initial"]) and len(spec["initial"]) <= 3:
            assert not (
                t_c in (100.0, 110.0, 120.0, 130.0)
                and float(spec["minutes"]) in (240.0, 120.0, 60.0, 30.0)
            ), f"{name} sits on a Yiltirak hold-out rung"


def test_the_fit_rows_are_unchanged_from_b2_2():
    """
    Amendment 9 says REFIT ON THE UNCHANGED DECLARED FIT ROWS. A wave that
    quietly adds a row is not a refit, and this is the mechanical check.
    """
    import scripts.generators.generate_kinetic_core_b2_2_fit as b22
    import scripts.generators.generate_kinetic_core_b2_3_fit as b23

    assert len(b23.ACTIVE_FIT_ROWS) == len(b22.ACTIVE_FIT_ROWS)
    assert ({r["id"] for r in b23.ACTIVE_FIT_ROWS}
            == {r["id"] for r in b22.ACTIVE_FIT_ROWS})
    for a, b in zip(b22.ACTIVE_FIT_ROWS, b23.ACTIVE_FIT_ROWS):
        assert a["target"] == b["target"], a["id"]
        assert a["sigma_log"] == b["sigma_log"], a["id"]
    assert set(b23.SYSTEMS) == set(b22.SYSTEMS)



# ===========================================================================
# 5. THE FROZEN FIT, once it exists
# ===========================================================================


def test_the_frozen_fit_adds_no_parameter():
    frozen = _frozen()["frozen_parameters"]
    n = (len(frozen["log10_k_ref_at_145C"]) + 1
         + len(frozen["decay_Ea_kJ_mol"]) + len(frozen["ph_drift"]))
    assert n == 48, (
        f"B2.3 must have B2.2's 48 free parameters, not {n}. A conservation "
        f"fix that costs a parameter is not a conservation fix.")


def test_the_kang_pot_no_longer_runs_alkaline():
    """
    THE DEFECT, PINNED SHUT. B2.2's diagnosis sec. 3 measured Kang's 100 and
    120 C pots finishing at a cooled pH of 11.42 against a measured 4.9. That
    number is the reason this wave exists, so a regression past pH 7 must break
    a test rather than a report.
    """
    import numpy as np

    import scripts.generators.generate_kinetic_core_b2_3_fit as b23

    frozen = _frozen()["frozen_parameters"]
    x = np.array(
        [frozen["log10_k_ref_at_145C"][k] for k in b23.PARAM_ORDER]
        + [frozen["lumped_formation_Ea_kJ_mol"]]
        + [frozen["decay_Ea_kJ_mol"][f] for f in b23.DECAY_FAMILY_ORDER]
        + [frozen["ph_drift"]["acid_yield_per_sink_event"],
           frozen["ph_drift"]["arp_secondary_ammonium_pKa"]])
    _f, _e, _d, drift = b23.unpack(x)
    run = S.integrate_sulfur(
        b23.build_parameters(x), 393.15,
        {"TTCA": 10.0, "OX": b23.OX_AMBIENT_MMOL_L},
        np.array([0.0, 120.0]), ph=7.0, buffer_spec=b23.BUFFER_NONE,
        ph_drift=drift, ph_nodes=41, ph_iterations=4, rtol=1e-7, atol=1e-14)
    assert run.metadata["ph_final_cooled"] < 7.0, (
        f"Kang's TTCA pot finishes at pH "
        f"{run.metadata['ph_final_cooled']:.2f}; the charge leak is back")
