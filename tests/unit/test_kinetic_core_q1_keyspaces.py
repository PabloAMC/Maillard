"""
Wave Q1 -- the key-space contract, and the report-flow bugs it was hiding.

WHY THIS FILE EXISTS
--------------------
Q1 found that three of the V1 API-feedback items were not merely ergonomic:
they were concealing live defects that NO test covered.

  * ``comparative_cli.predict_core`` looked up each row's OAV as
    ``per_species[species_key]``, but ``per_species`` is keyed by B4 KEY. The
    two spaces are identical for the trunk and sulfur lanes and DIFFERENT for
    the lipid one ("HEXANAL" vs "hexanal"), so every lipid compound's lookup
    missed and every report announced "no threshold" for two compounds that
    have a measured one.
  * Report rows read their interval out of the OAV table, which has no entry at
    all for a ``NO_B4_RECORD`` species -- so four of the lipid lane's seven
    products were printed as bare points, against the report layer's own rule
    that an absolute is never drawn without its interval.
  * ``engine.compare`` dropped its arms' OAV tables, so the HTML report rebuilt
    them by hand; the copy never learned B7's furanone band and drew narrower
    intervals than a predict of the identical arm.

Every one of those failures is SILENT: a key from the wrong space does not
raise, it returns "no measured threshold". So the tests below pin the
observable consequence, not the implementation.
"""

import pytest

from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram
from src.kinetic_core.engine import compare as core_compare
from src.kinetic_core.engine import predict
from src.kinetic_core.keyspaces import b4_key, keys_for, no_b4_reason
from src.kinetic_core.matrix_oav import compare_formulations
from src.kinetic_core.species_lipid import B4_COMPOUND_KEY, NO_B4_RECORD


def _lipid_spec(name: str = "carrier", cysteine: float = 5.0) -> FormulationSpec:
    """A CO-INTEGRATED lipid + sulfur charge: the case the bugs lived in."""
    return FormulationSpec(
        name=name,
        precursors={
            "pea protein isolate": 100.0,
            "L-Cysteine": cysteine,
            "D-Ribose": 5.0,
        },
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(140.0, 30.0),
            ph=6.5,
            water_activity=0.9,
            matrix="water",
        ),
    )


# ---------------------------------------------------------------------------
# The mapping itself
# ---------------------------------------------------------------------------


def test_the_two_key_spaces_actually_differ_on_the_lipid_lane():
    """If this ever became an identity map the bugs below would be untestable."""
    assert B4_COMPOUND_KEY["HEXANAL"] == "hexanal"
    assert B4_COMPOUND_KEY["DECADIENAL"] == "tt_2_4_decadienal"
    assert b4_key("HEXANAL") == "hexanal"
    # ...and an identity for every lane whose keys were never re-transcribed.
    assert b4_key("MFT") == "MFT"


def test_a_species_with_no_b4_record_has_NO_b4_key_and_a_recorded_reason():
    for species, reason in NO_B4_RECORD.items():
        assert b4_key(species) is None, (
            f"{species} has no B4 record; returning a key for it is exactly how "
            "a table lookup comes back as 'no measured threshold' with no reason"
        )
        assert no_b4_reason(species) == reason
        assert reason.strip(), "a NO_B4_RECORD entry must SAY why"


def test_keys_for_composes_both_hops():
    keys = keys_for("hexanal", {"hexanal": "HEXANAL"})
    assert (keys.display, keys.species, keys.b4) == ("hexanal", "HEXANAL", "hexanal")
    assert keys.has_b4_record is True
    assert keys.no_b4_reason is None

    keys = keys_for("pentane", {"pentane": "PENTANE"})
    assert keys.b4 is None and keys.has_b4_record is False
    assert "hydrocarbon" in keys.no_b4_reason


# ---------------------------------------------------------------------------
# The bug the mapping was hiding
# ---------------------------------------------------------------------------


def test_a_lipid_compound_with_a_measured_threshold_reports_its_OAV():
    """
    THE REGRESSION. hexanal and 2,4-decadienal both have measured B4 thresholds.
    Before Q1 every report said "no threshold" for both, because the row was
    looked up in the wrong key-space.
    """
    from src.comparative_cli import predict_core

    payload = predict_core(
        {
            "name": "lipid",
            "precursors": {
                "pea protein isolate": 100.0,
                "L-Cysteine": 5.0,
                "D-Ribose": 5.0,
            },
            "temp_C": 140.0,
            "time_min": 30.0,
            "ph": 6.5,
            "aw": 0.9,
            "matrix": "water",
        }
    )
    assert payload["answered"] is True
    by_compound = {row["compound"]: row for row in payload["rows"]}
    for compound in ("hexanal", "2,4-decadienal"):
        row = by_compound[compound]
        assert row["oav"] is not None, f"{compound} lost its OAV entirely"
        assert row["oav"]["available"] is True, (
            f"{compound} has a measured B4 threshold but the row reports "
            f"{row['oav'].get('reason')!r}"
        )
        assert row["oav"]["oav"] > 0.0
        assert row["oav"]["threshold_source"]


def test_a_compound_with_no_b4_record_still_carries_its_interval_and_its_reason():
    """
    A compound with no measured ODOUR THRESHOLD still has a perfectly
    well-defined CONCENTRATION INTERVAL. Refusing to print it was an accident
    of where the number was stored, not a statement about the evidence.
    """
    from src.comparative_cli import predict_core

    payload = predict_core(
        {
            "name": "lipid",
            "precursors": {"pea protein isolate": 100.0},
            "temp_C": 140.0,
            "time_min": 30.0,
            "ph": 6.5,
            "aw": 0.9,
            "matrix": "water",
        }
    )
    assert payload["answered"] is True
    rows = {row["compound"]: row for row in payload["rows"]}
    assert "pentane" in rows, "the NO_B4_RECORD species must still be reported"
    row = rows["pentane"]
    lo, hi = row["interval_ug_per_l"]
    assert lo is not None and hi is not None and hi > lo > 0.0
    assert row["band_x"] > 1.0
    # ...and the refusal to give it an OAV names its recorded reason.
    assert row["oav"]["available"] is False
    assert "hydrocarbon" in row["oav"]["reason"]


def test_every_answered_row_carries_its_own_interval():
    """The report layer must never have to join a row back onto another table."""
    run = predict(_lipid_spec(), ["hexanal", "pentane", "2-methyl-3-furanthiol (MFT)"])
    assert run.answered
    for row in run.interval_rows():
        lo, hi = row["interval_ug_per_l"]
        assert lo is not None and hi is not None, row["compound"]
        assert row["band_x"] is not None
        assert row["species_key"] and row["lane"]


# ---------------------------------------------------------------------------
# compare: the arms carry their own B4 output
# ---------------------------------------------------------------------------


def test_compare_emits_each_arm_s_own_oav_table_and_rows():
    payload = core_compare(
        _lipid_spec("high_cys", cysteine=20.0),
        _lipid_spec("control", cysteine=1.0),
        ["hexanal", "2-methyl-3-furanthiol (MFT)"],
    )
    assert payload["comparable"] is True
    for key in ("oav_table_a", "oav_table_b", "rows_a", "rows_b"):
        assert payload[key], f"compare dropped {key}"
    assert "hexanal" in payload["oav_table_a"]["per_species"]
    assert all(r["interval_ug_per_l"][1] is not None for r in payload["rows_a"])


def test_a_compare_arm_draws_the_same_interval_as_a_predict_of_that_arm():
    """
    The drift test. The report layer's hand-rolled copy of ``.oav()`` was never
    taught about B7's furanone declared-assumption band, so a compare drew
    NARROWER intervals than a predict of the identical arm. Both sides now come
    from the same method, and this pins that they agree.
    """
    spec = _lipid_spec("arm")
    targets = ["hexanal", "2-methyl-3-furanthiol (MFT)", "DMHF"]
    run = predict(spec, targets)
    payload = core_compare(spec, _lipid_spec("other", cysteine=1.0), targets)
    assert payload["comparable"] is True

    solo = {r["compound"]: r["interval_ug_per_l"] for r in run.interval_rows()}
    in_compare = {r["compound"]: r["interval_ug_per_l"] for r in payload["rows_a"]}
    assert solo and solo.keys() == in_compare.keys()
    for compound, interval in solo.items():
        assert in_compare[compound] == pytest.approx(interval, rel=1e-12), (
            f"{compound}: a compare arm and a predict of the same arm disagree "
            "on the interval"
        )


# ---------------------------------------------------------------------------
# n_resolved
# ---------------------------------------------------------------------------


def test_an_undefined_ratio_is_not_counted_as_resolved():
    """
    A ratio against zero is not a ratio, so it resolves nothing. It used to be
    counted in ``n_resolved`` because its row is built with
    ``within_reliability_band=False`` -- literally true, and read by every
    consumer as "resolved".
    """
    comparison = compare_formulations(
        {"MFT": 10.0, "FFT": 1.0, "ACTZ": 1.0},
        {"MFT": 0.0, "FFT": 1.0, "ACTZ": 100.0},
        label_a="A",
        label_b="B",
    )
    rows = {row["compound"]: row for row in comparison["rows"]}
    assert rows["MFT"]["direction"] == "undefined"
    assert comparison["n_compared"] == 3
    assert comparison["n_undefined"] == 1
    # ACTZ resolves (100x apart), FFT does not (identical), MFT is undefined.
    assert comparison["n_resolved"] == 1


def test_n_undefined_is_published_so_no_renderer_re_derives_the_predicate():
    comparison = compare_formulations({"MFT": 1.0}, {"MFT": 1.0})
    assert comparison["n_undefined"] == 0
    assert comparison["n_resolved"] == 0


# ---------------------------------------------------------------------------
# EnvelopeDeclaration.summary
# ---------------------------------------------------------------------------


def test_a_refusal_summarises_as_a_refusal_and_names_its_first_reason():
    run = predict(
        FormulationSpec(
            name="unnameable",
            precursors={"pea protein isolate": 100.0},
            process=ProcessSpec(
                thermal=ThermalProgram.isothermal(140.0, 30.0), ph=6.5,
                water_activity=0.9, matrix="water",
            ),
        ),
        ["1-hexanol"],
    )
    summary = run.declaration.summary
    assert summary.startswith("REFUSED")
    assert run.declaration.reasons[0] in summary
    assert run.declaration.as_dict()["summary"] == summary


def test_an_extrapolated_answer_never_summarises_as_a_plain_answer():
    run = predict(_lipid_spec(), ["hexanal"])
    assert run.answered
    assert run.declaration.state == "in_envelope_extrapolated"
    assert "EXTRAPOLATED" in run.declaration.summary


def test_the_summary_is_derived_and_cannot_disagree_with_its_fields():
    """No stored copy: there is no state in which the two can drift apart."""
    from dataclasses import fields

    from src.kinetic_core.engine import EnvelopeDeclaration

    assert "summary" not in {f.name for f in fields(EnvelopeDeclaration)}
    assert isinstance(EnvelopeDeclaration.summary, property)
