"""
tests/unit/test_kinetic_core_b5_cutover.py

Build Wave B5 -- THE PROPAGATOR CUTOVER. These tests guard the four properties
the cutover is only worth having if it keeps:

  1. ENVELOPE ENFORCEMENT. The core refuses what it cannot name, and a refusal
     carries no number. Asking a refused prediction for a value RAISES.
  2. NO FAST ABSOLUTES IN USER OUTPUT. The demoted screening lane's ppb cannot
     reach a terminal or a JSON payload, by any route the CLI offers.
  3. LABELS. Every FAST surface says ORDINAL SCREENING; every core surface says
     which lane answered.
  4. EXAM SCHEMA. The final-exam artifact keeps the shape the README and the
     validation contract quote from.

These are contract tests, not chemistry tests: they assert the SHAPE of the
refusals and the ABSENCE of leaked fields, and they pin no rate constant.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.kinetic_core.engine import (
    ACRYLAMIDE,
    SULFUR,
    TRUNK,
    EnvelopeDeclaration,
    FormulationSpec,
    OutOfEnvelope,
    ProcessSpec,
    ThermalProgram,
    declare_envelope,
    default_targets_for,
    engine_metadata,
    predict,
    resolve_lane,
    resolve_matrix,
)

ROOT = Path(__file__).resolve().parents[2]

MFT = "2-Methyl-3-furanthiol (MFT)"
FFT = "2-Furfurylthiol (FFT)"


def _sulfur_spec(temperature_c: float = 145.0, minutes: float = 20.0, ph: float = 5.0):
    return FormulationSpec(
        "sulfur",
        {"D-Ribose": 100.0, "L-Cysteine": 33.0},
        ProcessSpec(ThermalProgram.isothermal(temperature_c, minutes), ph=ph),
    )


def _acrylamide_spec():
    return FormulationSpec(
        "acrylamide",
        {"D-Glucose": 27.75, "L-Asparagine": 33.3},
        ProcessSpec(ThermalProgram.isothermal(180.0, 30.0), ph=6.0, water_activity=0.99),
    )


# ---------------------------------------------------------------------------
# 1. Envelope enforcement
# ---------------------------------------------------------------------------


def test_sulfur_request_is_in_envelope_and_routes_to_the_sulfur_lane():
    """
    B2.2 CHANGES THE EXPECTED STATE HERE, DELIBERATELY.

    The sulfur lane now carries a pH-TRAJECTORY state, and a pH trajectory
    needs a buffer. A request that does not declare one is answered -- but it
    is answered by EXTRAPOLATING the trajectory from water autoprotolysis and
    the charged solutes alone, which for a real buffered experiment
    over-predicts the drift. That is an extrapolation and the engine now says
    so, which is why the state is `in_envelope_extrapolated` rather than
    `in_envelope`. The diagnostic below is what B2.2's buffer probe measured:
    the same pot answered unbuffered and answered at Hofmann's declared 0.5 M
    phosphate differs by up to ~15x in predicted thiol.
    """
    run = predict(_sulfur_spec(), [MFT, FFT])
    assert run.answered
    assert run.declaration.lane == SULFUR
    assert run.declaration.state == "in_envelope_extrapolated"
    assert any("no buffer was declared" in w for w in run.declaration.warnings)
    assert set(run.concentrations_ug_per_l) == {MFT, FFT}
    assert all(v >= 0.0 for v in run.concentrations_ug_per_l.values())


def test_declaring_the_buffer_removes_the_extrapolation_flag():
    """
    The other half: the flag is not decoration. Supply the buffer the source
    actually states and the request is back inside the envelope.
    """
    from src.kinetic_core.engine import ProcessSpec, ThermalProgram
    from src.kinetic_core.ph_state import BufferSpec

    spec = FormulationSpec(
        "sulfur",
        {"D-Ribose": 100.0, "L-Cysteine": 33.0},
        ProcessSpec(
            ThermalProgram.isothermal(145.0, 20.0), ph=5.0,
            buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.5,
                              declared=True, source="Hofmann 1998 Methods"),
        ),
    )
    run = predict(spec, [MFT, FFT])
    assert run.answered
    assert not any("no buffer was declared" in w for w in run.declaration.warnings)
    assert run.declaration.state == "in_envelope"


def test_acrylamide_request_routes_to_the_acrylamide_lane():
    run = predict(_acrylamide_spec(), ["Acrylamide"])
    assert run.answered
    assert run.declaration.lane == ACRYLAMIDE
    assert run.concentrations_ug_per_l["Acrylamide"] > 0.0


@pytest.mark.parametrize(
    "compound, token",
    [
        # WAVE B7 REMOVED THESE TWO. 5-HMF and DMHF are now SPECIES with
        # routes -- see tests/unit/test_kinetic_core_b7.py -- and the compounds
        # that replace them here are the neighbours whose refusals B7 made
        # SHARPER rather than dropping. That is the same discipline B6 applied
        # to 1-hexanol and 2-pentylfuran two lines below, and the two B5
        # reasons that are retired are quoted verbatim in the B7 report so the
        # change of verdict is legible:
        #   "the hexose-dehydration route that forms it was never
        #    parameterised" -- NOW FALSE (Kocadagli's amine-free glucose
        #    system is ingested whole);
        #   "reporting NF as DMHF would be a species substitution the corpus
        #    does not license" -- STILL TRUE, and honoured by DMHF being its
        #    OWN species with its own route rather than an alias of NF.
        ("HEMF", "alanine"),
        ("2,5-dimethyl-4-hydroxy-3(2H)-thiophenone", "Haleva-Toledo"),
        # WAVE B6 updated these four rows. The lipid lane now EXISTS, so
        # "no lipid-oxidation path" is no longer true of any of them. What is
        # asserted instead is that each is still refused for a reason that is
        # SHARPER than the old one, and that the two whose branch fraction is
        # measured are answered elsewhere (see tests/unit/test_kinetic_core_b6.py).
        ("1-hexanol", "no aldehyde-reduction step"),
        ("2-pentylfuran", "not in Frankel"),
    ],
)
def test_unrepresented_compounds_are_declared_out_with_a_named_reason(compound, token):
    """The core must NAME why it cannot answer, not merely fail."""
    run = predict(_acrylamide_spec(), [compound])
    assert not run.answered
    assert run.declaration.state == "out_of_envelope"
    assert run.concentrations_ug_per_l == {}
    assert any(
        token.lower() in reason.lower() for reason in run.declaration.reasons
    ), (
        f"{compound}: no reason mentioned {token!r}; reasons were "
        f"{run.declaration.reasons}"
    )


def test_a_refused_prediction_raises_rather_than_returning_a_number():
    """
    The property the whole engine exists for: there is NO state in which a
    declined request hands back a float.
    """
    # B6: hexanal is a species now, but this charge declares NO LIPID CARRIER,
    # so the request is still refused -- and still refused with NO number.
    run = predict(_acrylamide_spec(), ["hexanal"])
    assert not run.answered
    with pytest.raises(OutOfEnvelope):
        run.require("hexanal")
    with pytest.raises(OutOfEnvelope):
        run.require("Acrylamide")


def test_an_intact_protein_is_not_a_precursor_the_core_can_charge():
    spec = FormulationSpec(
        "matrix",
        {"Soy Protein Isolate": 1000.0},
        ProcessSpec(ThermalProgram.isothermal(160.0, 25.0), ph=7.0),
    )
    # B6: a protein isolate is now recognised as a LIPID CARRIER, so a hexanal
    # request is answered -- but the half of the pre-B6 refusal that was RIGHT
    # survives: it charges NO Maillard network, and its declared "1000 mM" is
    # ignored rather than believed.
    declaration = declare_envelope(spec, ["hexanal"])
    assert declaration.mapped_precursors == {}
    assert declaration.lipid_carriers == ("soy_protein_isolate",)
    assert any("IGNORED" in w for w in declaration.warnings)

    # And a MAILLARD target from the same charge is still refused outright.
    run = predict(spec, ["melanoidins"])
    assert run.answered is False or run.concentrations_ug_per_l["melanoidins"] == 0.0


def test_the_lanes_do_not_compose_and_a_cross_lane_request_is_declared():
    """
    A request spanning the sulfur and acrylamide lanes is UNANSWERABLE, not
    hard: the acrylamide network deliberately omits every sulfur step.
    """
    lane, reasons = resolve_lane(["Cys", "Asn"], ["MFT", "ACR"])
    assert lane is None
    assert reasons and "LANE CONFLICT" in reasons[0]

    spec = FormulationSpec(
        "both",
        {"D-Glucose": 100.0, "L-Cysteine": 33.0, "L-Asparagine": 33.0},
        ProcessSpec(ThermalProgram.isothermal(160.0, 20.0), ph=6.0),
    )
    run = predict(spec, [MFT, "Acrylamide"])
    assert not run.answered
    assert any("LANE CONFLICT" in r for r in run.declaration.reasons)


def test_alanine_plus_furfural_is_the_schibilsky_declension():
    """The exam's furan-family refusal, pinned so it cannot regress to a number."""
    spec = FormulationSpec(
        "schibilsky",
        {"D-Glucose": 250.0, "L-Alanine": 250.0},
        ProcessSpec(ThermalProgram.isothermal(130.0, 120.0), ph=5.0),
    )
    run = predict(spec, ["Furfural"])
    assert not run.answered
    assert run.concentrations_ug_per_l == {}


def test_a_sulfur_target_without_a_sulfur_source_is_declared():
    spec = FormulationSpec(
        "no-sulfur",
        {"D-Glucose": 100.0, "Glycine": 100.0},
        ProcessSpec(ThermalProgram.isothermal(145.0, 20.0), ph=5.0),
    )
    run = predict(spec, [MFT])
    assert not run.answered
    assert any("NO sulfur source" in r for r in run.declaration.reasons)


def test_out_of_range_conditions_are_a_DECLARED_extrapolation_not_a_refusal():
    """
    An extrapolation still answers -- but it says so, and the warning travels
    with the number rather than being discovered later.
    """
    run = predict(_sulfur_spec(temperature_c=250.0, minutes=5.0), [MFT])
    assert run.answered
    assert run.declaration.state == "in_envelope_extrapolated"
    assert any("200 C" in w or "extrapolation" in w for w in run.declaration.warnings)


def test_the_pH_free_lanes_declare_that_they_ignore_pH():
    run = predict(_acrylamide_spec(), ["Acrylamide"])
    assert run.answered
    assert any("NO pH term" in w for w in run.declaration.warnings)


def test_declaration_serialises_with_its_reasons():
    declaration = declare_envelope(_acrylamide_spec(), ["1-hexanol"])
    payload = declaration.as_dict()
    assert payload["state"] == "out_of_envelope"
    assert payload["unrepresented_targets"]
    assert all(
        entry["compound"] and entry["reason"]
        for entry in payload["unrepresented_targets"]
    )


def test_default_targets_follow_the_resolved_lane():
    assert MFT.lower() in [t.lower() for t in default_targets_for({"D-Ribose": 1.0, "L-Cysteine": 1.0})]
    assert "acrylamide" in [t.lower() for t in default_targets_for({"D-Glucose": 1.0, "L-Asparagine": 1.0})]
    # B6: a protein isolate resolves to the LIPID lane, whose default targets
    # are Frankel's measured slate. Before B6 it resolved to nothing.
    assert "hexanal" in default_targets_for({"Soy Protein Isolate": 1.0})
    assert default_targets_for({"Whey Protein Hydrolysate": 1.0}) == ()


def test_aqueous_matrix_aliases_resolve_to_water_but_a_real_matrix_does_not():
    """
    The wiring bug B5 found: `protein_type: free` is not a threshold matrix, and
    a compound with a measured water threshold was reported as having none.
    """
    for descriptor in ("free", "water", "aqueous", "buffer", ""):
        assert resolve_matrix(descriptor) == "water"
    assert resolve_matrix("soy_paste_hong") == "soy_paste_hong"


def test_the_oav_table_is_reachable_and_finds_measured_thresholds():
    run = predict(_sulfur_spec(), [MFT, FFT])
    table = run.oav()
    per_species = table["per_species"]
    # MFT has a measured water threshold; a None here is the keying bug.
    assert per_species["MFT"]["threshold_ug_per_L"] is not None
    assert per_species["MFT"]["OAV_point"] is not None
    assert per_species["MFT"]["OAV_interval"][0] < per_species["MFT"]["OAV_interval"][1]


def test_engine_declares_that_it_fits_nothing():
    metadata = engine_metadata()
    assert metadata["fits_anything"] is False
    assert metadata["lanes_compose"] is False
    # B6 adds the LIPID lane. `lanes_compose` stays False because it describes
    # the three MAILLARD lanes; the lipid lane's separate co-integration ruling
    # is reported under its own key.
    assert set(metadata["lanes"]) == {TRUNK, SULFUR, ACRYLAMIDE, "lipid"}
    assert metadata["lipid_lane_cointegrates"]["rule"].startswith("direct sum")


# ---------------------------------------------------------------------------
# 2 + 3. CLI routing, FAST demotion, labels
# ---------------------------------------------------------------------------

_FORBIDDEN_FAST_FIELDS = ("a_ppb", "b_ppb", "predicted_ppb", "range_p5", "range_p95")


def test_screening_payload_strips_every_fast_absolute_field():
    from src.comparative_cli import SCREENING_LABEL, screening_payload

    raw = {
        "artifact": "maillard_predict",
        "rows": [
            {
                "compound": "MFT",
                "predicted_ppb": 1234.5,
                "range_p5": 1.0,
                "range_p95": 9.0,
                "range_available": True,
                "dominant_pathway": "x",
                "lane_reliability": None,
            }
        ],
        "caveats": {"absolute": "..."},
    }
    clean = screening_payload(raw)
    text = json.dumps(clean)
    for field in _FORBIDDEN_FAST_FIELDS:
        assert f'"{field}"' not in text, f"{field} survived the screening sanitiser"
    assert clean["lane_label"] == SCREENING_LABEL
    assert clean["absolutes_withheld"] is True
    assert clean["rows"][0]["absolute_ppb_withheld"] is True
    # The ORDERING -- the lane's measured product -- must survive.
    assert clean["rows"][0]["compound"] == "MFT"
    assert clean["rows"][0]["dominant_pathway"] == "x"


def test_screening_label_appears_in_the_rendered_fast_text():
    from src.comparative_cli import (
        SCREENING_LABEL,
        render_predict_text,
        screening_payload,
    )

    raw = {
        "artifact": "maillard_predict",
        "system": {"name": "demo", "spec": {}},
        "rows": [
            {
                "compound": "MFT",
                "predicted_ppb": 1234.5,
                "range_p5": 1.0,
                "range_p95": 9.0,
                "range_available": True,
                "dominant_pathway": "x",
                "lane_reliability": None,
                "sulfur": False,
            }
        ],
        "warnings": [],
        "caveats": {"absolute": "...", "sulfur": None, "no_range": "..."},
    }
    text = render_predict_text(screening_payload(raw))
    assert SCREENING_LABEL in text
    assert "1234" not in text, "a FAST absolute reached the rendered text"
    assert "ppb withheld" in text


def test_cli_defaults_to_the_core_lane():
    import scripts.maillard as cli

    args = cli.build_parser().parse_args(["predict", "spec.yml"])
    assert args.lane == "core"
    args = cli.build_parser().parse_args(["compare", "spec.yml"])
    assert args.lane == "core"


def test_cli_offers_both_lanes_and_no_others():
    import scripts.maillard as cli

    args = cli.build_parser().parse_args(["predict", "spec.yml", "--lane", "fast"])
    assert args.lane == "fast"
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["predict", "spec.yml", "--lane", "nonsense"])


def test_fast_lane_refuses_absolute_flag(tmp_path, capsys):
    """`--absolute` on the screening lane must be REFUSED, not quietly honoured."""
    import scripts.maillard as cli
    from src.comparative_cli import SPEC_TEMPLATE

    spec = tmp_path / "spec.yml"
    spec.write_text(SPEC_TEMPLATE)
    code = cli.main(["compare", str(spec), "--lane", "fast", "--absolute"])
    assert code == 2
    assert "not available on the screening lane" in capsys.readouterr().err


def test_core_predict_payload_carries_its_declaration_and_no_fast_fields():
    from src.comparative_cli import predict_core

    payload = predict_core(
        {
            "name": "sulfur",
            "precursors": {"D-Ribose": 100.0, "L-Cysteine": 33.0},
            "temp_C": 145.0,
            "time_min": 20.0,
            "ph": 5.0,
            "aw": 0.98,
            "protein_type": "free",
        }
    )
    assert payload["lane"] == "core"
    assert payload["answered"] is True
    assert payload["declaration"]["lane"] == SULFUR
    text = json.dumps(payload, default=str)
    for field in _FORBIDDEN_FAST_FIELDS:
        assert f'"{field}"' not in text
    assert payload["rows"], "the core lane returned no rows for an in-envelope spec"
    assert all("predicted_ug_per_l" in row for row in payload["rows"])


def test_core_predict_declares_rather_than_guessing_for_a_matrix_spec():
    from src.comparative_cli import predict_core

    payload = predict_core(
        {
            "name": "matrix",
            "precursors": {"Soy Protein Isolate": 1000.0},
            "temp_C": 160.0,
            "time_min": 25.0,
            "ph": 7.0,
            "aw": 0.35,
            "protein_type": "soy_iso",
        }
    )
    # WAVE B6. This spec is now ANSWERED, on the lipid lane -- but the thing
    # this test was written to protect is unchanged and is asserted directly:
    # the answer is declared an EXTRAPOLATION, it names the rate assumption,
    # and no number leaves without one.
    assert payload["declaration"]["state"] == "in_envelope_extrapolated"
    assert "lipid" in payload["declaration"]["lanes"]
    assert any(
        "ASSUMPTION, NOT A MEASUREMENT" in w
        for w in payload["declaration"]["warnings"]
    )


# ---------------------------------------------------------------------------
# 4. Exam-report schema
# ---------------------------------------------------------------------------

EXAM_JSON = ROOT / "results" / "validation" / "cutover_final_exam.json"
PREREG = ROOT / "results" / "validation" / "cutover_prereg.md"


def test_the_prereg_exists_and_names_the_out_of_envelope_families():
    assert PREREG.exists(), "the cutover pre-registration is missing"
    text = PREREG.read_text()
    for token in ("HMF", "DMHF", "lipid-oxidation", "OUT", "23 IN / 17 OUT"):
        assert token in text


@pytest.mark.skipif(not EXAM_JSON.exists(), reason="final exam not yet generated")
def test_exam_report_schema():
    payload = json.loads(EXAM_JSON.read_text())
    assert payload["artifact"] == "cutover_final_exam"
    assert payload["bundle_count"] == 21
    assert payload["pass_band_level_x"] == 3.0
    assert payload["pre_registration"].endswith("cutover_prereg.md")

    summary = payload["summary"]
    for key in ("core_answered", "core_declined", "core_scored", "core", "old_lane_all_points"):
        assert key in summary
    # The paired subset is the only apples-to-apples comparison and must be present.
    paired = summary["paired_subset"]
    assert paired["n"] > 0
    assert paired["core"]["median_fold_error"] is not None
    assert paired["old"]["median_fold_error"] is not None

    assert summary["core_answered"] + summary["core_declined"] == payload["point_count"]

    for row in payload["rows"]:
        for key in (
            "benchmark_id", "compound", "family", "envelope_state",
            "answered", "core_predicted", "core_fold_error", "old_fold_error",
        ):
            assert key in row, f"row missing {key}"
        # THE INVARIANT: a declined point carries no number.
        if not row["answered"]:
            assert row["core_predicted"] is None
            assert row["core_fold_error"] is None
            assert row["declaration_reasons"], "a declension with no reason"


@pytest.mark.skipif(not EXAM_JSON.exists(), reason="final exam not yet generated")
def test_exam_declines_every_lipid_and_hmf_point():
    """
    The pre-registered envelope, pinned against the artifact.

    WAVE B6 NARROWED THIS. HEXANAL is now answered: its branch fraction is
    measured (Frankel 1989, the Module 5 FIT column). Everything else on the
    original list is still declined, and for each one the REASON is what this
    test now pins:

      * nonanal      -- the oleate -> nonanal branch fraction is unmeasured;
      * 1-hexanol    -- no aldehyde-reduction step exists in the corpus;
      * 2-pentylfuran -- not in Frankel's six-product slate;
      * HMF / DMHF   -- NO LONGER ON THIS LIST. Wave B7 gives both a route,
                        and it is NOT a branch fraction: the HMF node's source
                        topology is adopted from four independent published
                        networks that all write it the same way, its seven
                        constants are INGESTED from one declared-FIT system
                        with nothing fitted, and the DMHF node has exactly ONE
                        fitted number, calibrated on a pentose system that
                        appears in no exam bundle. The limb SHARE is computed
                        from the pools and moves from 2 % to 96 % across
                        conditions, which is what makes it not a fraction.
                        tests/unit/test_kinetic_core_b7.py pins all of that.

    A future wave that quietly starts answering any of those has invented a
    branch fraction, and this test is where that shows up.
    """
    payload = json.loads(EXAM_JSON.read_text())
    still_declined = ("nonanal", "hexanol", "pentylfuran")
    for row in payload["rows"]:
        compound = row["compound"].lower()
        if compound == "hexanal":
            continue
        if any(token in compound for token in still_declined):
            assert not row["answered"], f"{row['compound']} should be declined"


@pytest.mark.skipif(not EXAM_JSON.exists(), reason="final exam not yet generated")
def test_exam_records_the_prereg_checks():
    payload = json.loads(EXAM_JSON.read_text())
    checks = payload["prereg_checks"]
    assert checks, "the exam must score its own pre-registration"
    outcomes = {c["outcome"] for c in checks}
    # A scorecard on which nothing can fail is not a scorecard.
    assert outcomes - {"HELD"}, "no pre-registered claim was falsified -- verify honestly"
