"""Wave B35 (2026-09-11): the benchmark provenance audit, and the silent-zero guard it found
(kinetic_core_b35_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import panel
from src.kinetic_core.engine import predict

BD = data_paths.BENCHMARKS_DIR
SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
ROASTED = BD / "external_validation" / "external_validation_bi_2020_roasted_pea_hexanal.json"
TRIK = BD / "pea_isolate_uht_140C_Trikusuma2019.json"


def test_a_pot_with_no_precursor_refuses_maillard_targets_instead_of_answering_zero():
    """THE GUARD, and the third instance of this bug family.

    A matrix-only charge declares a protein isolate, which is a LIPID CARRIER and deliberately not a
    precursor. The trunk, sulfur and acrylamide networks therefore integrate from an all-zero state
    and every species in them stays zero BY CONSTRUCTION. Until B35 that was reported as 0.0 with no
    refusal: furaneol scored 0.0 against a measured 2780 ug/kg, and 5-HMF would have done the same.
    """
    bench = json.loads(ROASTED.read_text())
    spec = panel.core_spec(bench)
    assert spec.precursors == {"Pea Protein Isolate": 1000.0}
    for compound in ("furaneol", "furfural", "5-HMF", "acrylamide"):
        run = predict(spec, [compound])
        assert not run.answered, f"{compound} answered on a pot with no precursor"
        reason = " ".join(run.declaration.reasons)
        assert "CHARGES NO PRECURSOR" in reason
        assert "THE CURE IS A CHARGE" in reason
    # the lipid lane is exempt: its charge IS the carrier
    assert predict(spec, ["hexanal"]).answered


def test_no_scored_row_anywhere_on_the_panel_is_exactly_zero():
    """A zero is the absence of a prediction, not a prediction. B28 said so and it stayed possible."""
    scores = json.loads(SCORES.read_text())
    zeros = [(b["benchmark_id"], c["compound"]) for b in scores["benchmarks"]
             for c in b["compounds"] if c["predicted"] == 0.0]
    assert zeros == [], f"scored rows with a prediction of exactly zero: {zeros}"


def test_the_two_stale_provenance_notes_are_corrected_and_their_claims_kept():
    for path, key in ((BD / "external_validation" / "maillard_path"
                       / "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021.json", "vessel"),
                      (BD / "external_validation"
                       / "external_validation_liu_2023_ppi_offnote_baseline.json", "buffer")):
        note = json.loads(path.read_text())["conditions"][key]["provenance_note"]
        assert note.startswith("WAVE B35") and "THE SOURCE IS ON DISK" in note
        assert "SUPERSEDED AND RETAINED AS THE AUDIT RECORD" in note
        assert note.index("THE SOURCE IS ON DISK") < note.index("SUPERSEDED")


def test_the_eight_audited_measurements_are_present_with_their_printed_values():
    roasted = json.loads(ROASTED.read_text())["measured_volatiles"]
    assert roasted["furaneol"]["conc_ppb"] == 2780.0
    assert roasted["2,5-dimethylpyrazine"]["conc_ppb"] == 5960.0
    assert roasted["furfural"]["conc_ppb"] == 327.0
    raw = json.loads((BD / "external_validation"
                      / "external_validation_bi_2020_raw_pea_hexanal.json").read_text())
    assert raw["measured_volatiles"]["nonanal"]["conc_ppb"] == 69.8
    t = json.loads(TRIK.read_text())
    for name, uht, ctrl in (("2,5-dimethylpyrazine", 2.29, 2.46), ("methional", 3.10, 0.55),
                            ("2-acetyl-1-pyrroline", 0.41, 0.29), ("(E,E)-2,4-decadienal", 46.9, 0.06)):
        assert t["measured_volatiles"][name]["conc_ppb"] == uht
        # every compound the source prints a control column for gets one declared (B31), not just three
        assert t["conditions"]["carried_volatiles"][name] == ctrl
    assert len(t["conditions"]["carried_volatiles"]) == 7


def test_the_refusals_are_the_ones_the_record_predicts():
    """T3: methional is B22's product and 2-acetyl-1-pyrroline is B24's, and both waves were refused
    with their steps left inert -- so the engine must refuse them by NAME rather than answer them."""
    scores = json.loads(SCORES.read_text())
    refused = {(r["benchmark_id"], r["compound"]): r["reason"] for r in scores["refused_compounds"]}
    for compound, wave in (("methional", "B22"), ("2-acetyl-1-pyrroline", "B24"),
                           ("2,5-dimethylpyrazine", "B18")):
        reason = refused[("pea_isolate_uht_140C_Trikusuma2019", compound)]
        assert wave.lower() in reason.lower(), (compound, reason[:120])
    # and the one added row that IS answerable is answered, on the lane that carries it
    bench = next(b for b in scores["benchmarks"] if b["benchmark_id"] == "pea_isolate_uht_140C_Trikusuma2019")
    dec = next(c for c in bench["compounds"] if c["compound"] == "(E,E)-2,4-decadienal")
    assert dec["lane"] == "lipid" and dec["within_band"]
