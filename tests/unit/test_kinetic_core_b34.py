"""Wave B34 (2026-09-11): the five observables already on disk, and the unit guard that should
have existed two waves ago (kinetic_core_b34_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import engine, panel
from src.kinetic_core.engine import TARGET_ALIASES, _REPORTED_IN_MMOL_PER_L
from src.kinetic_core.species_lipid import MOLECULAR_WEIGHT_G_PER_MOL as LIPID_MW
from src.kinetic_core.species_sulfur import MOLECULAR_WEIGHT_G_PER_MOL as MW

BUNDLE = (data_paths.BENCHMARKS_DIR / "external_validation" / "maillard_path"
          / "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021.json")
SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"


def test_every_reachable_target_can_be_weighed_or_is_a_declared_pool():
    """THE GUARD THAT WAS MISSING, and the reason this wave exists.

    `engine._concentrations` used to fall back to reporting mmol/L for any species with no molar
    mass. That is not a different unit, it is a wrong number: it reads as a prediction smaller than
    the truth by exactly the molar mass. It cost B28 a day (2-pentylfuran, diagnosed as a routing
    problem), and it was still live three days later when B34 asked for 3-deoxyglucosone (180 144x)
    and methylglyoxal (92 306x). A name may be reachable as a target only if it can be weighed, or
    if it is declared to be an accounting pool rather than a molecule.
    """
    unweighable = sorted(
        set(TARGET_ALIASES.values()) - set(MW) - set(LIPID_MW) - set(_REPORTED_IN_MMOL_PER_L) - {"ACR"}
    )
    assert unweighable == [], (
        "these names can be asked for and cannot be weighed, so they would be reported in mmol/L: "
        f"{unweighable}"
    )


def test_a_species_with_no_molar_mass_now_raises_instead_of_changing_the_unit():
    assert "MEL_N" in _REPORTED_IN_MMOL_PER_L        # a pool: mmol/L on purpose
    assert MW["TDG"] == 162.14 and MW["MGO"] == 72.06 and MW["ODG"] == 162.14
    # TDG and INT share a formula and are different molecules; both must be present and equal.
    assert MW["TDG"] == MW["INT"] == 162.14


def test_the_bundle_cites_the_right_authors_and_says_the_pdf_is_on_disk():
    bench = json.loads(BUNDLE.read_text())
    sm = bench["source_metadata"]
    assert sm["citation"].startswith("Leitzen, S.")
    assert "Steinhagen" not in sm["citation"]
    assert bench["source_doi"] == "10.3390/ph14111121"
    assert "Leitzen2021.pdf" in sm["pdf_on_disk"]
    # The id is deliberately NOT renamed; it must say so in its own provenance.
    assert bench["benchmark_id"].endswith("Steinhagen2021")
    assert "NOT renamed" in sm["citation_correction_2026_09_11"]
    # The buffer note's old "SOURCE PAPER NOT ON DISK" claim is kept as the audit record, and must
    # be LABELLED superseded rather than deleted -- the repo's standing practice for a corrected note.
    note = bench["conditions"]["buffer"]["provenance_note"]
    assert note.startswith("WAVE B34") and "THE SOURCE IS ON DISK" in note
    assert "SUPERSEDED 2026-09-11 AND RETAINED AS THE AUDIT RECORD" in note
    assert note.index("THE SOURCE IS ON DISK") < note.index("SUPERSEDED")


@pytest.mark.parametrize("compound,fold", [
    # RE-PINNED BY WAVE B41 (2026-09-11): 1.11 -> 1.27, 1.28 -> 1.01, 11.93 -> 9.31, 32.43 -> 6.63. The fed
    # 3-deoxy triangle and the formic-acid exit's pH term, fitted on pots this bundle is not in.
    ("3-deoxyglucosone", 1.27), ("methylglyoxal", 1.01), ("5-Hydroxymethylfurfural (HMF)", 9.31),
    ("3,4-dideoxyglucosone", 6.63), ("glyoxal", 34.83), ("glucosone", 62.99),
])
def test_the_six_observables_score_where_the_prereg_said(compound, fold):
    scores = json.loads(SCORES.read_text())
    bench = next(b for b in scores["benchmarks"] if b["benchmark_id"].endswith("Steinhagen2021"))
    row = next(c for c in bench["compounds"] if c["compound"] == compound)
    assert row["fold_error"] == pytest.approx(fold, rel=0.01)


def test_the_entry_is_right_and_the_step_after_it_is_not():
    """B33 IS REFUTED BY THIS ROW, and the deficit is localised instead.

    B33 forecast that the amine-free sugar entries are ~10x too slow in water. The entry lands at
    1.11x and the methylglyoxal route at 1.28x; it is the NEXT step, 3-DG -> 3,4-DGE, that is 32x
    low. Scaling the entries as B33 proposed would have broken two right answers to fix one wrong.
    """
    scores = json.loads(SCORES.read_text())
    bench = next(b for b in scores["benchmarks"] if b["benchmark_id"].endswith("Steinhagen2021"))
    f = {c["compound"]: c["fold_error"] for c in bench["compounds"]}
    # At B34 the step after the entry was 32x low. B41 (2026-09-11) fitted it on fed pots and it is now 6.6x:
    # still outside threefold, no longer the deficit this test was written about.
    assert f["3-deoxyglucosone"] < 1.5 and 3.0 < f["3,4-dideoxyglucosone"] < 10.0
    assert f["methylglyoxal"] < 1.5
    # and the two that are structurally absent in an amine-free pot, with a measured size
    assert f["glucosone"] > 30.0 and f["glyoxal"] > 30.0
