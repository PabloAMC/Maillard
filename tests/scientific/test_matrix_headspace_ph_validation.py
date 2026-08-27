import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark
from src.headspace import HeadspaceModel


# RE-PINNED 2026-08-27 (Wave J2). This constant is the model's OWN knob, not a measurement.
# src/headspace.py:241 computes factor = exp(log_slope * (reference_ph - pH)), clipped to
# [min_factor, max_factor]. With the shipped log_slope = 0.235 the pH 4.5 / pH 6.5 ratio is
# exp(0.235 * (6.0 - 4.5)) / exp(0.235 * (6.0 - 6.5)) = exp(0.235 * 2.0) = 1.6000 exactly,
# for every acid-sensitive compound and independent of its concentration.
_EXPECTED_PH_RELEASE_RATIO = 1.6000


def test_ph_release_surrogate_reproduces_its_own_log_slope_knob():
    """BEHAVIOURAL REGRESSION on a tuned knob. NOT a literature validation. See below.

    RENAMED AND RE-PINNED 2026-08-27 (Wave J2, red-team finding S5). The previous name was
    ``test_pouvreau_pea_ph_family_is_reproduced_as_an_acidic_release_trend`` and its
    docstring claimed:

        "Pouvreau 2021 reports roughly 55-65% higher headspace release at pH 4.5 vs pH 6.5
        for pea-isolate aldehydes/furans."

    THAT CITATION DOES NOT EXIST. Verified 2026-08-27 against CrossRef (author search,
    Pouvreau, 2018-2024) and against this repo's own
    ``results/validation/citation_verification_ledger.md``, where "Pouvreau" appears ZERO
    times -- it was never DOI-verified because it has no DOI. The claim is therefore removed
    rather than adjusted. Three separate defects were stacked in it:

    1. GHOST ATTRIBUTION. No 2021 Pouvreau paper on pea-protein aroma, headspace or pH
       exists. ``data/benchmarks/maillard_validation_benchmarks.md`` still carries it as
       PRIMARY with the self-declaring placeholder DOI
       "(Approx. DOI: 10.1016/j.foodchem.2021.xxx -- retrieve via Scopus)".
       The real underlying source, already repaired in ``data/lit`` on 2026-08-26 under
       ``doi_repair`` (record ``karolkowski_2021_ppi_hexanal_ph_release``, itself a
       misattribution), is:
           Fischer, E.; Cachon, R.; Cayot, N. "Effects of extraction pH on the volatile
           compounds from pea protein isolate." Food Res. Int. 2021, 150, 110760.
           DOI: 10.1016/j.foodres.2021.110760
    2. OVER-GENERALISATION. Fischer et al. report 59% -- for HEXANAL ONLY, and explicitly as
       "for example". Nothing in it supports a per-compound figure for nonanal or
       2-pentylfuran, so 2 of the 3 compounds this test asserted on had no quantitative
       backing at all. The repo's own record concedes this: its ``numeric_reference`` units
       are ``qualitative_release_order``, i.e. an ORDERING, not a percentage.
    3. WRONG PHYSICAL QUANTITY. Fischer varied the pH at which the isolate was EXTRACTED
       (a manufacturing variable that changes which volatiles the powder carries).
       ``predict_headspace(..., pH=...)`` varies pH as a RELEASE parameter at measurement
       time on a fixed matrix composition. These are not the same measurement.

    And the 1.5-1.7 band was not a literature band. exp(2 * 0.235) = 1.6000 sits dead-centre
    of it, so the test was pinning ``log_slope`` against itself and reporting the result as
    literature agreement. The three per-compound assertions were also the SAME assertion
    three times over: the factor in src/headspace.py depends only on pH and the knob, never
    on the compound's concentration, so the 205.0 / 75.0 / 310.0 inputs below could be any
    numbers at all without changing the outcome.

    WHAT THIS TEST NOW IS: a behavioural regression pin on an UNANCHORED surrogate. It
    catches a silent change to log_slope, to the clip bounds, or to the acid-sensitivity
    scope gate. It asserts nothing about the real world. Treat the knob as
    ``no_verifiable_source`` for any quantitative purpose: the underlying literature
    supports the DIRECTION (acidic release exceeds near-neutral release) and one 59% figure
    for one compound under a different experimental variable -- it does not support this
    magnitude, and no measurement in this repo constrains it.
    """
    model = HeadspaceModel()
    # These concentrations are inert: the pH factor is concentration-independent by
    # construction. They exist only so predict_headspace has a matrix to act on.
    matrix = {
        "Hexanal": 205.0,
        "Nonanal": 75.0,
        "2-Pentylfuran": 310.0,
        "2,5-Dimethylpyrazine": 10.0,
    }

    air_acid = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=4.5)
    air_less_acid = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=6.5)

    for compound in ["Hexanal", "Nonanal", "2-Pentylfuran"]:
        ratio = air_acid[compound] / air_less_acid[compound]
        # Tightened from the old 1.5-1.7 band to the knob's exact value: a band that wide
        # around a deterministic closed-form ratio only hides drift. If this fails, the
        # surrogate's parameters changed -- fix the pin and say why, do not widen it.
        assert ratio == pytest.approx(_EXPECTED_PH_RELEASE_RATIO, rel=1e-4), (
            f"{compound} acidic/near-neutral release ratio is {ratio:.4f}, not the "
            f"{_EXPECTED_PH_RELEASE_RATIO:.4f} implied by src/headspace.py's shipped "
            f"log_slope. The unanchored pH surrogate changed."
        )

    # Scope gate, not a chemistry result: src/headspace.py:229 selects acid-sensitive
    # compounds by the substring test ("anal" | "enal" | "furan"), which no pyrazine name
    # can match, so this is exactly 1.0 by construction. It is kept as a guard on the SCOPE
    # of the surrogate -- widening that gate to all volatiles would go red here -- and is
    # explicitly not evidence that pyrazine release is pH-independent in reality.
    pyrazine_ratio = air_acid["2,5-Dimethylpyrazine"] / air_less_acid["2,5-Dimethylpyrazine"]
    assert pyrazine_ratio == pytest.approx(1.0, rel=1e-9)


# RE-PINNED 2026-08-27 (Wave M). The old assertion was a blanket `ratio <= 1.01` on every
# compound of both Pratap-Singh benchmarks, which read as an accuracy claim but was really
# the signature of ALGEBRAIC RECOVERY: each compound's matrix observability factor had been
# back-solved from that compound's own measured ppb (see the comment block at the top of
# src/matrix_calibration_registry.py).
#
# The Wave K/M content correction broke that in the most informative possible way. Molecules
# 2021, 26, 4104 Table 1 (Europe PMC, PMC8271896) reports hexanal = 1138.00 ppb (pea) and
# 1621.71 ppb (soy); the repo carried 260 / 380. The 2-pentylfuran values (638 / 2492) were
# verified verbatim and are unchanged. So the factors still recover 2-pentylfuran exactly and
# now miss hexanal by exactly the size of the transcription error -- 1138/260.65 = 4.366x and
# 1621.71/379.88 = 4.269x. The observability factors were deliberately NOT refitted: that is
# an owner science decision (AUDIT.md, Round 3), and doing it in the same pass as the Wave N
# chemistry change would make any resulting agreement unattributable.
#
# The test is therefore split per compound rather than relaxed. The <=1.01 pin survives where
# it still means something (2-pentylfuran, still a genuine recovery), and hexanal is pinned
# two-sided at its honest miss so a silent drift in either direction is caught.
_PRATAP_SINGH_EXPECTED_RATIOS = {
    "pea_isolate_40C_PratapSingh2021": {"hexanal": 4.366, "2-pentylfuran": 1.0},
    "soy_isolate_40C_PratapSingh2021": {"hexanal": 4.269, "2-pentylfuran": 1.0},
}


def test_matrix_only_pratap_singh_baselines_remain_stable_at_reference_ph():
    pea = evaluate_benchmark(ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json")
    soy = evaluate_benchmark(ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json")

    for evaluation in [pea, soy]:
        assert evaluation.supported is True
        expected = _PRATAP_SINGH_EXPECTED_RATIOS[evaluation.benchmark_id]
        assert {c.compound for c in evaluation.comparisons} == set(expected), (
            f"{evaluation.benchmark_id} compound set changed; the hexanol rows were removed "
            "2026-08-27 because the paper reports n.d. for both matrices."
        )
        for comparison in evaluation.comparisons:
            if comparison.compound == "2-pentylfuran":
                # Still an exact algebraic recovery: this value WAS verified against the
                # paper, so the constant solved from it still reproduces it.
                assert comparison.ratio <= 1.01, (
                    f"{evaluation.benchmark_id} 2-pentylfuran recovery broke "
                    f"({comparison.ratio:.4f}x); its factor is solved from a verified value."
                )
            else:
                assert comparison.predicted_ppb < comparison.measured_ppb
                assert comparison.ratio == pytest.approx(
                    expected[comparison.compound], rel=0.01
                ), (
                    f"{evaluation.benchmark_id} {comparison.compound} ratio is "
                    f"{comparison.ratio:.3f}x, pinned at "
                    f"{expected[comparison.compound]}x -- the exact size of the 2026-08-27 "
                    "Wave K/M reference correction. If this moved, check whether the "
                    "observability factor was refitted (an owner decision that must be "
                    "recorded in the registry and AUDIT.md) or whether the benchmark was "
                    "reverted to its erroneous values."
                )