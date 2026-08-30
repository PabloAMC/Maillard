"""Ordinal sulfur-branch regression: pentose >> hexose as an MFT precursor.

This file replaces two tests that were deleted on 2026-08-26 together with their fixtures:

  * ``tests/scientific/test_mottram_coverage.py`` (fixture ``cys_ribose_150C_Mottram1994``)
  * ``tests/scientific/test_farmer_coverage.py`` (fixture ``cys_glucose_150C_Farmer1999``)

Both fixtures were fabricated -- source recovery established that no literature source
exists for either (full findings in ``data/benchmarks/quarantined/README.md``), so the files
were deleted rather than left quarantined. Their *reference values* were worthless, but the
assertions built on them were about reaction-network / SMIRKS coverage, which is real
regression value. That value is preserved here on synthetic formulations, so nothing in this
file depends on an unverifiable number.

The quantitative claim those files pretended to make is not recoverable from the literature:
the quantitative cysteine + sugar corpus of that era reports FD factors or mol % yields, never
absolute ppb triples. What the literature *does* support is an **ordinal** constraint, and that
is what is asserted below.

Literature anchor
-----------------
Hofmann, T.; Schieberle, P. "Identification of Potent Aroma Compounds in Thermally Treated
Mixtures of Glucose/Cysteine and Rhamnose/Cysteine Using Aroma Extract Dilution Techniques."
*J. Agric. Food Chem.* **1998**, 46, 235-241. DOI: ``10.1021/jf9705983``

Its abstract states verbatim that *"pentoses generated much higher amounts of FFT and MFT than
hexoses when heated in the presence of cysteine"*. 2-Methyl-3-furanthiol is a pentose marker:
it derives from the 4-hydroxy-5-methyl-3(2H)-furanone / norfuraneol branch reachable from a
pentose, whereas a hexose routes preferentially to HMF and 5-methylfurfural. The same paper is
why ``cys_glucose_150C_Farmer1999``'s headline "MFT from glucose" number was judged invented
rather than merely unlocatable.

The margin below is deliberately coarse (a factor, not a fitted ratio). The claim under test is
the *ordering and its scale separation*, which is what the source supports; it is not a
calibration and must never be reported as one.
"""

import sys
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark_payload  # noqa: E402


# Matched conditions for both sugars. Only the sugar identity differs between the two runs.
_MATCHED_CONDITIONS = {
    "temp_C": 150.0,
    "ph": 5.5,
    "water_activity": 0.95,
    "time_min": 60.0,
}

# The pentose system is the one where the full cysteine/ribose thiol branch is expected to be
# reachable, so it carries the three-analyte coverage assertion inherited from the deleted
# Mottram test.
_PENTOSE_ANALYTES = [
    "2-methyl-3-furanthiol",
    "bis(2-methyl-3-furyl) disulfide",
    "furfural",
]

# Inherited from the deleted Farmer test. Note that this asserts only that the network can
# REACH these species from a hexose + cysteine system, not that the yields are right -- and in
# particular not that MFT from glucose is chemically favoured, which the test below denies.
_HEXOSE_ANALYTES = [
    "2-methyl-3-furanthiol",
    "furfural",
    "pyrazine",
]

# Ordering margin. Hofmann & Schieberle report "much higher", not a number; 3x is a
# deliberately conservative reading of "much" that still fails loudly if the sugar-identity
# discrimination in the model ever collapses to parity.
_MIN_PENTOSE_HEXOSE_MFT_RATIO = 3.0


def _payload(benchmark_id: str, sugar: str, analytes: list[str]) -> dict:
    return {
        "benchmark_id": benchmark_id,
        "source_doi": "synthetic_ordinal_constraint_10.1021/jf9705983",
        "precursors": {
            "L-Cysteine": {"concentration_mM": 10.0},
            sugar: {"concentration_mM": 10.0},
        },
        "conditions": deepcopy(_MATCHED_CONDITIONS),
        "metadata": {
            "tier": "DIAGNOSTIC",
            "family": "free_aa_sugar_discrimination",
            "execution_path": "free_precursor",
            "notes": (
                "Synthetic formulation for an ordinal literature constraint; the conc_ppb "
                "placeholders below exist only so the evaluator produces comparison rows and "
                "are NOT reference values. Nothing here asserts agreement with them."
            ),
        },
        "validation_contract": {
            "scale_thresholds": {
                # Wide on purpose: this file never scores against the placeholders.
                "max_ratio": 1e9,
                "mean_abs_log10_error": 9.0,
            }
        },
        "measured_volatiles": {
            analyte: {"conc_ppb": 1.0, "uncertainty_pct": 100} for analyte in analytes
        },
    }


def _mft_predicted_ppb(evaluation) -> float:
    by_compound = {comparison.compound: comparison for comparison in evaluation.comparisons}
    return by_compound["2-methyl-3-furanthiol"].predicted_ppb


def test_pentose_cysteine_system_reaches_the_full_thiol_branch():
    """SMIRKS-coverage regression inherited from the deleted test_mottram_coverage.py."""
    evaluation = evaluate_benchmark_payload(
        _payload("synthetic_cys_ribose_coverage", "D-Ribose", _PENTOSE_ANALYTES)
    )

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0

    matched = {c.compound: c.matched_name for c in evaluation.comparisons}
    for analyte in _PENTOSE_ANALYTES:
        assert matched[analyte] is not None, f"cysteine + ribose cannot reach {analyte}"


def test_hexose_cysteine_system_reaches_its_claimed_products():
    """SMIRKS-coverage regression inherited from the deleted test_farmer_coverage.py."""
    evaluation = evaluate_benchmark_payload(
        _payload("synthetic_cys_glucose_coverage", "D-Glucose", _HEXOSE_ANALYTES)
    )

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0

    matched = {c.compound: c.matched_name for c in evaluation.comparisons}
    for analyte in _HEXOSE_ANALYTES:
        assert matched[analyte] is not None, f"cysteine + glucose cannot reach {analyte}"


def test_pentose_derived_mft_far_exceeds_hexose_derived_mft_at_matched_conditions():
    """Hofmann & Schieberle 1998 (10.1021/jf9705983): pentoses give *much* more MFT than hexoses.

    Same cysteine loading, same sugar molarity, same temperature, pH, water activity and time.
    The only difference between the two runs is ribose vs glucose, so any difference in
    predicted MFT is the model's sugar-identity discrimination and nothing else.
    """
    pentose = evaluate_benchmark_payload(
        _payload("synthetic_cys_ribose_ordinal", "D-Ribose", ["2-methyl-3-furanthiol"])
    )
    hexose = evaluate_benchmark_payload(
        _payload("synthetic_cys_glucose_ordinal", "D-Glucose", ["2-methyl-3-furanthiol"])
    )

    assert pentose.supported, pentose.reason
    assert hexose.supported, hexose.reason

    pentose_mft = _mft_predicted_ppb(pentose)
    hexose_mft = _mft_predicted_ppb(hexose)

    assert hexose_mft > 0.0, "degenerate comparison: hexose MFT prediction is zero"
    ratio = pentose_mft / hexose_mft
    assert ratio >= _MIN_PENTOSE_HEXOSE_MFT_RATIO, (
        "Pentose-derived MFT must far exceed hexose-derived MFT at matched conditions "
        f"(10.1021/jf9705983). Got ribose {pentose_mft:.4f} ppb vs glucose {hexose_mft:.4f} ppb "
        f"= {ratio:.2f}x, below the required {_MIN_PENTOSE_HEXOSE_MFT_RATIO:.1f}x."
    )
