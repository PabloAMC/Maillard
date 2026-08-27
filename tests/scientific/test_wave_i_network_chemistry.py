"""Wave I (2026-08-27) chemistry remediation regressions — scientific tier.

These exercise the whole enumeration and the benchmark prediction path, so they
are slower than the unit-tier file (tests/unit/test_wave_i_chemistry.py) that
covers the same wave's SMARTS, single-template and file-content assertions.

Every number below was MEASURED on the tree at the end of Wave I, not chosen.
The measurement conditions are stated with each pin, because a share or a ratio
is meaningless without the system it was measured on.
"""

import json
import sys
from collections import Counter
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.smirks_engine as smirks_engine  # noqa: E402
from src.benchmark_validation import evaluate_benchmark_payload  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from src.precursor_resolver import resolve_many  # noqa: E402
from src.smirks_engine import SmirksEngine  # noqa: E402


# ── the "core (no-lipid) network", fixed once so every share below is comparable

_CORE_PRECURSORS = [
    "D-Glucose",
    "D-Ribose",
    "Glycine",
    "L-Cysteine",
    "L-Lysine",
    "L-Asparagine",
]
_CORE_CONDITIONS = dict(pH=5.5, temperature_celsius=150.0, water_activity=0.95)
_CORE_GENERATIONS = 4


@pytest.fixture(scope="module")
def core_network():
    conditions = ReactionConditions(**_CORE_CONDITIONS)
    steps = SmirksEngine(conditions).enumerate(
        resolve_many(_CORE_PRECURSORS), max_generations=_CORE_GENERATIONS
    )
    return steps, Counter(s.reaction_family or "unknown" for s in steps)


# ── the ordinal pentose/hexose harness, matched to
#    tests/scientific/test_pentose_hexose_sulfur_ordering.py

_MATCHED_CONDITIONS = {
    "temp_C": 150.0,
    "ph": 5.5,
    "water_activity": 0.95,
    "time_min": 60.0,
}


def _ordinal_payload(benchmark_id: str, sugar: str) -> dict:
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
                "Synthetic formulation for a Wave I regression; the conc_ppb placeholder "
                "is NOT a reference value and nothing here scores against it."
            ),
        },
        "validation_contract": {
            "scale_thresholds": {"max_ratio": 1e9, "mean_abs_log10_error": 9.0}
        },
        "measured_volatiles": {
            "2-methyl-3-furanthiol": {"conc_ppb": 1.0, "uncertainty_pct": 100}
        },
    }


def _mft_ppb(evaluation) -> float:
    by_compound = {c.compound: c for c in evaluation.comparisons}
    return by_compound["2-methyl-3-furanthiol"].predicted_ppb


# ══ FIX 8 — MFT must not be a downstream dependent of pyrazine chemistry ═══

def test_mft_does_not_collapse_when_aminoketone_condensation_is_disabled(monkeypatch):
    """Red-team finding H4, the regression this wave exists for.

    The norfuraneol + H2S -> MFT step consumes the pool species `[HH]`, a 2[H]
    reducing-equivalent token. Before Wave I the token's only producer reachable
    from a ribose/cysteine system was `Aminoketone_Condensation` (dihydropyrazine
    -> pyrazine aromatisation), so switching the pyrazine lane off drove
    predicted MFT to EXACTLY 0.0 ppb -- a coupling nothing in the literature
    supports. `_thiol_reductant_pool` (2 cysteine -> cystine + 2[H]) gives the
    token a source that is a PRECURSOR of every system that can make MFT at all.

    MEASURED, ribose + cysteine, 150 C / pH 5.5 / aw 0.95 / 60 min:
        before Wave I, aminoketone condensation disabled:    0.00 ppb
        after  Wave I, aminoketone condensation disabled: 1038.80 ppb
    """
    monkeypatch.setattr(
        smirks_engine, "_aminoketone_condensation", lambda *args, **kwargs: []
    )
    evaluation = evaluate_benchmark_payload(
        _ordinal_payload("wave_i_mft_without_pyrazines", "D-Ribose")
    )
    assert evaluation.supported, evaluation.reason

    mft = _mft_ppb(evaluation)
    assert mft > 0.0, (
        "MFT collapsed to zero with the pyrazine lane disabled: the 2[H] token is "
        "coupled to aminoketone condensation again (red-team H4)."
    )
    # Order-of-magnitude pin, not a calibration: the point is that MFT survives
    # at a comparable level, not that it lands on a particular number.
    assert mft > 100.0, f"MFT is only {mft:.2f} ppb without the pyrazine lane"


def test_reducing_equivalent_token_has_a_non_pyrazine_producer(core_network):
    """The structural half of the same finding: at least one `[HH]` producer in
    a cysteine system must not be `Aminoketone_Condensation`."""
    steps, _families = core_network
    producers = {
        step.reaction_family
        for step in steps
        if any(p.smiles == "[HH]" for p in step.products)
    }
    assert producers, "nothing produces the 2[H] reducing-equivalent token at all"
    assert producers - {"Aminoketone_Condensation"}, (
        "the only producer of the 2[H] token is pyrazine chemistry again: "
        f"{sorted(producers)}"
    )
    assert "Thiol_Oxidation" in producers, (
        "the cysteine -> cystine reductant step (_thiol_reductant_pool) is gone"
    )


# ══ FIX 9 — Lipid_Schiff_Base enumeration share ════════════════════════════

# MEASURED on the core network defined above (no lipid precursor at all):
#   pre-Wave-I:            224 Lipid_Schiff_Base steps of 298 total = 75.2%
#   after fix 8 landed:    228 of 303                              = 75.2%
#   after fix 9 (SMARTS):   28 of 103                              = 27.2%
#   end of Wave I:          28 of 102                              = 27.5%
# The ceiling is set at 0.35 -- comfortably above the measured 0.275, so normal
# network growth does not trip it, and far below the 0.75 it was, so the defect
# cannot regrow. It is a RATCHET, not a target: if a later wave gets the share
# down further, lower this number rather than leaving the slack.
_LIPID_SCHIFF_SHARE_CEILING = 0.35


def test_lipid_schiff_base_is_not_the_majority_of_the_core_network(core_network):
    """Red-team finding H3. A family named for LIPID chemistry was 3 steps in
    every 4 of a network with no lipid in it."""
    steps, families = core_network
    total = len(steps)
    assert total > 0
    count = families.get("Lipid_Schiff_Base", 0)
    share = count / total
    assert share < _LIPID_SCHIFF_SHARE_CEILING, (
        f"Lipid_Schiff_Base is {count}/{total} = {share:.1%} of the core no-lipid "
        f"network, above the {_LIPID_SCHIFF_SHARE_CEILING:.0%} ceiling. It was 75.2% "
        "before Wave I fix 9 because the rule's SMARTS did not implement the "
        "alpha-hydroxy and C3+ exclusions its own comment claimed."
    )


def test_lipid_schiff_base_no_longer_consumes_sugar_derived_carbonyls(core_network):
    """The exclusion the comment always claimed: "excludes sugars"."""
    steps, _families = core_network
    offenders = []
    for step in steps:
        if step.reaction_family != "Lipid_Schiff_Base":
            continue
        for reactant in step.reactants:
            # alpha-hydroxy aldehydes: an OH on the carbon next to the CHO.
            if reactant.smiles in ("OCC=O", "O=CCO", "O=CC(O)CO"):
                offenders.append((reactant.label, reactant.smiles))
    assert not offenders, (
        f"Lipid_Schiff_Base is condensing sugar-derived alpha-hydroxy carbonyls: {offenders}"
    )


# ══ FIX 12 — the demoted MFT shortcut is hexose-only in a real enumeration ══

def test_pentose_system_does_not_use_the_legacy_mft_shortcut():
    conditions = ReactionConditions(**_CORE_CONDITIONS)
    steps = SmirksEngine(conditions).enumerate(
        resolve_many(["D-Ribose", "L-Cysteine"]), max_generations=_CORE_GENERATIONS
    )
    families = {s.reaction_family for s in steps}
    assert "Thiol_Addition_Legacy_Shortcut" not in families, (
        "the demoted 3-deoxyosone lump fired on a pentose again; pentoses reach "
        "MFT through the isotope-evidenced 1,4-dideoxyosone route"
    )
    # 2026-08-27 (Wave N) — RE-PINNED from Thiol_Addition_Norfuraneol. CAUSE:
    # route correction on isotope evidence (Cerny & Davidek 2003,
    # 10.1021/jf026123f; 2004, 10.1021/jf035265m).
    assert "Thiol_Addition_Pentodiulose" in families, (
        "the accepted pentose route disappeared along with the lump"
    )
    assert "Thiol_Addition_Norfuraneol" not in families, (
        "the retired norfuraneol->MFT step re-appeared"
    )


def test_hexose_system_keeps_its_lumped_path_to_mft():
    """Why fix 12 RESTRICTED rather than DELETED the lump."""
    conditions = ReactionConditions(**_CORE_CONDITIONS)
    steps = SmirksEngine(conditions).enumerate(
        resolve_many(["D-Glucose", "L-Cysteine"]), max_generations=_CORE_GENERATIONS
    )
    families = {s.reaction_family for s in steps}
    assert "Thiol_Addition_Hexose_Legacy_Shortcut" in families
    labels = {p.label for s in steps for p in s.products}
    assert "2-methyl-3-furanthiol" in labels, "MFT became unreachable from glucose"


def test_pentose_still_far_exceeds_hexose_as_an_mft_precursor():
    """Hofmann & Schieberle 1998 (10.1021/jf9705983): pentoses give *much* more
    MFT than hexoses. Wave I moved both sides; the ordering must survive.

    MEASURED at 150 C / pH 5.5 / aw 0.95 / 60 min, ribose vs glucose:
        pre-Wave-I:  315.26 / 44.74 ppb = 7.05x
        end of Wave I: 981.31 / 109.33 ppb = 8.98x
    The 3x floor is the same deliberately coarse reading of "much" used by
    tests/scientific/test_pentose_hexose_sulfur_ordering.py.
    """
    pentose = evaluate_benchmark_payload(_ordinal_payload("wave_i_rib", "D-Ribose"))
    hexose = evaluate_benchmark_payload(_ordinal_payload("wave_i_glc", "D-Glucose"))
    assert pentose.supported and hexose.supported

    pentose_mft, hexose_mft = _mft_ppb(pentose), _mft_ppb(hexose)
    assert hexose_mft > 0.0, "degenerate comparison: hexose MFT is zero"
    ratio = pentose_mft / hexose_mft
    assert ratio >= 3.0, (
        f"ribose {pentose_mft:.2f} ppb vs glucose {hexose_mft:.2f} ppb = {ratio:.2f}x"
    )


# ══ The flagship benchmark, recorded rather than asserted-into-agreement ════

def test_hofmann1998_sulfur_predictions_are_recorded_not_tuned():
    """cys_ribose_140C_Hofmann1998 is the ONLY surviving literature constraint on
    the sulfur branch, so what Wave I did to it must be visible.

    MEASURED (140 C / pH 5.0 / aw 0.98 / 30 min; reference MFT 342 ppb,
    FFT 200 ppb):
        pre-Wave-I        MFT  61.25 ppb (5.58x UNDER)   FFT  61.44 (3.26x UNDER)
        after fix 8       MFT 345.04 ppb (1.01x)         FFT 187.49 (1.07x UNDER)
        end of Wave I     MFT 235.32 ppb (1.45x UNDER)   FFT 219.96 (1.10x OVER)

    NOT A CALIBRATION. No barrier was touched in this wave; the movement is a
    consequence of decoupling the 2[H] token from pyrazine chemistry (fix 8) and
    of removing a lumped pentose shortcut that was double-producing MFT (fix 12).
    `thiol_addition_norfuraneol` is still the 26.85 that Wave H fitted against a
    network that no longer exists, so the remaining agreement is a coincidence of
    two independent choices -- see the Wave I note in src/barrier_constants.py.
    The bounds below are wide on purpose: they exist to catch a SILENT collapse
    or a silent blow-up, not to assert accuracy.
    """
    payload = json.loads(
        (ROOT / "data/benchmarks/cys_ribose_140C_Hofmann1998.json").read_text(
            encoding="utf-8"
        )
    )
    evaluation = evaluate_benchmark_payload(payload)
    assert evaluation.supported, evaluation.reason

    predicted = {c.compound: c.predicted_ppb for c in evaluation.comparisons}
    mft = predicted["2-methyl-3-furanthiol"]
    fft = predicted["2-furfurylthiol"]

    assert 100.0 < mft < 700.0, f"Hofmann1998 MFT = {mft:.2f} ppb (reference 342)"
    assert 80.0 < fft < 500.0, f"Hofmann1998 FFT = {fft:.2f} ppb (reference 200)"


# ══ FIX 18 — acrylamide in a real enumeration ══════════════════════════════

def _acrylamide_payload(benchmark_id: str, sugar: str, sugar_mM: float) -> dict:
    return {
        "benchmark_id": benchmark_id,
        "source_doi": "synthetic_wave_i_fix18",
        "precursors": {
            "L-Asparagine": {"concentration_mM": 10.0},
            sugar: {"concentration_mM": sugar_mM},
        },
        "conditions": {
            "temp_C": 160.0,
            "ph": 6.0,
            "water_activity": 0.6,
            "time_min": 20.0,
        },
        "metadata": {
            "tier": "DIAGNOSTIC",
            "family": "safety",
            "execution_path": "free_precursor",
            "notes": "Synthetic Wave I regression; the placeholder is not a reference value.",
        },
        "validation_contract": {
            "scale_thresholds": {"max_ratio": 1e9, "mean_abs_log10_error": 9.0}
        },
        "measured_volatiles": {"acrylamide": {"conc_ppb": 1.0, "uncertainty_pct": 100}},
    }


def test_predicted_acrylamide_scales_with_reducing_sugar_dose():
    """`src/safety.py::predict_acrylamide` is second order in [Asn]*[Sugar]; the
    template must not be able to make acrylamide out of a catalytic sugar.

    MEASURED (Asn 10 mM, glucose 1 / 10 / 100 mM, 160 C / 20 min / pH 6.0 /
    aw 0.6): 7.86 / 78.61 / 786.09 ppb -- exactly linear.
    """
    values = [
        {c.compound: c.predicted_ppb for c in
         evaluate_benchmark_payload(
             _acrylamide_payload(f"wave_i_acr_{dose}", "D-Glucose", dose)
         ).comparisons}["acrylamide"]
        for dose in (1.0, 10.0, 100.0)
    ]
    assert values[0] > 0.0
    assert values[1] > values[0] and values[2] > values[1], values
    assert values[1] / values[0] == pytest.approx(10.0, rel=0.05), values
    assert values[2] / values[1] == pytest.approx(10.0, rel=0.05), values


def test_acrylamide_step_consumes_the_sugar_in_a_real_enumeration(core_network):
    steps, _families = core_network
    acrylamide_steps = [
        s for s in steps if s.reaction_family == "Safety_Risk_Acrylamide"
    ]
    assert acrylamide_steps, "the acrylamide lane vanished from the core network"
    for step in acrylamide_steps:
        reactant_smiles = {r.smiles for r in step.reactants}
        product_smiles = {p.smiles for p in step.products}
        assert not (reactant_smiles & product_smiles), (
            "the sugar (or asparagine) is regenerated by the acrylamide step: "
            f"{sorted(reactant_smiles & product_smiles)}"
        )
