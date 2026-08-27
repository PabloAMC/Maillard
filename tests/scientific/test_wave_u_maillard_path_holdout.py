"""Wave U (2026-08-27) — the reaction network's FIRST out-of-sample test.

Wave S1 established that all four bundles in the legacy external hold-out run the
``matrix_only`` execution path and never reach ``Recommender.predict_from_steps``.
Three consecutive waves of reaction-network work moved in-panel predictions and
left all eight hold-out points bit-identical *for that structural reason*. The
chemistry this repository is about had never been scored on a system it had not
already seen.

Wave U harvested twelve content-verified free-precursor literature points into
``data/benchmarks/external_validation/maillard_path/`` and froze the CURRENT
model's predictions on them into
``results/validation/maillard_path_holdout_frozen_predictions.{json,md}`` BEFORE
any rate-calibration wave ran. That artifact is a PRE-REGISTRATION.

These tests guard the three things that make it one:

1. every bundle really does execute the network (``free_precursor``), because a
   bundle that silently drifted to ``matrix_only`` would stop testing the network
   while still appearing in the panel;
2. every bundle is labelled and located so ``scripts/ci/holdout_guard.py`` reaches
   it, and no bundle is named by any fit record;
3. the frozen artifact still says which git HEAD it was generated from, so a later
   wave cannot quietly regenerate it and compare it to itself.

They deliberately do NOT pin the scores. The scores are the thing under test in a
later wave; pinning them here would create exactly the "assert the number you just
produced" circularity this campaign spent three rounds removing. The one number
pinned is the bundle COUNT, so that silently dropping an inconvenient point breaks
the build.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
HOLDOUT_DIR = ROOT / "data" / "benchmarks" / "external_validation" / "maillard_path"
FROZEN_JSON = ROOT / "results" / "validation" / "maillard_path_holdout_frozen_predictions.json"

EXPECTED_BUNDLE_COUNT = 12


def _bundles() -> list[tuple[Path, dict]]:
    return [
        (path, json.loads(path.read_text(encoding="utf-8")))
        for path in sorted(HOLDOUT_DIR.glob("*.json"))
    ]


def test_the_maillard_path_holdout_has_not_shrunk():
    """A hold-out that loses its inconvenient points is not a hold-out.

    Wave U ingested 12 bundles. Dropping one is a defensible scientific act and an
    indefensible silent one, so it has to come with a deliberate edit here.
    """
    paths = sorted(HOLDOUT_DIR.glob("*.json"))
    assert len(paths) == EXPECTED_BUNDLE_COUNT, (
        f"{len(paths)} Maillard-path hold-out bundles found, Wave U froze "
        f"{EXPECTED_BUNDLE_COUNT}. If a point was removed for a stated reason, update this "
        "count in the same commit and record the reason in tasks/audit_remediation.md."
    )


@pytest.mark.parametrize("path,payload", _bundles(), ids=lambda item: getattr(item, "name", ""))
def test_every_bundle_executes_the_reaction_network(path: Path, payload: dict):
    """``free_precursor`` is the whole point: it is the path that calls predict_from_steps.

    ``matrix_only`` reads a lipid-oxidation load off a matrix profile and returns
    before the network is touched. A bundle that drifted to it would keep scoring,
    keep looking fine, and stop testing anything this wave exists to test.
    """
    assert payload.get("metadata", {}).get("execution_path") == "free_precursor", (
        f"{path.name}: execution_path is "
        f"{payload.get('metadata', {}).get('execution_path')!r}, must be 'free_precursor'."
    )
    assert payload.get("protein_type") == "free", (
        f"{path.name}: protein_type is {payload.get('protein_type')!r}. A matrix protein_type "
        "routes the prediction through the matrix observability layer instead of the network."
    )


@pytest.mark.parametrize("path,payload", _bundles(), ids=lambda item: getattr(item, "name", ""))
def test_every_bundle_is_labelled_so_the_guard_reaches_it(path: Path, payload: dict):
    assert payload.get("evidence_class") == "external_validation_only", (
        f"{path.name}: without this label the bundle is eligible for calibration."
    )
    assert payload.get("maillard_path_holdout") is True, (
        f"{path.name}: the flag is what holdout_guard check 4 keys on."
    )


@pytest.mark.parametrize("path,payload", _bundles(), ids=lambda item: getattr(item, "name", ""))
def test_every_value_carries_a_retrieved_quote_and_a_stated_basis(path: Path, payload: dict):
    """No value entered this hold-out without a quoted source. Keep it that way.

    This is the rule the 342/200 anchor broke: a number was carried as a literature
    measurement for four months while its only upstream was a table the repository
    had written itself. A target with no quote and no basis is that failure mode
    re-entering through the one lane that is supposed to be clean.
    """
    targets = payload.get("holdout_targets") or {}
    assert targets, f"{path.name}: no holdout_targets block."
    for compound, spec in targets.items():
        assert spec.get("target_quote"), f"{path.name} / {compound}: no target_quote."
        assert spec.get("target_value") is not None, f"{path.name} / {compound}: no target_value."
        assert spec.get("molar_basis") or spec.get("native_unit_and_conversion"), (
            f"{path.name} / {compound}: neither molar_basis nor native_unit_and_conversion. "
            "A number without a stated denominator is how an assumed basis becomes a free "
            "multiplicative parameter."
        )
    verification = payload.get("content_verification") or {}
    assert verification.get("access_route"), f"{path.name}: no access_route recorded."
    assert verification.get("retrieval_date"), f"{path.name}: no retrieval_date recorded."


def test_the_frozen_artifact_names_the_commit_it_was_generated_from():
    """Without the commit, the pre-registration is just a report.

    The artifact's whole function is that a later calibration wave is scored against
    predictions made at a NAMED earlier state. If the git block goes missing, the
    file can no longer prove which model produced its numbers.
    """
    if not FROZEN_JSON.is_file():
        pytest.skip(
            "results/validation/ is gitignored; regenerate with "
            "scripts/generators/generate_maillard_path_holdout_frozen_predictions.py"
        )
    payload = json.loads(FROZEN_JSON.read_text(encoding="utf-8"))
    git = payload.get("git") or {}
    assert git.get("commit") and git["commit"] != "unknown", "frozen artifact has no git commit"
    assert payload.get("pre_registration"), "frozen artifact lost its pre-registration statement"
    assert payload["summary"]["bundle_count"] == EXPECTED_BUNDLE_COUNT
    # Structural zeroes must stay visible. They have no finite fold error, so a
    # median that silently drops them reports a better model than exists.
    assert "structural_zero_count" in payload["summary"], (
        "the summary no longer counts structural zeroes; a compound the model returns "
        "NOTHING for would then disappear from the report rather than failing in it."
    )
