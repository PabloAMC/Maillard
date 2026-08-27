"""Regression tests for the Wave I audit-remediation tooling fixes.

Covers:
  * FIX 6  -- `src/family_sensitivity.py` offset-key routing (false zeros for
              schiff / enolisation / cysteine).
  * FIX 16 -- `src/plot_style.py` graceful degradation when LaTeX/dvipng is
              unavailable.
  * FIX 17 -- `src/matrix_calibration_optimizer.py` guardrail validating the
              baseline constants instead of the candidate, and not restoring the
              mutated globals on the success path.
"""

from __future__ import annotations

import ast
import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# FIX 6 -- family sensitivity offset-key routing
# ---------------------------------------------------------------------------


def _engine_family_labels_in_sources() -> set[str]:
    """AST-scan `src/` for the raw `reaction_family` strings the engine emits."""
    labels: set[str] = set()
    for path in sorted((ROOT / "src").glob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:  # pragma: no cover - defensive
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
            if name not in {"ElementaryStep", "_step"}:
                continue
            for keyword in node.keywords:
                if keyword.arg == "reaction_family" and isinstance(keyword.value, ast.Constant):
                    if isinstance(keyword.value.value, str):
                        labels.add(keyword.value.value)
            for arg in node.args:
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    labels.add(arg.value)
    return labels


def test_engine_family_labels_cover_sources():
    """The pinned label list must not fall behind the engine.

    A family whose raw label is missing here would silently route to a no-op
    offset and be reported as a zero sensitivity -- exactly the FIX 6 failure.
    """
    from src.family_sensitivity import ENGINE_FAMILY_LABELS

    missing = _engine_family_labels_in_sources() - set(ENGINE_FAMILY_LABELS)
    assert not missing, (
        "New reaction_family labels found in src/ that are not pinned in "
        f"ENGINE_FAMILY_LABELS: {sorted(missing)}"
    )


@pytest.mark.parametrize(
    "reaction_family,raw_label",
    [
        # These three were the false zeros: the short offset keys `schiff`,
        # `enol` and `cys` never reached the engine's own family labels because
        # `get_barrier()` matches offsets against the NORMALISED RAW label before
        # canonicalising it.
        ("schiff_condensation", "Schiff_Base_Formation"),
        ("1,2-enolisation", "Enolisation_1_2"),
        ("cysteine_thermolysis", "Cysteine_Degradation"),
        # Same class of bug, found while fixing the three above.
        ("retro_aldol", "Retro_Aldol_Fragmentation"),
        ("thiol_addition", "Thiol_Addition_H2"),
        # These two already routed; pinned so the fix cannot regress them.
        ("amadori_rearrangement", "Amadori_Rearrangement"),
        ("strecker_degradation", "Strecker_Degradation"),
    ],
)
def test_resolved_offset_keys_actually_shift_the_barrier(reaction_family, raw_label, monkeypatch):
    """The published offset keys must move `get_barrier()` for the raw label.

    This is the assertion that would have caught the false zero: the pre-fix
    sweep published only the short key, and `get_barrier(raw_label)` returned the
    UNCHANGED barrier, so every downstream sensitivity number was structurally
    0.00 rather than measured.
    """
    from src.barrier_constants import get_barrier
    from src.family_sensitivity import DEFAULT_FAMILY_OFFSET_KEYS, resolve_offset_keys

    offset_key = DEFAULT_FAMILY_OFFSET_KEYS[reaction_family]
    delta = -3.0

    monkeypatch.delenv("BARRIER_OFFSETS", raising=False)
    baseline_barrier, _ = get_barrier(raw_label)

    keys = resolve_offset_keys(reaction_family, offset_key)
    monkeypatch.setenv("BARRIER_OFFSETS", json.dumps({key: delta for key in keys}))
    shifted_barrier, _ = get_barrier(raw_label)

    assert shifted_barrier == pytest.approx(baseline_barrier + delta), (
        f"offset keys {keys} are a no-op for engine family {raw_label!r}: "
        f"barrier stayed at {baseline_barrier}"
    )


@pytest.mark.parametrize("reaction_family", ["schiff_condensation", "1,2-enolisation", "cysteine_thermolysis"])
def test_short_offset_key_alone_is_still_a_no_op(reaction_family):
    """Pins the ROOT CAUSE so the fix cannot be reverted to the short key.

    `get_barrier()` lives in `src/barrier_constants.py`, which Wave I did not
    touch. Documenting that the short key is inert is what makes the expansion in
    `resolve_offset_keys` load-bearing rather than decorative.

    If `get_barrier()` is ever fixed to canonicalise BEFORE the offset lookup,
    this test SHOULD start failing -- that is the signal that the expansion in
    `resolve_offset_keys` has become redundant and can be simplified, not a
    reason to delete the expansion without checking.
    """
    from src.barrier_constants import get_barrier
    from src.family_sensitivity import DEFAULT_FAMILY_OFFSET_KEYS

    raw_label = {
        "schiff_condensation": "Schiff_Base_Formation",
        "1,2-enolisation": "Enolisation_1_2",
        "cysteine_thermolysis": "Cysteine_Degradation",
    }[reaction_family]
    short_key = DEFAULT_FAMILY_OFFSET_KEYS[reaction_family]

    previous = os.environ.pop("BARRIER_OFFSETS", None)
    try:
        baseline_barrier, _ = get_barrier(raw_label)
        os.environ["BARRIER_OFFSETS"] = json.dumps({short_key: -3.0})
        short_key_barrier, _ = get_barrier(raw_label)
    finally:
        if previous is None:
            os.environ.pop("BARRIER_OFFSETS", None)
        else:
            os.environ["BARRIER_OFFSETS"] = previous

    assert short_key_barrier == pytest.approx(baseline_barrier)


def test_family_sensitivity_artifact_reports_the_expanded_keys():
    """The artifact must record what was actually written to BARRIER_OFFSETS."""
    from src.family_sensitivity import build_family_sensitivity_artifact

    payload = build_family_sensitivity_artifact(benchmark_files=[])
    assert payload["summary"]["evaluated_benchmark_count"] == 0

    # With no benchmarks the sweep short-circuits, so exercise the row shape via
    # the resolver directly on the tracked artifact when it exists.
    artifact = ROOT / "results" / "validation" / "family_sensitivity.json"
    if not artifact.exists():
        # LEGITIMATE SKIP, reason sharpened 2026-08-27 (Wave J2). The precondition is not
        # "someone forgot to run a generator" -- it is that `.gitignore` excludes
        # `results/validation/*` wholesale and family_sensitivity.json is NOT among the
        # force-added exceptions, so it is absent from every fresh clone and therefore from
        # CI. This test can only ever run in a working tree where
        # `scripts/generators/generate_family_sensitivity.py` has been run locally.
        #
        # CONSEQUENCE, STATED RATHER THAN HIDDEN: on CI this assertion is not enforced, so
        # the Wave I routing fix is NOT protected by a gate. Force-adding the artifact (as
        # was done for prediction_uncertainty.json and external_validation_report.json)
        # would fix that; .gitignore is outside this file's ownership, so it is recorded as
        # a handoff in tasks/audit_remediation.md instead of changed here.
        pytest.skip(
            "results/validation/family_sensitivity.json is gitignored (results/validation/* "
            "is excluded and this file is not force-added), so it is absent on CI and in "
            "fresh clones; generate it locally to exercise this check"
        )

    rows = json.loads(artifact.read_text(encoding="utf-8"))["families"]
    by_family = {row["reaction_family"]: row for row in rows}
    for family in ("schiff_condensation", "1,2-enolisation", "cysteine_thermolysis"):
        assert len(by_family[family]["offset_keys"]) > 1, (
            f"{family} published only its short offset key; the routing fix is not in the artifact"
        )


def test_generated_family_sensitivity_shows_real_enolisation_sensitivity():
    """1,2-enolisation was a false zero and is the sweep's top family once routed."""
    artifact = ROOT / "results" / "validation" / "family_sensitivity.json"
    if not artifact.exists():
        # LEGITIMATE SKIP, reason sharpened 2026-08-27 (Wave J2). The precondition is not
        # "someone forgot to run a generator" -- it is that `.gitignore` excludes
        # `results/validation/*` wholesale and family_sensitivity.json is NOT among the
        # force-added exceptions, so it is absent from every fresh clone and therefore from
        # CI. This test can only ever run in a working tree where
        # `scripts/generators/generate_family_sensitivity.py` has been run locally.
        #
        # CONSEQUENCE, STATED RATHER THAN HIDDEN: on CI this assertion is not enforced, so
        # the Wave I routing fix is NOT protected by a gate. Force-adding the artifact (as
        # was done for prediction_uncertainty.json and external_validation_report.json)
        # would fix that; .gitignore is outside this file's ownership, so it is recorded as
        # a handoff in tasks/audit_remediation.md instead of changed here.
        pytest.skip(
            "results/validation/family_sensitivity.json is gitignored (results/validation/* "
            "is excluded and this file is not force-added), so it is absent on CI and in "
            "fresh clones; generate it locally to exercise this check"
        )

    rows = json.loads(artifact.read_text(encoding="utf-8"))["families"]
    by_family = {row["reaction_family"]: row for row in rows}

    def max_abs_mae_shift(family: str) -> float:
        return max(
            max(abs(bench["mae_shift_ppb"]) for bench in scenario["benchmarks"])
            for scenario in by_family[family]["scenarios"]
        )

    # The weighted impact score is a composite that other chemistry work can move
    # around; the load-bearing assertion is that a +/-3 kcal/mol offset on these
    # two families reaches the predictions AT ALL. Pre-fix it was identically 0.
    assert max_abs_mae_shift("1,2-enolisation") > 0.0
    assert max_abs_mae_shift("schiff_condensation") > 0.0
    assert by_family["1,2-enolisation"]["max_weighted_impact_score"] > 0.0


# ---------------------------------------------------------------------------
# FIX 16 -- LaTeX / dvipng graceful degradation
# ---------------------------------------------------------------------------


def test_configure_science_plot_style_falls_back_to_mathtext(monkeypatch):
    import matplotlib

    import src.plot_style as plot_style

    monkeypatch.delenv(plot_style.STRICT_LATEX_ENV_VAR, raising=False)
    monkeypatch.setattr(plot_style, "latex_contract_failure", lambda: "simulated: dvipng missing")
    monkeypatch.setattr(plot_style, "_FALLBACK_WARNING_EMITTED", False)

    plot_style.configure_science_plot_style()
    assert matplotlib.rcParams["text.usetex"] is False


def test_configure_science_plot_style_still_raises_under_strict_flag(monkeypatch):
    import src.plot_style as plot_style

    monkeypatch.setenv(plot_style.STRICT_LATEX_ENV_VAR, "1")
    monkeypatch.setattr(plot_style, "latex_contract_failure", lambda: "simulated: dvipng missing")

    with pytest.raises(RuntimeError, match="simulated: dvipng missing"):
        plot_style.configure_science_plot_style()


def test_configure_science_plot_style_warns_once(monkeypatch, caplog):
    import src.plot_style as plot_style

    monkeypatch.delenv(plot_style.STRICT_LATEX_ENV_VAR, raising=False)
    monkeypatch.setattr(plot_style, "latex_contract_failure", lambda: "simulated: dvipng missing")
    monkeypatch.setattr(plot_style, "_FALLBACK_WARNING_EMITTED", False)

    with caplog.at_level("WARNING", logger="src.plot_style"):
        plot_style.configure_science_plot_style()
        plot_style.configure_science_plot_style()

    warnings = [record for record in caplog.records if "mathtext" in record.getMessage()]
    assert len(warnings) == 1


def test_configure_science_plot_style_runs_on_this_machine(monkeypatch):
    """Acceptance criterion: the documented conda path must not raise here."""
    import src.plot_style as plot_style

    monkeypatch.delenv(plot_style.STRICT_LATEX_ENV_VAR, raising=False)
    plot_style.configure_science_plot_style()


# ---------------------------------------------------------------------------
# FIX 17 -- guardrail must validate the candidate constants
# ---------------------------------------------------------------------------


_GUARDRAIL_PROBE_TEST = '''
def test_candidate_constants_are_visible_to_the_guardrail():
    from src.matrix_correction import MATRIX_CORRECTIONS, ProteinType

    retention = MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE].volatile_retention_denatured
    # The baseline is ~0.63. A candidate at 0.01 must be visible here and must
    # fail. Under the pre-Wave-I subprocess the guardrail re-imported the module
    # from disk and always saw the baseline, so this assertion passed and every
    # candidate was waved through.
    assert retention > 0.5, f"guardrail saw retention={retention}"
'''


@pytest.fixture()
def guardrail_probe(tmp_path: Path) -> str:
    probe = tmp_path / "test_guardrail_probe.py"
    probe.write_text(_GUARDRAIL_PROBE_TEST, encoding="utf-8")
    return str(probe)


def test_guardrail_validates_the_candidate_not_the_baseline(guardrail_probe):
    """A candidate that must fail the guardrail has to be reported as failing."""
    from src.matrix_calibration_optimizer import (
        _run_guardrail_tests,
        candidate_constants_payload,
    )
    from src.matrix_correction import ProteinType

    bad_candidate = candidate_constants_payload(
        ProteinType.SOY_ISOLATE, [0.01, 0.01, 0.01, 0.01]
    )
    assert _run_guardrail_tests(bad_candidate, test_paths=[guardrail_probe]) is False

    good_candidate = candidate_constants_payload(
        ProteinType.SOY_ISOLATE, [0.9, 0.9, 0.9, 0.9]
    )
    assert _run_guardrail_tests(good_candidate, test_paths=[guardrail_probe]) is True


def test_guardrail_does_not_mutate_the_calling_process_constants(guardrail_probe):
    from src.matrix_calibration_optimizer import (
        _run_guardrail_tests,
        candidate_constants_payload,
    )
    from src.matrix_correction import MATRIX_CORRECTIONS, ProteinType

    before = MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE].volatile_retention_denatured
    _run_guardrail_tests(
        candidate_constants_payload(ProteinType.SOY_ISOLATE, [0.01, 0.01, 0.01, 0.01]),
        test_paths=[guardrail_probe],
    )
    assert MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE].volatile_retention_denatured == before


def test_calibrate_matrix_constants_restores_globals_on_success(monkeypatch, tmp_path):
    """The success path used to leave the optimiser's last vector in the globals."""
    import numpy as np

    import src.matrix_calibration_optimizer as optimizer
    from src.matrix_correction import (
        MATRIX_CORRECTIONS,
        VOLATILE_CLASS_RETENTION_PROFILES,
        ProteinType,
    )

    before_corr = MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE]
    before_prof = VOLATILE_CLASS_RETENTION_PROFILES[ProteinType.SOY_ISOLATE]

    seen: dict[str, object] = {}

    # Cheap synthetic objective: strictly improves away from x0, and mutates the
    # globals exactly the way the real one does.
    def fake_error(experiments, protein_type, x):
        optimizer._monkey_patch_matrix_constants(
            protein_type,
            volatile_retention_denatured=float(x[0]),
            aldehyde_factor=float(x[1]),
            sulfur_factor=float(x[2]),
            pyrazine_factor=float(x[3]),
        )
        return float(np.sum((np.asarray(x) - 0.02) ** 2)) + 0.5

    def fake_guardrail(candidate=None, test_paths=None):
        seen["candidate"] = candidate
        seen["live_retention"] = MATRIX_CORRECTIONS[
            ProteinType.SOY_ISOLATE
        ].volatile_retention_denatured
        return True

    monkeypatch.setattr(optimizer, "_compute_prediction_error", fake_error)
    monkeypatch.setattr(optimizer, "_run_guardrail_tests", fake_guardrail)
    monkeypatch.setattr(optimizer, "CALIBRATION_HISTORY_DIR", tmp_path / "calibration_history")

    result = optimizer.calibrate_matrix_constants([{"benchmark_id": "synthetic"}], "soy_iso")

    assert result is not None
    # The guardrail must have been handed the candidate, and the live globals
    # must have been the candidate's at that moment.
    assert seen["candidate"]["volatile_retention_denatured"] == pytest.approx(
        result["parameters"]["volatile_retention_denatured"]["new"]
    )
    assert seen["live_retention"] == pytest.approx(
        result["parameters"]["volatile_retention_denatured"]["new"]
    )
    # ... and the globals must be back to baseline afterwards.
    assert MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE] is before_corr
    assert VOLATILE_CLASS_RETENTION_PROFILES[ProteinType.SOY_ISOLATE] is before_prof
