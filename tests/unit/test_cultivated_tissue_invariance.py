"""
The cultivated-tissue invariance test (results/validation/cultivated_tissue_invariance_prereg.md).

What these tests promise, section by section of the pre-registration:

* section 4 -- the module adds nothing to the engine: importing it leaves a reference prediction
  byte-identical, and it writes only under results/cultivated_tissue_invariance/;
* section 6 A4 -- the naive ranking is fractional deficit with the declared tie-break and drop rule;
* section 3 -- the verdict is read off mechanically, and each clause fires only on its threshold;
* section 6 A7 -- every range in either box is labelled and ordered; the frozen STUB_BOX is
  byte-for-byte the box the declared run's artifact echoes; and every non-stub range in the
  current BOX names an extraction dossier that exists on disk (section 7: a stub flips only
  with a dossier).

The declared run itself is not a unit test (minutes, not seconds); the artifact under
results/cultivated_tissue_invariance/ is its record.
"""
from __future__ import annotations

import importlib
import json

import pytest

from src import data_paths
from src.comparative_cli import spec_to_core, validate_spec
from src.kinetic_core import engine


_REFERENCE_SPEC = {
    "name": "reference",
    "precursors": {"ribose": 2.0, "cysteine": 0.2, "glucose": 5.0, "thiamine": 0.003},
    "temp_C": 100.0, "time_min": 20.0, "ph": 6.0, "aw": 0.98, "matrix": "water", "protein_type": "free",
    "buffer": {"kind": "phosphate", "phosphate_mol_l": 0.03, "source": "test"},
}
_REFERENCE_TARGETS = ["2-methyl-3-furanthiol", "2-furfurylthiol", "furfural"]


def _reference_prediction():
    run = engine.predict(spec_to_core(validate_spec(_REFERENCE_SPEC, label="reference")), _REFERENCE_TARGETS)
    return json.dumps([(r["compound"], r["predicted_ug_per_l"]) for r in run.interval_rows()], sort_keys=True)


def test_importing_the_module_does_not_change_the_engine():
    before = _reference_prediction()
    module = importlib.import_module("src.cultivated_tissue_invariance")
    importlib.reload(module)
    after = _reference_prediction()
    assert before == after


def _check_box_shape(box):
    from src.cultivated_tissue_invariance import CANDIDATES, LABELS, UNRANKABLE

    for tissue, rows in box.items():
        for name, r in rows.items():
            assert r.label in LABELS, (tissue, name)
            assert 0 < r.lo_mM < r.hi_mM, (tissue, name)
            assert r.note, (tissue, name)
    assert set(box["beef"]) == set(box["cultivated_muscle"])
    assert set(CANDIDATES) | set(UNRANKABLE) == set(box["beef"])


def test_the_frozen_stub_box_is_all_stub_and_matches_the_declared_runs_echo():
    import yaml

    from src.cultivated_tissue_invariance import BOX_ECHO, STUB_BOX, box_echo

    _check_box_shape(STUB_BOX)
    for tissue, rows in STUB_BOX.items():
        for name, r in rows.items():
            assert r.label == "stub" and r.dossiers == (), (tissue, name)
    if not BOX_ECHO.exists():
        pytest.skip("the declared run's echo is not on this checkout")
    echoed = yaml.safe_load(BOX_ECHO.read_text(encoding="utf-8"))
    frozen = box_echo(STUB_BOX)
    for tissue, rows in frozen.items():
        for name, r in rows.items():
            e = echoed[tissue][name]
            assert (e["lo_mM"], e["hi_mM"], e["label"]) == (r["lo_mM"], r["hi_mM"], r["label"]), (tissue, name)


def test_every_non_stub_range_in_the_current_box_names_a_dossier_on_disk():
    from src.cultivated_tissue_invariance import BOX, CANDIDATES

    _check_box_shape(BOX)
    dossiers = data_paths.EXTRACTION_DOSSIERS_DIR
    seen_non_stub = 0
    for tissue, rows in BOX.items():
        for name, r in rows.items():
            if r.label == "stub":
                assert r.dossiers == (), (tissue, name)
                assert "NO MEASUREMENT FOUND" in r.note or tissue == "beef", (tissue, name)
                continue
            seen_non_stub += 1
            assert r.dossiers, (tissue, name)
            for stem in r.dossiers:
                assert (dossiers / f"{stem}_extraction.md").exists(), (tissue, name, stem)
    assert seen_non_stub > 0
    # The read of 2026-09-14 found nothing for the four rankable precursors in cultured muscle.
    # If a later read does, this assertion is edited to name what it found.
    for name in CANDIDATES:
        assert BOX["cultivated_muscle"][name].label == "stub", name
    # ...and every beef-side range is measured or cited.
    for name in BOX["beef"]:
        assert BOX["beef"][name].label != "stub", name


def test_gap_map_lists_exactly_the_stubs_with_a_closing_measurement_for_each_rankable_one():
    from src.cultivated_tissue_invariance import BOX, CANDIDATES, gap_map

    gaps = gap_map(BOX)
    stubs = {(t, p) for t, rows in BOX.items() for p, r in rows.items() if r.label == "stub"}
    assert {(g["tissue"], g["precursor"]) for g in gaps} == stubs
    for g in gaps:
        if g["engine_rankable"]:
            assert g["closes_with"], g
    assert all(g["engine_rankable"] == (g["precursor"] in CANDIDATES) for g in gaps)


def test_naive_ranking_is_fractional_deficit_with_the_declared_tiebreak_and_drop_rule():
    from src.cultivated_tissue_invariance import naive_ranking

    beef = {"a": 10.0, "b": 1.0, "c": 5.0, "d": 2.0}
    cult = {"a": 9.0, "b": 0.1, "c": 5.0, "d": 4.0}
    # b: 90 % down; a: 10 % down; c: nothing to restore; d: above beef -> dropped
    assert naive_ranking(beef, cult, ["a", "b", "c", "d"]) == ["b", "a"]
    # tie on fractional deficit -> absolute deficit decides
    beef = {"a": 10.0, "b": 1.0}
    cult = {"a": 5.0, "b": 0.5}
    assert naive_ranking(beef, cult, ["a", "b"]) == ["a", "b"]


def test_kendall_tau_on_small_orderings():
    from src.cultivated_tissue_invariance import kendall_tau

    assert kendall_tau(["a", "b", "c", "d"], ["a", "b", "c", "d"]) == 1.0
    assert kendall_tau(["a", "b", "c", "d"], ["d", "c", "b", "a"]) == -1.0
    assert kendall_tau(["a", "b", "c", "d"], ["b", "a", "c", "d"]) == pytest.approx(4 / 6)
    assert kendall_tau(["a"], ["a"]) is None
    assert kendall_tau(["a", "b"], ["a", "c"]) is None


def _summary(**over):
    base = {
        "refusal_fraction": 0.0,
        "disagreement_fraction": 0.0,
        "top_agreement_fraction": 1.0,
        "mean_kendall_tau": 1.0,
        "dominant_pair": None,
        "dominant_pair_envelope_survival": None,
    }
    base.update(over)
    return base


def test_verdict_reads_section_3_mechanically():
    from src.cultivated_tissue_invariance import THRESHOLDS, verdict

    assert verdict(_summary())["verdict"] == "T2"
    assert verdict(_summary(refusal_fraction=0.5))["verdict"] == "T3"
    t1 = _summary(
        disagreement_fraction=0.25, top_agreement_fraction=0.75, mean_kendall_tau=0.6,
        dominant_pair={"pair": "ribose over thiamine", "count": 10, "share_of_disagreements": 0.6},
        dominant_pair_envelope_survival=0.85,
    )
    assert verdict(t1)["verdict"] == "T1"
    # each T1 clause alone is not enough
    assert verdict({**t1, "disagreement_fraction": 0.19})["verdict"] == "indeterminate"
    assert verdict({**t1, "dominant_pair": {**t1["dominant_pair"], "share_of_disagreements": 0.49}})["verdict"] == "indeterminate"
    assert verdict({**t1, "dominant_pair_envelope_survival": 0.79})["verdict"] == "indeterminate"
    # T2 needs both clauses
    assert verdict(_summary(mean_kendall_tau=0.7))["verdict"] == "indeterminate"
    assert verdict(_summary(top_agreement_fraction=0.89, disagreement_fraction=0.11))["verdict"] == "indeterminate"
    # T3 takes precedence over everything
    assert verdict({**t1, "refusal_fraction": 0.5})["verdict"] == "T3"
    assert THRESHOLDS["T1_min_disagreement_fraction"] == 0.20
    assert THRESHOLDS["T2_min_mean_tau"] == 0.75


def test_a_spec_from_the_box_validates_and_the_sulfur_arm_answers():
    from src.cultivated_tissue_invariance import BOX, run_arm

    mid = {p: (r.lo_mM * r.hi_mM) ** 0.5 for p, r in BOX["cultivated_muscle"].items()}
    arm = run_arm(mid, 100.0, 20.0, "midpoint")
    assert arm.answered, arm.refusals
    assert arm.metric is not None and arm.metric > 0
    assert set(arm.oav) == {"2-methyl-3-furanthiol", "2-furfurylthiol", "furfural"}


def test_write_targets_only_the_results_directory(tmp_path, monkeypatch):
    import src.cultivated_tissue_invariance as m

    monkeypatch.setattr(m, "OUTPUT_DIR", tmp_path / "out")
    monkeypatch.setattr(m, "OUTPUT_JSON", tmp_path / "out" / "x.json")
    monkeypatch.setattr(m, "BOX_ECHO", tmp_path / "out" / "box.yml")
    payload = {
        "provenance": {"generated_on": "today", "generated_by": "test"},
        "pre_registration": "prereg",
        "design": {"candidates": list(m.CANDIDATES), "targets": list(m.TARGETS), "metric_targets": list(m.METRIC_TARGETS),
                   "fixed_conditions": dict(m.FIXED_CONDITIONS), "n_draws": 1, "seed": 0, "n_envelope": 1,
                   "composition_box": m.box_echo(m.BOX),
                   "structural_refusals": [], "unrankable_declared_candidates": dict(m.UNRANKABLE)},
        "programmes": [],
        "overall": {"verdicts": {}},
    }
    json_path, md_path = m.write(payload)
    assert json_path.parent == tmp_path / "out" and md_path.parent == tmp_path / "out"
    assert (tmp_path / "out" / "box.yml").exists()
    # a stem never touches the declared run's files
    j2, m2 = m.write(payload, stem="probe")
    assert j2.name == "probe.json" and m2.name == "probe.md" and (tmp_path / "out" / "probe_box.yml").exists()
    assert json_path.exists() and json_path.name == "x.json"
    assert not any(p.is_relative_to(data_paths.DATA_ROOT) for p in (json_path, md_path))
