"""Wave G3 — front-door UX, CLI honesty, schema unification and determinism.

Covers the fixes landed on 2026-08-27 for the audit's "front-door smoke +
determinism" findings and module-hunt TIER 2 optimizer items. Each test names
the defect it pins so a regression is self-explaining.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.bayesian_optimizer import (  # noqa: E402
    CATEGORICAL_PARAM_KEYS,
    CONCENTRATION_PARAM_KEYS,
    amino_acid_concentration_param,
    build_molar_ratios,
    formulation_from_params,
    interventions_from_params,
    pre_processing_steps_from_params,
)
from src.formulation import CANONICAL_PIPELINE_CONDITION_KEYS, Formulation  # noqa: E402
from src.pipeline import MaillardPipeline, _extract_unresolved_names, _short_resolver_reason  # noqa: E402
from src.pre_processor import (  # noqa: E402
    KNOWN_INTERVENTIONS,
    PreProcessor,
    intervention_is_active,
    resolve_intervention,
)
from src.projection_utils import FINGERPRINT_EXCLUDED_INPUT_KEYS, fingerprint_inputs  # noqa: E402
from src.reporting import _md_cell, _render_glossary_markdown  # noqa: E402


def _load_script(name: str):
    """Import a top-level script from scripts/ as a module."""
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_script_{name}", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# ─────────────────────────────────────────────────────────────────────────────
# Fix 2 — optimizer trial-parameter mapping
# ─────────────────────────────────────────────────────────────────────────────

BEST_TRIAL_PARAMS = {
    "sugar_conc": 0.25,
    "aa_conc_sulfur": 0.75,
    "aa_conc_branched": 0.11,
    "aa_conc_basic": 0.42,
    "aa_conc_other": 0.06,
    "ph": 7.42,
    "temp": 118.2,
    "aw": 0.37,
    "time_minutes": 70.9,
    "intervention_agent": "rosemary_extract",
    "intervention_dose": 0.5,
    "pre_processing": "none",
}


def test_amino_acid_concentration_param_routes_by_class():
    assert amino_acid_concentration_param("Cysteine") == "aa_conc_sulfur"
    assert amino_acid_concentration_param("methionine") == "aa_conc_sulfur"
    assert amino_acid_concentration_param("leucine") == "aa_conc_branched"
    assert amino_acid_concentration_param("Lysine") == "aa_conc_basic"
    assert amino_acid_concentration_param("phenylalanine") == "aa_conc_other"


def test_build_molar_ratios_maps_class_knobs_onto_real_precursor_names():
    ratios = build_molar_ratios(
        BEST_TRIAL_PARAMS,
        ["ribose", "glucose"],
        ["cysteine", "leucine", "phenylalanine"],
        ["hexanal"],
    )
    assert ratios["ribose"] == pytest.approx(0.25)
    assert ratios["glucose"] == pytest.approx(0.25)
    assert ratios["cysteine"] == pytest.approx(0.75)
    assert ratios["leucine"] == pytest.approx(0.11)
    assert ratios["phenylalanine"] == pytest.approx(0.06)
    assert ratios["hexanal"] == pytest.approx(0.1)


def test_optimizer_formulation_does_not_leak_ph_or_categoricals_into_molar_ratios():
    """The old CLI passed `trial.params` straight in as `molar_ratios`.

    Two consequences, both pinned here: 'ph' substring-matched
    'L-Phenylalanine' (leaking pH 7.42 in as that amino acid's concentration),
    and float('rosemary_extract') crashed the run outright.
    """
    formulation = formulation_from_params(
        BEST_TRIAL_PARAMS,
        name="Best_Trial",
        fixed_sugars=["ribose"],
        fixed_amino_acids=["phenylalanine"],
        fixed_lipids=[],
        protein_type="pea_iso",
    )
    ratios = formulation["molar_ratios"]

    # No condition or categorical key survives into the concentration map.
    for leaked in ("ph", "temp", "aw", "time_minutes") + CATEGORICAL_PARAM_KEYS:
        assert leaked not in ratios
    # Every remaining value is a real number (the float() crash class).
    for value in ratios.values():
        assert isinstance(float(value), float)
    # And phenylalanine carries its own class concentration, not the pH.
    assert ratios["phenylalanine"] == pytest.approx(BEST_TRIAL_PARAMS["aa_conc_other"])
    assert ratios["phenylalanine"] != pytest.approx(BEST_TRIAL_PARAMS["ph"])

    # Conditions land where conditions belong.
    assert formulation["ph"] == pytest.approx(7.42)
    assert formulation["temp"] == pytest.approx(118.2)
    assert formulation["aw"] == pytest.approx(0.37)
    assert formulation["time_minutes"] == pytest.approx(70.9)


def test_optimizer_formulation_is_not_a_flat_ratio_map():
    """Unmapped names used to fall back to a flat 1:1:1 in the pipeline."""
    formulation = formulation_from_params(
        BEST_TRIAL_PARAMS,
        name="Best_Trial",
        fixed_sugars=["ribose"],
        fixed_amino_acids=["cysteine", "leucine"],
        fixed_lipids=[],
    )
    values = set(formulation["molar_ratios"].values())
    assert len(values) > 1
    assert values != {1.0}


def test_interventions_and_preprocessing_round_trip_from_params():
    assert interventions_from_params(BEST_TRIAL_PARAMS) == [
        {"name": "rosemary_extract", "dose": 0.5}
    ]
    assert interventions_from_params({"intervention_agent": "none"}) == []
    assert pre_processing_steps_from_params({"pre_processing": "both"}) == [
        "yeast_fermentation",
        "protease_hydrolysis",
    ]
    assert pre_processing_steps_from_params({"pre_processing": "none"}) == []
    assert pre_processing_steps_from_params({"pre_processing": "yeast_fermentation"}) == [
        "yeast_fermentation"
    ]


def test_optimizer_cli_prints_real_units_and_survives_categorical_params(capsys):
    optimize = _load_script("optimize_formulation")
    optimize.print_optimized_parameters(
        dict(BEST_TRIAL_PARAMS),
        sugars=["ribose"],
        aas=["cysteine", "phenylalanine"],
        lipids=[],
    )
    out = capsys.readouterr().out
    assert "mM" in out
    # The unit label was wrong (values are mM, printed as "M").
    assert " M\n" not in out
    # Categorical params render as themselves rather than crashing.
    assert "rosemary_extract" in out
    assert "pre_processing" in out


def test_optimizer_cli_flags_a_degenerate_objective():
    optimize = _load_script("optimize_formulation")
    flat = SimpleNamespace(
        trials=[SimpleNamespace(value=-1.5) for _ in range(5)],
    )
    spread = optimize.describe_objective_spread(flat)
    assert spread["degenerate"] is True
    assert spread["n"] == 5

    moving = SimpleNamespace(
        trials=[SimpleNamespace(value=v) for v in (-3.0, -1.5, -2.0)],
    )
    spread = optimize.describe_objective_spread(moving)
    assert spread["degenerate"] is False
    assert spread["spread"] == pytest.approx(1.5)


def test_wheat_gluten_is_not_offered_as_a_protein_type(monkeypatch, capsys):
    """`--protein-type wheat_gluten` died in ProteinType() on every run.

    ProteinType has no wheat_gluten member and the matrix layer has no
    wheat-gluten profile, so the CLIs must reject it up front rather than
    crashing deep in the matrix correction.
    """
    from src.matrix_correction import ProteinType

    assert "wheat_gluten" not in {member.value for member in ProteinType}

    run_pipeline = _load_script("run_pipeline")
    monkeypatch.setattr(
        sys, "argv", ["run_pipeline.py", "--sugars", "ribose", "--protein-type", "wheat_gluten"]
    )
    with pytest.raises(SystemExit) as excinfo:
        run_pipeline.parse_args()
    assert excinfo.value.code == 2
    assert "wheat_gluten" in capsys.readouterr().err

    optimize_source = (ROOT / "scripts" / "optimize_formulation.py").read_text()
    choices_line = next(
        line for line in optimize_source.splitlines() if "--protein-type" in line and "choices" in line
    )
    assert "wheat_gluten" not in choices_line.split("help=")[0]


# ─────────────────────────────────────────────────────────────────────────────
# Fix 4 — intervention token matching
# ─────────────────────────────────────────────────────────────────────────────


def test_known_interventions_are_the_ones_apply_can_simulate():
    assert set(KNOWN_INTERVENTIONS) == {"yeast_fermentation", "protease_hydrolysis"}


@pytest.mark.parametrize(
    "entry",
    [
        "no_yeast_fermentation",
        "skip_yeast_fermentation",
        {"yeast_fermentation": False},
        {"yeast_fermentation": None},
        {"yeast_fermentation": 0},
        {"yeast_fermentation": "false"},
        {"name": "protease_hydrolysis"},
        "protease_hydrolysis",
    ],
)
def test_yeast_fermentation_is_not_activated_by_these(entry):
    active, _ = resolve_intervention(entry, "yeast_fermentation")
    assert active is False


@pytest.mark.parametrize(
    "entry",
    [
        "yeast_fermentation",
        "  Yeast_Fermentation  ",
        {"yeast_fermentation": True},
        {"yeast_fermentation": {"time_hours": 5}},
        {"name": "yeast_fermentation", "dose": 0.4},
    ],
)
def test_yeast_fermentation_is_activated_by_these(entry):
    active, _ = resolve_intervention(entry, "yeast_fermentation")
    assert active is True


def test_disabling_yeast_fermentation_leaves_the_pool_untouched():
    """Regression: `{"yeast_fermentation": False}` used to turn it ON."""
    processor = PreProcessor()
    ratios = {"hexanal": 1.0, "nonanal": 0.5}

    off = processor.apply(dict(ratios), [{"yeast_fermentation": False}])
    assert off == ratios
    assert "hexanol" not in off

    negated = processor.apply(dict(ratios), ["no_yeast_fermentation"])
    assert negated == ratios

    on = processor.apply(dict(ratios), ["yeast_fermentation"])
    assert on["hexanal"] == pytest.approx(0.2)
    assert on["hexanol"] == pytest.approx(0.8)


def test_yeast_fermentation_time_hours_still_sets_efficiency():
    processor = PreProcessor()
    out = processor.apply({"hexanal": 1.0}, [{"yeast_fermentation": {"time_hours": 5}}])
    # eff = 1 - exp(-0.4 * 5) = 0.8647
    assert out["hexanol"] == pytest.approx(1.0 - 2.718281828459045 ** (-2.0), rel=1e-6)


def test_protease_hydrolysis_matches_dict_and_string_forms():
    processor = PreProcessor()
    for interventions in (["protease_hydrolysis"], [{"protease_hydrolysis": True}], [{"name": "protease_hydrolysis"}]):
        out = processor.apply({"cysteine": 1.0}, interventions)
        assert out["cysteine"] == pytest.approx(2.0)
    assert intervention_is_active(["no_protease_hydrolysis"], "protease_hydrolysis") is False


# ─────────────────────────────────────────────────────────────────────────────
# Fix 5 — Formulation schema unification
# ─────────────────────────────────────────────────────────────────────────────


def test_formulation_to_dict_emits_the_keys_the_pipeline_reads():
    """to_dict emitted temperature/water_activity; the pipeline reads temp/aw."""
    formulation = Formulation(name="f", temperature=137.0, water_activity=0.61, ph=5.2)
    payload = formulation.to_dict()
    for key in CANONICAL_PIPELINE_CONDITION_KEYS:
        assert key in payload, f"pipeline reads {key!r} but to_dict does not emit it"
    assert payload["temp"] == 137.0
    assert payload["aw"] == 0.61
    assert payload["ph"] == 5.2
    # Legacy dataclass names kept so nothing that read them breaks.
    assert payload["temperature"] == 137.0
    assert payload["water_activity"] == 0.61


def test_formulation_round_trips_through_to_dict():
    original = Formulation(
        name="round_trip",
        sugars=["ribose"],
        amino_acids=["cysteine"],
        molar_ratios={"ribose": 1.5},
        ph=5.5,
        temperature=105.0,
        water_activity=0.8,
        time_minutes=45.0,
        protein_type="pea_iso",
    )
    restored = Formulation.from_dict(original.to_dict())
    assert restored == original


def test_formulation_accepts_the_natural_pH_spelling():
    formulation = Formulation.from_dict({"name": "f", "pH": "6.4"})
    assert formulation.ph == 6.4
    assert formulation.get("pH") == 6.4
    assert formulation.get("temp") == formulation.temperature


# ─────────────────────────────────────────────────────────────────────────────
# Fix 3 — precursor-skip surfacing
# ─────────────────────────────────────────────────────────────────────────────


def test_resolver_error_is_reduced_to_the_offending_name():
    message = "Unknown precursor 'unobtanium'. Available: alanine, arginine, ribose"
    assert _extract_unresolved_names(message, ["ribose", "unobtanium"]) == ["unobtanium"]
    assert _short_resolver_reason(message) == "Unknown precursor 'unobtanium'."


def test_unresolvable_precursor_is_reported_not_silently_dropped():
    designer = MaillardPipeline("meaty", "beany")
    bogus = {
        "name": "Bogus_Formulation",
        "sugars": ["ribose"],
        "amino_acids": ["unobtanium"],
    }
    from src.conditions import ReactionConditions

    conditions = ReactionConditions(pH=6.0, temperature_celsius=120.0, water_activity=0.8)
    with pytest.raises(ValueError) as excinfo:
        designer.evaluate_single(bogus, conditions)

    # The error names the formulation and the offending precursor.
    assert "Bogus_Formulation" in str(excinfo.value)
    assert "unobtanium" in str(excinfo.value)

    skipped = designer.last_skipped_formulations
    assert len(skipped) == 1
    assert skipped[0]["name"] == "Bogus_Formulation"
    assert skipped[0]["unresolved_precursors"] == "unobtanium"
    # The giant catalogue of available names is not dumped into the payload.
    assert "Available:" not in skipped[0]["reason"]


# ─────────────────────────────────────────────────────────────────────────────
# Fix 6 — grid name handling
# ─────────────────────────────────────────────────────────────────────────────


def test_grid_name_resolution_is_tolerant_but_never_silent():
    designer = MaillardPipeline("meaty", "beany")
    exact = designer.resolve_grid_formulation("Cysteine Enrichment (Basic)")
    assert exact["name"] == "Cysteine Enrichment (Basic)"
    # Unique substring match (what the quickstart's shorthand needs).
    assert designer.resolve_grid_formulation("Cysteine Enrichment")["name"] == (
        "Cysteine Enrichment (Basic)"
    )
    # Case-insensitive.
    assert designer.resolve_grid_formulation("premium meaty mix")["name"] == "Premium Meaty Mix"

    with pytest.raises(ValueError) as excinfo:
        designer.resolve_grid_formulation("Baseline")
    # The documented-but-nonexistent name now yields the valid ones.
    assert "Available:" in str(excinfo.value)
    assert "Soy/Pea Base (Untreated)" in str(excinfo.value)


def test_quickstart_campaign_spec_names_exist_in_the_grid():
    import yaml

    spec_path = ROOT / "docs" / "examples" / "shareable_meaty_screen.yml"
    assert spec_path.exists(), "QUICKSTART references this campaign spec"
    spec = yaml.safe_load(spec_path.read_text())
    designer = MaillardPipeline(
        str(spec["campaign"]["target_tag"]), str(spec["campaign"]["minimize_tag"])
    )
    for entry in spec["formulations"]:
        assert designer.resolve_grid_formulation(entry["name"]) is not None


def test_quickstart_lane_uses_real_grid_names():
    shell = (ROOT / "scripts" / "docker_maillard.sh").read_text()
    assert "'Baseline,Cysteine Enrichment'" not in shell
    assert "Soy/Pea Base (Untreated),Cysteine Enrichment (Basic)" in shell


# ─────────────────────────────────────────────────────────────────────────────
# Fix 7 — markdown table integrity and glossary truthfulness
# ─────────────────────────────────────────────────────────────────────────────


def test_md_cell_keeps_a_pipe_joined_label_inside_one_cell():
    value = "static_class_profile | class_level | standard_matrix_support"
    rendered = _md_cell(value)
    row = f"| a | {rendered} | c |"
    # A markdown row splits on unescaped pipes only.
    assert row.count("|") - row.count("\\|") == 4


def test_glossary_describes_the_vocabularies_that_are_actually_emitted():
    glossary = _render_glossary_markdown()
    # The real per-compound tier values.
    for value in ("`high`", "`medium`", "`low`", "`exploratory`"):
        assert value in glossary
    # The real evidence-strength values.
    for value in ("literature_anchored", "class_anchored", "directional_transferred", "heuristic"):
        assert value in glossary
    # The separate literature-prior scale.
    assert "confidence_tier" in glossary
    assert "medium_high" in glossary
    # The fictional compound-level vocabulary is gone as a *tier* claim.
    assert "**Evidence tiers**" not in glossary


def test_readme_sample_output_matches_the_real_report_columns():
    readme = (ROOT / "README.md").read_text()
    assert "| Compound | Predicted | Tier | Score | Mode | Reachability |" in readme
    # The fictional per-compound confidence vocabulary is gone from the sample.
    sample = readme.split("## What the output looks like")[1].split("## How well calibrated")[0]
    assert "`bounded_calibration`" not in sample
    assert "`surrogate_family`" not in sample
    assert "exploratory" in sample


# ─────────────────────────────────────────────────────────────────────────────
# Fix 8 — determinism
# ─────────────────────────────────────────────────────────────────────────────


def test_input_fingerprint_ignores_output_dir_and_argv():
    a = {"sugars": "ribose", "ph": 5.5, "output_dir": "results/run_a", "report": True}
    b = {"sugars": "ribose", "ph": 5.5, "output_dir": "results/run_b", "report": True}
    assert fingerprint_inputs(a) == fingerprint_inputs(b)
    assert json.dumps(fingerprint_inputs(a), sort_keys=True) == json.dumps(
        fingerprint_inputs(b), sort_keys=True
    )
    assert "output_dir" in FINGERPRINT_EXCLUDED_INPUT_KEYS
    # A genuine scientific difference still moves the fingerprint.
    c = dict(a, ph=6.5)
    assert fingerprint_inputs(a) != fingerprint_inputs(c)


def test_build_artifact_provenance_fingerprint_is_output_dir_independent(tmp_path):
    from src.projection_utils import build_artifact_provenance

    dir_a = tmp_path / "a"
    dir_b = tmp_path / "b"
    dir_a.mkdir()
    dir_b.mkdir()
    inputs_a = {"ph": 5.5, "temp": 105.0, "output_dir": str(dir_a)}
    inputs_b = {"ph": 5.5, "temp": 105.0, "output_dir": str(dir_b)}
    prov_a = build_artifact_provenance("single_run_report", dir_a, inputs_a)
    prov_b = build_artifact_provenance("single_run_report", dir_b, inputs_b)
    assert prov_a["input_fingerprint_sha256"] == prov_b["input_fingerprint_sha256"]


def test_attribution_total_is_summed_in_a_stable_order():
    """The sum used to follow dict order, which follows PYTHONHASHSEED."""
    from src.usability_reports import _fallback_sensitivity_summary

    values = {"ribose": 0.1, "cysteine": 0.2, "leucine": 0.30000000000000004, "lysine": 0.7}
    forward = SimpleNamespace(
        precursor_contributions=dict(values),
        flagged_toxics=[],
        off_flavour_risk=0.0,
        target_score=0.0,
        safety_score=0.0,
        texture_risk=0.0,
    )
    reversed_order = SimpleNamespace(
        precursor_contributions={k: values[k] for k in reversed(list(values))},
        flagged_toxics=[],
        off_flavour_risk=0.0,
        target_score=0.0,
        safety_score=0.0,
        texture_risk=0.0,
    )
    assert json.dumps(_fallback_sensitivity_summary(forward), sort_keys=True) == json.dumps(
        _fallback_sensitivity_summary(reversed_order), sort_keys=True
    )


@pytest.mark.slow
def test_two_in_process_runs_are_bit_identical():
    from src.conditions import ReactionConditions

    formulation = {
        "name": "Determinism_Probe",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine", "lysine"],
        "lipids": ["hexanal"],
        "molar_ratios": {"ribose": 0.5, "cysteine": 0.2, "lysine": 0.1},
        "ph": 5.5,
        "temp": 105.0,
        "aw": 0.8,
        "time_minutes": 45.0,
        "protein_type": "pea_iso",
    }
    conditions = ReactionConditions(
        pH=5.5, temperature_celsius=105.0, water_activity=0.8, protein_type="pea_iso"
    )

    def snapshot():
        designer = MaillardPipeline("meaty", "beany")
        res = designer.evaluate_single(formulation, conditions)
        return json.dumps(
            {
                "predicted_ppb": {k: repr(v) for k, v in sorted(res.predicted_ppb.items())},
                "attribution": {
                    k: repr(v) for k, v in sorted(res.precursor_contributions.items())
                },
                "target_score": repr(res.target_score),
                "off_flavour_risk": repr(res.off_flavour_risk),
                "safety_score": repr(res.safety_score),
                "lysine_budget": repr(res.lysine_budget),
                "trapping_efficiency": repr(res.trapping_efficiency),
            },
            sort_keys=True,
        )

    assert snapshot() == snapshot()


# ─────────────────────────────────────────────────────────────────────────────
# Fix 1 — run_pipeline mode honesty
# ─────────────────────────────────────────────────────────────────────────────


def _args(**overrides) -> argparse.Namespace:
    base = dict(
        sugars="ribose",
        amino_acids="cysteine",
        additives="",
        lipids="",
        ratios="ribose:0.5,cysteine:0.2",
        ph=5.5,
        temp=105.0,
        aw=0.8,
        time_minutes=45.0,
        catalyst=None,
        sme_kj_per_kg=0.0,
        moisture_regime=None,
        sterilization_temp=None,
        sterilization_time_minutes=0.0,
        barrel_zones="",
        barrel_zone_time_fractions="",
        target="meaty",
        minimize="beany",
        xtb=False,
        protein_type="pea_iso",
        denaturation_state=0.5,
        report=False,
        output_dir=None,
        dry_run=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def test_forward_and_inverse_modes_build_the_same_formulation():
    run_pipeline = _load_script("run_pipeline")
    args = _args()
    forward = run_pipeline.build_formulation_from_args(args, "Single_Run")
    custom = run_pipeline.build_formulation_from_args(args, run_pipeline.CUSTOM_CANDIDATE_NAME)
    assert forward["molar_ratios"] == {"ribose": 0.5, "cysteine": 0.2}
    assert forward["time_minutes"] == 45.0
    assert {k: v for k, v in forward.items() if k != "name"} == {
        k: v for k, v in custom.items() if k != "name"
    }


def test_inverse_design_detects_a_user_supplied_formulation():
    run_pipeline = _load_script("run_pipeline")
    assert run_pipeline._user_supplied_precursors(_args()) is True
    assert (
        run_pipeline._user_supplied_precursors(_args(sugars="", amino_acids="", additives="", lipids=""))
        is False
    )


def test_inverse_design_report_inputs_describe_the_evaluated_candidate():
    """§1 used to echo the user's args on a report about a grid entry."""
    run_pipeline = _load_script("run_pipeline")
    winner = {
        "name": "Cysteine Enrichment (Basic)",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine", "glycine", "lysine"],
        "molar_ratios": {"ribose": 1.0, "cysteine": 2.0},
        "lipids": ["hexanal", "nonanal"],
    }
    payload = run_pipeline._report_inputs_for_inverse_design(
        _args(),
        winning_formulation=winner,
        candidate_count=16,
        custom_candidate_included=True,
        ignored_arguments=[],
    )
    assert payload["mode"] == "inverse_design"
    assert payload["selected_formulation"] == "Cysteine Enrichment (Basic)"
    assert payload["selected_sugars"] == "ribose"
    assert payload["selected_molar_ratios"] == {"ribose": 1.0, "cysteine": 2.0}
    assert payload["candidates_screened"] == 16
    assert payload["custom_candidate_included"] is True
    # The user's own --ratios are not presented as the evaluated ones.
    assert payload["selected_molar_ratios"] != {"ribose": 0.5, "cysteine": 0.2}


def test_inverse_design_names_arguments_it_cannot_apply():
    run_pipeline = _load_script("run_pipeline")
    payload = run_pipeline._report_inputs_for_inverse_design(
        _args(sugars="", amino_acids="", additives="", lipids=""),
        winning_formulation={"name": "Premium Meaty Mix"},
        candidate_count=15,
        custom_candidate_included=False,
        ignored_arguments=["--ratios"],
    )
    assert "--ratios" in payload["ignored_cli_arguments"]
    assert "forward mode" in payload["ignored_cli_arguments_note"]


def test_forward_mode_report_inputs_carry_no_output_plumbing():
    run_pipeline = _load_script("run_pipeline")
    args = _args(target=None, output_dir="results/somewhere", report=True)
    formulation = run_pipeline.build_formulation_from_args(args, "Single_Run")
    payload = run_pipeline._report_inputs_for_forward_mode(args, formulation)
    assert payload["mode"] == "forward_single_formulation"
    assert "output_dir" not in payload
    assert "report" not in payload
    assert payload["sugars"] == "ribose"
    assert payload["molar_ratios"] == {"ribose": 0.5, "cysteine": 0.2}


# ─────────────────────────────────────────────────────────────────────────────
# Fix 6 — ingest CLI defaults
# ─────────────────────────────────────────────────────────────────────────────


def test_ingest_preview_defaults_away_from_the_validation_directory():
    from src import data_ingest

    parser = data_ingest._build_parser()
    args = parser.parse_args(["--file", "x.csv"])
    assert args.output_dir is None
    assert args.confirm is False
    assert data_ingest.DEFAULT_PREVIEW_OUTPUT_DIR != data_ingest.DEFAULT_OUTPUT_DIR
    assert data_ingest.DEFAULT_PREVIEW_OUTPUT_DIR.name == "ingest_previews"


def test_ingest_confirm_help_does_not_claim_artifact_regeneration():
    from src import data_ingest

    parser = data_ingest._build_parser()
    confirm = next(a for a in parser._actions if a.dest == "confirm")
    assert "regenerate" in confirm.help
    assert "does not rebuild" in confirm.help


def test_readme_no_longer_claims_ingest_regenerates_artifacts():
    readme = (ROOT / "README.md").read_text()
    assert "updates the benchmark database and regenerates validation artifacts" not in readme
    assert "--water-activity" in readme
    assert "--time-min" in readme


def test_validated_envelope_generator_exposes_an_output_dir():
    source = (ROOT / "scripts" / "generators" / "generate_validated_envelope_report.py").read_text()
    assert '"--output-dir"' in source
    # The LaTeX contract must not run at import time (it broke --help).
    body = source.split("def main()")[0]
    assert "configure_science_plot_style()" not in body
