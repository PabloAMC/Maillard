"""
tests/test_pipeline.py — Test suite for Phase 7.1 end-to-end wiring.

Verifies:
1. `precursor_resolver` correctly maps string names to SMILES.
2. `Recommender.predict_from_steps` correctly finds the min-max bottleneck
   barrier via graph traversal on generated steps.
"""

import pytest
import time
from pathlib import Path

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.conditions import ReactionConditions
from src.formulation import Formulation
from src.pipeline import MaillardPipeline
from src.precursor_resolver import resolve, resolve_many  # noqa: E402
from src.pathway_extractor import Species, ElementaryStep  # noqa: E402
from src.recommend import Recommender  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]

def test_resolver_exact_match():
    sp = resolve("L-Cysteine")
    assert sp.label == "L-Cysteine"
    assert "S" in sp.smiles

def test_resolver_fuzzy_match():
    sp1 = resolve("ribose")
    sp2 = resolve("d-ribose")
    assert sp1.label == "D-Ribose"
    assert sp1.smiles == sp2.smiles

def test_resolver_multiple():
    sps = resolve_many(["glucose", "glycine", "thiamine"])
    assert len(sps) == 3
    assert sps[0].label == "D-Glucose"
    assert sps[1].label == "Glycine"
    assert sps[2].label == "Thiamine (Vitamin B1)"

def test_resolver_unknown():
    with pytest.raises(ValueError, match="Unknown precursor"):
        resolve("nonexistent_sugar")

def test_recommender_predict_from_steps():
    # We mock an empty results file since we rely on dynamic barriers passed in
    mock_results = Path("mock_results.json")
    if not mock_results.exists():
        mock_results.write_text("[]")
        
    rec = Recommender(mock_results)
    
    # We want to form 2-Methyl-3-furanthiol (MFT), canonical SMILES matches desirable targets
    # Let's use the exact SMILES from the desirable_targets.yml: "Cc1occc1S"
    MFT = Species("MFT", "Cc1occc1S")
    FURFURAL = Species("Furfural", "O=Cc1ccco1") # not a target, intermediate
    RIBOSE = Species("Ribose", "OC[C@H]1OC(O)[C@H](O)[C@@H]1O")
    CYS = Species("Cysteine", "N[C@@H](CS)C(=O)O")
    
    # Pathway 1: Ribose -> Furfural -> MFT
    # Barriers: Ribose->Furfural = 30.0, Furfural->MFT = 15.0
    # Min-max barrier required to reach MFT = max(30, 15) = 30.0
    
    step1 = ElementaryStep(reactants=[RIBOSE], products=[FURFURAL])
    step2 = ElementaryStep(reactants=[FURFURAL, CYS], products=[MFT])
    
    # Build dictionary
    s1_key = f"{RIBOSE.smiles}->{FURFURAL.smiles}"
    s2_key = f"{CYS.smiles}+{FURFURAL.smiles}->{MFT.smiles}"
    # Python sorting makes C comes before O, we must sort the smiles:
    s2_key = f"{sorted([CYS.smiles, FURFURAL.smiles])[0]}+{sorted([CYS.smiles, FURFURAL.smiles])[1]}->{MFT.smiles}"
    
    barriers = {
        s1_key: 30.0,
        s2_key: 15.0
    }
    
    # Also add a slower pathway to a competing target like hexanal to test sorting
    HEXANAL = Species("Hexanal", "CCCCCC=O")
    LINOLEIC = Species("Linoleic", "CCCCC=CCC=CCCCCCCCC(=O)O") # Mock pre-cursor
    step3 = ElementaryStep(reactants=[LINOLEIC], products=[HEXANAL])
    s3_key = f"{LINOLEIC.smiles}->{HEXANAL.smiles}"
    barriers[s3_key] = 45.0
    
    steps = [step1, step2, step3]
    initial_pool = {RIBOSE.smiles: 1.0, CYS.smiles: 1.0, LINOLEIC.smiles: 1.0}
    
    predictions = rec.predict_from_steps(steps, barriers, initial_pool)["targets"]
    
    # Updated count: 3 (MFT, Furfural, and Hexanal)
    assert len(predictions) == 3
    
    # Barriers are: MFT (30.0), Furfural (30.0), Hexanal (45.0)
    # Sorting: Furfural seems to come before 2-Methyl-3-furanthiol (MFT) in this environment
    p1 = predictions[0]
    p2 = predictions[1]
    p3 = predictions[2]
    
    assert p1["name"] == "Furfural"
    assert p1["span"] == pytest.approx(30.0, abs=1e-5)
    
    assert p2["name"] == "2-Methyl-3-furanthiol (MFT)"
    assert p2["span"] == pytest.approx(30.0, abs=1e-5)
    
    assert p3["name"] == "Hexanal"
    assert p3["span"] == pytest.approx(45.0, abs=1e-5)

    # Cleanup mock if needed
    if mock_results.exists():
        mock_results.unlink()


def test_recommender_predict_from_steps_kinetic_mode_returns_metadata():
    rec = Recommender()

    furfural = Species("Furfural", "O=Cc1ccco1")
    ribose = Species("Ribose", "OC[C@H]1OC(O)[C@H](O)[C@@H]1O")
    step = ElementaryStep(
        reactants=[ribose],
        products=[furfural],
        reaction_family="furfural_formation",
    )
    step_key = f"{ribose.smiles}->{furfural.smiles}"

    result = rec.predict_from_steps(
        [step],
        {step_key: 12.0},
        {ribose.smiles: 1.0},
        time_minutes=5.0,
        prediction_mode="kinetic",
        reaction_conditions=ReactionConditions(prediction_mode="kinetic"),
    )

    assert result["kinetic_metadata"]["prediction_engine"] == "time_resolved_microkinetic"
    assert result["kinetic_metadata"]["successful"] is True
    assert result["projection_context"]["prediction_engine"] == "time_resolved_microkinetic"
    assert result["kinetic_metadata"]["reaction_count"] == 1
    assert result["kinetic_metadata"]["observable_surface"] == "end_state_default"


def test_pipeline_evaluate_single_kinetic_mode_sets_confidence_metadata():
    pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    formulation = Formulation(
        name="kinetic-smoke",
        sugars=["ribose"],
        amino_acids=["cysteine"],
        molar_ratios={"ribose": 1.0, "cysteine": 1.0},
        temperature=140.0,
        time_minutes=10.0,
    )

    result = pipeline.evaluate_single(
        formulation,
        ReactionConditions(pH=6.0, temperature_celsius=140.0, prediction_mode="kinetic"),
    )

    assert result.confidence_metadata["prediction_engine"] == "time_resolved_microkinetic"
    assert result.confidence_metadata["kinetics"]["successful"] is True


def test_pipeline_matrix_benchmark_supports_kinetic_mode():
    pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    bench = load_benchmark(ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json")
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)
    conditions.prediction_mode = "kinetic"

    result = pipeline.evaluate_single(formulation, conditions)

    assert result.confidence_metadata["prediction_engine"] == "time_resolved_microkinetic"
    assert result.confidence_metadata["kinetics"]["reaction_count"] >= 1
    assert result.predicted_ppb


def test_projection_mode_remains_unchanged_when_explicitly_requested():
    rec = Recommender()
    furfural = Species("Furfural", "O=Cc1ccco1")
    ribose = Species("Ribose", "OC[C@H]1OC(O)[C@H](O)[C@@H]1O")
    step = ElementaryStep(
        reactants=[ribose],
        products=[furfural],
        reaction_family="furfural_formation",
    )
    step_key = f"{ribose.smiles}->{furfural.smiles}"

    baseline = rec.predict_from_steps([step], {step_key: 18.0}, {ribose.smiles: 1.0}, time_minutes=5.0)
    explicit_projection = rec.predict_from_steps(
        [step],
        {step_key: 18.0},
        {ribose.smiles: 1.0},
        time_minutes=5.0,
        prediction_mode="projection",
    )

    assert explicit_projection["predicted_ppb"] == baseline["predicted_ppb"]
    assert explicit_projection["predicted_proxy_ppb"] == baseline["predicted_proxy_ppb"]


def test_pipeline_kinetic_mode_small_run_completes_quickly():
    pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    formulation = Formulation(
        name="kinetic-performance-smoke",
        sugars=["ribose"],
        amino_acids=["cysteine"],
        molar_ratios={"ribose": 1.0, "cysteine": 1.0},
        temperature=140.0,
        time_minutes=5.0,
    )

    start = time.perf_counter()
    result = pipeline.evaluate_single(
        formulation,
        ReactionConditions(pH=6.0, temperature_celsius=140.0, prediction_mode="kinetic"),
    )
    elapsed = time.perf_counter() - start

    assert result.confidence_metadata["kinetics"]["successful"] is True
    assert elapsed < 10.0
