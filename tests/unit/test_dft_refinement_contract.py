from src.dft_refinement_contract import (
    DEFER,
    REJECT,
    RUN_NOW,
    DFTTargetCandidate,
    build_offline_dft_job,
    triage_dft_candidate,
)


def test_triage_rejects_low_value_candidates():
    candidate = DFTTargetCandidate(
        reaction_id="low_value_step",
        reaction_family="background",
        sensitivity_score=0.2,
        benchmark_impact_score=0.3,
        evidence_gap_score=0.5,
        current_barrier_source="heuristic",
        rationale="No benchmark-visible effect.",
    )

    decision = triage_dft_candidate(candidate)

    assert decision.decision == REJECT
    assert decision.priority_tier == "not_applicable"


def test_triage_defers_watchlist_candidates():
    candidate = DFTTargetCandidate(
        reaction_id="watchlist_step",
        reaction_family="sulfur",
        sensitivity_score=0.55,
        benchmark_impact_score=0.5,
        evidence_gap_score=0.55,
        current_barrier_source="heuristic",
        rationale="Interesting but not yet decisive.",
    )

    decision = triage_dft_candidate(candidate)

    assert decision.decision == DEFER
    assert decision.priority_tier == "watchlist"


def test_triage_runs_tier1_for_high_impact_barrier_refinement():
    candidate = DFTTargetCandidate(
        reaction_id="mft_forming_step",
        reaction_family="sulfur",
        sensitivity_score=0.9,
        benchmark_impact_score=0.85,
        evidence_gap_score=0.8,
        current_barrier_source="xtb",
        rationale="Dominates MFT ranking uncertainty.",
    )

    decision = triage_dft_candidate(candidate)

    assert decision.decision == RUN_NOW
    assert decision.priority_tier == "tier1_barrier_refinement"
    assert "surrogate_patch.json" in decision.required_outputs


def test_triage_runs_tier2_when_explicit_solvent_is_required():
    candidate = DFTTargetCandidate(
        reaction_id="water_mediated_pt",
        reaction_family="proton_transfer",
        sensitivity_score=0.92,
        benchmark_impact_score=0.78,
        evidence_gap_score=0.85,
        current_barrier_source="xtb",
        rationale="Water-mediated proton transfer flagged by xTB limitations.",
        requires_explicit_solvent=True,
        requires_irc=True,
    )

    decision = triage_dft_candidate(candidate)

    assert decision.decision == RUN_NOW
    assert decision.priority_tier == "tier2_explicit_solvent"
    assert "irc_path.xyz" in decision.required_outputs


def test_build_offline_dft_job_materializes_execution_contract():
    candidate = DFTTargetCandidate(
        reaction_id="fft_forming_step",
        reaction_family="sulfur",
        sensitivity_score=0.88,
        benchmark_impact_score=0.83,
        evidence_gap_score=0.75,
        current_barrier_source="heuristic",
        rationale="Expected to move FFT ranking.",
        charge=0,
        spin=0,
    )

    job = build_offline_dft_job(candidate, artifact_prefix="fft_batch")

    assert job.reaction_id == "fft_forming_step"
    assert job.priority_tier == "tier1_barrier_refinement"
    assert job.optimize_transition_state is True
    assert job.generate_quasi_harmonic is True
    assert job.artifact_prefix == "fft_batch"