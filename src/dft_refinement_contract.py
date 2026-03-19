from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


RUN_NOW = "run_now"
DEFER = "defer"
REJECT = "reject"


@dataclass(frozen=True)
class DFTTargetCandidate:
    reaction_id: str
    reaction_family: str
    sensitivity_score: float
    benchmark_impact_score: float
    evidence_gap_score: float
    current_barrier_source: str
    rationale: str
    requires_explicit_solvent: bool = False
    requires_irc: bool = False
    current_confidence: str = "heuristic"
    charge: int = 0
    spin: int = 0


@dataclass(frozen=True)
class DFTTriageDecision:
    reaction_id: str
    decision: str
    priority_tier: str
    rationale: str
    required_outputs: List[str]
    surrogate_plan: str


@dataclass(frozen=True)
class OfflineDFTJob:
    reaction_id: str
    priority_tier: str
    solvent_name: str
    temperature_k: float
    charge: int
    spin: int
    use_explicit_solvent: bool
    n_water: int
    optimize_transition_state: bool
    generate_irc: bool
    generate_quasi_harmonic: bool
    expected_outputs: List[str] = field(default_factory=list)
    artifact_prefix: str = "dft_batch"


def triage_dft_candidate(candidate: DFTTargetCandidate) -> DFTTriageDecision:
    if candidate.sensitivity_score < 0.45 or candidate.benchmark_impact_score < 0.45:
        return DFTTriageDecision(
            reaction_id=candidate.reaction_id,
            decision=REJECT,
            priority_tier="not_applicable",
            rationale="Low sensitivity or low benchmark impact. Keep this reaction on the cheap-first surrogate path.",
            required_outputs=[],
            surrogate_plan="No DFT. Preserve current surrogate until a benchmark-visible gap appears.",
        )

    weighted_priority = (
        0.4 * candidate.sensitivity_score
        + 0.4 * candidate.benchmark_impact_score
        + 0.2 * candidate.evidence_gap_score
    )

    if weighted_priority < 0.6:
        return DFTTriageDecision(
            reaction_id=candidate.reaction_id,
            decision=DEFER,
            priority_tier="watchlist",
            rationale="Potentially relevant, but not strong enough to justify expensive QM ahead of matrix priors or new experiments.",
            required_outputs=["priority_note.md"],
            surrogate_plan="Track in the watchlist and revisit after benchmark or sensitivity updates.",
        )

    priority_tier = "tier2_explicit_solvent" if candidate.requires_explicit_solvent else "tier1_barrier_refinement"
    outputs = [
        "optimized_reactant.xyz",
        "optimized_ts.xyz",
        "barrier_summary.json",
        "quasi_harmonic_thermo.json",
        "surrogate_patch.json",
    ]
    if candidate.requires_irc:
        outputs.append("irc_path.xyz")

    if candidate.requires_explicit_solvent:
        surrogate_plan = "Run a very small explicit-water refinement, then compress the result into a reusable correction delta rather than requiring repeat QM."
    else:
        surrogate_plan = "Refine only the decisive barrier and write the result back as a cached correction artifact for the cheap daily workflow."

    return DFTTriageDecision(
        reaction_id=candidate.reaction_id,
        decision=RUN_NOW,
        priority_tier=priority_tier,
        rationale=(
            "High-sensitivity, benchmark-relevant reaction with insufficient literature support. "
            "Eligible for selective offline DFT because it can plausibly move observable ranking or calibration."
        ),
        required_outputs=outputs,
        surrogate_plan=surrogate_plan,
    )


def build_offline_dft_job(
    candidate: DFTTargetCandidate,
    *,
    artifact_prefix: str,
    solvent_name: str = "water",
    temperature_k: float = 423.15,
) -> OfflineDFTJob:
    decision = triage_dft_candidate(candidate)
    if decision.decision != RUN_NOW:
        raise ValueError(f"Candidate {candidate.reaction_id} is not approved for DFT execution")

    return OfflineDFTJob(
        reaction_id=candidate.reaction_id,
        priority_tier=decision.priority_tier,
        solvent_name=solvent_name,
        temperature_k=float(temperature_k),
        charge=int(candidate.charge),
        spin=int(candidate.spin),
        use_explicit_solvent=bool(candidate.requires_explicit_solvent),
        n_water=3 if candidate.requires_explicit_solvent else 0,
        optimize_transition_state=True,
        generate_irc=bool(candidate.requires_irc),
        generate_quasi_harmonic=True,
        expected_outputs=decision.required_outputs,
        artifact_prefix=artifact_prefix,
    )