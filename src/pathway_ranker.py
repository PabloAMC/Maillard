"""
src/pathway_ranker.py — Parallel pathway screening execution and ranking.

Takes a list of reaction pathways (lists of elementary steps), evaluates them in
parallel via XTBScreener, and ranks them.

Two rankings are produced and explicitly labelled (2026-08-27, audit item 3.3):

* `by_energetic_span` — intrinsic, condition-free: the Kozuch-Shaik energetic
  span, lowest first.
* `by_conditioned_rate` — the bottleneck conditioned rate (pH x Arrhenius x
  water-activity multipliers), fastest first.  Only available when the ranker
  was given ReactionConditions.

The old single sort used the span alone, so the pH / temperature / water
activity multipliers were computed for every step and then thrown away: two
pathways with the same span ranked identically at pH 4 and pH 9.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from multiprocessing import Pool

from src.pathway_extractor import ElementaryStep  # noqa: E402
from src.xtb_screener import XTBScreener  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402

# Sentinel used when a step cannot be evaluated at all (matches the value the
# sibling failure branch in evaluate_single_step already returned).
FAILED_EVALUATION_SENTINEL_KCAL = 999.0

@dataclass
class PathwayProfile:
    pathway_name: str
    steps: List[ElementaryStep]
    deltaE_kcal_list: List[float]
    barrier_kcal_list: List[float]
    scaled_rates: Optional[List[float]] = None

    @property
    def has_steps(self) -> bool:
        return bool(self.barrier_kcal_list)

    @property
    def status(self) -> str:
        """Explicit evaluation status.

        An empty pathway used to report span 0.0 and therefore sorted FIRST —
        i.e. "no steps at all" was ranked as the fastest possible route.
        """
        if not self.has_steps:
            return "no_steps"
        if any(barrier >= FAILED_EVALUATION_SENTINEL_KCAL for barrier in self.barrier_kcal_list):
            return "evaluation_failed"
        return "ok"

    @property
    def bottleneck_scaled_rate(self) -> Optional[float]:
        """Slowest conditioned step rate, i.e. the pathway's rate bottleneck.

        None when no conditions were supplied (nothing to condition on).
        """
        if not self.scaled_rates:
            return None
        return min(float(rate) for rate in self.scaled_rates)

    @property
    def rate_limiting_barrier(self) -> float:
        if not self.barrier_kcal_list:
            return 0.0
        return max(self.barrier_kcal_list)

    @property
    def energetic_span(self) -> float:
        if not self.barrier_kcal_list:
            return 0.0
            
        max_span = 0.0
        current_energy = 0.0
        lowest_intermediate = 0.0
        
        for i in range(len(self.barrier_kcal_list)):
            # TS energy relative to start
            ts_energy = current_energy + self.barrier_kcal_list[i]
            span = ts_energy - lowest_intermediate
            if span > max_span:
                max_span = span
                
            # Update intermediate baseline
            current_energy += self.deltaE_kcal_list[i]
            if current_energy < lowest_intermediate:
                lowest_intermediate = current_energy
                
        return max_span
        
    @property
    def overall_thermodynamics(self) -> float:
        return sum(self.deltaE_kcal_list)
        
    def __str__(self) -> str:
        if not self.has_steps:
            return f"{self.pathway_name} | status: no_steps (not ranked as a viable route)"
        s = f"{self.pathway_name} | Max ΔE‡: {self.rate_limiting_barrier:.1f} kcal/mol | Energetic Span: {self.energetic_span:.1f} kcal/mol | ΔErxn: {self.overall_thermodynamics:.1f} kcal/mol"
        if self.scaled_rates:
            try:
                min_rate = min(r for r in self.scaled_rates if r > 0)
                s += f"\n  (Bottleneck scaled rate: {min_rate:.2e})"
            except ValueError:
                pass
        return s


def evaluate_single_step(step: ElementaryStep) -> Tuple[float, float]:
    """Evaluates a single elementary step via xTB (used for parallel mapping)."""
    # Steps whose species could not be resolved from the mechanism source are
    # not evaluable: computing an energy from a partial reactant set would be
    # arithmetic on a reaction that does not exist.
    if getattr(step, "source_quality", "") == "unresolved_species":
        return (FAILED_EVALUATION_SENTINEL_KCAL, FAILED_EVALUATION_SENTINEL_KCAL)

    if getattr(step, "barrier_kcal_mol", None) is not None:
        try:
            screener = XTBScreener()
            # Fast thermodynamic calculation (only optimize separate species, no NEB)
            reactants_total_E = sum(screener.optimize_species(r.smiles).energy_hartree for r in step.reactants)
            products_total_E = sum(screener.optimize_species(p.smiles).energy_hartree for p in step.products)
            delta_E = (products_total_E - reactants_total_E) * 627.509
        except Exception:
            # A failed ΔE is NOT a thermoneutral step: 0.0 silently claimed the
            # most favourable possible thermodynamics. Use the same honest
            # penalty sentinel as the branch below.
            delta_E = FAILED_EVALUATION_SENTINEL_KCAL
        return (delta_E, step.barrier_kcal_mol)

    screener = XTBScreener()
    try:
        return screener.compute_reaction_energy(step)
    except Exception:
        # Penalize drastically if evaluation fails (usually RDKit embedding failure for tricky intermediates)
        return (FAILED_EVALUATION_SENTINEL_KCAL, FAILED_EVALUATION_SENTINEL_KCAL)

def _evaluate_pathway(args: Tuple[str, List[ElementaryStep], Optional[ReactionConditions]]) -> PathwayProfile:
    """Helper to evaluate a full pathway sequence."""
    name, steps, conditions = args
    dEs, barriers, scaled_rates = [], [], []
    for step in steps:
        dE, bar = evaluate_single_step(step)
        dEs.append(dE)
        barriers.append(bar)
        
        if conditions:
            ph_mult = conditions.get_ph_multiplier(step.reaction_family or "")
            arrh_mult = conditions.get_arrhenius_multiplier(bar)
            aw_mult = conditions.get_water_activity_multiplier()
            scaled_rates.append(ph_mult * arrh_mult * aw_mult)
        
    return PathwayProfile(
        pathway_name=name,
        steps=steps,
        deltaE_kcal_list=dEs,
        barrier_kcal_list=barriers,
        scaled_rates=scaled_rates
    )

class PathwayRanker:
    """Evaluates multiple pathways and ranks them."""
    
    def __init__(self, n_cores: int = 4, conditions: Optional[ReactionConditions] = None):
        self.n_cores = n_cores
        self.conditions = conditions
        
    def evaluate_pathways(self, pathways: dict[str, List[ElementaryStep]]) -> List[PathwayProfile]:
        """Evaluate pathways (in parallel when n_cores > 1), unranked."""
        jobs = [(name, steps, self.conditions) for name, steps in pathways.items()]

        if self.n_cores > 1:
            with Pool(processes=self.n_cores) as pool:
                return list(pool.map(_evaluate_pathway, jobs))
        return [_evaluate_pathway(job) for job in jobs]

    @staticmethod
    def rank_by_energetic_span(profiles: List[PathwayProfile]) -> List[PathwayProfile]:
        """Intrinsic ranking: lowest energetic span first. Empty pathways last."""
        return sorted(
            profiles,
            key=lambda p: (
                0 if p.has_steps else 1,
                p.energetic_span if p.energetic_span > 0 else p.rate_limiting_barrier,
                p.pathway_name,
            ),
        )

    @staticmethod
    def rank_by_conditioned_rate(profiles: List[PathwayProfile]) -> List[PathwayProfile]:
        """Conditioned ranking: highest bottleneck scaled rate (fastest) first.

        Profiles with no conditioned rates (no conditions supplied, or no steps)
        rank after every profile that has one — they cannot be compared on rate.
        """
        def key(profile: PathwayProfile):
            rate = profile.bottleneck_scaled_rate
            if not profile.has_steps or rate is None:
                return (1, 0.0, profile.pathway_name)
            # Negated so that a larger rate sorts first.
            return (0, -float(rate), profile.pathway_name)

        return sorted(profiles, key=key)

    def rank_pathways(self, pathways: dict[str, List[ElementaryStep]]) -> Dict[str, List[PathwayProfile]]:
        """Return BOTH explicitly-labelled rankings for the same evaluation."""
        profiles = self.evaluate_pathways(pathways)
        return {
            "by_energetic_span": self.rank_by_energetic_span(profiles),
            "by_conditioned_rate": self.rank_by_conditioned_rate(profiles),
        }

    def screen_pathways(self, pathways: dict[str, List[ElementaryStep]]) -> List[PathwayProfile]:
        """
        Evaluate multiple pathways and return them ranked fastest-first.

        With ReactionConditions supplied the ranking is the CONDITIONED one, so
        pH / temperature / water activity actually change the order (previously
        the scaled rates were computed and then ignored). Without conditions
        there is nothing to condition on and the intrinsic energetic-span
        ranking is returned. Pathways with no steps always rank last, with
        `status == "no_steps"` (they used to rank first on a span of 0.0).
        """
        profiles = self.evaluate_pathways(pathways)
        if self.conditions is not None:
            return self.rank_by_conditioned_rate(profiles)
        return self.rank_by_energetic_span(profiles)
