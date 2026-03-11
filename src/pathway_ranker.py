"""
src/pathway_ranker.py — Parallel pathway screening execution and ranking.

Takes a list of reaction pathways (lists of elementary steps),
evaluates them in parallel via XTBScreener, and ranks them by 
maximum rate-limiting barrier height.
"""

from typing import List, Tuple, Optional
from dataclasses import dataclass
from multiprocessing import Pool

from src.pathway_extractor import ElementaryStep  # noqa: E402
from src.xtb_screener import XTBScreener  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402

@dataclass
class PathwayProfile:
    pathway_name: str
    steps: List[ElementaryStep]
    deltaE_kcal_list: List[float]
    barrier_kcal_list: List[float]
    scaled_rates: Optional[List[float]] = None
    
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
    screener = XTBScreener()
    try:
        return screener.compute_reaction_energy(step)
    except Exception:
        # Penalize drastically if evaluation fails (usually RDKit embedding failure for tricky intermediates)
        return (999.0, 999.0)

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
        
    def screen_pathways(self, pathways: dict[str, List[ElementaryStep]]) -> List[PathwayProfile]:
        """
        Evaluate multiple pathways in parallel and return them sorted
        from lowest rate-limiting barrier/span (fastest) to highest (slowest).
        """
        jobs = [(name, steps, self.conditions) for name, steps in pathways.items()]
        
        if self.n_cores > 1:
            with Pool(processes=self.n_cores) as pool:
                profiles = pool.map(_evaluate_pathway, jobs)
        else:
            profiles = [_evaluate_pathway(job) for job in jobs]
            
        # Sort by the energetic span, fallback to old max barrier
        profiles.sort(key=lambda p: p.energetic_span if p.energetic_span > 0 else p.rate_limiting_barrier)
        
        return profiles
