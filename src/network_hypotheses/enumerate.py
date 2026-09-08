"""Depth-limited enumeration of what the rules propose from a charge of engine species."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Set, Tuple

from src.network_hypotheses import rules as R
from src.network_hypotheses import structures as S

#: A product with more heavy atoms than this is not a small volatile or intermediate; it is dropped.
MAX_HEAVY_ATOMS = 40


@dataclass(frozen=True)
class Step:
    rule_id: str
    reactants: Tuple[str, ...]   # canonical SMILES
    products: Tuple[str, ...]    # canonical SMILES
    depth: int


def _heavy(smiles: str) -> int:
    Chem = S._rdkit()
    m = Chem.MolFromSmiles(smiles)
    return m.GetNumHeavyAtoms() if m else 10**6


def enumerate_from(charge: Sequence[str], table: Mapping[str, S.Structure], rules: Sequence[R.Rule],
                   depth: int = 2, max_new_products: int = 400,
                   may_expand=None) -> Tuple[List[Step], Dict[str, int]]:
    """Apply every rule to every reactant tuple drawn from the current pool, `depth` times.

    `charge` lists engine species keys; lumps are skipped (no structure). The charge and its
    first-generation products always react on; products found at depth 2 or deeper seed further steps
    only when `may_expand(canonical_smiles)` is true (default: every product), so a run stays on the
    known map (engine species and registry compounds) instead of walking into products no source names. Returns the steps found and the depth at which each
    canonical SMILES first appeared (0 for the charge).
    """
    may_expand = may_expand or (lambda canon: True)
    pool: Dict[str, int] = {}
    terminal: Set[str] = set()      # products of terminal (sink / adduct) rules: recorded, never expanded
    for key in charge:
        st = table[key]
        if st.kind == "molecule":
            pool.setdefault(st.canonical, 0)
    steps: List[Step] = []
    seen_steps: Set[Tuple[str, Tuple[str, ...]]] = set()
    for d in range(1, depth + 1):
        # the charge and its first-generation products always react on; deeper generations only
        # when `may_expand` says the compound is known (an engine species or a registry compound)
        current = [s for s, first in pool.items() if first < d and s not in terminal and (first <= 1 or may_expand(s))]
        new: Dict[str, int] = {}
        for rule in rules:
            if rule.arity == 1:
                tuples = [(a,) for a in current]
            else:
                tuples = [(a, b) for i, a in enumerate(current) for b in current[i:]]
            for reactants in tuples:
                if (rule.id, reactants) in seen_steps:
                    continue
                seen_steps.add((rule.id, reactants))
                for products in R.apply(rule, reactants):
                    products = tuple(p for p in products if _heavy(p) <= MAX_HEAVY_ATOMS)
                    if not products:
                        continue
                    steps.append(Step(rule.id, reactants, products, d))
                    for p in products:
                        if p not in pool and p not in new:
                            new[p] = d
                        if rule.terminal and p not in pool:
                            terminal.add(p)
        if len(pool) + len(new) > max_new_products:
            raise RuntimeError(f"enumeration exceeded {max_new_products} products at depth {d}; narrow the charge or the rules")
        pool.update(new)
        if not new:
            break
    return steps, pool
