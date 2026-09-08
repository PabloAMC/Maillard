"""Reaction rules: data/lit/reaction_rules.yml loaded, compiled and run, with their controls."""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Mapping, Optional, Sequence, Set, Tuple

from src import data_access, data_paths
from src.network_hypotheses import structures as S


@dataclass(frozen=True)
class Rule:
    id: str
    name: str
    klass: str
    status: str
    smirks: str
    source: Mapping[str, str]
    conditions: str
    controls: Mapping[str, list]
    arity: int
    terminal: bool = False


def load() -> List[Rule]:
    raw = data_access.load_yaml(data_paths.REACTION_RULES)["rules"]
    Chem = S._rdkit()
    from rdkit.Chem import AllChem

    out = []
    for r in raw:
        rxn = AllChem.ReactionFromSmarts(r["smirks"])
        out.append(Rule(r["id"], r["name"], r["class"], r["status"], r["smirks"], dict(r["source"]),
                        r.get("conditions", ""), dict(r.get("controls", {})), rxn.GetNumReactantTemplates(),
                        bool(r.get("terminal", False))))
    return out


def _compiled(rule: Rule):
    from rdkit.Chem import AllChem

    return AllChem.ReactionFromSmarts(rule.smirks)


def apply(rule: Rule, reactant_smiles: Sequence[str]) -> Set[Tuple[str, ...]]:
    """Every distinct product set (canonical SMILES, no stereo, sorted) the rule yields on these
    reactants, over every ordering of them. Products RDKit cannot sanitise are dropped: a rule that
    writes a broken molecule is a broken rule, and the controls catch it."""
    Chem = S._rdkit()
    rxn = _compiled(rule)
    if rxn.GetNumReactantTemplates() != len(reactant_smiles):
        return set()
    mols = [Chem.MolFromSmiles(s) for s in reactant_smiles]
    if any(m is None for m in mols):
        return set()
    out: Set[Tuple[str, ...]] = set()
    for perm in set(itertools.permutations(range(len(mols)))):
        try:
            product_sets = rxn.RunReactants(tuple(mols[i] for i in perm))
        except Exception:
            continue
        for pset in product_sets:
            names = []
            ok = True
            for p in pset:
                try:
                    Chem.SanitizeMol(p)
                except Exception:
                    ok = False
                    break
                names.append(Chem.MolToSmiles(p, isomericSmiles=False))
            if ok and names:
                out.add(tuple(sorted(names)))
    return out


def run_controls(rule: Rule, table: Mapping[str, S.Structure]) -> List[str]:
    """Failures of the rule's own positive and negative controls (empty when the rule is sound)."""
    failures: List[str] = []

    def smiles_of(key: str) -> str:
        st = table.get(key)
        if st is None or st.kind != "molecule":
            raise KeyError(f"{rule.id}: control names {key!r}, which is not a structured species")
        return st.smiles

    for ctl in rule.controls.get("positive", []):
        got = apply(rule, [smiles_of(k) for k in ctl["reactants"]])
        flat = {p for ps in got for p in ps}
        if ctl.get("fires"):
            if not got:
                failures.append(f"{rule.id}: positive control {ctl['reactants']} did not fire")
            continue
        want = {table[k].canonical for k in ctl.get("products", [])}
        want |= {S.canonical(s) for s in ctl.get("products_smiles", [])}
        missing = want - flat
        if missing:
            failures.append(f"{rule.id}: positive control {ctl['reactants']} lacks {sorted(missing)}; got {sorted(flat)}")
    for ctl in rule.controls.get("negative", []):
        got = apply(rule, [smiles_of(k) for k in ctl["reactants"]])
        if got:
            failures.append(f"{rule.id}: negative control {ctl['reactants']} fired: {sorted(got)[:3]}")
    return failures
