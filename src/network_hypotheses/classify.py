"""Place each proposed step and product against the engine's network and the compound registry."""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from src.network_hypotheses import structures as S
from src.network_hypotheses.enumerate import Step

#: Engine species that stand for bookkeeping pools rather than reagents; a proposed step is matched to
#: an engine reaction ignoring these on the engine side.
BOOKKEEPING = {"OX", "OXV", "OXR", "FRAG_C", "FRAG_N", "FRAG_S", "ACID", "CBX", "MEL_C", "MEL_N", "LIPID_FRAG_C"}

PLACEMENTS = {
    "modelled": "the engine has a reaction with these reactants and this product, with a rate",
    "mechanism_known": "the rule's source draws this step; the engine has no reaction for it",
    "proposed": "analogous to a cited rule, on reactants no source shows; no rate, no source",
}


def engine_reactions() -> List[Tuple[str, str, Set[str], Set[str]]]:
    """(lane, key, reactant keys, product keys) for every reaction of every lane, read-only."""
    from src.kinetic_core import acrylamide, lipid, network, sulfur

    out = []
    seen = set()
    for lane, reactions in (("trunk", getattr(network, "TRUNK_REACTIONS", network.REACTIONS)),
                            ("sulfur", sulfur.SULFUR_REACTIONS),
                            ("acrylamide", getattr(acrylamide, "FULL_ACRYLAMIDE_REACTIONS", acrylamide.ACRYLAMIDE_REACTIONS))):
        for r in reactions:
            if r.key in seen:
                continue
            seen.add(r.key)
            out.append((lane, r.key, set(r.reactants) - BOOKKEEPING, set(r.products) - BOOKKEEPING))
    # The lipid lane is a branch model, not a reaction list: each hydroperoxide pool splits into the
    # products its position can make (species_lipid.POSITION_PRODUCTS). Represented as one pseudo-step
    # per pool so the layer can recognise the lane's own products.
    from src.kinetic_core import species_lipid

    for pool, (position, _geometry) in species_lipid.LOOH_POOLS.items():
        out.append(("lipid", f"lipid_scission_{pool}", {pool}, set(species_lipid.POSITION_PRODUCTS[position])))
    return out


def species_by_canonical(table: Mapping[str, S.Structure]) -> Dict[str, str]:
    return {st.canonical: key for key, st in table.items() if st.kind == "molecule"}


def place_product(canon: str, by_canon: Mapping[str, str], registry: Mapping[str, str]) -> Tuple[str, Optional[str]]:
    """('species', key) | ('registry', compound id) | ('new', None)."""
    if canon in by_canon:
        return "species", by_canon[canon]
    if canon in registry:
        return "registry", registry[canon]
    return "new", None


def place_step(step: Step, rule_status: str, by_canon: Mapping[str, str],
               reactions: Sequence[Tuple[str, str, Set[str], Set[str]]]) -> Tuple[str, List[str]]:
    """The placement and the engine reaction keys it matches (empty unless modelled)."""
    r_keys = {by_canon[c] for c in step.reactants if c in by_canon}
    p_keys = {by_canon[c] for c in step.products if c in by_canon}
    if all(c in by_canon for c in step.reactants) and p_keys:
        matches = [key for lane, key, rs, ps in reactions if r_keys <= rs and (p_keys & ps)]
        if matches:
            return "modelled", sorted(matches)
        # a rule that lumps two engine steps (e.g. the Amadori compound through the Schiff base):
        # reactants into step A, one of A's products into step B, B makes the product
        via = []
        for lane_a, key_a, rs_a, ps_a in reactions:
            if not r_keys <= rs_a:
                continue
            for lane_b, key_b, rs_b, ps_b in reactions:
                if key_b != key_a and (ps_a & rs_b) and (p_keys & ps_b):
                    via.append(f"{key_a} > {key_b}")
        if via:
            return "modelled", sorted(via)
    return ("mechanism_known" if rule_status in ("established", "net") else "proposed"), []
