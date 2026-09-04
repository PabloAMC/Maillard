"""The data wishlist: what new measurement would improve this model most, and what each would unlock.

2026-09-04 (usability pass). Until now the answer to "what should we measure next?" was spread
over four artifacts -- the slice profile (which fitted coordinates the data do not identify),
the scorecard (which rows the engine declares not evaluable and which targets it refuses), the
directional scorecard (which axes are thin) and the value-of-information ranking (where the
envelope misses most). This module reads all four and writes ONE artifact,
``results/validation/data_wishlist.{json,md}``, that a scientist can read top to bottom.

Nothing here is a measurement or a fit. Every line is derived from a tracked artifact and says
which artifact; the "unlocks" sentences are structural (which observables sit downstream of a
step in the network, which rows leave the not-evaluable list), not forecasts of accuracy.
"""
from __future__ import annotations

import re
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from src import data_access, data_paths, provenance
from src.directional_reliability import COIN, MIN_EVALUABLE_FOR_TRUST, TRUST_MIN_RATE, VERDICT_TRUST, verdict_for, wilson_lower

PROFILE_PATH = data_paths.VALIDATION_DIR / "kinetic_core_b9_profile.json"
LAPLACE_PATH = data_paths.VALIDATION_DIR / "kinetic_core_b9_laplace_covariance.json"
RANKING_PATH = data_paths.VALIDATION_DIR / "experiment_value_ranking.json"
SULFUR_MODULE = data_paths.REPO_ROOT / "src" / "kinetic_core" / "sulfur.py"
OUTPUT_JSON = data_paths.VALIDATION_DIR / "data_wishlist.json"

#: Sulfur-lane species a scientist can measure (the panel's targets), in the lane's keys.
SULFUR_OBSERVABLES = ("MFT", "FFT", "FUR", "MFTD", "MESH", "ACTZ", "NF", "H2S")
#: Charged reactants, in the lane's keys, so "from a charge containing X" can be said.
SULFUR_CHARGES = {"Glc": "glucose", "Fru": "fructose", "PENT": "a pentose (ribose / xylose)",
                  "Cys": "cysteine", "THI": "thiamine", "H2S": "hydrogen sulfide"}
#: Coordinates whose story is already told elsewhere in the repository; the wishlist quotes it.
KNOWN_STORIES = {
    "k_glc_ha": "The hexose fragmentation entry (the primary-evidence refit left it on its band floor; the engine declares "
                "'HEXOSE ENTRY UNIDENTIFIED' on every hexose-only thiol request).",
    "Ea_decay_thiol_sink": "Pressed against Gigl 2021's measured covalent-capture ceiling (7-102 kJ/mol); "
                           "a measurement of the sink's barrier at two temperatures decides whether the "
                           "band or the model is wrong.",
    "ph_acid_yield_per_sink_event": "The pH-trajectory model's one calibrated fraction, at its floor: the "
                                    "acid yield of the pentose deoxyosone sink has never been measured.",
}


# ---------------------------------------------------------------------------- the network
_REACTION_RX = re.compile(
    r'Reaction\(\s*"(?P<name>\w+)",\s*(?P<react>\{[^}]*\}),\s*(?P<prod>\{[^}]*\}),\s*"(?P<key>k_\w+)"'
    r'(?:,\s*"(?P<doc>(?:[^"\\]|\\.)*)")?', re.S,
)
#: Bookkeeping pools a bench scientist does not quantify; dropped from measurement sentences.
_POOLS = {"FRAG_C", "FRAG_S", "FRAG_N", "CBX", "ACID", "OLG", "OX"}
#: Decay-family activation energies are shared by every step of the family, not one reaction.
_DECAY_FAMILY_TEXT = {
    "thiol_sink": "every step that removes a thiol into the matrix sink (the k_thiol_decay / k_dimer_decay steps)",
    "carbonyl_sink": "every step that removes a carbonyl intermediate into the sink",
}


def _parse_dict(text: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for m in re.finditer(r'"(\w+)":\s*([\d.]+)', text):
        out[m.group(1)] = int(float(m.group(2)))
    return out


def sulfur_reactions(source: Optional[str] = None) -> List[Dict[str, Any]]:
    """The sulfur lane's reaction table, parsed from its source so the wishlist cannot drift
    from the network: ``[{name, reactants, products, key, doc}]``."""
    text = source if source is not None else SULFUR_MODULE.read_text(encoding="utf-8")
    out = []
    for m in _REACTION_RX.finditer(text):
        doc = (m.group("doc") or "").replace('"\n        "', "").strip()
        out.append({
            "name": m.group("name"), "reactants": _parse_dict(m.group("react")),
            "products": _parse_dict(m.group("prod")), "key": m.group("key"),
            "doc": re.sub(r"\s+", " ", doc)[:200],
        })
    return out


def downstream_observables(reactions: Sequence[Mapping[str, Any]], species: Sequence[str]) -> List[str]:
    """Observables reachable from ``species`` through the reaction graph (products of products)."""
    edges: Dict[str, Set[str]] = defaultdict(set)
    for r in reactions:
        for a in r["reactants"]:
            edges[a].update(r["products"])
    seen: Set[str] = set(species)
    queue = deque(species)
    while queue:
        node = queue.popleft()
        for nxt in edges.get(node, ()):
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return [o for o in SULFUR_OBSERVABLES if o in seen]


def upstream_charges(reactions: Sequence[Mapping[str, Any]], species: Sequence[str]) -> List[str]:
    """Charged precursors from which ``species`` can be reached."""
    edges: Dict[str, Set[str]] = defaultdict(set)
    for r in reactions:
        for p in r["products"]:
            edges[p].update(r["reactants"])
    seen: Set[str] = set(species)
    queue = deque(species)
    while queue:
        node = queue.popleft()
        for nxt in edges.get(node, ()):
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return [c for c in SULFUR_CHARGES if c in seen]


# ---------------------------------------------------------------------------- the sections
def _measurement_for(kind: str, reaction: Mapping[str, Any], key: str = "") -> str:
    reactants = " + ".join(SULFUR_CHARGES.get(k, k) for k in reaction["reactants"]) or "the step's reactants"
    named = [k for k in reaction["products"] if k not in _POOLS]
    products = " + ".join(named) if named else "its loss"
    if key.startswith("Ea_decay_"):
        family = key[len("Ea_decay_"):]
        return (f"the rate of {_DECAY_FAMILY_TEXT.get(family, f'the {family} decay family')} at two "
                f"temperatures in the same matrix (e.g. FFT loss at 30 and 80 C, Hofmann 2002-style), so the "
                f"barrier is measured rather than banded")
    if kind == "Ea":
        return (f"the same fed-intermediate experiment at two temperatures (e.g. 100 and 145 C): "
                f"{reactants} -> {products}, timed, so the barrier is measured rather than banded")
    if kind == "acid_yield":
        return "titrate the acid formed per mole of pentose deoxyosone consumed (pH before/after a timed heat step)"
    if kind == "pka":
        return "a measured pKa of the Amadori secondary ammonium at reaction temperature"
    return (f"a fed-intermediate rate: {reactants} heated alone, {products} quantified against time "
            f"(SIDA where a labelled standard exists), at the lane's reference 145 C")


def unidentified_coordinates(profile: Mapping[str, Any], laplace: Mapping[str, Any],
                             reactions: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """Fitted coordinates the primary evidence does not pin: flat or bound-limited in the profile,
    or unsampled in the Laplace covariance."""
    unsampled = {c["key"]: c for c in laplace.get("coordinates", []) if c.get("sampled") is False}
    out = []
    for row in profile.get("rows", []):
        verdict = row.get("verdict")
        if verdict not in ("flat", "bound_limited") and row["key"] not in unsampled:
            continue
        steps = [r for r in reactions if r["key"] == row["key"]]
        if not steps:
            reaction: Dict[str, Any] = {"name": None, "reactants": {}, "products": {}, "doc": ""}
        elif len(steps) == 1:
            reaction = dict(steps[0])
        else:
            # one constant shared by several steps (e.g. k_thiol_decay): the measurement is of the
            # family, so every reactant and product is named and the step label lists them all
            reaction = {
                "name": f"{len(steps)} steps sharing the constant: " + ", ".join(s["name"] for s in steps),
                "reactants": {k: 1 for s in steps for k in s["reactants"]},
                "products": {k: 1 for s in steps for k in s["products"]},
                "doc": next((s["doc"] for s in steps if s["doc"]), ""),
            }
        family_steps = [r for r in reactions if row["key"].startswith("Ea_decay_") and r["key"].endswith("_decay")]
        if reaction["name"]:
            obs = downstream_observables(reactions, list(reaction["products"]))
            charges = upstream_charges(reactions, list(reaction["reactants"]))
        elif family_steps:
            consumed = sorted({k for r in family_steps for k in r["reactants"]})
            obs = [o for o in SULFUR_OBSERVABLES if o in consumed] or ["every thiol level outside 100-145 C"]
            charges = upstream_charges(reactions, consumed)
        else:
            obs, charges = [], []
        out.append({
            "key": row["key"], "kind": row.get("kind"), "profile_verdict": verdict,
            "sigma_used": row.get("sigma_used"), "slice_min_k": row.get("slice_min_k"),
            "reaction": reaction["name"],
            "reactants": list(reaction["reactants"]), "products": list(reaction["products"]),
            "what_the_code_says": reaction["doc"][:200],
            "story": KNOWN_STORIES.get(row["key"], ""),
            "measurement": _measurement_for(str(row.get("kind")), reaction, row["key"]),
            "unlocks_observables": obs,
            "unlocks_from_charges": [SULFUR_CHARGES[c] for c in charges],
        })
    order = {"flat": 0, "bound_limited": 1}
    out.sort(key=lambda r: (order.get(r["profile_verdict"], 2), r["key"]))
    return out


def not_evaluable_rows(scorecard: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Rows the engine answers but declares as running through an unidentified route."""
    rows = [r for r in scorecard.get("refused_compounds", []) if r.get("unidentified_route")]
    return [{"benchmark_id": r["benchmark_id"], "compound": r["compound"], "panel": r["panel"]} for r in rows]


def refused_targets(scorecard: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """What the panel asks for that no lane can represent, grouped by the engine's own reason."""
    groups: Dict[str, Dict[str, Any]] = {}
    for r in scorecard.get("refused_compounds", []):
        if r.get("unidentified_route"):
            continue
        reason = str(r.get("reason", ""))
        head = reason.split(" -- ")[0][:120]
        g = groups.setdefault(head, {"reason": reason[:400], "rows": []})
        g["rows"].append({"benchmark_id": r["benchmark_id"], "compound": r["compound"]})
    out = [{"what": k, "reason": v["reason"], "rows": v["rows"], "row_count": len(v["rows"])} for k, v in groups.items()]
    out.sort(key=lambda g: -g["row_count"])
    return out


def _claims_to_trust(agree: int, evaluable: int) -> Optional[int]:
    """Smallest number of ADDITIONAL agreeing independent claims that would make the axis 'trust'."""
    for extra in range(0, 60):
        if verdict_for(agree + extra, evaluable + extra) == VERDICT_TRUST:
            return extra
    return None


def thin_axes(directional: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Directional axes that are not 'trust': how many claims exist, how many are not evaluable
    and why, and how many more agreeing claims would be needed."""
    ind = directional["summary"]["independent"]
    out = []
    for axis, b in ind["by_category"].items():
        agree, evaluable = int(b["agree"]), int(b["evaluable"])
        verdict = verdict_for(agree, evaluable) if evaluable else "no evaluable claim"
        if verdict == VERDICT_TRUST:
            continue
        out.append({
            "axis": axis, "agree": agree, "evaluable": evaluable,
            "not_evaluable": int(b.get("not_evaluable", 0)), "verdict": verdict,
            "wilson_lower_bound": round(wilson_lower(agree, evaluable), 3) if evaluable else None,
            "additional_agreeing_claims_to_trust": _claims_to_trust(agree, evaluable),
            "structural_block": (
                "no lane carries a water-activity term: the engine refuses every a_w comparison, so "
                "claims cannot help until a moisture-dependent step is measured and fitted"
                if axis == "moisture_aw" else None
            ),
        })
    out.sort(key=lambda r: (r["evaluable"], r["axis"]))
    return out


def largest_misses(ranking: Optional[Mapping[str, Any]], top: int = 10) -> List[Dict[str, Any]]:
    if not ranking:
        return []
    rows = ranking.get("candidates") or []
    return [{"benchmark_id": r.get("benchmark_id"), "compound": r.get("compound"),
             "voi_score": r.get("voi_score"), "envelope_miss_log10": r.get("envelope_miss_log10"),
             "template": r.get("suggested_doe_template")} for r in rows[:top]]


def build(*, scorecard_path: Path = data_paths.CORE_PANEL_SCORES, profile_path: Path = PROFILE_PATH,
          laplace_path: Path = LAPLACE_PATH, directional_path: Path = data_paths.CORE_DIRECTIONAL_SCORES,
          ranking_path: Path = RANKING_PATH, sulfur_source: Optional[str] = None) -> Dict[str, Any]:
    scorecard = data_access.load_json(scorecard_path)
    profile = data_access.load_json(profile_path)
    laplace = data_access.load_json(laplace_path)
    directional = data_access.load_json(directional_path)
    ranking = data_access.load_json(ranking_path, missing_ok=True)
    reactions = sulfur_reactions(sulfur_source)
    coords = unidentified_coordinates(profile, laplace, reactions)
    nev = not_evaluable_rows(scorecard)
    refused = refused_targets(scorecard)
    axes = thin_axes(directional)
    misses = largest_misses(ranking)
    unlocked = []
    for c in coords:
        if c["unlocks_observables"]:
            unlocked.append({
                "if_measured": c["key"],
                "you_could_predict": (
                    f"absolute {', '.join(c['unlocks_observables'])} from charges containing "
                    f"{', '.join(c['unlocks_from_charges']) or 'its reactants'}, as a fitted value instead of a band artefact"
                    + (f"; {len(nev)} panel rows now NOT EVALUABLE would return to the scorecard" if c["key"] == "k_glc_ha" and nev else "")
                ),
            })
    for a in axes:
        if a["additional_agreeing_claims_to_trust"] is not None and not a["structural_block"]:
            unlocked.append({
                "if_measured": f"{a['additional_agreeing_claims_to_trust']} more independent ordering claims on {a['axis']} (all agreeing)",
                "you_could_predict": f"the direction of a change along {a['axis']} with a 'trust' tag instead of '{a['verdict']}'",
            })
    payload = {
        "artifact": "data_wishlist",
        "provenance": provenance.provenance_block(
            "data_wishlist", generated_by="src/data_wishlist.py",
            inputs=[scorecard_path, profile_path, laplace_path, directional_path, ranking_path, SULFUR_MODULE],
        ),
        "how_to_read": (
            "Sections in order of leverage. 1: fitted coordinates the primary evidence does not pin -- one "
            "fed-intermediate measurement each turns a band artefact into a fitted value. 2: panel rows the "
            "engine answers but declares not evaluable. 3: what the panel asks for that no lane represents. "
            "4: directional axes below 'trust' and how many agreeing claims would lift them. 5: where the "
            "envelope misses most (value of information). 6: what each measurement would let you predict."
        ),
        "unidentified_coordinates": coords,
        "not_evaluable_rows": nev,
        "refused_targets": refused,
        "thin_axes": axes,
        "largest_misses": misses,
        "what_you_could_predict": unlocked,
        "summary": {
            "unidentified_coordinates": len(coords),
            "not_evaluable_rows": len(nev),
            "refused_row_count": sum(g["row_count"] for g in refused),
            "axes_below_trust": len(axes),
            "trust_rule": f"rate >= {TRUST_MIN_RATE} on >= {MIN_EVALUABLE_FOR_TRUST} claims and Wilson lower bound > {COIN}",
        },
    }
    return payload


def render_markdown(payload: Mapping[str, Any]) -> str:
    s = payload["summary"]
    lines = [
        "# Data wishlist — what to measure next, and what it would unlock",
        "",
        f"> Generated by `src/data_wishlist.py` from the tracked scorecard, slice profile, Laplace covariance, "
        f"directional scorecard and value-of-information ranking ({payload['provenance']['generated_on']}). "
        "Regenerate with `./scripts/docker_maillard.sh wishlist`; read with `python scripts/maillard.py wishlist`.",
        "",
        payload["how_to_read"],
        "",
        f"**At a glance:** {s['unidentified_coordinates']} fitted coordinates unidentified by primary evidence · "
        f"{s['not_evaluable_rows']} panel rows answered but not evaluable · {s['refused_row_count']} panel rows refused "
        f"(targets no lane represents) · {s['axes_below_trust']} directional axes below 'trust' ({s['trust_rule']}).",
        "",
        "## 1. Fitted coordinates the data do not identify",
        "",
        "| coordinate | profile | step | one measurement that would identify it | unlocks |",
        "|---|---|---|---|---|",
    ]
    for c in payload["unidentified_coordinates"]:
        step = ((f"{c['reaction']}" if "steps sharing" in str(c["reaction"]) else f"`{c['reaction']}`: {' + '.join(c['reactants'])} → {' + '.join(c['products'])}") if c["reaction"]
                else ("shared barrier of a decay family" if c["key"].startswith("Ea_decay_") else "(no single step)"))
        unlocks = ", ".join(c["unlocks_observables"]) or "—"
        story = f" {c['story']}" if c["story"] else ""
        lines.append(f"| `{c['key']}` | {c['profile_verdict']} | {step} | {c['measurement']}.{story} | {unlocks} |")
    lines += ["", "## 2. Panel rows answered but not evaluable", ""]
    if payload["not_evaluable_rows"]:
        lines.append("The engine returns a number but declares the route unidentified; the rows leave the absolute count until section 1's `k_glc_ha` is measured:")
        lines.append("")
        for r in payload["not_evaluable_rows"]:
            lines.append(f"- `{r['benchmark_id']}` / {r['compound']} ({r['panel']})")
    else:
        lines.append("None.")
    lines += ["", "## 3. What the panel asks for that no lane represents", "", "| what | rows | the engine's reason |", "|---|---|---|"]
    for g in payload["refused_targets"]:
        lines.append(f"| {g['what']} | {g['row_count']} | {g['reason'][:260]} |")
    lines += ["", "## 4. Directional axes below 'trust'", "", "| axis | agree / evaluable | not evaluable | verdict | Wilson lower | agreeing claims to reach trust |", "|---|---|---|---|---|---|"]
    for a in payload["thin_axes"]:
        need = "blocked: " + a["structural_block"] if a["structural_block"] else (str(a["additional_agreeing_claims_to_trust"]) if a["additional_agreeing_claims_to_trust"] is not None else "—")
        lines.append(f"| {a['axis']} | {a['agree']} / {a['evaluable']} | {a['not_evaluable']} | {a['verdict']} | {a['wilson_lower_bound'] if a['wilson_lower_bound'] is not None else '—'} | {need} |")
    lines += ["", "## 5. Where the envelope misses most (value of information)", ""]
    if payload["largest_misses"]:
        lines += ["| # | benchmark | compound | VoI | miss (log10) | template |", "|---|---|---|---|---|---|"]
        for i, m in enumerate(payload["largest_misses"], 1):
            lines.append(f"| {i} | `{m['benchmark_id']}` | {m['compound']} | {m['voi_score']:.2f} | {m['envelope_miss_log10']:.2f} | {m['template']} |")
    else:
        lines.append("The value-of-information ranking is not generated (`./scripts/docker_maillard.sh experiment-value-ranking`).")
    lines += ["", "## 6. What you could predict if you had it", ""]
    for u in payload["what_you_could_predict"]:
        lines.append(f"- **If** `{u['if_measured']}` were measured: {u['you_could_predict']}.")
    lines += ["", "Nothing above is a forecast of accuracy: 'unlocks' means the observable sits downstream of the step in the network, so a measured rate would replace a band artefact with a fitted value. Whether the fitted value lands within 3x of a measurement is what the next pre-registered wave finds out (`scripts/generators/WAVES.md`)."]
    return "\n".join(lines) + "\n"


__all__ = ["OUTPUT_JSON", "build", "render_markdown", "sulfur_reactions", "downstream_observables", "upstream_charges"]
