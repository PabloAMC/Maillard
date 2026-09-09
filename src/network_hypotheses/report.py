"""The hypothesis artifact: what the rules propose from each lane's reference charge, placed against the engine."""
from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Mapping

from src import data_paths
from src.network_hypotheses import classify, rules as R, structures as S
from src.network_hypotheses.enumerate import enumerate_from

OUTPUT_JSON = data_paths.VALIDATION_DIR / "network_hypotheses.json"

#: The charges the layer is run on: each lane's reference pot, plus a probe of the thiols in the
#: company of the intermediates the pot holds, so every thiol rule gets a chance to fire.
CHARGES: Dict[str, Dict[str, Any]] = {
    "sugar_glycine": {"species": ["Glc", "Gly"], "lane": "trunk", "why": "the trunk's fitted pot (Martins 2005)"},
    "pentose_cysteine": {"species": ["PENT", "Cys"], "lane": "sulfur", "why": "the sulfur lane's reference pot (Hofmann 1998)"},
    "pentose_cysteine_thiamine": {"species": ["PENT", "Cys", "THI"], "lane": "sulfur", "why": "the thiamine route (Hofmann 1998 Table 8)"},
    "asparagine_glucose": {"species": ["Asn", "Glc"], "lane": "acrylamide", "why": "the acrylamide lane's pot (De Vleeschouwer)"},
    "linoleate_hydroperoxides": {"species": ["LOOH_13_ct", "LOOH_9_ct", "LOOH_10"], "lane": "lipid",
                                 "why": "the lipid lane's hydroperoxide pool (Frankel 1989) and the 10-hydroperoxide it lacks (Miyazaki 2023)"},
    "lipid_maillard_cross": {
        "species": ["DECADIENAL", "HEXANAL", "AMMONIA", "H2S", "Cys"],
        "lane": "lipid",
        "why": "the fatty aldehydes an isolate's lipid makes, with the ammonia and hydrogen sulfide the Maillard side supplies: the cross products no lane names (Zamora 2020, Zhou 2000, Du 2023)",
        "depth": 1,
    },
    "lipid_maillard_thiazoles": {
        "species": ["HEXANAL", "ACETOL", "ACETOIN", "AMMONIA", "H2S", "DECADIENAL"],
        "lane": "lipid",
        "why": "Elmore 1997's pot: a lipid alkanal, a Maillard hydroxyketone, ammonia and hydrogen sulfide give the 2-alkyl-3-thiazolines and, oxidised, the registry's 2-alkyl-4-methylthiazoles; the dienal in the same charge shows where the H2S goes instead (Farmer 1990, Mottram 2002)",
        "depth": 3,
    },
    "oleate_hydroperoxides": {
        "species": ["OL_8_OOH", "OL_9_OOH", "OL_10_OOH", "OL_11_OOH"],
        "lane": "lipid",
        "why": "the four oleate hydroperoxides the lipid lane lumps as LOOH_OL with no edge (nonanal is a declared hold-out): Cao 2020's routes",
    },
    "strecker_to_pyrazines": {
        "species": ["MGO", "GO", "DA", "Ala", "Cys"],
        "lane": "trunk",
        "why": "the small dicarbonyls with an amino acid: the Strecker aldehydes, the aminoketones and the pyrazines they condense to (no lane has a pyrazine)",
    },
    "thiol_sink_probe": {
        "species": ["MFT", "FFT", "MESH", "Cys", "H2S", "PENT", "NF", "FUR", "HMF", "MGO", "GO", "DA", "DECADIENAL", "HEXANAL", "ACR"],
        "lane": "sulfur",
        "why": "the two thiols with every carbonyl and thiol partner the pots hold: what could remove them",
        "depth": 1,
    },
}


def build() -> Dict[str, Any]:
    table = S.load()
    rules = R.load()
    by_canon = classify.species_by_canonical(table)
    registry = S.registry_by_canonical()
    reactions = classify.engine_reactions()
    rule_by_id = {r.id: r for r in rules}
    label: Dict[str, str] = {}
    for key, st in table.items():
        if st.kind == "molecule":
            label[st.canonical] = f"{label[st.canonical]}/{key}" if st.canonical in label else key
    known = lambda canon: canon in by_canon or canon in registry
    charges_out = []
    for name, spec in CHARGES.items():
        steps, pool = enumerate_from(spec["species"], table, rules, depth=spec.get("depth", 2), may_expand=known)
        rows = []
        for step in steps:
            placement, matches = classify.place_step(step, rule_by_id[step.rule_id].status, by_canon, reactions)
            rows.append({
                "rule": step.rule_id,
                "depth": step.depth,
                "reactants": [label.get(c, c) for c in step.reactants],
                "products": [label.get(c, c) for c in step.products],
                "placement": placement,
                "engine_reactions": matches,
            })
        products = []
        for canon, depth in pool.items():
            if depth == 0:
                continue
            kind, ident = classify.place_product(canon, by_canon, registry)
            products.append({"smiles": canon, "first_depth": depth, "kind": kind, "id": ident})
        charges_out.append({
            "charge": name, "lane": spec["lane"], "why": spec["why"], "species": spec["species"],
            "depth": spec.get("depth", 2),
            "steps": rows,
            "products": sorted(products, key=lambda p: (p["first_depth"], p["kind"], p["smiles"])),
            "summary": {
                "steps": dict(Counter(r["placement"] for r in rows)),
                "products": dict(Counter(p["kind"] for p in products)),
            },
        })
    return {
        "artifact": "network_hypotheses",
        "what_this_is": (
            "What the literature's reaction rules (data/lit/reaction_rules.yml) propose from each lane's "
            "reference charge, placed against the engine's own reactions. Steps and products only: no rate, "
            "no concentration, and nothing here is read by the engine. Beyond the first step only products "
            "that are engine species or registry compounds react further, so the walk stays on the known map."
        ),
        "placements": dict(classify.PLACEMENTS),
        "rules": [{"id": r.id, "name": r.name, "class": r.klass, "status": r.status, "source": dict(r.source)} for r in rules],
        "charges": charges_out,
        "summary": {
            "rules": len(rules),
            "steps": sum(len(c["steps"]) for c in charges_out),
            "steps_by_placement": dict(Counter(r["placement"] for c in charges_out for r in c["steps"])),
            "products_by_kind": dict(Counter(p["kind"] for c in charges_out for p in c["products"])),
        },
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    out = ["# What the reaction rules propose, against what the engine models", "",
           f"*{payload['what_this_is']}*", ""]
    s = payload["summary"]
    out += [f"{s['rules']} rules; {s['steps']} proposed steps: " +
            ", ".join(f"{n} {k.replace('_', ' ')}" for k, n in sorted(s["steps_by_placement"].items())) +
            "; products: " + ", ".join(f"{n} {k}" for k, n in sorted(s["products_by_kind"].items())) + ".", ""]
    out += ["| placement | meaning |", "|---|---|"] + [f"| {k} | {v} |" for k, v in payload["placements"].items()] + [""]
    for c in payload["charges"]:
        out += [f"## {c['charge'].replace('_', ' ')} ({c['lane']} lane, depth {c['depth']})", "", f"*{c['why']}.* Charge: {', '.join(c['species'])}.", ""]
        out += ["| rule | reactants | products | placement | engine reaction |", "|---|---|---|---|---|"]
        for r in sorted(c["steps"], key=lambda r: ({"modelled": 0, "mechanism_known": 1, "proposed": 2}[r["placement"]], r["rule"])):
            out.append(f"| {r['rule']} | {' + '.join(r['reactants'])} | {' + '.join(r['products'])} | {r['placement'].replace('_', ' ')} | {', '.join(r['engine_reactions']) or ''} |")
        new = [p for p in c["products"] if p["kind"] != "species"]
        if new:
            out += ["", "Products that are not engine species: " +
                    "; ".join(f"`{p['smiles']}`" + (f" ({p['id']})" if p["id"] else "") for p in new) + "."]
        out.append("")
    out += ["## Rules", "", "| rule | what | status | source |", "|---|---|---|---|"]
    out += [f"| {r['id']} | {r['name']} | {r['status']} | {r['source']['dossier']}: {r['source']['anchor']} |" for r in payload["rules"]]
    return "\n".join(out) + "\n"
