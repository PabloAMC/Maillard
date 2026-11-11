"""Ingredients in % w/w mapped to the engine's precursors (STATUS_AND_PLAN sec. 0 step 3).

A spec may carry ``ingredients: {name: percent}``; :func:`expand` turns them into precursor millimolar charges (added
to any the spec states itself), a protein matrix and its grams per litre, and a report of what was mapped, from
which source, with which band, and what was NOT mapped. The table is ``data/species/ingredients.yml``; every value in
it names its dossier. The pot is taken at 1 kg/L (1 % w/w = 10 g/L). Nothing here moves a shipped constant.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Tuple

from src import data_access, data_paths

G_PER_L_PER_PERCENT = 10.0

_TABLE: Dict[str, Any] = {}


class IngredientError(ValueError):
    pass


def table() -> Mapping[str, Any]:
    if not _TABLE:
        _TABLE.update(data_access.load_yaml(data_paths.INGREDIENTS)["ingredients"])
    return _TABLE


def _key(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_").replace("-", "_")


def expand(spec: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """(resolved spec, report). The resolved spec has ``precursors`` (the spec's own plus the ingredients'),
    and, for a protein isolate, ``matrix`` and ``protein_g_per_l``. A spec without ``ingredients`` is returned
    unchanged with an empty report."""
    ingredients = spec.get("ingredients") or {}
    if not ingredients:
        return dict(spec), {}
    tab = table()
    out = dict(spec)
    precursors: Dict[str, float] = {str(k): float(v) for k, v in (spec.get("precursors") or {}).items()}
    report: Dict[str, Any] = {"basis": "pot at 1 kg/L: 1 % w/w = 10 g/L", "ingredients": {}}
    matrices = []
    for raw_name, pct in ingredients.items():
        key = _key(raw_name)
        if key not in tab:
            raise IngredientError(
                f"ingredient {raw_name!r} is not in data/species/ingredients.yml (known: {', '.join(sorted(tab))}); "
                "give its precursors in mM instead, or add it with a sourced composition")
        entry = tab[key]
        g_per_l = float(pct) * G_PER_L_PER_PERCENT
        rec: Dict[str, Any] = {"percent_w_w": float(pct), "g_per_l": g_per_l, "mapped": {},
                               "not_mapped": dict(entry.get("not_mapped") or {})}
        for name, row in (entry.get("precursors") or {}).items():
            mm = g_per_l * float(row["central"])
            lo, hi = (g_per_l * float(b) for b in row["band"])
            precursors[name] = precursors.get(name, 0.0) + mm
            rec["mapped"][name] = {"mM": mm, "band_mM": [lo, hi], "source": row["source"],
                                   **({"caveat": row["caveat"]} if row.get("caveat") else {})}
        if entry.get("protein_fraction"):
            pf = entry["protein_fraction"]
            protein = g_per_l * float(pf["central"])
            rec["protein_g_per_l"] = protein
            rec["protein_band_g_per_l"] = [g_per_l * float(b) for b in pf["band"]]
            rec["protein_source"] = pf["source"]
            matrices.append((entry["matrix"], protein))
        report["ingredients"][key] = rec
    if len(matrices) > 1:
        raise IngredientError("two protein isolates in one pot: the engine charges one matrix; give protein_sites "
                              "explicitly instead")
    if matrices:
        if spec.get("matrix") not in (None, "water", matrices[0][0]):
            raise IngredientError(f"the spec's matrix {spec.get('matrix')!r} conflicts with the isolate's "
                                  f"{matrices[0][0]!r}")
        out["matrix"], out["protein_g_per_l"] = matrices[0][0], matrices[0][1] + float(spec.get("protein_g_per_l") or 0)
    if not precursors:
        raise IngredientError("the ingredients map to no precursor (a protein isolate alone carries only its "
                              "sites); add a sugar or a sulfur source in mM")
    out["precursors"] = precursors
    out.pop("ingredients", None)
    return out, report
