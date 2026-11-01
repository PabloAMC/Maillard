#!/usr/bin/env python
"""
check_engine_family_map.py -- the engine side of the pathway-condition atlas, checked against the code (2026-11-01).

Asserts that docs/reference/pathway_atlas/engine_family_map.yml
  * names every Reaction key of the three mass-action lanes (network.TRUNK_REACTIONS, sulfur.FULL_REACTIONS,
    acrylamide.FULL_ACRYLAMIDE_REACTIONS) with exactly the lanes whose tuple carries it, and the lipid lane's analytic
    steps;
  * gives each a primary family id from families.yml (or a declared pseudo id) and only valid secondary ids;
  * says `shipped` truthfully: k != 0 at 100 C (pH 7 on the sulfur lane, the reference elsewhere) and not in
    sulfur.inactive_reactions() -- rate constants only, no integration;
  * lists, per family, exactly the reactions that name it, and an engine_gaps.no_reaction list that matches;
then prints the coverage per family. Exit status 1 on any failure.

    ./scripts/docker_maillard.sh run "python -m scripts.generators.check_engine_family_map"
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import yaml

from src.kinetic_core import acrylamide, engine, network, sulfur
from src.kinetic_core.parameters_lipid import COVALENT_SINK, DECADIENAL_BREAKAGE

ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "docs" / "reference" / "pathway_atlas"
MAP_PATH = ATLAS / "engine_family_map.yml"
FAMILIES_PATH = ATLAS / "families.yml"

LANE_TUPLES = {
    "trunk": network.TRUNK_REACTIONS,
    "sulfur": sulfur.FULL_REACTIONS,
    "acrylamide": acrylamide.FULL_ACRYLAMIDE_REACTIONS,
}
#: The lipid lane's analytic steps (lipid.integrate_lipid) and their shipped state.
LIPID_STEPS = {
    "lipid.looh_decomposition": True,
    "lipid.slate_partition": True,
    "lipid.pentylfuran_from_linoleate": True,
    "lipid.nonanal_from_oleate": True,
    "lipid.decadienal_breakage": bool(DECADIENAL_BREAKAGE.enabled),
    "lipid.covalent_aldehyde_lysine_sink": bool(COVALENT_SINK.enabled),
}
REQUIRED_REACTION_FIELDS = ("lanes", "family", "shipped", "ph", "T")
FAMILY_BINS = ("ph", "T", "o2", "amine", "buffer", "water")
T_100C = 373.15


def shipped_rate_constants() -> dict:
    """lane -> {reaction key: k at 100 C} in the shipped configuration (no integration)."""
    return {
        "trunk": network.rate_constants_at(engine.core_parameters("trunk"), T_100C),
        "sulfur": sulfur.sulfur_rate_constants_at(engine.core_parameters("sulfur"), T_100C, 7.0),
        "acrylamide": acrylamide.acrylamide_rate_constants_at(engine.core_parameters("acrylamide"), T_100C),
    }


def main() -> int:
    families_doc = yaml.safe_load(FAMILIES_PATH.read_text())
    family_ids = [f["id"] for f in families_doc["families"]]
    doc = yaml.safe_load(MAP_PATH.read_text())
    pseudo = list(doc.get("pseudo_ids") or [])
    valid = set(family_ids) | set(pseudo)
    reactions = doc["reactions"]
    errors: list[str] = []

    # 1. every Reaction key of the three lanes, with its true lane membership
    membership: dict[str, set] = defaultdict(set)
    for lane, tup in LANE_TUPLES.items():
        keys = [r.key for r in tup]
        if len(keys) != len(set(keys)):
            errors.append(f"{lane}: duplicate reaction keys in the tuple")
        for key in keys:
            membership[key].add(lane)
    for key in LIPID_STEPS:
        membership[key].add("lipid")
    for key, lanes in sorted(membership.items()):
        entry = reactions.get(key)
        if entry is None:
            errors.append(f"{key}: missing from the map (lanes {sorted(lanes)})")
            continue
        if set(entry.get("lanes") or ()) != lanes:
            errors.append(f"{key}: map lanes {entry.get('lanes')} != tuple membership {sorted(lanes)}")
    for key in reactions:
        if key not in membership:
            errors.append(f"{key}: in the map but in no lane tuple and not a lipid step")

    # 2. fields and ids
    for key, entry in reactions.items():
        for field in REQUIRED_REACTION_FIELDS:
            if field not in entry:
                errors.append(f"{key}: missing field {field!r}")
        if entry.get("family") not in valid:
            errors.append(f"{key}: primary {entry.get('family')!r} is not a family id or a declared pseudo id")
        for sec in entry.get("secondary") or ():
            if sec not in family_ids:
                errors.append(f"{key}: secondary {sec!r} is not a family id")
            if sec == entry.get("family"):
                errors.append(f"{key}: secondary repeats the primary")
        if entry.get("shipped") is False and not entry.get("off_reason"):
            errors.append(f"{key}: shipped: false without off_reason")

    # 3. `on` against the shipped rate constants and the sulfur wave switches
    k = shipped_rate_constants()
    inactive = sulfur.inactive_reactions()
    for key, entry in reactions.items():
        lanes = membership.get(key, set())
        if "lipid" in lanes:
            truth = LIPID_STEPS[key]
        else:
            truth = all(k[lane][key] != 0.0 for lane in lanes) and key not in inactive
        if bool(entry.get("shipped")) != truth:
            errors.append(f"{key}: map says shipped={entry.get('shipped')}, shipped configuration says {truth}")

    # 4. the per-family lists and the gap list
    derived_primary: dict[str, list] = defaultdict(list)
    derived_secondary: dict[str, list] = defaultdict(list)
    for key, entry in reactions.items():
        derived_primary[entry["family"]].append(key)
        for sec in entry.get("secondary") or ():
            derived_secondary[sec].append(key)
    fam_block = doc.get("families") or {}
    for fid in family_ids:
        block = fam_block.get(fid)
        if block is None:
            errors.append(f"families: {fid} missing")
            continue
        for field in FAMILY_BINS:
            if field not in block:
                errors.append(f"families.{fid}: missing bin {field!r}")
        if set(block.get("primary") or ()) != set(derived_primary.get(fid, ())):
            errors.append(f"families.{fid}.primary {sorted(block.get('primary') or ())} != "
                          f"reactions naming it {sorted(derived_primary.get(fid, ()))}")
        if set(block.get("secondary") or ()) != set(derived_secondary.get(fid, ())):
            errors.append(f"families.{fid}.secondary {sorted(block.get('secondary') or ())} != "
                          f"reactions naming it {sorted(derived_secondary.get(fid, ()))}")
    for fid in fam_block:
        if fid not in family_ids:
            errors.append(f"families: {fid} is not in families.yml")
    no_reaction = {f for f in family_ids if not derived_primary.get(f) and not derived_secondary.get(f)}
    gaps = doc.get("engine_gaps") or {}
    if set(gaps.get("no_reaction") or ()) != no_reaction:
        errors.append(f"engine_gaps.no_reaction {sorted(gaps.get('no_reaction') or ())} != {sorted(no_reaction)}")
    on_by_family = {
        f: [r for r in derived_primary.get(f, []) + derived_secondary.get(f, []) if reactions[r]["shipped"]]
        for f in family_ids
    }
    no_shipped = {f for f in family_ids if f not in no_reaction and not on_by_family[f]}
    if set(gaps.get("no_shipped_reaction") or ()) != no_shipped:
        errors.append(f"engine_gaps.no_shipped_reaction {sorted(gaps.get('no_shipped_reaction') or ())} "
                      f"!= {sorted(no_shipped)}")
    for kind, block in (gaps.get("condition_not_expressible") or {}).items():
        for fid in block.get("families") or ():
            if fid not in family_ids:
                errors.append(f"engine_gaps.condition_not_expressible.{kind}: {fid!r} is not a family id")

    # 5. coverage report
    print("ENGINE FAMILY MAP -- coverage per family (P = primary, S = secondary, on = shipped with k != 0)")
    print(f"{'family':34s} {'P':>3s} {'S':>3s} {'P on':>5s} {'P+S on':>7s}  lanes")
    for fid in family_ids:
        prim = derived_primary.get(fid, [])
        sec = derived_secondary.get(fid, [])
        lanes = sorted({lane for r in prim + sec for lane in membership[r]})
        p_on = sum(1 for r in prim if reactions[r]["shipped"])
        flag = "  ABSENT" if fid in no_reaction else ("  NONE SHIPPED" if fid in no_shipped else "")
        print(f"{fid:34s} {len(prim):3d} {len(sec):3d} {p_on:5d} {len(on_by_family[fid]):7d}  "
              f"{','.join(lanes)}{flag}")
    for pid in pseudo:
        prim = derived_primary.get(pid, [])
        print(f"{pid:34s} {len(prim):3d} {'-':>3s} {sum(1 for r in prim if reactions[r]['shipped']):5d}")
    print()
    for lane, tup in LANE_TUPLES.items():
        keys = [r.key for r in tup]
        n_on = sum(1 for kk in keys if reactions.get(kk, {}).get("shipped"))
        n_fam = sum(1 for kk in keys if reactions.get(kk, {}).get("family") in family_ids)
        print(f"lane {lane:10s}: {len(keys):3d} reactions, {n_on:3d} on, {n_fam:3d} with a family primary, "
              f"{len(keys) - n_fam:3d} pseudo")
    print(f"lane lipid     : {len(LIPID_STEPS):3d} analytic steps, {sum(LIPID_STEPS.values()):3d} on")
    covered = [f for f in family_ids if f not in no_reaction]
    shipped = [f for f in covered if f not in no_shipped]
    print(f"\nunique keys mapped: {len(reactions)}; families with any engine reaction: {len(covered)}/{len(family_ids)}; "
          f"with a shipped reaction: {len(shipped)}/{len(family_ids)}")
    print(f"absent: {sorted(no_reaction)}")
    print(f"none shipped: {sorted(no_shipped)}")

    if errors:
        print(f"\nFAIL: {len(errors)} problem(s)")
        for e in errors:
            print("  -", e)
        return 1
    print("\nOK: every Reaction key of the three lanes is mapped with a valid family id; lists and gaps consistent.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
