#!/usr/bin/env python3
"""Build ``data/keys/compounds.yml`` -- the compound key registry.

Why this exists (cleaning branch, Phase 3, 2026-09-01)
-------------------------------------------------------
Every data file in this repository names compounds by free-text display name, and no
two files agree on the spelling: 148 distinct strings for roughly 85 real compounds
(``Hexanal``/``hexanal``, ``2-Methyl-3-furanthiol (MFT)``/``2-methyl-3-furanthiol``/
``2methyl3furanthiol``/``mft``, ``Furan-2-aldehyde (FA)`` for furfural ...). Five
independent alias tables in ``src/`` papered over that at runtime, one of them ending
in a ``difflib`` similarity match at ratio 0.75.

This script produces ONE registry, keyed on identity rather than spelling:

* ``id``            -- stable snake-case key used in code and, eventually, in data;
* ``inchikey``      -- computed with RDKit from ``smiles`` (molecules only);
* ``aliases``       -- every spelling seen in the repository, gathered automatically
                       from the species YAML, the matrix decision panel, the benchmark
                       payloads, the flavor/safety payloads and the Henry table, and
                       attached to an entry by normalised-name match;
* ``kind``          -- ``molecule`` | ``compound_class`` | ``marker_class`` | ``process_lever``;
* ``class``         -- the decision-panel target class where one exists.

The SEED table below is the hand-curated part. Species-YAML compounds bring their own
SMILES/InChI/CAS; entries added here by hand carry
``identity_source: "structure entered by hand 2026-09-01 (standard chemical identity)"``
so nobody mistakes them for measured data. Any repository spelling that resolves to no
seed makes the script fail loudly -- add a seed, do not let a name float free.

Run inside Docker (needs RDKit):
    python scripts/generators/build_compound_registry.py            # write data/keys/compounds.yml
    python scripts/generators/build_compound_registry.py --check    # exit 1 if the file is stale
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import yaml  # noqa: E402

from src import data_access, data_paths  # noqa: E402
from src.text_utils import normalize_compound_name as norm  # noqa: E402

HAND = "structure entered by hand 2026-09-01 (standard chemical identity, not a measurement)"

# --------------------------------------------------------------------------- seeds
# (id, display_name, kind, smiles, extra) -- aliases beyond the automatic ones go in extra.
SEEDS: List[Dict[str, Any]] = [
    # --- molecules that the species YAML does not carry -------------------------
    dict(id="2_3_butanedione", display_name="2,3-Butanedione", kind="molecule", smiles="CC(=O)C(C)=O", aliases=["diacetyl"]),
    dict(id="2_acetylfuran", display_name="2-Acetylfuran", kind="molecule", smiles="CC(=O)c1ccco1"),
    dict(id="2_methyltetrahydrofuran_3_one", display_name="2-Methyltetrahydrofuran-3-one", kind="molecule", smiles="CC1OCCC1=O"),
    dict(id="3_isobutyl_2_methoxypyrazine", display_name="3-Isobutyl-2-methoxypyrazine (IBMP)", kind="molecule", smiles="COc1nccnc1CC(C)C", aliases=["ibmp"], class_="methoxypyrazines"),
    dict(id="4_vinylguaiacol", display_name="4-Vinylguaiacol", kind="molecule", smiles="C=Cc1ccc(O)c(OC)c1"),
    dict(id="acetaldehyde", display_name="Acetaldehyde", kind="molecule", smiles="CC=O"),
    dict(id="acetoin", display_name="Acetoin", kind="molecule", smiles="CC(O)C(C)=O"),
    dict(id="acetone", display_name="Acetone", kind="molecule", smiles="CC(C)=O"),
    dict(id="acrolein", display_name="Acrolein", kind="molecule", smiles="C=CC=O"),
    dict(id="benzaldehyde", display_name="Benzaldehyde", kind="molecule", smiles="O=Cc1ccccc1"),
    dict(id="beta_ionone", display_name="beta-Ionone", kind="molecule", smiles="CC(=O)/C=C/C1=C(C)CCCC1(C)C"),
    dict(id="chlorogenic_acid", display_name="Chlorogenic acid", kind="molecule", smiles="O=C(OC1CC(O)(C(=O)O)CC(O)C1O)/C=C/c1ccc(O)c(O)c1", identity_note="ring stereocentres omitted; identity at constitution level only"),
    dict(id="e_2_octenal", display_name="(E)-2-Octenal", kind="molecule", smiles="CCCCC/C=C/C=O"),
    dict(id="e_e_2_4_heptadienal", display_name="(E,E)-2,4-Heptadienal", kind="molecule", smiles="CC/C=C/C=C/C=O"),
    dict(id="furan", display_name="Furan", kind="molecule", smiles="c1ccoc1"),
    dict(id="furosine", display_name="Furosine", kind="molecule", smiles="NC(CCCCNCC(=O)c1ccco1)C(=O)O", identity_note="alpha-carbon stereocentre omitted"),
    dict(id="heptanal", display_name="Heptanal", kind="molecule", smiles="CCCCCCC=O"),
    dict(id="methylpyrazine", display_name="2-Methylpyrazine", kind="molecule", smiles="Cc1cnccn1", class_="pyrazines"),
    dict(id="tetramethylpyrazine", display_name="Tetramethylpyrazine", kind="molecule", smiles="Cc1nc(C)c(C)nc1C", class_="pyrazines"),
    dict(id="2_6_dimethylpyrazine", display_name="2,6-Dimethylpyrazine", kind="molecule", smiles="Cc1cncc(C)n1", class_="pyrazines"),
    dict(id="2_ethylpyrazine", display_name="2-Ethylpyrazine", kind="molecule", smiles="CCc1cnccn1", class_="pyrazines"),
    dict(id="trimethylpyrazine", display_name="Trimethylpyrazine", kind="molecule", smiles="Cc1cnc(C)c(C)n1", class_="pyrazines"),
    dict(id="mercapto_2_propanone", display_name="Mercapto-2-propanone", kind="molecule", smiles="CC(=O)CS", aliases=["1-mercapto-2-propanone", "mercaptoacetone"]),
    dict(id="methanethiol", display_name="Methanethiol", kind="molecule", smiles="CS"),
    dict(id="norfuraneol", display_name="Norfuraneol (NF)", kind="molecule", smiles="CC1=C(O)C(=O)CO1", aliases=["4-hydroxy-5-methyl-3(2H)-furanone", "nf"]),
    dict(id="imp", display_name="Inosine 5'-monophosphate (IMP)", kind="molecule", smiles=None, cas="131-99-7", aliases=["imp", "inosine monophosphate", "inosinate"], class_="umami_support_markers"),
    dict(id="gmp", display_name="Guanosine 5'-monophosphate (GMP)", kind="molecule", smiles=None, cas="85-32-5", aliases=["gmp", "guanosine monophosphate", "guanylate"], class_="umami_support_markers"),
    # --- classes and markers -----------------------------------------------------
    dict(id="pyrazines", display_name="Pyrazines", kind="compound_class", aliases=["pyrazine", "pyrazine_family", "dimethylpyrazine"], members_of_class="pyrazines", identity_note="the bare word 'pyrazine' in benchmark and payload rows means the family, as the retired BENCHMARK_NAME_ALIASES table asserted"),
    dict(id="methoxypyrazines", display_name="Methoxypyrazines", kind="compound_class", aliases=["methoxypyrazines"], members_of_class="methoxypyrazines"),
    dict(id="cml_cel", display_name="CML + CEL (advanced glycation end-products)", kind="marker_class", aliases=["CML/CEL", "CML_CEL", "cml_cel"], members=["cml", "cel"]),
    dict(id="cml_cel_acrylamide", display_name="CML + CEL + acrylamide (damage marker set)", kind="marker_class", aliases=["cml_cel_acrylamide"], members=["cml", "cel", "acrylamide"]),
    dict(id="furosine_cml", display_name="Furosine + CML (early/advanced glycation pair)", kind="marker_class", aliases=["furosine_cml"], members=["furosine", "cml"]),
    dict(id="browning_intermediates", display_name="Browning intermediates (unspecified)", kind="marker_class", aliases=["browning_intermediates"]),
    dict(id="reactive_lysine", display_name="Reactive (available) lysine", kind="marker_class", aliases=["reactive_lysine"]),
    dict(id="heterocyclic_amines", display_name="Heterocyclic amines (HCAs)", kind="marker_class", aliases=["HCAs", "hcas"], members=["phip", "meiqx"]),
    # --- process levers that the decision panel scores like compounds ------------
    dict(id="thiamine_availability", display_name="Thiamine availability", kind="process_lever", class_="pretreatment_state_markers"),
    dict(id="nucleotide_enrichment", display_name="Nucleotide enrichment", kind="process_lever", class_="pretreatment_state_markers"),
    dict(id="free_amino_acid_enrichment", display_name="Free amino acid enrichment", kind="process_lever", class_="pretreatment_state_markers"),
    dict(id="pretreatment_ph_shift", display_name="Pretreatment pH shift", kind="process_lever", class_="pretreatment_state_markers"),
    dict(id="caramelization_severity", display_name="Caramelization severity", kind="process_lever", class_="severity_markers"),
]

# species-YAML display names -> registry id (the YAML stays the source of SMILES/CAS/InChI)
SPECIES_IDS = {
    "2-Methyl-3-furanthiol (MFT)": "2_methyl_3_furanthiol",
    "2-Furfurylthiol (FFT)": "2_furfurylthiol",
    "Methional (3-(methylthio)propanal)": "methional",
    "2-Methylbutanal": "2_methylbutanal",
    "3-Methylbutanal": "3_methylbutanal",
    "2-Methylpropanal (isobutyraldehyde)": "2_methylpropanal",
    "Phenylacetaldehyde": "phenylacetaldehyde",
    "2,3-Dimethylpyrazine": "2_3_dimethylpyrazine",
    "2,5-Dimethylpyrazine": "2_5_dimethylpyrazine",
    "2-Ethyl-3,5-dimethylpyrazine": "2_ethyl_3_5_dimethylpyrazine",
    "Dimethyl disulfide": "dimethyl_disulfide",
    "Bis(2-methyl-3-furyl) disulfide": "bis_2_methyl_3_furyl_disulfide",
    "Dimethyl trisulfide": "dimethyl_trisulfide",
    "2-Methylthiophene": "2_methylthiophene",
    "4,5-Dihydro-2-methylthiazole": "4_5_dihydro_2_methylthiazole",
    "2-pentyl-4-methylthiazole": "2_pentyl_4_methylthiazole",
    "2-hexyl-4-methylthiazole": "2_hexyl_4_methylthiazole",
    "Hydrogen Sulfide": "hydrogen_sulfide",
    "HEMF": "hemf",
    "DMHF": "hdmf",
    "Hexanal": "hexanal",
    "Nonanal": "nonanal",
    "1-Octen-3-ol": "1_octen_3_ol",
    "2-Pentylfuran": "2_pentylfuran",
    "1-Hexanol": "1_hexanol",
    "Furfural": "furfural",
    "Nε-(Carboxymethyl)lysine (CML)": "cml",
    "Nε-(Carboxyethyl)lysine (CEL)": "cel",
    "5-Hydroxymethylfurfural (HMF)": "hmf",
    "Lysinoalanine (LAL)": "lysinoalanine",
    "Lanthionine (LAN)": "lanthionine",
    "Acrylamide": "acrylamide",
    "2-Amino-1-methyl-6-phenylimidazo[4,5-b]pyridine (PhIP)": "phip",
    "2-Amino-3-methylimidazo[4,5-f]quinoxaline (MeIQx)": "meiqx",
}
# Extra spellings for species compounds that the automatic gather cannot derive.
EXTRA_ALIASES = {
    "2_methyl_3_furanthiol": ["2-methyl-3-furylthiol", "2-methylfuran-3-thiol", "2methylfuran3thiol", "2methyl3furylthiol"],
    "2_furfurylthiol": ["2-furylmethanethiol", "2furylmethanethiol", "furfuryl mercaptan"],
    "bis_2_methyl_3_furyl_disulfide": ["2-methyl-3-furyl 2-methyl-3-furyl disulfide", "mft disulfide", "bis(2-methyl-3-furyl) disulphide"],
    "methional": ["3-(methylthio)propanal", "3methylthiopropanal"],
    "3_methylbutanal": ["isovaleraldehyde"],
    "2_methylbutanal": ["2-methylbutyraldehyde"],
    "acrylamide": ["2-propenamide", "2propenamide"],
    "furfural": ["Furan-2-aldehyde (FA)", "furan-2-aldehyde", "furan-2-carbaldehyde", "fa"],
    "cml": ["carboxymethyl lysine", "nε-(carboxymethyl)lysine (cml)"],
    "cel": ["carboxyethyl lysine", "nε-(carboxyethyl)lysine (cel)"],
    "hdmf": ["furaneol", "hdmf", "dmhf", "4-hydroxy-2,5-dimethyl-3(2H)-furanone"],
    "hydrogen_sulfide": ["h2s"],
}
PANEL_CLASS_OVERRIDES = {"hdmf": "furans_furanones", "2_methyltetrahydrofuran_3_one": "furans_furanones"}


def _gather_spellings() -> Dict[str, set]:
    """Every compound spelling in the repository, with where it was seen."""
    seen: Dict[str, set] = {}

    def add(name: Any, where: str) -> None:
        if isinstance(name, str) and name.strip():
            seen.setdefault(name.strip(), set()).add(where)

    for f in glob.glob(str(data_paths.BENCHMARKS_DIR / "**" / "*.json"), recursive=True):
        d = json.load(open(f, encoding="utf-8"))
        for key in ("measured_volatiles", "reference_volatiles", "holdout_targets"):
            for name in (d.get(key) or {}):
                add(name, f"benchmarks:{key}")
    for p in (data_paths.DESIRABLE_TARGETS, data_paths.OFF_FLAVOUR_TARGETS, data_paths.TOXIC_MARKERS):
        for c in (data_access.load_yaml(p) or {}).get("compounds", []):
            add(c.get("name"), f"species:{p.name}")
    for c in (data_access.load_yaml(data_paths.HENRY_CONSTANTS) or {}).get("constants", []):
        add(c.get("name"), "henry_constants")
    for k, v in data_access.load_json(data_paths.FLAVOR_REFERENCE_PAYLOADS).items():
        if isinstance(v, list):
            for r in v:
                if isinstance(r, dict):
                    add(r.get("compound"), f"flavor:{k}")
    for r in data_access.load_json(data_paths.SAFETY_REFERENCE_PAYLOADS).get("entries", []):
        add(r.get("analyte"), "safety:analyte")
    for k, v in data_access.load_json(data_paths.MATRIX_DECISION_PANEL).get("entries", {}).items():
        add(k, "panel:key")
        add(v.get("display_name"), "panel:display_name")
        for a in v.get("aliases", []):
            add(a, "panel:alias")
    return seen


def _species_entries() -> List[Dict[str, Any]]:
    out = []
    for p, default_class in (
        (data_paths.DESIRABLE_TARGETS, None),
        (data_paths.OFF_FLAVOUR_TARGETS, None),
        (data_paths.TOXIC_MARKERS, "safety_markers"),
    ):
        for c in (data_access.load_yaml(p) or {}).get("compounds", []):
            name = c["name"]
            cid = SPECIES_IDS.get(name)
            if cid is None:
                raise SystemExit(f"species compound {name!r} in {p.name} has no id in SPECIES_IDS")
            aliases = {name, re.sub(r"\s*\([^()]*\)$", "", name).strip()}
            aliases.update(t.strip() for t in re.findall(r"\(([^()]*)\)$", name))
            aliases.update(EXTRA_ALIASES.get(cid, []))
            out.append(
                dict(
                    id=cid,
                    display_name=name,
                    kind="molecule",
                    smiles=c.get("smiles"),
                    inchi=c.get("inchi"),
                    cas=c.get("cas"),
                    class_=default_class,
                    aliases=sorted(aliases),
                    identity_source=data_paths.rel(p),
                )
            )
    return out


def _panel_classes() -> Dict[str, str]:
    panel = data_access.load_json(data_paths.MATRIX_DECISION_PANEL).get("entries", {})
    return {norm(k): v.get("target_class") for k, v in panel.items()} | {
        norm(v.get("display_name", k)): v.get("target_class") for k, v in panel.items()
    }


def build() -> Dict[str, Any]:
    try:
        from rdkit import Chem  # type: ignore
        from rdkit import RDLogger  # type: ignore

        RDLogger.DisableLog("rdApp.*")
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("RDKit is required (run inside the Docker container)") from exc

    entries: Dict[str, Dict[str, Any]] = {}
    for seed in _species_entries() + [dict(s) for s in SEEDS]:
        cid = seed["id"]
        if cid in entries:
            raise SystemExit(f"duplicate id {cid}")
        entries[cid] = seed

    panel_classes = _panel_classes()
    alias_index: Dict[str, str] = {}
    for cid, e in entries.items():
        names = set(e.get("aliases") or [])
        names.add(e["display_name"])
        names.add(re.sub(r"\s*\([^()]*\)$", "", e["display_name"]).strip())
        names.update(t.strip() for t in re.findall(r"\(([^()]*)\)$", e["display_name"]))
        names.add(cid)
        e["_norms"] = {norm(n) for n in names if norm(n)}
        for n in e["_norms"]:
            if n in alias_index and alias_index[n] != cid:
                raise SystemExit(f"alias {n!r} claimed by both {alias_index[n]} and {cid}")
            alias_index[n] = cid

    # attach every repository spelling; fail on anything unresolved
    spellings = _gather_spellings()
    unresolved = []
    for spelling, where in sorted(spellings.items()):
        cid = alias_index.get(norm(spelling))
        if cid is None:
            unresolved.append((spelling, sorted(where)))
            continue
        entries[cid].setdefault("aliases", [])
        if spelling not in entries[cid]["aliases"]:
            entries[cid]["aliases"].append(spelling)
        entries[cid].setdefault("_seen_in", set()).update(where)
    if unresolved:
        for s, w in unresolved:
            print(f"UNRESOLVED: {s!r} seen in {w}", file=sys.stderr)
        raise SystemExit(f"{len(unresolved)} repository spellings resolve to no registry entry; add seeds")

    # inchikeys + classes
    out = []
    for cid, e in entries.items():
        rec: Dict[str, Any] = {"id": cid, "display_name": e["display_name"], "kind": e["kind"]}
        smiles = e.get("smiles")
        if e["kind"] == "molecule":
            rec["smiles"] = smiles
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise SystemExit(f"{cid}: RDKit cannot parse SMILES {smiles!r}")
                rec["inchikey"] = Chem.MolToInchiKey(mol)
                rec["canonical_smiles"] = Chem.MolToSmiles(mol)
            else:
                rec["inchikey"] = None
            if e.get("inchi"):
                rec["inchi"] = e["inchi"]
            if e.get("cas"):
                rec["cas"] = e["cas"]
        cls = e.get("class_") or PANEL_CLASS_OVERRIDES.get(cid)
        if cls is None:
            for n in e["_norms"]:
                if n in panel_classes and panel_classes[n]:
                    cls = panel_classes[n]
                    break
        if cls:
            rec["class"] = cls
        if e.get("members"):
            rec["members"] = e["members"]
        if e.get("members_of_class"):
            rec["members_of_class"] = e["members_of_class"]
        rec["aliases"] = sorted({a for a in e.get("aliases", []) if a != e["display_name"]}, key=lambda a: (a.lower(), a))
        rec["seen_in"] = sorted(e.get("_seen_in", set()))
        rec["identity_source"] = e.get("identity_source") or (HAND if e["kind"] == "molecule" else "decision panel / safety registry vocabulary")
        if e.get("identity_note"):
            rec["identity_note"] = e["identity_note"]
        out.append(rec)

    # InChIKey uniqueness across molecules
    by_key: Dict[str, str] = {}
    for rec in out:
        key = rec.get("inchikey")
        if key:
            if key in by_key:
                raise SystemExit(f"{rec['id']} and {by_key[key]} share InChIKey {key}")
            by_key[key] = rec["id"]

    out.sort(key=lambda r: (r["kind"], r["id"]))
    for rec in out:
        if rec["kind"] == "molecule" and rec.get("members"):
            del rec["members"]
    return {
        "schema_version": 1,
        "generated_by": "scripts/generators/build_compound_registry.py",
        "note": (
            "One row per compound identity. Molecules are keyed by InChIKey (computed with "
            "RDKit from `smiles`); classes, marker sets and process levers are keyed by id. "
            "`aliases` is every spelling that appears anywhere in the repository; every "
            "compound string in data/ must resolve here (tests/unit/test_compound_keys.py). "
            "Do not hand-edit: change the seeds in the generator and re-run it."
        ),
        "compounds": out,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true", help="exit 1 if data/keys/compounds.yml is stale")
    args = parser.parse_args(argv)
    payload = build()
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=True, width=110)
    target = data_paths.COMPOUND_REGISTRY
    if args.check:
        current = target.read_text(encoding="utf-8") if target.exists() else ""
        if current != text:
            print(f"{data_paths.rel(target)} is stale; re-run without --check", file=sys.stderr)
            return 1
        print(f"{data_paths.rel(target)} is current ({len(payload['compounds'])} entries)")
        return 0
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    kinds = {}
    for rec in payload["compounds"]:
        kinds[rec["kind"]] = kinds.get(rec["kind"], 0) + 1
    print(f"wrote {data_paths.rel(target)}: {len(payload['compounds'])} entries {kinds}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
