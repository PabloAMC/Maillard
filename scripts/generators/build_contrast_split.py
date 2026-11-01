#!/usr/bin/env python
"""
build_contrast_split -- the frozen, hashed DEVELOPMENT / TEST split of the contrast benchmark (v1, 2026-11-01).

Plan from 2026-11-01, rule (i) (docs/guides/STATUS_AND_PLAN.md): no circular scoring. Before anything is scored:
  - DEVELOPMENT = every atlas cell cited in docs/reference/pathway_atlas/ENGINE_VS_LITERATURE.md secs. 2 - 5 (by file +
    index, or by naming its study or laboratory) or in the two pre-registrations written from it
    (results/validation/kinetic_core_pyrazine_strecker_step_ph_prereg.md, kinetic_core_mft_route_barriers_prereg.md),
    PLUS every cell of any laboratory that has a cited cell;
  - TEST = the cells of the laboratories with no cited cell.
Laboratories are the canonical ids of docs/reference/contrast_benchmark/lab_ids.yml.

The candidate contrasts are
  - `paired_arms`: cells with evidence `within_study_factor` whose contrast / quantity carries a numeric effect (two or
    more arm values, a fold, or a ratio) and does not compare against another study;
  - `route_share`: `labelling` cells that print a numeric route share (a % or an n : m split).
Each gets a stable id `F25#9@meynier1995` (family, 0-based cell index, first dossier stem).

Outputs docs/reference/contrast_benchmark/split.yml. Its `split:` block (membership only) is hashed with sha256 over
json.dumps(split, sort_keys=True, separators=(",", ":"), ensure_ascii=True); that hash is frozen below. NO ENGINE IS RUN
and nothing is scored here: the per-error counts are availability (power) only.

    python scripts/generators/build_contrast_split.py            # check: rebuild in memory, assert hash + file
    python scripts/generators/build_contrast_split.py --write    # (re)write split.yml; refuses if membership moved

Changing the split needs a NEW, dated, reasoned version (a new VERSION / FROZEN_SHA256 here, the old split kept as
split_<version>.yml, a README entry), never an edit of this one.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "docs" / "reference" / "pathway_atlas"
BENCH = ROOT / "docs" / "reference" / "contrast_benchmark"
LAB_IDS = BENCH / "lab_ids.yml"
SPLIT = BENCH / "split.yml"
EVL = ATLAS / "ENGINE_VS_LITERATURE.md"
REGISTER = ROOT / "docs" / "reference" / "conflicts_register.yml"
VALIDATION = ROOT / "results" / "validation"
CITED_PREREGS = (
    VALIDATION / "kinetic_core_pyrazine_strecker_step_ph_prereg.md",
    VALIDATION / "kinetic_core_mft_route_barriers_prereg.md",
)

VERSION = "contrast_split_v1"
FROZEN_ON = "2026-11-01"
FROZEN_SHA256 = "022d0cab5a53ececb27a1fc6837ff275fbab177cf2161ae1f8a4ed0341f273a7"   # v1; a membership change needs a new version
# The atlas cells the split covers: per family, indices 0 .. n-1 at the freeze. Cells appended later (rule (i)'s
# "newly keyed rows", plan track P1) are outside v1: they are listed under `annotations.post_freeze_cells` and enter
# the benchmark only through a new version. (The atlas README: new cells are added at the END of their file.)
FROZEN_EXTENT = {
    "F01": 10, "F02": 17, "F03": 17, "F04": 27, "F05": 25, "F06": 18, "F07": 21, "F08": 10, "F09": 9, "F10": 17,
    "F11": 17, "F12": 10, "F13": 16, "F14": 14, "F15": 10, "F16": 18, "F17": 11, "F18": 17, "F19": 24, "F20": 14,
    "F21": 14, "F22": 17, "F23": 23, "F24": 13, "F25": 24, "F26": 26, "F27": 27, "F28": 27, "F29": 26, "F30": 21,
    "F31": 18, "F32": 4,
}

# ---------------------------------------------------------------------------------------------- fixed vocabularies
# Odorant class of each family (what the family's cells measure). Declared, not inferred.
ODORANT_CLASS = {
    **{f: "precursor_non_odorant" for f in ("F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10",
                                             "F11", "F12", "F14", "F15", "F31", "F32")},
    "F13": "furans_furanones",
    **{f: "strecker_aldehyde_acid" for f in ("F16", "F17", "F18")},
    **{f: "pyrazine_n_heterocycle" for f in ("F19", "F20", "F21")},
    **{f: "thiol_h2s_sulfur" for f in ("F22", "F23", "F24", "F25", "F26", "F27", "F28")},
    "F29": "acrylamide",
    "F30": "lipid_crossover",
}

# The factor a paired-arms contrast varies, read from the contrast's head (the text before its first ': ').
# Heuristic, multi-label; a factor counts only when its vocabulary sits within WINDOW characters of a comparator.
COMPARATOR = re.compile(r"\bvs\.?|->|\+/-|\+-|\bover\b|\bagainst\b|\bwithout\b|\bor\b|\bladder\b|\bseries\b|"
                        r"\b\d[\d.]*\s*/\s*~?\d|\b\d[\d.]*\s*[-\u2013]\s*\d|\bno [a-z]+ /|\bwith(?:out)?\b|"
                        r"\balone\b|\bonly\b")
WINDOW = 40
FACTOR_VOCAB = {
    "ph": r"\bpH\b",
    "temperature": r"temperature|\b\d+(?:\.\d+)?\s*(?:/|vs\.?|->|-|to|and)\s*\d+(?:\.\d+)?\s*C\b|\bT\b",
    "buffer": r"phosphate|buffer|citrate|malonate|bicarbonate|acetate|borate|unbuffered",
    "oxygen": r"\bair\b|argon|\bN2\b|\bO2\b|aerat|degass|deaerat|\binert\b|oxygen|anaerob|glovebox",
    "metals": r"chelat|DTPA|EDTA|DETAPAC|Cu\(II\)|\bCu\b|\bFe\b|metal-free",
    "water": r"a_w|\bWATER\b|\bwater\b|\bdry\b|moisture|glycerol|solvent-free|\bpowder|\bsolid\b",
    "lipid": r"lipid|\bfat\b|\boil\b|lecithin|linole|linolen|18:[23]|PUFA|dienal|hexanal|nonanal|nonenal|"
             r"fatty acid|phospholipid|triglycer|\bester\b",
    "sugar": r"\bSUGAR\b|sugar|glucose|ribose|xylose|fructose|pentose|hexose|lactose|sucrose|maltose|rhamnose|"
             r"arabinose|\bGlc\b|\bXyl\b|\bRib\b|\bFru\b",
    "amine": r"amino acid|glycine|lysine|cysteine|cystine|\bCys\b|alanine|\bAla\b|\bGly\b|\bLys\b|asparagine|"
             r"\bAsn\b|peptide|\bamine\b|protein|leucine|\bLeu\b|phenylalanine|\bPhe\b|glutathione|\bGSH\b|"
             r"methionine|\bMet\b|ammonia|NH3|thiamin|\bARP\b|Amadori",
    "oxidant": r"oxidant|persulfate|K2S2O8|peroxide",
    "time": r"\btime\b|time course|\bhold",
    "concentration": r"concentration|loading|\bratio\b|halving|dilut|charge|\bdose",
}

# The scope of each open engine error E1 - E20 (ENGINE_VS_LITERATURE sec. 3 and the binding sec. 5): which candidate
# contrasts could test it. `families` restricts the family; `factors` (any of) restricts the varied factor;
# `kinds` restricts the contrast kind (rule (ii): route-share changes, E6 / E11, are judged on LABELLING cells);
# `text` (regex on contrast + quantity) and `sugar_tag` narrow further. Declared 2026-11-01 from secs. 3 and 5; after
# the first build E16's text was narrowed from (norfuraneol | NF | HDMF | furaneol) to norfuraneol only, because E16 is
# norfuraneol's decay and furaneol (HDMF) is another compound (the change removed one TEST contrast, F13#13).
ERROR_SCOPES = {
    "E1":  {"title": "hexose MFT has no intact-C5 route", "families": ["F23", "F24", "F25"],
            "sugar_tag": "hexose", "text": r"MFT|FFT|C5|glucose"},
    "E2":  {"title": "hexose 1,2- vs 2,3-enolisation split frozen in pH", "families": ["F06", "F07"],
            "factors": ["ph"]},
    "E3":  {"title": "no Strecker-aldehyde lane", "families": ["F16"]},
    "E4":  {"title": "pyrazine pH law on the Strecker step", "families": ["F16", "F19", "F20"], "factors": ["ph"]},
    "E5":  {"title": "no buffer-species catalysis", "factors": ["buffer"]},
    "E6":  {"title": "MFT via free norfuraneol (route share)", "families": ["F23", "F26"], "kinds": ["route_share"],
            "text": r"\bNF\b|norfuraneol|mercapto-3-pentanone|DDP|dideoxypentosone"},
    "E7":  {"title": "lipid lane inert to Maillard odorants", "families": ["F30", "F23", "F25", "F26", "F28"],
            "factors": ["lipid"], "all_of_family": ["F30"]},
    "E8":  {"title": "FFT / furfural supply pH law too steep", "families": ["F25", "F13"], "factors": ["ph"],
            "text": r"FFT|furfural|furfurylthiol"},
    "E9":  {"title": "acrylamide: glucose-only, phosphate law, irreversible adducts", "families": ["F29"]},
    "E10": {"title": "TTCA formation off, release too fast and pH-flat", "families": ["F04"]},
    "E11": {"title": "one formation barrier for every sulfur route (route mix vs T)",
            "families": ["F23", "F24", "F25", "F26"], "kinds": ["route_share"]},
    "E12": {"title": "thiol losses: disulfide terminal, one decay, fixed O2", "families": ["F28"]},
    "E13": {"title": "CF-077 methylpyrazine cross-pair factor", "families": ["F19", "F20"],
            "text": r"methylpyrazine|\bMPZ\b|cross"},
    "E14": {"title": "aminoketone condensation instant; no pyrrole sink / extra MPZ route",
            "families": ["F19", "F20", "F21"]},
    "E15": {"title": "no Strecker acid / Amadori O2 route", "families": ["F17"]},
    "E16": {"title": "norfuraneol decay lacks base / O2 catalysis", "families": ["F13", "F15"],
            "factors": ["ph", "oxygen"], "text": r"norfuraneol|\bNF\b"},
    "E17": {"title": "pentose ARP decays return no amine", "families": ["F03"]},
    "E18": {"title": "3-mercapto-2-pentanone only from thiamine; thiamine pH", "families": ["F26", "F27"],
            "text": r"mercapto-2-pentanone|3-MP|thiamin"},
    "E19": {"title": "cysteine-Strecker H2S off", "families": ["F22"]},
    "E20": {"title": "no lipid Strecker; no dienal + NH3 -> pentylpyridine", "families": ["F18", "F21", "F30"],
            "text_unless_family": {"F18": None, "F21": r"pyridine", "F30": r"pyridine"}},
}

# Hand-reviewed factor tags for every TEST contrast (2026-11-01), replacing the heuristic for the held-out set, where the
# per-error counts are read. The generator fails if a TEST contrast lacks a review or a review names no TEST contrast.
# DEVELOPMENT contrasts keep the heuristic tags (their counts are context, not power). Reading notes as comments.
FACTOR_REVIEW = {
    "F02#1@potman1989": ["buffer"],                        # phosphate 0 -> 0.5 M at start pH 5.6
    "F04#26@mizutani1960": ["amine"],                      # free Cys vs aldehyde-bound (thiazolidine) Cys
    "F05#18@warmbier1976": ["water", "temperature"],       # a_w ladder; 35 -> 45 C
    "F05#19@warmbier1976": ["water"],
    "F05#20@pereyragonzales2010": ["water"],
    "F05#21@pereyragonzales2010": ["water"],               # a_w 0.33 vs 0.43 across the glass transition
    "F06#16@lee2022": ["amine"],                           # glucose + Leu vs glucose only
    "F08#9@lee2022": ["amine", "temperature"],
    "F10#14@huyghuesdespointes1996": [],                   # label position; one condition
    "F11#15@novotny2008": ["oxidant"],                     # NaOH vs + persulfate
    "F11#16@huyghuesdespointes1996": [],
    "F12#9@novotny2008": ["oxidant"],
    "F13#12@haleva1999": ["ph", "amine"],                  # pH 3 vs 5 in ONE buffer (citrate-phosphate): not buffer
    "F13#13@haleva1999": ["ph"],
    "F14#7@garrett1969": ["ph", "temperature", "sugar"],   # HCl (a_H+) ladder, 60 - 80 C, four pentoses
    "F14#11@haleva1999": ["ph", "amine"],
    "F14#12@haleva1999": ["ph", "amine"],
    "F15#3@huyghuesdespointes1996": [],
    "F16#14@huang2017": ["temperature", "sugar", "amine"],
    "F19#17@shibamoto1977": ["oxygen"],
    "F20#13@shibamoto1977": ["oxygen", "metals"],
    "F22#5@mizutani1960": ["amine"],                       # thiazolidine-bound vs free Cys
    "F22#6@obata1960": ["amine", "concentration"],         # cystine vs cysteine; dilution
    "F27#18@thomas2014": ["amine", "concentration"],       # thiamine dose vs cysteine dose
    "F28#3@kumazawa2003": ["ph"],                          # buffer species changes along the ladder: confounded
    "F28#18@stack2018": ["temperature"],
    "F29#11@quan2020": ["temperature", "amine"],
    "F30#17@meynier2004": ["lipid"],                       # hexanal vs (E)-2-hexenal with protein
    "F31#0@fang2009": ["water"],
    "F31#1@fang2009": ["water", "ph"],
    "F31#2@fang2010": [],
    "F31#7@cammerer1994": ["sugar"],
}

# ------------------------------------------------------------------------------------------------- small helpers
def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def canonical_hash(obj) -> str:
    return sha256_text(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True))


def stem_of(dossier: str) -> str:
    return re.sub(r"_(extraction|reconciliation)\.md$", "", dossier)


def cell_sort_key(cell_id: str) -> Tuple[str, int]:
    fam, idx = cell_id.split("#")
    return fam, int(idx)


# ------------------------------------------------------------------------------------------------------- inputs
def load_atlas() -> Dict[str, dict]:
    cells: Dict[str, dict] = {}
    for path in sorted(ATLAS.glob("F[0-9][0-9]_*.yml")):
        fam = path.name[:3]
        doc = yaml.safe_load(path.read_text())
        for i, c in enumerate(doc["cells"]):
            cells[f"{fam}#{i}"] = {**c, "_family": fam, "_index": i, "_file": path.name}
    return cells


def load_lab_ids(cells: Dict[str, dict]):
    doc = yaml.safe_load(LAB_IDS.read_text())
    string_to_lab: Dict[str, str] = {}
    aliases: Dict[str, List[str]] = {}
    for lab in doc["labs"]:
        aliases[lab["id"]] = lab.get("text_aliases") or []
        for s in lab["strings"]:
            if s in string_to_lab:
                raise SystemExit(f"lab_ids.yml: string mapped twice: {s!r}")
            string_to_lab[s] = lab["id"]
    atlas_strings = {c["lab"] for c in cells.values()}
    unmapped = sorted(atlas_strings - set(string_to_lab))
    unused = sorted(set(string_to_lab) - atlas_strings)
    if unmapped:
        raise SystemExit("lab_ids.yml: atlas lab strings with no canonical id:\n  " + "\n  ".join(unmapped))
    if unused:
        raise SystemExit("lab_ids.yml: strings no atlas cell uses (stale map):\n  " + "\n  ".join(unused))
    lab_set = set(aliases)
    collab = doc.get("collaborations") or {}
    for s, extra in collab.items():
        if s not in string_to_lab or not set(extra) <= lab_set:
            raise SystemExit(f"lab_ids.yml: bad collaboration entry {s!r}")
    copies = doc.get("data_copies") or []
    register = {e["id"]: e for e in yaml.safe_load(REGISTER.read_text())["entries"]}
    for cp in copies:
        e = register.get(cp["register"])
        if e is None or e.get("kind") != "duplicate":
            raise SystemExit(f"lab_ids.yml: {cp['register']} is not a duplicate entry of the register")
        srcs = {s.lower() for s in e["sources"]}
        if cp["copy"].lower() not in srcs or cp["original"].lower() not in srcs:
            raise SystemExit(f"lab_ids.yml: {cp['register']} does not name {cp['copy']} / {cp['original']}")
        for cid in cp["cells"]:
            stems = {stem_of(d).lower() for d in cells[cid]["dossiers"]}
            if not stems & {cp["copy"].lower(), cp["original"].lower()}:
                raise SystemExit(f"lab_ids.yml: {cid} cites neither {cp['copy']} nor {cp['original']}")
        if cp.get("adds_lab") and cp["adds_lab"] not in lab_set:
            raise SystemExit(f"lab_ids.yml: unknown lab {cp['adds_lab']}")
    return doc, string_to_lab, aliases, collab, copies


def cell_labs(cells, string_to_lab, collab, copies) -> Dict[str, List[str]]:
    """Every lab id a cell carries: lead lab first, then collaboration partners and original-data owners."""
    extra: Dict[str, Set[str]] = defaultdict(set)
    for cp in copies:
        if cp.get("adds_lab"):
            for cid in cp["cells"]:
                extra[cid].add(cp["adds_lab"])
    out = {}
    for cid, c in cells.items():
        lead = string_to_lab[c["lab"]]
        rest = (set(collab.get(c["lab"], [])) | extra[cid]) - {lead}
        out[cid] = [lead] + sorted(rest)
    return out


def cited_sections(text: str) -> str:
    """ENGINE_VS_LITERATURE secs. 2 - 5: every `##` / `###` block numbered 2.x - 5.x (wherever it sits in the file;
    sec. 5.1 was appended after sec. 6) plus Appendix A (the verification probe behind sec. 2's V*p* rows)."""
    keep, out = False, []
    for line in text.split("\n"):
        m = re.match(r"^(#{1,3})\s+(.*)$", line)
        if m:
            title = m.group(2)
            num = re.match(r"(\d+)(?:\.\d+)*\.?\s", title + " ")
            if len(m.group(1)) == 1:
                keep = False
            elif num:
                keep = num.group(1) in {"2", "3", "4", "5"}
            elif title.startswith("Appendix A"):
                keep = True
            elif len(m.group(1)) == 2:
                keep = False
        if keep:
            out.append(line)
    return "\n".join(out)


# ------------------------------------------------------------------------------------------------- citations
CELL_REF = re.compile(r"(F\d{2})?#(\d+(?:\s*(?:,|-|–)\s*\d+)*)")


def parse_cell_refs(text: str, source: str, cells) -> List[dict]:
    hits, last_fam = [], None
    for m in CELL_REF.finditer(text):
        fam = m.group(1)
        if fam:
            last_fam = fam
        elif last_fam is None or text[max(0, m.start() - 3):m.start()].strip()[-1:] not in {",", "/", ";", "("}:
            continue          # a bare '#n' only continues a list right after a family reference
        else:
            fam = last_fam
        idxs: List[int] = []
        for part in re.split(r"\s*,\s*", m.group(2)):
            r = re.match(r"(\d+)\s*[-–]\s*(\d+)$", part)
            idxs.extend(range(int(r.group(1)), int(r.group(2)) + 1) if r else [int(part)])
        line = text.count("\n", 0, m.start()) + 1
        for i in idxs:
            cid = f"{fam}#{i}"
            if cid not in cells:
                raise SystemExit(f"{source}: cites {cid}, which is not an atlas cell")
            hits.append({"source": source, "line": line, "kind": "cell", "target": cid})
    return hits


def author_year_patterns(stem: str) -> List[re.Pattern]:
    m = re.match(r"([A-Za-z]+?)(\d{4})([a-z]*)(.*)$", stem)
    name, year, suffix, rest = m.groups()
    pats = [re.compile(r"(?<![A-Za-z0-9])" + re.escape(stem) + r"(?![A-Za-z0-9])", re.I)]
    if not rest:
        pats.append(re.compile(r"\b" + re.escape(name) + r"(?: et al\.?| (?:&|and) [A-Z][a-zA-Z\-]+)?,? \(?"
                               + year + suffix + r"(?![a-z0-9])", re.I))
    return pats


def parse_named(text: str, source: str, aliases, stem_labs) -> List[dict]:
    hits = []
    for lab, pats in aliases.items():
        for p in pats:
            for m in re.finditer(r"(?<![A-Za-z])" + p + r"(?![a-z])", text):
                hits.append({"source": source, "line": text.count("\n", 0, m.start()) + 1, "kind": "lab_name",
                             "target": lab, "match": m.group(0)})
    for stem, labs in stem_labs.items():
        for p in author_year_patterns(stem):
            for m in p.finditer(text):
                hits.append({"source": source, "line": text.count("\n", 0, m.start()) + 1, "kind": "study",
                             "target": stem, "labs": sorted(labs), "match": m.group(0)})
    return hits


# ---------------------------------------------------------------------------------------- candidate contrasts
STRIP = [
    r"\bpH\s*~?\d+(?:\.\d+)?(?:\s*(?:/|->|-|–|vs\.?|to|,|and|or)\s*~?\d+(?:\.\d+)?)*",
    r"\ba_w\s*~?\d+(?:\.\d+)?(?:\s*(?:/|->|-|–|vs\.?|to|,|and|or)\s*~?\d+(?:\.\d+)?)*",
    r"~?\d+(?:\.\d+)?(?:\s*(?:/|->|-|–|vs\.?|to|,|and|or)\s*~?\d+(?:\.\d+)?)*\s*(?:C|mM|M|mol/L|mmol/L|min|h|d|days|"
    r"s|kDa|mL|kJ/mol)\b",
    r"\b(?:Table|Tables|Fig|Figs|sec|p|pp)\.?\s*[\w.\-]+", r"F\d\d#\d+", r"CF-\d+", r"\b(?:19|20)\d\d[a-z]?\b",
    r"\[?\d+[A-Z][a-z]?\d*\]?", r"\bC\d\b", r"\b\d+(?:,\d+)*-(?=[A-Za-z])", r"\b[A-Z]\+\d\b",
]
NUMBER = re.compile(r"(?<![A-Za-z\d.])~?\d+(?:\.\d+)?")
FOLD = re.compile(r"(?<![A-Za-z])x\s?~?\d|\d(?:\.\d+)?\s?x\b|-fold\b|\bfold\b|\bratio\b[^.;]{0,20}\d|\d\s*%")
CROSS_STUDY = re.compile(r"second laborator|against [A-Z][\w\-]+(?: (?:&|and) [A-Z][\w\-]+)?'s|extrapolat|"
                         r"literature value|other stud", re.I)
SHARE = re.compile(r"\d\s*%|\b\d+(?:\.\d+)?\s*:\s*\d+(?:\.\d+)?\b")


def effect_numbers(s: str) -> Tuple[int, bool]:
    t = s
    for p in STRIP:
        t = re.sub(p, " ", t)
    return len(NUMBER.findall(t)), bool(FOLD.search(t))


def factors_of(contrast: str) -> List[str]:
    head = contrast.split(": ", 1)[0]
    comps = [m.start() for m in COMPARATOR.finditer(head)]
    out = []
    for f, vocab in FACTOR_VOCAB.items():
        for m in re.finditer(vocab, head):
            if any(abs(m.start() - c) <= WINDOW for c in comps):
                out.append(f)
                break
    return out


def classify(cells, labs_of) -> Tuple[List[dict], List[dict]]:
    contrasts, excluded = [], []
    for cid in sorted(cells, key=cell_sort_key):
        c = cells[cid]
        con, qty = c.get("contrast") or "", c.get("quantity") or ""
        kind, why = None, None
        if c["evidence"] == "within_study_factor":
            n_c, f_c = effect_numbers(con)
            n_q, f_q = effect_numbers(qty)
            if CROSS_STUDY.search(con):
                why = "compares against another study, not paired arms inside one"
            elif n_c >= 2 or f_c or f_q or n_q >= 2:
                kind = "paired_arms"
            else:
                why = "no numeric effect in contrast / quantity"
        elif c["evidence"] == "labelling":
            if SHARE.search(con + " " + qty):
                kind = "route_share"
            else:
                why = "labelling without a numeric route share"
        else:
            continue
        study = stem_of(c["dossiers"][0])
        rec = {"id": f"{cid}@{study}", "cell": cid, "family": c["_family"], "study": study,
               "lab": labs_of[cid][0], "labs": labs_of[cid], "kind": kind,
               "odorant_class": ODORANT_CLASS[c["_family"]], "status": c["status"]}
        if kind:
            rec["factors"] = factors_of(con)
            contrasts.append(rec)
        else:
            rec["reason"] = why
            excluded.append(rec)
    return contrasts, excluded


def in_scope(con: dict, cell: dict, scope: dict) -> bool:
    """Does candidate contrast `con` (atlas cell `cell`) fall in an error's declared scope (ERROR_SCOPES)?"""
    fam = con["family"]
    if scope.get("families") and fam not in scope["families"]:
        return False
    if scope.get("kinds") and con["kind"] not in scope["kinds"]:
        return False
    if fam in (scope.get("all_of_family") or []):
        return True
    text = " ".join(str(cell.get(k) or "") for k in ("contrast", "quantity"))
    pat = (scope.get("text_unless_family") or {}).get(fam)
    if pat and not re.search(pat, text):
        return False
    if scope.get("factors") and not set(scope["factors"]) & set(con["factors"]):
        return False
    if scope.get("text") and not re.search(scope["text"], text):
        return False
    if scope.get("sugar_tag"):
        sugar = cell["tags"].get("sugar")
        sugars = sugar if isinstance(sugar, list) else [sugar]
        if not any(str(x).startswith(scope["sugar_tag"]) for x in sugars):
            return False
    return True


# ------------------------------------------------------------------------------------------------------- build
def build() -> Tuple[dict, dict]:
    atlas = load_atlas()
    cells = {cid: c for cid, c in atlas.items() if c["_index"] < FROZEN_EXTENT.get(c["_family"], 0)}
    post_freeze = sorted((cid for cid in atlas if cid not in cells), key=cell_sort_key)
    short = {f: n for f, n in FROZEN_EXTENT.items()
             if sum(1 for c in atlas.values() if c["_family"] == f) < n}
    if short:
        raise SystemExit(f"atlas cells deleted below the frozen extent: {short} (indices are stable; cut a new version)")
    lab_doc, string_to_lab, aliases, collab, copies = load_lab_ids(cells)
    labs_of = cell_labs(cells, string_to_lab, collab, copies)
    stem_labs: Dict[str, Set[str]] = defaultdict(set)
    for cid, c in cells.items():
        for d in c["dossiers"]:
            stem_labs[stem_of(d)].update(labs_of[cid])

    sources = {"ENGINE_VS_LITERATURE.md secs. 2 - 5": cited_sections(EVL.read_text())}
    missing_preregs = []
    for p in CITED_PREREGS:
        if p.exists():
            sources[p.name] = p.read_text()
        else:
            missing_preregs.append(p.name)
    if missing_preregs:
        raise SystemExit("cited prereg(s) missing: " + ", ".join(missing_preregs))

    citations: List[dict] = []
    for name, text in sources.items():
        citations += parse_cell_refs(text, name, cells)
        citations += parse_named(text, name, aliases, stem_labs)

    cited_cells = {h["target"] for h in citations if h["kind"] == "cell"}
    for h in citations:
        if h["kind"] == "study":
            cited_cells |= {cid for cid, c in cells.items() if h["target"] in {stem_of(d) for d in c["dossiers"]}}
    dev_labs: Dict[str, Set[str]] = defaultdict(set)    # lab -> reasons
    for cid in cited_cells:
        for lab in labs_of[cid]:
            dev_labs[lab].add(f"cell {cid}")
    for h in citations:
        if h["kind"] == "lab_name":
            dev_labs[h["target"]].add(f"named '{h['match']}' in {h['source']} l. {h['line']}")
        elif h["kind"] == "study":
            for lab in h["labs"]:
                dev_labs[lab].add(f"study '{h['match']}' in {h['source']} l. {h['line']}")
    # A cell-index citation reason names its source too.
    cell_sources: Dict[str, Set[str]] = defaultdict(set)
    for h in citations:
        if h["kind"] == "cell":
            cell_sources[h["target"]].add(h["source"])

    all_labs = sorted({l for ls in labs_of.values() for l in ls})
    dev_set = set(dev_labs)
    dev_cells = sorted((cid for cid in cells if cid in cited_cells or set(labs_of[cid]) & dev_set), key=cell_sort_key)
    test_cells = sorted((cid for cid in cells if cid not in set(dev_cells)), key=cell_sort_key)
    test_labs = sorted({labs_of[cid][0] for cid in test_cells})
    if set(test_labs) & dev_set:
        raise SystemExit("internal: a TEST cell's lead lab is a development lab")

    contrasts, excluded = classify(cells, labs_of)
    dev_cell_set = set(dev_cells)
    for con in contrasts + excluded:
        con["partition"] = "development" if con["cell"] in dev_cell_set else "test"
    test_ids = {c["id"] for c in contrasts if c["partition"] == "test"}
    if set(FACTOR_REVIEW) != test_ids:
        raise SystemExit("FACTOR_REVIEW out of step with the TEST contrasts:\n  missing "
                         + ", ".join(sorted(test_ids - set(FACTOR_REVIEW))) + "\n  stale "
                         + ", ".join(sorted(set(FACTOR_REVIEW) - test_ids)))
    for con in contrasts:
        if con["partition"] == "test":
            con["factors_heuristic"] = con["factors"]
            con["factors"] = FACTOR_REVIEW[con["id"]]

    extent = Counter(c["_family"] for c in cells.values())
    split = {
        "version": VERSION,
        "atlas_extent": dict(sorted(extent.items())),
        "development_labs": sorted(dev_set),
        "test_labs": test_labs,
        "development_cells": dev_cells,
        "test_cells": test_cells,
        "contrasts": [{"id": c["id"], "kind": c["kind"], "lab": c["lab"], "partition": c["partition"]}
                      for c in contrasts],
    }

    # ---- counts (derived; deterministic) --------------------------------------------------------------------
    def tally(items, key):
        out: Dict[str, Dict[str, int]] = {}
        for it in items:
            k = key(it)
            out.setdefault(k, {"development": 0, "test": 0})
            out[k][it["partition"]] += 1
        return dict(sorted(out.items()))

    cell_items = [{"partition": "development" if cid in dev_cell_set else "test", "family": cells[cid]["_family"],
                   "odorant_class": ODORANT_CLASS[cells[cid]["_family"]]} for cid in cells]
    test_contrasts = [c for c in contrasts if c["partition"] == "test"]
    counts = {
        "lab_strings_in_atlas": len({c["lab"] for c in cells.values()}),
        "canonical_labs": len(all_labs),
        "canonical_labs_leading_a_cell": len({labs_of[cid][0] for cid in cells}),
        "development_labs": len(dev_set),
        "test_labs": len(test_labs),
        "cells": {"total": len(cells), "development": len(dev_cells), "test": len(test_cells)},
        "contrasts": {"total": len(contrasts),
                      "development": len(contrasts) - len(test_contrasts), "test": len(test_contrasts),
                      "by_kind": tally(contrasts, lambda c: c["kind"]),
                      "test_labs_with_a_contrast": len({c["lab"] for c in test_contrasts})},
        "cells_by_family": tally(cell_items, lambda c: c["family"]),
        "cells_by_odorant_class": tally(cell_items, lambda c: c["odorant_class"]),
        "contrasts_by_family": tally(contrasts, lambda c: c["family"]),
        "contrasts_by_odorant_class": tally(contrasts, lambda c: c["odorant_class"]),
        "test_labs_by_odorant_class": {k: len({c["lab"] for c in test_contrasts if c["odorant_class"] == k})
                                       for k in sorted({c["odorant_class"] for c in contrasts})},
    }

    # ---- per-error availability (power before score; no engine run) ----------------------------------------
    errors = {}
    for eid, scope in ERROR_SCOPES.items():
        sel = [c for c in contrasts if in_scope(c, cells[c["cell"]], scope)]
        tst = [c for c in sel if c["partition"] == "test"]
        dev = [c for c in sel if c["partition"] == "development"]
        tlabs = sorted({c["lab"] for c in tst})
        errors[eid] = {
            "title": scope["title"],
            "development_contrasts": len(dev),
            "development_labs": len({c["lab"] for c in dev}),
            "test_contrasts": len(tst),
            "test_labs": len(tlabs),
            "test_lab_ids": tlabs,
            "test_contrast_ids": [c["id"] for c in tst],
            "verdict": "testable" if len(tlabs) >= 2 else "UNTESTABLE_on_held_out",
        }

    # ---- annotations (NOT part of the frozen identity; may change as other preregs are written) ------------
    exposure: Dict[str, Set[str]] = defaultdict(set)
    cited = {p.name for p in CITED_PREREGS}
    test_stems = {stem_of(cells[cid]["dossiers"][0]) for cid in test_cells} | {
        stem_of(d) for cid in test_cells for d in cells[cid]["dossiers"]}
    for p in sorted(VALIDATION.glob("*prereg*.md")):
        if p.name in cited:
            continue
        t = p.read_text()
        for stem in sorted(test_stems):
            if any(pt.search(t) for pt in author_year_patterns(stem)):
                for lab in stem_labs[stem]:
                    if lab in test_labs:
                        exposure[lab].add(p.name)
    annotations = {
        "post_freeze_cells": [{"cell": cid, "lab": atlas[cid]["lab"]} for cid in post_freeze],
        "note": ("TEST labs whose papers another pre-registration (any results/validation/*prereg*.md other than the "
                 "two cited ones) names: the engine may have been fitted or checked on them. Not a split criterion "
                 "in v1 (rule (i) counts CELLS cited by the two preregs); a scorer should report TEST results with "
                 "and without these labs."),
        "test_labs_named_by_other_preregs": {lab: sorted(v) for lab, v in sorted(exposure.items())},
        "test_labs_named_by_no_prereg": sorted(set(test_labs) - set(exposure)),
        "errors_test_labs_excluding_prereg_named": {
            eid: sorted(set(e["test_lab_ids"]) - set(exposure)) for eid, e in errors.items()},
    }

    doc = {
        "version": VERSION,
        "frozen_on": FROZEN_ON,
        "sha256": canonical_hash(split),
        "hash_rule": "sha256 of json.dumps(split, sort_keys=True, separators=(',', ':'), ensure_ascii=True)",
        "rule": ("DEVELOPMENT = every atlas cell cited (by file + index, or by naming its study or laboratory) in "
                 "ENGINE_VS_LITERATURE.md secs. 2 - 5 (incl. Appendix A and sec. 5.1) or in "
                 + " / ".join(p.name for p in CITED_PREREGS)
                 + ", plus every cell of any laboratory (lab_ids.yml) that has a cited cell; a cell carrying a "
                   "collaborating or original-data lab is DEVELOPMENT if any of its labs is. TEST = all other cells."),
        "inputs": {
            "cited_text_sha256": {k: sha256_text(v) for k, v in sources.items()},
            "lab_ids_version": lab_doc["version"],
        },
        "split": split,
        "counts": counts,
        "errors": errors,
        "contrasts": contrasts,
        "excluded_candidates": excluded,
        "cells_cited_by_index": {cid: sorted(cell_sources[cid]) for cid in sorted(cell_sources, key=cell_sort_key)},
        "development_lab_reasons": {lab: sorted(r)[:6] + ([f"... {len(r) - 6} more"] if len(r) > 6 else [])
                                    for lab, r in sorted(dev_labs.items())},
        "annotations": annotations,
    }
    return doc, split


HEADER = """# The contrast benchmark's DEVELOPMENT / TEST split. GENERATED by scripts/generators/build_contrast_split.py from
# the pathway atlas, lab_ids.yml, ENGINE_VS_LITERATURE.md secs. 2 - 5 and the two cited preregs. DO NOT HAND-EDIT.
# FROZEN: `sha256` is over the `split:` block (membership); the generator asserts it on every run. A change of the split
# needs a new, dated, reasoned version (see README.md), never an edit of this file.
"""


def dump(doc: dict) -> str:
    return HEADER + yaml.safe_dump(doc, sort_keys=False, width=120, default_flow_style=None, allow_unicode=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--write", action="store_true", help="write split.yml (refuses if membership moved)")
    args = ap.parse_args()
    doc, split = build()
    sha = doc["sha256"]
    frozen = FROZEN_SHA256 if FROZEN_SHA256 != "TO_BE_SET" else None
    text = dump(doc)
    if args.write:
        if frozen and sha != frozen:
            print(f"REFUSED: membership hash {sha} != frozen {frozen}. Cut a new dated version instead.")
            return 1
        SPLIT.write_text(text)
        print(f"wrote {SPLIT.relative_to(ROOT)}  sha256 {sha}")
        if not frozen:
            print("FROZEN_SHA256 is not set: put this hash into the script to freeze v1.")
        return 0
    ok = True
    if not frozen:
        print("FAIL: FROZEN_SHA256 not set")
        ok = False
    elif sha != frozen:
        print(f"FAIL: the split's membership drifted: rebuilt {sha} vs frozen {frozen}")
        ok = False
    if not SPLIT.exists():
        print("FAIL: split.yml missing")
        return 1
    disk = yaml.safe_load(SPLIT.read_text())
    if canonical_hash(disk["split"]) != disk["sha256"] or (frozen and disk["sha256"] != frozen):
        print("FAIL: split.yml's split block does not hash to its recorded / frozen sha256 (hand edit?)")
        ok = False
    if disk.get("inputs", {}).get("cited_text_sha256") != doc["inputs"]["cited_text_sha256"]:
        print("FAIL: the cited text (ENGINE_VS_LITERATURE secs. 2 - 5 or a cited prereg) changed since the freeze: "
              "re-review lab_ids.yml text_aliases against it and cut a new version if membership must change")
        ok = False
    for k in ("counts", "errors", "contrasts", "excluded_candidates", "cells_cited_by_index", "development_lab_reasons"):
        if disk.get(k) != doc[k]:
            print(f"FAIL: split.yml block `{k}` differs from a rebuild (regenerate with --write; membership unchanged)")
            ok = False
    if disk.get("annotations") != doc["annotations"]:
        print("note: annotations changed (other preregs moved); not part of the frozen identity")
    print(("OK" if ok else "FAILED") + f": {VERSION} sha256 {sha}; DEV cells {len(split['development_cells'])}, "
          f"TEST cells {len(split['test_cells'])}; contrasts DEV "
          f"{sum(c['partition'] == 'development' for c in split['contrasts'])} / TEST "
          f"{sum(c['partition'] == 'test' for c in split['contrasts'])}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
