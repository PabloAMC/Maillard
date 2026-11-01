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
from functools import lru_cache
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


@lru_cache(maxsize=None)
def author_year_patterns(stem: str) -> Tuple[re.Pattern, ...]:
    m = re.match(r"([A-Za-z]+?)(\d{4})([a-z]*)(.*)$", stem)
    name, year, suffix, rest = m.groups()
    pats = [re.compile(r"(?<![A-Za-z0-9])" + re.escape(stem) + r"(?![A-Za-z0-9])", re.I)]
    if not rest:
        pats.append(re.compile(r"\b" + re.escape(name) + r"(?: et al\.?| (?:&|and) [A-Z][a-zA-Z\-]+)?,? \(?"
                               + year + suffix + r"(?![a-z0-9])", re.I))
    return tuple(pats)


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


# ================================================================================================================
# v2: PER-CHANGE hold-out (coordinator, 2026-11-01, declared before any scoring). v1 above stays frozen and untouched.
#   DEV(E)  = the cells cited AS EVIDENCE FOR E in ENGINE_VS_LITERATURE secs. 2 - 5 (E's mismatch rows in sec. 2.1, its
#             rows in secs. 3 / 3.1 / 5 / 5.1, its recommended-order item, its Appendix-A probes, its sec. 2.0
#             correction), plus E's own prereg cells, plus every cell of those cells' / named laboratories;
#   TEST(E) = the candidate contrasts in E's scope (ERROR_SCOPES_V2) from laboratories not in DEV(E), EXCLUDING any
#             contrast whose paper feeds the shipped engine (a calibration input: CALIBRATION_* below).
# ================================================================================================================
VERSION_V2 = "contrast_split_v2_per_change"
SPLIT_V2 = BENCH / "split_v2_per_change.yml"
FROZEN_V2_SHA256 = "dbe32a2b3af83f3bbb603efb391dd8477999263a6b24769ea946fcc0e0fe2e12"

# Which sec. 2.1 mismatch rows are evidence for which error: sec. 3's first column, written out (sec. 3 abbreviates
# "F16-a, -b"). F28-d (named in E7's sec. 3 row) has no sec. 2.1 row. F21-a carries both the pyrrole (E14) and the
# pyridine (E20). F06-d (the latent dpo_share split) belongs to E10: sec. 5 makes B59b's dpo_share fix E10's lever.
MISMATCH_ERRORS = {
    "F23-a": ["E1"], "F24-a": ["E1"], "F06-a": ["E2"], "F06-c": ["E2"], "F16-a": ["E3"], "F16-b": ["E3"],
    "F16-c": ["E4"], "F01-b": ["E5"], "F03-b": ["E5"], "F13-b": ["E5"], "F14-a": ["E5"], "F19-c": ["E5"],
    "F23-c": ["E5"], "F29-b": ["E5", "E9"], "F23-b": ["E6"], "F30-a": ["E7"], "F26-c": ["E7"], "F25-a": ["E8"],
    "F25-b": ["E8"], "F06-b": ["E8"], "F29-a": ["E9"], "F29-c": ["E9"], "F04-a": ["E10"], "F04-b": ["E10"],
    "F06-d": ["E10"], "F23-d": ["E11"], "F28-a": ["E12"], "F28-b": ["E12"], "F28-c": ["E12"], "F19-a": ["E13"],
    "F19-b": ["E14"], "F21-a": ["E14", "E20"], "F20-a": ["E14"], "F17-a": ["E15"], "F13-a": ["E16"],
    "F03-a": ["E17"], "F26-a": ["E18"], "F27-a": ["E18"], "F27-b": ["E18"], "F22-a": ["E19"], "F18-a": ["E20"],
}
# Sec. 2.0's two corrections, and Appendix A's probes (A: NF share; B: FFT vs pH; C, E: hexose MFT; D: MFT / FFT).
SEC20_ITEM_ERRORS = {1: ["E8"], 2: ["E10"]}
PROBE_ERRORS = {"A": ["E6"], "B": ["E8"], "C": ["E1"], "D": ["E6", "E8"], "E": ["E1"]}
PREREG_ERRORS = {CITED_PREREGS[0].name: ["E4"], CITED_PREREGS[1].name: ["E6", "E11"]}
# Evidence the cited text gives without a cell index or a name. Sec. 2.1 row F03-a: "hexose ARP 71 - 96 % (3 labs)"
# = F03#0 (Nestle, 0.74 - 0.88), F03#3 (Wageningen, 71 - 96 %), F03#11 (Ho network, 88 %). Sec. 5's E5 "8+ labs"
# names none and cannot be resolved (README caveat).
UNNAMED_EVIDENCE = {"E17": ["F03#0", "F03#3", "F03#11"]}

# v2 scopes = v1's, with two widenings that make each error's TEST set match its own prereg's search (declared now,
# before any v2 count was read): E4 also F17 and F21 (the E4 prereg searched F16, F17, F19, F20, F21 for pH
# contrasts); E6 every F23 route share (the E6 prereg scores intact-C5 / xylose shares, e.g. Wang 2012's, which the
# change moves), not only those whose text names NF.
ERROR_SCOPES_V2 = {**ERROR_SCOPES,
                   "E4": {**ERROR_SCOPES["E4"], "families": ["F16", "F17", "F19", "F20", "F21"]},
                   "E6": {"title": ERROR_SCOPES["E6"]["title"], "families": ["F23", "F26"], "kinds": ["route_share"],
                          "text_unless_family": {"F23": None, "F26": ERROR_SCOPES["E6"]["text"]}}}

# Calibration inputs (a test must not be one). Paper level (conservative: one shipped row excludes the paper):
#   (1) measurement-ledger entries `in_model_fit` without a not-shipped flag (data/lit/measurement_ledger/*.yml);
#   (2) the shipped sulfur fit B60's 52 objective rows, by their row-id prefix (the ledger predates part of B60);
#   (3) papers whose measured constant the engine ships FIXED (SRC_FIXED_CONSTANTS, read in src/kinetic_core).
# Recall net for (3): every atlas paper the engine's source cites (src/kinetic_core/*.py) that would enter a TEST(E)
# set must be reviewed into SRC_FIXED_CONSTANTS or SRC_MENTION_ONLY; the generator fails until it is.
LEDGER_DIR = ROOT / "data" / "lit" / "measurement_ledger"
B60_REPORT = VALIDATION / "kinetic_core_b60_fit_report.json"
B60_ROW_PREFIX_STEMS = {
    "hofmann_": ["hofmann1998"], "fed_": ["hofmann1998"], "zhou_": ["zhou2023"], "whitfield_": ["whitfield1999",
    "whitfield2001"], "vanseeventer_": ["vanseeventer2001"], "zhang_fig1_": ["Zhang2024"], "cerny_": ["cerny2007b"],
    "kang_": ["kang2026", "kang2026_SI"], "kumazawa_": ["kumazawa2003"], "feng_": ["feng2022"],
    "zhai_": ["zhai2023foodchem", "zhai2023jafc", "zhai2021"],
}   # zhai_*_fold_120_over_100 are Zhai 2023 Food Chem Table 1; its 100 C arm reprints Zhai 2021, as do JAFC's (CF-066)
ENGINE_SRC = ROOT / "src" / "kinetic_core"
SRC_FIXED_CONSTANTS = {
    "chan1994b": "2-acetyl-1-pyrroline Ea 60.2 kJ/mol (parameters_proline.py) and the measured sulfur-class Ea prior "
                 "(ea_class_prior.py, Table I)",
    "arnold1971": "B50: the share of degraded thiamine released as H2S vs pH (parameters_sulfur.py ~l. 402, sulfur.py)",
    "stack2018": "the thiol / o-quinone K(T) and forward Ea, implemented (parameters_sulfur.py ~l. 1204 - 1706)",
}
SRC_MENTION_ONLY = {   # cited in src/kinetic_core, read 2026-11-01: no shipped constant or fit row comes from it
    "du2023": "documents why HEXOSE_C5_ENTRY exists; the switch ships OFF (k = 0)",
    "zhao2019b": "documents why HEXOSE_C5_ENTRY exists; the switch ships OFF (k = 0)",
    "farmer1990": "cited as the lipid quench behind a REFUSED saturable-thioether constant (parameters_sulfur.py)",
    "whitfield1988": "cited as the lipid quench behind a REFUSED saturable-thioether constant (parameters_sulfur.py)",
    "meynier1995": "B51 records the buffer (pyrophosphate) its pots used: a scoring condition, not a fitted value",
    "mundt2004": "C/N 7.64 quoted as a comparison against the fixed MEL_C 8 : MEL_N 1 (species.py)",
    "quan2020": "its 88 rate constants are REFUSED (parameters_acrylamide.py); its glyoxal level is a check "
                "(ledger in_model_check)",
    "zamora2015": "the src 'Zamora 2015 (Food Chem 174:89)' is zamora2015b, a different paper (DECADIENAL_BREAKAGE)",
}

# Hand-reviewed factor tags (2026-11-01) for contrasts that v1 left on the heuristic and that could enter a factor-
# scoped v2 TEST set (E2 / E4 / E8 pH, E5 buffer, E7 lipid, E16 pH / O2). The generator lists, for each such error,
# every in-family candidate from a non-DEV(E), non-calibration lab whose text carries the factor's vocabulary ANYWHERE
# (recall), and fails until each has a review here; the review then decides (precision). Notes as comments.
FACTOR_REVIEW_V2: Dict[str, List[str]] = {
    "F04#2@kallen1971b": ["buffer"],                       # phosphate 0 -> 0.3 M at pH ~6.1 - 6.6, saturating
    "F04#19@zhang2025b": ["amine"],                        # lysine dose; the phosphate is held
    "F05#23@thornalley1999": ["amine", "buffer", "metals"],  # +/- HSA; phosphate 10 -> 100 mM; DETAPAC
    "F06#11@thornalley1999": ["buffer", "metals", "amine"],
    "F08#2@thornalley1999": ["metals", "buffer"],
    "F15#9@cerny2003": [],                                 # one condition, a label pattern
    "F16#6@weenen2001": ["ph"],                            # pH 7 phosphate vs pH 3.2 citrate: buffer co-varies
    "F16#11@monforte2021": ["ph", "metals", "sugar", "temperature"],  # pH 7 > 5 > 3.4; metals; glucose; 40 - 80 C
    "F22#7@li2022": ["time"],
    "F22#12@shu1985": ["ph"],                              # unbuffered pH 2.2 / 5.1 / 7.1
    "F25#17@us8920862": ["lipid"],                         # structured lipid vs aqueous buffer (patent)
    "F25#18@grosch1992": ["concentration", "temperature"],  # meat-like dilution at 100 C vs Farmer-style at 140 C
    "F27#7@dwivedi1972": ["ph"],                           # one buffer (0.1 M phosphate) across the pH ladder
    "F27#15@grosch1992": ["amine", "concentration"],       # thiamine vs ribose + Cys; Cys dose
    "F27#16@grosch1992": ["amine"],                        # TPP vs free thiamine (precursor identity)
    "F28#12@grosch1992": [],                               # reflux vs SDE: a method contrast
    "F28#20@xu2011": ["lipid"],                            # no / fresh / oxidised lard
    "F28#25@wang2022": ["temperature"],
    "F30#5@xu2011": ["lipid"],
    "F30#7@du2025": ["lipid"],
    "F30#11@kim1996b": ["amine"],                          # amino-acid identity with a dienal
    "F30#15@zamora2015": ["oxygen", "lipid"],              # N2 vs air; lipid identity; one buffer (citrate pH 3)
    "F31#3@mundt2004": ["sugar"],                          # glucose vs maltose, label positions; acetate held
}


def _segments_by_error(text: str) -> Dict[str, str]:
    """Split ENGINE_VS_LITERATURE secs. 2 - 5 (as cited_sections() returns it) into the text that is evidence for each
    error. Lines that belong to no error (headings, the mismatch-kind legend, the lower-ranked leads, sec. 4) drop out."""
    out: Dict[str, List[str]] = defaultdict(list)
    sec, item_errs = None, []
    for line in text.split("\n"):
        h = re.match(r"^#{2,3}\s+(Appendix A|\d+(?:\.\d+)?)", line)
        if h:
            sec, item_errs = h.group(1), []
            continue
        errs: List[str] = []
        if sec == "2.0":
            m = re.match(r"^(\d+)\.\s", line)
            if m:
                item_errs = SEC20_ITEM_ERRORS.get(int(m.group(1)), [])
            errs = item_errs
        elif sec == "2.1":
            m = re.match(r"^\| (F\d\d-[a-z]) \|", line)
            errs = MISMATCH_ERRORS.get(m.group(1), []) if m else []
        elif sec == "3":
            m = re.match(r"^\| \d+ \| (E\d+)\b", line)
            errs = [m.group(1)] if m else []
        elif sec == "3.1":
            m = re.match(r"^\| (E\d+)(?: - (E\d+))? \|", line)
            if m:
                a = int(m.group(1)[1:])
                b = int(m.group(2)[1:]) if m.group(2) else a
                errs = [f"E{i}" for i in range(a, b + 1)]
        elif sec == "Appendix A":
            m = re.match(r"^\| ([A-E]) \|", line)
            errs = PROBE_ERRORS[m.group(1)] if m else []
        elif sec == "5":
            m = re.match(r"^\| (E\d+) ", line)
            o = re.match(r"^(\d+)\.\s+\*\*(.*?)\*\*", line)
            if m:
                errs, item_errs = [m.group(1)], []
            elif o:
                item_errs = re.findall(r"E\d+", o.group(2))
                errs = item_errs
            elif line.startswith("   ") and item_errs:
                errs = item_errs
            else:
                item_errs = []
        elif sec == "5.1":
            errs = ["E4"]
        for e in errs:
            out[e].append(line)
    unmapped = set(MISMATCH_ERRORS) - set(re.findall(r"^\| (F\d\d-[a-z]) \|", text, re.M))
    if unmapped:
        raise SystemExit(f"MISMATCH_ERRORS names rows absent from sec. 2.1: {sorted(unmapped)}")
    return {e: "\n".join(v) for e, v in out.items()}


def _prereg_parts(text: str) -> Tuple[str, str, str]:
    """(everything but sec. 4.1, sec. 4.1 'Held-out', sec. 4.2 'Development checks')."""
    m = re.search(r"^### 4\.1 .*?$(.*?)(?=^### )", text, re.M | re.S)
    d = re.search(r"^### 4\.2 .*?$(.*?)(?=^##)", text, re.M | re.S)
    if not m or not d:
        raise SystemExit("cited prereg lacks its '### 4.1' held-out or '### 4.2' development-check section")
    return text[:m.start()] + text[m.end():], m.group(0), d.group(0)


def _calibration_inputs(stems: Set[str]) -> Dict[str, List[str]]:
    """atlas dossier stem -> the calibration sources that make it a calibration input."""
    out: Dict[str, Set[str]] = defaultdict(set)
    lower = {s.lower(): s for s in stems}
    for path in sorted(LEDGER_DIR.glob("*.yml")):
        for e in yaml.safe_load(path.read_text())["entries"]:
            if e["status"] != "in_model_fit":
                continue
            if any(f.startswith("fit_wave_not_shipped") or f == "fit_not_forwarded_route_stage1" for f in e["flags"]):
                continue
            s = lower.get(stem_of(e["dossier"]).lower())
            if s:
                out[s].add(f"ledger {path.name}: {e['id']}")
    rows = json.loads(B60_REPORT.read_text())["members"][0]["row_ids"]
    for r in rows:
        hit = [p for p in B60_ROW_PREFIX_STEMS if r.startswith(p)]
        if not hit:
            raise SystemExit(f"B60 row {r} has no declared source paper (B60_ROW_PREFIX_STEMS)")
        for s in B60_ROW_PREFIX_STEMS[hit[0]]:
            if s in stems:
                out[s].add(f"B60 fit row {hit[0]}*")
    for s, why in SRC_FIXED_CONSTANTS.items():
        if s not in stems:
            raise SystemExit(f"SRC_FIXED_CONSTANTS names {s}, not an atlas paper")
        out[s].add(f"fixed constant: {why}")
    # A data set published twice is one calibration input: a register copy of a calibration paper (or the original
    # of one) is a calibration input too (lab_ids.yml data_copies; e.g. Schieberle 1998 Table IV's aqueous arm is
    # Hofmann 1998 T1, a B60 fit row).
    norm = lambda x: lower.get(re.sub(r"_(extraction|reconciliation)$", "", x).lower())
    for cp in yaml.safe_load(LAB_IDS.read_text()).get("data_copies") or []:
        a, b = norm(cp["copy"]), norm(cp["original"])
        for x, y in ((a, b), (b, a)):
            if x and y and y in out and not any(t.startswith("data copy") for t in out[y]):
                out[x].add(f"data copy ({cp['register']}) of calibration input {y}")
    return {s: sorted(v) for s, v in sorted(out.items())}


def _src_cited(stems: Set[str]) -> Set[str]:
    src = "\n".join(p.read_text() for p in sorted(ENGINE_SRC.glob("*.py")))
    return {s for s in stems if any(pt.search(src) for pt in author_year_patterns(s))}


_CALIB_CACHE: Dict[frozenset, tuple] = {}


def build_v2(variant: str = "declared") -> dict:
    """variant: declared | ames_separate | kallen_not_cited | prereg_literal | calibration_narrow | calibration_broad."""
    atlas = load_atlas()
    cells = {cid: c for cid, c in atlas.items() if c["_index"] < FROZEN_EXTENT.get(c["_family"], 0)}
    lab_doc, string_to_lab, aliases, collab, copies = load_lab_ids(cells)
    aliases = {k: list(v) for k, v in aliases.items()}
    if variant == "ames_separate":
        for s in ("Univ. Reading / Dijon (Ames)", "University of Reading (Ames)"):
            string_to_lab[s] = "reading_ames"
        aliases["reading_mottram_parker"].remove("Ames")
        aliases["reading_ames"] = ["Ames"]
    if variant == "kallen_not_cited":
        aliases["kallen_pennsylvania"] = []
    labs_of = cell_labs(cells, string_to_lab, collab, copies)
    stem_labs: Dict[str, Set[str]] = defaultdict(set)
    for cid, c in cells.items():
        for d in c["dossiers"]:
            stem_labs[stem_of(d)].update(labs_of[cid])
    stem_cells: Dict[str, Set[str]] = defaultdict(set)
    for cid, c in cells.items():
        for d in c["dossiers"]:
            stem_cells[stem_of(d)].add(cid)

    evl = cited_sections(EVL.read_text())
    segments = _segments_by_error(evl)
    texts: Dict[str, Dict[str, str]] = defaultdict(dict)    # error -> source -> text
    for e, t in segments.items():
        texts[e]["ENGINE_VS_LITERATURE secs. 2 - 5"] = t
    prereg_parts: Dict[str, Tuple[str, str, str]] = {}
    for p in CITED_PREREGS:
        body, held, devchk = _prereg_parts(p.read_text())
        prereg_parts[p.name] = (body, held, devchk)
        for e in PREREG_ERRORS[p.name]:
            texts[e][p.name] = p.read_text() if variant == "prereg_literal" else body

    def key(h):
        return h["kind"], h["target"]

    # A prereg's declared held-out targets (cited in its sec. 4.1 and not among its sec. 4.2 development checks) do
    # not become its development, wherever else the prereg mentions them (summary, baselines, expectations).
    exempt: Dict[str, Set[Tuple[str, str]]] = {}
    for name, (body, held, devchk) in prereg_parts.items():
        hk = {key(h) for h in parse_cell_refs(held, name, cells) + parse_named(held, name, aliases, stem_labs)}
        dk = {key(h) for h in parse_cell_refs(devchk, name, cells) + parse_named(devchk, name, aliases, stem_labs)}
        exempt[name] = set() if variant == "prereg_literal" else hk - dk

    def targets(text: str, source: str):
        hits = parse_cell_refs(text, source, cells) + parse_named(text, source, aliases, stem_labs)
        hits = [h for h in hits if key(h) not in exempt.get(source, set())]
        tcells, tlabs = set(), set()
        for h in hits:
            if h["kind"] == "cell":
                tcells.add(h["target"])
            elif h["kind"] == "study":
                tcells |= stem_cells[h["target"]]
            else:
                tlabs.add(h["target"])
        return tcells, tlabs, hits

    contrasts, _ = classify(cells, labs_of)
    v1_test = {c["id"] for c in contrasts} & set(FACTOR_REVIEW)
    for c in contrasts:
        c["factors"] = (FACTOR_REVIEW.get(c["id"]) if c["id"] in v1_test else
                        FACTOR_REVIEW_V2.get(c["id"], c["factors"]))
        c["factors_source"] = ("v1_review" if c["id"] in v1_test else
                               "v2_review" if c["id"] in FACTOR_REVIEW_V2 else "heuristic")

    all_stems = {stem_of(d) for c in cells.values() for d in c["dossiers"]}
    key_ = frozenset(all_stems)
    if key_ not in _CALIB_CACHE:
        _CALIB_CACHE[key_] = (_calibration_inputs(all_stems), _src_cited(all_stems))
    calib, src_cited = _CALIB_CACHE[key_]
    calib = {k: list(v) for k, v in calib.items()}
    if variant == "calibration_narrow":
        calib = {s: [x for x in v if not x.startswith("fixed constant")] for s, v in calib.items()}
        calib = {s: v for s, v in calib.items() if v}
    if variant == "calibration_broad":
        for s_ in src_cited:
            calib.setdefault(s_, ["cited in src/kinetic_core"])

    errors, review_needed = {}, {}
    for eid, scope in ERROR_SCOPES_V2.items():
        cited_cells, named_labs, prov = set(), set(), []
        for src_name, t in texts.get(eid, {}).items():
            tc, tl, hits = targets(t, src_name)
            cited_cells |= tc
            named_labs |= tl
            prov += [f"{h['kind']} {h.get('match', h['target'])} ({src_name})" for h in hits]
        cited_cells |= set(UNNAMED_EVIDENCE.get(eid, []))
        held_out_declared = sorted({f"{k}:{t}" for p, es in PREREG_ERRORS.items() if eid in es
                                    for k, t in exempt[p]})
        dev_labs = {l for cid in cited_cells for l in labs_of[cid]} | named_labs
        dev_cells = {cid for cid in cells if cid in cited_cells or set(labs_of[cid]) & dev_labs}
        in_scope_all = [c for c in contrasts if in_scope(c, cells[c["cell"]], scope)]
        test, excluded_cal = [], []
        for c in contrasts:
            if c["cell"] in dev_cells:
                continue
            stems = {stem_of(d) for d in cells[c["cell"]]["dossiers"]}
            cal = sorted(s for s in stems if s in calib)
            if c in in_scope_all:
                (excluded_cal if cal else test).append({**c, "calibration_papers": cal})
                unrev = sorted(st for st in stems if st in src_cited and st not in calib
                               and st not in SRC_MENTION_ONLY)
                if not cal and unrev and variant == "declared":
                    review_needed.setdefault(f"SRC:{','.join(unrev)}", set()).add(eid)
            # recall net for factor-scoped errors: any family member whose text carries the factor vocabulary
            if scope.get("factors") and not cal and c["id"] not in v1_test and c["id"] not in FACTOR_REVIEW_V2:
                fam_ok = not scope.get("families") or c["family"] in scope["families"]
                txt = " ".join(str(cells[c["cell"]].get(k) or "") for k in ("contrast", "quantity"))
                if fam_ok and any(re.search(FACTOR_VOCAB[f], txt) for f in scope["factors"]):
                    review_needed.setdefault(c["id"], set()).add(eid)
        tlabs = sorted({c["lab"] for c in test})
        errors[eid] = {
            "title": scope["title"],
            "development_labs": sorted(dev_labs),
            "development_cells_cited": sorted(cited_cells, key=cell_sort_key),
            "prereg_declared_heldout_not_dev": held_out_declared,
            "test_contrasts": [c["id"] for c in test],
            "test_labs": tlabs,
            "n_test_contrasts": len(test),
            "n_test_labs": len(tlabs),
            "excluded_as_calibration_input": [f"{c['id']} ({', '.join(c['calibration_papers'])})"
                                              for c in excluded_cal],
            "verdict": "testable" if len(tlabs) >= 2 else "UNTESTABLE_on_held_out",
        }
    return {"errors": errors, "calibration_inputs": calib, "review_needed": review_needed,
            "contrasts": {c["id"]: c for c in contrasts}, "cells": cells}


HEADER_V2 = """# The contrast benchmark's v2 PER-CHANGE hold-out (one DEV / TEST split per open engine error E1 - E20). GENERATED by
# scripts/generators/build_contrast_split.py --v2 --write. DO NOT HAND-EDIT. v1 (split.yml) stays frozen beside it.
# FROZEN: `sha256` is over the `split:` block; `--v2` asserts it. A change needs a new, dated, reasoned version.
"""


def v2_document() -> dict:
    main_run = build_v2("declared")
    if main_run["review_needed"]:
        lines = [f"  {cid}  ({', '.join(sorted(es))})" for cid, es in sorted(main_run["review_needed"].items())]
        raise SystemExit("FACTOR_REVIEW_V2 needs these contrasts reviewed:\n" + "\n".join(lines))
    errs = main_run["errors"]
    v1_doc = yaml.safe_load(SPLIT.read_text())
    v1_test = [c["id"] for c in v1_doc["contrasts"] if c["partition"] == "test"]
    v2_test = sorted({cid for e in errs.values() for cid in e["test_contrasts"]},
                     key=lambda x: cell_sort_key(x.split("@")[0]))
    reread = sorted(set(v1_test) | set(v2_test), key=lambda x: cell_sort_key(x.split("@")[0]))
    split = {
        "version": VERSION_V2,
        "atlas_extent": dict(sorted(FROZEN_EXTENT.items())),
        "calibration_input_papers": sorted(main_run["calibration_inputs"]),
        "errors": {eid: {"development_labs": e["development_labs"], "test_contrasts": e["test_contrasts"]}
                   for eid, e in errs.items()},
    }
    alternatives = {}
    for var in ("ames_separate", "kallen_not_cited", "prereg_literal", "calibration_narrow", "calibration_broad"):
        r = build_v2(var)["errors"]
        alternatives[var] = {eid: {"n_test_labs": e["n_test_labs"], "n_test_contrasts": e["n_test_contrasts"],
                                   "test_labs": e["test_labs"], "verdict": e["verdict"]} for eid, e in r.items()}
    contrasts = main_run["contrasts"]
    return {
        "version": VERSION_V2,
        "frozen_on": FROZEN_ON,
        "sha256": canonical_hash(split),
        "hash_rule": "sha256 of json.dumps(split, sort_keys=True, separators=(',', ':'), ensure_ascii=True)",
        "rule": ("Per error E: DEV(E) = the cells cited as evidence for E in ENGINE_VS_LITERATURE.md secs. 2 - 5 "
                 "(its sec. 2.1 mismatch rows, its rows in secs. 3, 3.1, 5 and 5.1, its recommended-order item, its "
                 "Appendix-A probes, its sec. 2.0 correction) and in E's own prereg outside that prereg's sec. 4.1 "
                 "held-out list, plus every cell of those cells' and named laboratories (lab_ids.yml). TEST(E) = "
                 "the candidate contrasts in E's scope from laboratories not in DEV(E), excluding any contrast "
                 "whose paper is a calibration input of the shipped engine."),
        "inputs": {"cited_text_sha256": yaml.safe_load(SPLIT.read_text())["inputs"]["cited_text_sha256"],
                   "v1_sha256": v1_doc["sha256"], "lab_ids_version": "lab_ids_v1"},
        "split": split,
        "summary": {
            "testable_errors": [eid for eid, e in errs.items() if e["verdict"] == "testable"],
            "per_error": {eid: f"{e['n_test_labs']} labs / {e['n_test_contrasts']} contrasts ({e['verdict']})"
                          for eid, e in errs.items()},
            "v2_test_contrasts_union": len(v2_test),
            "reread_needed": {"v1_test": len(v1_test), "v2_new": len(set(v2_test) - set(v1_test)),
                              "total": len(reread)},
        },
        "errors": errs,
        "alternatives": alternatives,
        "reread_list": [{"id": cid, "in_v1_test": cid in v1_test,
                         "errors": [eid for eid, e in errs.items() if cid in e["test_contrasts"]],
                         "lab": contrasts[cid]["lab"] if cid in contrasts else None} for cid in reread],
        "test_contrast_factors": {cid: {"factors": contrasts[cid]["factors"], "source": contrasts[cid]["factors_source"]}
                                  for cid in v2_test},
        "calibration_inputs": main_run["calibration_inputs"],
    }


def main_v2(write: bool) -> int:
    doc = v2_document()
    sha = doc["sha256"]
    frozen = FROZEN_V2_SHA256 if FROZEN_V2_SHA256 != "TO_BE_SET" else None
    text = HEADER_V2 + yaml.safe_dump(doc, sort_keys=False, width=120, default_flow_style=None, allow_unicode=False)
    if write:
        if frozen and sha != frozen:
            print(f"REFUSED: v2 hash {sha} != frozen {frozen}. Cut a new dated version instead.")
            return 1
        SPLIT_V2.write_text(text)
        print(f"wrote {SPLIT_V2.relative_to(ROOT)}  sha256 {sha}" + ("" if frozen else "  (FROZEN_V2_SHA256 unset)"))
        return 0
    ok = bool(frozen) and sha == frozen
    if not ok:
        print(f"FAIL: v2 rebuilt {sha} vs frozen {frozen}")
    if not SPLIT_V2.exists():
        print("FAIL: split_v2_per_change.yml missing")
        return 1
    disk = yaml.safe_load(SPLIT_V2.read_text())
    if canonical_hash(disk["split"]) != disk["sha256"] or disk["sha256"] != frozen:
        print("FAIL: split_v2_per_change.yml's split block does not hash to its recorded / frozen sha256")
        ok = False
    for k in ("inputs", "summary", "errors", "alternatives", "reread_list", "test_contrast_factors",
              "calibration_inputs"):
        if disk.get(k) != doc[k]:
            print(f"FAIL: block `{k}` differs from a rebuild")
            ok = False
    print(("OK" if ok else "FAILED") + f": {VERSION_V2} sha256 {sha}; testable {doc['summary']['testable_errors']}")
    return 0 if ok else 1


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
    ap.add_argument("--v2", action="store_true", help="the v2 per-change split (split_v2_per_change.yml) instead")
    args = ap.parse_args()
    if args.v2:
        return main_v2(args.write)
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
