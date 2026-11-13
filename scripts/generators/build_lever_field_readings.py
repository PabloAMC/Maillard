#!/usr/bin/env python
"""
The PREREG amendment A6 audit (STATUS_AND_PLAN sec. 0, step 8; spec: docs/reference/lever_table/STEP8_FIELDS_SPEC.md).

Derives three facts per contrast or arm that the keyed free text holds only implicitly, overrides them with the
resolved reader checks in ``docs/reference/lever_table/field_checks/``, and reports what the frozen candidate rules
A6-R1 - R3 would read on the sources ALREADY read. An audit: no reading of the lever table moves, and nothing here is
evidence for the rules (they are scored only on sources keyed after the A6 commit).

- F1 ``ph_control`` per pH contrast: buffered | adjusted_unbuffered | unadjusted | unstated (| mixed, if arms differ);
- F2 ``final_ph`` per arm: arms keyed ``ph_basis: final``, else the checks;
- F3 the base mixture per arm, as compound classes (amine, sulfur, carbonyl, lipid, other) of the reactants other
  than the arm's own lever level: the checks' ``present`` list, else the keyed ``precursors``.

Writes ``results/validation/lever_field_readings.{json,md}``. The lever-table generator is imported, not edited (it
is an input of frozen records).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import build_lever_table as LT  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "lever_field_readings.json"
CHECKS_DIR = data_paths.LEVER_TABLE_DIR / "field_checks"
MAIN_TABLE = data_paths.VALIDATION_DIR / "lever_table.json"
AUDIT = data_paths.VALIDATION_DIR / "lever_disagreement_audit.json"
DELTA_FINAL_MIN = 0.5          # A6-R2: a slope on final pH only where the arms' final pH spans this much
UNBUFFERED = ("adjusted_unbuffered", "unadjusted")
CONTROL_VALUES = ("buffered", "adjusted_unbuffered", "unadjusted", "unstated")


# ------------------------------------------------------------------------------------------------- F1, F2, F3
def control_from_text(text: Any) -> str:
    """F1 from one arm's keyed `buffer` text (STEP8_FIELDS_SPEC sec. 2). A bare 'none' cannot tell titrated from
    untouched and reads 'unstated' unless a check resolves it."""
    t = str(text or "").strip().lower()
    if t in ("", "unstated", "not stated", "none stated"):
        return "unstated"
    if t.startswith(("none", "no buffer", "water")):
        if "adjust" in t:
            return "adjusted_unbuffered"
        return "unstated" if t == "none" else "unadjusted"
    return "buffered"


def contrast_control(arms: Sequence[Mapping[str, Any]], check: Optional[Mapping[str, Any]], cid: str) -> Tuple[str, str]:
    """(F1 value, where it came from)."""
    if check and "F1" in (check.get("fields") or []):
        v = (check.get("ph_control_overrides") or {}).get(cid) or check.get("ph_control")
        if v in CONTROL_VALUES:
            return v, "check"
    states = {control_from_text((a.get("conditions") or {}).get("buffer")) for a in arms}
    return (states.pop() if len(states) == 1 else "mixed"), "derived"


def arm_final_ph(arm: Mapping[str, Any], check: Optional[Mapping[str, Any]]) -> Optional[float]:
    cond = arm.get("conditions") or {}
    if cond.get("ph_basis") == "final" and LT._num(cond.get("ph")) is not None:
        return LT._num(cond.get("ph"))
    if check and "F2" in (check.get("fields") or []):
        rec = (check.get("arms") or {}).get(str(arm.get("label")))
        if isinstance(rec, Mapping):
            return LT._num(rec.get("final_ph"))
    return None


_CLASS_WORDS = (
    ("sulfur", ("cystein", "cystine", "glutathione", "gsh", "thiamin", "methionin", "sulfide", "h2s", "sulfur",
                "sulphur", "thiol")),
    ("lipid", ("lipid", "oil", "fat", "linole", "linolen", "oleic", "phospholipid", "triglycerid", "lecithin",
               "phosphatidyl", "hexanal", "enal")),
    ("carbonyl", ("ose", "sugar", "dextrin", "starch", "glyoxal", "dicarbonyl", "dione", "furfural",
                  "monophosphate", "imp", "gmp", "aldehyde")),
    ("amine", ("ine", "amino", "protein", "peptide", "hydrolys", "amadori", "ammoni", "amine", "isolate", "gluten",
               "casein", "whey", "meal", "powder", "yeast")),
)


def classes(name: str) -> set:
    """Compound classes of one reactant name, by keyword. Amadori products and glutathione carry two or three."""
    n = _expand(name)
    out = {c for c, words in _CLASS_WORDS if any(w in n for w in words)}
    if "amadori" in n or "arp" in re.findall(r"[a-z]+", n):
        out |= {"amine", "carbonyl"}
    if n in ("water", "buffer", "phosphate", "") or n.startswith(("buffer", "water", "ethanol", "acetonitrile")):
        return set()
    return out or {"other"}


#: abbreviations the keyed precursor names use (normalised before matching and classing)
ABBREVIATIONS = {"gsh": "glutathione", "cys": "cysteine", "gcys": "cystine", "cys-cys": "cystine", "vb1": "thiamine",
                 "xyl": "xylose", "glc": "glucose", "rib": "ribose", "met": "methionine", "gly": "glycine",
                 "tpp": "thiamine pyrophosphate", "arp": "amadori product", "rg-arp": "amadori product",
                 "nf": "norfuraneol", "imp": "inosine 5'-monophosphate", "gmp": "guanosine 5'-monophosphate"}


def _expand(s: str) -> str:
    t = str(s).strip().lower()
    for pre in ("l-", "d-"):
        if t.startswith(pre):
            t = t[len(pre):]
    head = t.split(" ")[0]
    return ABBREVIATIONS[head] + t[len(head):] if head in ABBREVIATIONS else t


def _core(s: str) -> str:
    """A reactant or level name without amounts, parentheses, salt forms or an L-/D- prefix."""
    t = LT._level(re.sub(r"\(.*?\)", " ", _expand(s)))
    t = re.sub(r"\b(hcl|hydrochloride|monohydrate)\b", " ", t)
    t = re.sub(r"-[ld]-", "", t).replace("-", " ")
    t = re.sub(r"[0-9.]+\s*(mmol|mol|mm|m|g|mg|%)\b", " ", t)
    return " ".join(t.split())


def _same_reactant(name: str, level: str) -> bool:
    """Whether a reactant is (part of) the arm's own lever level: 'glycine' is the level 'glycine (2 mmol)'; 'lysine'
    is part of 'free amino acids (lysine)'."""
    n, lv = _core(name), _core(level)
    full = re.sub(r"-[ld]-", "", LT._level(_expand(level))).replace("-", " ")
    parts = {" ".join(x.split()) for x in re.split(r"[+,()/]+", full)} | set(full.split())
    return bool(n) and (n == lv or n in parts or (bool(lv) and n.startswith(lv + " ")))


def _zero_amount(v: Any) -> bool:
    """A keyed precursor printed with amount zero (a zero-dose arm) is not in the pot."""
    x = LT._num(v)
    return x is not None and x == 0.0


def arm_base(arm: Mapping[str, Any], check: Optional[Mapping[str, Any]]) -> Tuple[frozenset, str]:
    """F3: the compound classes of the arm's reactants other than its own lever level."""
    level = LT._level(arm.get("lever_value") or arm.get("label"))
    names = None
    src = "derived"
    if check and "F3" in (check.get("fields") or []):
        rec = (check.get("arms") or {}).get(str(arm.get("label")))
        if isinstance(rec, Mapping) and rec.get("present") is not None:
            names, src = [str(x) for x in rec.get("present") or []], "check"
    if names is None:
        prec = (arm.get("conditions") or {}).get("precursors") or {}
        names = ([str(k) for k, v in prec.items() if not _zero_amount(v)] if isinstance(prec, Mapping) else [])
    kept = [n for n in names if not _same_reactant(n, level)] if level != "none" else names
    out = set()
    for n in kept:
        out |= classes(n)
    return frozenset(out), src


# ------------------------------------------------------------------------------------------------- loading
def load_checks(directory: Path = CHECKS_DIR) -> Dict[str, Dict[str, Any]]:
    """Resolved check files, keyed by source; a `<source>.F3.yml` file is merged into `<source>.yml`."""
    out: Dict[str, Dict[str, Any]] = {}
    if not directory.exists():
        return out
    for f in sorted(directory.glob("*.yml")):
        if ".reader_" in f.name:
            continue
        doc = yaml.safe_load(f.read_text()) or {}
        if doc.get("reader") != "resolved":
            raise SystemExit(f"{f.name}: not a resolved check file (reader: {doc.get('reader')!r})")
        src = str(doc["source"])
        cur = out.setdefault(src, {"source": src, "fields": [], "arms": {}})
        cur["fields"] = sorted(set(cur["fields"]) | set(doc.get("fields") or []))
        for k in ("ph_control", "ph_control_overrides"):
            if doc.get(k) is not None:
                cur[k] = doc[k]
        for label, rec in (doc.get("arms") or {}).items():
            cur["arms"].setdefault(str(label), {}).update(rec or {})
    return out


def observations_with_contrasts(sources, labs) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    """The lever table's observations, each carrying the contrast it came from. Built one resolved FILE at a time:
    contrast ids repeat across the files of one source (xu2007.yml and xu2007__families.yml both hold xu2007-3)."""
    obs, skipped = [], Counter()
    for doc in sources:
        by_id = {str(c.get("id")): c for c in doc["contrasts"] or []}
        o_doc, sk = LT.observations([doc], labs)
        for o in o_doc:
            o["_contrast"] = by_id[str(o["contrast"])]
        obs += o_doc
        skipped.update(sk)
    return obs, dict(skipped)


# ------------------------------------------------------------------------------------------------- readings
def annotate(obs, checks) -> Dict[str, Any]:
    """Attach F1 / F2 / F3 to each observation (in place); return coverage counts."""
    cov = {"ph_contrasts": 0, "ph_control": Counter(), "ph_control_source": Counter(), "final_all_arms": 0,
           "final_delta_ge_min": 0, "categorical_pairs": 0, "pairs_base_mismatch": 0, "base_source": Counter()}
    for o in obs:
        c = o["_contrast"]
        arms = c.get("arms") or []
        chk = checks.get(o["source"])
        if o["lever_class"] == "pH":
            cov["ph_contrasts"] += 1
            o["ph_control"], how = contrast_control(arms, chk, str(c.get("id")))
            cov["ph_control"][o["ph_control"]] += 1
            cov["ph_control_source"][how] += 1
            fin = [arm_final_ph(a, chk) for a in arms]
            o["final_ph"] = fin
            if fin and all(f is not None for f in fin):
                cov["final_all_arms"] += 1
                o["delta_final"] = max(fin) - min(fin)
                if o["delta_final"] >= DELTA_FINAL_MIN:
                    o["final_slope"] = final_ph_slope(c, arms, fin)
                    cov["final_delta_ge_min"] += int(o["final_slope"] is not None)
        if o.get("kind") == "categorical":
            bases = {}
            for a in arms:
                b, how = arm_base(a, chk)
                bases[LT._level(a.get("lever_value") or a.get("label"))] = b
                cov["base_source"][how] += 1
            o["base_by_level"] = {k: sorted(v) for k, v in bases.items()}
            for p in o["pairs"]:
                cov["categorical_pairs"] += 1
                p["base_matched"] = bases.get(p["level"]) == bases.get(p["against"])
                cov["pairs_base_mismatch"] += int(not p["base_matched"])
    return {k: (dict(v) if isinstance(v, Counter) else v) for k, v in cov.items()}


def final_ph_slope(c, arms, fin) -> Optional[Tuple[float, float]]:
    """A6-R2: log10 amount on the final pH (never the step of the main table), over the longest monotone run as
    PREREG sec. 4 takes it on the keyed basis."""
    pts = []
    for a, x in zip(arms, fin):
        ly, var, censor = LT.arm_log_and_var(a, bool(c.get("digitised")))
        if ly is not None and x is not None:
            pts.append((x, ly, var))
    pts.sort()
    ds = [q[1] - p[1] for p, q in zip(pts, pts[1:])]
    monotone = all(d >= 0 for d in ds) or all(d <= 0 for d in ds)
    run = pts if monotone else LT._longest_monotone(pts)    # PREREG sec. 4's run, as the main table
    if not run or max(p[0] for p in run) - min(p[0] for p in run) < DELTA_FINAL_MIN:
        return None          # A6-R2: the spread is checked on the arms the slope uses, not on the whole ladder
    return LT._wls_slope(run)


def main_readings() -> Dict[Tuple[str, str, str], Dict[str, Any]]:
    table = json.loads(MAIN_TABLE.read_text(encoding="utf-8"))
    out = {}
    for cell in table["cells"]:
        key = (cell["compound"], cell["lever_class"], cell["band"])
        if cell["lever_class"] in LT.NUMERIC:
            out[key] = {"reading": cell["reading"], "flags": cell.get("fragility", {}).get("flags", []),
                        "estimate": cell.get("fit", {}).get("estimate"), "labs": cell["labs"]}
        else:
            out[key] = {"levels": {lv: s.get("reading") for lv, s in cell["fit"]["levels"].items()},
                        "labs": cell["labs"]}
    return out


def ph_cell_readings(rows, main) -> Dict[str, Any]:
    """A6-R1 (bands by keyed control) and A6-R2 (final-pH slopes where available) for one pH cell."""
    out = {"contrasts": len(rows), "labs": len({r["lab"] for r in rows}),
           "ph_control": dict(Counter(r.get("ph_control") for r in rows)),
           "free_text_buffer_state": dict(Counter(r.get("buffer_state") for r in rows)),
           "main_reading": main.get("reading"), "main_flags": main.get("flags", [])}
    out["free_text_buffered_only"] = LT._alt_reading([r for r in rows if r.get("buffer_state") == "buffered"])
    out["R1_same_set_as_free_text"] = ({(r["source"], r["contrast"]) for r in rows if r.get("ph_control") == "buffered"}
                                       == {(r["source"], r["contrast"]) for r in rows
                                           if r.get("buffer_state") == "buffered"})
    out["R1"] = {"buffered": LT._alt_reading([r for r in rows if r.get("ph_control") == "buffered"]),
                 "unbuffered": LT._alt_reading([r for r in rows if r.get("ph_control") in UNBUFFERED])}
    swapped = [{**r, "y": r["final_slope"][0], "v": r["final_slope"][1]} if r.get("final_slope") else r for r in rows]
    n_sw = sum(1 for r in rows if r.get("final_slope"))
    out["R2"] = {"contrasts_on_final_ph": n_sw,
                 "reading": LT._alt_reading(swapped) if n_sw else None,
                 "per_contrast": [{"source": r["source"], "contrast": r["contrast"], "lab": r["lab"],
                                   "initial_range": r.get("range"), "final_ph": r.get("final_ph"),
                                   "y_keyed": r["y"], "y_final": (r.get("final_slope") or [None])[0]}
                                  for r in rows if r.get("final_ph") and any(f is not None for f in r["final_ph"])]}
    changes = []
    for band in ("buffered", "unbuffered"):
        rd = out["R1"][band]["reading"]
        if out["labs"] >= LT.MIN_LABS and rd != out["main_reading"]:
            changes.append(f"R1 {band} band: {rd}")
    if out["R2"]["reading"] and out["labs"] >= LT.MIN_LABS and out["R2"]["reading"]["reading"] != out["main_reading"]:
        changes.append(f"R2 final pH: {out['R2']['reading']['reading']}")
    out["changes"] = changes
    return out


def categorical_cell_readings(rows, lc, main) -> Dict[str, Any]:
    """A6-R3: each level's reading without the pairs whose base mixture differs beyond the lever."""
    ref = LT.CATEGORICAL_REFERENCE[lc]
    n_mis = sum(1 for r in rows for p in r["pairs"] if not p.get("base_matched", True))
    out = {"pairs": sum(len(r["pairs"]) for r in rows), "pairs_base_mismatch": n_mis,
           "reference_pot_by_lab": {}, "changes": []}
    for r in rows:
        b = r.get("base_by_level") or {}
        for p in r["pairs"]:
            out["reference_pot_by_lab"].setdefault(r["lab"], set()).add(
                f"{p['against']}: {'+'.join(b.get(p['against'], [])) or 'nothing else'}")
    out["reference_pot_by_lab"] = {k: sorted(v) for k, v in out["reference_pot_by_lab"].items()}
    if not n_mis:
        return out
    kept = [{**r, "pairs": [p for p in r["pairs"] if p.get("base_matched", True)]} for r in rows]
    kept = [r for r in kept if r["pairs"]]
    alt = LT.categorical_fit(kept, ref, LT.PRIOR_SCALE)["levels"] if kept else {}
    for lv, m in (main.get("levels") or {}).items():
        s = alt.get(lv, {})
        r2 = LT.reading(s["labs"], s) if "estimate" in s else "too few laboratories"
        if m and m != "not connected" and r2 != m:
            out["changes"].append(f"R3 {lv}: {m} -> {r2}")
    return out


def watched_cells() -> Dict[Tuple[str, str, str], str]:
    """Cells the step-8 record must report whatever happens: the answered and fragile cells of the main table and
    the audit's conflicts and within-laboratory sign changes."""
    out = {}
    table = json.loads(MAIN_TABLE.read_text(encoding="utf-8"))
    for cell in table["cells"]:
        key = (cell["compound"], cell["lever_class"], cell["band"])
        if cell.get("reading") == "answered by the literature":
            out[key] = "answered"
        if cell.get("fragility", {}).get("flags"):
            out[key] = "answered, fragile" if key in out else "fragile"
        for lv, s in (cell.get("fit", {}).get("levels") or {}).items():
            if s.get("reading") == "answered by the literature":
                out[key] = "answered (level)"
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    for c in audit["cells"]:
        if c["class"] in ("conflict between laboratories", "sign changes inside one laboratory"):
            out.setdefault((c["compound"], c["lever_class"], "all" if c["lever_class"] in LT.CATEGORICAL_REFERENCE
                            else c["band"]), c["class"])
    return out


# ------------------------------------------------------------------------------------------------- build
def build() -> Dict[str, Any]:
    sources = LT.load_sources()
    labs, _ = LT.lab_partition()
    obs, _ = observations_with_contrasts(sources, labs)
    checks = load_checks()
    coverage = annotate(obs, checks)
    main = main_readings()
    cells = defaultdict(list)
    for o in obs:
        cells[(o["compound"], o["lever_class"], o["band"])].append(o)
    watched = watched_cells()
    ph, cat = [], []
    for key, rows in sorted(cells.items()):
        compound, lc, band = key
        if lc == "pH" and len({r["lab"] for r in rows}) >= 2:
            ph.append({"compound": compound, "band": band, "watched": watched.get(key),
                       **ph_cell_readings(rows, main.get(key, {}))})
        elif lc in LT.CATEGORICAL_REFERENCE and len({r["lab"] for r in rows}) >= 2:
            r3 = categorical_cell_readings(rows, lc, main.get(key, {}))
            if r3["pairs_base_mismatch"] or key in watched:
                cat.append({"compound": compound, "lever_class": lc, "watched": watched.get(key), **r3})
    return {"artifact": "lever_field_readings", "status": "audit on sources already read; no reading of the lever "
            "table moves; PREREG A6 rules are scored only on sources keyed after the A6 commit",
            "spec": "docs/reference/lever_table/STEP8_FIELDS_SPEC.md", "checked_sources": sorted(checks),
            "coverage": coverage, "ph_cells": ph, "categorical_cells": cat}


def _f(x):
    return "-" if x is None else f"{x:+.2f}"


def render(p: Mapping[str, Any]) -> str:
    cov = p["coverage"]
    L = ["# PREREG A6 audit: the keyed pH control, final pH and base mixture (STATUS_AND_PLAN sec. 0, step 8)", "",
         f"*{p['status']}. Spec: [{Path(p['spec']).name}](../../{p['spec']}).*", "",
         "## Coverage", "",
         f"- pH contrasts: {cov['ph_contrasts']}; `ph_control`: "
         + ", ".join(f"{k} {v}" for k, v in sorted(cov['ph_control'].items()))
         + f" ({', '.join(f'{k} {v}' for k, v in sorted(cov['ph_control_source'].items()))}).",
         f"- final pH on every arm: {cov['final_all_arms']} contrasts; final-pH spread of at least {DELTA_FINAL_MIN} "
         f"pH unit over the arms the slope uses (A6-R2 slope reported): {cov['final_delta_ge_min']}.",
         f"- categorical pairs: {cov['categorical_pairs']}; base mixture differs beyond the lever (A6-R3 drops "
         f"them): {cov['pairs_base_mismatch']}; arms' base from "
         + ", ".join(f"{k} {v}" for k, v in sorted(cov['base_source'].items())) + ".",
         f"- reader-checked sources: {', '.join(p['checked_sources']) or 'none'}.", "",
         "## pH cells with two or more laboratories", "",
         "| odorant | band | watched | labs | control (keyed) | main reading | R1 buffered | buffered set = free text's? "
         "| R1 unbuffered | R2 final pH (contrasts) | changes |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in p["ph_cells"]:
        r2 = c["R2"]["reading"]
        L.append(f"| {c['compound']} | {c['band']} | {c['watched'] or ''} | {c['labs']} | "
                 + ", ".join(f"{k} {v}" for k, v in sorted(c['ph_control'].items())) + f" | {c['main_reading']} | "
                 f"{c['R1']['buffered']['reading']} ({c['R1']['buffered']['labs']} labs) | "
                 f"{'yes' if c['R1_same_set_as_free_text'] else 'NO'} | "
                 f"{c['R1']['unbuffered']['reading']} ({c['R1']['unbuffered']['labs']} labs) | "
                 f"{(r2 or {}).get('reading', '-')} ({c['R2']['contrasts_on_final_ph']}) | "
                 f"{'; '.join(c['changes']) or 'none'} |")
    L += ["", "### Contrasts that carry a final pH (initial range against the pH the pots ended at)", "",
          "| odorant | band | source | laboratory | initial or label pH | final pH by arm | slope, keyed basis | "
          "slope on final pH |", "|---|---|---|---|---|---|---|---|"]
    for c in p["ph_cells"]:
        for r in c["R2"]["per_contrast"]:
            fin = ", ".join("-" if f is None else f"{f:g}" for f in r["final_ph"])
            rng = r["initial_range"]
            L.append(f"| {c['compound']} | {c['band']} | {r['source']} | {r['lab']} | "
                     f"{rng[0]:g} - {rng[1]:g} | {fin} | {_f(r['y_keyed'])} | {_f(r['y_final'])} |")
    L += ["", "## Categorical cells: base mixture (A6-R3)", "",
          "| odorant | lever | watched | pairs | base differs | reference pot by laboratory | changes |",
          "|---|---|---|---|---|---|---|"]
    for c in p["categorical_cells"]:
        pots = "; ".join(f"{lab}: {', '.join(v)}" for lab, v in sorted(c["reference_pot_by_lab"].items()))
        L.append(f"| {c['compound']} | {c['lever_class']} | {c['watched'] or ''} | {c['pairs']} | "
                 f"{c['pairs_base_mismatch']} | {pots} | {'; '.join(c['changes']) or 'none'} |")
    L += ["", "R1 - R3 are PREREG amendment A6's frozen candidate rules, read here on contrasts already seen: a change "
          "is a finding about this table's sensitivity, not evidence for the rule. A reading needs three "
          "laboratories (PREREG sec. 5); 'changes' lists only readings that differ from the main table's."]
    return "\n".join(L) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.parse_args(argv)
    payload = build()
    inputs = ([p for p in sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml"))
               if p.name != "TEMPLATE.yml" and ".reader_" not in p.name]
              + [p for p in sorted(CHECKS_DIR.glob("*.yml")) if ".reader_" not in p.name]
              + [data_paths.CONTRAST_SPLIT, MAIN_TABLE, AUDIT])
    payload["provenance"] = provenance.provenance_block(
        OUT.stem, generated_by="scripts/generators/build_lever_field_readings.py", inputs=inputs)
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(payload['ph_cells'])} pH cells, "
          f"{len(payload['categorical_cells'])} categorical cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
