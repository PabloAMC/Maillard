#!/usr/bin/env python
"""
The experiment map (STATUS_AND_PLAN sec. 0, step 4): for every lever class x odorant family a formulator asks about,
one reading, and where the answer is missing, the smallest experiment that would supply it.

Readings (PREREG sec. 5 and STATUS_AND_PLAN sec. 0):
  * "answered by the literature"  -- the lever table answers the cell (>= 3 laboratories, the direction shown);
  * "answered by the engine"      -- left blank in v1: step 5 alone adds it;
  * "needs experiment"            -- everything else, with the reason: direction unclear across laboratories, too
                                     few laboratories, keyed but no contrast usable, or never measured on disk.

The design line is one line, no simulation (sec. 0 step 4): the lever and its levels, the replicate number the
cell's own spread implies (two-arm, alpha 0.05 one-sided, power 0.8, against the measurement SD of PREREG sec. 5
for the effect the cell suggests, else 0.3 dex), in a NEW laboratory (one laboratory is the binding constraint, not
replicates), and in a plant-protein matrix wherever the cell has no plant-protein contrast.

Inputs: results/validation/lever_table.json (step 2), docs/reference/contrast_benchmark/split.yml (how many atlas
contrasts exist on disk for the odorant families not yet keyed). Writes results/validation/experiment_map.{json,md}.
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "experiment_map.json"
LEVER_TABLE = data_paths.VALIDATION_DIR / "lever_table.json"

#: The formulator's levers (STATUS_AND_PLAN sec. 1 and W5.4) and the odorant families of the plan.
LEVERS = ["pH", "temperature", "time", "precursor dose", "sugar identity", "sulfur source", "yeast extract",
          "thiamine", "lipid", "water", "buffer", "heme / iron", "matrix", "cooking method"]
FAMILIES = {
    "meaty thiols": ["2-Methyl-3-furanthiol (MFT)", "2-Furfurylthiol (FFT)", "Hydrogen Sulfide",
                     "Bis(2-methyl-3-furyl) disulfide", "Bis(2-furfuryl) disulfide",
                     "2-Methyl-3-(2-furfuryldithio)furan"],
    "Strecker aldehydes, methional": [],
    "pyrazines": [],
    "lipid aldehydes (off-notes)": [],
    "furanones, furfural": [],
}
#: the split's odorant classes -> the families above (for "keyed later" counts)
SPLIT_CLASS = {"thiol_h2s_sulfur": "meaty thiols", "strecker_aldehyde_acid": "Strecker aldehydes, methional",
               "pyrazine_n_heterocycle": "pyrazines", "lipid_crossover": "lipid aldehydes (off-notes)",
               "furans_furanones": "furanones, furfural"}
SPLIT_FACTOR = {"ph": "pH", "temperature": "temperature", "time": "time", "concentration": "precursor dose",
                "sugar": "sugar identity", "amine": "sulfur source", "lipid": "lipid", "water": "water",
                "buffer": "buffer", "metals": "heme / iron", "oxygen": None}
#: levels a design line proposes per lever (the formulator's range; a convention, stated as such)
LEVELS = {"pH": "pH 5.0 / 6.0 / 7.0", "temperature": "120 / 140 / 160 C", "time": "5 / 10 / 20 / 40 min",
          "precursor dose": "0 / 0.5x / 1x / 2x of the base dose", "sugar identity": "ribose vs glucose (vs xylose)",
          "sulfur source": "cysteine vs glutathione vs thiamine vs none", "yeast extract": "0 / 0.5 / 1 / 2 % w/w",
          "thiamine": "0 / 0.05 / 0.1 % w/w", "lipid": "0 / 2 / 5 % sunflower oil", "water": "a_w 0.95 / 0.85 / 0.75",
          "buffer": "phosphate 0 / 0.1 / 0.2 M", "heme / iron": "0 / 0.2 / 0.5 % leghemoglobin or FeSO4 equivalent",
          "matrix": "buffer vs 15 % pea isolate", "cooking method": "sealed 140 C vs pan 180 C surface"}
ARM_SD = math.log10(1.2) / math.sqrt(2)   # PREREG sec. 5 per-arm measurement SD (20 % CV, n = 2)
Z = 1.645 + 0.842                          # one-sided alpha 0.05, power 0.8


def replicates(effect_dex: float) -> int:
    d = max(abs(effect_dex), 0.1)
    return max(2, math.ceil(2 * (Z * ARM_SD / d) ** 2 * 2))   # x2: the between-replicate pot SD ~ measurement SD


def _plant(matrices) -> bool:
    return any(any(k in str(m).lower() for k in ("pea", "soy", "plant", "isolate", "gluten", "chickpea"))
               for m in matrices)


def lever_of(cell: Mapping[str, Any]) -> str:
    lc = cell["lever_class"]
    return {"buffer concentration": "buffer"}.get(lc, lc)


def build() -> Dict[str, Any]:
    lt = json.loads(LEVER_TABLE.read_text())
    split = yaml.safe_load(data_paths.CONTRAST_SPLIT.read_text())
    atlas = defaultdict(lambda: {"contrasts": 0, "labs": set()})
    for c in split["contrasts"]:
        fam = SPLIT_CLASS.get(c.get("odorant_class"))
        if fam and c["kind"] == "paired_arms":
            for f in c.get("factors") or []:
                lv = SPLIT_FACTOR.get(f)
                if lv:
                    atlas[(fam, lv)]["contrasts"] += 1
                    atlas[(fam, lv)]["labs"].add(c["lab"])

    # the lever table's evidence per (family, lever): the best reading among its odorants / bands / levels
    rank = {"answered by the literature": 3, "direction unclear": 2, "too few laboratories": 1}
    keyed: Dict = defaultdict(list)
    for cell in lt["cells"]:
        fam = next((f for f, names in FAMILIES.items() if cell["compound"] in names), None)
        if not fam:
            continue
        lv = lever_of(cell)
        if "step" in cell:
            keyed[(fam, lv)].append({"what": f"{cell['compound']}, {cell['band']}", "reading": cell["reading"],
                                     "labs": cell["labs"], "estimate": cell["fit"]["estimate"],
                                     "new_lab_pi90": cell["fit"]["new_lab_pi90"], "plant": _plant(cell["matrices"]),
                                     "rows": [r["y"] for r in cell["rows"]]})
        else:
            for level, s in cell["fit"]["levels"].items():
                if "estimate" not in s:
                    continue
                # a thiamine level of the sulfur-source cell is the formulator's thiamine lever too
                for lv2 in ([lv, "thiamine"] if "thiamine" in level else [lv]):
                    keyed[(fam, lv2)].append({"what": f"{cell['compound']}, {level} vs {cell['reference']}",
                                         "reading": s["reading"], "labs": s["labs"], "estimate": s["estimate"],
                                          "new_lab_pi90": s["new_lab_pi90"], "plant": _plant(cell["matrices"]),
                                          "rows": []})
    # keyed thiamine contrasts that gave no estimate (censored "none" arm, or no chain to cysteine)
    keyed_thiamine = []
    for f in sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml")):
        if ".reader_" in f.name or f.name == "TEMPLATE.yml":
            continue
        doc = yaml.safe_load(f.read_text()) or {}
        for c in doc.get("contrasts") or []:
            if any("thiamin" in str(a.get("lever_value") or "").lower() for a in c.get("arms") or []):
                keyed_thiamine.append(doc["source"])
    rows = []
    for fam in FAMILIES:
        for lv in LEVERS:
            ev = keyed.get((fam, lv), [])
            at = atlas.get((fam, lv), {"contrasts": 0, "labs": set()})
            answered = [e for e in ev if e["reading"] == "answered by the literature"]
            best = max(ev, key=lambda e: (rank.get(e["reading"], 0), e["labs"]), default=None)
            plant = any(e["plant"] for e in ev)
            if answered:
                reading = "answered by the literature"
                why = "; ".join(f"{e['what']}: {e['estimate']:+.2f} dex (new laboratory {e['new_lab_pi90'][0]:+.2f} "
                                f"to {e['new_lab_pi90'][1]:+.2f})" for e in answered)
            elif ev:
                reading = "needs experiment"
                why = (f"keyed: best cell '{best['what']}' reads {best['reading']} ({best['labs']} laboratories)")
            elif lv == "thiamine" and fam == "meaty thiols" and keyed_thiamine:
                reading = "needs experiment"
                why = (f"keyed in {', '.join(sorted(set(keyed_thiamine)))} but no usable estimate: the no-thiamine arm "
                       "is 'not detected' without a limit, or thiamine is compared only with thiamine pyrophosphate "
                       "and never connects to cysteine")
            elif fam != "meaty thiols" and at["contrasts"]:
                reading = "not yet keyed"
                why = f"{at['contrasts']} atlas contrasts from {len(at['labs'])} laboratories await keying (step 2)"
            else:
                reading = "needs experiment"
                why = ("no atlas contrast on disk (the measurement ledger may hold unkeyed ones)" if not at["contrasts"]
                       else f"{at['contrasts']} atlas contrasts, none keyed as usable for this lever")
            design = None
            if reading == "needs experiment":
                eff = best["estimate"] if best else 0.3
                labs_needed = max(1, 3 - (best["labs"] if best else 0))
                design = (f"{lv}: {LEVELS[lv]}; {replicates(eff)} replicate pots per arm; "
                          f"{'in a 15 % pea-isolate matrix' if not plant else 'matrix as in the keyed studies'}; "
                          f"from {labs_needed} new laboratory(ies) to reach 3; measure "
                          f"{'MFT, FFT, H2S and the disulfides by SIDA' if fam == 'meaty thiols' else fam}")
            rows.append({"family": fam, "lever": lv, "reading": reading, "why": why, "plant_protein_evidence": plant,
                         "design": design, "engine": "blank in v1 (step 5 adds)"})
    counts = defaultdict(int)
    for r in rows:
        counts[r["reading"]] += 1
    return {"artifact": "experiment_map", "lever_table": str(LEVER_TABLE.relative_to(ROOT)), "rows": rows,
            "counts": dict(counts)}


def render(p: Mapping[str, Any]) -> str:
    L = ["# The experiment map", "",
         "*Generated by `scripts/generators/build_experiment_map.py` (STATUS_AND_PLAN sec. 0, step 4) from the lever "
         "table (`lever_table.md`). For each lever a formulator asks about and each odorant family: what the "
         "literature answers, and where it does not, the smallest experiment that would. The engine's column is "
         "blank in v1; step 5 alone fills it.*", "",
         "**Read this first.** " + ", ".join(f"{v} cells {k}" for k, v in sorted(p["counts"].items())) + ". "
         "No cell has a heated plant-protein contrast on the keyed record, so every design below names a pea-isolate "
         "matrix: that, and laboratories, are the binding gaps, not replicates.", "",
         "**One study covers most of the meaty-thiol column.** A single factorial in a 15 % pea-isolate matrix "
         "(sugar ribose vs glucose x sulfur source cysteine vs thiamine x yeast extract 0 vs 1 % x pH 5.5 vs 6.5, "
         "2 replicates, 32 pots, MFT / FFT / H2S / disulfides by SIDA) adds one laboratory to eight thiol cells at "
         "once, and is the matrix the keyed record lacks. That is the plan's P3 study (STATUS_AND_PLAN sec. 4.1); "
         "it brings each cell one laboratory closer, so cells keyed at two laboratories reach the three-laboratory "
         "rule.", ""]
    for fam in dict.fromkeys(r["family"] for r in p["rows"]):
        L += [f"## {fam}", "", "| lever | reading | evidence | smallest experiment |", "|---|---|---|---|"]
        for r in p["rows"]:
            if r["family"] == fam:
                L.append(f"| {r['lever']} | {r['reading']} | {r['why']} | {r['design'] or '-'} |")
        L.append("")
    L += ["Design conventions (stated, not fitted): levels are the formulator's usual range; replicates from a "
          "two-arm one-sided test at alpha 0.05 and power 0.8 against the PREREG measurement SD for the effect the "
          "cell suggests (0.1 dex floor; 0.3 dex where nothing is keyed); 'laboratories' counts toward the lever "
          "table's three-laboratory rule. Sending any request is the owner's action.", ""]
    return "\n".join(L)


def main() -> int:
    payload = build()
    payload["provenance"] = provenance.provenance_block(
        "experiment_map", generated_by="scripts/generators/build_experiment_map.py",
        inputs=[LEVER_TABLE, data_paths.CONTRAST_SPLIT])
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)}: {payload['counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
