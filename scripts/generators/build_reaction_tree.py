#!/usr/bin/env python
"""
The reaction tree the model carries, drawn from the code (2026-09-07), for docs/guides/REACTION_TREES.md.

One figure per lane (sugar, pentose-cysteine, acrylamide). Species are boxes; every step in the lane's
reaction list is an arrow, coloured by HOW ITS RATE CONSTANT IS KNOWN:

  measured directly          a published rate for this step (evidence_class measured_rate / measured_activation_energy)
  fitted, identified         a fitted coordinate the data pin down (Laplace / fit-report standard error)
  fitted, not identified     a fitted coordinate the data leave on its band (a band artefact)
  fitted earlier, frozen     fitted in an earlier calibration and carried unchanged
  bounded from a timescale   only a bracket is known (evidence_class bounded_from_a_timescale_bracket)
  no constant                an instantaneous or lumped step (parameter_key None) or a structural constant

Species boxes are coloured by HOW WELL THE PANEL MEASURES THEM: the median fold error over the
scorecard rows for that compound (green <= 3x, amber 3-10x, red > 10x); grey = never measured on the
panel. Bookkeeping pools (fragments, carboxylate, acid, oxidant, oligomer) are hidden.

    python scripts/generators/build_reaction_tree.py     # writes docs/assets/thiol_sink/1[0-2]_tree_*.png + 13_steps_by_status.png
"""
from __future__ import annotations

import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
from src import data_paths  # noqa: E402
from src.kinetic_core import acrylamide as acr_mod, network, species, species_acrylamide, species_sulfur, sulfur as sul_mod  # noqa: E402
from src.kinetic_core.engine import ACRYLAMIDE, SULFUR, TARGET_ALIASES, TRUNK, core_parameters  # noqa: E402

OUT = ROOT / "docs" / "assets" / "thiol_sink"
V = data_paths.VALIDATION_DIR
INK, MUTED = "#1E2A2C", "#5E6B6E"
HIDDEN = {"FRAG_C", "FRAG_N", "FRAG_S", "CBX", "ACID", "OLG", "OX", "OXR", "OXV", "MELE", "PROT_SS", "PRB", "MEL_N", "MEL"}
STATUS_STYLE = {   # the same scale as the field scheme in build_thiol_sink_figures.FIELD_STATUS
    "rate known at several temperatures (measured, or fitted with its barrier)": ("#2B5DA8", "-", 2.0),
    "rate known at one temperature only (fitted at 145 °C; no measured barrier)": ("#178F6E", "-", 2.0),
    "carried from an earlier calibration, not re-examined": ("#9AA6A3", "-", 1.6),
    "only a band or a bracket is known": ("#D9822B", "--", 2.0),
    "no rate constant (instantaneous or lumped)": ("#C9D0CD", "-", 1.0),
}
NODE_FILL = {"good": "#D9EFE3", "mid": "#FBE9D0", "bad": "#F6D9D9", "none": "#F2F3F1"}
NODE_EDGE = {"good": "#178F6E", "mid": "#D9822B", "bad": "#B23A3A", "none": "#9AA6A3"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150})


def _read(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# status of every rate constant
# ---------------------------------------------------------------------------
def constant_status(lane: str) -> Dict[str, str]:
    params = core_parameters(lane)
    priors = _read(V / "core_prediction_uncertainty.json")["priors"]
    fit_free, identified = set(), {}
    if lane == SULFUR:
        rep = _read(V / "kinetic_core_b9_fit_report.json")
        fit_free = set(rep["free_set"]["keys"])
        lap = _read(V / "kinetic_core_b9_laplace_covariance.json")
        identified = {(c["key"] or c["block"]): bool(i) for c, i in zip(lap["coordinates"], lap["identified"])}
    reason: Dict[str, Dict[str, str]] = defaultdict(dict)      # key -> {coordinate kind: reason}
    for p in priors:
        parts = p["key"].split(".")
        if len(parts) >= 3:
            reason[parts[1]][parts[2]] = p["reason"].lower()
    S = list(STATUS_STYLE)
    out = {}
    for key, param in params.items():
        ec = getattr(param, "evidence_class", None)
        # a step whose constant the shipped model carries at exactly zero (the inert defaults of the
        # refused sink variants: the disulfide release, the pot-made electrophile sites) carries no
        # flux and is not drawn; the tree is the network the model integrates, not the one it could
        if getattr(param, "k_ref", None) == 0.0 and ec == "derived_from_fit_data":
            out[key] = "inert"
            continue
        if ec in ("measured_rate", "measured_activation_energy"):
            out[key] = S[0]
        elif "fitted_wave_b18" in (getattr(param, "flags", None) or ()):
            # the pyrazine step's two Strecker constants: rate and barrier fitted together on a
            # fed-dicarbonyl ladder at three temperatures (Laplace sigma 0.08 dex); the envelope
            # does not sample them yet, so the priors table has no row for them
            out[key] = S[0]
        elif ec == "bounded_from_a_timescale_bracket":
            out[key] = S[3]
        elif ec == "derived_from_fit_data":
            if lane == SULFUR:
                if key in fit_free:
                    out[key] = S[1] if identified.get(key, False) else S[3]
                else:
                    out[key] = S[2]
            else:
                r = reason.get(key, {})
                rate = next((v for k, v in r.items() if "log10_k" in k), "")
                ea = next((v for k, v in r.items() if "ea" in k), "")
                if not r:
                    out[key] = S[2]
                elif "unidentified" in rate or "drawn" in rate:
                    out[key] = S[3]
                elif ea and ("unidentified" in ea or "drawn" in ea):
                    out[key] = S[1]
                else:
                    out[key] = S[0]
        else:
            out[key] = S[4]
    return out


# ---------------------------------------------------------------------------
# how well each species is measured on the panel
# ---------------------------------------------------------------------------
def species_quality() -> Dict[str, str]:
    sc = _read(V / "core_panel_scores.json")
    folds = defaultdict(list)
    for b in sc["benchmarks"]:
        for c in b["compounds"]:
            if c.get("fold_error") and c.get("predicted"):
                folds[c["compound"].lower()].append(float(c["fold_error"]))
    inv = defaultdict(list)
    for name, key in TARGET_ALIASES.items():
        inv[key].append(name.lower())
    out = {}
    for key, names in inv.items():
        vals = []
        for comp, f in folds.items():
            if any(n in comp for n in names):
                vals += f
        if vals:
            m = statistics.median(vals)
            out[key] = "good" if m <= 3 else ("mid" if m <= 10 else "bad")
    return out


# ---------------------------------------------------------------------------
# layout
# ---------------------------------------------------------------------------
#: Curated left-to-right stages per lane: reactants -> first intermediates -> fragments -> aroma compounds -> removal.
STAGES = {
    SULFUR: [
        ["PENT", "Cys", "THI", "Glc", "Gly"],
        ["TTCA", "ARP", "H2S", "MESH", "HMF"],
        ["DPO", "TDP", "HMP", "HMFAD"],
        ["DDP", "NF", "PTR", "FUR", "HA", "MP", "AF", "MGO", "DMHF"],
        ["MFT", "FFT", "MP3P", "MP2P", "ACTZ", "DMHFS"],
        ["MFTD", "FFTD", "MMFT", "BND", "BND_F", "SINK"],
    ],
    TRUNK: [
        ["Glc", "Gly"], ["SB", "Fru"], ["AMA"], ["TDG", "ODG", "DDG", "G"],
        ["HMF", "DMHF", "MGO", "GO", "DA", "AF", "FA", "AA", "HA"], ["MEL_C", "MEL", "MEL_N", "BROWN", "SINK"],
    ],
    ACRYLAMIDE: [["Asn", "Glc", "Cys", "Gln", "Lys", "Ala"], ["SBA"], ["INT1", "Asp"], ["ACR"], ["ACRCYS", "SINK"]],
}
LABELS = {
    "PENT": "ribose / xylose", "Cys": "cysteine", "THI": "thiamine", "Glc": "glucose", "Gly": "glycine", "TTCA": "ring intermediate (TTCA)",
    "ARP": "pentose Amadori compound", "H2S": "hydrogen sulfide", "MESH": "methanethiol", "DPO": "1-deoxypentosone", "TDP": "3-deoxypentosone",
    "HMP": "5-hydroxy-3-mercapto-2-pentanone", "DDP": "1,4-dideoxypentosone", "NF": "norfuraneol", "PTR": "2,3,4-pentanetrione",
    "FUR": "furfural", "HA": "hydroxyacetaldehyde", "MP": "1-mercapto-2-propanone", "MFT": "MFT (meaty thiol)", "FFT": "FFT (roasted thiol)",
    "MP3P": "2-mercapto-3-pentanone", "MP2P": "3-mercapto-2-pentanone", "ACTZ": "2-acetylthiazole", "MFTD": "MFT disulfide", "FFTD": "FFT disulfide",
    "MMFT": "MFT-methanethiol adduct", "BND": "MFT bound to the matrix", "BND_F": "FFT bound to the matrix", "HMF": "HMF", "HMFAD": "HMF-cysteine adducts",
    "DMHF": "caramel furanone (DMHF)", "DMHFS": "thio-furanone", "AF": "acetylformoin", "MGO": "methylglyoxal", "GO": "glyoxal", "DA": "diacetyl",
    "G": "glucosone", "TDG": "3-deoxyglucosone", "ODG": "1-deoxyglucosone", "DDG": "3,4-dideoxyglucosone", "SB": "Schiff base", "Fru": "fructose",
    "AMA": "Amadori compound", "FA": "formic acid", "AA": "acetic acid", "MEL": "melanoidins (brown)", "MEL_N": "melanoidins (brown)", "BROWN": "browning",
    "Asn": "asparagine", "SBA": "Schiff base", "INT1": "Maillard intermediate", "Asp": "aspartic acid", "ACR": "acrylamide", "ACRCYS": "acrylamide-cysteine adduct",
    "Gln": "glutamine", "Lys": "lysine", "Ala": "alanine", "MEL_C": "melanoidins (brown)", "CYC": "cyclic intermediate",
    "SINK": "removed into the matrix\n(sink steps)",
}


def staged_positions(G: nx.DiGraph, lane: str) -> Dict[str, Tuple[float, float]]:
    stages = [list(s) for s in STAGES[lane]]
    placed = {n for s in stages for n in s}
    # anything not curated goes to the stage after its last placed predecessor
    for n in list(G.nodes):
        if n in placed:
            continue
        preds = [p for p in G.predecessors(n) if p in placed]
        idx = max((i for i, s in enumerate(stages) for p in preds if p in s), default=0) + 1
        idx = min(idx, len(stages) - 1)
        stages[idx].append(n)
        placed.add(n)
    pos = {}
    tallest = max(len(s) for s in stages)
    for i, s in enumerate(stages):
        s = [n for n in s if n in G.nodes]
        for j, n in enumerate(s):
            pos[n] = (i * 2.6, -(j + 0.5) * tallest / max(len(s), 1))
    return pos


def short_label(sp) -> str:
    lab = sp.label.split(" (")[0]
    return lab if len(lab) <= 26 else lab[:24] + "…"


def draw_lane(reactions, species_list, lane: str, title: str, fname: str, figsize=(16, 10)) -> Dict[str, int]:
    status = constant_status(lane)
    quality = species_quality()
    labels = {s.key: LABELS.get(s.key, short_label(s)) for s in species_list}
    labels.update({k: v for k, v in LABELS.items() if k not in labels})
    G = nx.DiGraph()
    edges = []
    for r in reactions:
        none = list(STATUS_STYLE)[4]
        st = status.get(r.parameter_key, none) if r.parameter_key else none
        if st == "inert":
            continue
        visible_products = [b for b in r.products if b not in HIDDEN]
        for a in r.reactants:
            if a in HIDDEN:
                continue
            if not visible_products:
                # a removal step: everything it makes is a bookkeeping pool -> draw it to the sink box
                G.add_edge(a, "SINK")
                edges.append((a, "SINK", st, r.key))
                continue
            for b in visible_products:
                if a == b:
                    continue
                G.add_edge(a, b)
                edges.append((a, b, st, r.key))
    pos = staged_positions(G, lane)
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 1.3, max(xs) + 1.3)
    ax.set_ylim(min(ys) - 2.6, max(ys) + 0.8)
    # edges first (dedupe by pair, keep the "best-known" status for drawing)
    rank = list(STATUS_STYLE)
    best: Dict[Tuple[str, str], str] = {}
    for a, b, st, _ in edges:
        if (a, b) not in best or rank.index(st) < rank.index(best[(a, b)]):
            best[(a, b)] = st
    for (a, b), st in best.items():
        colour, ls, lw = STATUS_STYLE[st]
        dx = pos[b][0] - pos[a][0]
        rad = 0.18 if abs(dx) < 1e-9 else (0.08 if dx > 3 else 0.0)
        if dx < 0:
            rad = 0.25   # a backward step is drawn as a visible loop
        ax.add_patch(FancyArrowPatch(pos[a], pos[b], arrowstyle="-|>", mutation_scale=10, color=colour, lw=lw, linestyle=ls,
                                     connectionstyle=f"arc3,rad={rad}", shrinkA=16, shrinkB=16, alpha=0.85, zorder=1))
    bw, bh = 1.9, 0.5
    for n, (x, y) in pos.items():
        q = quality.get(n, "none")
        ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    fc=NODE_FILL[q], ec=NODE_EDGE[q], lw=1.3 if q != "none" else 0.8, zorder=2))
        ax.text(x, y, labels.get(n, n), ha="center", va="center", fontsize=8.2, color=INK, zorder=3)
    counts = defaultdict(int)
    for _, _, st, key in {(a, b, st, key) for a, b, st, key in edges}:
        pass
    seen_keys = set()
    for a, b, st, key in edges:
        if key not in seen_keys:
            seen_keys.add(key)
            counts[st] += 1
    handles = [plt.Line2D([0], [0], color=c, ls=ls, lw=lw, label=f"{s} ({counts.get(s, 0)})") for s, (c, ls, lw) in STATUS_STYLE.items()]
    handles += [plt.Rectangle((0, 0), 1, 1, fc=NODE_FILL[q], ec=NODE_EDGE[q], label=t)
                for q, t in (("good", "predicted within 3x where measured"), ("mid", "3-10x off"), ("bad", "more than 10x off"), ("none", "not measured on the panel"))]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False, ncol=2, title="arrows: how the step's rate constant is known · boxes: how well the panel's measurements of the species are predicted", title_fontsize=8)
    ax.set_title(title, loc="left", fontsize=12, color=INK)
    fig.tight_layout()
    fig.savefig(OUT / fname, bbox_inches="tight")
    plt.close(fig)
    return dict(counts)


def fig_summary(counts_by_lane: Dict[str, Dict[str, int]]) -> None:
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    lanes = list(counts_by_lane)
    left = [0.0] * len(lanes)
    for st, (colour, _, _) in STATUS_STYLE.items():
        vals = [counts_by_lane[l].get(st, 0) for l in lanes]
        ax.barh(lanes, vals, left=left, color=colour, label=st, height=0.6)
        for i, v in enumerate(vals):
            if v:
                ax.text(left[i] + v / 2, i, str(v), ha="center", va="center", fontsize=8, color="white" if st != "no constant" else INK)
        left = [a + b for a, b in zip(left, vals)]
    ax.invert_yaxis()
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("reaction steps")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8, frameon=False, ncol=2)
    ax.set_title("How each path's steps are known (same colours as the chemistry figure)", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "13_steps_by_status.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {}
    counts["sugar and amino acid"] = draw_lane(network.TRUNK_REACTIONS, species.SPECIES, TRUNK,
                                               "The sugar path: glucose and glycine to browning, furans and the small dicarbonyls",
                                               "10_tree_sugar.png", figsize=(16, 9))
    counts["pentose and cysteine"] = draw_lane(sul_mod.SULFUR_REACTIONS, species_sulfur.SULFUR_SPECIES, SULFUR,
                                               "The pentose-cysteine path: from ribose or xylose and cysteine to the meaty thiols and their removal",
                                               "11_tree_sulfur.png", figsize=(19, 11))
    counts["asparagine and glucose"] = draw_lane(acr_mod.ACRYLAMIDE_REACTIONS, species_acrylamide.ACRYLAMIDE_SPECIES, ACRYLAMIDE,
                                                 "The acrylamide path: asparagine and glucose to acrylamide and its elimination",
                                                 "12_tree_acrylamide.png", figsize=(13, 7))
    fig_summary(counts)
    from figure_manifest import ENGINE_SOURCES, record

    record("scripts/generators/build_reaction_tree.py",
           [V / "core_panel_scores.json", V / "core_prediction_uncertainty.json", V / "kinetic_core_b9_fit_report.json",
            V / "kinetic_core_b9_laplace_covariance.json", ROOT / "scripts" / "generators" / "build_reaction_tree.py", *ENGINE_SOURCES],
           ["10_tree_sugar.png", "11_tree_sulfur.png", "12_tree_acrylamide.png", "13_steps_by_status.png"])
    print("wrote 3 trees + summary:", {k: dict(v) for k, v in counts.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
