#!/usr/bin/env python
"""
Figures for docs/guides/INTRODUCTION.md and REACTION_TREES.md (2026-09-07).

Six plots of MEASURED against MODEL for the sulfur lane's thiol-sink diagnosis. Model values are
read from the frozen artifacts (the B16 ship rule, the directional scorecard); measured values are
literals with their source anchors, re-typed from the extraction dossiers named beside them.

    python scripts/generators/build_thiol_sink_figures.py      # writes docs/assets/thiol_sink/*.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from src import data_paths  # noqa: E402

OUT = ROOT / "docs" / "assets" / "thiol_sink"
V = data_paths.VALIDATION_DIR

MEASURED = "#2B5DA8"     # what the pot measured
SHIPPED = "#178F6E"      # the calibration the tool ships today (wave B9)
RETUNED = "#D9822B"      # the 2026-09-07 re-tuning attempt, sink barrier ceiling kept (wave B16)
LIFTED = "#B5468A"       # the same attempt with the ceiling lifted (wave B16, lift variant)
INK = "#1E2A2C"
MUTED = "#5E6B6E"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#D6DBD8", "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK, "axes.titleweight": "bold",
    "axes.titlesize": 11.5, "legend.frameon": False, "legend.fontsize": 9, "figure.dpi": 150,
})


def _read(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def _style(ax, ylabel: str, xlabel: str = "") -> None:
    ax.grid(True, axis="y", color="#E7EBE9", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)


def fig_schieberle(ship) -> None:
    """The fit's own pot at 100 C: Schieberle, Hofmann & Muench 2000 Table IV (SIDA), ug per 100 mL -> ug/L."""
    t = [30, 60, 360, 720]
    meas = {"MFT": [45, 138, 1560, 1790], "FFT": [20, 31, 1100, 1320]}     # schieberle2000_extraction.md sec. 2
    scores = _read(V / "core_directional_scores.json")
    sch = next(c for c in scores["claims"] if c["claim_id"] == "SCH-T-01")
    shipped_mft = sch["values_ug_per_l"]                                     # the shipped lane, scored on the panel
    b16 = ship["b16"]["T1_shape"]
    lift = ship["b16_lift"]["T1_shape"] if ship.get("b16_lift") else None
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
    for ax, sp in zip(axes, ("MFT", "FFT")):
        ax.plot(t, meas[sp], "-o", color=MEASURED, lw=2, ms=5, label="measured in the pot")
        if sp == "MFT":
            ax.plot(t, shipped_mft, "-o", color=SHIPPED, lw=2, ms=5, label="model as shipped")
        ax.plot(t, b16["mft_ug_per_l" if sp == "MFT" else "fft_ug_per_l"], "-o", color=RETUNED, lw=2, ms=5,
                label="model re-tuned on this series")
        if lift:
            ax.plot(t, lift["mft_ug_per_l" if sp == "MFT" else "fft_ug_per_l"], "--o", color=LIFTED, lw=2, ms=5,
                    label="re-tuned, sink barrier ceiling lifted")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(t)
        ax.set_xticklabels(["30 min", "1 h", "6 h", "12 h"])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(f"{'2-methyl-3-furanthiol (MFT)' if sp == 'MFT' else '2-furfurylthiol (FFT)'}")
        _style(ax, "µg per litre (log scale)" if sp == "MFT" else "", "time at 100 °C")
    axes[0].legend(loc="lower right")
    fig.suptitle("Ribose + cysteine, phosphate buffer pH 5, held at 100 °C", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "01_hofmann_pot_100C.png")
    plt.close(fig)


def fig_wang(scores) -> None:
    """Wang 2022 (FFJ): the model's own curves; the pot's shape as the paper states it."""
    by = {c["claim_id"]: c for c in scores["claims"]}
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot([30, 90, 180], by["WANG22-T-01"]["values_ug_per_l"], "-o", color=SHIPPED, lw=2, ms=5, label="model, 100 °C")
    ax.plot([30, 60, 120, 180], by["WANG22-T-03"]["values_ug_per_l"], "--o", color=SHIPPED, lw=2, ms=5, label="model, 140 °C")
    ax.set_yscale("log")
    ax.set_ylim(1, 20000)
    ax.annotate("what the pot does at 140 °C: rises to 60 min,\nthen declines gently to 180 min",
                xy=(60, 693), xytext=(75, 6000), color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.annotate("what the pot does at 100 °C:\nstill rising at 180 min", xy=(180, 391), xytext=(110, 40), color=INK, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_xticks([30, 60, 90, 120, 180])
    _style(ax, "2-methyl-3-furanthiol, µg per litre (log)", "minutes")
    ax.set_title("Cysteine + xylose in buffer (Wang 2022): the model at 100 and 140 °C")
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "02_wang2022_shapes.png")
    plt.close(fig)


def fig_liu(ship) -> None:
    """Liu 2023 (LWT) at 168 C, indexed to the 20-min point."""
    t = [20, 40, 60]
    meas_mft = [1253, 760, 628]    # liu2023b_extraction.md sec. 2, ng per vial
    meas_fft = [175, 198, 170]
    b16 = ship["b16"]["T6_liu2023_168C"]["mft_ug_per_l"]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(t, [v / meas_mft[0] for v in meas_mft], "-o", color=MEASURED, lw=2, ms=5, label="measured MFT")
    ax.plot(t, [v / meas_fft[0] for v in meas_fft], "--o", color=MEASURED, lw=2, ms=5, label="measured FFT")
    ax.plot(t, [v / b16[0] for v in b16], "-o", color=RETUNED, lw=2, ms=5, label="re-tuned model, MFT")
    ax.axhline(1.0, color="#D6DBD8", lw=0.8)
    ax.set_ylim(0.3, 1.3)
    ax.set_xticks(t)
    _style(ax, "relative to the 20-minute value", "minutes at 168 °C")
    ax.set_title("Cysteine-rich ribose pot at 168 °C (Liu 2023)")
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "03_liu2023_168C.png")
    plt.close(fig)


def fig_yiltirak(ship) -> None:
    rows9 = ship["b16"]["T3_yiltirak_levels"]["b9"]
    rows16 = ship["b16"]["T3_yiltirak_levels"]["b16"]
    labels = ["MFT\n100 °C, 4 h", "FFT\n100 °C, 4 h", "MFT\n110 °C, 2 h", "FFT\n110 °C, 2 h"]
    meas = [r["measured"] for r in rows9]
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    x = range(len(labels))
    w = 0.26
    ax.bar([i - w for i in x], meas, w, color=MEASURED, label="measured (Reading, 2026)")
    ax.bar(list(x), [r["predicted"] for r in rows9], w, color=SHIPPED, label="model as shipped")
    ax.bar([i + w for i in x], [r["predicted"] for r in rows16], w, color=RETUNED, label="model with the low-temperature sink weakened")
    for i, r in enumerate(rows9):
        ax.text(i, r["predicted"] * 1.15, f"{r['fold']:.0f}×", ha="center", fontsize=8.5, color=INK)
    for i, r in enumerate(rows16):
        ax.text(i + w, r["predicted"] * 1.15, f"{r['fold']:.0f}×", ha="center", fontsize=8.5, color=INK)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 30000)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    _style(ax, "µg per litre (log scale)")
    ax.set_title("Ribose + cysteine at 100 and 110 °C (Reading, 2026): a pot no model version was tuned on")
    ax.legend(loc="upper left", ncol=1)
    fig.tight_layout()
    fig.savefig(OUT / "04_yiltirak_ladder.png")
    plt.close(fig)


def fig_ttca() -> None:
    temps = ["100 °C", "120 °C", "140 °C"]
    meas = [10.331 - 0.0271 * 60, 9.9718 - 0.0651 * 60, 9.3375 - 0.0813 * 60]   # zhai2021_extraction.md sec. 2
    model = [1.58, 0.051, 0.001]                                                   # kinetic_core_b16_prereg.md sec. 6 probe
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    x = range(3)
    w = 0.34
    ax.bar([i - w / 2 for i in x], meas, w, color=MEASURED, label="measured (Zhai 2021)")
    ax.bar([i + w / 2 for i in x], model, w, color=SHIPPED, label="model as shipped")
    for i in x:
        ax.text(i - w / 2, meas[i] + 0.2, f"{meas[i]:.1f}", ha="center", fontsize=8.5, color=INK)
        ax.text(i + w / 2, max(model[i], 0.02) + 0.2, f"{model[i]:.2f}" if model[i] >= 0.01 else "≈0", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(list(x))
    ax.set_xticklabels(temps)
    ax.set_ylim(0, 11)
    _style(ax, "mmol per litre left after 60 min (of 10)")
    ax.set_title("The xylose–cysteine intermediate (TTCA), heated alone for 60 min")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "05_ttca_decay.png")
    plt.close(fig)


def fig_dicarbonyls(scores) -> None:
    dic = next(c for c in scores["claims"] if c["claim_id"] == "DIC-03")
    order = dic["observables"]                     # 3-deoxyglucosone, glucosone, glyoxal, methylglyoxal
    model = dict(zip(order, dic["values_ug_per_l"]))
    measured = {"3-deoxyglucosone": 52.2, "glucosone": 7.5, "glyoxal": 5.6, "methylglyoxal": 2.6}   # leitzen2021_extraction.md, ug/mL
    names = ["glucosone", "glyoxal", "methylglyoxal"]
    m_rel = [measured[n] / measured["3-deoxyglucosone"] for n in names]
    mod_rel = [model[n] / model["3-deoxyglucosone"] for n in names]
    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    x = range(3)
    w = 0.34
    ax.bar([i - w / 2 for i in x], m_rel, w, color=MEASURED, label="measured (Leitzen 2021, 121 °C)")
    ax.bar([i + w / 2 for i in x], mod_rel, w, color=SHIPPED, label="model (sugar lane)")
    ax.axhline(1.0, color="#D6DBD8", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 100)
    ax.set_xticks(list(x))
    ax.set_xticklabels(names)
    _style(ax, "amount relative to 3-deoxyglucosone (log)")
    ax.set_title("Glucose heated alone in water, 121 °C (Leitzen 2021)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "06_dicarbonyls_water.png")
    plt.close(fig)


def fig_map() -> None:
    """The reaction paths the model carries, coloured by how well each is predicted (docs guide, sec. 'The map')."""
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    GOOD, MID, BAD, NONE = ("#D9EFE3", "#178F6E"), ("#FBE9D0", "#D9822B"), ("#F6D9D9", "#B23A3A"), ("#EEEEEE", "#9AA6A3")
    nodes = {
        # key: (x, y, label, status)
        "S": (0.0, 4.2, "sugar +\namino acid", GOOD), "A": (1.6, 4.2, "Amadori\ncompound", GOOD),
        "D": (3.2, 4.2, "deoxyosones,\nsmall dicarbonyls", BAD), "B2": (4.9, 4.7, "HMF (2-12x off)", MID),
        "B": (4.9, 4.05, "brown colour", GOOD), "B3": (4.9, 3.4, "caramel furanone\n(50-270x off)", BAD),
        "P": (0.0, 2.6, "pentose sugar +\ncysteine", MID), "T": (1.6, 2.6, "ring intermediate\n(TTCA)", MID),
        "F": (3.2, 2.6, "furanones, furfural\n+ hydrogen sulfide", MID), "M": (4.9, 2.6, "meaty thiols\nMFT and FFT", BAD),
        "X": (6.5, 2.6, "thiol removal:\ndisulfides, adducts", BAD),
        "H": (0.0, 1.6, "hexose sugar +\ncysteine", NONE),
        "N": (0.0, 0.5, "asparagine +\nglucose", MID), "Y": (1.6, 0.5, "acrylamide", MID), "Z": (3.2, 0.5, "acrylamide\nelimination", MID),
        "L": (4.9, 0.5, "unsaturated fat", MID), "O": (6.5, 0.5, "hydroperoxides,\nhexanal + aldehydes", MID),
    }
    edges = [("S", "A"), ("A", "D"), ("D", "B"), ("D", "B2"), ("D", "B3"), ("P", "T"), ("T", "F"), ("F", "M"), ("M", "X"), ("N", "Y"), ("Y", "Z"), ("L", "O")]
    fig, ax = plt.subplots(figsize=(11, 5.4))
    ax.set_xlim(-0.8, 7.4)
    ax.set_ylim(-0.55, 5.1)
    ax.axis("off")
    bw, bh = 1.25, 0.62
    for key, (x, y, label, (fill, edge)) in nodes.items():
        dashed = key == "H"
        box = FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fill, ec=edge, lw=1.4,
                             ls="--" if dashed else "-")
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center", fontsize=9, color=INK if not dashed else MUTED)
    for a, b in edges:
        xa, ya = nodes[a][0], nodes[a][1]
        xb, yb = nodes[b][0], nodes[b][1]
        ax.add_patch(FancyArrowPatch((xa + bw / 2, ya), (xb - bw / 2, yb), arrowstyle="-|>", mutation_scale=12, color=MUTED, lw=1.2))
    ax.add_patch(FancyArrowPatch((nodes["H"][0] + bw / 2, nodes["H"][1]), (nodes["M"][0] - 0.1, nodes["M"][1] - bh / 2),
                                 arrowstyle="-|>", mutation_scale=12, color=MUTED, lw=1.2, ls="--", connectionstyle="arc3,rad=0.25"))
    ax.text(1.35, 1.0, "no route in the model", fontsize=8.5, color=MUTED, style="italic")
    legend = [("predicts held-out data within about 1.5x", GOOD), ("right shape; within 3x inside the source lab only", MID),
              ("wrong by 10x or more, or wrong in direction", BAD), ("no route exists", NONE)]
    # legend: two rows under the diagram
    for i, (txt, (fill, edge)) in enumerate(legend):
        x0 = -0.6 + (i % 2) * 3.9
        y0 = -0.12 - (i // 2) * 0.34
        ax.add_patch(FancyBboxPatch((x0, y0), 0.22, 0.18, boxstyle="round,pad=0.01", fc=fill, ec=edge, lw=1.2))
        ax.text(x0 + 0.3, y0 + 0.09, txt, fontsize=8.5, color=INK, va="center")
    ax.set_title("The reaction paths the model carries, coloured by how well each one is predicted", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "00_map.png")
    plt.close(fig)


#: The papers whose MEASUREMENTS became rate constants or fit rows, by path (hand-curated from the
#: parameter registries and the frozen generators' row anchors, 2026-09-07).
CONSTANT_SOURCES = {
    "sugar and amino acid": ["Martins & van Boekel 2005", "Martins & van Boekel 2003", "Kocadagli & Gokmen 2016", "Pereyra Gonzales 2010",
                             "Bell 1995", "Shu 1988", "Hamzalioglu 2018", "Poisson 2019", "Wang 2008"],
    "pentose and cysteine": ["Hofmann & Schieberle 1998", "Hofmann & Schieberle 2002", "Kumazawa 2003", "Cerny 2007", "Whitfield 1999",
                             "van Seeventer 2001", "Yaghmur 2005", "Zhou 2023", "Zhang 2024", "Kang 2026", "Feng 2022", "Zhai 2023",
                             "Charles-Bernard 2005", "Gigl 2021"],
    "asparagine and glucose": ["De Vleeschouwer 2006", "De Vleeschouwer 2007", "De Vleeschouwer 2008", "De Vleeschouwer 2009 I",
                               "De Vleeschouwer 2009 II", "Claeys 2005", "Knol 2005", "Knol 2009", "Knol 2010"],
    "fat oxidation": ["Frankel 1989", "Schroen 2022"],
}
PATH_COLOURS = {"sugar and amino acid": "#178F6E", "pentose and cysteine": "#B23A3A", "asparagine and glucose": "#D9822B", "fat oxidation": "#8A6BBF"}


def fig_funnel() -> None:
    """How much of the registered literature is inside the model, and in what role."""
    import yaml

    reg = yaml.safe_load((ROOT / "data" / "keys" / "papers.yml").read_text(encoding="utf-8"))["papers"]
    n_reg = len(reg)
    n_intake = sum(1 for p in reg if p["intake_ids"])
    n_dossier = sum(1 for p in reg if p["dossier"])
    n_bench = sum(1 for p in reg if any(f.startswith("data/benchmarks") for f in p["record_ids"]))
    panel = yaml.safe_load((ROOT / "docs" / "validation" / "directional_claims_panel.yml").read_text(encoding="utf-8"))
    n_claims_src = len(panel["sources"])
    n_const = sum(len(v) for v in CONSTANT_SOURCES.values())
    stages = [("registered in the corpus", n_reg, MUTED), ("screened and indexed", n_intake, MUTED), ("read in full (extraction dossier)", n_dossier, MUTED),
              ("supply a claim about direction (87 claims)", n_claims_src, "#4F80D0"), ("supply a validation measurement (a level)", n_bench, "#2B5DA8"),
              ("supply a rate constant or a fit row", n_const, "#1E2A2C")]
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    ys = list(range(len(stages)))[::-1]
    for y, (label, n, colour) in zip(ys, stages):
        if label.startswith("supply a rate"):
            x = 0
            for path, names in CONSTANT_SOURCES.items():
                ax.barh(y, len(names), left=x, color=PATH_COLOURS[path], height=0.62)
                ax.text(x + len(names) / 2, y, str(len(names)), ha="center", va="center", fontsize=8.5, color="white")
                x += len(names)
        else:
            ax.barh(y, n, color=colour, height=0.62)
        ax.text(n + 3, y, f"{n}", va="center", fontsize=10, color=INK, fontweight="bold")
        ax.text(-4, y, label, va="center", ha="right", fontsize=9.5, color=INK)
    ax.set_yticks([])
    ax.set_xlim(0, 320)
    ax.set_xlabel("papers")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis="x", color="#E7EBE9", linewidth=0.8)
    ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in PATH_COLOURS.values()]
    ax.legend(handles, [f"{k} ({len(v)})" for k, v in CONSTANT_SOURCES.items()], loc="lower right", title="rate constants, by path", fontsize=8.5, title_fontsize=8.5)
    ax.set_title("How much of the literature is inside the model, and in what role", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "07_literature_funnel.png", bbox_inches="tight")
    plt.close(fig)


def fig_scorecard() -> None:
    """The path-by-path scorecard as an image: what we have, how it does, what we lack."""
    rows = [
        ("sugar + amino acid\n-> brown colour", "9", "one glucose-glycine study at 3 temperatures;\nwater-activity and pH ratios",
         "browning within 1.5x (held out); HMF 2-12x;\ncaramel furanone 50-270x off", "small dicarbonyls: constants from a\nsugar glass, wrong order in water", "#FBE9D0"),
        ("pentose + cysteine\n-> meaty thiols", "14", "every step at 145 C from one lab's\nfed-intermediate experiments",
         "2-7x in that lab at 145 C and pH 5;\n20-140x at pH 3 or 7; 10-500x elsewhere", "how fast a thiol is REMOVED, at more\nthan one temperature; pH on formation", "#F6D9D9"),
        ("hexose + cysteine\n-> meaty thiols", "0", "nothing at step level", "declares 'unknown' (no route)", "the furfural / furfuryl-alcohol route;\none 168 C time series waits as its test", "#EEEEEE"),
        ("asparagine + glucose\n-> acrylamide", "9", "formation, elimination, pH and\nwater-activity effects, one lab, 120-200 C",
         "other labs' 180 C pots 2.5-220x (median 9x);\nextrusion in real food 10,000x", "a second laboratory's constants;\nreal-food matrices", "#FBE9D0"),
        ("unsaturated fat\n-> aldehydes", "2", "six products and their split\nfrom one 1989 study", "cooked rows 4-34x; 40 C storage rows\nnot comparable (model starts from zero)",
         "nonanal, 2-pentylfuran (no branch);\na storage baseline", "#FBE9D0"),
    ]
    cols = ["path", "papers behind\nits constants", "what we have", "how it does", "what we lack"]
    fig, ax = plt.subplots(figsize=(12.5, 5.2))
    ax.axis("off")
    widths = [0.15, 0.09, 0.28, 0.24, 0.24]
    x0 = [sum(widths[:i]) for i in range(len(widths))]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(rows) + 1)
    for i, c in enumerate(cols):
        ax.text(x0[i] + 0.01, len(rows) + 0.5, c, fontsize=9.5, color=MUTED, va="center", fontweight="bold")
    ax.plot([0, 1], [len(rows) + 0.1, len(rows) + 0.1], color="#1E2A2C", lw=1)
    for r, row in enumerate(rows):
        y = len(rows) - r - 0.5
        ax.add_patch(plt.Rectangle((0, y - 0.5), 1, 1, color=row[-1], alpha=0.55, lw=0))
        for i, cell in enumerate(row[:-1]):
            ax.text(x0[i] + 0.01, y, cell, fontsize=9.2 if i else 9.8, color=INK, va="center", fontweight="bold" if i == 0 else "normal",
                    ha="center" if i == 1 else "left", transform=ax.transData if i != 1 else ax.transData)
            if i == 1:
                ax.texts[-1].set_x(x0[1] + widths[1] / 2)
        ax.plot([0, 1], [y - 0.5, y - 0.5], color="#D6DBD8", lw=0.8)
    ax.set_title("Path by path: what the model has, how it does, what it lacks", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "08_path_scorecard.png", bbox_inches="tight")
    plt.close(fig)


def fig_repo_flow() -> None:
    """What the repository is made of, as a flow from papers to guides."""
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    boxes = [
        (0.0, "published papers", "PDFs in data/articles\n297 registered"),
        (1.0, "extraction dossiers", "data/lit/extraction_dossiers\n58 papers read in full,\nevery table re-typed"),
        (2.0, "three kinds of evidence", "rate constants (34 papers)\nbenchmark pots (53 files)\ndirectional claims (92)"),
        (3.0, "the kinetic model", "src/kinetic_core\n4 paths, 159 steps;\nevery calibration\npre-registered"),
        (4.0, "scorecards", "results/validation\nlevels, directions,\nintervals, wishlist"),
        (5.0, "guides and the tool", "docs, maillard.py\nthis guide; predict,\ncompare, rank, score"),
    ]
    fig, ax = plt.subplots(figsize=(14, 3.9))
    ax.set_xlim(-0.55, 5.55)
    ax.set_ylim(-0.95, 0.95)
    ax.axis("off")
    for x, head, sub in boxes:
        ax.add_patch(FancyBboxPatch((x - 0.46, -0.72), 0.92, 1.44, boxstyle="round,pad=0.02,rounding_size=0.06", fc="#F2F3F1", ec="#9AA6A3", lw=1.2))
        ax.text(x, 0.45, head, ha="center", va="center", fontsize=9.5, color=INK, fontweight="bold")
        ax.text(x, -0.12, sub, ha="center", va="center", fontsize=8, color=MUTED, linespacing=1.4)
    for x in range(5):
        ax.add_patch(FancyArrowPatch((x + 0.47, 0.0), (x + 0.53, 0.0), arrowstyle="-|>", mutation_scale=12, color=MUTED, lw=1.2))
    ax.set_title("What the repository is made of: from papers to predictions", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "09_repository_flow.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# What the FIELD knows (docs/guides/INTRODUCTION.md, sec. 1-2): the accepted scheme, annotated by
# how well each part has been measured in the published literature the repository has read.
# ---------------------------------------------------------------------------
FIELD_STATUS = {   # the same scale the reaction trees use (build_reaction_tree.STATUS_STYLE)
    "rate known at several temperatures": ("#2B5DA8", "-", 2.4),
    "rate or yield known at one temperature": ("#178F6E", "-", 2.0),
    "mechanism known, no rate (labelling, products)": ("#9AA6A3", "-", 1.6),
    "open: one measurement at 121 °C, no temperature dependence": ("#B23A3A", "--", 2.0),
}


def fig_field_scheme() -> None:
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    S = list(FIELD_STATUS)
    nodes = {
        "sug": (0.0, 5.0, "reducing sugar\n+ amino acid"), "ama": (1.7, 5.0, "Amadori / Heyns\ncompound"),
        "dox": (3.4, 5.6, "3-deoxyosone\n(1,2-enolisation)"), "dox2": (3.4, 4.4, "1-deoxyosone\n(2,3-enolisation)"),
        "frag": (3.4, 3.2, "sugar fragments:\nglyoxal, methylglyoxal,\ndiacetyl"),
        "hmf": (5.1, 5.6, "HMF, furfural"), "fur": (5.1, 4.4, "furanones\n(caramel, norfuraneol)"),
        "str": (6.8, 3.2, "Strecker aldehydes\n+ amino-ketones"), "mel": (6.8, 5.0, "melanoidins\n(brown colour)"),
        "pyr": (8.5, 3.2, "pyrazines, pyrroles"),
        "cys": (0.0, 1.6, "cysteine"), "h2s": (1.7, 1.6, "H2S, NH3,\nacetaldehyde"),
        "thiol": (5.1, 1.6, "meaty thiols\nMFT, FFT"), "sink": (6.8, 1.6, "disulfides, adducts,\npolymers"),
        "ttca": (3.4, 1.6, "ring intermediate\n(TTCA), deoxypentosones"),
        "asn": (0.0, -0.5, "asparagine\n+ sugar"), "acr": (3.4, -0.5, "acrylamide"), "acr2": (5.1, -0.5, "acrylamide\nelimination"),
        "lip": (0.0, -1.8, "unsaturated fat"), "ald": (3.4, -1.8, "hydroperoxides,\naldehydes"), "lm": (5.1, -1.8, "lipid-Maillard:\nalkylthiophenes"),
    }
    edges = [   # (from, to, status, label, rad, label dy)
        ("sug", "ama", S[0], "Martins 2005, 3 T", 0.0, 0.48), ("ama", "dox", S[0], "Martins 2005", 0.12, 0.0), ("ama", "dox2", S[0], "", 0.12, 0.0),
        ("ama", "frag", S[0], "Martins 2005; Kocadagli 2016 (glass)", 0.12, 0.0),
        ("dox", "hmf", S[0], "Martins 2005; Kocadagli 2016", 0.0, 0.48), ("dox2", "fur", S[0], "Kocadagli 2016: 160-200 C, glass only", 0.0, 0.48),
        ("frag", "str", S[1], "Hofmann 2000: yields at 98 C", 0.0, 0.48, -0.6), ("dox2", "frag", S[0], "", 0.12, 0.0),
        ("frag", "mel", S[1], "Martins 2005: from methylglyoxal, lumped", 0.18, 0.55, 1.2), ("str", "pyr", S[2], "mechanism", 0.0, 0.48),
        ("cys", "h2s", S[1], "Hofmann 1998, 145 C", 0.0, 0.48), ("cys", "ttca", S[0], "Zhai 2021: 100-140 C, zero order", 0.35, -0.72, -0.9),
        ("ttca", "thiol", S[1], "Kang 2026, Zhai 2023: 100-140 C levels", 0.0, 0.48),
        ("h2s", "thiol", S[1], "Whitfield 1999 / 2001: 140 C, pH 4.5 and 6.5", 0.32, -1.05, 0.6),
        ("fur", "thiol", S[1], "norfuraneol + H2S, 145 C", 0.12, -0.55, 0.95), ("hmf", "thiol", S[1], "furfural + H2S, 145 C", 0.12, 0.35, 0.95),
        ("thiol", "sink", S[3], "Kumazawa 2003: 121 C only; disulfides never quantified", 0.0, 0.48), ("frag", "h2s", S[2], "Strecker of cysteine (mechanism)", -0.25, -0.35, -0.5),
        ("asn", "acr", S[0], "De Vleeschouwer 2006-09, Knol: 120-200 C, pH, a_w", 0.0, 0.48), ("acr", "acr2", S[0], "same series", 0.0, 0.48),
        ("lip", "ald", S[1], "Frankel 1989", 0.0, 0.48), ("ald", "lm", S[2], "Wang 2022, products", 0.0, 0.48),
    ]
    fig, ax = plt.subplots(figsize=(15, 8.6))
    ax.set_xlim(-0.9, 9.5)
    ax.set_ylim(-3.2, 6.3)
    ax.axis("off")
    bw, bh = 1.35, 0.72
    for a, b, st, lab, rad, dy, *dx in edges:
        dx = dx[0] if dx else 0.0
        colour, ls, lw = FIELD_STATUS[st]
        (xa, ya, _), (xb, yb, _) = nodes[a], nodes[b]
        ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb), arrowstyle="-|>", mutation_scale=12, color=colour, lw=lw, linestyle=ls,
                                     connectionstyle=f"arc3,rad={rad}", shrinkA=26, shrinkB=26, zorder=1))
        if lab:
            xm, ym = (xa + xb) / 2 + dx, (ya + yb) / 2 + dy
            ax.text(xm, ym, lab, fontsize=6.8, color=colour, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=4)
    for key, (x, y, label) in nodes.items():
        ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.08", fc="#F2F3F1", ec="#9AA6A3", lw=1.1, zorder=2))
        ax.text(x, y, label, ha="center", va="center", fontsize=8.6, color=INK, zorder=3)
    for x, y, txt in ((-0.85, 6.05, "SUGAR AND AMINO ACID (the Hodge scheme)"), (-0.85, 2.25, "SULFUR: CYSTEINE AND A PENTOSE"),
                      (-0.85, 0.15, "ASPARAGINE"), (-0.85, -1.15, "FAT")):
        ax.text(x, y, txt, fontsize=8.5, color=MUTED, fontweight="bold", ha="left")
    handles = [plt.Line2D([0], [0], color=c, ls=ls, lw=lw, label=s) for s, (c, ls, lw) in FIELD_STATUS.items()]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.0), fontsize=8.5, frameon=False, ncol=2, title="how well the published literature has measured the step", title_fontsize=8.5)
    ax.set_title("The Maillard reaction as the field draws it, and how well each part has been measured", loc="left", fontsize=11.5)
    fig.tight_layout()
    fig.savefig(OUT / "14_field_scheme.png", bbox_inches="tight")
    plt.close(fig)


def fig_field_coverage() -> None:
    """Where the quantitative measurements the repository has read actually sit: path x temperature, by matrix."""
    rows = [
        # (path, T_lo, T_hi, matrix, label)
        ("sugar + amino acid", 100, 120, "water", "Martins & van Boekel 2005 (glucose-glycine, rates)"),
        ("sugar + amino acid", 120, 120, "water", "Brands & van Boekel 2001 (sugar-casein)"),
        ("sugar + amino acid", 111, 121, "water", "Leitzen 2021 (glucose alone, dicarbonyls)"),
        ("sugar + amino acid", 90, 110, "water", "Zhang 2021 (glucose-glutamate, dicarbonyls)"),
        ("sugar + amino acid", 160, 200, "dry glass", "Kocadagli & Gokmen 2016 (glucose glass)"),
        ("sugar + amino acid", 150, 170, "dry food", "Goncuoglu Tas 2016 (hazelnut)"),
        ("sugar + amino acid", 37, 60, "powder", "Pereyra Gonzales 2010 (milk powder, water activity)"),
        ("sugar + amino acid", 98, 100, "water", "Hofmann 2000 / 2000b (Strecker yields, air vs argon)"),
        ("pentose + cysteine", 145, 145, "water", "Hofmann & Schieberle 1998 (fed intermediates)"),
        ("pentose + cysteine", 140, 140, "water", "Whitfield & Mottram 1999 / 2001 (fed norfuraneol, pH 4.5 / 6.5)"),
        ("pentose + cysteine", 100, 100, "water", "Schieberle 2000 (time series)"),
        ("pentose + cysteine", 100, 140, "water", "Zhai 2021 / 2023, Kang 2026 (TTCA ladders)"),
        ("pentose + cysteine", 100, 140, "water", "Wang 2022 (five-temperature grid, figures only)"),
        ("pentose + cysteine", 100, 130, "water", "Yiltirak 2026 (ladder, stated vessel)"),
        ("pentose + cysteine", 121, 121, "water", "Kumazawa 2003 (thiol loss grid)"),
        ("pentose + cysteine", 80, 80, "brew", "Hofmann 2002 (FFT loss in coffee brew)"),
        ("pentose + cysteine", 168, 168, "water", "Liu 2023 (decline)"),
        ("asparagine + glucose", 120, 200, "water / powder", "De Vleeschouwer 2006-2009, Knol 2005-2010, Claeys 2005"),
        ("fat", 25, 100, "oil", "Frankel 1989, Schroen 2022"),
    ]
    paths = ["sugar + amino acid", "pentose + cysteine", "asparagine + glucose", "fat"]
    colours = {"water": "#2B5DA8", "dry glass": "#D9822B", "dry food": "#B5468A", "powder": "#8A6BBF", "water / powder": "#178F6E", "brew": "#5E6B6E", "oil": "#9AA6A3"}
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ypos = {p: i for i, p in enumerate(paths)}
    n_of = {p: sum(1 for r in rows if r[0] == p) for p in paths}
    stack = {p: 0 for p in paths}
    for path, lo, hi, matrix, label in rows:
        y = ypos[path] + (stack[path] - (n_of[path] - 1) / 2) * 0.1
        stack[path] += 1
        ax.plot([lo, hi], [y, y], "-", color=colours[matrix], lw=3.5, solid_capstyle="round", alpha=0.9)
        if lo == hi:
            ax.plot([lo], [y], "o", color=colours[matrix], ms=6)
        ax.text(hi + 3, y, label, fontsize=7.4, va="center", color=INK)
    ax.set_yticks(range(len(paths)))
    ax.set_yticklabels(paths)
    ax.invert_yaxis()
    ax.set_xlim(20, 320)
    ax.set_xlabel("temperature, °C (each bar: the temperatures one study measured at)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis="x", color="#E7EBE9", linewidth=0.8)
    ax.axvspan(100, 145, color="#F3E9D2", alpha=0.5, lw=0, zorder=0)
    ax.set_ylim(len(paths) - 0.4, -0.95)
    ax.text(122, -0.78, "cooking window with most data", ha="center", fontsize=8, color=MUTED)
    handles = [plt.Line2D([0], [0], color=c, lw=3.5, label=m) for m, c in colours.items()]
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.36), fontsize=8, frameon=False, title="matrix", title_fontsize=8)
    ax.set_title("Where the quantitative measurements sit: temperature and matrix, by path", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "15_field_coverage.png", bbox_inches="tight")
    plt.close(fig)


def fig_how_a_model_works() -> None:
    """The two ideas a kinetic model rests on: steps give time courses; rate constants rise with temperature."""
    import numpy as np

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.9))
    t = np.linspace(0, 120, 400)
    k1, k2 = 0.05, 0.02                          # per minute: A -> B -> C
    A = np.exp(-k1 * t)
    B = k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    C = 1 - A - B
    a1.plot(t, A, color="#9AA6A3", lw=2, label="starting material")
    a1.plot(t, B, color="#2B5DA8", lw=2.4, label="aroma compound (formed, then removed)")
    a1.plot(t, C, color="#B23A3A", lw=2, label="removal product")
    a1.set_xlabel("minutes at one temperature")
    a1.set_ylabel("fraction of the starting amount")
    a1.set_title("Two steps in a row: formation, then removal")
    a1.legend(loc="center right", fontsize=8)
    _style(a1, "fraction of the starting amount", "minutes at one temperature")
    T = np.linspace(80, 180, 200)
    R = 8.314e-3
    for ea, colour, lab in ((60, "#178F6E", "activation energy 60 kJ/mol"), (100, "#D9822B", "100 kJ/mol"), (160, "#B23A3A", "160 kJ/mol")):
        k = np.exp(-ea / R * (1 / (T + 273.15) - 1 / (145 + 273.15)))
        a2.plot(T, k, color=colour, lw=2.2, label=lab)
    a2.axvline(145, color="#D6DBD8", lw=1)
    a2.text(146, 0.02, "measured here", fontsize=8, color=MUTED)
    a2.set_yscale("log")
    a2.set_ylim(0.005, 20)
    a2.set_title("The same step at other temperatures, relative to 145 °C")
    a2.legend(loc="upper left", fontsize=8)
    _style(a2, "rate relative to the rate at 145 °C (log)", "temperature, °C")
    fig.tight_layout()
    fig.savefig(OUT / "16_how_a_kinetic_model_works.png", bbox_inches="tight")
    plt.close(fig)


def fig_papers_weight() -> None:
    """Which papers the model takes the most from: fit rows, benchmark pots and directional claims per paper,
    counted from the frozen generators' row anchors, the paper registry and the claims panel."""
    import re
    from collections import Counter

    import yaml

    pat = re.compile(r"\b([A-Z][a-zA-Z\-]+)(?: et al\.?| & (?:van |de |De )?[A-Z][a-zA-Z]+| and [A-Z][a-zA-Z]+)? ((?:19|20)\d{2})\b")
    rows: Counter = Counter()
    seen_anchors = set()       # the wave generators re-import earlier rows: count each anchor string once
    for f in sorted((ROOT / "scripts" / "generators").glob("generate_kinetic_core_b*_fit.py")):
        text = f.read_text(encoding="utf-8")
        for m in re.finditer(r'anchor\s*=\s*\(?((?:\s*f?"[^"]*"\s*\+?)+)', text):
            s = " ".join(re.findall(r'"([^"]*)"', m.group(1)))
            if s in seen_anchors:
                continue
            seen_anchors.add(s)
            for a, y in set(pat.findall(s)):
                rows[f"{a} {y}"] += 1
    panel = yaml.safe_load((ROOT / "docs" / "validation" / "directional_claims_panel.yml").read_text(encoding="utf-8"))
    claims: Counter = Counter()
    for c in panel["panel"]:
        ref = (c.get("source") or {}).get("ref") or ""
        m = re.match(r"([a-z\-]+)(\d{4})", ref)
        if m:
            claims[f"{m.group(1).capitalize()} {m.group(2)}"] += 1
    reg = yaml.safe_load((ROOT / "data" / "keys" / "papers.yml").read_text(encoding="utf-8"))["papers"]
    pots: Counter = Counter()
    for p in reg:
        n = sum(len(v) for k, v in p["record_ids"].items() if k.startswith("data/benchmarks"))
        if not n:
            continue
        cit = p.get("citation") or ""
        m = pat.search(cit) or re.match(r"([a-z]+)_?(?:[a-z]+_)*?(\d{4})", p["paper_id"])
        if m:
            pots[f"{m.group(1).capitalize()} {m.group(2)}"] += n
    consts: Counter = Counter()      # rate constants whose registry entry names the paper in its source string
    for f in ("parameters.py", "parameters_furanic.py", "parameters_dicarbonyl.py", "parameters_sulfur.py", "parameters_acrylamide.py",
              "parameters_lipid.py", "trunk_conditions.py", "acrylamide_conditions.py"):
        path = ROOT / "src" / "kinetic_core" / f
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for m in re.finditer(r'source(?:_anchor)?\s*=\s*\(?\s*((?:f?"[^"]*"\s*\+?\s*)+)', text):
            s = " ".join(re.findall(r'"([^"]*)"', m.group(1)))
            for a, y in set(pat.findall(s)):
                consts[f"{a} {y}"] += 1
    total = lambda k: rows[k] + claims[k] + pots[k] + consts[k]
    papers = sorted(set(rows) | set(claims) | set(pots) | set(consts), key=lambda k: -total(k))[:12]
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    y = range(len(papers))
    left = [0] * len(papers)
    for label, counter, colour in (("rate constants taken from the paper", consts, "#8A6BBF"), ("rows the model was tuned on", rows, "#1E2A2C"),
                                   ("benchmark pots it is scored on", pots, "#2B5DA8"), ("statements of direction it is scored on", claims, "#4F80D0")):
        vals = [counter[p] for p in papers]
        ax.barh(list(y), vals, left=left, color=colour, height=0.62, label=label)
        left = [a + b for a, b in zip(left, vals)]
    for i, p in enumerate(papers):
        ax.text(left[i] + 0.6, i, str(left[i]), va="center", fontsize=8.5, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels(papers)
    ax.invert_yaxis()
    ax.set_xlabel("records in the repository that cite the paper")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis="x", color="#E7EBE9", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax.set_title("The papers this model rests on most, counted from its own records", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "17_papers_by_weight.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fig_map()
    fig_repo_flow()
    fig_field_scheme()
    fig_field_coverage()
    fig_how_a_model_works()
    fig_papers_weight()
    fig_funnel()
    fig_scorecard()
    ship = _read(V / "kinetic_core_b16_ship_rule.json")
    scores = _read(V / "core_directional_scores.json")
    fig_schieberle(ship)
    fig_wang(scores)
    fig_liu(ship)
    fig_yiltirak(ship)
    fig_ttca()
    fig_dicarbonyls(scores)
    print(f"wrote 14 figures to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
