#!/usr/bin/env python
"""
Figures for docs/guides/WHERE_THE_THIOLS_GO.md (2026-09-07).

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


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ship = _read(V / "kinetic_core_b16_ship_rule.json")
    scores = _read(V / "core_directional_scores.json")
    fig_schieberle(ship)
    fig_wang(scores)
    fig_liu(ship)
    fig_yiltirak(ship)
    fig_ttca()
    fig_dicarbonyls(scores)
    print(f"wrote 6 figures to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
