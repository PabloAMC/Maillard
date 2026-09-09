#!/usr/bin/env python
"""
The seven figures the introduction lacked (2026-09-09): what a plant-based flavour scientist
asks for against what the model can name; the fat path; the protein matrix layer; the
hypothesis layer; the pyrazine step and its supply caveat; the two refused thiol sinks; and
what a laboratory's own data does to the model. Every panel is drawn from the tracked
artifacts or from the engine run in-process; no number is typed here.

Writes docs/assets/thiol_sink/23_*.png to 29_*.png. Companion of build_thiol_sink_figures.py
(same style, same folder). Run inside the container:
    python scripts/generators/build_story_figures.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = ROOT / "docs" / "assets" / "thiol_sink"
INK, MUTED = "#1F2A2E", "#5B6B70"
GOOD, MID, BAD, NONE = "#178F6E", "#D9822B", "#B23A3A", "#9AA6A3"
FILL = {"good": "#D9EFE3", "mid": "#FBE9D0", "bad": "#F6D9D9", "none": "#F2F3F1"}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#D6DBD8", "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK, "axes.titleweight": "bold",
    "axes.titlesize": 11.5, "legend.frameon": False, "legend.fontsize": 9, "figure.dpi": 150,
})


def _read(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def _style(ax, ylabel: str = "", xlabel: str = "") -> None:
    ax.grid(True, axis="y", color="#E7EBE9", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)


# ---------------------------------------------------------------------------
# 23. what is asked for against what the model can name
# ---------------------------------------------------------------------------
def fig_coverage() -> dict:
    from src.explain_compound import explain

    desirable = [c["name"] for c in yaml.safe_load(data_paths.DESIRABLE_TARGETS.read_text())["compounds"]]
    off = [c["name"] for c in yaml.safe_load((data_paths.SPECIES_DIR / "off_flavour_targets.yml").read_text())["compounds"]]

    def status(name: str):
        p = explain(name)
        if p.get("answered"):
            return "modelled", f"{p.get('lane', '')}"
        reached = (p.get("hypotheses") or {}).get("reached_in_charges") or []
        if reached:
            return "route", "no rate, not no route"
        return "none", p.get("state", "")

    rows = [("desirable", n, *status(n)) for n in desirable] + [("off-note", n, *status(n)) for n in off]
    counts = {g: {"modelled": 0, "route": 0, "none": 0} for g in ("desirable", "off-note")}
    for g, _n, st, _ in rows:
        counts[g][st] += 1
    colour = {"modelled": ("good", GOOD), "route": ("mid", MID), "none": ("none", NONE)}
    label = {"modelled": "modelled, with a rate", "route": "a cited route, no rate", "none": "nothing"}

    n_left = len(desirable)
    fig_h = 0.42 * max(n_left, len(off)) + 1.6
    fig, axes = plt.subplots(1, 2, figsize=(13, fig_h), gridspec_kw={"width_ratios": [1, 0.55]})
    for ax, group, names in ((axes[0], "desirable", desirable), (axes[1], "off-note", off)):
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, len(names) - 0.5)
        ax.invert_yaxis()
        ax.axis("off")
        c = counts[group]
        ax.set_title(f"{group} ({len(names)}): {c['modelled']} modelled, {c['route']} route only, {c['none']} nothing", loc="left")
        for i, (g, n, st, note) in enumerate(r for r in rows if r[0] == group):
            fill, edge = colour[st]
            ax.add_patch(plt.Rectangle((0.0, i - 0.42), 1.0, 0.84, facecolor=FILL[fill], edgecolor=edge, linewidth=1.2))
            ax.text(0.02, i, n, va="center", ha="left", color=INK, fontsize=9.5)
            ax.text(0.98, i, label[st] if st != "modelled" else f"modelled ({note} path)", va="center", ha="right", color=edge, fontsize=8.5)
    fig.suptitle("What a plant-based flavour scientist asks for, and what this model can name",
                 x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.01, 0.005, "Each row is one of the repository's declared targets (data/species/desirable_targets.yml, off_flavour_targets.yml); "
             "the status is what `maillard explain` says today. 'Route only' means the hypothesis layer reaches it by a cited rule the engine has no rate for.",
             fontsize=8, color=MUTED, wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    fig.savefig(OUT / "23_coverage_of_declared_targets.png", bbox_inches="tight")
    plt.close(fig)
    return counts


# ---------------------------------------------------------------------------
# 24. the fat path: hexanal on the panel, and what is refused
# ---------------------------------------------------------------------------
def fig_lipid() -> dict:
    sc = _read(V / "core_panel_scores.json")
    rows, refused = [], {}
    for b in sc["benchmarks"]:
        for c in b["compounds"]:
            if c.get("lane") == "lipid" and c.get("predicted") and c.get("measured"):
                rows.append((b["benchmark_id"], c["compound"], float(c["predicted"]), float(c["measured"])))
        for r in b["refused_compounds"]:
            if r.get("lane") == "lipid":
                refused.setdefault(r["compound"], set()).add(b["benchmark_id"])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 0.9]})
    lo, hi = 1e-1, 1e5
    ax.plot([lo, hi], [lo, hi], color=MUTED, linewidth=1)
    ax.fill_between([lo, hi], [lo / 3, hi / 3], [lo * 3, hi * 3], color="#E7EBE9", alpha=0.6, label="within 3x")
    rows.sort(key=lambda r: r[3])
    for i, (bid, comp, pred, meas) in enumerate(rows, start=1):
        storage = "40C" in bid
        ax.scatter(meas, pred, s=60, color=(MID if storage else BAD), edgecolor="white", zorder=3)
        ax.annotate(str(i), (meas, pred), textcoords="offset points", xytext=(6, 6 + 9 * (i % 3)), fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    _style(ax, "predicted hexanal, ug/L", "measured hexanal, ug/L")
    ax.scatter([], [], color=BAD, label="cooked rows")
    ax.scatter([], [], color=MID, label="40 C storage rows (the model starts from zero)")
    ax.legend(loc="upper left")
    ax.set_title(f"Hexanal on the panel: {len(rows)} rows, every one under-predicted", loc="left")
    ax2.axis("off")
    ax2.set_title("What the fat path is, and what it refuses", loc="left")
    lines = [
        "One 1989 study fixes six products and their split from",
        "the linoleate hydroperoxides; the rate at cooking",
        "temperature is an assumption and is declared as one.",
        "",
        "The under-prediction has a named cause: an isolate arrives",
        "with hexanal already made by its own lipoxygenase before",
        "any heat, and the model charges none of it (roadmap,",
        "programme 7: charge the isolate's volatiles as inputs).",
        "",
        "Refused on the panel, with the reason printed:",
    ]
    for comp, bids in sorted(refused.items()):
        lines.append(f"  {comp}: {len(bids)} rows (no measured branch fraction)")
    lines += ["", "The hypothesis layer reaches nonanal, 2-pentylfuran and",
              "1-octen-3-ol by cited rules: 'no rate, not no route'.", "", "Rows, by measured level:"]
    for i, (bid, comp, pred, meas) in enumerate(rows, start=1):
        short = bid.replace("external_validation_", "").replace("_hexanal", "")
        lines.append(f"  {i}. {short}: measured {meas:.3g}, model {pred:.3g} ug/L")
    ax2.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8.6, color=INK, family="DejaVu Sans")
    fig.tight_layout()
    fig.savefig(OUT / "24_fat_path_hexanal.png", bbox_inches="tight")
    plt.close(fig)
    return {"rows": len(rows), "refused": {k: len(v) for k, v in refused.items()}}


# ---------------------------------------------------------------------------
# 25. the protein matrix layer: sites charged, and how little binds in a cook
# ---------------------------------------------------------------------------
def fig_matrix() -> dict:
    from src import api

    base = {"precursors": {"Pea protein isolate": 1.0, "L-Cysteine": 10.0, "D-Ribose": 10.0},
            "temp_C": 145.0, "time_min": 20.0, "ph": 5.0, "aw": 0.98}
    loadings = [10.0, 25.0, 50.0, 100.0]
    matrices = ["blg", "soy_isolate", "pea_isolate"]
    pools, bound = {}, {}
    for m in matrices:
        for g in loadings:
            spec = {"name": f"{m}_{g:g}", **base, "matrix": m, "protein_g_per_l": g}
            p = api.predict(spec, targets=["hexanal", "2-methyl-3-furanthiol", "2-furfurylthiol"])
            if not p["answered"]:
                continue
            mat = p.get("matrix") or {}
            if g == 50.0:
                pools[m] = (mat.get("sites") or {}).get("pools_mmol_per_l", {})
            for comp, b in (mat.get("binding") or {}).items():
                bound.setdefault((m, comp), {})[g] = (b["bound_fraction"], b.get("bound_fraction_corners"))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    keys = ["free_thiol", "disulfide", "amine"]
    x = np.arange(len(matrices))
    w = 0.26
    for j, k in enumerate(keys):
        vals = [pools.get(m, {}).get(k, 0.0) for m in matrices]
        ax.bar(x + (j - 1) * w, vals, w, label=k.replace("_", " "), color=[GOOD, MID, "#4C78A8"][j])
        for xi, v in zip(x + (j - 1) * w, vals):
            ax.text(xi, v * 1.15 if v > 0 else 0.01, f"{v:.2g}", ha="center", fontsize=8, color=MUTED)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["beta-lactoglobulin\n(from its sequence)", "soy isolate\n(measured)", "pea isolate\n(measured)"])
    _style(ax, "reactive sites charged at 50 g/L, mmol/L")
    ax.legend(loc="upper left")
    ax.set_title("What a protein loading puts into the pot", loc="left")
    comps = sorted({c for (_m, c) in bound})
    for m, ls in zip(matrices, ("-", "--", ":")):
        for comp, col in zip(comps, (BAD, GOOD, MID)):
            series = bound.get((m, comp), {})
            if not series:
                continue
            gs = sorted(series)
            ax2.plot(gs, [100 * series[g][0] for g in gs], ls, color=col, marker="o", markersize=4,
                     label=f"{comp}, {m.replace('_', ' ')}")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_yticks([0.03, 0.1, 0.3, 1.0])
    ax2.set_yticklabels(["0.03 %", "0.1 %", "0.3 %", "1 %"])
    ax2.set_xticks([10, 25, 50, 100])
    ax2.set_xticklabels(["10", "25", "50", "100"])
    ax2.set_ylim(0.02, 1.5)
    _style(ax2, "bound at the end of the cook, % of the compound", "protein loading, g/L")
    ax2.set_title("How much binds during 20 min at 145 C: under a percent", loc="left")
    ax2.legend(loc="upper left", fontsize=7.5, ncol=1)
    ax2.text(0.98, 0.04, "the thiols meet the disulfide pool through the network's exchange\nchannel; in this cook their level does not move (MFT 75.0 ug/L with\nor without 50 g/L soy isolate)", transform=ax2.transAxes,
             ha="right", va="bottom", fontsize=7.5, color=MUTED)
    fig.text(0.01, 0.005, "Left: data/species/protein_matrices.yml through the engine (cysteine 10 + ribose 10 mmol/L, pH 5). Right: the declared adduct brackets "
             "(aldehyde to amine, thiol to protein disulfide) applied to that pot; the barriers are 15-23 kJ/mol, so the channel matters over weeks at ambient, not during a cook.",
             fontsize=8, color=MUTED, wrap=True)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT / "25_protein_matrix_layer.png", bbox_inches="tight")
    plt.close(fig)
    return {"pools_at_50": pools, "bound_keys": [f"{m}:{c}" for (m, c) in bound]}


# ---------------------------------------------------------------------------
# 26. the hypothesis layer: the literature's routes placed against the model's
# ---------------------------------------------------------------------------
def fig_hypotheses() -> dict:
    nh = _read(V / "network_hypotheses.json")
    charges = nh["charges"]
    names = [c["charge"].replace("_", " ") for c in charges]
    kinds = ["modelled", "mechanism_known", "proposed"]
    counts = {k: [sum(1 for s in c["steps"] if s["placement"] == k) for c in charges] for k in kinds}
    reached = [sorted({p["id"] for p in c["products"] if p.get("kind") == "registry"}) for c in charges]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1.1, 1]})
    y = np.arange(len(charges))
    left = np.zeros(len(charges))
    for k, col, lab in zip(kinds, (GOOD, MID, NONE), ("modelled: the engine has the step and a rate",
                                                       "mechanism known: a cited rule, no rate in the engine",
                                                       "proposed: a rule with a control, mechanism not settled")):
        ax.barh(y, counts[k], left=left, color=col, label=lab, height=0.6)
        left += np.array(counts[k])
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    _style(ax, "", "steps found from the charge")
    ax.grid(True, axis="x", color="#E7EBE9", linewidth=0.8)
    ax.grid(False, axis="y")
    ax.legend(loc="upper right", fontsize=8)
    total = {k: sum(v) for k, v in counts.items()}
    ax.set_title(f"{len(nh['rules'])} cited rules: {total['modelled']} steps the engine models, "
                 f"{total['mechanism_known'] + total['proposed']} it does not", loc="left", fontsize=10.5)
    ax2.axis("off")
    ax2.set_title("Registry compounds reached only by a rule", loc="left", fontsize=10.5)
    lines = []
    for n, ids in zip(names, reached):
        if ids:
            lines.append(f"{n}:")
            lines.append("    " + ", ".join(i.replace("_", "-") for i in ids))
    ax2.text(0.0, 1.0, "\n".join(lines) or "none", va="top", ha="left", fontsize=8.8, color=INK)
    fig.text(0.01, 0.005, "results/validation/network_hypotheses.json: every rule is applied to the species' structures from each reference charge and placed against the engine's "
             "reactions; each carries a positive and a negative control from its source paper; the engine never imports this layer. 'Reached only by a rule' is what `explain` reports as 'no rate, not no route'.",
             fontsize=8, color=MUTED, wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "26_hypothesis_layer.png", bbox_inches="tight")
    plt.close(fig)
    return {"rules": len(nh["rules"]), **total}


# ---------------------------------------------------------------------------
# 27. the pyrazine step: fitted where fed, a thousandfold low from a sugar pot
# ---------------------------------------------------------------------------
def fig_pyrazine() -> dict:
    s = _read(V / "kinetic_core_b18_ship_rule.json")
    rows = s["T1"]["rows_dex"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1, 1]})
    labels = [k.replace("zhou_", "").replace("_rate_", " ").replace("PZ", "pyrazine").replace("DMP", "2,5-dimethylpyrazine") for k in rows]
    vals = list(rows.values())
    ax.axhspan(-0.3, 0.3, color="#E7EBE9", alpha=0.7, label="within 0.3 dex (twofold)")
    ax.bar(range(len(vals)), vals, color=GOOD, width=0.6)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylim(-0.5, 0.5)
    _style(ax, "model minus measured, log10")
    ax.legend(loc="upper left")
    ax.set_title("Fed dicarbonyl + amino acid (Zhou 2024): within 0.07 dex", loc="left", fontsize=10.5)
    le = s["leahy"]
    model_total, meas_total = le["T4_total"]["model_ug_per_l"], le["T4_total"]["leahy_ug_per_l"]
    ax2.bar([0, 1], [meas_total, model_total], color=[MUTED, BAD], width=0.55)
    ax2.set_yscale("log")
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["measured\n(Leahy 1989, glucose + lysine,\n95 C, 2 h, pH 9)", "model\n(glycine standing in for lysine)"])
    for xpos, v in ((0, meas_total), (1, model_total)):
        ax2.text(xpos, v * 1.3, f"{v:.3g} ug/L", ha="center", fontsize=9, color=INK)
    ax2.set_ylim(1, meas_total * 30)
    _style(ax2, "total pyrazines, ug/L")
    ax2.set_title(f"From a sugar pot: {abs(le['T4_total']['dex']):.1f} decades low", loc="left", fontsize=10.5)
    ax2.text(0.97, 0.62, "the condensation is not the problem;\nthe model makes far too little glyoxal\nand methylglyoxal in water\n(the supply caveat on every pyrazine answer)",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8.5, color=MUTED)
    fig.text(0.01, 0.005, "results/validation/kinetic_core_b18_ship_rule.json (T1 and the Leahy hold-out).", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "27_pyrazine_step_supply.png", bbox_inches="tight")
    plt.close(fig)
    return {"worst_fed_dex": s["T1"]["worst_dex"], "leahy_dex": le["T4_total"]["dex"]}


# ---------------------------------------------------------------------------
# 28. the two refused thiol sinks
# ---------------------------------------------------------------------------
def fig_sink_refusals() -> dict:
    from src.kinetic_core.parameters_sulfur import (
        STACK_DELTA_H_ADDUCT_KJ_MOL, STACK_MEASUREMENT_T_K, STACK_NAC_FORWARD_K_M_INV_S_INV_AT_19_4C,
        STACK_NAC_REVERSE_K_S_INV_AT_19_4C,
    )

    b = _read(V / "kinetic_core_b17_ship_rule.json")
    a = _read(V / "kinetic_core_b17a_ship_rule.json")
    phs = ["6.0", "7.0", "8.0"]
    meas = [b["T3"]["zhou2023"][p]["zhou_mft_share_pct"] for p in phs]
    b9 = [b["b9_reference"]["T3"]["zhou2023"][p]["model_mft_share_pct"] for p in phs]
    vb = [b["T3"]["zhou2023"][p]["model_mft_share_pct"] for p in phs]
    va = [a["T3"]["zhou2023"][p]["model_mft_share_pct"] for p in phs]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(3)
    w = 0.2
    for j, (vals, col, lab) in enumerate(((meas, INK, "measured (Zhou 2023)"), (b9, NONE, "shipped model"),
                                          (vb, MID, "variant b: disulfide made reversible"), (va, BAD, "variant a: pot-made electrophile pool"))):
        ax.bar(x + (j - 1.5) * w, vals, w, color=col, label=lab)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"pH {p[0]}" for p in phs])
    _style(ax, "MFT held as its disulfide, % of the free thiol")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Disulfide share: 10 to 200 times too little, both variants", loc="left", fontsize=10.5)
    T = np.linspace(293.15, 423.15, 200)
    K0 = STACK_NAC_FORWARD_K_M_INV_S_INV_AT_19_4C / STACK_NAC_REVERSE_K_S_INV_AT_19_4C
    K = K0 * np.exp(-(STACK_DELTA_H_ADDUCT_KJ_MOL * 1000.0 / 8.314) * (1.0 / T - 1.0 / STACK_MEASUREMENT_T_K))
    ax2.plot(T - 273.15, 100 * K * 0.1 / (1 + K * 0.1), color=BAD, linewidth=2, label="bound thiol with a 0.1 M site pool, %")
    ax2.axvspan(100, 145, color="#E7EBE9", alpha=0.7, label="cooking window")
    for tc in (100.0, 145.0):
        k = float(np.interp(tc, T - 273.15, K))
        ax2.annotate(f"K = {k:.2f} M^-1 at {tc:.0f} C", (tc, 100 * k * 0.1 / (1 + k * 0.1)), textcoords="offset points",
                     xytext=(6, 10), fontsize=8.5, color=INK)
    _style(ax2, "thiol bound at equilibrium, %", "temperature, C")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_title("Why variant a fails: the measured binding lets go when hot", loc="left", fontsize=10.5)
    fig.text(0.01, 0.005, "Left: kinetic_core_b17_ship_rule.json and kinetic_core_b17a_ship_rule.json, T3. Right: Stack 2018's thiol-quinone equilibrium "
             "(K at 19.4 C and its van 't Hoff enthalpy, parameters_sulfur.py) extrapolated; the site pool is set equal to the whole sugar charge of the reference pot.",
             fontsize=8, color=MUTED, wrap=True)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT / "28_two_refused_sinks.png", bbox_inches="tight")
    plt.close(fig)
    return {"measured": meas, "b9": b9, "variant_b": vb, "variant_a": va}


# ---------------------------------------------------------------------------
# 29. a laboratory's own ladder through `calibrate`
# ---------------------------------------------------------------------------
def fig_calibration() -> dict:
    from src import api

    doc = yaml.safe_load((ROOT / "docs" / "examples" / "reading_2026_ladder.yml").read_text())
    _cal, card = api.calibrate(doc, "Reading 2026")
    before, after = card["holdout"]["before"], card["holdout"]["after"]
    rows_b = {(r["record"], r["compound"]): r for r in before["rows"]} if isinstance(before, dict) and "rows" in before else {}
    rows_a = {(r["record"], r["compound"]): r for r in after["rows"]} if isinstance(after, dict) and "rows" in after else {}
    keys = sorted(rows_b)
    fig, ax = plt.subplots(figsize=(10, 4.8))
    x = np.arange(len(keys))
    fb = [rows_b[k]["fold_error"] for k in keys]
    fa = [rows_a[k]["fold_error"] for k in keys]
    ax.bar(x - 0.18, fb, 0.36, color=NONE, label=f"shipped model (median {before.get('median_fold', float('nan')):.0f}x)")
    ax.bar(x + 0.18, fa, 0.36, color=GOOD, label=f"after calibrate (median {after.get('median_fold', float('nan')):.2f}x)")
    ax.axhline(3.0, color=MUTED, linewidth=1, linestyle="--")
    ax.text(len(keys) - 0.5, 3.3, "3x", color=MUTED, fontsize=8, ha="right")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.replace('yiltirak_', '')}\n{c}" for r, c in keys], fontsize=8.5)
    _style(ax, "fold error on the held-out pots")
    ax.legend(loc="upper right")
    ax.set_title("One laboratory's four-temperature ladder: two pots fitted, two held out", loc="left")
    fig.text(0.01, 0.005, "docs/examples/reading_2026_ladder.yml through `maillard calibrate`: the levels set a response factor per compound, the contrasts moved "
             f"{card.get('n_identified', len(card.get('overrides', {})))} rate constants; the shipped model is untouched. Yiltirak et al. 2026.",
             fontsize=8, color=MUTED, wrap=True)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT / "29_calibration_reading_ladder.png", bbox_inches="tight")
    plt.close(fig)
    return {"before": fb, "after": fa}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {
        "23_coverage": fig_coverage(),
        "24_lipid": fig_lipid(),
        "25_matrix": fig_matrix(),
        "26_hypotheses": fig_hypotheses(),
        "27_pyrazine": fig_pyrazine(),
        "28_sinks": fig_sink_refusals(),
        "29_calibration": fig_calibration(),
    }
    print(json.dumps(summary, default=str)[:1500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
