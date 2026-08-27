#!/usr/bin/env python3
"""Generates a publication-quality reaction network visualization for Maillard pathways.

Constructs a directed graph of the curated Maillard, Strecker, off-flavor trapping,
and cross-linking pathways using NetworkX and Matplotlib, styling with the LaTeX
science theme if available. Saves the output to results/validation/reaction_network.png.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx  # type: ignore

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.curated_pathways import PATHWAYS  # noqa: E402
from src.plot_style import configure_science_plot_style  # noqa: E402

try:
    configure_science_plot_style()
except Exception as e:
    print(
        f"Warning: Could not configure LaTeX science style: {e}. Falling back to default."
    )
    plt.style.use("default")


# Edge colour per curated reaction family.
#
# 2026-08-27 (audit remediation B): Wave G1 renamed the reaction-family
# vocabulary (src/reaction_templates.py, mirrored by src/curated_pathways.py) and
# this map still carried the pre-G1 keys, so most curated steps fell through to
# the grey `#7f7f7f` default and the legend advertised families the engine never
# emits. The keys below are exactly the families `src.curated_pathways.PATHWAYS`
# emits today (see `curated_reaction_families()` and the coverage check in
# `main()`), so the map cannot silently drift again.
#
# Retired keys and what happened to them:
#   * `Enolisation`      -> split by G1 into `Enolisation_1_2`,
#                           `Enolisation_2_3_Amadori` and `Enolisation_Intermediate`.
#                           The original colour stays with `Enolisation_1_2`; the
#                           other two get neighbouring hues so the lane still
#                           reads as one family group.
#   * `Sugar_Dehydration` -> no longer emitted anywhere in the repo (neither
#                           curated nor engine templates). Its colour is
#                           repointed to `Furanone_Cyclisation`, the curated
#                           cyclisation/dehydration lane that replaced it.
FAMILY_EDGE_COLORS: dict[str, str] = {
    "Schiff_Base_Formation": "#4682B4",
    "Amadori_Rearrangement": "#DAA520",
    "Enolisation_1_2": "#CD5C5C",
    "Enolisation_2_3_Amadori": "#E9967A",
    "Enolisation_Intermediate": "#B03060",
    "Furanone_Cyclisation": "#8FBC8F",
    "Strecker_Degradation": "#BC8F8F",
    "Cysteine_Degradation": "#D2691E",
    "Aminoketone_Condensation": "#6A5ACD",
    "Thiol_Addition": "#48D1CC",
    "Thiol_Addition_Norfuraneol": "#20B2AA",
    "Thiol_Oxidation": "#008B8B",
    "Lipid_Schiff_Base": "#BA55D3",
    "Beta_Elimination": "#CD853F",
    "DHA_Crosslinking": "#B22222",
}

# Colour used when a family has no entry above. Reaching it means the map has
# drifted from the family vocabulary again.
DEFAULT_EDGE_COLOR = "#7f7f7f"


def curated_reaction_families() -> set[str]:
    """Every reaction family `src.curated_pathways.PATHWAYS` actually emits."""
    return {step.reaction_family for steps in PATHWAYS.values() for step in steps}


def uncoloured_reaction_families() -> set[str]:
    """Curated families that would fall through to `DEFAULT_EDGE_COLOR`."""
    return curated_reaction_families() - set(FAMILY_EDGE_COLORS)


def stale_colour_map_keys() -> set[str]:
    """Colour-map keys no curated pathway emits any more."""
    return set(FAMILY_EDGE_COLORS) - curated_reaction_families()


def main() -> int:
    # 1. Create a Directed Graph
    G = nx.DiGraph()

    # We will exclude small co-products/byproducts like water, CO2, H2S, ammonia,
    # and hydrogen to keep the diagram clean and focused on organic backbones.
    EXCLUDED_SPECIES = {"water", "CO2", "H2S", "ammonia", "hydrogen"}

    # Define color scheme (hex colors for premium aesthetics)
    # Precursors: Cool blue/green
    # Intermediates: Soft warm orange
    # Volatiles/Aromas: Vibrant teal
    # Off-flavor traps: Gentle lavender/purple
    # Safety/Damage markers: Bold crimson
    NODE_STYLES = {
        "D-glucose": {
            "color": "#1f77b4",
            "shape": "o",
            "size": 1800,
            "label": "D-Glucose",
        },
        "D-ribose": {
            "color": "#1f77b4",
            "shape": "o",
            "size": 1800,
            "label": "D-Ribose",
        },
        "glycine": {"color": "#2ca02c", "shape": "o", "size": 1400, "label": "Glycine"},
        "L-cysteine": {
            "color": "#2ca02c",
            "shape": "o",
            "size": 1400,
            "label": "L-Cysteine",
        },
        "L-leucine": {
            "color": "#2ca02c",
            "shape": "o",
            "size": 1400,
            "label": "L-Leucine",
        },
        "L-lysine": {
            "color": "#2ca02c",
            "shape": "o",
            "size": 1400,
            "label": "L-Lysine",
        },
        "hexanal": {
            "color": "#e377c2",
            "shape": "o",
            "size": 1400,
            "label": "Hexanal\n(Off-flavor)",
        },
        "glucose-glycine-Schiff-base": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Glucose-Gly\nSchiff Base",
        },
        "ribose-glycine-Schiff-base": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Ribose-Gly\nSchiff Base",
        },
        "glucose-glycine-Amadori": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Glucose-Gly\nAmadori",
        },
        "ribose-glycine-Amadori": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Ribose-Gly\nAmadori",
        },
        "glucose-3-deoxyosone": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Glucose-3-\nDeoxyosone",
        },
        "3-deoxyosone": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "3-Deoxyosone",
        },
        "dehydroalanine": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Dehydroalanine\n(DHA)",
        },
        "pyruvaldehyde": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Pyruvaldehyde",
        },
        "HMF": {
            "color": "#d62728",
            "shape": "D",
            "size": 1800,
            "label": "5-HMF\n(Damage/Desirable)",
        },
        "furfural": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "Furfural\n(Aroma)",
        },
        "2-furfurylthiol": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "2-Furfurylthiol\n(FFT / Meaty)",
        },
        "3-methylbutanal": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "3-Methylbutanal\n(Malty Aroma)",
        },
        "hexanal-glycine-Schiff-base": {
            "color": "#9467bd",
            "shape": "^",
            "size": 1600,
            "label": "Hexanal-Gly\nSchiff (Trapped)",
        },
        "hexanal-lysine-Schiff-base": {
            "color": "#9467bd",
            "shape": "^",
            "size": 1600,
            "label": "Hexanal-Lys\nSchiff (Trapped)",
        },
        "lysinoalanine": {
            "color": "#d62728",
            "shape": "h",
            "size": 1800,
            "label": "Lysinoalanine\n(LAL / Damage)",
        },
        "L-asparagine": {
            "color": "#2ca02c",
            "shape": "o",
            "size": 1400,
            "label": "L-Asparagine",
        },
        "glucose-asparagine-Schiff-base": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Glucose-Asn\nSchiff Base",
        },
        "acrylamide": {
            "color": "#d62728",
            "shape": "D",
            "size": 1800,
            "label": "Acrylamide\n(Safety Hazard)",
        },
        "2-methyl-3-furanthiol": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "2-Methyl-3-\nfuranthiol (MFT)",
        },
        "bis(2-methyl-3-furyl) disulfide": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "MFT Disulfide\n(Meaty Dimer)",
        },
        "2,5-dimethylpyrazine": {
            "color": "#17becf",
            "shape": "D",
            "size": 1800,
            "label": "2,5-Dimethyl-\npyrazine (Aroma)",
        },
        "aminoacetone": {
            "color": "#ff7f0e",
            "shape": "s",
            "size": 1600,
            "label": "Aminoacetone",
        },
    }

    # 2. Add edges from pathways
    family_edge_colors = FAMILY_EDGE_COLORS

    # Coverage check: every family the curated pathways emit must have a colour,
    # otherwise the diagram silently paints real chemistry in the grey default.
    missing = sorted(uncoloured_reaction_families())
    if missing:
        raise RuntimeError(
            "reaction-family colour map is out of date — no colour for: "
            + ", ".join(missing)
            + ". Add them to FAMILY_EDGE_COLORS (do not let them fall through to "
            f"{DEFAULT_EDGE_COLOR})."
        )
    stale = sorted(stale_colour_map_keys())
    if stale:
        print(
            "Warning: FAMILY_EDGE_COLORS carries key(s) no curated pathway emits: "
            + ", ".join(stale)
        )
    print(
        f"Family colour coverage: {len(curated_reaction_families())}/"
        f"{len(curated_reaction_families())} curated families coloured, 0 falling "
        "through to the default."
    )

    for path_name, steps in PATHWAYS.items():
        for step in steps:
            reactants = [
                r.label for r in step.reactants if r.label not in EXCLUDED_SPECIES
            ]
            products = [
                p.label for p in step.products if p.label not in EXCLUDED_SPECIES
            ]
            family = step.reaction_family

            for r in reactants:
                for p in products:
                    if r in NODE_STYLES and p in NODE_STYLES:
                        G.add_edge(
                            r,
                            p,
                            family=family,
                            color=family_edge_colors.get(family, DEFAULT_EDGE_COLOR),
                            source_quality=step.source_quality,
                        )

    # 3. Define Positions (flow from left to right)
    POSITIONS = {
        # Column 0: Sugar Precursors
        "D-glucose": (1.0, 3.5),
        "D-ribose": (1.0, 2.0),
        # Column 1: Amino/Amine Precursors
        "glycine": (2.2, 2.8),
        "L-leucine": (2.2, 0.9),
        "L-cysteine": (2.2, -0.6),
        "L-lysine": (2.2, -2.1),
        "L-asparagine": (2.2, -3.2),
        # Column 2: Off-flavor / Lipid
        "hexanal": (1.0, -3.2),
        # Column 3: Early Intermediates / Schiff Bases
        "glucose-glycine-Schiff-base": (3.6, 3.5),
        "ribose-glycine-Schiff-base": (3.6, 2.0),
        "dehydroalanine": (3.6, -0.6),
        "hexanal-glycine-Schiff-base": (3.6, -2.8),
        "hexanal-lysine-Schiff-base": (3.6, -3.8),
        # Column 4: Amadori / Pyruvaldehyde
        "glucose-glycine-Amadori": (5.2, 3.5),
        "ribose-glycine-Amadori": (5.2, 2.0),
        "pyruvaldehyde": (5.2, 0.5),
        "glucose-asparagine-Schiff-base": (5.2, -4.6),
        # Column 5: Deoxyosones
        "glucose-3-deoxyosone": (6.8, 3.5),
        "3-deoxyosone": (6.8, 2.0),
        "aminoacetone": (6.8, -0.6),
        # Column 6: Endpoints / Aromas & Damage
        "HMF": (8.4, 3.5),
        "furfural": (8.4, 2.0),
        "3-methylbutanal": (8.4, 0.5),
        "lysinoalanine": (8.4, -1.4),
        "2,5-dimethylpyrazine": (8.4, -0.6),
        "acrylamide": (8.4, -4.6),
        # Column 7: Secondary Additions
        "2-furfurylthiol": (10.0, 0.8),
        "2-methyl-3-furanthiol": (10.0, 2.0),
        "bis(2-methyl-3-furyl) disulfide": (11.2, 2.0),
    }

    # Verify all nodes in G have styles and positions
    for node in list(G.nodes):
        if node not in NODE_STYLES or node not in POSITIONS:
            print(f"Removing node {node} because it is not styled or positioned.")
            G.remove_node(node)

    # 4. Set up figure
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.2, 12.0)
    ax.set_ylim(-5.2, 5.2)

    nodes_by_shape: dict[str, list[str]] = {}
    for node in G.nodes:
        shape = str(NODE_STYLES[node]["shape"])
        if shape not in nodes_by_shape:
            nodes_by_shape[shape] = []
        nodes_by_shape[shape].append(node)

    # Draw nodes with decreased alpha for high label contrast
    for shape, nlist in nodes_by_shape.items():
        node_colors = [NODE_STYLES[n]["color"] for n in nlist]
        node_sizes = [NODE_STYLES[n]["size"] for n in nlist]
        nx.draw_networkx_nodes(
            G,
            POSITIONS,
            nodelist=nlist,
            node_color=node_colors,
            node_shape=shape,
            node_size=node_sizes,  # type: ignore
            edgecolors="#4c4c4c",
            linewidths=1.0,
            alpha=0.25,
            ax=ax,
        )

    # Draw node labels
    labels = {node: NODE_STYLES[node]["label"] for node in G.nodes}
    nx.draw_networkx_labels(
        G,
        POSITIONS,
        labels=labels,
        font_size=8,
        font_family="sans-serif",
        font_weight="bold",
        font_color="#1c1c1c",
        ax=ax,
    )

    # Draw edges group by group to indicate intensity / kinetic flux
    for u, v, data in G.edges(data=True):
        fam = data["family"]
        color = data["color"]
        quality = data.get("source_quality", "heuristic")

        # Map family to intensity attributes based on reaction barriers/rates
        if fam == "DHA_Crosslinking":
            width, alpha = 4.0, 0.95
        elif fam in {"Schiff_Base_Formation", "Lipid_Schiff_Base", "Thiol_Addition"}:
            width, alpha = 2.8, 0.85
        elif fam in {
            "Amadori_Rearrangement",
            "Strecker_Degradation",
            "Cysteine_Degradation",
        }:
            width, alpha = 1.8, 0.70
        else:
            # 2026-08-27: comment corrected to the post-Wave-G1 family names.
            # Everything else — the three Enolisation lanes, Beta_Elimination,
            # Furanone_Cyclisation, Aminoketone_Condensation, Thiol_Oxidation and
            # Thiol_Addition_Norfuraneol — carries a high barrier and draws thin.
            # Group membership is unchanged; only the names were stale.
            width, alpha = 1.0, 0.45

        # Style by confidence tier / source quality (Solid = Calibrated, Dashed = Heuristic)
        style = "solid" if quality == "literature" else "dashed"

        nx.draw_networkx_edges(
            G,
            POSITIONS,
            edgelist=[(u, v)],
            edge_color=color,
            width=width,
            alpha=alpha,
            style=style,
            arrowstyle="-|>",
            arrowsize=12,
            connectionstyle="arc3,rad=0.1",
            ax=ax,
        )

    # 5. Build Custom Legend
    node_legend_handles = [
        mpatches.Patch(facecolor="#1f77b4", label="Carbohydrate Precursors (F07/F09)"),
        mpatches.Patch(facecolor="#2ca02c", label="Amino Acid Precursors (F01/F06)"),
        mpatches.Patch(facecolor="#e377c2", label="Lipid Oxidation Precursors (F02)"),
        mpatches.Patch(facecolor="#ff7f0e", label="Reactive Intermediates / Amadori"),
        mpatches.Patch(facecolor="#17becf", label="Desirable Volatiles (Aromas)"),
        mpatches.Patch(facecolor="#9467bd", label="Trapped Off-flavor Compounds"),
        mpatches.Patch(
            facecolor="#d62728", label="Safety / Protein Damage Markers (F12)"
        ),
    ]

    # 2026-08-27: legend lists only the families that actually have an edge in
    # this figure, instead of every key in the colour map — the old legend
    # advertised lanes the drawing does not contain.
    drawn_families = {str(data["family"]) for _, _, data in G.edges(data=True)}
    edge_legend_handles = [
        plt.Line2D([0], [0], color=color, lw=2, label=fam.replace("_", " "))
        for fam, color in family_edge_colors.items()
        if fam in drawn_families
    ]

    # Add a legend entry to define transition intensity
    edge_legend_handles.extend(
        [
            plt.Line2D(
                [0],
                [0],
                color="0.4",
                lw=4.0,
                label="High Intensity (Fast Michael addition)",
            ),
            plt.Line2D(
                [0],
                [0],
                color="0.4",
                lw=2.8,
                label="Medium-High (Fast Schiff/thiol add)",
            ),
            plt.Line2D(
                [0],
                [0],
                color="0.4",
                lw=1.8,
                label="Medium-Low (Moderate rearrangements)",
            ),
            plt.Line2D(
                [0],
                [0],
                color="0.4",
                lw=1.0,
                ls="--",
                label="Low Intensity (Rate-limiting steps)",
            ),
            plt.Line2D(
                [0],
                [0],
                color="0.2",
                lw=2.0,
                ls="-",
                label="Solid: Calibrated (Prior/Literature)",
            ),
            plt.Line2D(
                [0],
                [0],
                color="0.2",
                lw=2.0,
                ls="--",
                label="Dashed: Heuristic (Uncalibrated)",
            ),
        ]
    )

    legend_nodes = ax.legend(
        handles=node_legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.15),
        title=r"\textbf{Species Classification}",
        fontsize=9,
        title_fontsize=10,
        framealpha=0.9,
        edgecolor="0.75",
    )
    ax.add_artist(legend_nodes)

    ax.legend(
        handles=edge_legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.15),
        title=r"\textbf{Reaction Families / Lanes \& Intensities}",
        fontsize=9,
        title_fontsize=10,
        framealpha=0.9,
        edgecolor="0.75",
        ncol=2,
    )

    # (Title removed to prevent overlap with legends and support clean embedding in documents)

    # 6. Save Plot
    output_dir = ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path_png = output_dir / "reaction_network.png"
    output_path_pdf = output_dir / "reaction_network.pdf"

    # Save the file with padding at the top for legends
    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)
    print(
        f"Successfully generated and saved Maillard reaction network plots to {output_path_png} and {output_path_pdf}"
    )

    # Copy to assets dir for docs
    docs_asset_dir = ROOT / "docs" / "assets"
    docs_asset_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(output_path_png, docs_asset_dir / "reaction_network.png")
    shutil.copyfile(output_path_pdf, docs_asset_dir / "reaction_network.pdf")
    print(f"Copied reaction network plots to {docs_asset_dir / 'reaction_network.png'} and {docs_asset_dir / 'reaction_network.pdf'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
