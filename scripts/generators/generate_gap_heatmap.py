"""S22.3 — render a benchmark × compound gap heatmap.

Reads `results/validation/experiment_value_ranking.json` (S22.1 output) and
draws a PNG heatmap where:

* rows are benchmarks (sorted by max VoI desc),
* columns are compounds (sorted by max VoI desc),
* cell colour = VoI score (yellow→red),
* cell annotation = "*" if measured value sits outside the 90% CI envelope.

The artifact lives at `results/validation/gap_heatmap.png` and is embedded
into the README "Where the next experiments matter most" section.

Run inside the Docker container:

    ./scripts/docker_maillard.sh run "python scripts/generators/generate_gap_heatmap.py"
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "results" / "validation" / "experiment_value_ranking.json"
DEFAULT_OUTPUT = ROOT / "results" / "validation" / "gap_heatmap.png"
DEFAULT_HTML_OUTPUT = ROOT / "results" / "validation" / "experiment_brief_cards.html"

_MEATY_KEYWORDS = ("furanthiol", "furfurylthiol", "methional", "thiazole", "mft", "fft")
_OFFNOTE_KEYWORDS = ("hexanal", "nonanal", "octenal", "pentylfuran", "hexanol")
_SAFETY_KEYWORDS = ("acrylamide", "cml", "cel", "furosine", "hmf")

# Detailed chemical pathway and reaction step mappings for each compound
REACTION_MAPPINGS = {
    "2-methyl-3-furanthiol": {
        "pathway_name": "Sulfur-Maillard (MFT) Pathway",
        "description": "Formed via thiolation of pentose/hexose dehydration intermediates (deoxyosones) with hydrogen sulfide.",
        "steps": [
            "<strong>Sugar Dehydration</strong>: D-Ribose/D-Glucose &rarr; 3-Deoxyosone + H<sub>2</sub>O",
            "<strong>Cysteine Degradation</strong>: L-Cysteine + H<sub>2</sub>O &rarr; H<sub>2</sub>S + Ammonia + Acetaldehyde + CO<sub>2</sub>",
            "<strong>Thiol Addition</strong>: 3-Deoxyosone + H<sub>2</sub>S &rarr; 2-Methyl-3-furanthiol + H<sub>2</sub>O"
        ]
    },
    "2-methyl-3-furanthiol (mft)": {
        "pathway_name": "Sulfur-Maillard (MFT) Pathway",
        "description": "Formed via thiolation of pentose/hexose dehydration intermediates (deoxyosones) with hydrogen sulfide.",
        "steps": [
            "<strong>Sugar Dehydration</strong>: D-Ribose/D-Glucose &rarr; 3-Deoxyosone + H<sub>2</sub>O",
            "<strong>Cysteine Degradation</strong>: L-Cysteine + H<sub>2</sub>O &rarr; H<sub>2</sub>S + Ammonia + Acetaldehyde + CO<sub>2</sub>",
            "<strong>Thiol Addition</strong>: 3-Deoxyosone + H<sub>2</sub>S &rarr; 2-Methyl-3-furanthiol + H<sub>2</sub>O"
        ]
    },
    "2-furfurylthiol": {
        "pathway_name": "Furanthiol (FFT) Pathway",
        "description": "Formed via nucleophilic addition of hydrogen sulfide to furfural under reducing conditions.",
        "steps": [
            "<strong>Sugar Dehydration</strong>: D-Ribose &rarr; Furfural + 3 H<sub>2</sub>O",
            "<strong>Cysteine Degradation</strong>: L-Cysteine + H<sub>2</sub>O &rarr; H<sub>2</sub>S + Ammonia + Acetaldehyde + CO<sub>2</sub>",
            "<strong>Thiol Addition (Reduction)</strong>: Furfural + H<sub>2</sub>S + [H] &rarr; 2-Furfurylthiol + H<sub>2</sub>O"
        ]
    },
    "2-furfurylthiol (fft)": {
        "pathway_name": "Furanthiol (FFT) Pathway",
        "description": "Formed via nucleophilic addition of hydrogen sulfide to furfural under reducing conditions.",
        "steps": [
            "<strong>Sugar Dehydration</strong>: D-Ribose &rarr; Furfural + 3 H<sub>2</sub>O",
            "<strong>Cysteine Degradation</strong>: L-Cysteine + H<sub>2</sub>O &rarr; H<sub>2</sub>S + Ammonia + Acetaldehyde + CO<sub>2</sub>",
            "<strong>Thiol Addition (Reduction)</strong>: Furfural + H<sub>2</sub>S + [H] &rarr; 2-Furfurylthiol + H<sub>2</sub>O"
        ]
    },
    "bis(2-methyl-3-furyl) disulfide": {
        "pathway_name": "MFT Disulfide Dimerization",
        "description": "Oxidative dimerization of 2-methyl-3-furanthiol (MFT) monomer in the presence of trace oxygen/catalysts.",
        "steps": [
            "<strong>MFT Generation</strong>: Pentose + Cysteine &rarr; 2-Methyl-3-furanthiol (MFT)",
            "<strong>Oxidative Dimerization</strong>: 2 MFT + [O] &rarr; Bis(2-methyl-3-furyl) disulfide + H<sub>2</sub>O"
        ]
    },
    "2,5-dimethylpyrazine": {
        "pathway_name": "Alkylpyrazine Condensation",
        "description": "Formed via self-condensation of Strecker-derived aminoacetone intermediates.",
        "steps": [
            "<strong>Sugar Dehydration</strong>: Sugar &rarr; 3-Deoxyosone &rarr; Pyruvaldehyde + Glycolaldehyde",
            "<strong>Strecker Degradation</strong>: Pyruvaldehyde + L-Alanine/Glycine &rarr; Aminoacetone + CO<sub>2</sub>",
            "<strong>Aminoketone Condensation</strong>: 2 Aminoacetone &rarr; 2,5-Dimethylpyrazine + 2 H<sub>2</sub>O + H<sub>2</sub>"
        ]
    },
    "methional": {
        "pathway_name": "Methionine Strecker Degradation",
        "description": "Formed directly via Strecker degradation of L-methionine with dicarbonyl cleavage agents.",
        "steps": [
            "<strong>Dicarbonyl Formation</strong>: Sugar &rarr; 3-Deoxyosone &rarr; Pyruvaldehyde",
            "<strong>Strecker Degradation</strong>: L-Methionine + Pyruvaldehyde &rarr; Methional + Aminoacetone + CO<sub>2</sub>"
        ]
    },
    "acrylamide": {
        "pathway_name": "Asparagine-Maillard Pathway",
        "description": "Formed from the condensation of L-asparagine with reducing sugars followed by decarboxylation.",
        "steps": [
            "<strong>Schiff Base Formation</strong>: L-Asparagine + Glucose &rarr; Schiff Base + H<sub>2</sub>O",
            "<strong>Decarboxylation & Elimination</strong>: Schiff Base &rarr; Acrylamide + CO<sub>2</sub> + other fragments"
        ]
    },
    "nε-(carboxymethyl)lysine (cml)": {
        "pathway_name": "Advanced Glycation End-product (AGE) Pathway",
        "description": "Formed via the oxidative cleavage of Amadori intermediates of lysine.",
        "steps": [
            "<strong>Amadori Formation</strong>: L-Lysine + D-Glucose/D-Ribose &rarr; Amadori Product",
            "<strong>Oxidative Cleavage</strong>: Amadori Product + O<sub>2</sub> &rarr; N&epsilon;-(Carboxymethyl)lysine + Erythronic acid"
        ]
    },
    "nε-(carboxyethyl)lysine (cel)": {
        "pathway_name": "Advanced Glycation End-product (AGE) Pathway",
        "description": "Formed by reaction of lysine with triose phosphates/methylglyoxal.",
        "steps": [
            "<strong>Dicarbonyl Generation</strong>: Sugar &rarr; Methylglyoxal",
            "<strong>Covalent Adduct Formation</strong>: L-Lysine + Methylglyoxal &rarr; N&epsilon;-(Carboxyethyl)lysine"
        ]
    },
    "furosine": {
        "pathway_name": "Furosine Formation (Glycation Index)",
        "description": "Analytical marker formed by acid hydrolysis of lysine-sugar Amadori products.",
        "steps": [
            "<strong>Amadori Formation</strong>: L-Lysine + D-Glucose/D-Galactose &rarr; Amadori Product",
            "<strong>Acid Hydrolysis</strong>: Amadori Product + HCl (Analytical heating) &rarr; Furosine (approx. 33% yield)"
        ]
    },
    "furfural": {
        "pathway_name": "Pentose Dehydration Pathway",
        "description": "Direct acid-catalyzed thermal dehydration of pentose sugars.",
        "steps": [
            "<strong>Pentose Enolization</strong>: D-Ribose/D-Xylose &rarr; 1,2-Enediol",
            "<strong>Triple Dehydration</strong>: 1,2-Enediol &rarr; 3-Deoxyosone &rarr; Furfural + 3 H<sub>2</sub>O"
        ]
    },
    "5-hydroxymethylfurfural (hmf)": {
        "pathway_name": "Hexose Dehydration Pathway",
        "description": "Direct thermal dehydration of hexose sugars (glucose, fructose).",
        "steps": [
            "<strong>Hexose Enolization</strong>: D-Glucose &rarr; 1,2-Enediol",
            "<strong>Triple Dehydration</strong>: 1,2-Enediol &rarr; 3-Deoxyosone &rarr; HMF + 3 H<sub>2</sub>O"
        ]
    },
    "hexanal": {
        "pathway_name": "Lipid Autoxidation Pathway",
        "description": "Formed via homolytic beta-cleavage of linoleic acid 13-hydroperoxide.",
        "steps": [
            "<strong>Lipid Peroxidation</strong>: Linoleic Acid + O<sub>2</sub> &rarr; Linoleic Acid 13-Hydroperoxide",
            "<strong>Homolytic Cleavage</strong>: 13-Hydroperoxide &rarr; Hexanal + 12-Oxo-9-dodecenoic acid"
        ]
    },
    "nonanal": {
        "pathway_name": "Lipid Autoxidation Pathway",
        "description": "Formed via cleavage of oleic acid 9-hydroperoxide.",
        "steps": [
            "<strong>Lipid Peroxidation</strong>: Oleic Acid + O<sub>2</sub> &rarr; Oleic Acid 9-Hydroperoxide",
            "<strong>Homolytic Cleavage</strong>: 9-Hydroperoxide &rarr; Nonanal + other cleavage products"
        ]
    },
    "pyrazine": {
        "pathway_name": "Alkylpyrazine Condensation",
        "description": "Formed via condensation of Strecker-derived aminoketone intermediates.",
        "steps": [
            "<strong>Sugar Fragmentation</strong>: Sugar &rarr; Dicarbonyl (Glyoxal / Methylglyoxal)",
            "<strong>Strecker-like Condensation</strong>: Dicarbonyl + Amino Acid &rarr; Aminoketone + CO<sub>2</sub>",
            "<strong>Dimerization</strong>: 2 Aminoketone &rarr; Pyrazine + 2 H<sub>2</sub>O"
        ]
    }
}

PLANNED_CAMPAIGNS = [
    {
        "id": "pea_iso_ribose_cys_95C_pH5p5_meaty",
        "title": "Pea Protein Isolate (PPI) Meaty-Positive Benchmark",
        "matrix": "pea_iso",
        "conditions": "pH 5.5, 95°C, 0 to 240 minutes",
        "precursors": "D-Ribose (1.0 mM) + L-Cysteine (1.0 mM)",
        "goal": "Quantify meaty volatile generation (MFT, FFT) and lipid oxidation off-flavours (hexanal, 2-pentylfuran) under moderate aqueous heating.",
        "protocol_link": "docs/protocols/pea_matrix_meaty_benchmark.md",
        "protocol_name": "pea_matrix_meaty_benchmark.md",
        "badge_class": "badge-planned"
    },
    {
        "id": "soy_iso_ribose_cys_120C_pH5p8_meaty",
        "title": "Soy Protein Isolate (SPI) High-Severity Meaty Benchmark",
        "matrix": "soy_iso",
        "conditions": "pH 5.8, 120°C, 0 to 240 minutes",
        "precursors": "D-Ribose (1.0 mM) + L-Cysteine (1.0 mM)",
        "goal": "Capture sulfur volatile yields and off-flavour suppression under high-severity conditions, plus safety marker (acrylamide) formation.",
        "protocol_link": "docs/protocols/soy_matrix_meaty_benchmark.md",
        "protocol_name": "soy_matrix_meaty_benchmark.md",
        "badge_class": "badge-planned"
    },
    {
        "id": "matrix_accessibility_assays",
        "title": "Precursor & Denaturation Accessibility Assays",
        "matrix": "pea_iso / soy_iso",
        "conditions": "Ellman (free -SH), OPA (amino exposure), DSC (denaturation)",
        "precursors": "Endogenous protein reactive sites",
        "goal": "Close the gap between physical accessibility kinetics (site exposure) and chemical reaction rates.",
        "protocol_link": "docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md",
        "protocol_name": "PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md",
        "badge_class": "badge-missing"
    }
]


def _short_compound(name: str) -> str:
    name = name.split("(")[0].strip()
    if len(name) > 32:
        name = name[:30] + "…"
    return name


def _short_benchmark(name: str) -> str:
    # Explicit clean formatting map for existing and planned/missing benchmarks
    clean_map = {
        "pea_iso_ribose_cys_95C_pH5p5_meaty": "[MISSING] Pea Isolate Meaty Benchmark (95°C)",
        "soy_iso_ribose_cys_120C_pH5p8_meaty": "[MISSING] Soy Isolate Meaty Benchmark (120°C)",
        "myco_ribose_cys_meaty": "[MISSING] Mycoprotein Meaty Benchmark",
        "wheat_gluten_ribose_cys_meaty": "[MISSING] Wheat Gluten Meaty Benchmark",
        "chickpea_ribose_cys_meaty": "[MISSING] Chickpea Meaty Benchmark",
        "cys_glucose_150C_Farmer1999": "Farmer 1999 (Cys + Glucose, 150°C)",
        "wheat_gluten_hvp_xylose_120C_PMC9905368": "PMC9905368 (Wheat Gluten + HVP + Xylose, 120°C)",
        "spi_hvp_xylose_120C_PMC9905368": "PMC9905368 (SPI + HVP + Xylose, 120°C)",
        "acrylamide_asparagine_glucose_Parker2012": "Parker 2012 (Asparagine + Glucose)",
        "acrylamide_spi_extrusion_130C_ACSRef3": "ACSRef3 (SPI Extrusion, 130°C)",
        "cml_cel_commercial_pbma_Foods2023": "Foods 2023 (Commercial PBMA)",
        "cys_ribose_140C_Hofmann1998": "Hofmann 1998 (Cys + Ribose, 140°C)",
        "cys_ribose_150C_Mottram1994": "Mottram 1994 (Cys + Ribose, 150°C)",
        "furosine_extrusion_crossover_140C_RamirezJimenez2000": "Ramirez-Jimenez 2000 (Extrusion, 140°C)",
        "pea_isolate_40C_PratapSingh2021": "Pratap Singh 2021 (Pea Isolate, 40°C)",
        "pea_isolate_ribose_cysteine_100C_45min_Internal2026": "Internal 2026 (Pea Isolate + Ribose + Cys, 100°C)",
        "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026": "Protocol Pilot 2026 (Pea Isolate + Ribose + Cys, 100°C)",
        "pea_isolate_uht_140C_Trikusuma2019": "Trikusuma 2019 (Pea Isolate UHT, 140°C)",
        "resconi_2023_pbma_beef_identity_benchmark": "Resconi 2023 (PBMA Beef Identity)",
        "soy_isolate_40C_PratapSingh2021": "Pratap Singh 2021 (Soy Isolate, 40°C)",
        "soy_isolate_ribose_cysteine_100C_45min_Internal2026": "Internal 2026 (Soy Isolate + Ribose + Cys, 100°C)",
        "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026": "Protocol Pilot 2026 (Soy Isolate + Ribose + Cys, 100°C)",
        "thiamine_cys_ribose_100C_Hofmann1996": "Hofmann 1996 (Thiamine + Cys + Ribose, 100°C)",
        "thiamine_cys_xylose_145C_Cerny2008": "Cerny 2008 (Thiamine + Cys + Xylose, 145°C)",
    }
    if name in clean_map:
        return clean_map[name]
    
    # Fallback parser for new/unlisted benchmarks
    import re
    # Remove '_benchmark' suffix if present
    clean_name = re.sub(r"_benchmark$", "", name)
    parts = clean_name.split("_")
    
    # Try to find reference/year info
    year_regex = re.compile(r"^(19|20)\d{2}$")
    temp_regex = re.compile(r"^\d+C$")
    
    year_part = ""
    author_part = ""
    temp_part = ""
    rest_parts = []
    
    # Standard mapping for common terms
    term_map = {
        "cys": "Cys",
        "cysteine": "Cys",
        "glucose": "Glucose",
        "ribose": "Ribose",
        "xylose": "Xylose",
        "spi": "SPI",
        "hvp": "HVP",
        "pbma": "PBMA",
        "pea": "Pea",
        "soy": "Soy",
        "isolate": "Isolate",
        "asparagine": "Asparagine",
        "wheat": "Wheat",
        "gluten": "Gluten",
        "thiamine": "Thiamine",
    }
    
    for part in parts:
        if temp_regex.match(part):
            temp_part = part.replace("C", "°C")
        elif year_regex.match(part):
            year_part = part
        else:
            # Check if it has year inside it, like Farmer1999
            m = re.match(r"^([a-zA-Z\-]+)((19|20)\d{2})$", part)
            if m:
                author_part = m.group(1).capitalize()
                year_part = m.group(2)
            else:
                lpart = part.lower()
                if lpart in term_map:
                    rest_parts.append(term_map[lpart])
                elif lpart not in ["min", "45min", "only", "validation"]:
                    rest_parts.append(part.capitalize())
                    
    citation = ""
    if author_part and year_part:
        citation = f"{author_part} {year_part}"
    elif year_part:
        # Look for potential author at the end
        if parts:
            potential_author = parts[-1]
            potential_author = re.sub(r"\d+", "", potential_author)
            if potential_author and potential_author.lower() not in ["benchmark", "validation", "only"]:
                citation = f"{potential_author.capitalize()} {year_part}"
            else:
                citation = year_part
        else:
            citation = year_part
    else:
        citation = name
        
    system = " + ".join(rest_parts)
    if temp_part:
        system = f"{system}, {temp_part}" if system else temp_part
        
    if system and citation != name:
        return f"{citation} ({system})"
    return citation


PLANNED_BENCHMARKS = [
    "pea_iso_ribose_cys_95C_pH5p5_meaty",
    "soy_iso_ribose_cys_120C_pH5p8_meaty",
    "myco_ribose_cys_meaty",
    "wheat_gluten_ribose_cys_meaty",
    "chickpea_ribose_cys_meaty",
]


def build_grid(payload: dict) -> Tuple[List[str], List[str], np.ndarray, np.ndarray, np.ndarray]:
    """Return (benchmark_labels, compound_labels, voi_grid, outside_grid, planned_mask)."""
    candidates = payload.get("candidates", [])
    if not candidates:
        raise ValueError("experiment_value_ranking.json has no candidates")

    bench_max: Dict[str, float] = {}
    comp_max: Dict[str, float] = {}
    for cand in candidates:
        b = str(cand["benchmark_id"])
        c = str(cand["compound"])
        v = float(cand.get("voi_score", 0.0))
        bench_max[b] = max(bench_max.get(b, 0.0), v)
        comp_max[c] = max(comp_max.get(c, 0.0), v)

    bench_order = PLANNED_BENCHMARKS + [b for b, _ in sorted(bench_max.items(), key=lambda kv: -kv[1])]
    comp_order = [c for c, _ in sorted(comp_max.items(), key=lambda kv: -kv[1])]

    n_benchmarks = len(bench_order)
    n_compounds = len(comp_order)

    voi = np.full((n_benchmarks, n_compounds), np.nan)
    outside = np.zeros_like(voi, dtype=bool)
    planned_mask = np.zeros_like(voi, dtype=bool)

    # Mark planned/missing rows
    for idx, b in enumerate(bench_order):
        if b in PLANNED_BENCHMARKS:
            planned_mask[idx, :] = True

    bidx = {b: i for i, b in enumerate(bench_order)}
    cidx = {c: i for i, c in enumerate(comp_order)}

    for cand in candidates:
        b = str(cand["benchmark_id"])
        c = str(cand["compound"])
        if b in bidx and c in cidx:
            i, j = bidx[b], cidx[c]
            voi[i, j] = float(cand.get("voi_score", 0.0))
            if not bool(cand.get("inside_ci", True)):
                outside[i, j] = True
    return bench_order, comp_order, voi, outside, planned_mask


def render(
    bench_labels: List[str],
    comp_labels: List[str],
    voi: np.ndarray,
    outside: np.ndarray,
    planned_mask: np.ndarray,
    output: Path,
) -> Path:
    n_rows, n_cols = voi.shape
    fig_w = max(7.5, 0.45 * n_cols + 4.5)
    fig_h = max(4.5, 0.32 * n_rows + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("YlOrRd")
    cmap.set_bad(color="#eeeeee")
    masked = np.ma.masked_invalid(voi)
    vmax = float(np.nanmax(voi)) if np.isfinite(np.nanmax(voi)) else 1.0
    im = ax.imshow(masked, cmap=cmap, vmin=0.0, vmax=vmax, aspect="auto")

    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))
    ax.set_xticklabels([_short_compound(c) for c in comp_labels], rotation=55, ha="right", fontsize=8)
    ax.set_yticklabels([_short_benchmark(b) for b in bench_labels], fontsize=8)

    # Draw hatch patterns for planned/missing benchmarks
    for i in range(n_rows):
        for j in range(n_cols):
            if planned_mask[i, j]:
                rect = plt.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    fill=False, hatch='//', edgecolor='#bbbbbb', linewidth=0.5
                )
                ax.add_patch(rect)
            elif outside[i, j]:
                ax.text(j, i, "*", ha="center", va="center",
                        color="black", fontsize=9, fontweight="bold")

    ax.set_title(
        "Benchmark × compound experiment-value gaps\n"
        "Cell colour: VoI score (yellow → red). '*' marks outside 90% CI.\n"
        "Hatched cells (//) represent missing high-value physical experiments.",
        fontsize=10,
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("VoI score", fontsize=9)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def _get_log_pcts(p5: float, p50: float, p95: float, measured: float) -> Tuple[float, float, float, float, int, int]:
    epsilon = 1e-6
    vals = [v for v in [p5, p50, p95, measured] if v > 0.0]
    if not vals:
        return 0.0, 0.0, 0.0, 0.0, -4, 4
    
    min_val = min(vals)
    max_val = max(vals)
    
    log_min = math.floor(math.log10(min_val)) - 1
    log_max = math.ceil(math.log10(max_val)) + 1
    
    if log_max - log_min < 4:
        log_max = log_min + 4
        
    log_range = log_max - log_min
    
    def to_pct(val: float) -> float:
        if val <= 0.0:
            return 0.0
        log_v = math.log10(val)
        pct = 100.0 * (log_v - log_min) / log_range
        return min(100.0, max(0.0, pct))
        
    return to_pct(p5), to_pct(p50), to_pct(p95), to_pct(measured), log_min, log_max


def generate_html_briefs(payload: dict, output_html_path: Path) -> Path:
    """Generate a premium glassmorphic HTML dashboard with interactive briefs cards."""
    candidates = payload.get("candidates", []) or []
    
    # Render interactive top calibration cards
    gap_cards_html = []
    for cand in candidates[:15]:  # Display top 15 highest-value gaps
        comp = cand.get("compound", "")
        comp_lower = comp.lower().strip()
        comp_clean = _short_compound(comp)
        bench_id = cand.get("benchmark_id", "")
        bench_clean = _short_benchmark(bench_id)
        
        # Determine compound class for filtering
        comp_class = "kinetics"
        if any(k in comp_lower for k in _MEATY_KEYWORDS):
            comp_class = "meaty"
        elif any(k in comp_lower for k in _OFFNOTE_KEYWORDS):
            comp_class = "off-notes"
        elif any(k in comp_lower for k in _SAFETY_KEYWORDS):
            comp_class = "safety"
            
        # Lookup chemistry pathway information
        reaction_info = REACTION_MAPPINGS.get(comp_lower)
        if not reaction_info:
            # Fallback for unrecognized compounds
            reaction_info = {
                "pathway_name": f"{comp_clean} Synthesis pathway",
                "description": f"Formed via thermal Maillard precursor interaction in the {cand.get('matrix_family', 'unknown')} matrix.",
                "steps": [
                    f"<strong>Precursor Activation</strong>: Matrix components and thermal energy trigger reactive cascades.",
                    f"<strong>Product Condensation</strong>: Intermediates react to yield target {comp_clean}."
                ]
            }
            
        steps_list_html = "".join([f"<li>{step}</li>" for step in reaction_info["steps"]])
        
        # Calculate visual slider percentages
        p5 = float(cand.get("predicted_p5", 0.0))
        p50 = float(cand.get("predicted_p50", 0.0))
        p95 = float(cand.get("predicted_p95", 0.0))
        measured = float(cand.get("measured_ppb", 0.0))
        inside_ci = bool(cand.get("inside_ci", True))
        
        p5_pct, p50_pct, p95_pct, measured_pct, log_min, log_max = _get_log_pcts(p5, p50, p95, measured)
        ci_width_pct = max(1.0, p95_pct - p5_pct)
        
        # Rationale & DoE template info
        rationale = cand.get("rationale", "")
        template = cand.get("suggested_doe_template", "missing_absolute_anchor")
        
        # Standard protocol mapping
        protocol_link = "../../docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md"
        protocol_name = "PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md"
        
        status_class = "inside" if inside_ci else "outside"
        status_badge = '<span class="badge badge-success">Inside 90% CI</span>' if inside_ci else '<span class="badge badge-danger">OUTSIDE 90% CI</span>'
        # Extract and format odour threshold (ODT) if available
        odt_val = cand.get("odour_threshold_ug_per_kg")
        if odt_val is not None:
            try:
                odt_val = float(odt_val)
            except (TypeError, ValueError):
                odt_val = None

        if odt_val is not None and odt_val > 0.0:
            odt_html = f"""
                    <div class="quant-item">
                        <span class="quant-label">Odour Threshold (ODT)</span>
                        <span class="quant-value">{odt_val:.3g} ppb</span>
                    </div>
            """
        else:
            odt_html = ""

        quant_grid_html = f"""
                <!-- Quantitative bounds grid -->
                <div class="quant-grid">
                    <div class="quant-item">
                        <span class="quant-label">Measured (Lit)</span>
                        <span class="quant-value">{measured:.3g} ppb</span>
                    </div>
                    <div class="quant-item">
                        <span class="quant-label">Model Median (P50)</span>
                        <span class="quant-value">{p50:.3g} ppb</span>
                    </div>
                    <div class="quant-item">
                        <span class="quant-label">90% Confidence Interval</span>
                        <span class="quant-value">{p5:.3g} &ndash; {p95:.3g} ppb</span>
                    </div>
                    {odt_html}
                </div>
        """

        card_html = f"""
        <div class="card gap-card" data-matrix="{cand.get('matrix_family')}" data-template="{template}" data-class="{comp_class}">
            <div class="card-header">
                <div class="card-rank">#{cand.get('rank')}</div>
                <div class="card-title-group">
                    <h3>{comp}</h3>
                    <div class="card-subtitle">{bench_clean}</div>
                </div>
                <div class="card-voi" title="Value of Information: how much running this experiment reduces model uncertainty, weighted by sensory relevance (ODT proximity). Higher = higher priority.">{cand.get('voi_score', 0.0):.2f} VoI</div>
            </div>
            
            <div class="card-body">
                <div class="badges-row">
                    <span class="badge badge-matrix">{cand.get('matrix_family', 'unknown')} matrix</span>
                    {status_badge}
                    <span class="badge badge-template" title="Gap type: characterises why this experiment is high-priority">{template.replace('_', ' ')}</span>
                </div>
                
                <p class="rationale-text"><strong>Rationale</strong>: {rationale}</p>
                
                <!-- Chemistry brief section -->
                <div class="chemistry-brief">
                    <h4 class="chem-title">{reaction_info['pathway_name']}</h4>
                    <p class="chem-desc">{reaction_info['description']}</p>
                    <ol class="chem-steps">
                        {steps_list_html}
                    </ol>
                </div>
                
                {quant_grid_html}
                
                <!-- Log scale visualization bar -->
                <div class="scale-section">
                    <div class="scale-header">
                        <span>Model uncertainty vs. observation (log₁₀ ppb) — band width reflects how well-constrained the prediction is</span>
                    </div>
                    <div class="scale-container">
                        <div class="scale-bounds">
                            <span>10<sup>{log_min}</sup></span>
                            <span>10<sup>{log_min + (log_max-log_min)//2}</sup></span>
                            <span>10<sup>{log_max}</sup></span>
                        </div>
                        <div class="scale-bar-track">
                            <div class="scale-bar-ci {status_class}" style="left: {p5_pct}%; width: {ci_width_pct}%;"></div>
                            <div class="scale-p50" style="left: {p50_pct}%;" title="Model Median (P50): {p50:.3g} ppb"></div>
                            <div class="scale-marker" style="left: {measured_pct}%;" title="Measured: {measured:.3g} ppb"></div>
                        </div>
                        <div class="scale-legend">
                            <span class="legend-item"><span class="legend-dot ci-dot {status_class}"></span>90% CI ({p5:.3g} - {p95:.3g} ppb)</span>
                            <span class="legend-item"><span class="legend-dot p50-dot"></span>Median ({p50:.3g} ppb)</span>
                            <span class="legend-item"><span class="legend-dot measured-dot"></span>Measured ({measured:.3g} ppb)</span>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="card-footer">
                <a href="{protocol_link}" class="btn-protocol" target="_blank">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="16" y1="13" x2="8" y2="13"></line><line x1="16" y1="17" x2="8" y2="17"></line><polyline points="10 9 9 9 8 9"></polyline></svg>
                    View {protocol_name}
                </a>
            </div>
        </div>
        """
        gap_cards_html.append(card_html)
        
    gap_cards_merged = "\n".join(gap_cards_html)
    
    # Render Planned Campaign Cards
    planned_cards_html = []
    for camp in PLANNED_CAMPAIGNS:
        steps_html = f"""
        <li><strong>Matrix</strong>: {camp['matrix']}</li>
        <li><strong>Conditions</strong>: {camp['conditions']}</li>
        <li><strong>Precursors</strong>: {camp['precursors']}</li>
        """
        
        card_html = f"""
        <div class="card campaign-card">
            <div class="card-header">
                <div class="card-rank planned-mark">&beta;</div>
                <div class="card-title-group">
                    <h3>{camp['title']}</h3>
                    <div class="card-subtitle">{camp['conditions']}</div>
                </div>
            </div>
            
            <div class="card-body">
                <div class="badges-row">
                    <span class="badge badge-matrix">{camp['matrix']} matrix</span>
                    <span class="badge {camp['badge_class']}">{camp['badge_class'].replace('badge-', '').upper()}</span>
                </div>
                
                <p class="rationale-text"><strong>Campaign Target</strong>: {camp['goal']}</p>
                
                <div class="chemistry-brief">
                    <h4 class="chem-title">Experimental Parameters</h4>
                    <ul class="chem-steps">
                        {steps_html}
                    </ul>
                </div>
            </div>
            
            <div class="card-footer">
                <a href="../../{camp['protocol_link']}" class="btn-protocol" target="_blank">
                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="16" y1="13" x2="8" y2="13"></line><line x1="16" y1="17" x2="8" y2="17"></line><polyline points="10 9 9 9 8 9"></polyline></svg>
                    View {camp['protocol_name']}
                </a>
            </div>
        </div>
        """
        planned_cards_html.append(card_html)
        
    planned_cards_merged = "\n".join(planned_cards_html)
    
    # HTML Layout & CSS styles (highly aesthetic glassmorphic dark theme)
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Maillard Chemistry — Experiment Brief Dashboard</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=Plus+Jakarta+Sans:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    
    <style>
        :root {{
            --bg-color: #fafafa;
            --container-bg: #ffffff;
            --text-primary: #111827;
            --text-secondary: #4b5563;
            --text-muted: #9ca3af;
            --text-accent: #0284c7;
            --accent-success: #059669;
            --accent-danger: #dc2626;
            --accent-warning: #d97706;
            --accent-blue: #2563eb;
            --font-family: 'Plus Jakarta Sans', system-ui, -apple-system, sans-serif;
            --font-header: 'Outfit', sans-serif;
            --card-border: #e5e7eb;
            --card-hover-border: #9ca3af;
            --shadow-subtle: 0 1px 3px rgba(0, 0, 0, 0.05);
            --shadow-hover: 0 4px 12px rgba(0, 0, 0, 0.05);
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            background: var(--bg-color);
            color: var(--text-primary);
            font-family: var(--font-family);
            min-height: 100vh;
            padding: 2.5rem 2rem;
            line-height: 1.5;
        }}

        .container {{
            max-width: 1300px;
            margin: 0 auto;
            background: var(--container-bg);
            padding: 2.5rem;
            border-radius: 12px;
            border: 1px solid var(--card-border);
            box-shadow: var(--shadow-subtle);
        }}

        /* Header Styles */
        header {{
            margin-bottom: 2.5rem;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 2rem;
        }}

        .header-top {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 1rem;
        }}

        h1 {{
            font-family: var(--font-header);
            font-size: 2.2rem;
            font-weight: 700;
            color: var(--text-primary);
            letter-spacing: -0.02em;
        }}

        .subtitle {{
            color: var(--text-secondary);
            font-size: 1rem;
            margin-top: 0.25rem;
        }}

        /* Stats badges in header */
        .header-stats {{
            display: flex;
            gap: 1rem;
        }}

        .stat-badge {{
            background: #f9fafb;
            border: 1px solid var(--card-border);
            border-radius: 8px;
            padding: 0.5rem 1rem;
            text-align: center;
            min-width: 100px;
        }}

        .stat-val {{
            font-family: var(--font-header);
            font-size: 1.2rem;
            font-weight: 700;
            color: var(--text-primary);
        }}

        .stat-lbl {{
            font-size: 0.7rem;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-top: 0.1rem;
        }}

        /* Interactive tab switcher */
        .tabs-row {{
            display: flex;
            gap: 1rem;
            margin-bottom: 2rem;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1rem;
        }}

        .tab-btn {{
            background: none;
            border: none;
            color: var(--text-secondary);
            font-family: var(--font-header);
            font-size: 1.1rem;
            font-weight: 500;
            padding: 0.5rem 1rem;
            cursor: pointer;
            position: relative;
            transition: color 0.2s ease;
        }}

        .tab-btn:hover {{
            color: var(--text-primary);
        }}

        .tab-btn.active {{
            color: var(--text-accent);
        }}

        .tab-btn.active::after {{
            content: '';
            position: absolute;
            bottom: -1rem;
            left: 0;
            width: 100%;
            height: 2px;
            background: var(--text-accent);
        }}

        /* Filters section */
        .filters-row {{
            display: flex;
            gap: 1rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            align-items: center;
        }}

        .filter-group {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
            background: #ffffff;
            border: 1px solid var(--card-border);
            padding: 0.4rem 0.85rem;
            border-radius: 8px;
            box-shadow: var(--shadow-subtle);
        }}

        .filter-group label {{
            font-size: 0.75rem;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.03em;
            font-weight: 600;
        }}

        .filter-group select {{
            background: none;
            border: none;
            color: var(--text-primary);
            font-family: var(--font-family);
            font-size: 0.9rem;
            font-weight: 500;
            cursor: pointer;
            outline: none;
        }}

        .filter-group select option {{
            background: #ffffff;
            color: var(--text-primary);
        }}

        /* Tab Content structures */
        .tab-content {{
            display: none;
        }}

        .tab-content.active {{
            display: block;
        }}

        /* Cards Grid */
        .cards-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(380px, 1fr));
            gap: 1.5rem;
        }}

        /* Card common styles */
        .card {{
            background: #ffffff;
            border: 1px solid var(--card-border);
            border-radius: 12px;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: var(--shadow-subtle);
            transition: transform 0.2s ease, 
                        border-color 0.2s ease, 
                        box-shadow 0.2s ease;
        }}

        .card:hover {{
            transform: translateY(-2px);
            border-color: var(--card-hover-border);
            box-shadow: var(--shadow-hover);
        }}

        .card-header {{
            padding: 1.25rem;
            border-bottom: 1px solid var(--card-border);
            display: flex;
            align-items: flex-start;
            gap: 1rem;
            background: #fafafa;
        }}

        .card-rank {{
            font-family: var(--font-header);
            font-weight: 700;
            font-size: 1.1rem;
            color: var(--text-primary);
            background: #e5e7eb;
            width: 32px;
            height: 32px;
            border-radius: 6px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
        }}

        .planned-mark {{
            color: #4f46e5;
            background: #e0e7ff;
        }}

        .card-title-group {{
            flex-grow: 1;
            min-width: 0;
        }}

        .card-title-group h3 {{
            font-family: var(--font-header);
            font-size: 1.1rem;
            font-weight: 600;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            color: var(--text-primary);
        }}

        .card-subtitle {{
            font-size: 0.8rem;
            color: var(--text-secondary);
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            margin-top: 0.15rem;
        }}

        .card-voi {{
            font-family: var(--font-header);
            font-weight: 600;
            font-size: 0.85rem;
            color: #92400e;
            background: #fef3c7;
            padding: 0.25rem 0.6rem;
            border-radius: 6px;
            border: 1px solid #fde68a;
            flex-shrink: 0;
        }}

        .card-body {{
            padding: 1.25rem;
            flex-grow: 1;
            display: flex;
            flex-direction: column;
            gap: 1rem;
        }}

        /* Badges Row */
        .badges-row {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
        }}

        .badge {{
            font-size: 0.7rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.03em;
            padding: 0.2rem 0.5rem;
            border-radius: 4px;
        }}

        .badge-matrix {{
            background: #f3f4f6;
            color: var(--text-secondary);
            border: 1px solid var(--card-border);
        }}

        .badge-success {{
            background: #ecfdf5;
            color: #047857;
            border: 1px solid #a7f3d0;
        }}

        .badge-danger {{
            background: #fef2f2;
            color: #b91c1c;
            border: 1px solid #fca5a5;
        }}

        .badge-template {{
            background: #f0f9ff;
            color: #0369a1;
            border: 1px solid #bae6fd;
        }}

        .badge-planned {{
            background: #e0e7ff;
            color: #4338ca;
            border: 1px solid #c7d2fe;
        }}

        .badge-missing {{
            background: #fffbeb;
            color: #b45309;
            border: 1px solid #fde68a;
        }}

        .rationale-text {{
            font-size: 0.85rem;
            color: var(--text-secondary);
        }}

        /* Chemistry Section inside Cards */
        .chemistry-brief {{
            background: #fafafa;
            border: 1px solid var(--card-border);
            border-radius: 8px;
            padding: 0.85rem;
        }}

        .chem-title {{
            font-size: 0.85rem;
            font-weight: 600;
            color: var(--text-accent);
            margin-bottom: 0.25rem;
        }}

        .chem-desc {{
            font-size: 0.78rem;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }}

        .chem-steps {{
            font-size: 0.78rem;
            color: var(--text-primary);
            padding-left: 1.1rem;
        }}

        .chem-steps li {{
            margin-bottom: 0.35rem;
        }}

        .chem-steps li:last-child {{
            margin-bottom: 0;
        }}

        /* Quantitative bounds grid */
        .quant-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 0.75rem;
            background: #fafafa;
            border: 1px solid var(--card-border);
            border-radius: 8px;
            padding: 0.85rem;
            font-size: 0.8rem;
        }}

        .quant-item {{
            display: flex;
            flex-direction: column;
            gap: 0.15rem;
        }}

        .quant-label {{
            font-size: 0.7rem;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.03em;
            font-weight: 600;
        }}

        .quant-value {{
            font-weight: 700;
            color: var(--text-primary);
            font-family: var(--font-header);
        }}

        /* Log scale predicted vs. measured bar */
        .scale-section {{
            margin-top: auto;
            border-top: 1px solid var(--card-border);
            padding-top: 1rem;
        }}

        .scale-header {{
            font-size: 0.75rem;
            font-weight: 600;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.03em;
            margin-bottom: 0.5rem;
        }}

        .scale-container {{
            background: #fafafa;
            border: 1px solid var(--card-border);
            border-radius: 6px;
            padding: 0.75rem;
        }}

        .scale-bounds {{
            display: flex;
            justify-content: space-between;
            font-size: 0.7rem;
            color: var(--text-muted);
            margin-bottom: 0.25rem;
        }}

        .scale-bar-track {{
            height: 8px;
            background: #e5e7eb;
            border-radius: 4px;
            position: relative;
            margin-bottom: 0.5rem;
        }}

        .scale-bar-ci {{
            height: 100%;
            border-radius: 4px;
            position: absolute;
        }}

        .scale-bar-ci.inside {{
            background: var(--accent-success);
            opacity: 0.3;
        }}

        .scale-bar-ci.outside {{
            background: var(--accent-danger);
            opacity: 0.3;
        }}

        .scale-p50 {{
            position: absolute;
            width: 3px;
            height: 14px;
            background: var(--text-primary);
            top: -3px;
            border-radius: 1px;
        }}

        .scale-marker {{
            position: absolute;
            width: 10px;
            height: 10px;
            background: var(--accent-blue);
            border: 2px solid #ffffff;
            border-radius: 50%;
            top: -1px;
            transform: translateX(-5px);
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.2);
            z-index: 2;
        }}

        .scale-legend {{
            display: flex;
            justify-content: space-between;
            font-size: 0.65rem;
            color: var(--text-secondary);
            flex-wrap: wrap;
            gap: 0.5rem;
        }}

        .legend-item {{
            display: flex;
            align-items: center;
            gap: 0.25rem;
        }}

        .legend-dot {{
            width: 6px;
            height: 6px;
            border-radius: 50%;
            display: inline-block;
        }}

        .ci-dot.inside {{ background: var(--accent-success); }}
        .ci-dot.outside {{ background: var(--accent-danger); }}
        .p50-dot {{ background: var(--text-primary); width: 2px; height: 8px; border-radius: 0; }}
        .measured-dot {{ background: var(--accent-blue); border: 1px solid #ffffff; }}

        .card-footer {{
            padding: 1rem 1.25rem;
            background: #fafafa;
            border-top: 1px solid var(--card-border);
            display: flex;
            justify-content: flex-end;
        }}

        .btn-protocol {{
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            color: var(--text-accent);
            text-decoration: none;
            font-size: 0.8rem;
            font-weight: 600;
            transition: color 0.2s ease;
        }}

        .btn-protocol:hover {{
            color: #0369a1;
        }}

        .btn-protocol svg {{
            flex-shrink: 0;
        }}

        /* Explainer legend row */
        .explainer-row {{
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            margin-top: 1.5rem;
            padding-top: 1.25rem;
            border-top: 1px solid var(--card-border);
        }}

        .explainer-item {{
            flex: 1 1 280px;
            display: flex;
            gap: 0.6rem;
            align-items: flex-start;
            background: #fafafa;
            border: 1px solid var(--card-border);
            border-radius: 8px;
            padding: 0.75rem 1rem;
            font-size: 0.8rem;
            color: var(--text-secondary);
            line-height: 1.5;
        }}

        .explainer-icon {{
            font-size: 1.1rem;
            flex-shrink: 0;
            margin-top: 0.05rem;
        }}

        .explainer-item strong {{
            color: var(--text-primary);
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="header-top">
                <div>
                    <h1>Maillard Chemistry</h1>
                    <div class="subtitle">High-Value Experiment Briefs — ranked by Value of Information (VoI)</div>
                </div>
                <div class="header-stats">
                    <div class="stat-badge">
                        <div class="stat-val">{len(candidates)}</div>
                        <div class="stat-lbl">Calibration Gaps Ranked</div>
                    </div>
                    <div class="stat-badge">
                        <div class="stat-val">{sum(1 for c in candidates if not c.get('inside_ci'))}</div>
                        <div class="stat-lbl">Mismatches (Outside 90% CI)</div>
                    </div>
                </div>
            </div>
            <div class="explainer-row">
                <div class="explainer-item">
                    <span class="explainer-icon">📊</span>
                    <span><strong>VoI (Value of Information)</strong> — how much running this experiment would reduce prediction uncertainty, weighted by the compound&#8217;s sensory relevance (proximity to odour detection threshold). Higher = more worth doing first.</span>
                </div>
                <div class="explainer-item">
                    <span class="explainer-icon">⚠️</span>
                    <span><strong>Benchmark gap</strong> — each card represents an existing literature benchmark study where the model has high uncertainty or mismatches the published measurement. The brief outlines the proposed action (using the suggested DoE template) to re-run or refine the chemistry to resolve this calibration gap.</span>
                </div>
                <div class="explainer-item">
                    <span class="explainer-icon">📈</span>
                    <span><strong>Uncertainty bar</strong> — shows the model&#8217;s predicted 90% CI (coloured band) and median (white line) vs. the known literature measurement (blue dot) on a log₁₀ ppb scale. Green = observation falls within CI; red = mismatch.</span>
                </div>
            </div>
        </header>



        <!-- Calibration Gaps Tab -->
        <div id="tab-gaps" class="tab-content active">
            <div class="filters-row">
                <div class="filter-group">
                    <label for="matrix-filter">Matrix:</label>
                    <select id="matrix-filter" onchange="filterCards()">
                        <option value="all">All Matrices</option>
                        <option value="free">Free Precursor</option>
                        <option value="pea_iso">Pea Isolate</option>
                        <option value="soy_iso">Soy Isolate</option>
                        <option value="wheat_gluten">Wheat Gluten</option>
                    </select>
                </div>
                <div class="filter-group">
                    <label for="template-filter">Gap Type:</label>
                    <select id="template-filter" onchange="filterCards()">
                        <option value="all">All Gap Types</option>
                        <option value="blocking_benchmark_gap">Wide CI — benchmark exists, model poorly constrained</option>
                        <option value="missing_absolute_anchor">Missing absolute concentration anchor</option>
                        <option value="missing_kinetic_dataset">Missing kinetic time-course data</option>
                        <option value="missing_positive_flavor_anchor">Missing positive flavour anchor</option>
                    </select>
                </div>
            </div>

            <div class="cards-grid">
                {gap_cards_merged}
            </div>
        </div>


    </div>

    <script>

        function filterCards() {{
            const matrix = document.getElementById('matrix-filter').value;
            const template = document.getElementById('template-filter').value;
            
            document.querySelectorAll('.gap-card').forEach(card => {{
                const cardMatrix = card.dataset.matrix;
                const cardTemplate = card.dataset.template;
                
                const matrixMatch = (matrix === 'all' || cardMatrix === matrix);
                const templateMatch = (template === 'all' || cardTemplate === template);
                
                if (matrixMatch && templateMatch) {{
                    card.style.display = 'flex';
                }} else {{
                    card.style.display = 'none';
                }}
            }});
        }}
    </script>
</body>
</html>
"""
    output_html_path.parent.mkdir(parents=True, exist_ok=True)
    output_html_path.write_text(html_content, encoding="utf-8")
    return output_html_path


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--output-html", type=Path, default=DEFAULT_HTML_OUTPUT)
    args = parser.parse_args(argv)

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    bench_labels, comp_labels, voi, outside, planned_mask = build_grid(payload)
    output = render(bench_labels, comp_labels, voi, outside, planned_mask, args.output)
    output_html = generate_html_briefs(payload, args.output_html)
    
    finite_voi = voi[np.isfinite(voi)]
    print(
        f"wrote {output.relative_to(ROOT)} | "
        f"{len(bench_labels)} benchmarks × {len(comp_labels)} compounds | "
        f"max VoI={float(np.nanmax(finite_voi)):.2f} | "
        f"outside-CI cells={int(outside.sum())}"
    )
    print(f"wrote {output_html.relative_to(ROOT)} (Interactive Experiment Briefs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
