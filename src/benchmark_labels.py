from __future__ import annotations

from typing import Iterable


_BENCHMARK_LABELS: dict[str, dict[str, str]] = {
    "acrylamide_asparagine_glucose_Parker2012": {
        "plain": "Asparagine + glucose, acrylamide benchmark (Parker, 2012)",
        "latex": r"Asparagine + glucose, acrylamide benchmark (Parker, 2012)",
    },
    "acrylamide_spi_extrusion_130C_ACSRef3": {
        "plain": "SPI extrusion, acrylamide benchmark, 130 C (ACS Ref. 3)",
        "latex": r"SPI extrusion, acrylamide benchmark, $130\,^{\circ}$C (ACS Ref. 3)",
    },
    "cml_cel_commercial_pbma_Foods2023": {
        "plain": "Commercial PBMA, CML/CEL benchmark (Foods, 2023)",
        "latex": r"Commercial PBMA, CML/CEL benchmark (Foods, 2023)",
    },
    "cys_glucose_150C_Farmer1999": {
        "plain": "Cysteine + glucose, 150 C (Farmer, 1999)",
        "latex": r"Cysteine + glucose, $150\,^{\circ}$C (Farmer, 1999)",
    },
    "cys_ribose_140C_Hofmann1998": {
        "plain": "Cysteine + ribose, 140 C (Hofmann, 1998)",
        "latex": r"Cysteine + ribose, $140\,^{\circ}$C (Hofmann, 1998)",
    },
    "cys_ribose_150C_Mottram1994": {
        "plain": "Cysteine + ribose, 150 C (Mottram, 1994)",
        "latex": r"Cysteine + ribose, $150\,^{\circ}$C (Mottram, 1994)",
    },
    "furosine_extrusion_crossover_140C_RamirezJimenez2000": {
        "plain": "Furosine crossover, 140 C (Ramirez-Jimenez, 2000)",
        "latex": r"Furosine crossover, $140\,^{\circ}$C (Ram\'irez-Jim\'enez, 2000)",
    },
    "pea_isolate_40C_PratapSingh2021": {
        "plain": "Pea isolate, 40 C (Pratap Singh, 2021)",
        "latex": r"Pea isolate, $40\,^{\circ}$C (Pratap Singh, 2021)",
    },
    "pea_isolate_ribose_cysteine_100C_45min_Internal2026": {
        "plain": "Pea isolate + ribose + cysteine, 100 C (Internal, 2026)",
        "latex": r"Pea isolate + ribose + cysteine, $100\,^{\circ}$C (Internal, 2026)",
    },
    "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026": {
        "plain": "Pea isolate + ribose + cysteine, 100 C (Protocol Pilot, 2026)",
        "latex": r"Pea isolate + ribose + cysteine, $100\,^{\circ}$C (Protocol Pilot, 2026)",
    },
    "pea_isolate_uht_140C_Trikusuma2019": {
        "plain": "Pea isolate UHT, 140 C (Trikusuma, 2019)",
        "latex": r"Pea isolate UHT, $140\,^{\circ}$C (Trikusuma, 2019)",
    },
    "soy_isolate_40C_PratapSingh2021": {
        "plain": "Soy isolate, 40 C (Pratap Singh, 2021)",
        "latex": r"Soy isolate, $40\,^{\circ}$C (Pratap Singh, 2021)",
    },
    "soy_isolate_ribose_cysteine_100C_45min_Internal2026": {
        "plain": "Soy isolate + ribose + cysteine, 100 C (Internal, 2026)",
        "latex": r"Soy isolate + ribose + cysteine, $100\,^{\circ}$C (Internal, 2026)",
    },
    "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026": {
        "plain": "Soy isolate + ribose + cysteine, 100 C (Protocol Pilot, 2026)",
        "latex": r"Soy isolate + ribose + cysteine, $100\,^{\circ}$C (Protocol Pilot, 2026)",
    },
    "spi_hvp_xylose_120C_PMC9905368": {
        "plain": "SPI + HVP + xylose, 120 C (PMC9905368, 2023)",
        "latex": r"SPI + HVP + xylose, $120\,^{\circ}$C (PMC9905368, 2023)",
    },
    "thiamine_cys_ribose_100C_Hofmann1996": {
        "plain": "Thiamine + cysteine + ribose, 100 C (Hofmann, 1996)",
        "latex": r"Thiamine + cysteine + ribose, $100\,^{\circ}$C (Hofmann, 1996)",
    },
    "thiamine_cys_xylose_145C_Cerny2008": {
        "plain": "Thiamine + cysteine + xylose, 145 C (Cerny, 2008) [reference anchor]",
        "latex": r"Thiamine + cysteine + xylose, $145\,^{\circ}$C (Cerny, 2008) [reference anchor]",
    },
    "wheat_gluten_hvp_xylose_120C_PMC9905368": {
        "plain": "Wheat gluten + HVP + xylose, 120 C (PMC9905368, 2023)",
        "latex": r"Wheat gluten + HVP + xylose, $120\,^{\circ}$C (PMC9905368, 2023)",
    },
}


def _fallback_plain(benchmark_id: str) -> str:
    return benchmark_id.replace("_", " ")


def _fallback_latex(benchmark_id: str) -> str:
    return benchmark_id.replace("_", r"\_")


def benchmark_label(benchmark_id: str, *, style: str = "plain") -> str:
    entry = _BENCHMARK_LABELS.get(benchmark_id, {})
    if style == "latex":
        return entry.get("latex", _fallback_latex(benchmark_id))
    return entry.get("plain", _fallback_plain(benchmark_id))


def format_benchmark_list(benchmark_ids: Iterable[str], *, style: str = "plain") -> str:
    labels = benchmark_label_list(benchmark_ids, style=style)
    return ", ".join(labels) if labels else "none"


def benchmark_label_list(benchmark_ids: Iterable[str], *, style: str = "plain") -> list[str]:
    return [benchmark_label(str(benchmark_id), style=style) for benchmark_id in benchmark_ids]
