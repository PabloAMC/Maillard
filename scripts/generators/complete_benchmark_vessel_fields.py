#!/usr/bin/env python3
"""
scripts/generators/complete_benchmark_vessel_fields.py

PROGRAMME STEP R1 (2026-09-06) -- THE VESSEL-FIELD COMPLETION OF EVERY PANEL BUNDLE.

WHY THIS EXISTS
===============
`tasks/data_restructure_plan.md`, "Reaction-modelling programme", finding: the largest
between-laboratory term in the sulfur lane's residuals is the oxygen the closed vessel held
above the sample, which no bundle recorded. This script writes exactly one new key per bundle,
``conditions.vessel``, from each bundle's SOURCE PAPER (PDF on disk) or from a verbatim methods
quotation the bundle already carries, with a per-field provenance class and note
(FIT_HOLDOUT_DECLARATION.md Amendment 19). It edits nothing else; the hash test in
tests/unit/test_kinetic_core_b2_3.py proves it.

PROVENANCE CLASSES (the buffer completion's, plus one)
=====================================================
  primary_source_pdf           read from a PDF in data/articles/ (named in the note)
  repo_verbatim_methods_quote  read from a verbatim quotation the bundle already carries;
                               source NOT ON DISK; second-hand, and the note says so
  unknown                      neither; recorded as unstated, never guessed
  not_applicable               no cook (unheated ingredient, commercial product, synthetic)

Idempotent: re-running rewrites the same block. ``--check`` reports drift and exits 1.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

# ---------------------------------------------------------------------------
# Shared texts
# ---------------------------------------------------------------------------

_HOFMANN_PDF = "data/articles/hofmann1998.pdf"
_HOFMANN_NOTE_100 = (
    "Hofmann & Schieberle 1998, JAFC 46:235, read from " + _HOFMANN_PDF + " (p. 236, EXPERIMENTAL "
    "PROCEDURES, 'Model Experiments'), verbatim: 'the reactants (amounts detailed in the tables) were "
    "dissolved in phosphate buffer (0.5 mol/L) at pH 3.0, 5.0, or 7.0 and thermally treated for 20 min at "
    "145 degC in a laboratory autoclave (200 mL; Type II; Roth, Germany)'. Table 1 footnote, verbatim: "
    "'A solution (100 mL) of cysteine (3.3 mmol) and the corresponding carbohydrate (10.0 mmol) was "
    "reacted'. Headspace = 200 - 100 = 100 mL. ATMOSPHERE: no inert-gas step appears anywhere in the "
    "procedures, so the headspace is taken as laboratory air at closure (air_by_default). Water source "
    "for the buffer is not stated. Stirring not stated."
)
_HOFMANN_NOTE_50 = (
    "Hofmann & Schieberle 1998, JAFC 46:235, read from " + _HOFMANN_PDF + " (p. 236, EXPERIMENTAL "
    "PROCEDURES, 'Model Experiments'): 'thermally treated for 20 min at 145 degC in a laboratory "
    "autoclave (200 mL; Type II; Roth, Germany)'. This bundle's own quoted table footnote gives the "
    "fed-intermediate system volume verbatim as 'An aqueous solution (50 mL)' of 1 mmol of each "
    "precursor. Headspace = 200 - 50 = 150 mL. ATMOSPHERE: no inert-gas step appears anywhere in the "
    "procedures, so the headspace is taken as laboratory air at closure (air_by_default). Water source "
    "not stated. Stirring not stated."
)
_YILTIRAK_NOTE = (
    "Yiltirak et al. 2026, Food Res. Int. 231:118600, read from data/articles/Yiltirak2026.pdf on "
    "2026-09-06. Sec. 2.4, verbatim: 'model systems (3 mL) ... were transferred to 20-mL Duran tubes with "
    "screw-caps lined with PTFE gaskets (Duran Wheaton Kimble, Wertheim, Germany) ... in a heating block "
    "equipped with a magnetic stirrer ... PTFE-covered magnetic stir bars (12 mm x 4.5 mm) were used ... "
    "After heating, the tubes were cooled immediately to room temperature in an ice bath ... Subsequently, "
    "they were flushed with Pureshield Argon (B.O.C., UK) to prevent further oxidation'. The argon comes "
    "AFTER the cook: the 17 mL headspace was laboratory air during heating (air_by_default). Sec. 2.1: "
    "the buffer was made 'in tap water'. Dossier: data/lit/extraction_dossiers/yiltirak2026_extraction.md."
)
_BOLTON_NOTE = (
    "Bolton, Reineccius & Liardon 1994 (ACS Symp. Ser. 543, ch. 22), read from data/articles/bolton1993.pdf "
    "on 2026-09-04; this bundle's content_verification.quoted_conditions, verbatim: 'Samples were prepared "
    "in triplicate using 10 g of the model system and 23.3 g of distilled water ... thermally processed at "
    "120 C for one hour in 125 mL glass vials sealed with Teflon-lined septa'. Fill taken as 33.3 mL at "
    "1 g/mL (a 30 %-solids slurry is denser, ~1.1 g/mL, so the true headspace is ~92-95 mL; the 91.7 mL "
    "used here errs by < 4 %). No inert-gas step: air_by_default. Water: distilled. Stirring not stated."
)
_SCHIBILSKY_NOTE = (
    "Schibilsky 2019 (dissertation), SOURCE NOT ON DISK; second-hand, from this bundle's own "
    "content_verification.quoted_method, verbatim (German): 'Die Thermolysen erfolgten in abgeschmolzenen "
    "10 mL Glasspiessampullen (Neolab) oder in 20 mL Headspace-Vials mit Boerdelkappe (Amchro) ... Das "
    "Probenvolumen in den Glasspiessampullen betrug 2 mL und in den Headspace-Vials 6 mL.' The thesis "
    "does not say WHICH vessel the glucose/alanine model (5.1.6.2) used, so fill and vessel are left null "
    "rather than chosen: the two candidates give 8 or 14 mL of air (0.07 or 0.12 mmol O2). Sealed under "
    "air in either case (air_by_default); 'bidestilliertem Wasser' = distilled."
)
_CHANG_NOTE = (
    "Chang et al. 2021, SOURCE NOT ON DISK; second-hand, from this bundle's own "
    "content_verification.quoted_method, verbatim: 'The pH value of each solution was first adjusted to 5.8 "
    "... and topped up with distilled water to 100 mL. Each solution was heated at 180 degC for 10, 20, and 30 "
    "min, followed by cooling in tap water.' The quotation names the solution volume (100 mL) but neither "
    "the vessel nor whether it was closed, so vessel volume is null and the atmosphere is unstated. Water: "
    "distilled."
)
_LIN_NOTE = (
    "Lin et al. 2022, Polymers 14:1565, SOURCE NOT ON DISK; second-hand, from this bundle's own "
    "content_verification.quoted_method, verbatim: 'Each solution was topped up with distilled water to 100 "
    "mL. Thereafter, each solution was placed in a dry-bath incubator (DB200-2 ...) maintained at 180 degC "
    "for 10, 20, and 30 min'. A dry-bath incubator holds tubes; the quotation names neither the tube nor a "
    "closure, so vessel volume is null and the atmosphere is unstated. Water: distilled."
)
_YE_NOTE = (
    "Ye et al. 2024, SOURCE NOT ON DISK; second-hand, from this bundle's own content_verification."
    "quoted_method, verbatim: 'an equimolar solution of glucose and asparagine was accurately prepared in "
    "phosphate buffer (0.1 M, pH 6.86), and 4 mL of the solution was transferred to a 25 mL thick-walled "
    "pressurized glass tube in an oil bath ... at 180 degC and kept for 30 min.' Headspace 21 mL. No "
    "inert-gas step in the quotation: air_by_default. Water source for the buffer not stated."
)
_STEINHAGEN_NOTE = (
    # WAVE B35 (2026-09-11): this note used to open "SOURCE NOT ON DISK". It is on disk, and was when
    # the note was written: the paper is Leitzen et al. 2021, data/articles/Leitzen2021.pdf since
    # 2026-09-07, dossier leitzen2021_extraction.md. B34 corrected this bundle's citation and its
    # BUFFER note and missed this one; the audit that followed it found the miss. The false claim is
    # kept below, labelled, because it is the audit record.
    "WAVE B35 (2026-09-11): THE SOURCE IS ON DISK -- Leitzen et al. 2021 (Pharmaceuticals 14:1121), "
    "data/articles/Leitzen2021.pdf, dossier leitzen2021_extraction.md. ||| PRIOR NOTE, SUPERSEDED AND "
    "RETAINED AS THE AUDIT RECORD: "
    "Steinhagen et al. 2021, SOURCE NOT ON DISK; second-hand, from this bundle's own content_verification."
    "quoted_method: 'Glucose solutions (10%, w/v) were autoclaved at 111 degC, 116 degC, and 121 degC for "
    "different durations.' The quotation gives neither the vessel, the fill nor the closure; autoclaving "
    "implies a vented or loosely closed container but that is an inference, so everything is unstated."
)
_ACS_NOTE = (
    "Ma et al. 2024, Int. J. Mol. Sci. 25:8668, twin-screw extrusion of SPI + corn starch at 30 % moisture "
    "(die-zone 130 C; the bundle's 25 s is the residence time). NOT ON DISK as a PDF; the repo carries its "
    "Europe PMC full text for the directional panel (second-hand). An extruder has no fixed headspace: "
    "oxygen access is process-specific (continuous_process). Water source not applicable to a melt."
)
_FOODS_NOTE = (
    "Fu et al. 2023 (Foods, PMC10217484): a survey of COMMERCIAL plant-based meat analogues; the bundle's "
    "conditions block is a proxy operating point, not a cook the paper ran. There is no vessel to record "
    "(not_applicable). Verified from the PMC full text on 2026-09-04; see content_verification."
)
_PRATAP_NOTE = (
    "Pratap-Singh et al. 2021, Molecules 26:4104: an UNHEATED protein powder; the bundle's 40 C / 10 min "
    "is the HS-SPME incubation (1 g powder in 7 mL water), not a thermal process, as the 2026-09-04 "
    "calibration diagnosis records. No cook, so no vessel (not_applicable). SOURCE NOT ON DISK."
)
_TRIKUSUMA_NOTE = (
    "Trikusuma, Paravisini & Peterson 2020, Food Chem. 312:126082, read from data/articles/trikusuma2020.pdf "
    "on 2026-09-04; this bundle's content_verification.quoted_conditions, verbatim: 'aseptic UHT processing "
    "... processed at a final temperature and pressure of 140 C and 0.55 MPa for 6 s'. A UHT line is a "
    "continuous closed process with no fixed headspace (continuous_process); dissolved oxygen in the "
    "beverage is not reported. Water source for the 3 % w/w isolate beverage not stated."
)
_RESCONI_NOTE = (
    "Hernandez, Woerner, Brooks & Legako 2023, Molecules 28:3151: a COMMERCIAL plant-based meat analogue "
    "cooked as a reference product; the bundle's 150 C / 60 min is a proxy for an unspecified commercial "
    "process (process_metadata.extrusion_history = commercial_pbma_unknown). No vessel to record "
    "(not_applicable). SOURCE NOT ON DISK."
)
_BI_RAW_NOTE = (
    "Bi et al. 2020, JAFC 68:2718, read from data/articles/bi2020.pdf on 2026-09-04: this bundle is RAW "
    "pea flour, never heated (the 40 C / 10 min block is the HS-SPME incubation, p. 2719). No cook, so no "
    "vessel (not_applicable). The roasted-pea twin bundle carries the oven."
)
_BI_ROAST_NOTE = (
    "Bi et al. 2020, JAFC 68:2718, read from data/articles/bi2020.pdf on 2026-09-04; this bundle's "
    "content_verification.quoted_conditions, verbatim: 'The roasted peas (500 g) were prepared using an "
    "Isotemp forced-air oven at 160 C for 30 min.' An open forced-air oven: oxygen is not limited by a "
    "headspace (open). Water source not applicable to dry roasting."
)
_LI_NOTE = (
    "Li et al. 2026, Foods 15(5):912, SOURCE NOT ON DISK; second-hand, from the PMC full text the repo read "
    "on 2026-09-03 (content_verification): a high-moisture extrusion melt (57 wt % moisture, 160 C). An "
    "extruder has no fixed headspace; oxygen access is process-specific (continuous_process). Water source "
    "not stated."
)
_LIU_NOTE = (
    "Liu, Cadwallader & Drake 2023, Food Chem. 406:134998: commercial pea protein rehydrated to 10 % solids "
    "in deionized water and NEVER HEATED; the bundle's 40 C / 10 min is the headspace equilibration, as its "
    "buffer block already records. No cook, so no vessel (not_applicable); the PDF is on disk "
    "(data/articles/liu2023.pdf, read 2026-09-04)."
)
_SYNTHETIC_NOTE = (
    "Synthetic snapshot (_Internal2026): legacy-model output, not a measurement, kept off the scored panel "
    "by src/kinetic_core/panel.is_scored_bundle. There is no experiment and no vessel (not_applicable)."
)


def _block(fill, vessel, closure, atmosphere, water, stirred, pclass, note) -> Dict[str, Any]:
    return {
        "fill_mL": fill, "vessel_mL": vessel, "closure": closure, "atmosphere": atmosphere,
        "water_source": water, "stirred": stirred, "provenance_class": pclass, "provenance_note": note,
    }


_AUTOCLAVE = "laboratory autoclave (200 mL; Type II; Roth, Germany), sealed"
HOF100 = _block(100.0, 200.0, _AUTOCLAVE, "air_by_default", "unstated", None, "primary_source_pdf", _HOFMANN_NOTE_100)
HOF50 = _block(50.0, 200.0, _AUTOCLAVE, "air_by_default", "unstated", None, "primary_source_pdf", _HOFMANN_NOTE_50)
YIL = _block(3.0, 20.0, "20 mL Duran tube, screw cap with PTFE gasket, magnetically stirred", "air_by_default", "tap", True, "primary_source_pdf", _YILTIRAK_NOTE)
NA = lambda note, pclass="not_applicable": _block(None, None, "no cook", "not_applicable", "not_applicable", None, pclass, note)  # noqa: E731

BLOCKS: Dict[str, Dict[str, Any]] = {
    # --- trust loop (data/benchmarks/*.json) ---
    "acrylamide_spi_extrusion_130C_ACSRef3": _block(None, None, "twin-screw extruder", "continuous_process", "not_applicable", None, "repo_verbatim_methods_quote", _ACS_NOTE),
    "cml_cel_commercial_pbma_Foods2023": NA(_FOODS_NOTE),
    "hofmann1998_c2c3_recombination_145C_20min_pH3": HOF50,
    "hofmann1998_c2c3_recombination_145C_20min_pH5": HOF50,
    "hofmann1998_c2c3_recombination_145C_20min_pH7": HOF50,
    "hofmann1998_fructose_cysteine_145C_20min_pH5": HOF100,
    "hofmann1998_furan2aldehyde_h2s_145C_20min_pH5": HOF50,
    "hofmann1998_glucose_cysteine_145C_20min_pH5": HOF100,
    "hofmann1998_norfuraneol_cysteine_145C_20min_pH5": HOF50,
    "hofmann1998_norfuraneol_h2s_145C_20min_pH5": HOF50,
    "hofmann1998_ribose_cysteine_145C_20min_pH5": HOF100,
    "pea_isolate_40C_PratapSingh2021": NA(_PRATAP_NOTE),
    "pea_isolate_ribose_cysteine_100C_45min_Internal2026": NA(_SYNTHETIC_NOTE),
    "pea_isolate_uht_140C_Trikusuma2019": _block(None, None, "aseptic UHT line, 140 C / 6 s", "continuous_process", "unstated", None, "primary_source_pdf", _TRIKUSUMA_NOTE),
    "resconi_2023_pbma_beef_identity_benchmark": NA(_RESCONI_NOTE),
    "soy_isolate_40C_PratapSingh2021": NA(_PRATAP_NOTE),
    "soy_isolate_ribose_cysteine_100C_45min_Internal2026": NA(_SYNTHETIC_NOTE),
    "thiamine_cys_glucose_120C_Bolton1994": _block(33.3, 125.0, "125 mL glass vial sealed with Teflon-lined septum", "air_by_default", "distilled", None, "primary_source_pdf", _BOLTON_NOTE),
    # --- external matrix ---
    "external_validation_bi_2020_raw_pea_hexanal": NA(_BI_RAW_NOTE, "primary_source_pdf"),
    "external_validation_bi_2020_roasted_pea_hexanal": _block(None, None, "forced-air oven, 160 C / 30 min, 500 g peas", "open", "not_applicable", None, "primary_source_pdf", _BI_ROAST_NOTE),
    "external_validation_li_2026_spi_wg_hme_control": _block(None, None, "high-moisture extruder", "continuous_process", "unstated", None, "repo_verbatim_methods_quote", _LI_NOTE),
    "external_validation_liu_2023_ppi_offnote_baseline": NA(_LIU_NOTE, "primary_source_pdf"),
    # --- maillard_path hold-outs ---
    "mp_holdout_fructose_asparagine_180C_Lin2022": _block(100.0, None, "not stated (dry-bath incubator)", "unstated", "distilled", None, "repo_verbatim_methods_quote", _LIN_NOTE),
    "mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019": _block(None, None, "sealed 10 mL glass ampoule (2 mL) or 20 mL crimp-capped headspace vial (6 mL); thesis does not say which", "air_by_default", "distilled", None, "repo_verbatim_methods_quote", _SCHIBILSKY_NOTE),
    "mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019": _block(None, None, "sealed 10 mL glass ampoule (2 mL) or 20 mL crimp-capped headspace vial (6 mL); thesis does not say which", "air_by_default", "distilled", None, "repo_verbatim_methods_quote", _SCHIBILSKY_NOTE),
    "mp_holdout_glucose_asparagine_180C_10min_Chang2021": _block(100.0, None, "not stated", "unstated", "distilled", None, "repo_verbatim_methods_quote", _CHANG_NOTE),
    "mp_holdout_glucose_asparagine_180C_30min_Chang2021": _block(100.0, None, "not stated", "unstated", "distilled", None, "repo_verbatim_methods_quote", _CHANG_NOTE),
    "mp_holdout_glucose_asparagine_180C_30min_water_Chang2021": _block(100.0, None, "not stated", "unstated", "distilled", None, "repo_verbatim_methods_quote", _CHANG_NOTE),
    "mp_holdout_glucose_asparagine_180C_Ye2024": _block(4.0, 25.0, "25 mL thick-walled pressurized glass tube, sealed, oil bath", "air_by_default", "unstated", None, "repo_verbatim_methods_quote", _YE_NOTE),
    "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021": _block(None, None, "not stated (autoclaved)", "unstated", "unstated", None, "repo_verbatim_methods_quote", _STEINHAGEN_NOTE),
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3": HOF100,
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7": HOF100,
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3": HOF100,
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7": HOF100,
    "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5": HOF100,
    "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026": YIL,
    "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026": YIL,
    "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026": YIL,
    "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026": YIL,
}


def panel_files():
    yield from sorted(data_paths.BENCHMARKS_DIR.glob("*.json"))
    yield from sorted(data_paths.EXTERNAL_VALIDATION_DIR.glob("*.json"))
    yield from sorted(data_paths.MAILLARD_PATH_HOLDOUT_DIR.glob("*.json"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[3])
    parser.add_argument("--check", action="store_true", help="report drift, write nothing")
    args = parser.parse_args()
    drift = []
    seen = set()
    for path in panel_files():
        payload = json.loads(path.read_text(encoding="utf-8"))
        key = str(payload.get("benchmark_id") or path.stem)
        block = BLOCKS.get(key)
        if block is None:
            raise SystemExit(f"{key}: no vessel block declared in this script -- add one deliberately")
        seen.add(key)
        current = (payload.get("conditions") or {}).get("vessel")
        if current == block:
            continue
        drift.append(key)
        if not args.check:
            payload.setdefault("conditions", {})["vessel"] = block
            path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    unused = sorted(set(BLOCKS) - seen)
    if unused:
        raise SystemExit(f"declared blocks with no bundle on the panel: {unused}")
    if args.check:
        if drift:
            print(f"vessel blocks stale on {len(drift)} bundle(s): {drift}", file=sys.stderr)
            return 1
        print("vessel blocks current on every panel bundle")
        return 0
    print(f"wrote vessel blocks on {len(drift)} bundle(s); {len(seen) - len(drift)} already current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
