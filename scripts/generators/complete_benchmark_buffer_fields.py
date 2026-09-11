#!/usr/bin/env python3
"""
scripts/generators/complete_benchmark_buffer_fields.py

BUILD WAVE B2.3, TASK 2 -- THE BUFFER-FIELD COMPLETION OF THE 21 FROZEN
EXTERNAL-VALIDATION BUNDLES.

WHY THIS EXISTS
===============
`results/validation/kinetic_core_b2_2_diagnosis.md` sec. 4 probed the frozen
bundles' condition SCHEMA (keys only, no value read) and found it to be exactly
`{ph, temp_C, time_min, water_activity}`. There is no buffer field. The
engine's declared default for an undeclared buffer is *unbuffered plus an
extrapolation warning*, so bundles literally named
`..._ribose_cysteine_buffer_...` were being integrated as if they were water --
and the same diagnosis sized that schema gap at an **11x swing in predicted
2-furfurylthiol** on exactly the system shape whose exam row is the worst point
in the whole report.

`docs/reference/FIT_HOLDOUT_DECLARATION.md` **Amendment 9 clause 2** authorises
this edit and fences it:

  * the CONDITION records may be completed from each bundle's SOURCE PAPER,
    with per-field provenance;
  * measured values, compound lists, `evidence_class` and roles remain
    BYTE-IDENTICAL;
  * the exam is thereafter reported BOTH WAYS -- buffer-completed and as-was --
    in one artifact, permanently.

WHAT THIS SCRIPT MAY AND MAY NOT TOUCH
======================================
It writes exactly one new key per bundle: ``conditions.buffer``. It does not
edit `conditions.ph`, and specifically not where the source disagrees with it
(four such disagreements are recorded below as `ph_disagreement` INSIDE the
buffer block, reported and not acted on -- editing an executable condition
would change a prediction, which is a different licence from the one Amendment
9 grants). `tests/unit/test_kinetic_core_b2_3.py` proves the restriction by
hashing every bundle with `conditions.buffer` deleted and comparing against
baselines taken from git HEAD before this script first ran.

PROVENANCE CLASSES, AND WHY THE DISTINCTION IS KEPT VISIBLE
===========================================================
Only TWO of the eight source papers are on disk in `data/articles/`
(`hofmann1998.pdf`, `bi2020.pdf`), covering 7 of the 21 bundles. For the rest
the buffer is read from a VERBATIM METHODS QUOTATION already carried inside the
bundle's own `content_verification.quoted_method` -- the repo's own prior
extraction, with its retrieval route named. That is evidence, not a guess, but
it is second-hand and every field says so:

  * ``primary_source_pdf``            -- read from the PDF this wave opened.
  * ``repo_verbatim_methods_quote``   -- read from a verbatim quotation the
                                         repo already carries, source paper NOT
                                         on disk.
  * ``unknown``                       -- neither exists. Recorded as
                                         ``buffer_unknown`` with a note. NOT
                                         guessed, NOT inferred from the
                                         benchmark's own name (a bundle called
                                         "..._buffer_..." is not evidence of
                                         which buffer).

Run:  python scripts/generators/complete_benchmark_buffer_fields.py [--check]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

BUNDLE_DIR = data_paths.EXTERNAL_VALIDATION_DIR

# ===========================================================================
# THE COMPLETED BUFFER RECORDS
# ===========================================================================
# species:          a buffer identity, "none", or "buffer_unknown"
# concentration_M:  TOTAL buffer molarity, or null
# initial_ph_scale: "bench_cooled" (a room-temperature reading or a nominal
#                   buffer pH set at the bench, which is what every pH number
#                   in this corpus is) | "at_temperature" | "unstated"
# provenance_note:  what the source says, with a page/section anchor, and the
#                   provenance class.

_HOFMANN = dict(
    species="phosphate",
    concentration_M=0.5,
    initial_ph_scale="bench_cooled",
    provenance_class="primary_source_pdf",
    provenance_note=(
        "Hofmann & Schieberle 1998, JAFC 46:235, p. 236, EXPERIMENTAL "
        "PROCEDURES / 'Model Experiments', verbatim: 'the reactants (amounts "
        "detailed in the tables) were dissolved in phosphate buffer "
        "(0.5 mol/L) at pH 3.0, 5.0, or 7.0 and thermally treated for 20 min "
        "at 145 degC in a laboratory autoclave (200 mL; Type II; Roth, "
        "Germany)'. Read from data/articles/hofmann1998.pdf. The pH is the "
        "buffer's NOMINAL bench value; the paper does not measure it at "
        "temperature or after heating, and the bundle's own "
        "precursor_concentration_provenance.residual_caveats (b) already says "
        "so. NOT THIS ARM: the same paper's SECOND series dry-heats on silica "
        "gel wetted with 2 mol/L buffer -- 4x stronger -- and none of the five "
        "frozen Hofmann bundles uses it."
    ),
)

_YILTIRAK = dict(
    species="potassium_phosphate",
    concentration_M=0.5,
    initial_ph_scale="bench_cooled",
    provenance_class="repo_verbatim_methods_quote",
    provenance_note=(
        "Yiltirak et al. 2026, Food Res. Int. 231:118600, sec. 2.3 (model "
        "preparation), verbatim as carried in this bundle's own "
        "content_verification.quoted_method: 'All buffer phases were prepared "
        "by dissolving ribose (25 mM) and cysteine (25 mM) in potassium "
        "phosphate buffer (0.5 M, pH 5.5).' The scored arm is model (iv), "
        "'buffer system (100% buffer)' -- oil-free and emulsifier-free. "
        "SOURCE PAPER IS NOT ON DISK in data/articles/; this is the repo's own "
        "verbatim extraction, second-hand. TWO CAVEATS THE SOURCE ITSELF "
        "RAISES, carried here rather than dropped: the buffer was made in TAP "
        "WATER, not deionised water, so trace-metal catalysis is uncontrolled; "
        "and the paper does not report the post-heating pH."
    ),
)

_CHANG_ACETATE = dict(
    species="acetate",
    concentration_M=0.17,
    initial_ph_scale="bench_cooled",
    provenance_class="repo_verbatim_methods_quote",
    provenance_note=(
        "Chang / Sung et al. 2021, Polymers 13(12):1901, sec. 2.2, verbatim as "
        "carried in this bundle's content_verification.quoted_method: 'MRPs "
        "were prepared by combining 0.5 g asparagine, 0.5 g chitosan, and "
        "0.5 g glucose with 1% acetic acid' and 'The pH value of each solution "
        "was first adjusted to 5.8 by adding 1 N sodium hydroxide (NaOH) and "
        "then 0.001 N NaOH to pH 6.0 and topped up with distilled water to "
        "100 mL.' Acetic acid back-titrated with NaOH to pH 6.0 IS an acetate "
        "buffer, which is why this is recorded as a buffer and not as water. "
        "SOURCE PAPER NOT ON DISK; second-hand. TWO AMBIGUITIES, REPORTED NOT "
        "RESOLVED: (i) '1%' carries no w/v-vs-v/v basis, so the total acetate "
        "is 0.167-0.175 M and 0.17 is the midpoint; (ii) the quoted sentence "
        "ties the 1% acetic acid to a formulation that ALSO contains 0.5 g "
        "chitosan, while this bundle's vehicle_note labels the same arm merely "
        "'1% acetic acid, back-titrated to pH 6.0 (Figure 3B arm)'. Whether "
        "the scored arm contains chitosan cannot be settled from the quoted "
        "text."
    ),
)


BUFFER_RECORDS: Dict[str, Dict[str, Any]] = {
    # ---- source paper ON DISK -------------------------------------------
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3": _HOFMANN,
    "mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7": _HOFMANN,
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3": _HOFMANN,
    "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7": _HOFMANN,
    "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5": _HOFMANN,
    "external_validation_bi_2020_raw_pea_hexanal": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="unstated",
        provenance_class="primary_source_pdf",
        provenance_note=(
            "Bi et al. 2020, JAFC 68:2718, p. 2719, 'Isolation of Volatile "
            "Compounds from Raw and Roasted Pea Flour', verbatim: 'Sodium "
            "chloride (1.8 g) was dissolved in a 5 mL aliquot of distilled "
            "water contained in a 20 mL vial ... raw or roasted pea flour "
            "(1.5 g) was added.' p. 2719 'Chemicals': 'de-ionized water was "
            "used throughout the experiment.' NO BUFFER: distilled water plus "
            "NaCl, the salt being an SPME salting-out aid rather than a "
            "buffer. Read from data/articles/bi2020.pdf. "
            "ph_disagreement: THE PAPER STATES NO pH ANYWHERE (every page "
            "searched; zero hits). This bundle's conditions.ph = 6.0 is "
            "therefore an ASSUMPTION with no source. It is REPORTED and NOT "
            "EDITED -- changing an executable condition is outside Amendment "
            "9's licence and would change a prediction."
        ),
        ph_disagreement=(
            "conditions.ph = 6.0 is unsourced: Bi 2020 states no pH."
        ),
    ),
    "external_validation_bi_2020_roasted_pea_hexanal": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="unstated",
        provenance_class="primary_source_pdf",
        provenance_note=(
            "As the raw-pea bundle: Bi et al. 2020, JAFC 68:2718, p. 2719, "
            "distilled water + NaCl for the headspace isolation, and the "
            "roasting itself is DRY (forced-air oven, p. 2719 'Sample "
            "Preparation and Roasting Procedures'), so there is no aqueous "
            "phase to buffer during the reaction at all. Read from "
            "data/articles/bi2020.pdf. "
            "ph_disagreement: the paper states no pH anywhere; conditions.ph "
            "= 6.0 is an unsourced assumption, reported and not edited."
        ),
        ph_disagreement=(
            "conditions.ph = 6.0 is unsourced: Bi 2020 states no pH."
        ),
    ),
    # ---- verbatim methods quotation carried by the repo ------------------
    "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026": _YILTIRAK,
    "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026": _YILTIRAK,
    "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026": _YILTIRAK,
    "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026": _YILTIRAK,
    "mp_holdout_glucose_asparagine_180C_10min_Chang2021": _CHANG_ACETATE,
    "mp_holdout_glucose_asparagine_180C_30min_Chang2021": _CHANG_ACETATE,
    "mp_holdout_glucose_asparagine_180C_30min_water_Chang2021": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "Chang / Sung et al. 2021, Polymers 13(12):1901, sec. 2.2, "
            "verbatim: 'Solutions containing 0.5 g asparagine, 0.5 g glucose, "
            "and 0.5 g hydroxymethylfurfural in distilled water were also "
            "prepared', and the same pH ladder ('adjusted to 5.8 by adding 1 N "
            "NaOH and then 0.001 N NaOH to pH 6.0'). This bundle's own "
            "vehicle_note names it 'Deionized water, pH adjusted to 6.0 "
            "(Figure 3A arm)'. POSITIVE FINDING OF NO BUFFER, not an absence "
            "of information -- and it is the WATER half of the two-arm pair "
            "whose arms the exam currently predicts identically. SOURCE PAPER "
            "NOT ON DISK; second-hand."
        ),
    ),
    "mp_holdout_fructose_asparagine_180C_Lin2022": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "Lin et al. 2022, Polymers 14(8):1565, sec. 2.2 'Preparation of "
            "Maillard Reaction Model System', verbatim as carried in this "
            "bundle's content_verification.quoted_method: 'We adjusted the pH "
            "of each solution to 5.8 by adding 1 N sodium hydroxide (NaOH) and "
            "then to 6.0 by adding 0.001 N NaOH' and 'Each solution was topped "
            "up with distilled water to 100 mL.' NO BUFFER: distilled water, "
            "pH set at the bench with NaOH and free to drift over the hold. "
            "SOURCE PAPER NOT ON DISK; second-hand."
        ),
    ),
    "mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "Schibilsky 2019, TU Berlin dissertation, sec. 5.1.6 / 5.1.6.2, "
            "verbatim as carried in this bundle's "
            "content_verification.quoted_method: 'mit bidestilliertem Wasser "
            "aufgefuellt' and 'Des Weiteren wurde der pH-Wert auf 5,00 +/- "
            "0,01 und 8,00 +/- 0,01 eingestellt.' NO BUFFER: bidistilled "
            "water, pH set at t=0 with dilute NaOH/HCl. The bundle's own "
            "content_verification.quantification_verdict already states the "
            "consequence in capitals -- 'the system is UNBUFFERED - pH is set "
            "at t=0 ... and will drift over 2 h at 130 degC, and the thesis "
            "does not report the final pH'. SOURCE NOT ON DISK; second-hand."
        ),
    ),
    "mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "As the pH-5.0 arm: Schibilsky 2019, sec. 5.1.6, bidistilled "
            "water, UNBUFFERED, pH set at t=0 ('auf 5,00 +/- 0,01 und 8,00 "
            "+/- 0,01 eingestellt'). This is the pH-8.0 arm. SOURCE NOT ON "
            "DISK; second-hand."
        ),
    ),
    "mp_holdout_glucose_asparagine_180C_Ye2024": dict(
        species="phosphate",
        concentration_M=0.1,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "Ye et al. 2024, Foods 13(17):2836, sec. 2.2 'Glc-Asn Thermal "
            "Model Reaction', verbatim as carried in this bundle's "
            "content_verification.quoted_method: 'an equimolar solution of "
            "glucose and asparagine was accurately prepared in phosphate "
            "buffer (0.1 M, pH 6.86), and 4 mL of the solution was "
            "transferred to a 25 mL thick-walled pressurized glass tube ... at "
            "180 degC and kept for 30 min.' The 6.86 is the buffer's nominal "
            "bench pH. SOURCE PAPER NOT ON DISK; second-hand. NOTE, already "
            "recorded by the bundle: the same sentence is the ONLY statement "
            "of the medium and the paper never states the reactant molarity."
        ),
    ),
    "external_validation_liu_2023_ppi_offnote_baseline": dict(
        species="none",
        concentration_M=None,
        initial_ph_scale="bench_cooled",
        provenance_class="repo_verbatim_methods_quote",
        provenance_note=(
            "Liu, Cadwallader & Drake -- thesis sec. 2.9 / Food Chemistry 406: "
            "134998, as carried in "
            "data/benchmarks/external_validation/intake/"
            "external_validation_liu_2023_ppi_offnote_baseline.yaml "
            "(matrix_format): 'commercial pea protein rehydrated to 10% solids "
            "(w/w) in deionized water'. NO BUFFER. The slurry is never heated "
            "-- 40 C is a headspace equilibration, not a reaction -- so the "
            "buffer field is close to inert here in any case. SOURCE PAPER NOT "
            "ON DISK; second-hand. "
            "ph_disagreement: the source reports the rehydrated 10%-solids "
            "slurries at pH 6.3-7.3 (mean ~6.8) while conditions.ph = 6.0. The "
            "protocol YAML already self-declares this ('the proxied pH 6.0 is "
            "NOT the source pH ... left unchanged because editing an "
            "executable condition would change the prediction'). Reported, not "
            "edited."
        ),
        ph_disagreement=(
            "conditions.ph = 6.0 against a source range of 6.3-7.3 "
            "(mean ~6.8); already self-declared in the protocol YAML."
        ),
    ),
    # ---- neither a PDF nor a verbatim quotation: UNKNOWN ------------------
    "external_validation_li_2026_spi_wg_hme_control": dict(
        species="buffer_unknown",
        concentration_M=None,
        initial_ph_scale="unstated",
        provenance_class="unknown",
        provenance_note=(
            "Li et al. 2026, Foods 15(5):912. SOURCE PAPER NOT ON DISK and the "
            "repo carries NO verbatim methods quotation for the medium. The "
            "scored system is a 57 wt%-moisture high-moisture-extrusion melt, "
            "not an aqueous solution, so 'buffer' may not even be the right "
            "question -- but that is a reason to record UNKNOWN, not a licence "
            "to record 'none'. The protocol's own note is explicit: 'The paper "
            "does not publish a final-blend pH or water-activity closure' "
            "(data/benchmarks/external_validation/intake/"
            "external_validation_li_2026_spi_wg_hme_control.yaml, "
            "benchmark_alignment.notes). "
            "ph_disagreement: the only pH in the whole record is 7.0 for the "
            "WHEAT-GLUTEN ENZYMATIC PRETREATMENT at 30 C -- a different unit "
            "operation -- and it has been transposed onto a 160 C extrusion "
            "run. conditions.ph = 7.0 is unsourced AS AN EXTRUSION CONDITION. "
            "Reported, not edited."
        ),
        ph_disagreement=(
            "conditions.ph = 7.0 is the wheat-gluten pretreatment pH at 30 C, "
            "not a stated extrusion condition."
        ),
    ),
    "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021": dict(
        species="buffer_unknown",
        concentration_M=None,
        initial_ph_scale="at_temperature",
        provenance_class="unknown",
        provenance_note=(
            "Steinhagen et al. 2021, Pharmaceuticals 14(11):1121. SOURCE PAPER "
            "NOT ON DISK and the repo carries no verbatim statement of the "
            "medium beyond '10% (w/v) glucose solution ... autoclaved at 111, "
            "116, and 121 degC'. An unbuffered infusion solution is the "
            "LIKELY reading and it is exactly the reading this field refuses "
            "to record, because likely is not stated. "
            "ph_disagreement, AND IT IS THE ONE THAT BREAKS THE SCALE "
            "CONVENTION: conditions.ph = 4.36 is the pH measured AFTER "
            "autoclaving (the source's Table 3). The initial / non-autoclaved "
            "value is 4.98 (Table 6). Every other bundle's pH is an INITIAL "
            "bench reading, so this one is on a different scale from all "
            "twenty others -- which is why initial_ph_scale is recorded here "
            "as 'at_temperature' rather than 'bench_cooled'. The bundle "
            "discloses the choice honestly but a model that treats 4.36 as an "
            "initial condition is being handed the answer to part of its own "
            "trajectory. Reported, not edited."
        ),
        ph_disagreement=(
            "conditions.ph = 4.36 is the POST-autoclave measurement; the "
            "initial value is 4.98. This bundle's pH is on a different scale "
            "from the other twenty."
        ),
    ),
}


# ---------------------------------------------------------------------------
# WAVE B36 (2026-09-11): seven sources arrived. Five buffer notes above open with a
# "SOURCE PAPER NOT ON DISK; second-hand" claim that was true when written and is false now.
# The repo's standing practice for a corrected note (B34, B35): the new reading first, the
# old claim kept after it, labelled superseded. Applied here so the records above stay the
# audit record as written.
# ---------------------------------------------------------------------------
_B36_SUPERSEDED = " ||| PRIOR NOTE, SUPERSEDED 2026-09-11 AND RETAINED AS THE AUDIT RECORD: "
_B36_ON_DISK = {
    "mp_holdout_glucose_asparagine_180C_10min_Chang2021": ("primary_source_pdf",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Lin, Chan, Kao & Sung 2021 (Polymers 13:1901), "
        "data/articles/lin2021.pdf, dossier lin2021_extraction.md. The quotation below is verbatim in sec. 2.2 "
        "of the print; the print says '1% acetic acid' and no more, so the w/v-versus-v/v ambiguity stands."),
    "mp_holdout_glucose_asparagine_180C_30min_Chang2021": ("primary_source_pdf",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Lin, Chan, Kao & Sung 2021 (Polymers 13:1901), "
        "data/articles/lin2021.pdf, dossier lin2021_extraction.md. The quotation below is verbatim in sec. 2.2 "
        "of the print; the print says '1% acetic acid' and no more, so the w/v-versus-v/v ambiguity stands."),
    "mp_holdout_glucose_asparagine_180C_30min_water_Chang2021": ("primary_source_pdf",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Lin, Chan, Kao & Sung 2021 (Polymers 13:1901), "
        "data/articles/lin2021.pdf, dossier lin2021_extraction.md. The distilled-water arm and its pH ladder "
        "are verbatim in sec. 2.2 of the print: a positive finding of no buffer, confirmed."),
    "mp_holdout_fructose_asparagine_180C_Lin2022": ("primary_source_pdf",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Lin, Ting, Ndraha, Hsiao & Sung 2022 (Polymers 14:1565), "
        "data/articles/lin2022.pdf, dossier lin2022_extraction.md. The pH ladder and 'topped up with distilled "
        "water to 100 mL' are verbatim in sec. 2.2 of the print: no buffer, confirmed."),
    "mp_holdout_glucose_asparagine_180C_Ye2024": ("primary_source_pdf",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Ye et al. 2024 (Foods 13:2836), data/articles/ye2024.pdf, "
        "dossier ye2024_extraction.md. The 0.1 M phosphate buffer at pH 6.86 is verbatim in sec. 2.2 of the print, "
        "which still never states the reactant molarity."),
    "external_validation_li_2026_spi_wg_hme_control": ("unknown",
        "WAVE B36 (2026-09-11): THE SOURCE IS ON DISK -- Li, Dai, Mao, An, Bai & Kaur 2026 (Foods 15:912), "
        "data/articles/li2026.pdf, dossier li2026_extraction.md -- AND THE PRINT DOES NOT STATE THE MEDIUM: sec. 2.3 "
        "gives the blend ratio, moisture (about 57 %), screw speed, feed rate and the six zone temperatures, and no "
        "pH or buffer for the extrusion blend. buffer_unknown stays, now because the paper is silent, not absent."),
}


# ---------------------------------------------------------------------------
# The records above were written once (Amendment 9, 2026-09-04) and later waves edited SEVEN
# bundles' buffer blocks in place rather than here: B34 (Steinhagen), B35 (liu_2023), the
# Yiltirak primary reading (four bundles) and li_2026. Re-running this script on 2026-09-11
# silently regressed all seven, which wave B36 caught before committing. The blocks as they
# stand on disk are carried here verbatim, field by field, so that the script reproduces the
# repository and ``--check`` reports drift instead of pretending there is none. Any future edit
# to a buffer block belongs HERE, not in the bundle.
# ---------------------------------------------------------------------------
_EDITED_IN_PLACE: Dict[str, Dict[str, Any]] = {'external_validation_li_2026_spi_wg_hme_control': {'provenance_note': 'Li et al. 2026, Foods '
                                                                       '15(5):912. SOURCE PAPER '
                                                                       'NOT ON DISK and the repo '
                                                                       'carries NO verbatim '
                                                                       'methods quotation for the '
                                                                       'medium. The scored system '
                                                                       'is a 57 wt%-moisture '
                                                                       'high-moisture-extrusion '
                                                                       'melt, not an aqueous '
                                                                       "solution, so 'buffer' may "
                                                                       'not even be the right '
                                                                       'question -- but that is a '
                                                                       'reason to record UNKNOWN, '
                                                                       'not a licence to record '
                                                                       "'none'. The protocol's own "
                                                                       "note is explicit: 'The "
                                                                       'paper does not publish a '
                                                                       'final-blend pH or '
                                                                       "water-activity closure' "
                                                                       '(data/protocols/external_validation/external_validation_li_2026_spi_wg_hme_control.yaml, '
                                                                       'benchmark_alignment.notes). '
                                                                       'ph_disagreement: the only '
                                                                       'pH in the whole record is '
                                                                       '7.0 for the WHEAT-GLUTEN '
                                                                       'ENZYMATIC PRETREATMENT at '
                                                                       '30 C -- a different unit '
                                                                       'operation -- and it has '
                                                                       'been transposed onto a 160 '
                                                                       'C extrusion run. '
                                                                       'conditions.ph = 7.0 is '
                                                                       'unsourced AS AN EXTRUSION '
                                                                       'CONDITION. Reported, not '
                                                                       'edited.'},
 'external_validation_liu_2023_ppi_offnote_baseline': {'provenance_note': 'WAVE B35 (2026-09-11): '
                                                                          'THE SOURCE IS ON DISK '
                                                                          '-- '
                                                                          'data/articles/liu2023.pdf, '
                                                                          '717 KB, read '
                                                                          '2026-09-04, as this '
                                                                          "bundle's OWN vessel "
                                                                          'note already said. This '
                                                                          'note and that one '
                                                                          'disagreed: the vessel '
                                                                          'block was updated when '
                                                                          'the paper arrived and '
                                                                          'this one was not. ||| '
                                                                          'PRIOR NOTE, SUPERSEDED '
                                                                          'AND RETAINED AS THE '
                                                                          'AUDIT RECORD, and its '
                                                                          "'SOURCE PAPER NOT ON "
                                                                          "DISK; second-hand' was "
                                                                          'already false when the '
                                                                          'vessel note was '
                                                                          'written: Liu, '
                                                                          'Cadwallader & Drake -- '
                                                                          'thesis sec. 2.9 / Food '
                                                                          'Chemistry 406: 134998, '
                                                                          'as carried in '
                                                                          'data/protocols/external_validation/external_validation_liu_2023_ppi_offnote_baseline.yaml '
                                                                          '(matrix_format): '
                                                                          "'commercial pea protein "
                                                                          'rehydrated to 10% '
                                                                          'solids (w/w) in '
                                                                          "deionized water'. NO "
                                                                          'BUFFER. The slurry is '
                                                                          'never heated -- 40 C is '
                                                                          'a headspace '
                                                                          'equilibration, not a '
                                                                          'reaction -- so the '
                                                                          'buffer field is close '
                                                                          'to inert here in any '
                                                                          'case. SOURCE PAPER NOT '
                                                                          'ON DISK; second-hand. '
                                                                          'ph_disagreement: the '
                                                                          'source reports the '
                                                                          'rehydrated 10%-solids '
                                                                          'slurries at pH 6.3-7.3 '
                                                                          '(mean ~6.8) while '
                                                                          'conditions.ph = 6.0. '
                                                                          'The protocol YAML '
                                                                          'already self-declares '
                                                                          "this ('the proxied pH "
                                                                          '6.0 is NOT the source '
                                                                          'pH ... left unchanged '
                                                                          'because editing an '
                                                                          'executable condition '
                                                                          'would change the '
                                                                          "prediction'). Reported, "
                                                                          'not edited.'},
 'mp_holdout_glucose_only_autoclave_121C_Steinhagen2021': {'ph_disagreement': 'conditions.ph = '
                                                                              '4.36 sits below '
                                                                              "Leitzen's printed "
                                                                              '4.4-5.2 '
                                                                              'post-autoclave '
                                                                              'range (4.98 '
                                                                              'before). Reported, '
                                                                              'not edited.',
                                                           'provenance_note': 'WAVE B34 '
                                                                              '(2026-09-11): THE '
                                                                              'SOURCE IS ON DISK '
                                                                              'and the previous '
                                                                              'note saying '
                                                                              'otherwise was '
                                                                              'wrong. Leitzen et '
                                                                              'al. 2021 '
                                                                              '(Pharmaceuticals '
                                                                              '14:1121, '
                                                                              'data/articles/Leitzen2021.pdf), '
                                                                              'Methods: 10 % (w/v) '
                                                                              'glucose in water in '
                                                                              'PP bottles, '
                                                                              'autoclaved; NO '
                                                                              'BUFFER is added, '
                                                                              'and the paper '
                                                                              'measures the pH '
                                                                              'itself -- 4.98 '
                                                                              'before autoclaving '
                                                                              'and 4.4-5.2 after '
                                                                              '(scheme A). This '
                                                                              "bundle's "
                                                                              'conditions.ph of '
                                                                              '4.36 is BELOW the '
                                                                              'range the paper '
                                                                              'prints and is left '
                                                                              'unedited, because '
                                                                              'changing an '
                                                                              'executable '
                                                                              'condition is '
                                                                              'outside the licence '
                                                                              'an extraction '
                                                                              'correction carries; '
                                                                              'it is reported here '
                                                                              'instead. ||| PRIOR '
                                                                              'NOTE, SUPERSEDED '
                                                                              '2026-09-11 AND '
                                                                              'RETAINED AS THE '
                                                                              'AUDIT RECORD -- its '
                                                                              'claim that the '
                                                                              'source is not on '
                                                                              'disk was false when '
                                                                              'written and is '
                                                                              'false now: '
                                                                              'Steinhagen et al. '
                                                                              '2021, '
                                                                              'Pharmaceuticals '
                                                                              '14(11):1121. SOURCE '
                                                                              'PAPER NOT ON DISK '
                                                                              'and the repo '
                                                                              'carries no verbatim '
                                                                              'statement of the '
                                                                              "medium beyond '10% "
                                                                              '(w/v) glucose '
                                                                              'solution ... '
                                                                              'autoclaved at 111, '
                                                                              "116, and 121 degC'. "
                                                                              'An unbuffered '
                                                                              'infusion solution '
                                                                              'is the LIKELY '
                                                                              'reading and it is '
                                                                              'exactly the reading '
                                                                              'this field refuses '
                                                                              'to record, because '
                                                                              'likely is not '
                                                                              'stated. '
                                                                              'ph_disagreement, '
                                                                              'AND IT IS THE ONE '
                                                                              'THAT BRE'},
 'mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026': {'provenance_class': 'primary_source_pdf',
                                                            'provenance_note': 'Yiltirak et al. '
                                                                               '2026, Food Res. '
                                                                               'Int. 231:118600, '
                                                                               'read from '
                                                                               'data/articles/Yiltirak2026.pdf '
                                                                               'on 2026-09-06. '
                                                                               'Sec. 2.1 (p. 2), '
                                                                               'verbatim: '
                                                                               "'Potassium "
                                                                               'phosphate buffer '
                                                                               '(0.5 M, pH 5.5) '
                                                                               'was made using '
                                                                               'potassium '
                                                                               'dihydrogen '
                                                                               'phosphate and '
                                                                               'dipotassium '
                                                                               'hydrogen phosphate '
                                                                               "in tap water.' "
                                                                               'Sec. 2.3 (pp. '
                                                                               '2-3), verbatim: '
                                                                               "'All buffer phases "
                                                                               'were prepared by '
                                                                               'dissolving ribose '
                                                                               '(25 mM) and '
                                                                               'cysteine (25 mM) '
                                                                               'in potassium '
                                                                               'phosphate buffer '
                                                                               "(0.5 M, pH 5.5).' "
                                                                               'The scored arm is '
                                                                               'model (iv), '
                                                                               "'buffer system "
                                                                               "(100% buffer)' -- "
                                                                               'oil-free and '
                                                                               'emulsifier-free. '
                                                                               'TWO CAVEATS THE '
                                                                               'SOURCE ITSELF '
                                                                               'RAISES, carried '
                                                                               'rather than '
                                                                               'dropped: the '
                                                                               'buffer was made in '
                                                                               'TAP WATER, not '
                                                                               'deionised water, '
                                                                               'so trace-metal '
                                                                               'catalysis is '
                                                                               'uncontrolled; and '
                                                                               'the paper reports '
                                                                               'pH before heating '
                                                                               'only (0.5 M '
                                                                               'phosphate; the '
                                                                               'drift over the '
                                                                               'cook is not '
                                                                               'measured). Until '
                                                                               '2026-09-06 this '
                                                                               'block was a '
                                                                               'second-hand '
                                                                               'reading of the '
                                                                               "repo's own "
                                                                               'verbatim quote; '
                                                                               'the PDF now on '
                                                                               'disk confirms it '
                                                                               'word for word.'},
 'mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026': {'provenance_class': 'primary_source_pdf',
                                                            'provenance_note': 'Yiltirak et al. '
                                                                               '2026, Food Res. '
                                                                               'Int. 231:118600, '
                                                                               'read from '
                                                                               'data/articles/Yiltirak2026.pdf '
                                                                               'on 2026-09-06. '
                                                                               'Sec. 2.1 (p. 2), '
                                                                               'verbatim: '
                                                                               "'Potassium "
                                                                               'phosphate buffer '
                                                                               '(0.5 M, pH 5.5) '
                                                                               'was made using '
                                                                               'potassium '
                                                                               'dihydrogen '
                                                                               'phosphate and '
                                                                               'dipotassium '
                                                                               'hydrogen phosphate '
                                                                               "in tap water.' "
                                                                               'Sec. 2.3 (pp. '
                                                                               '2-3), verbatim: '
                                                                               "'All buffer phases "
                                                                               'were prepared by '
                                                                               'dissolving ribose '
                                                                               '(25 mM) and '
                                                                               'cysteine (25 mM) '
                                                                               'in potassium '
                                                                               'phosphate buffer '
                                                                               "(0.5 M, pH 5.5).' "
                                                                               'The scored arm is '
                                                                               'model (iv), '
                                                                               "'buffer system "
                                                                               "(100% buffer)' -- "
                                                                               'oil-free and '
                                                                               'emulsifier-free. '
                                                                               'TWO CAVEATS THE '
                                                                               'SOURCE ITSELF '
                                                                               'RAISES, carried '
                                                                               'rather than '
                                                                               'dropped: the '
                                                                               'buffer was made in '
                                                                               'TAP WATER, not '
                                                                               'deionised water, '
                                                                               'so trace-metal '
                                                                               'catalysis is '
                                                                               'uncontrolled; and '
                                                                               'the paper reports '
                                                                               'pH before heating '
                                                                               'only (0.5 M '
                                                                               'phosphate; the '
                                                                               'drift over the '
                                                                               'cook is not '
                                                                               'measured). Until '
                                                                               '2026-09-06 this '
                                                                               'block was a '
                                                                               'second-hand '
                                                                               'reading of the '
                                                                               "repo's own "
                                                                               'verbatim quote; '
                                                                               'the PDF now on '
                                                                               'disk confirms it '
                                                                               'word for word.'},
 'mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026': {'provenance_class': 'primary_source_pdf',
                                                            'provenance_note': 'Yiltirak et al. '
                                                                               '2026, Food Res. '
                                                                               'Int. 231:118600, '
                                                                               'read from '
                                                                               'data/articles/Yiltirak2026.pdf '
                                                                               'on 2026-09-06. '
                                                                               'Sec. 2.1 (p. 2), '
                                                                               'verbatim: '
                                                                               "'Potassium "
                                                                               'phosphate buffer '
                                                                               '(0.5 M, pH 5.5) '
                                                                               'was made using '
                                                                               'potassium '
                                                                               'dihydrogen '
                                                                               'phosphate and '
                                                                               'dipotassium '
                                                                               'hydrogen phosphate '
                                                                               "in tap water.' "
                                                                               'Sec. 2.3 (pp. '
                                                                               '2-3), verbatim: '
                                                                               "'All buffer phases "
                                                                               'were prepared by '
                                                                               'dissolving ribose '
                                                                               '(25 mM) and '
                                                                               'cysteine (25 mM) '
                                                                               'in potassium '
                                                                               'phosphate buffer '
                                                                               "(0.5 M, pH 5.5).' "
                                                                               'The scored arm is '
                                                                               'model (iv), '
                                                                               "'buffer system "
                                                                               "(100% buffer)' -- "
                                                                               'oil-free and '
                                                                               'emulsifier-free. '
                                                                               'TWO CAVEATS THE '
                                                                               'SOURCE ITSELF '
                                                                               'RAISES, carried '
                                                                               'rather than '
                                                                               'dropped: the '
                                                                               'buffer was made in '
                                                                               'TAP WATER, not '
                                                                               'deionised water, '
                                                                               'so trace-metal '
                                                                               'catalysis is '
                                                                               'uncontrolled; and '
                                                                               'the paper reports '
                                                                               'pH before heating '
                                                                               'only (0.5 M '
                                                                               'phosphate; the '
                                                                               'drift over the '
                                                                               'cook is not '
                                                                               'measured). Until '
                                                                               '2026-09-06 this '
                                                                               'block was a '
                                                                               'second-hand '
                                                                               'reading of the '
                                                                               "repo's own "
                                                                               'verbatim quote; '
                                                                               'the PDF now on '
                                                                               'disk confirms it '
                                                                               'word for word.'},
 'mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026': {'provenance_class': 'primary_source_pdf',
                                                               'provenance_note': 'Yiltirak et al. '
                                                                                  '2026, Food Res. '
                                                                                  'Int. '
                                                                                  '231:118600, '
                                                                                  'read from '
                                                                                  'data/articles/Yiltirak2026.pdf '
                                                                                  'on 2026-09-06. '
                                                                                  'Sec. 2.1 (p. '
                                                                                  '2), verbatim: '
                                                                                  "'Potassium "
                                                                                  'phosphate '
                                                                                  'buffer (0.5 M, '
                                                                                  'pH 5.5) was '
                                                                                  'made using '
                                                                                  'potassium '
                                                                                  'dihydrogen '
                                                                                  'phosphate and '
                                                                                  'dipotassium '
                                                                                  'hydrogen '
                                                                                  'phosphate in '
                                                                                  "tap water.' "
                                                                                  'Sec. 2.3 (pp. '
                                                                                  '2-3), verbatim: '
                                                                                  "'All buffer "
                                                                                  'phases were '
                                                                                  'prepared by '
                                                                                  'dissolving '
                                                                                  'ribose (25 mM) '
                                                                                  'and cysteine '
                                                                                  '(25 mM) in '
                                                                                  'potassium '
                                                                                  'phosphate '
                                                                                  'buffer (0.5 M, '
                                                                                  "pH 5.5).' The "
                                                                                  'scored arm is '
                                                                                  'model (iv), '
                                                                                  "'buffer system "
                                                                                  "(100% buffer)' "
                                                                                  '-- oil-free and '
                                                                                  'emulsifier-free. '
                                                                                  'TWO CAVEATS THE '
                                                                                  'SOURCE ITSELF '
                                                                                  'RAISES, carried '
                                                                                  'rather than '
                                                                                  'dropped: the '
                                                                                  'buffer was made '
                                                                                  'in TAP WATER, '
                                                                                  'not deionised '
                                                                                  'water, so '
                                                                                  'trace-metal '
                                                                                  'catalysis is '
                                                                                  'uncontrolled; '
                                                                                  'and the paper '
                                                                                  'reports pH '
                                                                                  'before heating '
                                                                                  'only (0.5 M '
                                                                                  'phosphate; the '
                                                                                  'drift over the '
                                                                                  'cook is not '
                                                                                  'measured). '
                                                                                  'Until '
                                                                                  '2026-09-06 this '
                                                                                  'block was a '
                                                                                  'second-hand '
                                                                                  'reading of the '
                                                                                  "repo's own "
                                                                                  'verbatim quote; '
                                                                                  'the PDF now on '
                                                                                  'disk confirms '
                                                                                  'it word for '
                                                                                  'word.'}}


# ---------------------------------------------------------------------------
# WAVE B45 (2026-09-11): the tap-water caveat gets a floor and a consequence.
#
# These four bundles have said since 2026-09-06 that "the buffer was made in TAP WATER,
# not deionised water, so trace-metal catalysis is uncontrolled". That was the right
# thing to write and it was unquantified: nobody knew whether it mattered at 1 % or at
# 100 %. Baldus, Klie, De & Methner 2017 (JAFC, 10.1021/acs.jafc.6b05472,
# baldus2017_extraction.md) puts a floor under it, and the floor is low.
#
# They measured by ICP-OES that a 0.1 M acetate buffer made with ULTRAPURE Milli-Q water
# held every transition metal below the 0.08 uM detection limit -- and that adding 300 uM
# of highest-purity-grade cysteine RAISED copper to 0.26 uM, as a reagent impurity. That
# 0.26 uM sufficed to drive a statistically resolved thioether oxidation at 95 C / 180 min
# (DMS 15.32 +/- 0.15 uM with 13 uM cysteine, against a 16.73 +/- 0.11 uM control), and
# EDTA in MOLAR EXCESS OVER THE THIOL abolished it entirely (16.80 +/- 0.07, same Tukey
# group as the control). In the cysteine-plus-copper system itself, "Cys was almost
# completely degraded in 5 minutes during heating from 40-60 C".
#
# These are the only four bundles in the corpus with water_source: tap, and they are
# hold-out rows the thiol sink fails on. B38 found that sink unreachable by any refit --
# barrier and both dimerisation rates already on their ceilings -- and B45's probe P2
# measured the shape of the gap: the engine retains 99.49 % of charged cysteine after
# 5 min at 95 C and 83 % after three hours. The missing quantity is a CATALYTIC CHANNEL,
# not a larger constant.
#
# NOTHING IS EDITED AND NO TOLERANCE MOVES. There is no rate constant to adopt -- Baldus
# state outright that "there are no rate constants available" -- and inventing one to close
# the gap is the failure this repository exists to avoid. The clause below records what is
# uncontrolled, what sizes it, and what a replacement experiment must contain.
# ---------------------------------------------------------------------------
_B45_TAP_WATER_FLOOR = (
    " ||| WAVE B45 (2026-09-11), A FLOOR UNDER THE TAP-WATER CAVEAT ABOVE: Baldus, Klie, De & "
    "Methner 2017 (J. Agric. Food Chem., 10.1021/acs.jafc.6b05472; dossier baldus2017_extraction.md) "
    "measured by ICP-OES that a 0.1 M acetate buffer in ULTRAPURE Milli-Q water held every transition "
    "metal below the 0.08 uM detection limit, and that adding 300 uM of highest-purity-grade cysteine "
    "raised copper to 0.26 uM AS A REAGENT IMPURITY. That trace drove statistically resolved thioether "
    "oxidation at 95 C / 180 min (DMS 15.32 +/- 0.15 uM with 13 uM cysteine vs a 16.73 +/- 0.11 uM "
    "control, n=3, Tukey-Kramer HSD), and EDTA IN MOLAR EXCESS OVER THE THIOL abolished it (16.80 +/- "
    "0.07, the control's group). In their cysteine-plus-copper system, 'Cys was almost completely "
    "degraded in 5 minutes during heating from 40-60 C'. A FORTIORI, tap water is not a controlled "
    "medium for a thiol pot. THE SIZE OF THE EFFECT IN THIS POT IS UNKNOWN AND IS NOT ESTIMATED: this "
    "bundle is NOT edited, its values are NOT adjusted and its tolerance is NOT widened. What this "
    "clause records is that the row's thiol trajectory has an uncontrolled catalytic term the model "
    "does not carry -- B38 found the thiol sink unreachable by refit (barrier and both dimerisation "
    "rates on their ceilings) and B45's probe P2 measured the gap's shape (the engine keeps 99.49 % of "
    "charged cysteine after 5 min at 95 C, 83 % after 3 h). A REPLACEMENT EXPERIMENT MUST BE RUN IN "
    "WATER OF STATED PROVENANCE, PAIRED WITH AND WITHOUT A CHELATOR IN MOLAR EXCESS OVER THE THIOL, "
    "WITH THE DISULFIDE QUANTIFIED IN THE SAME RUN; Baldus show that a chelator at 1:1 with the metal "
    "and sub-stoichiometric to the thiol makes matters WORSE, not better (DMS loss ~40 % -> ~72 %), so "
    "an unchecked 'add EDTA' control can invert the intended comparison."
)
_B45_TAP_WATER = (
    "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026",
    "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026",
)


def current_record(key: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """The buffer block a bundle must carry: base record, in-place edits, the B36 reading, the B45 clause."""
    out = dict(record)
    out.update(_EDITED_IN_PLACE.get(key, {}))
    hit = _B36_ON_DISK.get(key)
    if hit is not None:
        pclass, prefix = hit
        out["provenance_class"] = pclass
        out["provenance_note"] = prefix + _B36_SUPERSEDED + out["provenance_note"]
    if key in _B45_TAP_WATER:
        out["provenance_note"] = out["provenance_note"] + _B45_TAP_WATER_FLOOR
    return out


def bundle_paths():
    return sorted(BUNDLE_DIR.rglob("*.json"))


def stable_hash_without_buffer(payload: Dict[str, Any]) -> str:
    """
    SHA-256 of the whole bundle with ``conditions.buffer`` removed.

    Hashing EVERYTHING-EXCEPT is stronger than hashing the measured block: it
    proves not only that the measured values are untouched but that nothing
    else moved either -- roles, compound lists, evidence_class, provenance,
    the other three condition fields.
    """
    clone = json.loads(json.dumps(payload))
    conditions = clone.get("conditions")
    if isinstance(conditions, dict):
        conditions.pop("buffer", None)
    return hashlib.sha256(
        json.dumps(clone, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check", action="store_true",
        help="report what would change and write nothing")
    parser.add_argument(
        "--emit-baselines", action="store_true",
        help="print the everything-except-buffer hashes, for the unit test")
    args = parser.parse_args()

    paths = bundle_paths()
    if len(paths) != 21:
        raise SystemExit(f"expected 21 frozen bundles, found {len(paths)}")

    summary = []
    baselines = {}
    drift = []
    for path in paths:
        payload = json.loads(path.read_text())
        key = payload.get("benchmark_id") or path.stem
        record = BUFFER_RECORDS.get(key)
        if record is None:
            raise SystemExit(
                f"{path.name}: no buffer record declared for benchmark_id "
                f"{key!r}. Every one of the 21 bundles must be decided, "
                f"including as buffer_unknown -- silence is not an option."
            )
        before = stable_hash_without_buffer(payload)
        conditions = payload.setdefault("conditions", {})
        wanted = current_record(key, record)
        if conditions.get("buffer") != wanted:
            drift.append(key)
        conditions["buffer"] = wanted
        after = stable_hash_without_buffer(payload)
        if before != after:
            raise SystemExit(
                f"{path.name}: the edit changed something other than "
                f"conditions.buffer. Refusing to write."
            )
        baselines[key] = before
        summary.append((key, record["species"], record["concentration_M"],
                        record["provenance_class"],
                        "ph_disagreement" in record))
        if not args.check:
            path.write_text(json.dumps(payload, indent=2) + "\n")

    if args.emit_baselines:
        print(json.dumps(baselines, indent=4, sort_keys=True))
        return 0

    print(f"{'benchmark_id':62s} {'buffer':22s} {'M':>7s}  provenance")
    for key, species, molarity, klass, disagreement in summary:
        print(f"{key[:62]:62s} {species:22s} "
              f"{'--' if molarity is None else molarity:>7} "
              f" {klass}{'  [pH DISAGREEMENT]' if disagreement else ''}")
    current = {k: current_record(k, r) for k, r in BUFFER_RECORDS.items()}
    n_unknown = sum(1 for r in current.values()
                    if r["species"] == "buffer_unknown")
    n_none = sum(1 for r in current.values() if r["species"] == "none")
    n_primary = sum(1 for r in current.values()
                    if r["provenance_class"] == "primary_source_pdf")
    n_disagree = sum(1 for r in current.values() if "ph_disagreement" in r)
    print(f"\n21 bundles: {n_primary} from a PDF on disk, "
          f"{21 - n_primary - n_unknown} from a repo verbatim quotation, "
          f"{n_unknown} buffer_unknown. "
          f"{n_none} state NO BUFFER positively. "
          f"{n_disagree} carry a pH disagreement with their own source "
          f"(reported, not edited).")
    if args.check:
        if drift:
            print("--check: buffer blocks DRIFT from this script on: " + ", ".join(drift))
            return 1
        print("--check: buffer blocks current on every bundle; nothing was written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
