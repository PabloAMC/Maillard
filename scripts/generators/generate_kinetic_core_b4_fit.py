#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b4_fit.py

THE FIT REPORT OF BUILD WAVE B4 (matrix / OAV output layer).

Fits what is fittable on FIT ROWS ONLY, freezes it, and -- in the same pass and
BEFORE any hold-out value is read -- writes the PRE-REGISTERED BLIND
PREDICTIONS for the Hong 2020 flagship, together with the failures this wave
expects and the reasons it expects them.

WHAT IS FITTED
--------------
Two things, and both are deterministic functions of the frozen registry rather
than optimiser output, so "freezing the fit" and "freezing the registry" are the
same act:

  1. per-CLASS per-gram reversible binding constants, pooled WITHIN METHOD
     FAMILY from the Meynier partition RATIOS and the Leksrisompong K RATIOS
     (Amendment 4 FIT) and the Damodaran / Andriot constants (D.6 Module 6 FIT);
  2. the alpha,beta-unsaturation penalty, from the Vega gelatin ladder (D.6
     Module 7 FIT) and Meynier's skim-milk contrast (Amendment 4 FIT).

WHAT IS NOT FITTED, AND WHY
---------------------------
  * Nothing is fitted to Brewer 1995 (HOLD-OUT, and reclassified
    `dose_added_pre_cook`), to Hong 2020 (GATING HOLD-OUT), to Leksrisompong's
    BETs (not fit-eligible under Amendment 4), or to Barallat-Perez's lupin and
    mucin constants (Module 6 STAR HOLD-OUT).
  * No absolute partition coefficient is fitted. Meynier's own measured air/water
    K sits 6.24x below the Henry values in its own Table I, systematically, and
    Leksrisompong lands 6-17x low for the same methodological reason. The offset
    cancels in a within-study ratio, which is exactly why Amendment 4 makes the
    RATIOS fit-eligible and the absolutes never.
  * No general matrix correction factor is fitted, because k2 sec. D.1 computed
    what one would cost.

Writes results/validation/kinetic_core_b4_fit_report.{json,md}.
Run with: ./scripts/docker_maillard.sh run \
  "python scripts/generators/generate_kinetic_core_b4_fit.py"
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src import data_paths  # noqa: E402
from src.kinetic_core.matrix_oav import (  # noqa: E402
    MATRIX_THRESHOLDS,
    SEALED_OR_REFUSED_MATRICES,
    UNIMPLEMENTED_CANDIDATE_TERMS,
    WATER_THRESHOLDS,
    absolute_concentration,
    compare_formulations,
    covalent_channel_state,
    fit_class_binding_constants,
    fit_unsaturation_penalty,
    layer_metadata,
    oav_table,
    predict_matrix_shift,
    select_threshold,
)
from src.kinetic_core.parameters_matrix import (  # noqa: E402
    COMPOUND_STRUCTURE,
    MATRIX_LOADING,
    SOURCE_CONTRADICTIONS,
    matrix_registry_metadata,
)

OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b4_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b4_fit_report.md"
FROZEN_PREDICTIONS = data_paths.VALIDATION_DIR / "kinetic_core_b4_frozen_predictions.json"

# The Hong 2020 panel, in the order the public manifest lists it. The manifest
# is firewalled: it carries structure only, no threshold, ratio, sign or
# direction. These ten keys are the ONLY thing this script takes from it.
HONG_PANEL = (
    "2_3_dimethylpyrazine",
    "ethyl_4_methylpentanoate",
    "2_pentylfuran",
    "4_vinylphenol",
    "hexanal",
    "3_methylbutanal",
    "2_methylbutanal",
    "butyric_acid",
    "4_ethylphenol",
    "dimethyl_disulfide",
)

# =========================================================================
# PRE-REGISTERED EXPECTED FAILURES -- written down BEFORE any score exists
# =========================================================================
PRE_REGISTERED = {
    "gating_criterion": (
        "Declaration Amendment 4: >=7/10 of the paired soy/water ratios within "
        "5x AND correct sign on all 10, including the ethyl-4-methylpentanoate "
        "inversion."),
    "overall_expectation": "FAIL, and by construction rather than by accident.",
    "expected_failures": [
        {
            "compound": "ethyl_4_methylpentanoate",
            "expect": "SIGN WRONG",
            "confidence": "very high (~95%)",
            "why": (
                "The declaration names this compound as carrying an INVERSION. "
                "The only term this layer has for an ester is Meynier's measured "
                "ester SUPPRESSION (1.07-1.33x at 33.9 g/L dairy protein), which "
                "can only push the ratio ABOVE 1. Nothing in the FIT column "
                "produces a sign reversal for an ester: the two measured "
                "enhancements in the corpus are a lactone and a furanone, not an "
                "ester. The one mechanism that could -- background-odour masking "
                "-- is unimplementable, because Hong measured no background "
                "volatiles in the soy matrix. So this single compound fails the "
                "'correct sign on all 10' clause on its own, and therefore the "
                "gate fails whatever the magnitudes do."),
        },
        {
            "compounds": ["2_3_dimethylpyrazine", "2_pentylfuran", "4_vinylphenol",
                          "butyric_acid", "4_ethylphenol", "dimethyl_disulfide"],
            "expect": "NO PREDICTION AT ALL -> scored as sign-fail",
            "confidence": "certain (structural)",
            "why": (
                "Six of the ten compounds belong to classes for which the entire "
                "FIT column contains NO per-gram binding constant: alkylpyrazine, "
                "alkylfuran, phenol (x2), carboxylic acid, disulfide. The layer "
                "reports `no_measured_constant_for_this_class` and emits 1.0, "
                "which is the ABSENCE of a prediction, not a prediction of no "
                "shift. Imputing a constant for them would be inventing a "
                "parameter, and a class-mean over unrelated chemistry would be "
                "the general correction factor k2 sec. D.1 refutes. This is the "
                "single most informative thing the hold-out will show: the "
                "corpus supports a matrix prediction for 4 of Hong's 10 "
                "compounds and no more."),
        },
        {
            "compounds": ["hexanal", "3_methylbutanal", "2_methylbutanal"],
            "expect": "SIGN RIGHT, MAGNITUDE FAR TOO SMALL",
            "confidence": "high (~85%)",
            "why": (
                "These three do get a term. But Amendment 6 ruling 2 already "
                "computed the answer: reversible binding explains ~25 % of a "
                "hexanal matrix log-shift and covalent ~0.06 %. A model whose "
                "only active term is reversible binding must therefore "
                "UNDER-PREDICT a real matrix shift by roughly the remaining "
                "0.75 of the log, which for a shift of order 100x is a fold "
                "error of order 20-30x -- far outside the 5x window. If any of "
                "these three lands inside 5x it will be because the true shift "
                "is small, not because the model captured a large one."),
        },
        {
            "compounds": ["2_methylbutanal", "3_methylbutanal"],
            "expect": "IDENTICAL PREDICTIONS FOR TWO DIFFERENT COMPOUNDS",
            "confidence": "certain (structural)",
            "why": (
                "Both are C5 branched alkanals reached by the same chain-length "
                "surrogate from the C6 n-alkanal. Branch position is unmeasured "
                "anywhere in the corpus, so the layer CANNOT distinguish them. "
                "Any difference Hong measured between the two is, for this "
                "layer, pure unexplained residual by construction."),
        },
        {
            "compound": "dimethyl_disulfide",
            "expect": "carries an unresolved SOURCE CONTRADICTION",
            "confidence": "certain",
            "why": (
                "anantharamkrishnan2020b's Table 2 and its own Results text "
                "disagree on whether DMDS forms a covalent adduct. Resolved "
                "conservatively in favour of the table (no covalent term) and "
                "reported on the row rather than decided silently. It changes "
                "nothing numerically here, because the covalent term is inert "
                "anyway."),
        },
    ],
    "what_would_falsify_the_expectation": (
        "If the ester's measured ratio turns out to be ABOVE 1 the sign call "
        "would pass and this pre-registration would be wrong in the model's "
        "favour. If several of the six no-constant compounds happen to sit "
        "within 5x of 1.0, the magnitude count could reach 7/10 while the sign "
        "clause still fails -- that would be a pass on one half of a criterion "
        "that requires both, and it must not be reported as a partial pass."),
    "declared_before": "any read of results/validation/holdout_frozen/"
                       "hong2020_paired_thresholds.json by this wave.",
}

# =========================================================================
# DISCLOSURE -- hold-out exposure incurred by the instructed reading
# =========================================================================
HOLDOUT_EXPOSURE_DISCLOSURE = {
    "policy": (
        "Amendment 2 set the precedent: exposure that occurred during directed "
        "reads is disclosed in the report, and its non-use is enforced by a "
        "literal-grep firewall test rather than asserted."),
    "sealed_successfully": [
        "results/validation/holdout_frozen/hong2020_paired_thresholds.json was "
        "produced by a custodian agent under an explicit no-outcome reporting "
        "rule and was NOT read by this wave while the layer was built.",
        "data/lit/extraction_dossiers/hong2020_extraction.md sec. 4 (the paired "
        "BET table) was NOT opened; only the firewalled public manifest was used, "
        "and it carries structure, method and matrix description only.",
        "data/benchmarks/external_validation/ was never opened.",
    ],
    "exposure_incurred": [
        {
            "what": "Hong 2020 AGGREGATE statistics and four per-compound facts",
            "where": ("data/lit/extraction_dossiers/"
                      "k4b_paired_thresholds_and_browning.md, which the build "
                      "brief instructed this wave to read"),
            "detail": (
                "k4b's summary sections print, outside the paired table: the "
                "range and geometric mean of the soy shifts; the log10 SD; the "
                "fact that one of ten inverts and that it is the ester; the "
                "hexanal soy/water value; and the largest and smallest elevated "
                "values with the chemical class of each. That is roughly four of "
                "the ten rows disclosed in whole or in part."),
            "consequence": (
                "The Hong hold-out is NOT fully blind for this wave and must not "
                "be described as such. It is scored TWO WAYS in the hold-out "
                "report: over all ten rows, and over the SIX rows for which "
                "nothing leaked. The six-row score is the one that carries "
                "evidential weight."),
            "mitigation": (
                "None of the exposed values enters a parameter, a bound, an "
                "initialisation, or a class assignment. Enforced by "
                "tests/unit/test_kinetic_core_b4.py::test_holdout_firewall, "
                "which greps the executable code of every B4 runtime and fit file "
                "for the exposed literals."),
        },
        {
            "what": "Brewer 1995's six beef thresholds and their ratios",
            "where": "k2_matrix_and_thresholds.md sec. A.1 and sec. A.8, also "
                     "an instructed read",
            "detail": "Brewer is D.6 Module 7 HOLD-OUT (reclassified "
                      "`dose_added_pre_cook`). Its numbers were seen.",
            "consequence": (
                "The beef alkenal/alkanal ratio (2.01x) is the observation that "
                "would pull the unsaturation penalty DOWN toward k2's stated "
                "2-3x band. It is EXCLUDED from the fit, which leaves the "
                "fitted penalty above that band. Excluding it is the "
                "conservative direction for the fit and the honest one for the "
                "hold-out."),
        },
    ],
}


def _round(value, digits=6):
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return str(value)
        return round(value, digits)
    return value


def build() -> dict:
    classes = fit_class_binding_constants()
    penalty = fit_unsaturation_penalty()

    # ---- FIT-row self-consistency: reproduce what was fitted -----------
    fit_rows = []
    for compound, matrix, measured in (
        ("hexanal", "skim_milk", 1.39),
        ("ethyl_pentanoate", "skim_milk", 1.33),
        ("amyl_acetate", "skim_milk", 1.20),
        ("isoamyl_acetate", "skim_milk", 1.07),
        ("diacetyl", "caseinate_1pct", 1.73),
        ("delta_decalactone", "caseinate_1pct", 0.77),
        ("furaneol", "caseinate_1pct", 0.38),
    ):
        prediction = predict_matrix_shift(compound, matrix,
                                          class_constants=classes,
                                          unsaturation=penalty)
        fit_rows.append({
            "compound": compound, "matrix": matrix,
            "measured_ratio": measured,
            "predicted_ratio": _round(prediction.predicted_ratio),
            "fold_error": _round(max(measured, prediction.predicted_ratio)
                                 / min(measured, prediction.predicted_ratio)),
            "sign_measured": "elevated" if measured > 1 else "inverted",
            "sign_predicted": prediction.as_dict()["predicted_sign"],
            "state": prediction.state,
            "role": "FIT",
        })

    # ---- The Vega gelatin ladder, reproduced as lookup entries ---------
    gelatin_rows = []
    for compound in ("pentanal", "hexanal", "heptanal", "t_2_hexenal",
                     "t_2_octenal", "tt_2_4_decadienal"):
        water = select_threshold(compound, "water")
        gel = select_threshold(compound, "gelatin_3pct", temperature_c=22.0)
        measured_ratio = gel.value_ug_per_l / water.value_ug_per_l
        prediction = predict_matrix_shift(compound, "gelatin_3pct",
                                          class_constants=classes,
                                          unsaturation=penalty)
        gelatin_rows.append({
            "compound": compound,
            "threshold_water_ug_per_L": water.value_ug_per_l,
            "threshold_gelatin_22C_ug_per_L": gel.value_ug_per_l,
            "measured_ratio": _round(measured_ratio),
            "predicted_ratio": _round(prediction.predicted_ratio),
            "state": prediction.state,
            "role": "FIT (lookup-table entries, D.6 Module 7)",
            "cross_method_ratio": True,
        })

    # ---- THE PRE-REGISTERED BLIND PREDICTIONS ---------------------------
    frozen = []
    for compound in HONG_PANEL:
        prediction = predict_matrix_shift(compound, "soy_paste_hong",
                                          ph=None,
                                          class_constants=classes,
                                          unsaturation=penalty)
        row = prediction.as_dict()
        row["compound_display"] = COMPOUND_STRUCTURE[compound].display
        row["binding_class"] = COMPOUND_STRUCTURE[compound].binding_class
        row["covalent_gate"] = covalent_channel_state(compound)["state"]
        for key in ("predicted_matrix_over_water_ratio",):
            row[key] = _round(row[key])
        row["predicted_interval"] = [_round(v) for v in row["predicted_interval"]]
        frozen.append(row)

    # ---- The OAV code path, exercised on the sulfur panel ---------------
    # Zhou SI Table S2 prints MFT OAV 3.18e5 and disulfide 3.21e5 at pH 7; the
    # DIMER'S OAV MATCHES THE MONOMER'S despite carrying <10 % of the mass.
    demo_conc = {"MFT": 1.59, "MFTD": 0.1027, "FFT": 0.5, "ACTZ": 20.0,
                 "MMFT": 0.4, "H2S": 1.0}
    oav_demo = oav_table(demo_conc, matrix="water")

    # ---- The ratio path, exercised on two formulations ------------------
    ratio_demo = compare_formulations(
        {"MFT": 1.59, "MFTD": 0.1027, "hexanal": 400.0},
        {"MFT": 0.21, "MFTD": 0.0090, "hexanal": 1600.0},
        "pea_high_cysteine", "pea_control")

    interval_demo = absolute_concentration(
        400.0, via_partition=True, provenance="illustrative hexanal ppb").as_dict()

    return {
        "wave": "B4 -- matrix / OAV output layer",
        "generated_by": "scripts/generators/generate_kinetic_core_b4_fit.py",
        "generated_on": str(date.today()),
        "declaration": {
            "document": "docs/reference/FIT_HOLDOUT_DECLARATION.md, amendments 1-6",
            "fit_rows_used": [
                "Meynier 2002 partition RATIOS + enthalpies (Amendment 4) -- "
                "NEVER its absolute K_aw",
                "Leksrisompong 2010 K RATIOS (Amendment 4)",
                "Vega 1994 gelatin ladder as lookup-table entries (D.6 Module 7)",
                "Damodaran 1981 soy K's, Andriot 2000 beta-lg K's (D.6 Module 6)",
                "Zhou 2023 SI Table S2 + Guadagni-via-Vega thresholds as "
                "REFERENCE INPUTS for OAV (thresholds are inputs, not fit targets)",
            ],
            "fit_rows_declared_but_not_used_by_this_layer": [
                "Pereyra Gonzales k set (browning kinetics; no matrix/OAV role)",
                "Miao Ea (browning; no matrix/OAV role)",
                "Lievonen pH-drift (relevant to the pH ADDUCT GATE as context -- "
                "browning drives pH down 1.26-4.43 units in unbuffered systems, "
                "which would push a real system TOWARD the pH-3 gate -- but no "
                "value is used, because this layer has no pH trajectory state)",
            ],
            "holdouts_respected": [
                "Hong 2020 paired thresholds (GATING)",
                "Brewer 1995 beef (reclassified dose_added_pre_cook)",
                "Zhou OAV arithmetic check",
                "Barallat-Perez lupin + mucin binding constants",
                "Leksrisompong 24 BETs",
                "the frozen matrix-path bundles under "
                "data/benchmarks/external_validation/ -- NEVER OPENED",
            ],
        },
        "layer": layer_metadata(),
        "registry": matrix_registry_metadata(),
        "fitted": {
            "class_binding_constants": classes,
            "unsaturation_penalty": penalty,
            "note": ("Both are deterministic functions of the frozen registry. "
                     "There is no optimiser in this file and no free parameter: "
                     "'freezing the fit' and 'freezing the registry' are the "
                     "same act, which is why the fit is exactly reproducible."),
        },
        "fit_rows": fit_rows,
        "gelatin_ladder_rows": gelatin_rows,
        "pre_registered_holdout_expectations": PRE_REGISTERED,
        "frozen_blind_predictions_hong2020": frozen,
        "holdout_exposure_disclosure": HOLDOUT_EXPOSURE_DISCLOSURE,
        "source_contradictions": dict(SOURCE_CONTRADICTIONS),
        "unimplemented_candidate_terms": dict(UNIMPLEMENTED_CANDIDATE_TERMS),
        "sealed_or_refused_matrices": dict(SEALED_OR_REFUSED_MATRICES),
        "oav_code_path_demo": oav_demo,
        "ratio_code_path_demo": ratio_demo,
        "absolute_interval_demo": interval_demo,
        "counts": {
            "water_thresholds": len(WATER_THRESHOLDS),
            "matrix_thresholds": len(MATRIX_THRESHOLDS),
            "binding_constants": len(matrix_registry_metadata()["binding_constants"]),
            "hong_compounds_with_a_term": sum(
                1 for r in frozen if r["state"] != "no_binding_constant_for_class"),
            "hong_compounds_without_a_term": sum(
                1 for r in frozen if r["state"] == "no_binding_constant_for_class"),
        },
    }


def to_markdown(report: dict) -> str:
    lines = []
    add = lines.append
    add("# Kinetic core B4 -- matrix / OAV output layer: FIT REPORT")
    add("")
    add(f"Generated {report['generated_on']} by `{report['generated_by']}`.")
    add("")
    add("## The four-line answer")
    add("")
    add("1. **Ratios lead.** The layer's primary output is formulation-vs-formulation")
    add("   ratios with a validity bound; absolutes are emitted only as intervals")
    add("   carrying the measured reliability band (HS-SPME same-sample dispersion")
    add("   10-23x, air/water K +/-0.5 decades).")
    add("2. **Three named terms, and nothing beyond them.** Reversible binding")
    add("   (capped at ~25 % of a log-shift), the alpha,beta-unsaturation penalty,")
    add("   and the covalent ceiling as an INERT bound. Everything else is reported")
    add("   as quantified unexplained residual.")
    add(f"3. **The corpus supports a matrix prediction for "
        f"{report['counts']['hong_compounds_with_a_term']} of the flagship")
    add(f"   hold-out's 10 compounds and no more.** "
        f"{report['counts']['hong_compounds_without_a_term']} belong to classes for")
    add("   which the entire FIT column contains no binding constant.")
    add("4. **The gate is expected to fail, by construction.** The declared criterion")
    add("   requires the correct sign on all ten including an inversion; nothing in")
    add("   the FIT column produces a sign reversal for an ester.")
    add("")
    add("## Fitted constants (FIT rows only)")
    add("")
    add("| class | K_g (L/g) | n FIT rows | members | sign | sources |")
    add("|---|---:|---:|---|---|---|")
    for name, entry in sorted(report["fitted"]["class_binding_constants"].items()):
        members = ", ".join(entry["members"]) or "(chain-length surrogate)"
        add(f"| {name} | {entry['k_g_l_per_g']:.4g} | {entry['n_fit_rows']} | "
            f"{members} | {entry['sign']} | {', '.join(entry['sources'])} |")
    add("")
    penalty = report["fitted"]["unsaturation_penalty"]
    add(f"**alpha,beta-unsaturation penalty: {penalty['penalty_x']:.3g}x** from "
        f"{penalty['n_fit_rows']} FIT rows "
        f"({', '.join(f'{k} {v}x' for k, v in penalty['observations'].items())}), "
        f"spread {penalty['spread_x']:.2f}x.")
    add("")
    add(f"> {penalty['caveat']}")
    add("")
    add("## FIT-row reproduction")
    add("")
    add("| compound | matrix | measured | predicted | fold | sign ok | state |")
    add("|---|---|---:|---:|---:|---|---|")
    for row in report["fit_rows"]:
        ok = "yes" if row["sign_measured"] == row["sign_predicted"] else "NO"
        add(f"| {row['compound']} | {row['matrix']} | {row['measured_ratio']} | "
            f"{row['predicted_ratio']} | {row['fold_error']} | {ok} | {row['state']} |")
    add("")
    add("## The Vega gelatin ladder as lookup-table entries")
    add("")
    add("| compound | water (ug/L) | gelatin 22 C (ug/L) | measured ratio | "
        "predicted | state |")
    add("|---|---:|---:|---:|---:|---|")
    for row in report["gelatin_ladder_rows"]:
        add(f"| {row['compound']} | {row['threshold_water_ug_per_L']} | "
            f"{row['threshold_gelatin_22C_ug_per_L']} | {row['measured_ratio']} | "
            f"{row['predicted_ratio']} | {row['state']} |")
    add("")
    add("## PRE-REGISTERED blind predictions -- Hong 2020 flagship")
    add("")
    add("| # | compound | class | predicted ratio | interval | predicted sign | state |")
    add("|---:|---|---|---:|---|---|---|")
    for i, row in enumerate(report["frozen_blind_predictions_hong2020"], 1):
        lo, hi = row["predicted_interval"]
        add(f"| {i} | {row['compound_display']} | {row['binding_class']} | "
            f"{row['predicted_matrix_over_water_ratio']} | [{lo}, {hi}] | "
            f"{row['predicted_sign']} | {row['state']} |")
    add("")
    add("## Pre-registered expected failures")
    add("")
    add(f"**Criterion:** {PRE_REGISTERED['gating_criterion']}")
    add("")
    add(f"**Overall expectation: {PRE_REGISTERED['overall_expectation']}**")
    add("")
    for item in PRE_REGISTERED["expected_failures"]:
        who = item.get("compound") or ", ".join(item.get("compounds", []))
        add(f"### {who} -- {item['expect']} ({item['confidence']})")
        add("")
        add(item["why"])
        add("")
    add("### What would falsify this expectation")
    add("")
    add(PRE_REGISTERED["what_would_falsify_the_expectation"])
    add("")
    add("## Hold-out exposure disclosure")
    add("")
    add(HOLDOUT_EXPOSURE_DISCLOSURE["policy"])
    add("")
    add("**Sealed successfully:**")
    add("")
    for item in HOLDOUT_EXPOSURE_DISCLOSURE["sealed_successfully"]:
        add(f"- {item}")
    add("")
    add("**Exposure incurred:**")
    add("")
    for item in HOLDOUT_EXPOSURE_DISCLOSURE["exposure_incurred"]:
        add(f"- **{item['what']}** ({item['where']}). {item['detail']}")
        add(f"  - *Consequence:* {item['consequence']}")
        add(f"  - *Mitigation:* {item.get('mitigation', 'see consequence')}")
    add("")
    add("## Named terms NOT implemented, and why")
    add("")
    for name, why in report["unimplemented_candidate_terms"].items():
        add(f"- **{name}** -- {why}")
    add("")
    add("## Matrices the layer refuses to tabulate")
    add("")
    for name, why in report["sealed_or_refused_matrices"].items():
        add(f"- **{name}** -- {why}")
    add("")
    add("## Source contradictions carried, not resolved")
    add("")
    for name, why in report["source_contradictions"].items():
        add(f"- **{name}** -- {why}")
    add("")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fit Build Wave B4's matrix / OAV output layer on fit rows only and "
            "freeze the blind Hong 2020 predictions; writes "
            "results/validation/kinetic_core_b4_fit_report.{json,md} and "
            "kinetic_core_b4_frozen_predictions.json."
        )
    )
    parser.parse_args(argv)

    report = build()
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2, default=str))
    OUT_MD.write_text(to_markdown(report))
    FROZEN_PREDICTIONS.write_text(json.dumps({
        "frozen_on": str(date.today()),
        "declaration": "Amendment 4 gating hold-out: Hong 2020 paired ratios.",
        "rule": ("These predictions were written BEFORE any read of "
                 "hong2020_paired_thresholds.json by this wave. The hold-out "
                 "scorer reads this file and changes nothing in it."),
        "pre_registered_expectations": PRE_REGISTERED,
        "predictions": report["frozen_blind_predictions_hong2020"],
        "class_binding_constants": report["fitted"]["class_binding_constants"],
        "unsaturation_penalty": report["fitted"]["unsaturation_penalty"],
    }, indent=2, default=str))
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    print(f"wrote {OUT_MD.relative_to(REPO)}")
    print(f"wrote {FROZEN_PREDICTIONS.relative_to(REPO)}")
    print(f"  {report['counts']['hong_compounds_with_a_term']} of 10 Hong compounds "
          f"have a term; {report['counts']['hong_compounds_without_a_term']} do not.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
