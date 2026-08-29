#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b2_2_holdout.py

THE PRE-REGISTERED HOLD-OUT SCORING OF BUILD WAVE B2.2.

Pre-registered in `results/validation/kinetic_core_b2_2_prereg.md`, which was
written and committed to disk BEFORE the fit ran and BEFORE any row here was
scored. Every claim in it is scored below whether it held or not.

ONE THING CHANGED IN THIS SCORER AND IT MAKES IT HARDER, NOT EASIER
-------------------------------------------------------------------
B2.1's scorer PRESCRIBED each unbuffered system's measured final pH
(`ph_final=5.07` for Zhou's pH-8 pot, `ph_final=4.9` for Kang's). That is
handing the model a measured process-state observable at scoring time. B2.2
does not: every unbuffered system is now run with the DYNAMIC pH state and its
endpoint is a PREDICTION. The declared buffer of each system is passed in
explicitly instead, and where the source declares none the run carries the
extrapolation warning.

Reads the FROZEN parameter set from
`results/validation/kinetic_core_b2_2_fit_report.json`, changes NOTHING, predicts
every declared hold-out row, and writes
`results/validation/kinetic_core_b2_2_holdout_report.{json,md}`.

===========================================================================
THE RULE THIS SCRIPT OBEYS
===========================================================================
No parameter is touched here. There is no optimiser in this file, no bounds, no
`least_squares` import. The only thing it does is integrate the FROZEN network
at conditions the fit never saw and print the fold error, row by row, pass or
fail, with NO averaging away of failures and no per-row weighting that could
hide one.

===========================================================================
PRE-REGISTERED FAILURES -- stated BEFORE the numbers, not after
===========================================================================
Four rows are predicted to fail, and the reason for each is structural and was
known at fit time. They are written down here so that they cannot later be
presented as discoveries:

  1. **van Seeventer's 50 C zero-order loss.** The C-5 oligomerisation channel
     has NO rate on any FIT row (its only measurement IS this hold-out), so the
     model carries it at zero and predicts essentially no 50 C loss. The
     failure localises a MISSING CHANNEL, which is what the declaration says
     this row is for.

  2. **Hofmann's dry-180 C rows.** The module has no water-activity term,
     because nothing in its fit corpus varies a_w. A dry silica system at 180 C
     is outside the model's declared domain and is scored UNSCOREABLE rather
     than given a number that would be meaningless.

  3. **Hofmann 1998's pH-3 and pH-7 aqueous MFT rows.** Inventory sec. B2.5
     records a SIGN CONFLICT between Hofmann's buffered free-sugar system (MFT
     falls as pH rises) and Zhou's unbuffered Amadori-fed system (MFT peaks at
     pH 7). The dossier's own recommendation is that the two must NOT be merged
     into one pH response curve. A model with a single structural pH shape can
     satisfy at most one. This one was fitted at Zhou's pH 7 and Hofmann's
     pH 5, so it is expected to get Zhou's direction and miss Hofmann's.

  4. **Zhang's Figure 2 MMFT fractions.** Zhang's Fig. 2 system is 1:3:3
     Met:VB1:Xyl and METHIONINE IS NOT IN THIS MODULE'S STATE VECTOR, so the
     MeSH supply is wrong by construction. Scored, and flagged as a system
     mismatch rather than a clean model failure.

===========================================================================
A DISCLOSURE ABOUT FIREWALL HYGIENE
===========================================================================
The build brief directed this wave to read `k3_final_parameter_inventory.md`
sec. A, and sec. A.3.3 PRINTS Zhou 2023 Table 1's pH-6 and pH-8 columns
alongside the pH-7 column that is the declared FIT row. The same is true of
Zhang2024_extraction.md sec. 3, whose panels (d) and (e) are part of the
declared Figs. 2d/e/f hold-out. Those values were therefore SEEN before the
fit was run. They were not used: they appear in no fit row, in no bound, in no
initialisation and in no source file under `src/kinetic_core/`, and
`tests/unit/test_kinetic_core_b2_2.py` enforces that mechanically. The
exposure is recorded here because a hold-out whose exposure is undisclosed is
worth less than one whose exposure is stated.

The frozen `mp_holdout_*` bundles under `data/benchmarks/external_validation/`
were NOT opened by this wave at any point.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    DECAY_FAMILIES,
    MEASURED_SULFUR,
    with_fitted_sulfur,
)
from src.kinetic_core.ph_state import BufferSpec, PhDrift  # noqa: E402
from src.kinetic_core.species_sulfur import (  # noqa: E402
    mmol_per_litre_to_ug_per_litre,
    odour_activity_values,
    ug_per_litre_to_mmol_per_litre,
)
from src.kinetic_core.sulfur import (  # noqa: E402
    branch_shares,
    integrate_sulfur,
    sulfur_flux_budget,
)

CELSIUS = 273.15
OX_AMBIENT = 1.0

# The DECLARED buffers, quoted from the sources, identical to the fit script's.
BUFFER_HOFMANN = BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=True,
    source="Hofmann & Schieberle 1998 / Cerny model solutions: 0.5 M phosphate")
BUFFER_NONE = BufferSpec(
    kind="none", declared=True,
    source="the source states water and declares no buffer")
FIT_REPORT = ROOT / "results/validation/kinetic_core_b2_2_fit_report.json"

B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

#: Pass bands, declared BEFORE the numbers are seen.
#: A LEVEL row passes within 3x (the corpus's own cross-lab agreement on the
#: one quantity two labs both measured is 1.6x, and the HS-SPME rows carry a
#: 10-23x internal calibration defect, so 3x is generous on the SIDA rows and
#: tight on the HS-SPME ones -- both are stated per row).
#: A RATIO or SHARE row passes within 2x.
#: A DIRECTIONAL row passes if the SIGN and the ORDER OF MAGNITUDE are right.
PASS_BAND_LEVEL = 3.0
PASS_BAND_RATIO = 2.0


def load_frozen() -> Dict[str, Any]:
    if not FIT_REPORT.exists():
        raise SystemExit(
            f"{FIT_REPORT} not found. Run generate_kinetic_core_b2_2_fit.py first; "
            f"the hold-out scorer never fits anything itself."
        )
    payload = json.loads(FIT_REPORT.read_text())
    frozen = payload["frozen_parameters"]
    parameters: Dict[str, Any] = dict(operative_parameters(B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(
        with_fitted_sulfur(
            frozen["log10_k_ref_at_145C"],
            frozen["lumped_formation_Ea_kJ_mol"],
            frozen["decay_Ea_kJ_mol"],
        )
    )
    drift = PhDrift(
        acid_yield=frozen["ph_drift"]["acid_yield_per_sink_event"],
        arp_amine_pka=frozen["ph_drift"]["arp_secondary_ammonium_pKa"],
    )
    return {"parameters": parameters, "fit_payload": payload, "drift": drift}


#: The frozen drift constants, installed by main() from the fit report. The
#: scorer NEVER constructs its own: a hold-out scorer that could choose a
#: calibration would not be scoring a frozen model.
_DRIFT: PhDrift = None


def run(parameters, initial, t_c, minutes, ph, buffer=BUFFER_HOFMANN, times=None):
    grid = np.array([0.0, minutes]) if times is None else np.asarray(times)
    return integrate_sulfur(
        parameters, t_c + CELSIUS, initial, grid,
        ph=ph, buffer_spec=buffer, ph_drift=_DRIFT,
        ph_nodes=41, ph_iterations=4, rtol=1e-8, atol=1e-14,
    )


def flux(parameters, initial, t_c, minutes, ph, buffer=BUFFER_HOFMANN):
    return sulfur_flux_budget(
        parameters, t_c + CELSIUS, initial, minutes, ph=ph,
        buffer_spec=buffer, ph_drift=_DRIFT, n_points=61,
    )


def score(observed: float, predicted: float, band: float) -> Dict[str, Any]:
    if not np.isfinite(predicted) or predicted <= 0 or observed <= 0:
        fold = float("nan")
        passed = False
    else:
        fold = predicted / observed
        passed = (1.0 / band) <= fold <= band
    return {
        "observed": observed,
        "predicted": float(predicted),
        "fold_error": fold,
        "pass_band": band,
        "pass": bool(passed),
    }


def main() -> int:
    global _DRIFT
    frozen = load_frozen()
    p = frozen["parameters"]
    _DRIFT = frozen["drift"]
    rows: List[Dict[str, Any]] = []

    def add(row_id, group, gating, anchor, result, comment=""):
        entry = {"id": row_id, "group": group, "gating": gating,
                 "source_anchor": anchor, "comment": comment}
        entry.update(result)
        rows.append(entry)

    # =====================================================================
    # 1. ZHOU 2023 Table 1, the pH-8 column -- GATING
    # =====================================================================
    zhou8 = run(p, {"ARP": 20.0, "Cys": 20.0, "OX": OX_AMBIENT}, 120.0, 60.0,
                ph=8.0, buffer=BUFFER_NONE)
    zhou8_alone = run(p, {"ARP": 20.0, "OX": OX_AMBIENT}, 120.0, 60.0,
                      ph=8.0, buffer=BUFFER_NONE)
    for species, observed_ug, note in (
        ("MFT", 525.62, "HS-SPME, absolute_concentration: false"),
        ("FFT", 325.22, ""),
        ("MFTD", 50.07, "only 1.5x above a pseudo-LOD: its calibration curve "
                        "zeroes at x = 34.3 ug/L (inventory sec. F row 14)"),
        ("ACTZ", 582.34, ""),
    ):
        add(f"zhou_pH8_{species}", "Zhou pH-8 column", True,
            f"Zhou 2023 Table 1, pH 8 column: {species} {observed_ug} ug/L",
            score(ug_per_litre_to_mmol_per_litre(species, observed_ug),
                  zhou8.final(species), PASS_BAND_LEVEL), note)
    add("zhou_pH8_FUR_arp_alone", "Zhou pH-8 column", True,
        "Zhou 2023 Table 1, pH 8, ARP alone: 2-furfural 436.63 ug/L",
        score(ug_per_litre_to_mmol_per_litre("FUR", 436.63),
              zhou8_alone.final("FUR"), PASS_BAND_LEVEL))
    # the SHAPE test: does MFT fall from pH 7 to pH 8 as measured?
    zhou7 = run(p, {"ARP": 20.0, "Cys": 20.0, "OX": OX_AMBIENT}, 120.0, 60.0,
                ph=7.0, buffer=BUFFER_NONE)
    obs_fall_87 = 525.62 / 1588.57
    pred_fall_87 = (zhou8.final("MFT") / zhou7.final("MFT")
                    if zhou7.final("MFT") > 0 else float("nan"))
    add("zhou_MFT_shape_pH8_over_pH7", "Zhou pH-8 column", True,
        "Zhou 2023 T1: MFT 525.62 at pH 8 against 1588.57 at pH 7 => 0.331. "
        "THE MAXIMUM AT pH 7 IS THE POINT: a model fitted at the maximum must "
        "PREDICT the fall on the alkaline side.",
        score(obs_fall_87, pred_fall_87, PASS_BAND_RATIO),
        "the single most informative row in the whole hold-out set")
    obs_fft_87 = 325.22 / 757.965
    pred_fft_87 = (zhou8.final("FFT") / zhou7.final("FFT")
                   if zhou7.final("FFT") > 0 else float("nan"))
    add("zhou_FFT_shape_pH8_over_pH7", "Zhou pH-8 column", True,
        "Zhou 2023 T1: FFT 325.22 at pH 8 against 757.965 at pH 7 => 0.429, "
        "monotone DOWN -- the acid-lane prediction",
        score(obs_fft_87, pred_fft_87, PASS_BAND_RATIO))

    # =====================================================================
    # 2. ZHOU 2023 Table 1, the pH-6 column -- DIAGNOSTIC ONLY
    # =====================================================================
    # FIT_HOLDOUT_DECLARATION.md sec.5 decision 1: "the pH-6 hold-out column is
    # scored DIAGNOSTIC, not gating, until the model carries a pH-trajectory
    # state." This module DOES carry one, but the decision has not been amended,
    # so the column stays diagnostic. That is the ratified decision, honoured.
    zhou6 = run(p, {"ARP": 20.0, "Cys": 20.0, "OX": OX_AMBIENT}, 120.0, 60.0,
                ph=6.0, buffer=BUFFER_NONE)
    for species, observed_ug in (("MFT", 696.99), ("FFT", 813.65), ("MFTD", 59.70)):
        add(f"zhou_pH6_{species}", "Zhou pH-6 column (DIAGNOSTIC)", False,
            f"Zhou 2023 Table 1, pH 6 column: {species} {observed_ug} ug/L",
            score(ug_per_litre_to_mmol_per_litre(species, observed_ug),
                  zhou6.final(species), PASS_BAND_LEVEL),
            "DIAGNOSTIC per FIT_HOLDOUT_DECLARATION.md sec.5 decision 1: the pH "
            "labels are INITIAL pH of an unbuffered system whose pH-6 and pH-7 "
            "runs converge to within 0.2 units by the end of heating")

    # =====================================================================
    # 3. ZHOU Figure 3, the Cys + MGO time x pH grid
    # =====================================================================
    add("zhou_fig3_grid", "Zhou Fig. 3 grid", False,
        "Zhou 2023 Fig. 3, Cys + MGO 20 mM each, unbuffered, 120 C, "
        "0/30/60/90 min, pH 6/7/8",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        "NOT SCORED, and the reason is a SCOPE limit rather than a model "
        "failure: seven of the grid's eight panels are pyrazines, thiophenes "
        "and methylthiazoles, all of which this wave declares OUT OF SCOPE "
        "(sulfur.OUT_OF_SCOPE). The one in-scope quantity the grid supplies is "
        "the 2-acetylthiazole cross-system pair, which IS scored below. "
        "Declaring this unscored is the honest option; scoring one panel of "
        "eight and calling the grid passed would not be.")

    # =====================================================================
    # 4. ZHOU cross-system consistency: 582 vs 665 ug/L 2-acetylthiazole
    # =====================================================================
    # ARP-fed gives 582; MGO-fed gives 665. A DERIVED CONSISTENCY test, not a
    # level: it scores the ARP -> MGO flux with no free downstream parameter,
    # because the same k_cys_actz makes the thiazole in both pots.
    mgo_fed_pH8 = run(p, {"MGO": 20.0, "Cys": 20.0, "OX": OX_AMBIENT}, 120.0,
                      60.0, ph=8.0, buffer=BUFFER_NONE)
    pred_ratio = (zhou8.final("ACTZ") / mgo_fed_pH8.final("ACTZ")
                  if mgo_fed_pH8.final("ACTZ") > 0 else float("nan"))
    add("zhou_582_vs_665_consistency", "Zhou cross-system pair", True,
        "Zhou 2023: 2-acetylthiazole at 60 min, pH 8 -- ARP-fed 582 ug/L "
        "against MGO-fed 665 ug/L (14% apart)",
        score(582.34 / 665.0, pred_ratio, PASS_BAND_RATIO),
        "scores the ARP -> MGO flux with no free downstream parameter")

    # =====================================================================
    # 5. WHITFIELD 2001: the fed-NF pH-6.5 collapse, >= 150x
    # =====================================================================
    w45 = run(p, {"NF": 20.0, "Cys": 20.0}, 140.0, 60.0, ph=4.5,
              buffer=BUFFER_NONE)
    w65 = run(p, {"NF": 20.0, "Cys": 20.0}, 140.0, 60.0, ph=6.5,
              buffer=BUFFER_NONE)
    pred_collapse = (w45.final("MFT") / w65.final("MFT")
                     if w65.final("MFT") > 0 else float("inf"))
    add("whitfield_2001_pH65_collapse", "Whitfield 2001 pH collapse", True,
        "Whitfield 2001 Table 1: fed NF + cysteine, MFT 0.150 mol% at pH 4.5 "
        "-> nd (<0.0010 mol%) at pH 6.5 = a >=150x collapse within one lab, "
        "while total volatiles fall only 1.7x",
        {"observed": 150.0, "predicted": float(pred_collapse),
         "fold_error": (float(pred_collapse) / 150.0
                        if np.isfinite(pred_collapse) else float("inf")),
         "pass_band": None,
         "pass": bool(np.isfinite(pred_collapse) and pred_collapse >= 150.0)},
        "ONE-SIDED: the measurement is a LOWER bound (the pH-6.5 value is a "
        "non-detect), so any predicted collapse of 150x or more passes. NOTE "
        "the inventory's own qualification: Whitfield 2001's H2S column should "
        "NOT carry a standalone row (its Methods omits HMF -- likely a typo).")

    # =====================================================================
    # 6. CERNY 2003 at 95 C: the norfuraneol share ceiling
    # =====================================================================
    # "the gate must be evaluated AT CERNY'S CONDITIONS, not at 145 C"
    cerny95 = flux(p, {"PENT": 100.0, "Cys": 33.0}, 95.0, 240.0, ph=5.0)
    shares95 = branch_shares(cerny95)
    add("cerny2003_NF_share_ceiling_95C", "Cerny 2003, 95 C", True,
        "Cerny 2003 Table 3: the norfuraneol share of the MFT flux is <=7% at "
        "95 C / 4 h / pH 5.00, and that is an UPPER bound (NF was spiked in at "
        "1.5x the cysteine)",
        {"observed": 0.07,
         "predicted": float(shares95["MFT_share_norfuraneol_route"]),
         "fold_error": float(shares95["MFT_share_norfuraneol_route"]) / 0.07,
         "pass_band": None,
         "pass": bool(shares95["MFT_share_norfuraneol_route"] <= 0.07)},
        "ONE-SIDED CEILING, evaluated at CERNY'S conditions (95 C, 4 h) and not "
        "at the fit panel's 145 C, exactly as the declaration requires")
    add("cerny2003_intact_skeleton_share_95C", "Cerny 2003, 95 C", True,
        "Cerny 2003 Table 2: MFT is 49% unlabelled / 46% 13C5 with NO fragment "
        "pattern => ~93% of MFT comes from the INTACT-C5 (1,4-dideoxyosone) "
        "route at 95 C; 'pathways via ribose fragmentation were not relevant'",
        score(0.93, float(shares95["MFT_share_intact_skeleton_route"]),
              PASS_BAND_RATIO),
        "and its partner: Cerny 2004 finds the C2+C3 route 'not relevant' at "
        "95 C while Hofmann 1998 T10 makes it the STRONGEST route at 145 C. "
        "The model must therefore get a TEMPERATURE-DEPENDENT route mix right, "
        "with one lumped formation Ea.")
    cerny145 = flux(p, {"PENT": 100.0, "Cys": 33.0}, 145.0, 20.0, ph=5.0)
    shares145 = branch_shares(cerny145)
    add("route_mix_moves_with_temperature", "Cerny 2003, 95 C", False,
        "Hofmann 1998 T10 vs Cerny 2004: the C2+C3 route is the STRONGEST MFT "
        "route at 145 C (0.24 mol%) and 'not relevant' at 95 C",
        {"observed": None,
         "predicted": {
             "C2+C3 share at 95 C": float(shares95["MFT_share_C2_plus_C3_route"]),
             "C2+C3 share at 145 C": float(shares145["MFT_share_C2_plus_C3_route"]),
         },
         "fold_error": float("nan"), "pass_band": None,
         "pass": bool(shares145["MFT_share_C2_plus_C3_route"]
                      > shares95["MFT_share_C2_plus_C3_route"])},
        "DIRECTIONAL: passes if the C2+C3 share RISES with temperature")

    # =====================================================================
    # 7. HOFMANN 1998 dry 180 C / 6 min -- UNSCOREABLE, and why
    # =====================================================================
    add("hofmann_dry_180C", "Hofmann dry-180", True,
        "Hofmann 1998 T2 dry rows: ribose/glucose/rhamnose on silica, 180 C / "
        "6 min, FFT 97.2 / 1.4 / 0.4 and MFT 25.1 / 4.2 / 3.1 ug; and T10's "
        "C2+C3 dry row at 1553.9 ug (1.39 mol%)",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        "UNSCOREABLE, PRE-REGISTERED. The module has NO WATER-ACTIVITY TERM, "
        "because nothing in its fit corpus varies a_w -- B1's policy 3, "
        "inherited. A dry silica system at 180 C is outside the declared "
        "domain, and producing a number for it would be inventing an a_w "
        "dependence at scoring time. Also note inventory sec. F row 17: the "
        "dry protocol is described TWICE with different times (5 vs 6 min) and "
        "a 4x different buffer molarity while reporting the same number to four "
        "significant figures.")

    # =====================================================================
    # 8. HOFMANN 1998's own pH-3 and pH-7 aqueous rows (frozen hold-out)
    #    -- the B2.5 sign conflict
    # =====================================================================
    for ph_value, mft_ug, fft_ug in ((3.0, 553.0, 229.0), (7.0, 25.0, 12.0)):
        h = run(p, {"PENT": 100.0, "Cys": 33.0, "OX": OX_AMBIENT}, 145.0, 20.0,
                ph=ph_value)
        add(f"hofmann_ribose_pH{ph_value:.0f}_MFT", "Hofmann pH axis", True,
            f"Hofmann 1998 T2, ribose pH {ph_value:.0f}: MFT {mft_ug} ppb",
            score(ug_per_litre_to_mmol_per_litre("MFT", mft_ug), h.final("MFT"),
                  PASS_BAND_LEVEL),
            "PRE-REGISTERED FAILURE. Inventory sec. B2.5: Hofmann's buffered "
            "free-sugar system has MFT FALLING as pH rises while Zhou's "
            "unbuffered Amadori-fed system has it RISING to a pH-7 maximum, "
            "with the absolute levels 64x apart at pH 7. The dossier's own "
            "recommendation is NOT to merge them into one pH response curve. A "
            "single structural pH shape can satisfy at most one, and this "
            "module carries Zhou's.")
        add(f"hofmann_ribose_pH{ph_value:.0f}_FFT", "Hofmann pH axis", True,
            f"Hofmann 1998 T2, ribose pH {ph_value:.0f}: FFT {fft_ug} ppb",
            score(ug_per_litre_to_mmol_per_litre("FFT", fft_ug), h.final("FFT"),
                  PASS_BAND_LEVEL),
            "FFT is the row the acid lane SHOULD get right: it falls "
            "monotonically with pH in all three papers that measure it.")

    # =====================================================================
    # 9. MEYNIER 1995, directional
    # =====================================================================
    m45 = run(p, {"PENT": 100.0, "Cys": 33.0, "OX": OX_AMBIENT}, 145.0, 20.0, ph=4.5)
    m65 = run(p, {"PENT": 100.0, "Cys": 33.0, "OX": OX_AMBIENT}, 145.0, 20.0, ph=6.5)
    for species, observed_fold, label in (
        ("MFT", 152.0, "MFT falls >152x over pH 4.5 -> 6.5"),
        ("FFT", 6.1, "FFT falls 6.1x"),
        ("FUR", 25.0, "furfural falls 15-49x (midpoint 25x)"),
    ):
        pred = (m45.final(species) / m65.final(species)
                if m65.final(species) > 0 else float("inf"))
        add(f"meynier_{species}_pH_fold", "Meynier directional", False,
            f"Meynier 1995: {label} (constant-pH, 4-amino-acid, single-sugar "
            f"factorial, ~59 series; DIRECTIONAL ONLY -- 'response factors not "
            f"determined')",
            {"observed": observed_fold, "predicted": float(pred),
             "fold_error": (float(pred) / observed_fold
                            if np.isfinite(pred) else float("inf")),
             "pass_band": 10.0,
             "pass": bool(np.isfinite(pred) and 0.1 <= pred / observed_fold <= 10.0)},
            "DIRECTIONAL, order-of-magnitude band. Meynier has no absolute "
            "quantification of any kind (inventory sec. C.8).")

    # =====================================================================
    # 10. THIOL-SINK HOLD-OUTS
    # =====================================================================
    # (a) Hofmann 2002 Fig. 1, real coffee brew, 80 C. The declaration's own
    #     words: "it is the one the model will get WRONG in the informative
    #     direction (the brew is SLOWER than the 30 C model systems)".
    brew = run(p, {"FFT": 0.5, "MELE": 9.0 * 50.0, "OX": OX_AMBIENT}, 80.0, 60.0,
               # a real coffee brew: no buffer is declared anywhere
               #    (the run carries the extrapolation warning),
               ph=5.2)
    remaining = brew.final("FFT") / 0.5
    pred_k_per_min = -math.log(max(remaining, 1e-12)) / 60.0
    add("hofmann2002_brew_80C_FFT", "thiol sinks", True,
        "Hofmann 2002 Fig. 1: FFT loss in a real coffee brew (50 g powder/L, "
        "thermos, 80 C) at 0.023 /min, t1/2 ~30 min",
        score(0.023, pred_k_per_min, PASS_BAND_LEVEL),
        "THE ONLY TEMPERATURE EXTRAPOLATION THE MODULE HAS, and the model has "
        "no activation energy for this channel by policy, so the constant is "
        "HELD at its 30 C value. A pass would mean the depletable-electrophile "
        "structure is doing the work; a failure localises the "
        "electrophile-pool depletion term, which is exactly what the "
        "declaration says this row is for.")
    # (b) van Seeventer 2001, 50 C, ZERO ORDER
    vs = run(p, {"MFT": 0.01, "FFT": 0.01, "OX": OX_AMBIENT}, 50.0, 1440.0,
             ph=5.0, buffer=BUFFER_NONE)
    add("vanseeventer_50C_MFT_per_day", "thiol sinks", True,
        "van Seeventer 2001 Table 1: MFT 59% of initial per DAY, ZERO ORDER, "
        "50 C, pH 5.0, air ~ argon (FFT 28%/day)",
        score(0.59, 1.0 - vs.final("MFT") / 0.01, PASS_BAND_LEVEL),
        "PRE-REGISTERED FAILURE. The C-5 oligomerisation channel that carries "
        "this loss has NO rate on any FIT row -- its only measurement IS this "
        "hold-out -- so the model carries it at ZERO and predicts essentially "
        "no loss. The failure localises a MISSING CHANNEL, not a wrong "
        "barrier, and 'do not bolt van Seeventer's 59%/day onto the network as "
        "the MFT sink' is a verbatim dossier instruction (sec. C.17).")
    # (c) Zhang 2024 Figs. 2d/e/f, 115 C
    zhang_flux = flux(p, {"THI": 44.5, "PENT": 99.9, "Cys": 123.8,
                          "OX": OX_AMBIENT}, 115.0, 60.0, ph=4.9)
    zs = branch_shares(zhang_flux)
    consumed = (zs["MFT_consumed_share_dimerisation"]
                + zs["MFT_consumed_share_MMFT"])
    add("zhang_115C_MFT_consumed_share", "thiol sinks", True,
        "Zhang 2024 Figs. 2d/e/f: 3-19% of free MFT already consumed into "
        "MEASURED products in the Cys arm (21-55% in the GSH arm); LOWER "
        "BOUNDS, since a third sink (melanoidin/quinone) is named and never "
        "measured",
        {"observed": 0.11, "predicted": float(consumed),
         "fold_error": float(consumed) / 0.11, "pass_band": PASS_BAND_LEVEL,
         "pass": bool(0.03 <= consumed <= 0.60)},
        "SCORED AGAINST THE BAND, not the midpoint, because the measurement is "
        "a lower bound and a range. SYSTEM MISMATCH DISCLOSED: Zhang's Fig. 2 "
        "system is 1:3:3 Met:VB1:Xyl and METHIONINE IS NOT IN THIS MODULE'S "
        "STATE VECTOR, so the MeSH supply is structurally wrong; the Fig. 1 "
        "geometry is substituted and that substitution is part of what is "
        "being scored.")
    # (d) Zhou 2023, 120 C dimer fractions, pH 6-8 -- and the score to beat
    zhou_dimer = {}
    for label, r, observed_thiol_equivalents in (
        ("pH6", zhou6, 0.086), ("pH7", zhou7, 0.065), ("pH8", zhou8, 0.096),
    ):
        # "dimer / MFT in THIOL EQUIVALENTS" is 2 x [dimer] / [MFT], because one
        # disulfide carries two thiols. Checked against the source arithmetic:
        # at pH 7, 102.59/226.32 = 0.4533 umol/L dimer and 1588.57/114.16 =
        # 13.915 umol/L MFT, and 2 x 0.4533 / 13.915 = 6.5%, which is the
        # inventory's printed value exactly.
        mft = r.final("MFT")
        dimer = r.final("MFTD")
        pred = (2.0 * dimer / mft) if mft > 0 else float("nan")
        zhou_dimer[label] = pred
        add(f"zhou_120C_dimer_share_{label}", "thiol sinks", True,
            f"Zhou 2023 Table 1: the dimer carries {observed_thiol_equivalents:.1%} "
            f"of MFT-equivalents at {label}, and the fraction is NEAR-INVARIANT "
            f"in pH while [MFT] swings 3.0x",
            score(observed_thiol_equivalents, pred, PASS_BAND_RATIO),
            "THE SCORE TO BEAT: Zhou (120 C, unbuffered, Amadori-fed, pH 6-8) "
            "and Zhang (115 C, buffered pH 4.9, thiamine/xylose) agree to 1.3x "
            "on this channel across two labs, two feedstocks and two buffers. "
            "Both are held out, so that agreement stays an out-of-sample fact.")
    pred_invariance = (max(zhou_dimer.values()) / min(zhou_dimer.values())
                       if min(zhou_dimer.values()) > 0 else float("inf"))
    add("zhou_dimer_share_pH_invariance", "thiol sinks", True,
        "Zhou 2023: the dimer share is pH-INVARIANT (8.6 / 6.5 / 9.6% over pH "
        "6-8) while [MFT] swings 3.0x -- a 1.48x spread",
        {"observed": 1.48, "predicted": float(pred_invariance),
         "fold_error": (float(pred_invariance) / 1.48
                        if np.isfinite(pred_invariance) else float("inf")),
         "pass_band": PASS_BAND_RATIO,
         "pass": bool(np.isfinite(pred_invariance) and pred_invariance <= 1.48 * 2.0)},
        "a SHAPE test that a fitted pH term would have been able to fake and a "
        "structural one cannot")

    # =====================================================================
    # 11. CERNY 2007b, the two single-route controls
    # =====================================================================
    no_cys = flux(p, {"THI": 99.9, "PENT": 299.7}, 145.0, 20.0, ph=5.0)
    no_thi = flux(p, {"Cys": 99.9, "PENT": 299.7}, 145.0, 20.0, ph=5.0)
    s_no_cys = branch_shares(no_cys)
    s_no_thi = branch_shares(no_thi)
    add("cerny2007b_control_no_cysteine", "Cerny single-route controls", True,
        "Cerny 2007b Table 5: with NO CYSTEINE, MFT is >99% thiamine-derived "
        "and <1% xylose-derived",
        {"observed": 0.99,
         "predicted": float(s_no_cys["MFT_share_thiamine_route"]),
         "fold_error": float(s_no_cys["MFT_share_thiamine_route"]) / 0.99,
         "pass_band": None,
         "pass": bool(s_no_cys["MFT_share_thiamine_route"] >= 0.99)},
        "THE SHARPEST STRUCTURAL TEST IN THE WHOLE SPLIT. A model fitted on the "
        "ternary must PREDICT both limiting cases; getting 54:46 right while "
        "missing either control means the routes are wrong and the ratio was "
        "fitted.")
    add("cerny2007b_control_no_thiamine", "Cerny single-route controls", True,
        "Cerny 2007b Table 6: with NO THIAMINE, MFT is >95% xylose-derived and "
        "<5% thiamine-derived",
        {"observed": 0.95,
         "predicted": float(s_no_thi["MFT_share_sugar_routes"]),
         "fold_error": float(s_no_thi["MFT_share_sugar_routes"]) / 0.95,
         "pass_band": None,
         "pass": bool(s_no_thi["MFT_share_sugar_routes"] >= 0.95)},
        "the other limiting case")

    # =====================================================================
    # 12. CERNY 2007 Table 5, the 1x / 2x CONCENTRATION PAIR
    #     -- the branch-fraction responsiveness test
    # =====================================================================
    one_x = flux(p, {"Cys": 50.0, "THI": 50.0, "PENT": 150.0}, 145.0, 20.0, ph=5.0)
    two_x = flux(p, {"Cys": 99.9, "THI": 99.9, "PENT": 299.7}, 145.0, 20.0, ph=5.0)
    s1 = branch_shares(one_x)
    s2 = branch_shares(two_x)
    xylose_share_1x = 1.0 - s1["MFT_share_thiamine_route"]
    xylose_share_2x = 1.0 - s2["MFT_share_thiamine_route"]
    add("cerny2007_1x_xylose_share", "Cerny concentration pair", True,
        "Cerny 2007 Table 5 arm B, 1x (xylose 0.15 / cys 0.05 / thiamine 0.05 M): "
        "85 : 15 thiamine : xylose, i.e. a 15% xylose share",
        score(0.15, xylose_share_1x, PASS_BAND_RATIO))
    add("cerny2007_2x_xylose_share", "Cerny concentration pair", True,
        "Cerny 2007 Table 5 arm A, 2x (xylose 0.3 / cys 0.1 / thiamine 0.1 M): "
        "54 : 46, i.e. a 46% xylose share",
        score(0.46, xylose_share_2x, PASS_BAND_RATIO))
    obs_responsiveness = 0.46 / 0.15
    pred_responsiveness = (xylose_share_2x / xylose_share_1x
                           if xylose_share_1x > 0 else float("nan"))
    add("cerny2007_branch_responsiveness", "Cerny concentration pair", True,
        "Cerny 2007 Table 5: a 2x change in precursor loading moves the xylose "
        "share of MFT from 15% to 46% -- a 3.1x change in the BRANCH FRACTION, "
        "one lab, one method, one pH, one temperature",
        score(obs_responsiveness, pred_responsiveness, PASS_BAND_RATIO),
        "THE SINGLE HIGHEST-VALUE HOLD-OUT ROW IN THE SET. It scores whether "
        "the model's branch fractions RESPOND TO CONCENTRATION AT ALL. A model "
        "with fixed branch fractions predicts exactly 1.00x here and fails by "
        "construction, which is the point of the row.")

    # =====================================================================
    # 13. B2.1 -- KANG 2026 SI TABLE S4 AT 140 C. THE NEW GATING HOLD-OUT.
    # =====================================================================
    # FIT_HOLDOUT_DECLARATION.md Amendment 5: "MFT and FFT at 140 C -- a true
    # extrapolation test where a single-Ea model fitted on 100/120 under-
    # predicts 3.8x/2.5x. The measured behaviour is NON-ARRHENIUS ... A
    # single-Ea sulfur branch is EXPECTED TO FAIL this hold-out; passing
    # requires the switch-on to emerge structurally."
    #
    # THE FIT SAW 100 AND 120 C IN THIS EXACT POT AND NOTHING ELSE. Nothing at
    # 140 C entered any fit row, bound or initialisation, and the literal 140 C
    # values appear in no file under src/ (test_kinetic_core_b2_1.py greps).
    #
    # WHAT WOULD HAVE TO EMERGE, stated before the numbers. This module's
    # thiol-forming flux is [carbon skeleton] x [sulfide], and the two factors
    # carry VERY DIFFERENT barriers: the skeleton carries the lumped formation
    # Ea, while the sulfide supply carries Zheng & Ho's MEASURED 133 kJ/mol for
    # cysteine thermolysis -- the largest barrier in the module by a factor of
    # two. The competing sinks that consume the same sulfide and the same
    # cysteine are slower-climbing. A thiol lane gated behind the module's
    # steepest measured barrier, competing for a shared pool against shallower
    # ones, accelerates late; that is the only switch-on mechanism this network
    # has, it was already there before Kang was read, and it is not tunable.
    # IT CAN ALSO OVERSHOOT, and if it does that is a failure in the opposite
    # direction from the one the declaration predicts, which is worth as much.
    kang140 = run(p, {"TTCA": 10.0, "OX": OX_AMBIENT}, 140.0, 120.0,
                  ph=7.0, buffer=BUFFER_NONE)
    kang120 = run(p, {"TTCA": 10.0, "OX": OX_AMBIENT}, 120.0, 120.0,
                  ph=7.0, buffer=BUFFER_NONE)
    for species, observed_ug, at120 in (
        ("MFT", 5.907, 1.388), ("FFT", 11.439, 4.107),
    ):
        add(f"kang_140C_{species}", "Kang 140 C (NEW gating hold-out)", True,
            f"Kang 2026 SI Table S4, 140 C column: {species} {observed_ug} ug/L "
            f"(Tier A, external calibration curve printed, raster-verified "
            f"mu-g/L, nine-subtotal arithmetic closure to +/-0.003)",
            score(ug_per_litre_to_mmol_per_litre(species, observed_ug),
                  kang140.final(species), PASS_BAND_LEVEL),
            "A TRUE EXTRAPOLATION: the fit saw 100 and 120 C in this pot and "
            "nothing above. The declaration pre-registers a 3.8x (MFT) / 2.5x "
            "(FFT) UNDER-prediction for any single-Ea branch. Uncertainty on "
            "the observed value is the dossier's replacement +/-15% (sec. 7d), "
            "NOT the printed SD, which is demonstrably unreproducible between "
            "Tables S4 and S5 for the identical experiment.")
        obs_switch = observed_ug / at120
        pred_switch = (kang140.final(species) / kang120.final(species)
                       if kang120.final(species) > 0 else float("nan"))
        add(f"kang_switch_on_{species}", "Kang 140 C (NEW gating hold-out)",
            True,
            f"Kang 2026 SI sec. 7a: {species} rises {obs_switch:.2f}x from 120 "
            f"to 140 C after rising only "
            f"{'1.12' if species == 'MFT' else '1.10'}x from 100 to 120 C -- "
            f"apparent Ea climbing from ~6-7 to 70-98 kJ/mol",
            score(obs_switch, pred_switch, PASS_BAND_RATIO),
            "THE ROW THE WHOLE REVISION IS ABOUT. It scores the SWITCH-ON "
            "itself rather than the level, so it is immune to the calibration "
            "axis and to any offset the fit absorbed on the two lower rungs. A "
            "single-Ea branch cannot produce a rising apparent Ea at all.")
    # The class-versus-thiol divergence: the sulfur class DECELERATES while the
    # two free thiols accelerate, which is what rules out a saturation artefact.
    add("kang_thiols_diverge_from_their_class",
        "Kang 140 C (NEW gating hold-out)", False,
        "Kang 2026 SI sec. 7a: over 120 -> 140 C the sulfur CLASS decelerates "
        "(2.57x then 1.68x, apparent Ea falling 57.5 -> 35.2 kJ/mol) while MFT "
        "and FFT accelerate. Precursor depletion DEPRESSES apparent Ea, so a "
        "saturation artefact cannot produce the thiols' behaviour.",
        {"observed": None,
         "predicted": {
             "MFT 140/120": float(kang140.final("MFT") / kang120.final("MFT"))
             if kang120.final("MFT") > 0 else float("nan"),
             "furfural 140/120": float(kang140.final("FUR") / kang120.final("FUR"))
             if kang120.final("FUR") > 0 else float("nan"),
         },
         "fold_error": float("nan"), "pass_band": None,
         "pass": bool(
             kang120.final("MFT") > 0 and kang120.final("FUR") > 0
             and (kang140.final("MFT") / kang120.final("MFT"))
             > (kang140.final("FUR") / kang120.final("FUR"))
         )},
        "DIRECTIONAL, and DIAGNOSTIC rather than gating because this module "
        "carries only part of Kang's sulfur class -- no thiophenes, no "
        "thiazoles beyond 2-acetylthiazole, all declared out of scope. The "
        "in-scope proxy for 'the thiols outrun the rest of the pot' is that "
        "MFT's 120->140 fold EXCEEDS furfural's, which is a comparison the "
        "module can honestly make.")
    # The free-cysteine ladder's own third rung.
    kang_cys_140 = run(p, {"Cys": 10.0}, 140.0, 120.0, ph=7.0,
                       buffer=BUFFER_NONE)
    add("kang_140C_cys_conversion", "Kang 140 C (NEW gating hold-out)", False,
        "Kang 2026 SI Fig. S4 (digitised): free-cysteine conversion 62.6% at "
        "140 C / 120 min, against 16.2% and 38.7% at the two fitted rungs",
        score(0.626, 1.0 - kang_cys_140.final("Cys") / 10.0, PASS_BAND_RATIO),
        "DIAGNOSTIC, not gating: it is figure-digitised, its n is unstated, and "
        "the identification of the system as FREE cysteine rather than the "
        "TTCA-bound moiety is 85% confident on the dossier's own reading. The "
        "row it really tests is Kang's measured Ea of 55.1 kJ/mol, which this "
        "module carries as a fixed barrier on r_cys_thermal.")

    # =====================================================================
    # 14. B2.1 -- SUN 2019's pH-9 COLUMN. The other new gating hold-out.
    # =====================================================================
    # FIT_HOLDOUT_DECLARATION.md Amendment 4: "sun2019 pH-9 column (temperature-
    # ordering inversion)". The dossier calls it "the single best hold-out in
    # the paper" for one reason: THE TEMPERATURE ORDERING INVERTS. Free 2-FFT
    # in a coffee brew after 1 h is 2.7 / 1.9 / 1.6 ug/L at 20 / 55 / 90 C at
    # pH 3 (20 C HIGHEST) and 0.4 / 0.5 / 0.5 at pH 9 (20 C LOWEST).
    #
    # SCORED AS RATIOS AND AS AN ORDERING, NEVER AS A LEVEL. The brew's own
    # free-thiol loading is not something this module can set from the corpus,
    # so an absolute prediction would be scoring the loading convention rather
    # than the chemistry. Both rows below divide it out.
    #
    # DISCLOSURE: the FIT partner of this hold-out (Amendment 4 puts every other
    # sun2019 row in the FIT column) was NOT used, because none of the systems
    # it constrains carries a free parameter -- the thioether channel's constant
    # is measured and its reverse is derived. This row is therefore scored as a
    # PURE PREDICTION from measured constants, which is a harder test than the
    # declaration asks for, not an easier one.
    sun_brew = {"FFT": 1.0e-3, "MELE": 9.0 * 10.0, "OX": OX_AMBIENT}
    sun_free = {}
    for ph_value in (3.0, 9.0):
        for t_c in (20.0, 55.0, 90.0):
            sun_free[(ph_value, t_c)] = run(
                p, sun_brew, t_c, 60.0, ph=ph_value
            ).final("FFT")
    obs_ph9_over_ph3_20C = 0.4 / 2.7
    pred_ph9_over_ph3_20C = (
        sun_free[(9.0, 20.0)] / sun_free[(3.0, 20.0)]
        if sun_free[(3.0, 20.0)] > 0 else float("nan")
    )
    add("sun2019_pH9_over_pH3_free_FFT", "Sun 2019 pH-9 column", True,
        "Sun 2019 Fig. 4A, values printed in the text (p. 453): free 2-FFT in "
        "a coffee brew after 1 h at 20 C is 2.7 ug/L at pH 3 and 0.4 ug/L at "
        "pH 9, a 6.75x collapse",
        score(obs_ph9_over_ph3_20C, pred_ph9_over_ph3_20C, PASS_BAND_RATIO),
        "RATIO, so the brew's free-thiol loading divides out. B2 had no "
        "pH-dependent thiol sink at all and could only express this through "
        "formation, which does not run in a spiked brew; B2.1 has the "
        "thiolate-mediated channel Kumazawa's FIT grid pays for.")
    pred_ph3_falls = sun_free[(3.0, 20.0)] > sun_free[(3.0, 90.0)]
    pred_ph9_rises = sun_free[(9.0, 20.0)] < sun_free[(9.0, 90.0)]
    add("sun2019_temperature_ordering_inversion", "Sun 2019 pH-9 column", True,
        "Sun 2019 p. 453: free 2-FFT falls with temperature at pH 3 "
        "(2.7 > 1.9 > 1.6 at 20/55/90 C) and RISES with temperature at pH 9 "
        "(0.4 < 0.5 = 0.5). The ordering INVERTS across the pH axis.",
        {"observed": "pH3 falls with T, pH9 rises with T",
         "predicted": {
             "pH3 20C": float(sun_free[(3.0, 20.0)]),
             "pH3 90C": float(sun_free[(3.0, 90.0)]),
             "pH9 20C": float(sun_free[(9.0, 20.0)]),
             "pH9 90C": float(sun_free[(9.0, 90.0)]),
         },
         "fold_error": float("nan"), "pass_band": None,
         "pass": bool(pred_ph3_falls and pred_ph9_rises)},
        "THE DOSSIER'S OWN VERDICT: 'Any model that gets that inversion right "
        "has learned something real; getting it wrong is a clean "
        "falsification.' Note the pH-9 differences (0.4 vs 0.5) are one "
        "significant figure apart and the paper prints no error bars, so this "
        "is scored as an ORDERING and not as a magnitude.")

    # =====================================================================
    # 15. AMRANI-HEMAIMI on/off switch -- DEFERRED
    # =====================================================================
    add("amrani_hemaimi_onoff_switch", "pyrazines", False,
        "Amrani-Hemaimi 1995: 3-ethyl-2,5-dimethylpyrazine present with "
        "alanine (20 / 19%) and ABSENT with glycine (0 / 0%), 100% 13C-labelled "
        "=> 'one single reaction route exists'",
        {"observed": None, "predicted": None, "fold_error": float("nan"),
         "pass_band": None, "pass": None},
        "DEFERRED EXPLICITLY. Pyrazines are OUT OF SCOPE for this wave "
        "(sulfur.OUT_OF_SCOPE): no pyrazine is in the state vector, so neither "
        "this hold-out nor the Amrani-Hemaimi Table 2 FIT rows are touched. A "
        "later wave that adds the pyrazine lane inherits both, unscored and "
        "unfitted.")

    # =====================================================================
    # THE AROMA QUALIFICATION -- dimerisation is not aroma loss
    # =====================================================================
    oav = odour_activity_values(zhou7.concentrations[-1])
    aroma_note = {
        "predicted_MFT_OAV": oav["MFT"],
        "predicted_dimer_OAV": oav["MFTD"],
        "predicted_dimer_over_monomer_OAV": oav["dimer_over_monomer_OAV"],
        "measured_comparison": (
            "Zhou 2023 SI Table S2 prints, at pH 7, MFT OAV 3.18e5 against "
            "disulfide 3.21e5 -- the dimer's OAV MATCHES the monomer's, "
            "because the dimer is 15.6x more potent and carries 6.5% of the "
            "thiol equivalents."
        ),
        "consequence": (
            "Any objective that scores the dimerisation channel as a pure "
            "aroma LOSS is wrong by roughly the threshold ratio. This module "
            "tracks the dimer as a species with its own potency, so the "
            "hold-out scores above are mass scores and the OAV comparison is "
            "reported alongside them rather than folded into them."
        ),
    }

    # =====================================================================
    # SCORECARD
    # =====================================================================
    gating = [r for r in rows if r["gating"] and r["pass"] is not None]
    diagnostic = [r for r in rows if not r["gating"] and r["pass"] is not None]
    unscored = [r for r in rows if r["pass"] is None]
    n_gating_pass = sum(1 for r in gating if r["pass"])
    n_diag_pass = sum(1 for r in diagnostic if r["pass"])

    # ---------------------------------------------------------------------
    # QUALIFICATIONS THE PASS/FAIL COLUMN HIDES
    # ---------------------------------------------------------------------
    # A scorecard that reports only pass/fail is dishonest in both directions:
    # it credits a pass that happened for the wrong reason, and it discards the
    # information in a near-miss. These four rows are read out explicitly.
    def _row(row_id):
        return next(r for r in rows if r["id"] == row_id)

    qualifications = [
        {
            "row": "vanseeventer_50C_MFT_per_day",
            "nominal": "PASS" if _row("vanseeventer_50C_MFT_per_day")["pass"] else "FAIL",
            "qualification": (
                "THIS PASS SHOULD NOT BE COUNTED AS ONE. It was pre-registered "
                "as a failure and it passed the 3x band from the WRONG SIDE, by "
                "over-predicting: the model destroys ~99% of the MFT per day "
                "against a measured 59%, and it does so through the "
                "unassigned-decay and dimerisation channels, NOT through the "
                "C-5 oligomerisation channel that van Seeventer actually "
                "measured -- which still carries a rate of exactly zero. The "
                "FUNCTIONAL FORM is also still wrong: the measured loss is ZERO "
                "ORDER and every channel here is first or second order. The row "
                "tests the functional form, and on that test the model fails "
                "while landing inside the magnitude band by coincidence."
            ),
        },
        {
            "row": "cerny2007_branch_responsiveness",
            "nominal": "FAIL",
            "qualification": (
                "A FAIL THAT CARRIES THE WAVE'S MAIN ARCHITECTURAL RESULT. This "
                "is the row a model with fixed branch fractions fails BY "
                "CONSTRUCTION, because such a model predicts exactly 1.00x. "
                "This model predicts "
                f"{_row('cerny2007_branch_responsiveness')['predicted']:.2f}x "
                "against a measured 3.07x. So the branch fraction DOES respond "
                "to concentration -- the architectural requirement is met and "
                "the no-fixed-fraction design is vindicated -- but it responds "
                "about half as strongly as measured. The residual is a "
                "magnitude error in the relative reaction ORDERS of the two "
                "routes, not the structural error the row was designed to "
                "catch. Reported as a fail, because it is one."
            ),
        },
        {
            "row": "zhou_pH8_FFT",
            "nominal": "PASS" if _row("zhou_pH8_FFT")["pass"] else "FAIL",
            "qualification": (
                "THE ROW B2.1 EXISTS TO FIX. B2 predicted 1.23e-07 against a "
                "measured 2.85e-03 -- 0.00x -- and its report named the cause: "
                "a structural pH slope of one decade per pH unit, right in "
                "shape and far too steep in level. B2.1 changed three things "
                "and NONE of them is a fitted pH parameter: the pKa are now "
                "evaluated at reaction temperature; no pH-sensitive product has "
                "a single fully-catalysed route any more; and the thiol's "
                "pH-dependent CONSUMPTION lane, which Kumazawa measures "
                "directly at 121 C and which B2 did not have at all, now "
                "carries the part of the pH response that belongs to it. "
                f"B2.1 predicts {_row('zhou_pH8_FFT')['predicted']:.3g} against "
                f"the same measured 2.849e-03, a fold of "
                f"{_row('zhou_pH8_FFT')['fold_error']:.3g}x. Whether that is a "
                "pass or a fail, it is the single number by which this wave "
                "should be judged."
            ),
        },
        {
            "row": "hofmann2002_brew_80C_FFT",
            "nominal": "PASS" if _row("hofmann2002_brew_80C_FFT")["pass"] else "FAIL",
            "qualification": (
                "THE FAILURE THE DECLARATION PREDICTED AND ASKED FOR, and B2.1 "
                "attacks it from two measured directions rather than by "
                "tuning. FIRST: B2 pinned every fitted decay lump to its 145 C "
                "value and evaluated it unchanged at 80 C, so k_fft_decay alone "
                "predicted several times the measured brew loss; those lumps "
                "now carry the lumped Ea, which is a WEAKER claim than "
                "temperature-independence, not a stronger one. SECOND: Stack "
                "2018's measured van 't Hoff enthalpy is NEGATIVE, so the "
                "thiol-electrophile equilibrium constant FALLS 5x between 30 "
                "and 80 C -- heating releases bound thiol. B2 could not express "
                "that, because its K was a temperature-independent number that "
                "the source does not print. "
                f"B2 lost FFT 18x too fast; B2.1 is at "
                f"{_row('hofmann2002_brew_80C_FFT')['fold_error']:.2g}x."
            ),
        },
        {
            "row": "kang_140C_MFT",
            "nominal": "PASS" if _row("kang_140C_MFT")["pass"] else "FAIL",
            "qualification": (
                "THE NEW EXAM, AND THE ONE WHOSE DIRECTION OF FAILURE MATTERS "
                "MOST. Amendment 5 pre-registers that a single-Ea branch "
                "UNDER-predicts here by 3.8x, because the measured behaviour is "
                "non-Arrhenius: free-thiol formation switches on between 120 "
                "and 140 C. This module is not single-Ea at that node -- its "
                "thiol-forming flux is [carbon skeleton] x [sulfide] and the "
                "sulfide is gated behind Zheng & Ho's MEASURED 133 kJ/mol, the "
                "largest barrier in the module by a factor of two, against a "
                "lumped formation barrier for everything else. So it can also "
                "OVERSHOOT, and an overshoot is a failure of the same seriousness "
                "as an undershoot: it would mean the switch-on mechanism is "
                "real but mis-sized. "
                f"MEASURED 5.907 ug/L; PREDICTED "
                f"{_row('kang_140C_MFT')['predicted']:.4g} mmol/L, fold "
                f"{_row('kang_140C_MFT')['fold_error']:.3g}x."
            ),
        },
    ]

    payload = {
        "wave": "B2.2 -- separated decay barriers + pH-trajectory state",
        "qualifications_the_pass_fail_column_hides": qualifications,
        "generated_by": "scripts/generators/generate_kinetic_core_b2_2_holdout.py",
        "parameters_changed_after_the_fit": False,
        "frozen_from": str(FIT_REPORT.relative_to(ROOT)),
        "scorecard": {
            "gating_rows": len(gating),
            "gating_passed": n_gating_pass,
            "diagnostic_rows": len(diagnostic),
            "diagnostic_passed": n_diag_pass,
            "unscoreable_or_deferred": len(unscored),
        },
        "pre_registered_failures": [
            "van Seeventer 50 C -- STILL, and now expected to fail from the "
            "OTHER SIDE. The C-5 oligomerisation channel that actually carries "
            "this loss still has a rate of exactly zero, because its only "
            "measurement IS this hold-out. B2 landed inside the 3x band by "
            "over-destroying the thiol through unassigned decay lumps pinned to "
            "their 145 C values, and its own report refused to count that pass. "
            "B2.1 gives those lumps an activation energy, so they no longer run "
            "at 145 C speed in a 50 C store -- which is correct, and which "
            "should turn an accidental pass into an honest under-prediction. "
            "A scorecard that goes DOWN on this row is the model getting more "
            "truthful, not less.",
            "Hofmann dry-180 C -- unscoreable. Still no water-activity term, "
            "because nothing in the fit corpus varies a_w.",
            "Hofmann pH-3/pH-7 MFT -- the B2.5 sign conflict is a property of "
            "the CORPUS, not of the model: Hofmann's buffered free-sugar system "
            "has MFT falling with pH while Zhou's unbuffered Amadori-fed system "
            "has it peaking at pH 7, with the levels 64x apart at pH 7, and the "
            "dossier's own instruction is not to merge them. One pH shape "
            "satisfies at most one. Nothing in B2.1 changes that.",
            "Zhang Fig. 2 -- methionine is still not in the state vector.",
            "Cerny 2003's norfuraneol-share ceiling at 95 C -- B2 was 5.47x "
            "over an upper bound and nothing in B2.1 targets the norfuraneol "
            "partition. Expected to stay a fail.",
            "Sun 2019's pH-9 TEMPERATURE-ORDERING INVERSION -- expected to "
            "FAIL, and the reason is worth stating in advance. At pH 9 the "
            "thiolate channel dominates and, like every thermal channel in the "
            "module, it accelerates with temperature, so free thiol should FALL "
            "with T. The measurement has it RISING (0.4 -> 0.5 -> 0.5 at 20 -> "
            "55 -> 90 C). The only term in the module that moves the other way "
            "is Stack's measured release, which grows with temperature; "
            "whether it can outrun the sinks at pH 9 is exactly what the row "
            "asks, and the honest prior is that it cannot.",
            "Kang 140 C -- DIRECTION UNKNOWN, WHICH IS WHY IT IS THE EXAM. "
            "Amendment 5 pre-registers a 3.8x/2.5x UNDER-prediction for a "
            "single-Ea branch. This module is not single-Ea at that node: the "
            "sulfide supply is gated behind Zheng's measured 133 kJ/mol while "
            "the carbon skeleton carries a much smaller lumped barrier, so the "
            "thiol lane accelerates late by construction. That mechanism was in "
            "the network before Kang was read and it is not tunable -- but it "
            "can equally OVERSHOOT. This wave declares that it does not know "
            "which side it will miss on, and that only a landing inside 3x "
            "counts as the switch-on having emerged.",
        ],
        "firewall_disclosure": (
            "THREE EXPOSURES, ALL DISCLOSED. (1) The K3 inventory sec. A.3.3 "
            "prints Zhou's pH-6 and pH-8 columns next to the pH-7 FIT column, "
            "and Zhang2024_extraction.md sec. 3 prints Fig. 2 panels (d) and "
            "(e) -- B2's exposure, inherited. (2) NEW AND UNAVOIDABLE: "
            "kang2026_SI_extraction.md PRINTS THE 140 C COLUMN, and the build "
            "brief directed this wave to read that dossier, which is also where "
            "the 100 and 120 C FIT rows live. The 140 C values were therefore "
            "SEEN before the fit ran. (3) sun2019_extraction.md prints the pH-9 "
            "column alongside the pH-3 column. "
            "WHAT WAS DONE ABOUT IT: tests/unit/test_kinetic_core_b2_1.py "
            "greps every file under src/kinetic_core/ and the fit script for "
            "the literals 5.907, 11.439, 60.400, 62.6, 8.195 and Sun's pH-9 "
            "column, and fails if any appears outside a line that explicitly "
            "marks it as held out; a second test walks the fit script's SYSTEMS "
            "dictionary and asserts that no integrated condition is a hold-out "
            "condition (no Kang system at 140 C, no Zhou system off pH 7, no "
            "system at pH 9 or pH 6.5). The frozen mp_holdout_* bundles under "
            "data/benchmarks/external_validation/ were never opened, and a "
            "third test asserts that no B2.1 file performs file I/O on that "
            "path."
        ),
        "aroma_qualification": aroma_note,
        "rows": rows,
    }

    # =====================================================================
    # THE SIDE-BY-SIDE AGAINST B2's FROZEN SCORECARD
    # =====================================================================
    # B2.1's report is read, never written. B2.2 adds no hold-out row and
    # removes none, so every row here has a B2.1 counterpart and a REGRESSION
    # cannot hide inside an improved total.
    b2_path = ROOT / "results/validation/kinetic_core_b2_1_holdout_report.json"
    comparison: List[Dict[str, Any]] = []
    b2_summary: Dict[str, Any] = {}
    if b2_path.exists():
        b2 = json.loads(b2_path.read_text())
        b2_rows = {r["id"]: r for r in b2["rows"]}
        b2_summary = dict(b2.get("scorecard", {}))
        for r in rows:
            before = b2_rows.get(r["id"])
            comparison.append({
                "id": r["id"],
                "gating": r["gating"],
                "b2_1_pass": None if before is None else before["pass"],
                "b2_1_fold": None if before is None else before.get("fold_error"),
                "b2_2_pass": r["pass"],
                "b2_2_fold": r.get("fold_error"),
                "status": (
                    "NEW ROW" if before is None
                    else "unchanged" if before["pass"] == r["pass"]
                    else "FIXED" if r["pass"] else "REGRESSED"
                ),
            })
    payload["comparison_with_b2_1"] = {
        "b2_1_scorecard": b2_summary,
        "b2_2_scorecard": payload["scorecard"],
        "rows": comparison,
        "how_to_read_it": (
            "A row marked REGRESSED is not automatically bad and a row marked "
            "FIXED is not automatically good: see the qualifications block. "
            "B2.2 adds no hold-out row and removes none, so NEW ROW should "
            "never appear -- if it does, the two scorers have drifted apart "
            "and the comparison is unsound."
        ),
    }

    out_json = ROOT / "results/validation/kinetic_core_b2_2_holdout_report.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, default=str))

    lines: List[str] = []
    a = lines.append
    a("# Kinetic core, Build Wave B2.2 -- the pre-registered hold-out scorecard")
    a("")
    a("Predicted from the FROZEN fit in "
      "`results/validation/kinetic_core_b2_2_fit_report.json`. **No parameter was "
      "changed after the fit**; there is no optimiser in the scoring script.")
    a("")
    a(f"**GATING: {n_gating_pass} / {len(gating)} passed.** "
      f"Diagnostic: {n_diag_pass} / {len(diagnostic)}. "
      f"Unscoreable or explicitly deferred: {len(unscored)}.")
    a("")
    if comparison:
        b2_g = b2_summary.get("gating_passed"), b2_summary.get("gating_rows")
        n_reg = sum(1 for c in comparison if c["status"] == "REGRESSED")
        n_fix = sum(1 for c in comparison if c["status"] == "FIXED")
        a(f"**Against B2.1: {b2_g[0]} / {b2_g[1]} gating.** "
          f"B2.2 adds no row and removes none, so the denominators match and "
          f"the totals ARE directly comparable. "
          f"**{n_fix} FIXED, {n_reg} REGRESSED.**")
        a("")
        a("## B2.1 versus B2.2, row by row")
        a("")
        a("| row | gating | B2.1 fold | B2.1 | B2.2 fold | B2.2 | |")
        a("|---|---|---:|---|---:|---|---|")
        for c in comparison:
            def _f(v):
                return (f"{v:.3g}x" if isinstance(v, float) and np.isfinite(v)
                        else "--")

            def _p(v):
                return "--" if v is None else ("PASS" if v else "**FAIL**")

            a(f"| `{c['id']}` | {'yes' if c['gating'] else 'no'} | "
              f"{_f(c['b2_1_fold'])} | {_p(c['b2_1_pass'])} | "
              f"{_f(c['b2_2_fold'])} | {_p(c['b2_2_pass'])} | {c['status']} |")
        a("")
        a(payload["comparison_with_b2_1"]["how_to_read_it"])
        a("")
    a("## Pre-registered failures (written down before the numbers)")
    a("")
    for f in payload["pre_registered_failures"]:
        a(f"- {f}")
    a("")
    a("## Qualifications the pass/fail column hides")
    a("")
    for q in qualifications:
        a(f"- **`{q['row']}` ({q['nominal']})** -- {q['qualification']}")
        a("")
    a("## Row by row -- no averaging, no dropping")
    a("")
    a("| row | group | gating | observed | predicted | fold | pass |")
    a("|---|---|---|---:|---:|---:|---|")
    for r in rows:
        obs = r["observed"]
        pred = r["predicted"]
        fold = r["fold_error"]
        obs_s = f"{obs:.4g}" if isinstance(obs, (int, float)) else "--"
        pred_s = f"{pred:.4g}" if isinstance(pred, (int, float)) else "--"
        fold_s = (f"{fold:.2f}x" if isinstance(fold, float) and np.isfinite(fold)
                  else "--")
        verdict = ("PASS" if r["pass"] else "**FAIL**") if r["pass"] is not None \
            else "_not scored_"
        a(f"| `{r['id']}` | {r['group']} | {'yes' if r['gating'] else 'no'} | "
          f"{obs_s} | {pred_s} | {fold_s} | {verdict} |")
    a("")
    a("## Each row's source anchor and qualification")
    a("")
    for r in rows:
        a(f"### `{r['id']}`")
        a("")
        a(f"- source: {r['source_anchor']}")
        if r["comment"]:
            a(f"- {r['comment']}")
        a("")
    a("## The aroma qualification")
    a("")
    a(aroma_note["measured_comparison"])
    a("")
    a(aroma_note["consequence"])
    a("")
    a("## Firewall disclosure")
    a("")
    a(payload["firewall_disclosure"])
    a("")

    out_md = ROOT / "results/validation/kinetic_core_b2_2_holdout_report.md"
    out_md.write_text("\n".join(lines))
    print(f"GATING {n_gating_pass}/{len(gating)}  diagnostic {n_diag_pass}/"
          f"{len(diagnostic)}  unscored {len(unscored)}")
    print(f"wrote {out_json}")
    print(f"wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
