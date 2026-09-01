from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from typing import Dict, Optional

from src import data_paths
from src import compound_keys

_RUNTIME_MULTIPLIER_ENV = "MAILLARD_MATRIX_CALIBRATION_MULTIPLIERS"


# 2026-08-27 (Wave I) — the `fitted_to_benchmark` evidence strength.
#
# The cold-start red team (forensic F3, scientific C3) found that several entries in this
# registry were labelled `literature_anchored` / `conditional_literature_anchored` when they
# were in fact BACK-SOLVED from the very benchmark they are then scored against. The tell was
# an 8-to-17-significant-figure constant sitting next to a benchmark row reporting
# Pearson R = 1.000 and max_ratio = 1.000.
#
# THE ARITHMETIC, so this is checkable rather than asserted:
#
#   Pea reference lane (`MATRIX_BENCHMARK_BASE_MARKER_YIELDS` in benchmark_validation.py)
#     Hexanal        0.205 x 1268.3 = 260 ppb  <- pea_isolate_40C_PratapSingh2021 measured
#     2-Pentylfuran  0.502 x 1270.9 = 638 ppb  <- same benchmark
#     1-Hexanol      0.063 x 1269.8 =  80 ppb  <- same benchmark
#   One common scale (1269 +/- 0.2%) recovers all three measured values exactly. Those
#   "yields" ARE the benchmark's own measurements, rescaled by a single constant.
#
#   Soy ambient lane (the `0.453 / 0.205`-style expressions below)
#     numerators 0.453 / 2.972 / 0.143 x 838.8 = 380 / 2492 / 120 ppb
#                                                <- soy_isolate_40C_PratapSingh2021 measured
#   Again one common scale (838.8 +/- 0.08%) across all three. Same construction.
#
#   Trikusuma pea-UHT heated lane: three 17-significant-figure constants whose benchmark row
#   scores max_ratio 1.000 / MALE 0.000 on all three analytes. That is what an exact
#   algebraic recovery looks like; no measurement is reported to 17 figures.
#
# None of this is wrong to DO -- a matrix observability factor has to be anchored to
# something, and the anchor is a real measurement from a real paper. What was wrong was the
# LABEL: calling it `literature_anchored` and then reporting the resulting agreement as
# validation. A lane whose constant was solved from a benchmark cannot also be evidence about
# that benchmark. Such lanes are now labelled `fitted_to_benchmark`, and the benchmark
# reporting layer surfaces them as "fit-recovery (not predictive)" instead of "pass".
#
# 2026-08-27 (Wave M) -- AND THE ANCHOR ITSELF WAS WRONG. The arithmetic above is unchanged
# and still describes exactly what these constants are, but the sentence "the anchor is a
# real measurement from a real paper" no longer holds for the two hexanal rows and does not
# hold at all for the 1-hexanol rows. Wave K read Molecules 2021, 26, 4104 Table 1 (Europe
# PMC, PMC8271896): the paper's pea hexanal is 1138.00 ppb and its soy hexanal is 1621.71
# ppb -- 4.38x and 4.27x above the 260 / 380 these factors were back-solved from -- and it
# reports n.d. for hexanol in BOTH matrices, so the 80 / 120 ppb had no source at all.
# The benchmark files are corrected (see their `content_correction_note`); THESE FACTORS ARE
# NOT, and that is deliberate: refitting them is a science decision for an owner, not a
# propagation step, and doing it in the same pass as a chemistry change would make the
# agreement unattributable.
#
# CONSEQUENCE, so nobody reads the current numbers as agreement: the pea and soy ambient
# hexanal lanes now UNDER-predict by exactly the correction factor (predicted 260.6 vs 1138
# measured, 379.9 vs 1621.71), because the factor still reproduces the erroneous value it
# was solved from. The 2-pentylfuran factors are unaffected -- 638 and 2492 were verified
# verbatim against the paper. The 1-hexanol entries now anchor to a compound the paper says
# was not detected; they no longer have a target at all.
#
# 2026-08-27 (Wave O) -- THE REFIT, OWNER-APPROVED. The two ambient hexanal factors have now
# been refitted against the VERIFIED values (1138.00 / 1621.71 ppb).
#   Generator / record: scripts/generators/refit_matrix_observability_pratap_singh.py ->
#                       results/validation/matrix_observability_refit_pratap_singh.{json,md}
#   Fitted: ONE parameter, a common multiplicative scale s = 4.317249 on the ambient-slurry
#           HEXANAL factors of BOTH lanes. Nothing else moved -- no barrier, no projection
#           constant, no marker yield, and NOT the 2-pentylfuran factors (638 / 2492 were
#           verified verbatim and are already recovered to 4e-4).
#   Why one parameter and not two: two factors against two rows is exactly determined, so its
#           zero residual is arithmetic and says nothing. One shared scale leaves the data a
#           degree of freedom to disagree with, and it very nearly does not: the pea lane
#           wants 4.36606x, the soy lane 4.26899x, and the shared 4.317249x satisfies both to
#           1.0113x. The correction is an ABSOLUTE-SCALE error; the pea-vs-soy release ratio
#           this registry encodes (2.2098) survives it (a per-lane fit would make it 2.1606).
#   Status UNCHANGED: both benchmarks stay `fitted_to_benchmark`, stay out of the honest
#           literature-coverage numerator and denominator, and are reported as fit-recovery.
#           A refit changes the SIZE of a recovery (4.37x/4.27x under -> 1.0113x), never its
#           evidential status.
#   HONEST CONSEQUENCE, stated here because it is the number that matters: the EXTERNAL
#           HOLD-OUT gets WORSE. The pea ambient lane carries two mutually contradictory
#           external measurements (Bi 2020 1260 ppb vs Liu 2023 band midpoint 51.96 ppb, a
#           24x spread at nominally identical conditions), and the erroneous 260 ppb constant
#           sat almost exactly at their geometric mean (sqrt(1260 x 51.96) = 255.9). Moving to
#           the verified anchor moves the lane onto the Bi side of that disagreement: Bi raw
#           pea 5.37x -> 1.24x, but Liu 2023 4.52x -> 19.5x and Li 2026 hexanal 21.6x -> 93x.
#           Hold-out median |log10| 1.1849 -> 1.6298 dex (15.31x -> 42.64x). See
#           tasks/audit_remediation.md, Wave O.
#   NOT REFITTED, because there is nothing left to fit them to: the 1-hexanol factors. The
#           paper reports n.d. in both matrices, so `0.143 / 0.063` is a ratio of two numbers
#           that appear nowhere in it. Left untouched and flagged; retiring it is a separate
#           science decision.
#
# 2026-08-27 (Wave T3) -- WHERE 260 / 380 / 80 / 120 ACTUALLY CAME FROM. Waves K, M and O all
# recorded these as "origin unknown, no derivation found". They are not unknown. All four are
# printed in `data/benchmarks/maillard_validation_benchmarks.md` section 3.1 as `~260 +/- 35`,
# `~380 +/- 42`, `~80` and `~120` -- the repository's own abstract-reconstructed internal
# brief, the same file and the same commit era as its section 1.3, which Wave S2b proved
# fabricated (the MFT 342 / FFT 200 case). Section 3.1 carries the identical forensic
# fingerprint, and it is now validated twice independently: BOLD or unhedged = transcribed
# from a real abstract, `~`-hedged = invented. Section 3.1's only two bold cells are the
# 2-pentylfuran `638 +/- 49` and `2492 +/- 199` that Wave K confirmed VERBATIM; its
# tilde-hedged cells are exactly the four Wave K found wrong or unsourced. The fingerprint
# predicts Wave K's result with no misses. So these constants were not mis-transcribed from a
# paper -- they were back-solved from numbers this repository generated about itself.
#
# WHAT THAT CHANGES, AND WHAT IT DOES NOT. It changes nothing about the two hexanal factors:
# Wave O already refitted them onto the verified 1138.00 / 1621.71 ppb, and that stands. It
# changes the STATUS of the two 1-hexanol factors, which Wave O correctly declined to refit:
# they are not merely "unanchored", they are solved from a fabricated pair, and they are LIVE
# -- `0.143/0.063` ships, and the li_2026_hme 1-hexanol hold-out row misses by 1117x. Both
# entries are now labelled `no_verifiable_source` in `source` and carry the full trail in
# `notes`. WAVE T3 DID NOT REFIT OR RETIRE THEM. Refitting needs a measurement that does not
# exist; retiring the lane changes predictions and the hold-out. Which of the two to do is an
# OWNER DECISION, carried as [P] in tasks/audit_remediation.md.
#
# 2026-08-28 (Wave Y) -- THE AMBIENT HEXANAL SCALE HAS MOVED OUT OF THIS FILE. Wave O's
# shared scale 4.317249 is not withdrawn and no parameter was added or removed; it now lives
# in `MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal']` (0.205 -> 0.885036) and every hexanal
# `observable_factor` here has been divided by it:
#
#     pea_iso / ambient_slurry   4.31725    ->  1.0                        (reference lane)
#     pea_iso / heated_matrix    0.228776   ->  0.0529912                  (propagation)
#     soy_iso / ambient_slurry   9.54007    ->  0.453 / 0.205              (the soy/pea ratio)
#     soy_iso / heated_matrix    2.80478    ->  (0.453/0.205)*(1-0.7060)   (composition rule)
#
# THE ARGUMENT IS A UNIT ARGUMENT (Wave S4 (b)). An `observable_factor` is the fraction of a
# total that a measurement sees, so it cannot exceed 1 -- and after Wave O this file shipped
# 4.32 and 9.54. A marker YIELD has no such bound: it multiplies `hydroperoxide_scale`
# (1.0e6, an arbitrary constant in data/lit/lipid_oxidation_calibration.json), so only the
# product is identified and a yield above 1 says nothing. The scale therefore belongs on the
# yield side, and Wave S4 (c) independently evidences the value this lane returns to:
# Pratap-Singh spike their standards INTO the slurry, i.e. matrix-matched quantification,
# i.e. the measurement reads the TOTAL and the observable fraction is 1.0.
#
# NOTHING PREDICTED MOVED. The product `Y * cal` is preserved on every lane to 6 significant
# figures, verified by re-running the decomposition before and after; the hold-out's eight
# points and the four synthetic snapshots are unchanged. What DOES move is the UNCALIBRATED
# tier (`_uncalibrated_prediction_ppb`), which reads the yield and never reads this file --
# see results/validation/matrix_sigma_residual_derivation.{json,md}. Wave O wrote that no
# refit of THESE constants could ever move that derivation; that remains true, and a yield
# refit is the case it does not cover.
#
# AND HALF OF WAVE S4's PREDICTION IS FALSIFIED, which is recorded here rather than in a
# report nobody opens. S4 expected the factors to come back under 1 once the yields were
# fixed. Two do. SIX DO NOT, and every one of them is a SOY factor: ambient hexanal 2.2098,
# ambient 2-pentylfuran 5.9203, ambient 1-hexanol 2.2698, ambient nonanal 1.0667, and the
# soy class anchors 2.209 (aldehyde) / 5.92 (furan). A marker yield is shared across
# matrices, so it can absorb a GLOBAL scale error and can never absorb a LANE one. Measured
# with observability pinned to its evidenced 1.0 on both ambient lanes, the soy-vs-pea
# required-yield ratio is 2.1606x on hexanal and 5.9221x on 2-pentylfuran -- and those two
# differ, so the deficit is COMPOUND-specific too and cannot be repaired by the soy lipid
# profile either. That is the next workstream. [P]
#
# Derivation, corpus, identifiability and the full before/after census:
#   scripts/generators/rederive_matrix_marker_yields.py
#   results/validation/matrix_marker_yield_rederivation.{json,md}
FITTED_TO_BENCHMARK = "fitted_to_benchmark"


@dataclass(frozen=True)
class MatrixCalibrationRecord:
    protein_type: str
    process_state: str
    compound: str
    observable_factor: float
    evidence_strength: str
    source: str
    fallback_mode: str
    notes: str = ""
    # Set when a shipped constant is re-expressed (e.g. rounded from fit residue to a
    # defensible precision). Retains the exact prior value so the change is auditable.
    previous_value: Optional[float] = None
    # When evidence_strength == FITTED_TO_BENCHMARK, the benchmark id the factor was
    # solved from. Consumed by the reporting layer to mark the row fit-recovery.
    fitted_from_benchmark: Optional[str] = None

@dataclass(frozen=True)
class MatrixCalibrationAnchor:
    protein_type: str
    target_class: str
    observable_factor: float
    evidence_strength: str
    source: str
    fallback_mode: str
    notes: str = ""
    previous_value: Optional[float] = None
    fitted_from_benchmark: Optional[str] = None

@dataclass(frozen=True)
class MatrixRuntimeCompositionRule:
    protein_type: str
    compound: str
    mode: str
    active_process_states: tuple[str, ...]
    source: str
    notes: str = ""


_MATRIX_CALIBRATION_RECORDS = (
    # --- Pea ambient reference lane -------------------------------------------------
    # 2026-08-27 (Wave I): relabelled literature_anchored -> fitted_to_benchmark.
    # These factors were 1.0 BY CONSTRUCTION: this lane is the reference against which the
    # others are expressed, and the yields it multiplies
    # (benchmark_validation.MATRIX_BENCHMARK_BASE_MARKER_YIELDS) are
    # pea_isolate_40C_PratapSingh2021's own measured ppb divided by one common scale of
    # 1268-1271. So `pea_isolate_40C_PratapSingh2021` scoring Pearson 1.000 / max_ratio 1.002
    # was the lane reproducing the numbers it was built from, not a prediction.
    # 2026-08-27 (Wave O): the hexanal factor is no longer 1.0 -- the "1268-1271" scale was
    # built from a hexanal value that is not in the paper, and the refit (see the module
    # header) corrects the absolute scale by 4.317249x. The 2-pentylfuran / 1-hexanol /
    # nonanal factors of this lane stay at 1.0, so the lane's INTERNAL reference role is
    # preserved for every compound whose anchor was verified or absent.
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        # Wave O refit, owner-approved: 1.0 -> 4.31725 (x4.317249, the shared ambient-hexanal
        # observability scale fitted against the CONTENT-VERIFIED 1138.00 ppb of Molecules
        # 2021, 26, 4104 Table 1 -- not the 260 ppb transcription error the old 1.0 recovered).
        # Record: results/validation/matrix_observability_refit_pratap_singh.{json,md};
        # generator: scripts/generators/refit_matrix_observability_pratap_singh.py (2026-08-27).
        #
        # 2026-08-28 (Wave Y) -- THE SCALE MOVED SIDES. 4.31725 -> 1.0, and the SAME constant
        # now lives in MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal'] (0.205 -> 0.885036).
        # Wave O's fit is not withdrawn and no parameter was added: one shared scale, the
        # same two verified anchors, the same 1.0113x residual. What changed is which side of
        # the product carries it, and the argument is a UNIT argument (Wave S4 (b)): an
        # observability factor is the fraction of a total that a measurement sees and cannot
        # exceed 1, so "4.32" was not an observability at all. A marker yield can exceed 1
        # because it multiplies an arbitrary `hydroperoxide_scale`. This lane's factor is
        # therefore back at the 1.0 that its OWN definition requires -- it is the reference
        # lane -- and, independently, at the value Wave S4 (c) evidenced from the paper's
        # verbatim methods: Pratap-Singh spike their standards into the slurry, i.e. the
        # measurement is MATRIX-MATCHED and reads the TOTAL, i.e. the observable fraction
        # is 1.0. The two arguments agree, which is why this is the value that ships.
        # Derivation: scripts/generators/rederive_matrix_marker_yields.py ->
        # results/validation/matrix_marker_yield_rederivation.{json,md}.
        observable_factor=1.0,
        previous_value=4.31725,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline (Wave O refit to the verified 1138.00 ppb)",
        fallback_mode="compound_specific",
        notes="Reference compound for the pea matrix-only intake/headspace lane. FITTED: the factor is a shared ambient-hexanal observability scale solved so this benchmark's own verified measurement is reproduced, so this lane is fit-recovery, not a prediction. Refitting changed the size of the recovery (4.37x under -> 1.0113x), not its evidential status.",
        fitted_from_benchmark="pea_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="2-pentylfuran",
        observable_factor=1.0,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
        notes="Reference furan marker for the pea matrix-only intake/headspace lane. Factor is 1.0 by construction; see the hexanal entry.",
        fitted_from_benchmark="pea_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="1-hexanol",
        observable_factor=1.0,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="no_verifiable_source (2026-08-27, Wave T3) -- back-solved from a fabricated value; was 'Pratap-Singh 2021 pea isolate ambient slurry baseline'",
        fallback_mode="compound_specific",
        notes=(
            "NO VERIFIABLE SOURCE (2026-08-27, Wave T3; findings T1-04/T1-05). This factor is "
            "1.0 by construction because it DEFINES the pea reference lane: the module header's "
            "arithmetic is 0.063 x 1269.8 = 80 ppb, and the soy 1-hexanol factor below is the "
            "ratio 0.143/0.063 built on the same 80 / 120 ppb pair. Pratap-Singh et al. "
            "(Molecules 2021, 26, 4104, Table 1) report n.d. for hexanol in BOTH matrices and "
            "state pea proteins 'contained no alcohol compounds', so this lane was never "
            "anchored to that paper. Wave T3 identified where 80 and 120 DID come from: "
            "data/benchmarks/maillard_validation_benchmarks.md section 3.1, tilde-hedged rows "
            "'~80' and '~120' -- the repository's own abstract-reconstructed brief, same table "
            "and same fingerprint as the section 1.3 fabrication Wave S2b settled. These are "
            "LIVE CONSTANTS and this lane carries the hold-out's worst miss (li_2026_hme "
            "1-hexanol, 1117x). REFITTING OR RETIRING THE 1-HEXANOL LANE IS AN OWNER DECISION "
            "[P]. Wave T3 deliberately did NOT refit: there is no measurement to fit to, and "
            "inventing one would repeat the defect that produced 80 and 120 in the first place."
        ),
        fitted_from_benchmark="pea_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="nonanal",
        observable_factor=1.0,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
        fitted_from_benchmark="pea_isolate_40C_PratapSingh2021",
    ),
    # --- Trikusuma pea-UHT heated lane -----------------------------------------------
    # 2026-08-27 (Wave I): relabelled conditional_literature_anchored ->
    # fitted_to_benchmark, and the 17-significant-figure fit residue rounded to 6
    # significant figures (previous_value retained). These three constants were BACK-SOLVED
    # so that pea_isolate_uht_140C_Trikusuma2019 reproduces its own measured 782 / 163 / 24
    # ppb; the benchmark row consequently reports Pearson 1.000, max_ratio 1.000,
    # MALE 0.000. That is fit recovery, not agreement.
    # Rounding changes each factor by < 1e-6 relative -- verified not to move any scored
    # benchmark row.
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="hexanal",
        # 2026-08-28 (Wave Y): 0.228776 -> 0.0529912. PROPAGATION, NOT A NEW FIT. The
        # back-solve that produced 0.228776 is unchanged in every respect except which side
        # of `Y * cal` carries the ambient hexanal scale; dividing by the same 4.317249 that
        # went into MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal'] preserves this lane's
        # recovery of its own 782 ppb exactly (measured: 782.0056 ppb before and after).
        # Trikusuma 2020 is still the last content-unverified pillar of the matrix lane and
        # this row is still fit recovery, not evidence -- neither of those changes here.
        observable_factor=0.0529912,
        previous_value=0.228776,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Trikusuma 2020 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT aldehyde anchor carried onto the matrix-only oxidation/headspace lane. BACK-SOLVED from this benchmark's own measured 782 ppb -- a process-state-specific observable correction, not a global oxidation law, and not independent evidence about this benchmark.",
        fitted_from_benchmark="pea_isolate_uht_140C_Trikusuma2019",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="2-pentylfuran",
        observable_factor=0.0194733,
        previous_value=0.019473307397293472,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Trikusuma 2020 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT furan anchor. BACK-SOLVED from this benchmark's own measured 163 ppb.",
        fitted_from_benchmark="pea_isolate_uht_140C_Trikusuma2019",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="nonanal",
        observable_factor=0.00959565,
        previous_value=0.009595650239086601,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Trikusuma 2020 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        # 2026-08-27 (Wave P item 4) — THIS CONSTANT NO LONGER RECOVERS ITS OWN ANCHOR,
        # AND IT WAS NOT REFITTED. It was back-solved while nonanal was (wrongly) scaled
        # off the LINOLEATE hydroperoxide pool. Nonanal is the C9 fragment of the OLEATE
        # double bond and is now scaled off `oleic_acid_pct` (see
        # src/lipid_oxidation.MARKER_HYDROPEROXIDE_POOL; Miyazaki 2023 10.1093/bbb/zbac189
        # finds nonanal in NEITHER linoleate hydroperoxide isomer's product list). The
        # prediction consequently fell from 24.00 ppb to 10.56 ppb, i.e. 2.2727x under —
        # which is EXACTLY 1 / (oleic 22.0 / linoleic 50.0) = 1/0.44. That exact
        # arithmetic identity is the fingerprint: the factor was absorbing a
        # substrate-assignment error, nothing else. Refitting it here would re-absorb the
        # correction into the same constant and make the fix invisible, so it is left at
        # its old value and the miss is reported. Trikusuma 2020 is also still the last
        # content-unverified pillar of the matrix lane (Wave O [P] item 5), so there is no
        # verified anchor to refit against even if refitting were wanted.
        notes="Heated pea UHT nonanal anchor. BACK-SOLVED from this benchmark's own measured 24 ppb, against the pre-Wave-P linoleate-pool nonanal. NOT refitted after the Wave P oleate substrate correction: the row now reads 2.2727x under, which is exactly the oleic/linoleic ratio.",
        fitted_from_benchmark="pea_isolate_uht_140C_Trikusuma2019",
    ),
    # --- Soy ambient lane ------------------------------------------------------------
    # 2026-08-27 (Wave I): relabelled literature_anchored -> fitted_to_benchmark.
    # The expressions read like literature ratios, and the DENOMINATORS are the pea
    # reference yields -- but the NUMERATORS (0.453 / 2.972 / 0.143 / 0.160) are
    # soy_isolate_40C_PratapSingh2021's own measured ppb (380 / 2492 / 120) divided by one
    # common scale of 838.5-839.2. So both halves of the ratio are the two benchmarks these
    # lanes are then scored against, which is why soy scores max_ratio 1.001 / Pearson 1.000.
    # The expressions are LEFT AS EXPRESSIONS on purpose: they show the construction.
    # 2026-08-27 (Wave O): the hexanal expression is GONE, because its numerator 0.453 was
    # built from the 380 ppb transcription error. It is replaced by an explicit fitted
    # constant; the other three keep their expressions, which are still exactly what they
    # were.
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        # Wave O refit, owner-approved: 0.453/0.205 = 2.2097561 -> 9.54007 (x4.317249, the
        # SAME shared ambient-hexanal observability scale as the pea lane), fitted against the
        # CONTENT-VERIFIED 1621.71 ppb of Molecules 2021, 26, 4104 Table 1. Because one scale
        # serves both lanes, the soy-vs-pea release ratio is unchanged at 2.2097561 -- the
        # correction was to the absolute scale, not to the relative structure.
        # Record: results/validation/matrix_observability_refit_pratap_singh.{json,md}.
        #
        # 2026-08-28 (Wave Y) -- 9.54007 -> 0.453/0.205, the shared scale relocated into
        # MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal'] (see the pea ambient entry above and
        # results/validation/matrix_marker_yield_rederivation.md). The expression is restored
        # verbatim BECAUSE IT IS SELF-DOCUMENTING: this constant is the soy-vs-pea ratio and
        # nothing else, and writing it as a ratio is what makes that checkable.
        #
        # AND IT IS STILL ABOVE 1, WHICH FALSIFIES HALF OF WAVE S4's DIAGNOSIS. S4 predicted
        # the factors would come back under 1 once the yields were fixed. The pea lane does
        # (4.31725 -> 1.0) and the soy heated lane does (2.80478 -> 0.6497); this one does
        # not, and neither do the four other soy factors, which the relocation never touched.
        # The reason is structural and is measured in the Wave Y record: a marker yield is
        # shared across matrices, so it can absorb a GLOBAL scale error and can never absorb
        # a LANE one. With observability pinned to its evidenced 1.0 on both ambient lanes,
        # the soy/pea required-yield ratio is 2.1606x on hexanal and 5.9221x on
        # 2-pentylfuran -- and because those two differ from each other, the residual is
        # COMPOUND-specific as well, so it cannot be repaired by the soy lipid profile
        # either (that would move both linoleate markers by one factor). The soy-vs-pea
        # marker-ratio deficit is the next workstream, and it is not in this one. [P]
        observable_factor=0.453 / 0.205,
        previous_value=9.54007,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 soy ambient slurry hexanal (Wave O refit to the verified 1621.71 ppb)",
        fallback_mode="compound_specific",
        notes="Soy ambient hexanal observability. FITTED to this benchmark's own verified measurement via a scale shared with the pea lane -- fit-recovery, not a prediction. Residual after the shared-scale fit is 1.0113x, and that residual exists only because the second degree of freedom was declined.",
        fitted_from_benchmark="soy_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="2-pentylfuran",
        observable_factor=2.972 / 0.502,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
        fitted_from_benchmark="soy_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="1-hexanol",
        observable_factor=0.143 / 0.063,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="no_verifiable_source (2026-08-27, Wave T3) -- back-solved from a fabricated value; was 'Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio'",
        fallback_mode="compound_specific",
        notes=(
            "NO VERIFIABLE SOURCE (2026-08-27, Wave T3; findings T1-04/T1-05). 0.143/0.063 is "
            "the ratio of two numbers that appear NOWHERE in Pratap-Singh et al. (Molecules "
            "2021, 26, 4104): its Table 1 reports n.d. for hexanol in both matrices, and the "
            "paper's entire soy alcohol fraction is 40 +/- 9 ppb of 1-octen-3-ol -- one third "
            "of the 120 ppb this numerator was solved from. The 0.143 x 838.8 = 120 and 0.063 x "
            "1269.8 = 80 arithmetic is in the module header. Wave T3 traced 80 and 120 to "
            "data/benchmarks/maillard_validation_benchmarks.md section 3.1, tilde-hedged rows "
            "'~80' / '~120', in the repository's own abstract-reconstructed brief -- the same "
            "table and the same bold-vs-tilde fingerprint as the section 1.3 fabrication. The "
            "ratio is therefore a ratio of two INVENTED numbers, and it ships. REFIT OR RETIRE "
            "IS AN OWNER DECISION [P]; Wave T3 did not refit, because there is nothing to fit "
            "to and substituting a plausible number would be the original defect again."
        ),
        fitted_from_benchmark="soy_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="nonanal",
        observable_factor=0.160 / 0.150,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
        notes="Nonanal is NOT reported for soy in the ambient benchmark; 0.160 is a lane-internal value carried on the same construction as the other three, so it is neither literature-anchored nor independently checkable.",
        fitted_from_benchmark="soy_isolate_40C_PratapSingh2021",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="heated_matrix",
        compound="hexanal",
        # 2026-08-27 (Wave O): PROPAGATED, not fitted. This constant is DEFINED as the soy
        # ambient hexanal baseline times the Shu 2024 attenuation, so when the baseline was
        # refitted the composition had to follow -- freezing it would leave a corrected fit
        # composed with a stale baseline, which is the exact defect Wave O exists to remove.
        # 0.649668 -> 2.80478 (the same x4.317249). No benchmark constrains this lane; the
        # movement is a propagation of someone else's fit, and it is the mechanism by which
        # the Li 2026 hold-out hexanal point degrades from 21.6x to 93x over.
        #
        # 2026-08-28 (Wave Y): 2.80478 -> 0.6496683, following its baseline back down as the
        # shared scale relocated into the marker yield. The composition rule is unchanged --
        # this constant is still DEFINED as (soy ambient baseline) x (1 - 0.7060) -- so the
        # product that reaches a prediction is preserved and the Li 2026 hold-out hexanal
        # point does not move. Note what that means: this factor was above 1 and now is not,
        # so the Wave S4 unit objection is answered on THIS lane by arithmetic that changed
        # no prediction at all. See results/validation/matrix_marker_yield_rederivation.md.
        observable_factor=(0.453 / 0.205) * (1.0 - 0.7060),
        previous_value=2.8047805800000005,
        evidence_strength="conditional_literature_anchored",
        source="Shu 2024 heated soy off-flavour attenuation carried onto the Pratap-Singh soy ambient baseline",
        fallback_mode="compound_specific_process_state",
        notes="High-severity soy treatment prior for heated matrix states. Useful for reliability and directional accuracy, but not a meaty benchmark anchor. 2026-08-27 (Wave I): kept as conditional_literature_anchored -- the ATTENUATION (1 - 0.7060) is a real Shu 2024 literature figure and no panel benchmark exists for heated soy, so this lane is not fit-recovery. But note its BASELINE factor is fit-recovery (see the soy ambient block), so the anchoring is only as strong as one literature attenuation applied to a back-solved base. 2026-08-27 (Wave O): baseline refitted 2.2097561 -> 9.54007, so this composed value moved with it.",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="heated_matrix",
        compound="2-pentylfuran",
        observable_factor=(2.972 / 0.502) * 0.03,
        evidence_strength="conditional_literature_anchored",
        source="Shu 2024 heated soy 2-pentylfuran below-detection carryover onto the Pratap-Singh soy ambient baseline",
        fallback_mode="compound_specific_process_state",
        notes="Below-detection is carried as a small non-zero censoring surrogate so heated soy ranking stays numerically stable while honoring the reported severe attenuation. 2026-08-27 (Wave I): same reading as the heated soy hexanal entry -- literature attenuation on a fit-recovery baseline.",
    ),
)


def fitted_to_benchmark_lanes() -> Dict[str, tuple[str, ...]]:
    """Benchmark id -> the compounds whose observable factor was solved FROM it.

    2026-08-27 (Wave I).  A benchmark listed here cannot be scored as evidence about the
    model: the constants under test were derived from its own measured values, so any
    agreement is algebraic recovery.  Consumers (`benchmark_validation`, the summary
    generators, `scripts/ci/fit_target_gate.py`) use this to label such rows
    ``fit_recovery`` instead of ``pass`` and to exclude them from literature-coverage
    numerators and denominators.
    """
    lanes: Dict[str, list[str]] = {}
    for record in _MATRIX_CALIBRATION_RECORDS:
        if record.evidence_strength != FITTED_TO_BENCHMARK:
            continue
        if not record.fitted_from_benchmark:
            continue
        lanes.setdefault(record.fitted_from_benchmark, []).append(record.compound)
    for anchor in _MATRIX_CLASS_ANCHORS:
        if anchor.evidence_strength != FITTED_TO_BENCHMARK:
            continue
        if not anchor.fitted_from_benchmark:
            continue
        lanes.setdefault(anchor.fitted_from_benchmark, []).append(anchor.target_class)
    return {key: tuple(sorted(set(value))) for key, value in sorted(lanes.items())}


def is_fit_recovery_benchmark(benchmark_id: Optional[str]) -> bool:
    """True when every scored row of this benchmark is algebraic recovery of a fit.

    See `fitted_to_benchmark_lanes`.  Used to keep fit-recovery rows out of the honest
    literature-coverage count rather than letting a Pearson of exactly 1.000 read as a pass.
    """
    if not benchmark_id:
        return False
    return str(benchmark_id) in fitted_to_benchmark_lanes()


_MATRIX_CLASS_ANCHORS = (
    MatrixCalibrationAnchor(
        protein_type="pea_iso",
        target_class="aldehyde",
        observable_factor=1.0,
        evidence_strength="class_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline (generic aldehyde transfer)",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="pea_iso",
        target_class="furan",
        observable_factor=1.0,
        evidence_strength="class_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer)",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="pea_iso",
        target_class="sulfur",
        observable_factor=1.0,
        evidence_strength="directional_transferred",
        source="Interpolated base sulfur yield matching internal benchmark limits",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="pea_iso",
        target_class="pyrazine",
        observable_factor=1.0,
        evidence_strength="directional_transferred",
        source="Interpolated base pyrazine yield matching internal benchmark limits",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="soy_iso",
        target_class="aldehyde",
        observable_factor=2.209,
        evidence_strength="class_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio applied over aldehyde class",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="soy_iso",
        target_class="furan",
        observable_factor=5.92,
        evidence_strength="class_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio applied over furan class",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="soy_iso",
        target_class="sulfur",
        observable_factor=1.0,
        evidence_strength="directional_transferred",
        source="Interpolated base sulfur yield matching internal benchmark limits",
        fallback_mode="class_level",
    ),
    MatrixCalibrationAnchor(
        protein_type="soy_iso",
        target_class="pyrazine",
        observable_factor=1.0,
        evidence_strength="directional_transferred",
        source="Interpolated base pyrazine yield matching internal benchmark limits",
        fallback_mode="class_level",
    ),
    # --- Surrogate-only placeholders (Lane G, 2026-05-10b). No matrix-specific
    # calibration data exists for these protein types; class-level fallback only
    # so they are recognised as targets and the scope-check (Lane E) flags them
    # as out of scope instead of silently inheriting pea_iso/soy_iso factors.
    MatrixCalibrationAnchor(
        protein_type="wheat_gluten",
        target_class="aldehyde",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no wheat_gluten matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="wheat_gluten",
        target_class="furan",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no wheat_gluten matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="wheat_gluten",
        target_class="sulfur",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no wheat_gluten matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="wheat_gluten",
        target_class="pyrazine",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no wheat_gluten matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="myco",
        target_class="aldehyde",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no mycoprotein matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="myco",
        target_class="furan",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no mycoprotein matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="myco",
        target_class="sulfur",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no mycoprotein matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
    MatrixCalibrationAnchor(
        protein_type="myco",
        target_class="pyrazine",
        observable_factor=1.0,
        evidence_strength="surrogate_family",
        source="Placeholder: no mycoprotein matrix calibration data yet",
        fallback_mode="class_level",
        notes="No matrix-specific calibration data; flagged out of scope by scope-check.",
    ),
)

_MATRIX_RUNTIME_COMPOSITION_RULES = (
    MatrixRuntimeCompositionRule(
        protein_type="soy_iso",
        compound="hexanal",
        mode="compose_dynamic_retention",
        active_process_states=("intermediate_matrix", "heated_matrix", "aqueous_pre_extrusion_model", "extrusion_structured"),
        source="Ince 2024 reversible soy hexanal binding plus Xu 2023 thermal attenuation prior",
        notes="Ambient slurry remains frozen to preserve the historical Pratap-Singh benchmark calibration.",
    ),
)


def _load_persisted_matrix_multipliers() -> Dict[str, Dict[str, object]]:
    if not data_paths.MATRIX_CALIBRATION_OFFSETS.exists():
        return {}
    try:
        payload = json.loads(data_paths.MATRIX_CALIBRATION_OFFSETS.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    entries = payload.get("entries", {}) if isinstance(payload, dict) else {}
    if not isinstance(entries, dict):
        return {}
    normalized: Dict[str, Dict[str, object]] = {}
    for protein_type, row in entries.items():
        if not isinstance(row, dict):
            continue
        normalized[str(protein_type)] = dict(row)
    return normalized


def _runtime_matrix_multipliers() -> Dict[str, float]:
    raw = os.environ.get(_RUNTIME_MULTIPLIER_ENV)
    if not raw:
        return {}
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    if not isinstance(payload, dict):
        return {}
    multipliers: Dict[str, float] = {}
    for protein_type, value in payload.items():
        try:
            multipliers[str(protein_type)] = float(value)
        except (TypeError, ValueError):
            continue
    return multipliers


def get_matrix_observable_multiplier(protein_type: Optional[str]) -> Dict[str, object]:
    normalized = str(protein_type or "").strip()
    runtime = _runtime_matrix_multipliers()
    if normalized in runtime:
        return {
            "multiplier": float(runtime[normalized]),
            "source": "runtime_override",
        }

    persisted = _load_persisted_matrix_multipliers().get(normalized, {})
    if persisted:
        try:
            multiplier = float(persisted.get("observable_factor_multiplier", 1.0) or 1.0)
        except (TypeError, ValueError):
            multiplier = 1.0
        return {
            "multiplier": multiplier,
            "source": str(persisted.get("source", "matrix_calibration_offsets")),
            "provenance": str(persisted.get("provenance", "matrix_calibration_offsets")),
        }

    return {
        "multiplier": 1.0,
        "source": "baseline_registry",
    }


def _apply_observable_multiplier_to_record(record: MatrixCalibrationRecord) -> MatrixCalibrationRecord:
    multiplier_info = get_matrix_observable_multiplier(record.protein_type)
    multiplier = float(multiplier_info.get("multiplier", 1.0) or 1.0)
    if math.isclose(multiplier, 1.0, rel_tol=1e-9, abs_tol=1e-9):
        return record
    extra_note = f" Observable multiplier {multiplier:.3f} from {multiplier_info.get('source', 'unknown')}."
    return MatrixCalibrationRecord(
        protein_type=record.protein_type,
        process_state=record.process_state,
        compound=record.compound,
        observable_factor=float(record.observable_factor) * multiplier,
        evidence_strength=record.evidence_strength,
        source=record.source,
        fallback_mode=record.fallback_mode,
        notes=f"{record.notes}{extra_note}".strip(),
    )


def _normalize_compound(name: str) -> str:
    """Compound identity for record matching: the registry id, or the lower-cased spelling
    when the name is not in data/keys/compounds.yml (2026-09-02)."""
    key = compound_keys.resolve(str(name))
    return key.id if key else str(name).strip().lower()


def _process_state_fallback_order(requested_state: str) -> tuple[str, ...]:
    requested = str(requested_state or "ambient_slurry")
    if requested == "extrusion_structured":
        return (requested, "aqueous_pre_extrusion_model", "heated_matrix", "intermediate_matrix", "ambient_slurry")
    if requested == "aqueous_pre_extrusion_model":
        return (requested, "heated_matrix", "intermediate_matrix", "ambient_slurry")
    if requested == "heated_matrix":
        return (requested, "intermediate_matrix", "ambient_slurry")
    if requested == "intermediate_matrix":
        return (requested, "ambient_slurry")
    return (requested,)


def determine_matrix_process_state(*, temperature_celsius: float, time_minutes: float, water_activity: Optional[float] = None) -> str:
    if water_activity is not None:
        aw = float(water_activity)
        if temperature_celsius >= 160.0 and aw <= 0.45:
            return "extrusion_structured"
        if temperature_celsius >= 140.0 and aw <= 0.65:
            return "aqueous_pre_extrusion_model"

    if temperature_celsius <= 55.0 and time_minutes <= 30.0:
        return "ambient_slurry"
    if temperature_celsius >= 110.0 or time_minutes >= 90.0:
        return "heated_matrix"
    return "intermediate_matrix"


# (protein_type, process_state) pairs for which we hold compound-specific
# observable calibration anchors. Used by is_formulation_in_calibration_scope().
_CALIBRATED_PROTEIN_PROCESS_PAIRS = frozenset(
    (record.protein_type, record.process_state) for record in _MATRIX_CALIBRATION_RECORDS
)


@dataclass(frozen=True)
class ScopeAssessment:
    """Result of comparing a formulation against the calibration convex hull.

    Attributes
    ----------
    in_scope:
        True iff the (protein_type, process_state) pair is present in the
        compound-specific calibration anchor set AND the temperature/pH/time
        envelope falls within the calibrated range for that pair.
    reasons:
        Human-readable list of *why* the formulation was flagged out of scope.
        Empty when ``in_scope`` is True.
    nearest_calibrated:
        The closest calibrated (protein_type, process_state) pair, used by the
        report banner to suggest a comparable in-scope formulation.
    """

    in_scope: bool
    reasons: tuple[str, ...]
    nearest_calibrated: Dict[str, str]

    def to_dict(self) -> Dict[str, object]:
        return {
            "in_scope": bool(self.in_scope),
            "reasons": list(self.reasons),
            "nearest_calibrated": dict(self.nearest_calibrated),
        }


# Calibrated envelope per (protein_type, process_state). Derived from the
# benchmark + anchor coverage actually used by the matrix panel today; kept
# explicit (rather than inferred at runtime) so the scope check is auditable.
_CALIBRATED_ENVELOPES: Dict[tuple[str, str], Dict[str, tuple[float, float]]] = {
    ("pea_iso", "ambient_slurry"): {"temperature_c": (4.0, 55.0), "pH": (5.5, 8.5), "time_min": (0.0, 60.0)},
    ("pea_iso", "heated_matrix"): {"temperature_c": (90.0, 145.0), "pH": (5.5, 8.5), "time_min": (5.0, 240.0)},
    ("soy_iso", "ambient_slurry"): {"temperature_c": (4.0, 55.0), "pH": (5.5, 8.5), "time_min": (0.0, 60.0)},
    ("soy_iso", "heated_matrix"): {"temperature_c": (90.0, 145.0), "pH": (5.5, 8.5), "time_min": (5.0, 240.0)},
}


def is_formulation_in_calibration_scope(
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
    temperature_celsius: Optional[float] = None,
    time_minutes: Optional[float] = None,
    pH: Optional[float] = None,
) -> ScopeAssessment:
    """Return whether a formulation falls inside the calibration convex hull.

    The hull is intentionally conservative: only (protein_type, process_state)
    pairs with explicit `MatrixCalibrationRecord` anchors count as in-scope, and
    only when temperature/pH/time fall inside the envelope observed during
    calibration. Anything else is flagged so the report can downgrade tiers
    (Lane E, sprint 2026-05-10b).
    """
    reasons: list[str] = []
    nearest: Dict[str, str] = {
        "protein_type": str(protein_type or "unknown"),
        "process_state": str(process_state or "unknown"),
    }

    if not protein_type:
        return ScopeAssessment(in_scope=False, reasons=("protein_type is missing",), nearest_calibrated=nearest)
    if not process_state:
        return ScopeAssessment(in_scope=False, reasons=("process_state is missing",), nearest_calibrated=nearest)

    pair = (str(protein_type), str(process_state))
    if pair not in _CALIBRATED_PROTEIN_PROCESS_PAIRS:
        # Suggest the nearest calibrated process_state for the same protein type.
        same_protein = [ps for (pt, ps) in _CALIBRATED_PROTEIN_PROCESS_PAIRS if pt == protein_type]
        if same_protein:
            nearest = {"protein_type": str(protein_type), "process_state": same_protein[0]}
            reasons.append(
                f"No calibration anchor for ({protein_type}, {process_state}); nearest calibrated process_state is {same_protein[0]}"
            )
        else:
            reasons.append(
                f"No calibration anchor for protein_type '{protein_type}' (calibrated set: pea_iso, soy_iso)"
            )
        return ScopeAssessment(in_scope=False, reasons=tuple(reasons), nearest_calibrated=nearest)

    envelope = _CALIBRATED_ENVELOPES.get(pair, {})
    if temperature_celsius is not None and "temperature_c" in envelope:
        lo, hi = envelope["temperature_c"]
        if not (lo <= float(temperature_celsius) <= hi):
            reasons.append(f"temperature {temperature_celsius:.1f} °C outside calibrated range [{lo:.0f}, {hi:.0f}]")
    if time_minutes is not None and "time_min" in envelope:
        lo, hi = envelope["time_min"]
        if not (lo <= float(time_minutes) <= hi):
            reasons.append(f"time {time_minutes:.1f} min outside calibrated range [{lo:.0f}, {hi:.0f}]")
    if pH is not None and "pH" in envelope:
        lo, hi = envelope["pH"]
        if not (lo <= float(pH) <= hi):
            reasons.append(f"pH {pH:.2f} outside calibrated range [{lo:.1f}, {hi:.1f}]")

    return ScopeAssessment(in_scope=not reasons, reasons=tuple(reasons), nearest_calibrated=nearest)


def get_matrix_calibration_record(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Optional[MatrixCalibrationRecord]:
    if not protein_type:
        return None
    normalized = _normalize_compound(compound)
    requested_state = process_state or "ambient_slurry"

    for candidate_state in _process_state_fallback_order(requested_state):
        for record in _MATRIX_CALIBRATION_RECORDS:
            if record.protein_type != protein_type:
                continue
            if record.process_state != candidate_state:
                continue
            if _normalize_compound(record.compound) != normalized:
                continue
            if candidate_state == requested_state:
                return _apply_observable_multiplier_to_record(record)
            return _apply_observable_multiplier_to_record(MatrixCalibrationRecord(
                protein_type=record.protein_type,
                process_state=requested_state,
                compound=record.compound,
                observable_factor=record.observable_factor,
                evidence_strength="process_state_mismatch",
                source=record.source,
                fallback_mode="nearest_process_state",
                notes=f"Requested process state '{requested_state}' falls back to nearest calibrated state '{candidate_state}'.",
            ))
    return None


def get_matrix_runtime_composition_policy(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Dict[str, str]:
    if not protein_type:
        return {
            "mode": "static_observable_calibration",
            "source": "none",
            "notes": "No protein-type-specific runtime composition policy is registered.",
        }

    normalized = _normalize_compound(compound)
    requested_state = process_state or "ambient_slurry"
    for rule in _MATRIX_RUNTIME_COMPOSITION_RULES:
        if rule.protein_type != protein_type:
            continue
        if _normalize_compound(rule.compound) != normalized:
            continue
        if requested_state in rule.active_process_states:
            return {
                "mode": rule.mode,
                "source": rule.source,
                "notes": rule.notes,
            }

    return {
        "mode": "static_observable_calibration",
        "source": "historical_calibration_default",
        "notes": "Observable calibration is used as-is for this compound/process-state pair.",
    }


def describe_matrix_calibration(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Dict[str, object]:
    multiplier_info = get_matrix_observable_multiplier(protein_type)
    multiplier = float(multiplier_info.get("multiplier", 1.0) or 1.0)
    record = get_matrix_calibration_record(
        compound,
        protein_type=protein_type,
        process_state=process_state,
    )
    if record is None:
        from src.matrix_correction import classify_volatile_matrix_family
        target_class = classify_volatile_matrix_family(compound)
        for anchor in _MATRIX_CLASS_ANCHORS:
            if anchor.protein_type == protein_type and anchor.target_class == target_class:
                return {
                    "calibration_source": anchor.source,
                    "calibration_process_state": process_state or "unknown",
                    "calibration_evidence_strength": anchor.evidence_strength,
                    "calibration_fallback_mode": anchor.fallback_mode,
                    "calibration_observable_factor": float(anchor.observable_factor) * multiplier,
                    "calibration_observable_multiplier": multiplier,
                    "calibration_observable_multiplier_source": multiplier_info.get("source", "baseline_registry"),
                    "calibration_notes": anchor.notes,
                }
        return {
            "calibration_source": "class_fallback",
            "calibration_process_state": process_state or "unknown",
            "calibration_evidence_strength": "heuristic",
            "calibration_fallback_mode": "class_level",
            "calibration_observable_factor": None,
            "calibration_observable_multiplier": multiplier,
            "calibration_observable_multiplier_source": multiplier_info.get("source", "baseline_registry"),
            "calibration_notes": "No compound-specific matrix calibration is registered for this compound/process state.",
        }
    return {
        "calibration_source": record.source,
        "calibration_process_state": record.process_state,
        "calibration_evidence_strength": record.evidence_strength,
        "calibration_fallback_mode": record.fallback_mode,
        "calibration_observable_factor": float(record.observable_factor),
        "calibration_observable_multiplier": multiplier,
        "calibration_observable_multiplier_source": multiplier_info.get("source", "baseline_registry"),
        "calibration_notes": record.notes,
    }
