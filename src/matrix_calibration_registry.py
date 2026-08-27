from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional


ROOT = Path(__file__).resolve().parents[1]
MATRIX_CALIBRATION_OFFSETS_PATH = ROOT / "data" / "lit" / "matrix_calibration_offsets.json"
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
    # These factors are 1.0 BY CONSTRUCTION: this lane is the reference against which the
    # others are expressed, and the yields it multiplies
    # (benchmark_validation.MATRIX_BENCHMARK_BASE_MARKER_YIELDS) are
    # pea_isolate_40C_PratapSingh2021's own measured ppb divided by one common scale of
    # 1268-1271. So `pea_isolate_40C_PratapSingh2021` scoring Pearson 1.000 / max_ratio 1.002
    # is the lane reproducing the numbers it was built from, not a prediction. The factor is
    # still 1.0 -- nothing changes at runtime -- but the label now says what it is.
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        observable_factor=1.0,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
        notes="Reference compound for the pea matrix-only intake/headspace lane. Factor is 1.0 by construction; the underlying base marker yield is this benchmark's own measurement rescaled, so this lane is fit-recovery, not a prediction.",
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
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
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
        observable_factor=0.228776,
        previous_value=0.22877612093571738,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
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
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
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
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT nonanal anchor. BACK-SOLVED from this benchmark's own measured 24 ppb.",
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
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        observable_factor=0.453 / 0.205,
        evidence_strength=FITTED_TO_BENCHMARK,
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
        notes="Soy release ratio carried relative to the pea reference intake lane. Both numerator and denominator are rescaled measurements from the two benchmarks this lane is scored against -- fit-recovery, not a prediction.",
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
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
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
        observable_factor=(0.453 / 0.205) * (1.0 - 0.7060),
        evidence_strength="conditional_literature_anchored",
        source="Shu 2024 heated soy off-flavour attenuation carried onto the Pratap-Singh soy ambient baseline",
        fallback_mode="compound_specific_process_state",
        notes="High-severity soy treatment prior for heated matrix states. Useful for reliability and directional accuracy, but not a meaty benchmark anchor. 2026-08-27 (Wave I): kept as conditional_literature_anchored -- the ATTENUATION (1 - 0.7060) is a real Shu 2024 literature figure and no panel benchmark exists for heated soy, so this lane is not fit-recovery. But note its BASELINE factor (0.453/0.205) is fit-recovery (see the soy ambient block), so the anchoring is only as strong as one literature attenuation applied to a back-solved base.",
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
    if not MATRIX_CALIBRATION_OFFSETS_PATH.exists():
        return {}
    try:
        payload = json.loads(MATRIX_CALIBRATION_OFFSETS_PATH.read_text(encoding="utf-8"))
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
    return str(name).strip().lower()


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
