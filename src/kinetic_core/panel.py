"""
src/kinetic_core/panel.py -- THE ONE bundle -> core-spec mapping (retirement step B2).

Before B2 the mapping from a frozen benchmark bundle to a ``FormulationSpec``
lived inside ``scripts/generators/generate_cutover_final_exam.py`` (deleted at B5b). The core
Monte-Carlo envelope needs exactly the same mapping, and two copies of a
schema reader are how the matrix_path bundles came to be scored with
``measured = None`` for a whole wave (see :func:`measured_value`). So the
mapping is lifted HERE, and both the exam and the envelope import it.

Nothing in this module integrates anything or reads a fit report. It answers
three questions about a bundle -- what is charged, under what process, and
what was measured -- and one question about the panel: which bundles are on
it, and under which tag.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from src import data_paths

# ---------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------

#: The ratio units a bundle's ``holdout_targets`` may declare, and the factor
#: that turns a mole fraction (mol product / mol limiting precursor) into them.
RATIO_UNIT_FACTORS: Mapping[str, float] = {
    "mol_percent": 100.0,
    "umol_per_mol_limiting_precursor": 1.0e6,
}


# ---------------------------------------------------------------------------
# Panel membership
# ---------------------------------------------------------------------------

PANEL_TRUST_LOOP = "trust_loop"
PANEL_MAILLARD_PATH_HOLDOUT = "maillard_path_holdout"
PANEL_EXTERNAL_MATRIX = "external_matrix"

#: Bundles that are LEGACY-MODEL OUTPUT, not measurements. They stay on disk
#: (tests reference them) but leave the scored panel; the plan is to regenerate
#: them later as core frozen predictions.
SYNTHETIC_SNAPSHOT_SUFFIX = "_Internal2026"

#: ``evidence_class`` values that keep a bundle OFF the scored panel.
UNSCORED_EVIDENCE_CLASSES: Tuple[str, ...] = ("diagnostic_only",)


#: WAVE B2.4 -- THE FOUR SHARED ROWS, DECLARED.
#:
#: D1's reconciliation sec. 5 established that four of the cutover exam's
#: answered points are not merely analogous to four rows of the B2.1/B2.2/B2.3
#: hold-out panel -- they are THE SAME MEASUREMENTS. Hofmann & Schieberle
#: 1998's ribose + cysteine pots at 145 C / 20 min, pH 3 and pH 7, scored once
#: in the exam in fold error and once in the panel against a 3x band. The
#: panel's `hofmann_ribose_pH7_FFT` reads 0.00200x and the exam reads 498.99x:
#: the same number, inverted.
#:
#: Consequence, and the reason this table exists: THE EXAM AND THE PANEL ARE
#: NOT INDEPENDENT EVIDENCE. Anyone who reads "the exam says X and the panel
#: agrees" is, on these four points, reading one measurement twice. It is
#: declared here, in the exam artifact, in the panel artifact and in the core
#: envelope artifact so the double-counting is visible everywhere it matters.
SHARED_WITH_HOLDOUT_PANEL: Dict[Tuple[str, str], str] = {
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3",
     "2-Methyl-3-furanthiol (MFT)"): "hofmann_ribose_pH3_MFT",
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3",
     "2-Furfurylthiol (FFT)"): "hofmann_ribose_pH3_FFT",
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7",
     "2-Methyl-3-furanthiol (MFT)"): "hofmann_ribose_pH7_MFT",
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7",
     "2-Furfurylthiol (FFT)"): "hofmann_ribose_pH7_FFT",
}


def load_bundle(path: Path | str) -> Dict[str, Any]:
    """Read one bundle. Plain JSON; no schema massaging happens here."""
    return json.loads(Path(path).read_text())


def is_scored_bundle(bench: Mapping[str, Any], path: Path | str) -> Tuple[bool, str]:
    """
    Whether a bundle belongs on the SCORED panel, and if not, why.

    The two ``*_Internal2026`` snapshots are legacy-model output and any bundle
    self-labelled ``evidence_class: diagnostic_only`` is a diagnostic, not a
    measurement. Both are recorded as skipped rather than silently dropped.
    """
    stem = Path(path).stem
    if stem.endswith(SYNTHETIC_SNAPSHOT_SUFFIX):
        return False, (
            f"synthetic snapshot ({SYNTHETIC_SNAPSHOT_SUFFIX}): legacy-model "
            "output, not a measurement"
        )
    evidence_class = str(bench.get("evidence_class") or "")
    if evidence_class in UNSCORED_EVIDENCE_CLASSES:
        return False, f"evidence_class: {evidence_class}"
    return True, ""


def panel_bundles(
    *,
    benchmarks_dir: Path = data_paths.BENCHMARKS_DIR,
    maillard_path_dir: Path = data_paths.MAILLARD_PATH_HOLDOUT_DIR,
    external_dir: Path = data_paths.EXTERNAL_VALIDATION_DIR,
) -> Tuple[List[Tuple[Path, str]], List[Dict[str, str]]]:
    """
    The union panel: ``[(path, panel_tag), ...]`` plus the skipped bundles.

    * ``data/benchmarks/*.json`` (NON-recursive; the quarantined and
      step-level-unreachable subdirectories are physically separated and
      stay out)                                    -> ``trust_loop``
    * ``data/benchmarks/external_validation/maillard_path/*.json`` -> ``maillard_path_holdout``
    * ``data/benchmarks/external_validation/*.json``               -> ``external_matrix``
    """
    groups = (
        (benchmarks_dir, PANEL_TRUST_LOOP),
        (maillard_path_dir, PANEL_MAILLARD_PATH_HOLDOUT),
        (external_dir, PANEL_EXTERNAL_MATRIX),
    )
    scored: List[Tuple[Path, str]] = []
    skipped: List[Dict[str, str]] = []
    for directory, tag in groups:
        if not directory.exists():
            continue
        for path in sorted(directory.glob("*.json")):
            bench = load_bundle(path)
            keep, reason = is_scored_bundle(bench, path)
            if keep:
                scored.append((path, tag))
            else:
                skipped.append(
                    {
                        "benchmark_id": str(bench.get("benchmark_id") or path.stem),
                        "bench_file": data_paths.rel(path),
                        "panel": tag,
                        "reason": reason,
                    }
                )
    return scored, skipped


# ---------------------------------------------------------------------------
# Reading a bundle
# ---------------------------------------------------------------------------


def bundle_targets(bench: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    The compounds a bundle asks to be scored on, ``{compound: target_spec}``.

    ``holdout_targets`` (the maillard_path schema) wins over
    ``measured_volatiles`` (the trust-loop and matrix schemas). A bundle with
    neither returns an empty dict -- and is then reported as carrying no
    target, never as having been answered.
    """
    targets = bench.get("holdout_targets") or bench.get("measured_volatiles") or {}
    return {
        str(compound): (spec if isinstance(spec, dict) else {})
        for compound, spec in targets.items()
    }


def limiting_precursor_molar(
    bench: Mapping[str, Any],
) -> Tuple[Optional[str], Optional[float]]:
    """The smallest non-zero precursor charge, as ``(name, mol/L)``."""
    items = [
        (str(name), float(data.get("concentration_mM", 0.0)) / 1000.0)
        for name, data in (bench.get("precursors") or {}).items()
        if float(data.get("concentration_mM", 0.0)) > 0.0
    ]
    if not items:
        return None, None
    return min(items, key=lambda kv: kv[1])


def measured_value(
    bench: Mapping[str, Any], compound: str, target_spec: Mapping[str, Any]
) -> Optional[float]:
    """
    The measured value for one point, from wherever this bundle's schema keeps it.

    WAVE B6 FIX -- A PRE-EXISTING DEFECT, FOUND BY THE LIPID LANE.
    ------------------------------------------------------------
    The exam's ``score_core`` and ``score_old_lane`` both took ``target_value``
    and then fell back to ``reference_volatiles[compound]["conc_ppb"]``. The
    seventeen maillard_path bundles carry ``holdout_targets`` with
    ``target_value``, so they scored. The FOUR matrix_path bundles do not:
    they carry the value in ``measured_volatiles[compound]["conc_ppb"]``, and
    ``reference_volatiles`` does not exist in them at all. Every matrix_path
    point therefore came back ``measured = None`` and was reported with NO
    fold error -- in BOTH lanes, since the two functions shared the bug.

    That went unnoticed because the core REFUSED all four bundles before B6,
    so a missing measured value cost nothing visible. B6 answers them, and a
    prediction with no referee is exactly what the exam exists to prevent.

    This is a schema-reading fix, defensible with no reference to any
    hold-out value -- the same standing as Amendment 9 clause 2's buffer-field
    completion. It changed the OLD lane's reported numbers too, which is the
    honest consequence of a shared bug: those points were never scored, and
    now they are.
    """
    for candidate in (
        target_spec.get("target_value"),
        (bench.get("reference_volatiles") or {}).get(compound, {}).get("conc_ppb"),
        target_spec.get("conc_ppb"),
        (bench.get("measured_volatiles") or {}).get(compound, {}).get("conc_ppb"),
    ):
        if candidate is not None:
            return float(candidate)
    return None


# ---------------------------------------------------------------------------
# B2.3: THE BUFFER FIELD, AND WHY THE EXAM IS REPORTED BOTH WAYS
# ---------------------------------------------------------------------------
# `results/validation/kinetic_core_b2_2_diagnosis.md` sec. 4 found that the
# frozen bundles record NO buffer, so bundles literally named
# `..._ribose_cysteine_buffer_...` were being integrated as WATER -- and sized
# that schema gap at an 11x swing in predicted 2-furfurylthiol on exactly the
# system shape that is the worst point in the exam report.
#
# `docs/reference/FIT_HOLDOUT_DECLARATION.md` Amendment 9 clause 2 completes
# the buffer field AND requires the exam to be reported BOTH WAYS --
# buffer-completed and as-was -- in one artifact, PERMANENTLY. The reason the
# as-was column is permanent rather than transitional: every earlier number
# this repo has published was computed as-was, and a report that silently
# replaces them makes its own history unreadable. Two columns cost one extra
# integration per row and settle the question of how much of the Yiltirak
# failure was chemistry and how much was a data-schema gap.
#
# NOTE ON REACH: only the SULFUR lane carries a pH state, so the buffer can
# only move a sulfur-lane row. The acrylamide-lane bundles (Chang, Lin, Ye,
# Schibilsky) are IDENTICAL in both columns and that identity is a reported
# result, not an omission -- it is the same gap X-7 names.


def buffer_from_bundle(bench: Mapping[str, Any]):
    """
    Translate a bundle's completed ``conditions.buffer`` block into a
    ``BufferSpec``. Returns ``None`` when the bundle has no block at all.

    ``buffer_unknown`` maps to the engine's DECLARED DEFAULT -- unbuffered,
    ``declared=False``, which raises the extrapolation warning. That is the
    correct handling and it is the point of recording unknown rather than
    guessing: an unknown buffer produces a WARNED extrapolation, not a silent
    assumption in either direction.
    """
    from src.kinetic_core.ph_state import BufferSpec

    block = (bench.get("conditions") or {}).get("buffer")
    if not isinstance(block, dict):
        return None
    species = str(block.get("species", "buffer_unknown"))
    molarity = block.get("concentration_M")
    note = str(block.get("provenance_note", ""))[:200]
    if species == "buffer_unknown":
        return BufferSpec(
            kind="none", declared=False,
            source=f"BUFFER UNKNOWN, recorded as such rather than guessed: {note}")
    if species == "none":
        return BufferSpec(
            kind="none", declared=True,
            source=f"the source states NO BUFFER: {note}")
    if species in ("phosphate", "potassium_phosphate"):
        return BufferSpec(
            kind="phosphate", phosphate_mol_l=float(molarity), declared=True,
            source=note)
    if species == "acetate":
        return BufferSpec(
            kind="acetate", acetate_mol_l=float(molarity), declared=True,
            source=note)
    if species == "citrate_phosphate":
        return BufferSpec(
            kind="citrate_phosphate", phosphate_mol_l=float(molarity) * 2.0 / 3.0,
            citrate_mol_l=float(molarity) / 3.0, declared=True, source=note)
    raise SystemExit(
        f"{bench.get('benchmark_id')}: buffer species {species!r} has no "
        f"BufferSpec translation. Add one deliberately -- do not fall through "
        f"to water, which is the defect this whole block exists to remove.")


def core_spec(bench: Mapping[str, Any], *, use_buffer: bool = True):
    """
    A bundle as a ``FormulationSpec``: precursors in mM, an isothermal hold,
    the declared pH and water activity, the protein type as the matrix
    descriptor, and -- unless ``use_buffer=False`` -- the completed buffer.

    ``use_buffer=False`` reproduces every pre-B2.3 number exactly: no buffer
    supplied means the engine's declared default, which is unbuffered plus an
    extrapolation warning.
    """
    from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram

    conditions = bench.get("conditions") or {}
    return FormulationSpec(
        name=str(bench.get("benchmark_id")),
        precursors={
            str(name): float(data.get("concentration_mM", 0.0))
            for name, data in (bench.get("precursors") or {}).items()
        },
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(
                float(conditions.get("temp_C", 0.0)),
                float(conditions.get("time_min", 0.0)),
            ),
            ph=float(conditions.get("ph", 7.0)),
            water_activity=(
                float(conditions["water_activity"])
                if conditions.get("water_activity") is not None
                else None
            ),
            matrix=str(bench.get("protein_type") or "water"),
            buffer=buffer_from_bundle(bench) if use_buffer else None,
        ),
    )


def core_native_value(
    run, compound: str, unit: str, limiting_molar: Optional[float]
) -> Optional[float]:
    """
    The core's prediction in the TARGET's own unit.

    ppb is direct (the core reports ug/L and every scored bundle is aqueous at
    ~1 kg/L, which is the bundles' own stated basis). The ratio units are
    computed from the core's MOLAR state, so no molecular weight and no assumed
    concentration basis enters -- the 342/200 lesson.
    """
    if unit == "ppb":
        return run.concentrations_ug_per_l.get(compound)
    if unit in RATIO_UNIT_FACTORS:
        key = run.declaration.mapped_targets.get(compound)
        if key is None or not limiting_molar:
            return None
        mmol_per_l = float(run.species_mmol_per_l.get(key, 0.0))
        mol_per_l = mmol_per_l / 1000.0
        return RATIO_UNIT_FACTORS[unit] * (mol_per_l / limiting_molar)
    return None


# ---------------------------------------------------------------------------
# How the measurement was made -- decides which observable bands apply
# ---------------------------------------------------------------------------

#: The bundle quantified the analyte in the HEADSPACE (HS-SPME, static
#: headspace): the core's declared K_aw band and the HS-SPME same-sample
#: dispersion are measurement facts of THIS number.
QUANTIFICATION_HEADSPACE = "headspace"
#: The bundle quantified the analyte in the LIQUID/EXTRACT (SIDA after solvent
#: extraction, HPLC-UV, LC-MS/MS, internal-standard GC/MS of an extract): no
#: air/water partition step and no SPME fibre enters the number, so neither
#: headspace band applies to it.
QUANTIFICATION_EXTRACTION = "extraction"
#: The bundle declares no ``quantification_class``. The envelope then keeps
#: the core's own convention (``matrix_oav.absolute_concentration``: every
#: absolute ppb carries the HS-SPME dispersion) and says so on the row.
QUANTIFICATION_UNDECLARED = "undeclared"

_HEADSPACE_MARKERS = ("spme", "headspace", "hs-gc", "hs_gc", "dhs", "dynamic_headspace")
_EXTRACTION_MARKERS = (
    "isotope_dilution", "sida", "hplc", "lcms", "lc-ms", "internal_standard",
    "external_standard", "response_factor", "calibration_curve",
)


def _first_value(node: Any, key: str) -> Optional[str]:
    """Depth-first first occurrence of ``key`` in a nested bundle."""
    if isinstance(node, Mapping):
        if key in node and isinstance(node[key], str):
            return node[key]
        for value in node.values():
            found = _first_value(value, key)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = _first_value(value, key)
            if found is not None:
                return found
    return None


def quantification_family(bench: Mapping[str, Any]) -> Tuple[str, str]:
    """
    ``(family, why)`` for a bundle: one of :data:`QUANTIFICATION_HEADSPACE`,
    :data:`QUANTIFICATION_EXTRACTION`, :data:`QUANTIFICATION_UNDECLARED`.

    Read from the bundle's own ``quantification_class`` (written during the
    primary-source verification waves under ``content_verification``); the
    headspace markers are checked first because a headspace class can also
    name its calibration ("... internal standard ..., SPME-GC-MS").
    """
    declared = _first_value(bench, "quantification_class")
    if not declared:
        return QUANTIFICATION_UNDECLARED, "no quantification_class in the bundle"
    lowered = declared.lower()
    if any(marker in lowered for marker in _HEADSPACE_MARKERS):
        return QUANTIFICATION_HEADSPACE, f"quantification_class={declared[:80]!r}"
    if any(marker in lowered for marker in _EXTRACTION_MARKERS):
        return QUANTIFICATION_EXTRACTION, f"quantification_class={declared[:80]!r}"
    return QUANTIFICATION_UNDECLARED, f"quantification_class={declared[:80]!r} matches no known family"


__all__ = [
    "PANEL_EXTERNAL_MATRIX",
    "PANEL_MAILLARD_PATH_HOLDOUT",
    "PANEL_TRUST_LOOP",
    "QUANTIFICATION_EXTRACTION",
    "QUANTIFICATION_HEADSPACE",
    "QUANTIFICATION_UNDECLARED",
    "RATIO_UNIT_FACTORS",
    "SHARED_WITH_HOLDOUT_PANEL",
    "SYNTHETIC_SNAPSHOT_SUFFIX",
    "UNSCORED_EVIDENCE_CLASSES",
    "buffer_from_bundle",
    "bundle_targets",
    "core_native_value",
    "core_spec",
    "is_scored_bundle",
    "limiting_precursor_molar",
    "load_bundle",
    "measured_value",
    "panel_bundles",
    "quantification_family",
]
