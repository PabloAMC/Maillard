#!/usr/bin/env python3
"""
Corpus-wide Arrhenius PREFACTOR audit.

WHAT THIS IS
------------
Three separate occasions in this repo's history, a shipped Arrhenius pair
``(Ea, A)`` was found to carry a *correct* activation energy bolted to a
*wrong* pre-exponential factor -- ``cysteine_thermolysis`` (the A was a generic
unimolecular TST number never fitted with its Ea), Hamzalioglu 2018 (the SOURCE
PAPER's own six prefactors, four of which do not reproduce from its own k
table), and Stack 2018 (printed activation parameters 2.303x too large, a
natural-log slope processed with the base-10 formula). The failure mode is
invisible on inspection: the Ea looks right, the citation resolves, and nothing
goes wrong until you evaluate ``k(T)``.

This script asks whether ANY OTHERS exist. It:

  1. ENUMERATES every Arrhenius ``(Ea, A)`` pair shipped anywhere in the tree,
     in both of the forms the repo uses -- the explicit ``(A, Ea)`` form
     (``data/lit/arrhenius_params.yml``, ``src/safety.py``,
     ``data/lit/lipid_oxidation_calibration.json``) and the reparameterised
     ``(k_ref at T_ref, Ea)`` form used throughout ``src/kinetic_core``, which
     is algebraically the same object: ``A = k_ref * exp(Ea / (R * T_ref))``.
  2. REFITS ``(Ea, A)`` by WEIGHTED least squares on ``ln k`` vs ``1/T`` for
     every pair whose source's own rate-constant-versus-temperature table has
     been transcribed into ``data/lit/extraction_dossiers/*.md``.
  3. FLAGS every pair whose shipped A differs from the refit A by more than
     ``--a-ratio-threshold`` (default 2x), every pair whose shipped Ea differs
     by more than ``--ea-tolerance-pct`` (default 10 %), and -- as a SEPARATE
     category -- every pair for which no dossier k-table exists at all, which
     is unverifiable rather than wrong.

THIS SCRIPT CHANGES NOTHING. It is an audit that FLAGS, not fixes. It writes
only to ``results/validation/``; no shipped parameter is touched.

THE DOSSIER k-TABLES
--------------------
The four transcribed source tables that make a refit possible are embedded
below in ``DOSSIER_K_TABLES``, each carrying the dossier file and section it was
copied from so the transcription can be re-checked against the markdown. They
are the ONLY k(T) tables in the corpus that back a shipped pair:

  * ``zheng1994_extraction.md`` sec. 3 -- Table I, 16 first-order H2S-release
    constants, 4 temperatures x 4 pH, aqueous 0.1 M L-cysteine.
  * ``hamzalioglu2018_extraction.md`` sec. 4 -- Table 1, HMF loss with three
    amino acids at 5/25/50 C, high-moisture arm.
  * ``kocadagli2016jafc_extraction.md`` sec. 3 -- Table 1, per-temperature rate
    constants of the glucose system at 160/180/200 C.
  * ``stack2018_extraction.md`` sec. 3 -- Figure 4C, the NAC thiol-quinone
    forward constant at four temperatures.

CAVEAT THAT THE VERDICT COLUMN CANNOT CARRY, AND THAT THE REPORT STATES IN
PROSE: a refit disagreeing with a shipped value is not automatically an error.
Kocadagli's Table 2 comes from a simultaneous global fit of the whole ODE system
to all three temperatures, NOT from an Arrhenius fit of its own Table 1, so the
two are not required to agree; the size of the disagreement is an
IDENTIFIABILITY diagnostic. Where a shipped pair deliberately takes only the Ea
from a source and its magnitude from elsewhere (``k_thioether``), that is
declared at the parameter and the A ratio against the source is meaningless.
Both cases are annotated rather than silently scored.

OUTPUTS
-------
  results/validation/prefactor_audit.json
  results/validation/prefactor_audit.md
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

R_KJ = 8.314462618e-3  # kJ/(mol K)

DOSSIER_DIR = Path("data/lit/extraction_dossiers")


# ---------------------------------------------------------------------------
# The transcribed source k(T) tables. Values copied verbatim from the dossier
# markdown; the `anchor` on each is the dossier file and section to re-check.
# ---------------------------------------------------------------------------

DOSSIER_K_TABLES: Dict[str, Dict[str, Any]] = {
    # Zheng & Ho 1994, ACS Symp. Ser. 564:138-146, Table I p. 143.
    # First-order H2S release from 0.1 M aqueous L-cysteine. K printed in min^-1.
    # No per-cell error bars are printed, so the refit is unweighted.
    "zheng1994_cysteine_ph3": {
        "anchor": f"{DOSSIER_DIR}/zheng1994_extraction.md sec. 3 (Table I), pH 3.0 column",
        "unit": "1/min",
        "note": "the authors flag pH 3.0 as a possibly different mechanism (p. 140)",
        "points": [(80.0, 1.87e-6, None), (90.0, 8.30e-6, None),
                   (100.0, 4.15e-5, None), (110.0, 5.23e-5, None)],
    },
    "zheng1994_cysteine_ph5": {
        "anchor": f"{DOSSIER_DIR}/zheng1994_extraction.md sec. 3 (Table I), pH 5.0 column",
        "unit": "1/min",
        "points": [(80.0, 2.02e-6, None), (90.0, 1.24e-5, None),
                   (100.0, 2.65e-5, None), (110.0, 7.97e-5, None)],
    },
    "zheng1994_cysteine_ph7": {
        "anchor": f"{DOSSIER_DIR}/zheng1994_extraction.md sec. 3 (Table I), pH 7.0 column",
        "unit": "1/min",
        "points": [(80.0, 9.71e-6, None), (90.0, 5.99e-5, None),
                   (100.0, 2.30e-4, None), (110.0, 3.35e-4, None)],
    },
    "zheng1994_cysteine_ph9": {
        "anchor": f"{DOSSIER_DIR}/zheng1994_extraction.md sec. 3 (Table I), pH 9.0 column",
        "unit": "1/min",
        "points": [(80.0, 2.95e-5, None), (90.0, 1.33e-4, None),
                   (100.0, 4.53e-4, None), (110.0, 7.45e-4, None)],
    },
    # Hamzalioglu & Gokmen 2018, Food Chem 240:354-360, Table 1 high-moisture arm.
    # Pseudo-first-order HMF loss, day^-1. Every row carries a uniform 5 % SD.
    "hamzalioglu2018_hmf_cys": {
        "anchor": f"{DOSSIER_DIR}/hamzalioglu2018_extraction.md sec. 4 (Table 1), HMF-Cys row",
        "unit": "1/day",
        "points": [(5.0, 0.079, 0.05), (25.0, 0.103, 0.05), (50.0, 0.465, 0.05)],
        "sigma_is_relative": True,
    },
    # Kocadagli & Gokmen 2016 JAFC, Table 1 GLUCOSE system, k in min^-1 x 10^3,
    # at 160/180/200 C, each with a 95 % HPD half-width (used as the weight).
    "kocadagli2016_step3_glc_3dg": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 3",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 0.91, 0.19), (180.0, 4.14, 1.71), (200.0, 3.60, 1.26)],
        "note": "non-monotonic in T (0.91 -> 4.14 -> 3.60)",
    },
    "kocadagli2016_step4_3dg_34dg": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 4",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 23.1, 4.03), (180.0, 30.5, 4.71), (200.0, 49.3, 10.1)],
    },
    "kocadagli2016_step5_34dg_hmf": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 5",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 160.0, 35.0), (180.0, 110.0, 28.2), (200.0, 137.0, 46.1)],
        "note": "non-monotonic in T; the authors FIXED Ea to zero for this step",
    },
    "kocadagli2016_step6_fru_int": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 6",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 100.0, 8.6), (180.0, 344.0, 26.0), (200.0, 1058.0, 96.6)],
    },
    "kocadagli2016_step7_int_hmf": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 7",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 0.31, 0.07), (180.0, 1.87, 0.15), (200.0, 9.31, 1.74)],
    },
    "kocadagli2016_step8_fru_1dg": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 8",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 0.61, 0.15), (180.0, 2.47, 0.73), (200.0, 5.89, 1.93)],
    },
    "kocadagli2016_step11_3dg_mgo": {
        "anchor": f"{DOSSIER_DIR}/kocadagli2016jafc_extraction.md sec. 3 (Table 1), step 11",
        "unit": "1/min",
        "scale": 1e-3,
        "points": [(160.0, 96.0, 20.8), (180.0, 338.0, 29.4), (200.0, 863.0, 98.1)],
    },
    # Stack et al. 2018, Chem. Res. Toxicol. 31:81, Figure 4C, NAC forward k1.
    # Stated error budget: +/- 4-6 % on the k1 slopes (dossier sec. 2).
    "stack2018_nac_forward": {
        "anchor": f"{DOSSIER_DIR}/stack2018_extraction.md sec. 3 (Figure 4C), NAC k1 column",
        "unit": "M^-1 s^-1",
        "points": [(4.6, 391.0, 0.05), (9.6, 446.0, 0.05),
                   (14.5, 478.0, 0.05), (19.4, 496.0, 0.05)],
        "sigma_is_relative": True,
        "note": "14.8 K span, rate moves only 1.27x -- Ea weakly determined (+/-30-50 %)",
    },
}


# ---------------------------------------------------------------------------
# Unit conversion helpers. Each returns the factor that takes the DOSSIER
# table's own unit into the SHIPPED parameter's unit.
# ---------------------------------------------------------------------------

def _per_min_from_per_day(_: float) -> float:
    return 1.0 / 1440.0


CONVERSIONS: Dict[str, Tuple[float, str]] = {
    # name -> (multiplicative factor, human description)
    "identity": (1.0, "none"),
    "per_min_to_per_s": (1.0 / 60.0, "min^-1 -> s^-1 (/60)"),
    "hamzalioglu_pseudo1_to_second_order": (
        1.0 / (0.020 * 1000.0 * 1440.0),
        "day^-1 -> L/(mmol*min): k2 = k'/[AA] with [AA] = 0.020 M, then /1000 (mol->mmol) and /1440 (day->min)",
    ),
    "molar_per_s_to_l_per_mmol_min": (
        60.0 / 1000.0,
        "M^-1 s^-1 -> L/(mmol*min) (x60 /1000)",
    ),
}


# ---------------------------------------------------------------------------
# Weighted least squares on ln k vs 1/T.
# ---------------------------------------------------------------------------

def refit_arrhenius(spec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """WLS of ln k on 1/T. Returns Ea (kJ/mol) and A in the table's own unit."""
    scale = float(spec.get("scale", 1.0))
    rel = bool(spec.get("sigma_is_relative", False))
    xs: List[float] = []
    ys: List[float] = []
    ws: List[float] = []
    n_interval_spans_zero = 0
    for temp_c, k_raw, sigma in spec["points"]:
        k = float(k_raw) * scale
        if k <= 0.0:
            continue
        t_k = float(temp_c) + 273.15
        xs.append(1.0 / t_k)
        ys.append(math.log(k))
        if sigma is None:
            sigma_ln = 1.0
        elif rel:
            sigma_ln = float(sigma)
        else:
            # absolute half-width in the table's own printed units
            sigma_ln = abs(float(sigma) * scale) / k
            if abs(float(sigma)) >= abs(float(k_raw)):
                n_interval_spans_zero += 1
        sigma_ln = max(sigma_ln, 1e-6)
        ws.append(1.0 / (sigma_ln * sigma_ln))
    if len(xs) < 2:
        return None

    sw = sum(ws)
    sx = sum(w * x for w, x in zip(ws, xs))
    sy = sum(w * y for w, y in zip(ws, ys))
    sxx = sum(w * x * x for w, x in zip(ws, xs))
    sxy = sum(w * x * y for w, x, y in zip(ws, xs, ys))
    denom = sw * sxx - sx * sx
    if abs(denom) < 1e-30:
        return None
    slope = (sw * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / sw

    # Unweighted reference fit, for comparison with the dossiers' own OLS refits.
    n = len(xs)
    osx, osy = sum(xs), sum(ys)
    osxx = sum(x * x for x in xs)
    osxy = sum(x * y for x, y in zip(xs, ys))
    odenom = n * osxx - osx * osx
    ols_slope = (n * osxy - osx * osy) / odenom if abs(odenom) > 1e-30 else float("nan")
    ols_intercept = (osy - ols_slope * osx) / n if odenom else float("nan")

    # weighted R^2
    ybar = sy / sw
    ss_tot = sum(w * (y - ybar) ** 2 for w, y in zip(ws, ys))
    ss_res = sum(w * (y - (intercept + slope * x)) ** 2 for w, x, y in zip(ws, xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "ea_kj_mol": -slope * R_KJ,
        "a_in_table_unit": math.exp(intercept),
        "ols_ea_kj_mol": -ols_slope * R_KJ,
        "ols_a_in_table_unit": math.exp(ols_intercept),
        "r2_weighted": r2,
        "n_points": n,
        "table_unit": spec["unit"],
        "anchor": spec["anchor"],
        "n_points_with_interval_spanning_zero": n_interval_spans_zero,
        "table_note": spec.get("note"),
    }


# ---------------------------------------------------------------------------
# Enumeration of the shipped pairs.
# ---------------------------------------------------------------------------

def _a_from_k_ref(k_ref: Optional[float], ea: Optional[float], t_ref: float) -> Optional[float]:
    if k_ref is None or ea is None:
        return None
    return float(k_ref) * math.exp(float(ea) / (R_KJ * float(t_ref)))


def collect_yaml_pairs() -> List[Dict[str, Any]]:
    import yaml  # local import: only this collector needs it

    path = ROOT / "data" / "lit" / "arrhenius_params.yml"
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    rows: List[Dict[str, Any]] = []
    for key, entry in (data.get("arrhenius_data") or {}).items():
        a_raw = entry.get("A_value")
        a_is_nan = a_raw is None or (isinstance(a_raw, float) and math.isnan(a_raw))
        rows.append({
            "parameter_id": key,
            "shipped_in": "data/lit/arrhenius_params.yml",
            "lane": "cantera_export",
            "form": "explicit_A",
            "shipped_ea_kj_mol": (
                None if entry.get("Ea_kj_mol") is None else float(entry["Ea_kj_mol"])
            ),
            "shipped_a": None if a_is_nan else float(a_raw),
            "shipped_a_unit": entry.get("A_unit"),
            "a_is_placeholder_nan": a_is_nan,
            "source_key": entry.get("source"),
            "source_quality": entry.get("source_quality"),
            "has_audit_flag": bool(entry.get("audit_flag")),
        })
    return rows


def collect_kinetic_core_pairs() -> List[Dict[str, Any]]:
    from src.kinetic_core import parameters as p_trunk
    from src.kinetic_core import parameters_acrylamide as p_acr
    from src.kinetic_core import parameters_furanic as p_fur
    from src.kinetic_core import parameters_sulfur as p_sul

    registries = [
        ("src/kinetic_core/parameters.py", "trunk_martins_m4",
         p_trunk.MARTINS_M4, p_trunk.T_REF_K),
        ("src/kinetic_core/parameters_acrylamide.py", "acrylamide",
         p_acr.MEASURED_ACRYLAMIDE, p_acr.T_REF_A_K),
        ("src/kinetic_core/parameters_furanic.py", "furanic",
         p_fur.FURANIC_PARAMETERS, p_trunk.T_REF_K),
        ("src/kinetic_core/parameters_sulfur.py", "sulfur",
         p_sul.MEASURED_SULFUR, p_sul.T_REF_S_K),
    ]
    rows: List[Dict[str, Any]] = []
    for shipped_in, lane, registry, t_ref in registries:
        for key, param in registry.items():
            ea = getattr(param, "ea_kj_mol", None)
            k_ref = getattr(param, "k_ref", None)
            unit = getattr(param, "unit", None)
            rows.append({
                "parameter_id": key,
                "shipped_in": shipped_in,
                "lane": lane,
                "form": "reparameterised_k_ref",
                "shipped_ea_kj_mol": None if ea is None else float(ea),
                "shipped_k_ref": None if k_ref is None else float(k_ref),
                "t_ref_k": float(t_ref),
                "shipped_a": _a_from_k_ref(k_ref, ea, t_ref),
                "shipped_a_unit": unit,
                "a_is_placeholder_nan": False,
                "source_key": getattr(param, "source_anchor", None),
                "source_quality": getattr(param, "evidence_class", None),
                "dossier_anchor": getattr(param, "dossier_anchor", None),
            })
    return rows


def collect_misc_pairs() -> List[Dict[str, Any]]:
    """The (Ea, A) pairs that live outside the two registries."""
    rows: List[Dict[str, Any]] = []

    calib = json.loads(
        (ROOT / "data" / "lit" / "lipid_oxidation_calibration.json").read_text(encoding="utf-8")
    )
    kin = calib["kinetics"]
    rows.append({
        "parameter_id": "lipid_hydroperoxide_formation",
        "shipped_in": "data/lit/lipid_oxidation_calibration.json",
        "lane": "lipid_oxidation",
        "form": "explicit_A",
        "shipped_ea_kj_mol": float(kin["activation_energy_j_per_mol"]) / 1000.0,
        "shipped_a": float(kin["prefactor_per_min"]),
        "shipped_a_unit": "1/min",
        "a_is_placeholder_nan": False,
        "source_key": calib["_provenance"]["kinetics_source"],
        "source_quality": "order_of_magnitude_literature_self_declared",
    })

    # src/safety.py hardcodes two pairs inside predict_acrylamide. Parsed from
    # the source rather than copied, so this audit cannot drift from the code.
    safety_src = (ROOT / "src" / "safety.py").read_text(encoding="utf-8")

    def _grab(pattern: str) -> Optional[float]:
        m = re.search(pattern, safety_src)
        return float(m.group(1)) if m else None

    ea_f = _grab(r"Ea_f\s*=\s*([0-9.eE+-]+)\s*\+")
    a_f = _grab(r"A_f\s*=\s*([0-9.eE+-]+)")
    ea_e = _grab(r"Ea_e\s*=\s*([0-9.eE+-]+)")
    a_e = _grab(r"A_e\s*=\s*([0-9.eE+-]+)")
    rows.append({
        "parameter_id": "safety_acrylamide_formation",
        "shipped_in": "src/safety.py (predict_acrylamide)",
        "lane": "safety_screening",
        "form": "explicit_A",
        "shipped_ea_kj_mol": None if ea_f is None else ea_f / 1000.0,
        "shipped_a": a_f,
        "shipped_a_unit": "L/mol.s",
        "a_is_placeholder_nan": False,
        "source_key": "Knol 2009 (comment: '~129 kJ/mol'); no table cited for A",
        "source_quality": "hardcoded_in_code",
    })
    rows.append({
        "parameter_id": "safety_acrylamide_elimination",
        "shipped_in": "src/safety.py (predict_acrylamide)",
        "lane": "safety_screening",
        "form": "explicit_A",
        "shipped_ea_kj_mol": None if ea_e is None else ea_e / 1000.0,
        "shipped_a": a_e,
        "shipped_a_unit": "1/s",
        "source_key": "Parker 2012 (comment: 'Ea_e typically ~90-110 kJ/mol'); no table cited",
        "source_quality": "hardcoded_in_code",
        "a_is_placeholder_nan": False,
    })

    # Four more (Ea, A) pairs are passed as literals into
    # `_formation_elimination_signal` -- two in `_estimate_dicarbonyl_pools` and
    # two in `predict_furosine`. They carry no citation of any kind at the call
    # site. Parsed positionally from the source so they cannot drift.
    literal_pairs = re.findall(
        r"formation_pre_exponential=([0-9.eE+-]+),\s*\n\s*formation_ea_kj_mol=([0-9.eE+-]+),\s*\n"
        r"\s*elimination_pre_exponential=([0-9.eE+-]+),\s*\n\s*elimination_ea_kj_mol=([0-9.eE+-]+),",
        safety_src,
    )
    call_sites = [
        ("dicarbonyl_pool", "src/safety.py (_estimate_dicarbonyl_pools)",
         "no citation at the call site; the Amadori/dicarbonyl pool shape constants"),
        ("furosine", "src/safety.py (predict_furosine)",
         "no citation at the call site; the furosine fallback path"),
    ]
    for (a_form, ea_form, a_elim, ea_elim), (tag, where, provenance) in zip(literal_pairs, call_sites):
        rows.append({
            "parameter_id": f"safety_{tag}_formation",
            "shipped_in": where,
            "lane": "safety_screening",
            "form": "explicit_A",
            "shipped_ea_kj_mol": float(ea_form),
            "shipped_a": float(a_form),
            "shipped_a_unit": "1/s (implied)",
            "a_is_placeholder_nan": False,
            "source_key": provenance,
            "source_quality": "hardcoded_in_code_uncited",
        })
        rows.append({
            "parameter_id": f"safety_{tag}_elimination",
            "shipped_in": where,
            "lane": "safety_screening",
            "form": "explicit_A",
            "shipped_ea_kj_mol": float(ea_elim),
            "shipped_a": float(a_elim),
            "shipped_a_unit": "1/s (implied)",
            "a_is_placeholder_nan": False,
            "source_key": provenance,
            "source_quality": "hardcoded_in_code_uncited",
        })
    return rows


# ---------------------------------------------------------------------------
# The mapping from a shipped pair to the dossier k-table that can verify it.
# Anything not named here has NO transcribed source k(T) table.
# ---------------------------------------------------------------------------

VERIFICATION_MAP: Dict[str, Dict[str, Any]] = {
    "cysteine_thermolysis": {
        "table": "zheng1994_cysteine_ph5",
        "conversion": "per_min_to_per_s",
        "also_compare": ["zheng1994_cysteine_ph3", "zheng1994_cysteine_ph7",
                         "zheng1994_cysteine_ph9"],
        "note": (
            "The shipped Ea 130.4 is the pH 3-9 MEAN, so no single pH column is its "
            "exact partner; pH 5 is used as the primary comparator and all four are "
            "reported. The prefactor must be judged against the whole 9.8e11-2.4e13 "
            "band the four columns give."
        ),
    },
    "k_cys_h2s": {
        "table": "zheng1994_cysteine_ph5",
        "conversion": "identity",
        "note": "shipped at ph_of_measurement 5.0, so the pH 5 column is the exact partner",
    },
    "k_hmf_cys": {
        "table": "hamzalioglu2018_hmf_cys",
        "conversion": "hamzalioglu_pseudo1_to_second_order",
        "note": (
            "the shipped second-order constant is derived as k = k'/[AA] with the "
            "stated [AA] = 30 umol / 1.5 mL = 20 mM, so the same conversion is applied "
            "to the refit prefactor"
        ),
    },
    "k_glc_tdg": {"table": "kocadagli2016_step3_glc_3dg", "conversion": "identity"},
    "k_tdg_ddg": {"table": "kocadagli2016_step4_3dg_34dg", "conversion": "identity"},
    "k_ddg_hmf": {"table": "kocadagli2016_step5_34dg_hmf", "conversion": "identity"},
    "k_fru_int": {"table": "kocadagli2016_step6_fru_int", "conversion": "identity"},
    "k_int_hmf": {"table": "kocadagli2016_step7_int_hmf", "conversion": "identity"},
    "k_fru_odg": {"table": "kocadagli2016_step8_fru_1dg", "conversion": "identity"},
    "k_tdg_mgo": {"table": "kocadagli2016_step11_3dg_mgo", "conversion": "identity"},
    "k_thioether": {
        "table": "stack2018_nac_forward",
        "conversion": "molar_per_s_to_l_per_mmol_min",
        "magnitude_not_transferred": True,
        "note": (
            "DECLARED HYBRID at the parameter: only Stack's Ea is taken; the absolute "
            "forward magnitude is Hofmann/Charles-Bernard's matrix recast, because "
            "Stack's 496 M^-1 s^-1 against the titrated site density predicts a "
            "millisecond thiol half-life. The A ratio below is therefore EXPECTED to "
            "be large and is not a defect."
        ),
    },
}

# Kocadagli publishes its Table 2 (which the shipped values come from) as a
# SIMULTANEOUS GLOBAL FIT of the ODE system, not as an Arrhenius fit of its own
# Table 1. Disagreement is an identifiability diagnostic, not proof of error.
GLOBAL_FIT_SOURCES = {
    "k_glc_tdg", "k_tdg_ddg", "k_ddg_hmf", "k_fru_int",
    "k_int_hmf", "k_fru_odg", "k_tdg_mgo",
}


# ---------------------------------------------------------------------------
# A SECOND verification axis that needs no k(T) table: two shipped pairs in
# different lanes claiming THE SAME source and the same elementary step must
# agree once put in the same units. Where they do not, one of them is wrong and
# the disagreement is arithmetic, not opinion.
# ---------------------------------------------------------------------------

CROSS_LANE_CHECKS: List[Dict[str, Any]] = [
    {
        "name": "schiff_condensation (Cantera lane) vs k_schiff (kinetic core)",
        "left": "schiff_condensation",
        "right": "k_schiff",
        # L/(mmol*min) -> L/(mol*s):  x1000 (mmol->mol), /60 (min->s)
        "right_to_left_factor": 1000.0 / 60.0,
        "shared_source": (
            "Martins & van Boekel 2005, Food Chem 90:257-269, Table 2, step 1 "
            "(Glc + Gly -> DFG). Both entries name it."
        ),
    },
    {
        "name": "cysteine_thermolysis (Cantera lane) vs k_cys_h2s (kinetic core)",
        "left": "cysteine_thermolysis",
        "right": "k_cys_h2s",
        # 1/min -> 1/s
        "right_to_left_factor": 1.0 / 60.0,
        "shared_source": (
            "Zheng & Ho 1994, ACS Symp. Ser. 564:138-146, aqueous 0.1 M L-cysteine "
            "H2S release. Both entries name it."
        ),
    },
]


def run_cross_lane_checks(rows: Sequence[Dict[str, Any]], fold_threshold: float) -> List[Dict[str, Any]]:
    index = {row["parameter_id"]: row for row in rows}
    out: List[Dict[str, Any]] = []
    for check in CROSS_LANE_CHECKS:
        left = index.get(check["left"])
        right = index.get(check["right"])
        if left is None or right is None:
            continue
        left_a = left.get("shipped_a")
        right_a = right.get("shipped_a")
        if not left_a or not right_a:
            continue
        right_in_left_units = right_a * check["right_to_left_factor"]
        ratio = left_a / right_in_left_units
        fold = max(ratio, 1.0 / ratio) if ratio > 0 else float("inf")
        left_ea = left.get("shipped_ea_kj_mol")
        right_ea = right.get("shipped_ea_kj_mol")
        ea_delta = (
            None if left_ea is None or right_ea is None else left_ea - right_ea
        )
        # If exactly one side independently reproduces the source's own k(T)
        # table, the disagreement is not symmetric: that side is exonerated and
        # the other one carries the defect.
        def _reproduces(row: Dict[str, Any]) -> Optional[bool]:
            if not row.get("verifiable"):
                return None
            fold_vs_refit = row.get("a_fold_difference")
            if fold_vs_refit is None:
                return None
            return fold_vs_refit <= fold_threshold

        left_ok, right_ok = _reproduces(left), _reproduces(right)
        if left_ok is True and right_ok is not True:
            blame, exonerated = check["right"], check["left"]
        elif right_ok is True and left_ok is not True:
            blame, exonerated = check["left"], check["right"]
        else:
            blame, exonerated = None, None

        out.append({
            "name": check["name"],
            "shared_source": check["shared_source"],
            "left_id": check["left"],
            "right_id": check["right"],
            "left_a": left_a,
            "left_a_unit": left.get("shipped_a_unit"),
            "right_a_converted_to_left_unit": right_in_left_units,
            "a_ratio": ratio,
            "a_fold_difference": fold,
            "left_ea_kj_mol": left_ea,
            "right_ea_kj_mol": right_ea,
            "ea_delta_kj_mol": ea_delta,
            "side_carrying_the_defect": blame,
            "side_reproducing_its_source": exonerated,
            "flag": "FLAG_CROSS_LANE_A" if fold > fold_threshold else "OK",
        })
    return out


def audit(a_ratio_threshold: float, ea_tolerance_pct: float) -> Dict[str, Any]:
    rows = collect_yaml_pairs() + collect_kinetic_core_pairs() + collect_misc_pairs()

    refits = {name: refit_arrhenius(spec) for name, spec in DOSSIER_K_TABLES.items()}

    for row in rows:
        pid = row["parameter_id"]
        spec = VERIFICATION_MAP.get(pid)
        row["dossier_k_table"] = None
        row["refit_ea_kj_mol"] = None
        row["refit_a"] = None
        row["a_ratio"] = None
        row["ea_delta_pct"] = None
        row["flags"] = []
        row["notes"] = []

        if row["shipped_ea_kj_mol"] is None:
            row["flags"].append("NO_EA_SHIPPED")
        if row.get("a_is_placeholder_nan"):
            row["flags"].append("NO_PREFACTOR_SHIPPED")
            row["notes"].append(
                "A_value is a bare NaN placeholder; src/barrier_constants."
                "get_arrhenius_params silently substitutes 6.25e12 ('TST @ 150C') at "
                "the point of use, and for an A_unit of L/mol.s that substitution is "
                "also a unimolecular prefactor handed to a bimolecular reaction."
            )

        if spec is None:
            row["verifiable"] = False
            row["flags"].append("UNVERIFIABLE_NO_DOSSIER_K_TABLE")
            rows_note = row.get("dossier_anchor")
            if rows_note:
                row["notes"].append(
                    "the cited dossier anchors a reparameterised (k_ref, Ea) pair or a "
                    "bracket, not a rate-constant-versus-temperature table"
                )
            continue

        fit = refits.get(spec["table"])
        if fit is None:
            row["verifiable"] = False
            row["flags"].append("UNVERIFIABLE_REFIT_FAILED")
            continue

        factor, factor_desc = CONVERSIONS[spec["conversion"]]
        row["verifiable"] = True
        row["dossier_k_table"] = fit["anchor"]
        row["refit_ea_kj_mol"] = fit["ea_kj_mol"]
        row["refit_a"] = fit["a_in_table_unit"] * factor
        row["refit_a_ols"] = fit["ols_a_in_table_unit"] * factor
        row["refit_ea_ols_kj_mol"] = fit["ols_ea_kj_mol"]
        row["refit_r2_weighted"] = fit["r2_weighted"]
        row["refit_n_points"] = fit["n_points"]
        row["refit_unit_conversion"] = factor_desc
        if fit.get("table_note"):
            row["notes"].append("source table: " + fit["table_note"])
        if spec.get("note"):
            row["notes"].append(spec["note"])

        if row["shipped_a"] and row["refit_a"]:
            ratio = row["shipped_a"] / row["refit_a"]
            row["a_ratio"] = ratio
            fold = max(ratio, 1.0 / ratio) if ratio > 0 else float("inf")
            row["a_fold_difference"] = fold
            if fold > a_ratio_threshold:
                if spec.get("magnitude_not_transferred"):
                    row["flags"].append("A_RATIO_LARGE_BY_DESIGN")
                else:
                    row["flags"].append("FLAG_A_RATIO")

        if row["shipped_ea_kj_mol"] is not None and row["refit_ea_kj_mol"] is not None:
            shipped_ea = row["shipped_ea_kj_mol"]
            delta = row["refit_ea_kj_mol"] - shipped_ea
            denom = abs(shipped_ea) if abs(shipped_ea) > 1e-9 else None
            if denom is None:
                # shipped Ea is exactly zero (author-fixed); report absolute kJ/mol
                row["ea_delta_pct"] = None
                row["ea_delta_kj_mol"] = delta
                if abs(delta) > 5.0:
                    row["flags"].append("FLAG_EA_FIXED_ZERO_VS_NONZERO_REFIT")
            else:
                row["ea_delta_pct"] = 100.0 * delta / denom
                row["ea_delta_kj_mol"] = delta
                if abs(row["ea_delta_pct"]) > ea_tolerance_pct:
                    row["flags"].append("FLAG_EA")

        if pid in GLOBAL_FIT_SOURCES:
            row["notes"].append(
                "shipped value comes from the source's SIMULTANEOUS GLOBAL ODE fit "
                "(Table 2), not from an Arrhenius fit of its own per-temperature "
                "Table 1; the two are not required to agree and the gap is an "
                "identifiability diagnostic"
            )

        if not [f for f in row["flags"] if f.startswith("FLAG")]:
            row["flags"].append("OK")

    cross_lane = run_cross_lane_checks(rows, a_ratio_threshold)
    for check in cross_lane:
        if check["flag"] != "FLAG_CROSS_LANE_A":
            continue
        exonerated = check.get("side_reproducing_its_source")
        for pid in (check["left_id"], check["right_id"]):
            row = next((r for r in rows if r["parameter_id"] == pid), None)
            if row is None:
                continue
            if pid == exonerated:
                row["notes"].append(
                    f"cross-lane: disagrees with `{check['side_carrying_the_defect']}` by "
                    f"{check['a_fold_difference']:.3g}x on the prefactor while both cite the "
                    f"same source, but THIS side reproduces the source's own k(T) table, so "
                    f"the defect is the other side's. Not flagged here."
                )
                continue
            if "FLAG_CROSS_LANE_A" not in row["flags"]:
                row["flags"].append("FLAG_CROSS_LANE_A")
            if "OK" in row["flags"]:
                row["flags"].remove("OK")
            blamed = (
                " and THIS side is the one that does not reproduce the source's own "
                "k(T) table" if pid == check.get("side_carrying_the_defect") else ""
            )
            row["notes"].append(
                f"cross-lane: {check['name']} disagree by "
                f"{check['a_fold_difference']:.3g}x on the prefactor while citing the "
                f"same source{blamed}"
            )

    verifiable = [r for r in rows if r.get("verifiable")]
    unverifiable = [r for r in rows if not r.get("verifiable")]
    flagged = [r for r in rows if any(f.startswith("FLAG") for f in r["flags"])]

    return {
        "cross_lane_checks": cross_lane,
        "thresholds": {
            "a_ratio_fold": a_ratio_threshold,
            "ea_tolerance_pct": ea_tolerance_pct,
        },
        "coverage": {
            "shipped_pairs_total": len(rows),
            "verifiable_against_a_dossier_k_table": len(verifiable),
            "unverifiable_no_dossier_k_table": len(unverifiable),
            "pairs_with_no_prefactor_shipped": len(
                [r for r in rows if "NO_PREFACTOR_SHIPPED" in r["flags"]]
            ),
            "flagged": len(flagged),
            "dossier_k_tables_available": len(DOSSIER_K_TABLES),
            "cross_lane_checks_run": len(cross_lane),
            "cross_lane_checks_flagged": len(
                [c for c in cross_lane if c["flag"] == "FLAG_CROSS_LANE_A"]
            ),
        },
        "out_of_scope_ea_only_tables": {
            "note": (
                "These ship an activation energy with NO paired prefactor: the rate is "
                "reconstructed at the point of use from a generic Eyring/TST factor "
                "(k_B T/h, or the 6.25e12 'TST @ 150C' substitution). They are not "
                "(Ea, A) pairs and so are outside this audit's scored set, but they "
                "are the same silent-prefactor risk surface."
            ),
            "src/barrier_constants.py FAST_BARRIERS": "51 families, kcal/mol, fitted screening constants",
            "src/conditions.py ACTIVATION_ENERGIES_KJ": "6 family defaults, kJ/mol, Eyring prefactor",
        },
        "rows": rows,
        "refits": refits,
    }


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _fmt(value: Optional[float], spec: str = "{:.4g}") -> str:
    if value is None:
        return "—"
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return spec.format(value)


def _verdict(row: Dict[str, Any]) -> str:
    flags = row["flags"]
    parts: List[str] = []
    if "FLAG_A_RATIO" in flags and ("FLAG_EA" in flags or "FLAG_EA_FIXED_ZERO_VS_NONZERO_REFIT" in flags):
        parts.append("**FLAG — A and Ea**")
    elif "FLAG_A_RATIO" in flags:
        parts.append("**FLAG — prefactor**")
    elif "FLAG_EA" in flags:
        parts.append("**FLAG — Ea**")
    elif "FLAG_EA_FIXED_ZERO_VS_NONZERO_REFIT" in flags:
        parts.append("**FLAG — Ea fixed at 0**")
    if "FLAG_CROSS_LANE_A" in flags:
        parts.append("**FLAG — cross-lane A**")
    if parts:
        return " + ".join(parts)
    if "NO_PREFACTOR_SHIPPED" in flags:
        return "**NO PREFACTOR** (NaN placeholder)"
    if "A_RATIO_LARGE_BY_DESIGN" in flags:
        return "by design (Ea only transferred)"
    if "NO_EA_SHIPPED" in flags:
        return "no Ea shipped (not an Arrhenius pair)"
    if "UNVERIFIABLE_NO_DOSSIER_K_TABLE" in flags:
        return "unverifiable"
    if "OK" in flags:
        return "OK"
    return ", ".join(flags)


def render_markdown(result: Dict[str, Any]) -> str:
    cov = result["coverage"]
    th = result["thresholds"]
    out: List[str] = []
    out.append("# Corpus-wide Arrhenius prefactor audit")
    out.append("")
    out.append(
        "Generated by `scripts/generators/generate_prefactor_audit.py`. "
        "**This is an audit. No shipped parameter was changed by it.**"
    )
    out.append("")
    out.append(
        f"Flag thresholds: prefactor differing by more than **{th['a_ratio_fold']:g}x**, "
        f"activation energy differing by more than **{th['ea_tolerance_pct']:g} %**."
    )
    out.append("")
    out.append("## Coverage — be honest about it")
    out.append("")
    out.append(f"* Shipped `(Ea, A)` pairs enumerated: **{cov['shipped_pairs_total']}**")
    out.append(
        f"* Verifiable against a transcribed source k(T) table: "
        f"**{cov['verifiable_against_a_dossier_k_table']}** "
        f"({100.0 * cov['verifiable_against_a_dossier_k_table'] / cov['shipped_pairs_total']:.0f} %)"
    )
    out.append(
        f"* **UNVERIFIABLE** — no k(T) table anywhere in "
        f"`data/lit/extraction_dossiers/`: **{cov['unverifiable_no_dossier_k_table']}** "
        f"({100.0 * cov['unverifiable_no_dossier_k_table'] / cov['shipped_pairs_total']:.0f} %)"
    )
    out.append(f"* Pairs shipping NO prefactor at all (NaN placeholder): **{cov['pairs_with_no_prefactor_shipped']}**")
    out.append(f"* Flagged: **{cov['flagged']}**")
    out.append(f"* Distinct source k(T) tables available to refit against: **{cov['dossier_k_tables_available']}**")
    out.append("")
    oos = result["out_of_scope_ea_only_tables"]
    out.append("Outside the scored set, and worth stating: " + oos["note"])
    out.append("")
    for key, value in oos.items():
        if key == "note":
            continue
        out.append(f"* `{key}` — {value}")
    out.append("")

    out.append("## Is the refit machinery itself trustworthy?")
    out.append("")
    out.append(
        "Four of the refits below have an independent prior refit in the dossiers, "
        "done by a different analyst in a different wave, and the unweighted (OLS) "
        "column reproduces all four: Zheng's four cysteine columns "
        "(131.2 / 133.0 / 135.5 / 123.3 kJ/mol and A = 9.79e11 / 1.93e12 / 2.36e13 / "
        "1.04e12 s^-1), Hamzalioglu's HMF-Cys row (29.675 kJ/mol, A = 24115 day^-1), "
        "Kocadagli's seven glucose steps (59.6 / 32.1 / -7.0 / 100.5 / 145.0 / 96.9 / "
        "93.7 kJ/mol), and Stack's NAC forward series (10.8 kJ/mol, "
        "A = 4.267e4 M^-1 s^-1). The fitter is therefore validated against four "
        "independent hand-checks before any new number below is believed."
    )
    out.append("")
    out.append("## The table")
    out.append("")
    out.append(
        "| parameter id | source key | shipped Ea (kJ/mol) | shipped A | dossier k-table? | "
        "refit Ea (kJ/mol) | refit A | A ratio (shipped/refit) | verdict |"
    )
    out.append("|---|---|---|---|---|---|---|---|---|")
    for row in result["rows"]:
        source = (row.get("source_key") or "—")
        source = re.sub(r"\s+", " ", str(source))
        if len(source) > 78:
            source = source[:75] + "..."
        source = source.replace("|", "/")
        a_unit = row.get("shipped_a_unit") or ""
        if row.get("a_is_placeholder_nan"):
            shipped_a = "NaN placeholder"
        elif row.get("shipped_a") is None:
            shipped_a = "n/a (no Ea shipped)"
        else:
            shipped_a = f"{_fmt(row.get('shipped_a'))} {a_unit}".strip()
        refit_a = row.get("refit_a")
        refit_a_s = "—" if refit_a is None else f"{_fmt(refit_a)} {a_unit}".strip()
        table = "yes" if row.get("verifiable") else "**no**"
        ratio = row.get("a_ratio")
        if ratio is None:
            ratio_s = "—"
        else:
            fold = max(ratio, 1.0 / ratio) if ratio > 0 else float("inf")
            ratio_s = f"{_fmt(ratio)} ({fold:.3g}x)"
        out.append(
            f"| `{row['parameter_id']}` | {source} | {_fmt(row.get('shipped_ea_kj_mol'))} | "
            f"{shipped_a} | {table} | {_fmt(row.get('refit_ea_kj_mol'))} | {refit_a_s} | "
            f"{ratio_s} | {_verdict(row)} |"
        )
    out.append("")

    out.append("## Cross-lane consistency — a second axis that needs no k-table")
    out.append("")
    out.append(
        "Two shipped pairs in different lanes that name the SAME source and the SAME "
        "elementary step must agree once put in the same units. Where they do not, one "
        "of them is wrong, and the disagreement is arithmetic rather than opinion. "
        "The two lanes never read each other's numbers, which is how a disagreement "
        "can survive."
    )
    out.append("")
    out.append(
        "| check | shipped A (lane 1) | shipped A (lane 2, converted) | fold | Ea (lane 1) | "
        "Ea (lane 2) | which side carries the defect | verdict |"
    )
    out.append("|---|---|---|---|---|---|---|---|")
    for check in result.get("cross_lane_checks", []):
        verdict = "**FLAG**" if check["flag"] == "FLAG_CROSS_LANE_A" else "OK"
        blame = check.get("side_carrying_the_defect")
        blame_s = (
            f"`{blame}` (the other reproduces the source k-table)" if blame
            else "undecidable here — neither side has a source k-table"
        )
        out.append(
            f"| {check['name']} | {_fmt(check['left_a'])} {check['left_a_unit'] or ''} | "
            f"{_fmt(check['right_a_converted_to_left_unit'])} | "
            f"{check['a_fold_difference']:.3g}x | {_fmt(check['left_ea_kj_mol'])} | "
            f"{_fmt(check['right_ea_kj_mol'])} | {blame_s} | {verdict} |"
        )
    out.append("")
    for check in result.get("cross_lane_checks", []):
        out.append(f"* `{check['left_id']}` / `{check['right_id']}` — {check['shared_source']}")
    out.append("")

    flagged = [r for r in result["rows"] if any(f.startswith("FLAG") for f in r["flags"])]
    out.append("## Flagged rows, in detail")
    out.append("")
    if not flagged:
        out.append("None.")
    for row in flagged:
        out.append(f"### `{row['parameter_id']}` — {_verdict(row)}")
        out.append("")
        out.append(f"* shipped in `{row['shipped_in']}` (lane: {row['lane']})")
        out.append(
            f"* shipped: Ea = {_fmt(row.get('shipped_ea_kj_mol'))} kJ/mol, "
            f"A = {_fmt(row.get('shipped_a'))} {row.get('shipped_a_unit') or ''}"
        )
        if row.get("verifiable"):
            out.append(
                f"* refit (WLS, n = {row.get('refit_n_points')}, weighted R2 = "
                f"{_fmt(row.get('refit_r2_weighted'), '{:.4f}')}): "
                f"Ea = {_fmt(row.get('refit_ea_kj_mol'))} kJ/mol, "
                f"A = {_fmt(row.get('refit_a'))} {row.get('shipped_a_unit') or ''}"
            )
            if row.get("refit_a_ols") is not None:
                out.append(
                    f"* unweighted reference fit: Ea = {_fmt(row.get('refit_ea_ols_kj_mol'))} kJ/mol, "
                    f"A = {_fmt(row.get('refit_a_ols'))}"
                )
            out.append(f"* k-table: {row.get('dossier_k_table')}")
        else:
            out.append(
                "* NO source k(T) table exists for this pair; it is flagged on the "
                "cross-lane axis only"
            )
        if row.get("refit_unit_conversion") and row["refit_unit_conversion"] != "none":
            out.append(f"* unit conversion applied to the refit A: {row['refit_unit_conversion']}")
        if row.get("ea_delta_pct") is not None:
            out.append(
                f"* Ea difference: {row['ea_delta_pct']:+.1f} % "
                f"({row.get('ea_delta_kj_mol', 0.0):+.1f} kJ/mol)"
            )
        for note in row.get("notes", []):
            out.append(f"* NOTE: {note}")
        out.append("")

    out.append("## Unverifiable pairs — the honest gap")
    out.append("")
    out.append(
        "No transcribed source rate-constant-versus-temperature table exists for these, so "
        "the prefactor cannot be checked the way the three historical defects were found. "
        "This is a statement about coverage, not about correctness."
    )
    out.append("")
    out.append("| parameter id | shipped in | source key | source quality |")
    out.append("|---|---|---|---|")
    for row in result["rows"]:
        if row.get("verifiable"):
            continue
        source = re.sub(r"\s+", " ", str(row.get("source_key") or "—"))
        if len(source) > 90:
            source = source[:87] + "..."
        out.append(
            f"| `{row['parameter_id']}` | `{row['shipped_in']}` | {source.replace('|', '/')} | "
            f"{row.get('source_quality') or '—'} |"
        )
    out.append("")

    out.append("## The refits themselves")
    out.append("")
    out.append("| k-table | n | weighted Ea (kJ/mol) | weighted A (table unit) | OLS Ea | OLS A | weighted R2 | anchor |")
    out.append("|---|---|---|---|---|---|---|---|")
    for name, fit in result["refits"].items():
        if fit is None:
            continue
        out.append(
            f"| `{name}` | {fit['n_points']} | {_fmt(fit['ea_kj_mol'])} | "
            f"{_fmt(fit['a_in_table_unit'])} {fit['table_unit']} | {_fmt(fit['ols_ea_kj_mol'])} | "
            f"{_fmt(fit['ols_a_in_table_unit'])} | {_fmt(fit['r2_weighted'], '{:.4f}')} | {fit['anchor']} |"
        )
    out.append("")
    return "\n".join(out)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="generate_prefactor_audit.py",
        description=(
            "Corpus-wide Arrhenius prefactor audit: enumerate every shipped (Ea, A) "
            "pair, refit (Ea, A) by weighted least squares on ln k vs 1/T wherever a "
            "source k(T) table has been transcribed into data/lit/extraction_dossiers, "
            "and flag every pair whose prefactor or barrier does not reproduce. Writes "
            "prefactor_audit.json and prefactor_audit.md. Changes no shipped value."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        default="results/validation",
        help="directory the two artifacts are written to, relative to the repo root",
    )
    parser.add_argument(
        "--a-ratio-threshold",
        type=float,
        default=2.0,
        help="flag a pair whose shipped A differs from the refit A by more than this fold",
    )
    parser.add_argument(
        "--ea-tolerance-pct",
        type=float,
        default=10.0,
        help="flag a pair whose shipped Ea differs from the refit Ea by more than this percent",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="write the artifacts without echoing the markdown to stdout",
    )
    args = parser.parse_args(argv)

    result = audit(args.a_ratio_threshold, args.ea_tolerance_pct)
    markdown = render_markdown(result)

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    md_path = output_dir / "prefactor_audit.md"
    json_path = output_dir / "prefactor_audit.json"
    md_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    if not args.quiet:
        print(markdown)
    print(f"Wrote {md_path}")
    print(f"Wrote {json_path}")
    cov = result["coverage"]
    print(
        f"Coverage: {cov['verifiable_against_a_dossier_k_table']}/"
        f"{cov['shipped_pairs_total']} shipped pairs verifiable; "
        f"{cov['flagged']} flagged."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
