#!/usr/bin/env python
"""
BUILD WAVE B7 -- THE FURANIC CHANNEL'S FIT (2026-08-29).

Pre-registered in ``results/validation/kinetic_core_b7_prereg.md``, written
before this file existed. Read it first: it states the objective, the two
seeded starts, the expected residual size and the four things that would
falsify the build.

WHAT IS ACTUALLY FITTED HERE: ONE NUMBER.
-----------------------------------------
``k_dpo_af`` -- the 1-deoxypentosone + glycine -> acetylformoin edge -- against
THREE cells of Blank, Fay, Lakner & Schlosser 1997 Table 1, a declared FIT set
(FIT_HOLDOUT_DECLARATION Amendment 12). Everything else in the channel is
INGESTED (Kocadagli's seven glucose constants, Hamzalioglu's sink), DERIVED
(the hexose transfer), a DIGITISED PRIOR (Edge B's level) or EXACTLY ZERO
(Edge C).

This script therefore also does four things that are not fits and are reported
beside it, because a module with one fitted number has to defend the ingested
ones some other way:

  1. an INGESTION CHECK on Hamzalioglu's sink: the refit Arrhenius line is
     evaluated at 5 / 25 / 50 C and compared with the three second-order
     constants derived from the paper's own table;
  2. an INGESTION CHECK on Kocadagli's seven constants: each is re-referenced
     from T_b = 180 C to 100 C and back, and must return the published k_b;
  3. the SENSITIVITY SWEEP the pre-registration demands on ``k_af_dmhf``, the
     one constant in the channel with no source at all;
  4. a BRANCH-SHARE DEMONSTRATION: the same model run on two different sugar
     charges must give two different HMF limb shares, because there is no
     branch fraction to hard-code.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUTPUT_DIR = ROOT / "results" / "validation"
BASENAME = "kinetic_core_b7_fit_report"

SEED = 20260829
SEEDED_STARTS_LOG10 = (-8.0, -4.0)
START_AGREEMENT_TOLERANCE = 1.0e-6


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"

    return {
        "commit": _run("rev-parse", "HEAD"),
        "short": _run("rev-parse", "--short", "HEAD"),
        "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": "yes" if _run("status", "--porcelain") else "no",
    }


# ---------------------------------------------------------------------------
# The Blank 1997 forward model
# ---------------------------------------------------------------------------


def _blank_prediction_ug_per_mmol(k_dpo_af: float) -> float:
    """
    Predicted HDMF in Blank's own unit, for one value of the fitted constant.

    ONE prediction for all three sugars, and that is a DECLARED limitation:
    the sulfur lane carries one generic aldopentose, so the 1.96x
    arabinose > ribose > xylose spread Blank measures cannot be reproduced and
    shows up as a fit residual instead. Same treatment B2 gives Hofmann's 1.38x
    ribose/xylose gap.
    """
    from src.kinetic_core.engine import (
        SULFUR,
        FormulationSpec,
        ProcessSpec,
        ThermalProgram,
        core_parameters,
        predict,
    )
    from src.kinetic_core.furanic import ug_per_mmol_from_state
    from src.kinetic_core.parameters_furanic import (
        BLANK_AMINE_LOADING_MMOL_PER_L,
        BLANK_PH,
        BLANK_SUGAR_LOADING_MMOL_PER_L,
        BLANK_TEMPERATURE_C,
        BLANK_TIME_MIN,
        with_fitted_furanic,
    )
    from src.kinetic_core.ph_state import BufferSpec
    from src.kinetic_core.species_sulfur import MOLECULAR_WEIGHT_G_PER_MOL

    parameters = dict(core_parameters(SULFUR))
    parameters.update(with_fitted_furanic(float(k_dpo_af)))

    spec = FormulationSpec(
        name="blank1997_pentose_glycine_90C_1h_pH6",
        precursors={
            "xylose": BLANK_SUGAR_LOADING_MMOL_PER_L,
            "glycine": BLANK_AMINE_LOADING_MMOL_PER_L,
        },
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(BLANK_TEMPERATURE_C, BLANK_TIME_MIN),
            ph=BLANK_PH,
            buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                               source="Blank 1997 sec. 2: 0.2 M phosphate"),
        ),
    )
    run = predict(spec, ["DMHF"], parameters=parameters)
    if not run.answered:
        raise SystemExit(
            "the Blank 1997 fit condition is out of envelope: "
            + " | ".join(run.declaration.reasons)
        )
    return ug_per_mmol_from_state(
        float(run.species_mmol_per_l.get("DMHF", 0.0)),
        BLANK_SUGAR_LOADING_MMOL_PER_L,
        MOLECULAR_WEIGHT_G_PER_MOL["DMHF"],
    )


def _objective(log10_k: float, targets: List[float], sigma_log10: float) -> float:
    predicted = _blank_prediction_ug_per_mmol(10.0 ** float(log10_k))
    if predicted <= 0.0:
        return 1.0e9
    return sum(
        ((math.log10(predicted) - math.log10(t)) / sigma_log10) ** 2 for t in targets
    )


def _nelder_mead_1d(
    f, x0: float, step: float = 0.5, tol: float = 1.0e-9, max_iter: int = 200
) -> Tuple[float, float]:
    """A deterministic 1-D simplex. No RNG, so the 'seed' is the start pair."""
    a, b = x0, x0 + step
    fa, fb = f(a), f(b)
    if fb < fa:
        a, b, fa, fb = b, a, fb, fa
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        r = a + (a - b)          # reflect
        fr = f(r)
        if fr < fa:
            e = a + 2.0 * (a - b)   # expand
            fe = f(e)
            b, fb = (e, fe) if fe < fr else (r, fr)
        else:
            c = a + 0.5 * (b - a)   # contract
            fc = f(c)
            b, fb = (c, fc)
        if fb < fa:
            a, b, fa, fb = b, a, fb, fa
    return a, fa


def fit_k_dpo_af() -> Dict[str, Any]:
    from src.kinetic_core.furanic import blank1997_fit_cells

    cells = blank1997_fit_cells("fit")
    targets = [c.ug_per_mmol_sugar for c in cells]
    sigma_log10 = math.log10(1.10)

    def f(log10_k: float) -> float:
        return _objective(log10_k, targets, sigma_log10)

    starts = []
    for start in SEEDED_STARTS_LOG10:
        x, fx = _nelder_mead_1d(f, start)
        starts.append({"start_log10_k": start, "log10_k": x, "objective": fx})

    spread = abs(starts[0]["log10_k"] - starts[1]["log10_k"])
    best = min(starts, key=lambda s: s["objective"])
    k = 10.0 ** best["log10_k"]
    predicted = _blank_prediction_ug_per_mmol(k)

    residuals = {
        f"{c.sugar}/glycine HDMF": math.log10(predicted / c.ug_per_mmol_sugar)
        for c in cells
    }
    rms = math.sqrt(sum(v * v for v in residuals.values()) / len(residuals))
    return {
        "fitted_parameter": "k_dpo_af",
        "unit": "L/(mmol*min)",
        "value": k,
        "log10_value": best["log10_k"],
        "objective": best["objective"],
        "sigma_log10": sigma_log10,
        "sigma_basis": "Blank 1997's own stated maximum SD of 10 %",
        "seed": SEED,
        "starts": starts,
        "start_agreement_log10": spread,
        "starts_agree": spread <= START_AGREEMENT_TOLERANCE,
        "predicted_ug_per_mmol": predicted,
        "targets_ug_per_mmol": {
            f"{c.sugar}/glycine HDMF": c.ug_per_mmol_sugar for c in cells
        },
        "residual_log10_fold": residuals,
        "rms_residual_log10": rms,
        "rms_residual_fold": 10.0 ** rms,
        "residual_is_the_sugar_axis": (
            "The model carries ONE generic aldopentose, so it emits one number "
            "for three sugars and the residual IS Blank's 1.96x sugar spread. "
            "The pre-registration predicted ~0.12 decades on exactly this "
            "ground; a materially smaller RMS would mean something beyond the "
            "one declared parameter had been fitted."
        ),
    }


# ---------------------------------------------------------------------------
# The three ingestion checks and the sweep
# ---------------------------------------------------------------------------


def hamzalioglu_ingestion_check() -> Dict[str, Any]:
    """
    Does the refit Arrhenius line reproduce the three measured second-order
    constants? It is NOT expected to reproduce them exactly -- the refit R^2 is
    0.874 -- and the size of the disagreement is the diagnostic.
    """
    from src.kinetic_core.parameters_furanic import (
        MEASURED_FURANIC,
        m_inv_day_to_core_units,
    )

    measured = {5.0: 3.95, 25.0: 5.15, 50.0: 23.3}   # M^-1 day^-1, K5a 7.1 #1
    parameter = MEASURED_FURANIC["k_hmf_cys"]
    rows = []
    for t_c, k_m_inv_day in measured.items():
        predicted_core = parameter.k_at(t_c + 273.15)
        measured_core = m_inv_day_to_core_units(k_m_inv_day)
        rows.append({
            "temperature_C": t_c,
            "measured_M_inv_day": k_m_inv_day,
            "measured_L_per_mmol_min": measured_core,
            "line_L_per_mmol_min": predicted_core,
            "fold": max(predicted_core, measured_core)
            / min(predicted_core, measured_core),
        })
    return {
        "what": (
            "the DOSSIER'S REFIT prefactor (24115.1 day^-1 pseudo-first-order, "
            "/0.020 M) and Ea (29.675 kJ/mol), evaluated against the three "
            "second-order constants derived from Hamzalioglu's own Table 1"
        ),
        "why_the_refit_and_not_the_published_A": (
            "All six of that paper's activation energies reproduce from its own "
            "Table 1 to four decimal places and only TWO of its six "
            "pre-exponentials do -- a sign flip on every negative intercept and "
            "a SWAP of the Coffee-Cys/Coffee-Lys pair. It is the third audited "
            "case in the corpus of a correct Ea bolted to a defective "
            "prefactor. HMF-Cys happens to be one of the two that pass "
            "(published 23980.59 vs refit 24115.1, 0.56 % apart); the refit is "
            "used anyway, because using the published value for the rows that "
            "pass and the refit for the rows that fail would be a selection "
            "nobody could audit. Amendment 12 makes it binding."
        ),
        "rows": rows,
        "max_fold": max(r["fold"] for r in rows),
        "expected": (
            "NOT 1.00. The three-point refit has R^2 = 0.874 and its slope "
            "rests on the single 50 C point at which the authors declare the "
            "pseudo-first-order assumption compromised."
        ),
    }


def kocadagli_ingestion_check() -> Dict[str, Any]:
    """Every ingested Kocadagli constant, round-tripped 180 C -> 100 C -> 180 C."""
    from src.kinetic_core.parameters_furanic import MEASURED_FURANIC

    published_k_b_x1000 = {
        "k_glc_tdg": 4.19, "k_tdg_ddg": 30.5, "k_ddg_hmf": 119.0,
        "k_fru_int": 330.0, "k_int_hmf": 1.84, "k_fru_odg": 2.11,
        "k_tdg_mgo": 304.0,
    }
    rows = []
    for key, k_b in published_k_b_x1000.items():
        parameter = MEASURED_FURANIC[key]
        recovered = parameter.k_at(453.15) * 1000.0
        rows.append({
            "parameter": key,
            "published_k_b_min_inv_x1000": k_b,
            "recovered_k_b_min_inv_x1000": recovered,
            "relative_error": abs(recovered - k_b) / k_b,
            "ea_kj_mol": parameter.ea_kj_mol,
        })
    return {
        "what": (
            "each ingested Kocadagli glucose-system constant, re-referenced "
            "from the paper's T_b = 180 C to the core's 100 C at import and "
            "evaluated back at 180 C here. A pure arithmetic round trip, which "
            "is exactly why it is worth doing: the exponent's sign is the one "
            "place this ingestion could silently invert."
        ),
        "rows": rows,
        "max_relative_error": max(r["relative_error"] for r in rows),
    }


def af_sensitivity_sweep() -> Dict[str, Any]:
    """The pre-registered three-decade sweep on the one sourceless constant."""
    from dataclasses import replace

    from src.kinetic_core.engine import (
        SULFUR, FormulationSpec, ProcessSpec, ThermalProgram,
        core_parameters, predict,
    )
    from src.kinetic_core.parameters_furanic import (
        BLANK_AMINE_LOADING_MMOL_PER_L, BLANK_PH,
        BLANK_SUGAR_LOADING_MMOL_PER_L, BLANK_TEMPERATURE_C, BLANK_TIME_MIN,
        with_fitted_furanic,
    )
    from src.kinetic_core.ph_state import BufferSpec

    base = dict(core_parameters(SULFUR))
    from src.kinetic_core.parameters_furanic import FROZEN_K_DPO_AF_L_PER_MMOL_MIN
    base.update(with_fitted_furanic(FROZEN_K_DPO_AF_L_PER_MMOL_MIN))
    reference = base["k_af_dmhf"].k_ref

    spec = FormulationSpec(
        name="blank1997_sweep",
        precursors={
            "xylose": BLANK_SUGAR_LOADING_MMOL_PER_L,
            "glycine": BLANK_AMINE_LOADING_MMOL_PER_L,
        },
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(BLANK_TEMPERATURE_C, BLANK_TIME_MIN),
            ph=BLANK_PH,
            buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True,
                               source="Blank 1997 sec. 2: 0.2 M phosphate"),
        ),
    )
    rows = []
    for factor in (0.1, 1.0, 10.0, 100.0):
        parameters = dict(base)
        parameters["k_af_dmhf"] = replace(
            base["k_af_dmhf"], k_ref=reference * factor
        )
        run = predict(spec, ["DMHF"], parameters=parameters)
        rows.append({
            "factor_on_k_af_dmhf": factor,
            "dmhf_mmol_per_l": float(run.species_mmol_per_l.get("DMHF", 0.0)),
        })
    values = [r["dmhf_mmol_per_l"] for r in rows if r["dmhf_mmol_per_l"] > 0.0]
    span = (max(values) / min(values)) if len(values) > 1 else float("nan")
    return {
        "what": (
            "``k_af_dmhf`` (acetylformoin -> DMHF) has NO SOURCE OF ANY KIND. "
            "It encodes the assumption 'acetylformoin does not accumulate'. "
            "The pre-registration requires that sweeping it over three decades "
            "move the DMHF prediction by less than 1.2x; if it moves more, the "
            "constant is rate-limiting after all and the assumption is false."
        ),
        "rows": rows,
        "span_fold": span,
        "prereg_threshold_fold": 1.2,
        "passes": bool(span == span and span < 1.2),
    }


def branch_share_demonstration() -> Dict[str, Any]:
    """
    THE NO-FIXED-BRANCH-FRACTION DEMONSTRATION.

    Two charges, one model, two different HMF limb shares -- because the share
    is a consequence of the pools and not a stored number. K5a MUST-NOT #1.
    """
    from src.kinetic_core.engine import (
        TRUNK, FormulationSpec, ProcessSpec, ThermalProgram,
        core_parameters, predict,
    )
    from src.kinetic_core.furanic import hmf_limb_shares
    from src.kinetic_core.network import rate_constants_at

    parameters = core_parameters(TRUNK)
    rows = []
    for name, charge in (
        ("glucose only, 160 C", {"glucose": 300.0}),
        ("fructose only, 160 C", {"fructose": 300.0}),
        ("glucose + glycine, 160 C", {"glucose": 300.0, "glycine": 300.0}),
        ("glucose only, 120 C", {"glucose": 300.0}),
    ):
        temperature = 120.0 if "120" in name else 160.0
        spec = FormulationSpec(
            name=name,
            precursors=charge,
            process=ProcessSpec(
                thermal=ThermalProgram.isothermal(temperature, 30.0)
            ),
        )
        run = predict(spec, ["5-HMF"], parameters=parameters)
        k = rate_constants_at(parameters, temperature + 273.15)
        shares = hmf_limb_shares(run.species_mmol_per_l, k)
        rows.append({
            "system": name,
            "hmf_ug_per_l": run.concentrations_ug_per_l.get("5-HMF"),
            "fructose_limb_share": shares["fructose_limb"],
            "three_deoxyglucosone_limb_share": shares["three_deoxyglucosone_limb"],
        })
    shares = [r["fructose_limb_share"] for r in rows]
    return {
        "what": (
            "K5a MUST-NOT #1: the HMF node may not carry a fixed branch "
            "fraction. It does not -- the share below is computed from the "
            "Fru and 3,4-DG pools at the end of each run, and it MOVES with "
            "sugar identity, amine presence and temperature."
        ),
        "rows": rows,
        "share_range": [min(shares), max(shares)],
        "share_spread": max(shares) - min(shares),
        "moves": bool(max(shares) - min(shares) > 1.0e-6),
        "corroboration": (
            "Six papers, six matrices, verdicts spanning 'fructose limb "
            "dominant' to '3-DG limb dominant by infinity', and EVERY paper "
            "that names why its limb wins names a SUPPLY reason -- pool size, "
            "a starved 3-DG source, a drained cation pool -- never a terminal "
            "rate constant. In four of six published comparisons the terminal "
            "constant points the other way."
        ),
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def build_report() -> Dict[str, Any]:
    from src.kinetic_core.furanic import (
        DECLARED_GAPS, blank1997_conditions, blank1997_fit_cells,
        blank1997_structural_summary, describe_furanic,
    )
    from src.kinetic_core.parameters_furanic import furanic_registry_metadata

    fit = fit_k_dpo_af()
    return {
        "wave": "B7",
        "artifact": "kinetic_core_b7_fit_report",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "pre_registration": "results/validation/kinetic_core_b7_prereg.md",
        "declaration": (
            "docs/reference/FIT_HOLDOUT_DECLARATION.md Amendments 8 and 12"
        ),
        "fit": fit,
        "frozen_parameters": {"k_dpo_af": fit["value"]},
        "derived_parameters": {
            "k_odg_af": {
                "value": fit["value"] * 1000.0,
                "rule": (
                    "k_odg_af [1/min] = k_dpo_af [L/(mmol*min)] x 1000 mmol/L, "
                    "the pentose edge's pseudo-first-order constant at Blank's "
                    "1 M amine loading"
                ),
                "status": (
                    "DECLARED TRANSFER, and the weakest link in the DMHF node. "
                    "There is NO absolute hexose DMHF yield in any of the five "
                    "papers of the cluster; only the intact-C6 STRUCTURE is "
                    "measured, and it is measured twice."
                ),
            },
        },
        "fit_rows": [
            {
                "table": c.table, "sugar": c.sugar, "amine": c.amine,
                "compound": c.compound, "ug_per_mmol_sugar": c.ug_per_mmol_sugar,
                "role": c.role, "reason": c.reason,
            }
            for c in blank1997_fit_cells()
        ],
        "fit_conditions": blank1997_conditions(),
        "blank1997_structural_summary": blank1997_structural_summary(),
        "ingestion_check_hamzalioglu": hamzalioglu_ingestion_check(),
        "ingestion_check_kocadagli": kocadagli_ingestion_check(),
        "af_sensitivity_sweep": af_sensitivity_sweep(),
        "branch_share_demonstration": branch_share_demonstration(),
        "channel": describe_furanic(),
        "registry": furanic_registry_metadata(),
        "declared_gaps": list(DECLARED_GAPS),
    }


def render_markdown(payload: Dict[str, Any]) -> str:
    fit = payload["fit"]
    lines: List[str] = []
    add = lines.append
    add("# Build Wave B7 — the furanic channel — FIT REPORT")
    add("")
    add(f"Generated {payload['generated_on']} at `{payload['git']['short']}` "
        f"on `{payload['git']['branch']}`.")
    add(f"Pre-registered in `{payload['pre_registration']}`, written before this "
        f"script existed. Roles: {payload['declaration']}.")
    add("")
    add("## The whole fit, in one line")
    add("")
    add(f"**One free parameter.** `k_dpo_af` = **{fit['value']:.6e}** "
        f"{fit['unit']}, fitted to **three** cells of Blank 1997 Table 1.")
    add("")
    add("Everything else in the channel is ingested (Kocadagli's seven glucose "
        "constants, Hamzalioglu's sink), derived (the hexose transfer), a "
        "digitised prior (Edge B's level) or **exactly zero** (Edge C).")
    add("")
    add("### The two seeded starts")
    add("")
    add("| start `log10 k` | converged `log10 k` | objective |")
    add("|---|---:|---:|")
    for s in fit["starts"]:
        add(f"| {s['start_log10_k']:.1f} | {s['log10_k']:.9f} | "
            f"{s['objective']:.6g} |")
    add("")
    add(f"Agreement: **{fit['start_agreement_log10']:.2e}** decades "
        f"({'IDENTIFIED' if fit['starts_agree'] else 'NOT IDENTIFIED'}; "
        f"the pre-registration's threshold is 1e-6).")
    add("")
    add("### The three fit rows")
    add("")
    add("| row | measured µg/mmol | predicted µg/mmol | residual (decades) |")
    add("|---|---:|---:|---:|")
    for key, value in fit["targets_ug_per_mmol"].items():
        add(f"| {key} | {value:.2f} | {fit['predicted_ug_per_mmol']:.3f} | "
            f"{fit['residual_log10_fold'][key]:+.4f} |")
    add("")
    add(f"RMS residual **{fit['rms_residual_log10']:.4f} decades = "
        f"{fit['rms_residual_fold']:.3f}×**. The pre-registration predicted "
        f"≈0.12 decades ≈1.33× and said a materially SMALLER value would mean "
        f"something had been fitted beyond the one declared parameter.")
    add("")
    add("> " + fit["residual_is_the_sugar_axis"])
    add("")
    add("### The nine declared-FIT cells the core cannot represent")
    add("")
    add("| table | system | compound | µg/mmol | why not represented |")
    add("|---|---|---|---:|---|")
    for row in payload["fit_rows"]:
        if row["role"] != "fit_unrepresentable":
            continue
        add(f"| {row['table']} | {row['sugar']}/{row['amine']} | "
            f"{row['compound']} | {row['ug_per_mmol_sugar']:.1f} | "
            f"{row['reason']} |")
    add("")
    add("## The derived hexose transfer")
    add("")
    derived = payload["derived_parameters"]["k_odg_af"]
    add(f"`k_odg_af` = **{derived['value']:.6e}** 1/min. Rule: {derived['rule']}.")
    add("")
    add(f"> {derived['status']}")
    add("")
    add("## Ingestion check 1 — Hamzalioglu's HMF + cysteine sink")
    add("")
    ham = payload["ingestion_check_hamzalioglu"]
    add(ham["what"] + ".")
    add("")
    add("| T (°C) | measured (M⁻¹ day⁻¹) | measured (core units) | "
        "refit line (core units) | fold |")
    add("|---:|---:|---:|---:|---:|")
    for r in ham["rows"]:
        add(f"| {r['temperature_C']:.0f} | {r['measured_M_inv_day']:.2f} | "
            f"{r['measured_L_per_mmol_min']:.4e} | "
            f"{r['line_L_per_mmol_min']:.4e} | {r['fold']:.3f}× |")
    add("")
    add(f"Max fold **{ham['max_fold']:.3f}×**. {ham['expected']}")
    add("")
    add("**Why the refit prefactor and not the published one:** "
        + ham["why_the_refit_and_not_the_published_A"])
    add("")
    add("## Ingestion check 2 — Kocadagli's seven glucose constants")
    add("")
    koc = payload["ingestion_check_kocadagli"]
    add(koc["what"])
    add("")
    add("| parameter | published k_b (min⁻¹×10³) | recovered | rel. error | Ea |")
    add("|---|---:|---:|---:|---:|")
    for r in koc["rows"]:
        add(f"| `{r['parameter']}` | {r['published_k_b_min_inv_x1000']:.2f} | "
            f"{r['recovered_k_b_min_inv_x1000']:.4f} | "
            f"{r['relative_error']:.2e} | {r['ea_kj_mol']:.1f} |")
    add("")
    add(f"Max relative error **{koc['max_relative_error']:.2e}**.")
    add("")
    add("## The sourceless constant, swept")
    add("")
    sweep = payload["af_sensitivity_sweep"]
    add(sweep["what"])
    add("")
    add("| ×`k_af_dmhf` | DMHF (mmol/L) |")
    add("|---:|---:|")
    for r in sweep["rows"]:
        add(f"| {r['factor_on_k_af_dmhf']:g} | {r['dmhf_mmol_per_l']:.6e} |")
    add("")
    add(f"Span over three decades: **{sweep['span_fold']:.4f}×** against a "
        f"pre-registered threshold of {sweep['prereg_threshold_fold']}× — "
        f"**{'PASS' if sweep['passes'] else 'FAIL'}**.")
    add("")
    add("## No fixed branch fraction — demonstrated, not asserted")
    add("")
    branch = payload["branch_share_demonstration"]
    add(branch["what"])
    add("")
    add("| system | HMF (µg/L) | fructose-limb share | 3-DG-limb share |")
    add("|---|---:|---:|---:|")
    for r in branch["rows"]:
        hmf = r["hmf_ug_per_l"]
        add(f"| {r['system']} | {hmf:.4g} | "
            f"{r['fructose_limb_share']:.4f} | "
            f"{r['three_deoxyglucosone_limb_share']:.4f} |")
    add("")
    add(f"Share range **{branch['share_range'][0]:.4f} – "
        f"{branch['share_range'][1]:.4f}**; it MOVES: "
        f"**{branch['moves']}**.")
    add("")
    add("> " + branch["corroboration"])
    add("")
    add("## Declared gaps carried out of this wave")
    add("")
    for gap in payload["declared_gaps"]:
        add(f"* {gap}")
    add("")
    return "\n".join(lines) + "\n"


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))
    args = parser.parse_args(argv)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = build_report()
    (out / f"{BASENAME}.json").write_text(json.dumps(payload, indent=2) + "\n")
    (out / f"{BASENAME}.md").write_text(render_markdown(payload))
    print(f"k_dpo_af = {payload['fit']['value']:.10e} L/(mmol*min)")
    print(f"wrote {out / BASENAME}.json and .md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
