#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b6_fit.py

THE FIT REPORT OF BUILD WAVE B6 (the lipid-oxidation module).

WHAT IS FITTED
--------------
Four 3-way branch simplexes (13-OOH and 9-OOH, x cis,trans and trans,trans
geometry) and four hydroperoxide-composition parameters, to EIGHTEEN numbers:
Frankel 1989's three ZERO-ADDITIVE columns. That is the whole fit column
declared for Module 5 in ``docs/reference/FIT_HOLDOUT_DECLARATION.md`` D.6.

WHAT IS NOT FITTED, AND WHY
---------------------------
  * The RATE. No rate constant for linoleate hydroperoxide -> hexanal at
    cooking temperature exists anywhere (declared gap F.3, re-affirming k3
    sec. C.9). The rate is a bounded INPUT with a labelled Q10 assumption.
  * The Q10. It is the REPO'S assumption, exposed with a band.
  * Any activation energy. Frankel is ONE TEMPERATURE (180 C).
  * The hydrogen-donor suppression ``d``. It has no stored value at all: every
    donor claim this module makes is a monotonicity theorem over ``d in (0,1)``.
  * Anything at all against the alpha-tocopherol arms, the 1,4-cyclohexadiene
    arms, the nonanal absence, or the two frozen Bi 2020 bundles.

Writes results/validation/kinetic_core_b6_fit_report.{json,md}.
Run with: ./scripts/docker_maillard.sh run \
  "python scripts/generators/generate_kinetic_core_b6_fit.py"
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.kinetic_core.lipid import (  # noqa: E402
    REFERENCE_AUTOXIDATION_SYSTEM,
    Y_HEXANAL_PER_LOOH,
    describe_lipid,
    fit_branch_model,
    slate_yields,
    validate_lipid_structure,
)
from src.kinetic_core.parameters_lipid import (  # noqa: E402
    FRANKEL_ANCHOR,
    FRANKEL_DOSSIER,
    FRANKEL_ZERO_ADDITIVE,
    FRANKEL_ZERO_ADDITIVE_TOTAL_AREA,
    LIPID_SOURCE_CONTRADICTIONS,
    PROHIBITED_DERIVATIONS,
    SHIPPED_VALUES_REFUTED,
    lipid_registry_metadata,
)
from src.kinetic_core.species_lipid import (  # noqa: E402
    FRANKEL_SLATE,
    NAMED_UNQUANTIFIED_COPRODUCTS,
)

OUT_JSON = REPO / "results/validation/kinetic_core_b6_fit_report.json"
OUT_MD = REPO / "results/validation/kinetic_core_b6_fit_report.md"
PREREG = REPO / "results/validation/kinetic_core_b6_prereg.md"

#: The pre-registered fit tolerances, quoted from the prereg so the report
#: cannot grade itself leniently.
PREREG_MEDIAN_PP = 3.0
PREREG_WORST_PP = 8.0


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=REPO, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"

    return {"commit": _run("rev-parse", "HEAD"),
            "branch": _run("rev-parse", "--abbrev-ref", "HEAD")}


def main() -> int:
    if not PREREG.exists():
        raise SystemExit(
            f"{PREREG} is missing. The pre-registration is written BEFORE the "
            f"fit, and this script refuses to run without it."
        )

    fit = fit_branch_model()
    branch = fit["branch_model"]
    compositions = fit["compositions"]

    # --- the predicted FIT table, side by side with the measured one -------
    table: List[Dict[str, Any]] = []
    for system, observed in FRANKEL_ZERO_ADDITIVE.items():
        predicted = branch.slate_shares(compositions[system], donor=0.0)
        scale = sum(observed.values())
        for product in FRANKEL_SLATE:
            table.append({
                "system": system,
                "product": product,
                "measured_relative_percent": float(observed[product]),
                "predicted_relative_percent": predicted[product] * scale,
                "residual_pp": predicted[product] * scale - float(observed[product]),
            })

    # --- F-2: the shipped 0.37, tested like for like ----------------------
    reference_composition = compositions[REFERENCE_AUTOXIDATION_SYSTEM]
    hexanal_slate_share = branch.slate_shares(reference_composition)["HEXANAL"]
    measured_range = [
        min(v["HEXANAL"] for v in FRANKEL_ZERO_ADDITIVE.values()),
        max(v["HEXANAL"] for v in FRANKEL_ZERO_ADDITIVE.values()),
    ]

    # --- the 8.5x pairing anomaly, computed rather than asserted ----------
    pairing = {
        system: float(v["PENTANE"]) / float(v["ME_13_OXO_TRIDECADIENOATE"])
        for system, v in FRANKEL_ZERO_ADDITIVE.items()
    }
    pairing_spread = max(pairing.values()) / min(pairing.values())

    payload: Dict[str, Any] = {
        "artifact": "kinetic_core_b6_fit_report",
        "wave": "B6 -- the lipid-oxidation module",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "pre_registration": str(PREREG.relative_to(REPO)),
        "fit_column": {
            "source": FRANKEL_ANCHOR,
            "dossier": FRANKEL_DOSSIER,
            "systems": list(FRANKEL_ZERO_ADDITIVE),
            "values": {k: dict(v) for k, v in FRANKEL_ZERO_ADDITIVE.items()},
            "n_values": 18,
            "what_was_read": (
                "the three ZERO-ADDITIVE columns and nothing else. "
                "fit_residuals() asserts the array shape before it computes."
            ),
        },
        "frozen_parameters": {
            "branch_model": branch.as_dict(),
            "system_compositions": {
                system: composition.as_dict()
                for system, composition in compositions.items()
            },
            "default_pool_composition": compositions[
                REFERENCE_AUTOXIDATION_SYSTEM
            ].as_dict(),
            "default_pool_composition_provenance": (
                "Frankel 1989 Table 1 is the AUTOXIDATION mixture, so its "
                "fitted composition is the closest the corpus comes to 'what an "
                "oxidising food lipid's hydroperoxide pool looks like'. It is a "
                "FIT quantity, not an assumption."
            ),
            "Y_hexanal_per_LOOH": Y_HEXANAL_PER_LOOH,
            "slate_yields_at_reference": slate_yields(branch, reference_composition),
        },
        "fit_quality": {
            "median_abs_residual_pp": fit["median_abs_residual_pp"],
            "worst_abs_residual_pp": fit["worst_abs_residual_pp"],
            "sum_squared_residuals": fit["sum_squared_residuals"],
            "n_fit_values": fit["n_fit_values"],
            "n_free_parameters": fit["n_free_parameters"],
            "degrees_of_freedom": fit["degrees_of_freedom"],
            "prereg_median_bound_pp": PREREG_MEDIAN_PP,
            "prereg_worst_bound_pp": PREREG_WORST_PP,
            "prereg_F1_met": bool(
                fit["median_abs_residual_pp"] <= PREREG_MEDIAN_PP
                and fit["worst_abs_residual_pp"] <= PREREG_WORST_PP
            ),
        },
        "fit_table": table,
        "prereg_F2_shipped_0_37_refuted": {
            "fitted_hexanal_share_of_six_product_slate": hexanal_slate_share,
            "frankel_measured_zero_additive_range_percent": measured_range,
            "shipped_fast_lane_value": 0.37,
            "verdict": (
                "REFUTED"
                if hexanal_slate_share < 0.37 else "NOT REFUTED"
            ),
            "note": (
                "Compared LIKE FOR LIKE: both are a hexanal share of the "
                "six-product slate. The per-isomer simplex values are NOT "
                "comparable to 0.37 -- they are shares of a three-product "
                "sub-slate -- and quoting one of those instead would be the "
                "kind of denominator swap this repository has already been "
                "caught by."
            ),
        },
        "prereg_F3_nonanal_refuted_structurally": validate_lipid_structure(branch),
        "pentane_me13oxo_pairing_anomaly": {
            "ratio_by_system": pairing,
            "spread_x": pairing_spread,
            "expected_if_1to1": 1.0,
            "finding": LIPID_SOURCE_CONTRADICTIONS["pentane_vs_me13oxo_pairing"],
        },
        "total_peak_areas_NOT_fitted": {
            "values": dict(FRANKEL_ZERO_ADDITIVE_TOTAL_AREA),
            "why_not": (
                "relative to an internal standard at ONE temperature, and the "
                "three preparations were injected at different loadings "
                "(13.1 mg for the mixture; Table 2's totals are ~7x smaller). "
                "They carry no absolute yield and are not fitted."
            ),
        },
        "registry": lipid_registry_metadata(),
        "network": describe_lipid(),
        "named_unquantified_coproducts": dict(NAMED_UNQUANTIFIED_COPRODUCTS),
        "prohibited_derivations": dict(PROHIBITED_DERIVATIONS),
        "source_contradictions": dict(LIPID_SOURCE_CONTRADICTIONS),
        "shipped_values_refuted": dict(SHIPPED_VALUES_REFUTED),
        "holdout_exposure_disclosure": (
            "THE BUILDER OF THIS WAVE SAW THE ALPHA-TOCOPHEROL HOLD-OUT "
            "COLUMNS. Frankel 1989 prints the zero-additive column (FIT) and "
            "the tocopherol columns (HOLD-OUT) in the SAME TABLE ROWS, and the "
            "abstract states the hold-out result in prose. Reading the fit "
            "column without seeing the hold-out column is impossible from that "
            "PDF. Consequence, adopted in the pre-registration sec. 0 and "
            "binding: the tocopherol hold-out is scored seen_diagnostic, never "
            "gating (the Amendment 9 clause 1 precedent), and the mitigation is "
            "structural -- the donor term has NO FITTED PARAMETER, so there is "
            "nothing to tune toward the seen numbers. The nonanal absence stays "
            "GATING because it is a structural zero fixed by molecular "
            "topology. The two frozen Bi 2020 bundles were never opened."
        ),
    }

    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n")

    lines: List[str] = []
    lines.append("# Kinetic core B6 -- FIT report (the lipid-oxidation module)\n")
    lines.append(f"Generated {payload['generated_on']} at "
                 f"`{payload['git']['commit'][:12]}` "
                 f"(branch `{payload['git']['branch']}`).\n")
    lines.append(f"Pre-registration: `{payload['pre_registration']}` -- "
                 f"written before this file existed.\n")

    lines.append("## The one thing to read first\n")
    lines.append(
        "**The branch DISTRIBUTION is measured. The absolute RATE is not.** "
        "The rate enters as a bounded input anchored at Schroen & "
        "Berton-Carabin 2022's 25 C, hand-fitted, no-standard-error, lumped "
        "`k4 = 6e-3 /h`, with a Q10 that is the repo's own assumption carried "
        "with a [2, 3] band and a mandatory warning on every prediction. No "
        "Q10 number is baked into any stored constant.\n")

    lines.append("## Hold-out exposure disclosure\n")
    lines.append(payload["holdout_exposure_disclosure"] + "\n")

    lines.append("## The fit\n")
    quality = payload["fit_quality"]
    lines.append(
        f"| statistic | value | pre-registered bound |\n|---|---:|---:|\n"
        f"| median abs residual (percentage points) | "
        f"{quality['median_abs_residual_pp']:.2f} | {PREREG_MEDIAN_PP:.1f} |\n"
        f"| worst abs residual (percentage points) | "
        f"{quality['worst_abs_residual_pp']:.2f} | {PREREG_WORST_PP:.1f} |\n"
        f"| fit values / free parameters / df | "
        f"{quality['n_fit_values']} / {quality['n_free_parameters']} / "
        f"{quality['degrees_of_freedom']} | -- |\n"
        f"| **F-1 met?** | **{'YES' if quality['prereg_F1_met'] else 'NO'}** | -- |\n")

    lines.append("\n### Measured vs predicted, all eighteen FIT numbers\n")
    lines.append("| system | product | measured % | predicted % | residual pp |")
    lines.append("|---|---|---:|---:|---:|")
    for row in table:
        lines.append(
            f"| {row['system']} | {row['product']} | "
            f"{row['measured_relative_percent']:.1f} | "
            f"{row['predicted_relative_percent']:.1f} | "
            f"{row['residual_pp']:+.2f} |"
        )

    lines.append("\n### Frozen branch simplexes\n")
    lines.append("| hydroperoxide | product | share |")
    lines.append("|---|---|---:|")
    for cell, shares in sorted(payload["frozen_parameters"]["branch_model"]["simplexes"].items()):
        for product, share in shares.items():
            lines.append(f"| {cell} | {product} | {share:.4f} |")

    lines.append("\n### Fitted hydroperoxide-pool compositions\n")
    lines.append("| system | 13-ct | 13-tt | 9-ct | 9-tt |")
    lines.append("|---|---:|---:|---:|---:|")
    for system, cells in payload["frozen_parameters"]["system_compositions"].items():
        lines.append(
            f"| {system} | {cells['LOOH_13_ct']:.3f} | {cells['LOOH_13_tt']:.3f} "
            f"| {cells['LOOH_9_ct']:.3f} | {cells['LOOH_9_tt']:.3f} |"
        )

    f2 = payload["prereg_F2_shipped_0_37_refuted"]
    lines.append("\n## F-2 -- the shipped `hexanal 0.37`\n")
    lines.append(
        f"Fitted hexanal share of the six-product slate at the autoxidation "
        f"reference composition: **{f2['fitted_hexanal_share_of_six_product_slate']:.4f}**. "
        f"Frankel's own zero-additive range: "
        f"{f2['frankel_measured_zero_additive_range_percent'][0]:.0f}-"
        f"{f2['frankel_measured_zero_additive_range_percent'][1]:.0f} %. "
        f"Shipped FAST-lane value: 0.37. **{f2['verdict']}.**\n")
    lines.append(f2["note"] + "\n")

    lines.append("## F-3 -- the shipped `nonanal 0.15`\n")
    lines.append(payload["prereg_F3_nonanal_refuted_structurally"]["nonanal"] + "\n")

    anomaly = payload["pentane_me13oxo_pairing_anomaly"]
    lines.append("## The falsified 1:1 pairing (pre-registered as F-1's reason)\n")
    lines.append(
        "| system | pentane : methyl 13-oxo-tridecadienoate |\n|---|---:|\n"
        + "\n".join(f"| {k} | {v:.3f} |" for k, v in anomaly["ratio_by_system"].items())
        + f"\n\nSpread: **{anomaly['spread_x']:.1f}x** against an expected 1.0. "
        + anomaly["finding"] + "\n")

    lines.append("## What was refused\n")
    for what, why in payload["prohibited_derivations"].items():
        lines.append(f"* **{what}** -- {why}")

    lines.append("\n## Contradictions reported, not resolved\n")
    for what, why in payload["source_contradictions"].items():
        lines.append(f"* **{what}** -- {why}")

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    print(f"wrote {OUT_MD.relative_to(REPO)}")
    print(f"median residual {quality['median_abs_residual_pp']:.2f} pp, "
          f"worst {quality['worst_abs_residual_pp']:.2f} pp, "
          f"F-1 {'MET' if quality['prereg_F1_met'] else 'MISSED'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
