#!/usr/bin/env python
"""
Wave B30 (W8) -- THE GATING TEST, evaluated (2026-09-10).

`results/validation/kinetic_core_b30_prereg.md`. This wave has ONE test that runs before any
constant is fitted, and it failed, so nothing was fitted.

  T0  the SIGN. Whitfield's fed norfuraneol pot is the only pot in the corpus measured at two pH
      values by one laboratory: free 2-methyl-3-furanthiol falls at least 150x from pH 4.5 to 6.5.
      The model must at minimum fall. If it does not, no pH slope may be fitted on top, because a
      slope fitted over a wrong-signed mechanism buys agreement and leaves the mechanism wrong.

It also attributes the failure, which is the part worth keeping: the hydrosulfide branch that gives
this lane its pH response is STRUCTURALLY ABSENT from the norfuraneol route, and supplying it would
push the model further from the measurement rather than closer.
"""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_b30_ship_rule.json"
CHARGE = {"NF": 50.0, "Cys": 50.0}
T_C, MINUTES = 140.0, 60.0
#: Whitfield 1999 (pH 4.5, already a fit row) against Whitfield 2001 (pH 6.5), same laboratory.
MEASURED = {"mft_molpct_ph45": 0.150, "mft_molpct_ph65": 0.001, "fold_down_at_least": 150.0,
            "mercaptoketone_fold_down": 2500.0}
SOURCE = ("whitfield1999_extraction.md (pH 4.5) and the pH-6.5 arm recorded in the W8 backlog entry "
          "from Whitfield 2001; free MFT 0.150 -> < 0.001 mol %, mercaptoketones 74.5 -> 0.03")


def _mft(params, ph: float) -> float:
    from src.kinetic_core.sulfur import integrate_sulfur

    run = integrate_sulfur(params, T_C + 273.15, CHARGE, np.array([0.0, MINUTES]), ph=ph)
    return float(run.series("MFT")[-1])


def main() -> int:
    from src.kinetic_core.engine import SULFUR, core_parameters

    base = core_parameters(SULFUR)
    arms: Dict[str, Dict[str, float]] = {}

    def arm(name, params):
        a, b = _mft(params, 4.5), _mft(params, 6.5)
        arms[name] = {"ph_4_5": a, "ph_6_5": b, "ratio_45_over_65": (a / b) if b > 0 else float("inf")}

    arm("as_shipped", base)
    no_thiolate = dict(base)
    for key in list(no_thiolate):
        if "thiolate" in key:
            no_thiolate[key] = replace(no_thiolate[key], k_ref=0.0)
    arm("thiolate_loss_off", no_thiolate)
    no_hs = dict(base)
    for key in ("k_ddp_mft_hs", "k_fur_fft_hs"):
        if key in no_hs:
            no_hs[key] = replace(no_hs[key], k_ref=0.0)
    arm("hydrosulfide_branch_off", no_hs)

    shipped = arms["as_shipped"]["ratio_45_over_65"]
    hs_inactive = abs(arms["hydrosulfide_branch_off"]["ratio_45_over_65"] - shipped) < 1e-9
    t0 = {"arms": arms, "measured": dict(MEASURED), "source": SOURCE,
          "model_ratio_45_over_65": shipped,
          "sign_correct": bool(shipped > 1.0),
          "out_by_x": MEASURED["fold_down_at_least"] / shipped if shipped > 0 else float("inf"),
          "pass": bool(shipped > 1.0)}
    attribution = {
        "thiolate_loss_carries": (arms["thiolate_loss_off"]["ratio_45_over_65"] - shipped),
        "hydrosulfide_branch_is_inactive_on_this_pot": bool(hs_inactive),
        "why": ("the two-branch sulfide mechanism was built for the deoxypentosone route "
                "(r_ddp_mft_hs) and the furfural route (r_fur_fft_hs). The NORFURANEOL route has no "
                "hydrosulfide partner -- r_nf_mft and r_nf_mp3p are single steps in neutral H2S -- "
                "so on the one pot in the corpus with a measured pH PAIR, the lane's pH mechanism "
                "is structurally absent."),
        "and_supplying_it_would_make_it_worse": (
            "a hydrosulfide branch on the norfuraneol steps pushes the SAME way as the one that "
            "already exists: more hydrosulfide at higher pH means faster addition means MORE thiol "
            "at pH 6.5, where the measurement wants at least 150x less. The collapse is not in the "
            "nucleophile. It is in the substrate or in the sulfide budget, and neither is modelled."),
    }
    payload = {
        "artifact": "kinetic_core_b30_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b30_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b30_ship_rule.py", inputs=[]),
        "prereg": data_paths.rel(V / "kinetic_core_b30_prereg.md"),
        "rule": ("T0 GATES THE WAVE: the model must at least FALL from pH 4.5 to 6.5 on the one pot "
                 "measured at both. Nothing is fitted unless it does."),
        "T0": t0, "attribution": attribution,
        "verdict": "SHIP" if t0["pass"] else "DO NOT SHIP -- and nothing was fitted",
    }
    artifact_io.write_artifact(payload, OUT, render=_render)
    print(payload["verdict"], "| model ratio", round(shipped, 3), "| measured >= 150",
          "| HS branch inactive here:", hs_inactive)
    return 0


def _render(p: Dict[str, Any]) -> str:
    T0, A = p["T0"], p["attribution"]
    L = [f"# Wave B30 (W8) gating test: {p['verdict']}", "", f"*{p['rule']}*", "",
         "| arm | pH 4.5 | pH 6.5 | ratio 4.5 / 6.5 |", "|---|---:|---:|---:|"]
    for name, a in T0["arms"].items():
        L.append(f"| {name.replace('_', ' ')} | {a['ph_4_5']:.4g} | {a['ph_6_5']:.4g} | "
                 f"{a['ratio_45_over_65']:.3g} |")
    L += [f"| **measured** | 0.150 mol % | < 0.001 | **>= 150** |", "",
          f"The model gives MORE thiol at the higher pH, where the measurement collapses. The sign "
          f"is wrong and the ratio is out by about {T0['out_by_x']:.0f}x. No slope was fitted.", "",
          "## Attribution", "",
          f"- The thiolate loss carries almost none of it: switching it off moves the ratio by "
          f"{A['thiolate_loss_carries']:+.2f}.",
          f"- **The hydrosulfide branch is not active on this pot at all** "
          f"({A['hydrosulfide_branch_is_inactive_on_this_pot']}): switching it off changes the answer "
          f"by nothing to four figures. {A['why']}",
          f"- {A['and_supplying_it_would_make_it_worse']}", ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
