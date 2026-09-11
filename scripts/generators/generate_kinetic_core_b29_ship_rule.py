#!/usr/bin/env python
"""
Wave B29 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-10).

`results/validation/kinetic_core_b29_prereg.md` section 4.

  T1  the two air/argon ratios within 0.2 dex on a SINGLE f(argon) -- decisive, and the whole point:
      one multiplier has to explain 9.2 and 3.5
  T2  the copper arm, likewise
  T3  air is EXACTLY 1: a pot that declares nothing is bit-for-bit what it was before this axis
  T4  identification
Ship rule: SHIP if T1, T2, T3 and T4 hold.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
REPORT = V / "kinetic_core_b29_fit_report.json"
OUT = V / "kinetic_core_b29_ship_rule.json"


def _oxidative_share() -> Dict[str, float]:
    """How much of the Strecker aldehyde comes through the oxidative entries, per pot."""
    import numpy as np
    from types import SimpleNamespace

    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    import generate_kinetic_core_b29_fit as G
    from src.kinetic_core import trunk_conditions
    from src.kinetic_core.integrate import integrate

    out = {}
    for tag, init in (("fed_amadori", {"AMA": G.CHARGE_MMOL_L, "Gly": G.CHARGE_MMOL_L}),
                      ("glucose_glycine", {"Glc": G.CHARGE_MMOL_L, "Gly": G.CHARGE_MMOL_L})):
        vals = []
        for f in (1.0, 1e-12):
            p, _ = trunk_conditions.apply(
                dict(G.BASE), SimpleNamespace(ph=G.PH, water_activity=None, atmosphere="argon"),
                atmosphere_factors={"air": 1.0, "argon": f})
            r = integrate(p, G.T_C + 273.15, init, np.array([0.0, G.MINUTES]), rtol=1e-8, atol=1e-16)
            vals.append(sum(float(r.series(k)[-1]) for k in G.OBSERVABLE))
        out[tag] = 1.0 - (vals[1] / vals[0] if vals[0] > 0 else 0.0)
    return out


def main() -> int:
    rep = json.loads(REPORT.read_text(encoding="utf-8"))
    resid = rep["residual_by_row_dex"]
    pred = rep["predicted"]
    t1 = {"rows": {k: {"printed": rep["objective"]["targets"][k], "model": pred[k], "dex": resid[k]}
                   for k in ("arp_air_over_argon", "glc_air_over_argon")},
          "model_ratio_of_the_two_pots": pred["arp_air_over_argon"] / pred["glc_air_over_argon"],
          "printed_ratio_of_the_two_pots": 9.2 / 3.5,
          "pass": bool(max(abs(resid["arp_air_over_argon"]), abs(resid["glc_air_over_argon"])) <= 0.2)}
    t2 = {"rows": {k: {"printed": rep["objective"]["targets"][k], "model": pred[k], "dex": resid[k]}
                   for k in ("arp_aircu_over_air", "glc_aircu_over_air")},
          "pass": bool(max(abs(resid["arp_aircu_over_air"]), abs(resid["glc_aircu_over_air"])) <= 0.2)}
    t3 = dict(rep["air_is_exactly_one"])
    lap = rep["laplace"]
    t4 = {"sigma": lap["sigma"], "identified": lap["identified"], "on_bound": lap["on_bound"],
          "pass": bool(all(lap["identified"].values()) and not any(lap["on_bound"].values()))}
    share = _oxidative_share()
    diagnosis = {
        "oxidative_share_of_the_strecker_aldehyde": share,
        "what_it_means": (
            "CORRECTED ON REVIEW, 2026-09-10. The first run measured AKG alone -- glyoxal's Strecker "
            "product -- found it 100 % oxidative in both pots, and concluded that a non-oxidative "
            "route to the Strecker aldehyde was missing. That was an artefact of the observable: the "
            "trunk HAS a non-oxidative route (Amadori -> 1-deoxyosone -> methylglyoxal -> AKM), and "
            "summed over both Strecker products the non-oxidative share is 45 % in the Amadori pot "
            "and 36 % in the sugar pot. THE REAL FINDING IS THE ORDER. Hofmann's Amadori pot is the "
            "MORE oxygen-sensitive (9.2x against 3.5x), so its oxidative share must be the larger. "
            "The model's is the SMALLER. That is a statement about the trunk's branching between "
            "the oxidative route to glucosone and the non-oxidative routes to the deoxyosones, and "
            "it is the reason one multiplier cannot serve both pots."),
    }
    ships = bool(t1["pass"] and t2["pass"] and t3["pass"] and t4["pass"])
    payload = {
        "artifact": "kinetic_core_b29_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b29_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b29_ship_rule.py", inputs=[REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b29_prereg.md"),
        "rule": "SHIP if T1 (both air/argon ratios on one multiplier), T2, T3 (air exactly 1) and T4 hold",
        "T1": t1, "T2": t2, "T3": t3, "T4": t4, "diagnosis": diagnosis,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=_render)
    print(payload["verdict"], "| T1", t1["pass"], "| T2", t2["pass"], "| T3", t3["pass"],
          "| T4", t4["pass"], "| oxidative share", {k: f"{v:.1%}" for k, v in share.items()})
    return 0


def _render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4, D = p["T1"], p["T2"], p["T3"], p["T4"], p["diagnosis"]
    L = [f"# Wave B29 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}.*", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 the two air/argon ratios | model separates the pots by "
         f"{T1['model_ratio_of_the_two_pots']:.2f}x against a printed "
         f"{T1['printed_ratio_of_the_two_pots']:.2f}x | {T1['pass']} |",
         f"| T2 the copper arm | worst "
         f"{max(abs(r['dex']) for r in T2['rows'].values()):+.2f} dex | {T2['pass']} |",
         f"| T3 air is exactly 1 | parameters identical {T3['parameters_identical']}, observable "
         f"identical {T3['observable_identical']} | {T3['pass']} |",
         f"| T4 identification | sigma {json.dumps({k: round(v, 2) for k, v in T4['sigma'].items()})} | {T4['pass']} |",
         "", "| ratio | printed | model | dex |", "|---|---:|---:|---:|"]
    for block in (T1["rows"], T2["rows"]):
        for k, r in block.items():
            L.append(f"| {k} | {r['printed']:.3g} | {r['model']:.4g} | {r['dex']:+.3f} |")
    L += ["", "## What the axis exposed", "",
          "| pot | share of the Strecker aldehyde made through the oxidative entries |", "|---|---:|"]
    for k, v in D["oxidative_share_of_the_strecker_aldehyde"].items():
        L.append(f"| {k} | {v:.1%} |")
    L += ["", f"> {D['what_it_means']}", ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
