#!/usr/bin/env python
"""
probe_pentose_amine_entry.py -- the pentose + amine entry's cheap probe (kinetic_core_pentose_amine_entry_prereg.md
sec. 3, step 1; 2026-10-30). Nothing is fitted.

On the joint thiol fit's s042 engine state, with PENTOSE_AMINE_ENTRY set on a 3 x 3 grid of the entry rate and the
acetic-acid rate (both pH slopes at a declared 0.3), score every pot that the entry or the acetic step touches:
Davidek 2008 (xylose + glycine, acetic acid and residual xylose, pH 6 / 8), Desclaux 2006 (xylose + glycine,
norfuraneol and 3-deoxypentosone), Davidek 2006a glucose + glycine acetic acid (the hexose side: must not move), and
Zhou 2023's final pH (ARP-fed: the acetic step, not the entry). Also: ARP's derived enolisation rates at Desclaux's
conditions against Desclaux's measured glycine-xylose constants (the prereg's declared-assumption check).

    python scripts/generators/probe_pentose_amine_entry.py
Writes results/validation/kinetic_core_pentose_amine_entry_probe.{json,md}.
"""
from __future__ import annotations

import dataclasses
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, fit_rows, parameters_sulfur, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402

BASE_FIT = data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit_s042.json"
OUT = data_paths.VALIDATION_DIR / "kinetic_core_pentose_amine_entry_probe"
#: the grid (prereg sec. 3): Davidek's xylose half-life implies ~2e-4 L/(mmol min) at 100 C / pH 6 (6e-4 at 120 C)
#: second pass (after the first found acetic acid still low at 3e-2 and xylose overshot at 6e-4): wider in acetic acid
K_ENTRY = (3e-5, 7e-5, 2e-4)
K_DPO_AA = (1e-2, 3e-2, 1e-1, 3e-1)
SLOPE = 0.3
#: (table, pot-id prefix) of every pot scored
POTS = (("ph_acid_davidek_acids", "davidek"), ("deoxyosone_submodel", "desclaux2006"), ("ph_acid_submodel", "zhou2023"))
#: Desclaux 2006 (plan sec. 7): ARP -> 1-DP 2.79e-2 /min, -> 3-DP 5.50e-4 /min, xylose + glycine, pH 6, 100 C
DESCLAUX = {"r_pama_dpo": 2.79e-2, "r_pama_tdp": 5.50e-4}


FACTORS = {"pama_dpo_factor": 1.0, "pama_tdp_factor": 1.0}


def _install(job: Dict[str, Any], k_entry: float, k_aa: float) -> None:
    joint.sm3.install(job)
    parameters_sulfur.STRECKER_SULFUR = dataclasses.replace(parameters_sulfur.STRECKER_SULFUR,
                                                            h2s_share=float(job["strecker_share"]))
    parameters_sulfur.PENTOSE_AMINE_ENTRY = parameters_sulfur.PentoseAmineEntry(
        k_entry_l_per_mmol_min_100c_ph6=k_entry, entry_ph_slope=SLOPE, k_dpo_aa_per_min_100c_ph6=k_aa,
        aa_ph_slope=SLOPE, **FACTORS)


def main() -> int:
    fit = json.loads(BASE_FIT.read_text())
    theta = fit["fit"]["theta"]
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    sulfur.PH_COUPLING = "inline"
    jobs: List[tuple] = []
    for table, prefix in POTS:
        pots, rows = fit_rows.load(table)
        for j in probe.plan(pots, rows):
            if j["pot"].startswith(prefix) and "skip" not in j:
                jobs.append((table, j, [r for r in rows if r.pot == j["pot"]]))

    # the declared-assumption check first (prereg sec. 2): PAMA at Desclaux's measured level when ARP's is > x3 off
    _install({"mgo_share": trial["mgo_share"], "routes": trial["routes"], "strecker_share": fit["strecker_share"]},
             1e-4, 0.0)
    k100 = sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), 373.15, 6.0, ph_label=6.0)
    check = {k: {"lane": float(k100[k]), "desclaux": v, "log10_ratio": math.log10(float(k100[k]) / v)}
             for k, v in DESCLAUX.items()}
    for k, name in (("r_pama_dpo", "pama_dpo_factor"), ("r_pama_tdp", "pama_tdp_factor")):
        if abs(check[k]["log10_ratio"]) > math.log10(3.0):
            FACTORS[name] = check[k]["desclaux"] / check[k]["lane"]

    grid: Dict[str, Any] = {}
    for k_entry in (0.0,) + K_ENTRY:
        for k_aa in ((0.0,) if k_entry == 0.0 else K_DPO_AA):
            label = "off" if k_entry == 0.0 else f"entry {k_entry:g} / aa {k_aa:g}"
            cells: Dict[str, Dict[str, List[float]]] = {}
            for table, j, rows in jobs:
                job = dict(j, **trial, strecker_share=float(fit["strecker_share"]))
                _install(job, k_entry, k_aa)
                result = probe.run_pot_limited(job, joint.POT_SECONDS)
                for r in rows:
                    if r.species not in probe.SULFUR_STATE_KEYS and r.species != "pH":
                        continue
                    p = probe.predict(r, job, result)
                    d = (p - r.value) if r.species == "pH" else (math.log10(p / r.value) if p > 0 else float("nan"))
                    cells.setdefault(f"{r.ledger_id} | {r.species}", {}).setdefault("res", []).append(d)
            grid[label] = {k: {"n": len(v["res"]), "mean": float(sum(v["res"]) / len(v["res"])),
                               "rms": float(math.sqrt(sum(x * x for x in v["res"]) / len(v["res"])))}
                           for k, v in cells.items()}
            print(f"  {label}: done", flush=True)

    payload = {"artifact": OUT.name,
               "provenance": provenance.provenance_block(
                   OUT.name, generated_by="scripts/generators/probe_pentose_amine_entry.py",
                   inputs=[BASE_FIT, ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py"]
                   + [fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t, _ in POTS]),
               "base_fit": data_paths.rel(BASE_FIT), "ph_slopes": SLOPE, "pama_factors": FACTORS, "grid": grid,
               "desclaux_constant_check": check}
    OUT.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    series = sorted({k for g in grid.values() for k in g})
    lines = ["# Pentose + amine entry: cheap probe (nothing fitted)", "",
             f"*Generated by `scripts/generators/probe_pentose_amine_entry.py` on `{data_paths.rel(BASE_FIT)}`'s engine "
             f"state; both pH slopes {SLOPE}. Residual = predicted minus measured: dex, or pH units for pH rows. "
             "Rates at 100 C / pH 6: entry L/(mmol min), acetic 1/min.*", "",
             "## ARP's enolisation rates against Desclaux 2006's measured glycine-xylose constants (100 C, pH 6)", ""]
    lines += [f"- {k}: lane {v['lane']:.3g} /min, Desclaux {v['desclaux']:.3g} /min, log10 ratio {v['log10_ratio']:+.2f}"
              for k, v in check.items()]
    lines += [f"- applied (ARP's > x3 off): " + ", ".join(f"{k} {v:.3g}" for k, v in FACTORS.items())]
    lines += ["", "## Mean / RMS per data set x observable", "",
              "| setting | " + " | ".join(series) + " |", "|---" * (len(series) + 1) + "|"]
    for label, g in grid.items():
        lines.append(f"| {label} | " + " | ".join(
            f"{g[s]['mean']:+.2f} / {g[s]['rms']:.2f}" if s in g else "" for s in series) + " |")
    OUT.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(OUT.with_suffix('.md'))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
