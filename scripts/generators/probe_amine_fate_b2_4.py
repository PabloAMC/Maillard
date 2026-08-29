"""
B2.4 -- THE AMINE-FATE PROBE, REBUILT AGAINST THE CURRENT SPECIES SET.

WHY THIS EXISTS
---------------
`ph_state.AMINE_FATE_BASIS` -- the B2.3 wave's self-declared "strongest single
assumption" -- cites a pre-freeze probe, `scratch/b23_encoding_probe.py`, as
the evidence that sized it, and the B2.3 pre-registration sec. 0 prints that
probe's three-row table. D1 sec. 7 found the probe **cannot run on the shipped
tree**: it raises `KeyError: 'AMN'`, because `AMN` is not a species anywhere in
`src/kinetic_core/`. The shipped encoding does not DESTROY an amine pool; it
HAS no amine pool. Two of that probe's three "encodings" therefore collapse
onto the same unpatched code path and the amine-fate axis cannot be re-probed
at all.

The most likely history is benign -- an `AMN` species existed in the working
tree when the probe ran and was removed once the choice was made. That is an
honest sequence with a reproducibility defect, and the defect is what this file
repairs.

HOW THE AXIS IS RECONSTRUCTED WITHOUT RE-ADDING A SPECIES
---------------------------------------------------------
The liberated amino nitrogen is not a species, but it is DERIVABLE. The amine
centres in the sulfur lane are exactly those declared in
`ph_state.TITRATABLE_CENTRES` -- one each on `Cys`, `ARP` and `TTCA`. The
amount of amino nitrogen released by time t is therefore

    released(t) = (Cys + ARP + TTCA at t=0) - (Cys + ARP + TTCA at t)

computable from the concentration vector alone, with NO new species, NO new
parameter and NO change to any reaction. The probe patches
`ph_state.titratable_inventory` to add that quantity back as an ammonium centre
(pKa 9.25, the value `AMINE_FATE_BASIS` itself names), which realises the
"carry both centres" encoding, and compares it against the shipped encoding
("declare the amine non-titratable at the point of release").

WHAT IS AND IS NOT CLAIMED. This reconstruction is NOT the object the B2.3
probe ran against -- that object no longer exists, so its digits are not
reproducible and this file does not pretend to reproduce them. What is testable
is the DIRECTION: does carrying the released amine as ammonium make Amendment
7's three declared FIT drift anchors unreachable, as AMINE_FATE_BASIS asserts?
Pre-registered in `results/validation/kinetic_core_b2_4_prereg.md` sec. 5
before this file ran.

ANCHORS. Amendment 7's THREE DECLARED FIT drift endpoints, and Kang's 4.9 as
the out-of-sample DIAGNOSTIC B2.2's diagnosis sec. 3 already disclosed
(declared in no column, never a scored row). NO HOLD-OUT VALUE IS READ and
`data/benchmarks/external_validation/` is never opened.

Parameters are B2.2's OWN frozen vector, exactly as the original probe used --
so the comparison is of ENCODINGS at fixed chemistry, not of fits.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b2_3_fit as F  # noqa: E402
from src.kinetic_core import ph_state as PH  # noqa: E402
from src.kinetic_core import sulfur as S  # noqa: E402

#: The released amino nitrogen, carried as an ammonium centre. pKa 9.25 is the
#: value `AMINE_FATE_BASIS` itself names for surviving ammonium; dH_ion 51
#: kJ/mol is [GKL]'s primary-ammonium class value.
RELEASED_AMINE = PH.TitratableGroup(
    "released_amine", "liberated amino nitrogen carried as ammonium (RNH3+)",
    (9.25,), (51.0,), 1.0,
    "pKa 9.25, the value ph_state.AMINE_FATE_BASIS itself names for surviving "
    "ammonium; dH_ion 51 kJ/mol [GKL, primary-ammonium class]. Used ONLY by "
    "this probe, to realise the counterfactual encoding. Not imported by the "
    "engine and operative nowhere.",
)

#: The species carrying an amine centre, per ph_state.TITRATABLE_CENTRES.
AMINE_BEARING = tuple(
    key for key, centres in PH.TITRATABLE_CENTRES.items()
    if centres.get("amine")
)

CASES = (
    # name, initial, T_C, minutes, label pH, anchor, anchor class
    ("zhou_pH6", {"ARP": 20.0, "Cys": 20.0, "OX": F.OX_AMBIENT_MMOL_L},
     120.0, 60.0, 6.0, 3.22, "Amendment 7 DECLARED FIT"),
    ("zhou_pH7", {"ARP": 20.0, "Cys": 20.0, "OX": F.OX_AMBIENT_MMOL_L},
     120.0, 60.0, 7.0, 3.42, "Amendment 7 DECLARED FIT"),
    ("zhou_pH8", {"ARP": 20.0, "Cys": 20.0, "OX": F.OX_AMBIENT_MMOL_L},
     120.0, 60.0, 8.0, 5.07, "Amendment 7 DECLARED FIT"),
    ("kang_100C", {"TTCA": 10.0, "OX": F.OX_AMBIENT_MMOL_L},
     100.0, 120.0, 7.0, 4.9, "OUT-OF-SAMPLE DIAGNOSTIC (declared in no column)"),
    ("kang_120C", {"TTCA": 10.0, "OX": F.OX_AMBIENT_MMOL_L},
     120.0, 120.0, 7.0, 4.9, "OUT-OF-SAMPLE DIAGNOSTIC (declared in no column)"),
)

_ORIGINAL_INVENTORY = PH.titratable_inventory


def _patched_inventory(mode: str, initial_amine_mmol: float):
    """
    Three encodings, and unlike the B2.3 probe all three are DISTINCT on the
    current tree:

      shipped               -- the tree as it stands: the carboxyl is carried
                               into CBX, the amine is not represented.
      b22_no_carried_acid   -- CBX zeroed, i.e. B2.2's delete-the-centre
                               encoding for the carboxyl too.
      carry_both            -- CBX carried AND the released amino nitrogen
                               added back as ammonium at pKa 9.25.
    """
    def inventory(concentrations, drift, buffer_spec):
        c = dict(concentrations)
        if mode == "b22_no_carried_acid":
            c["CBX"] = 0.0
        out = list(_ORIGINAL_INVENTORY(c, drift, buffer_spec))
        if mode == "carry_both":
            present = sum(float(c.get(k, 0.0)) for k in AMINE_BEARING)
            released = max(0.0, initial_amine_mmol - present)
            if released > 0.0:
                out.append((RELEASED_AMINE, released))
        return tuple(out)
    return inventory


def run_probe() -> Dict[str, Any]:
    report = json.loads(
        (ROOT / "results/validation/kinetic_core_b2_2_fit_report.json").read_text())
    fz = report["frozen_parameters"]
    x = np.array(
        [fz["log10_k_ref_at_145C"][k] for k in F.PARAM_ORDER]
        + [fz["lumped_formation_Ea_kJ_mol"]]
        + [fz["decay_Ea_kJ_mol"][f] for f in F.DECAY_FAMILY_ORDER]
        + [fz["ph_drift"]["acid_yield_per_sink_event"],
           fz["ph_drift"]["arp_secondary_ammonium_pKa"]])
    parameters = F.build_parameters(x)
    _f, _e, _d, drift = F.unpack(x)

    rows: List[Dict[str, Any]] = []
    try:
        for mode in ("shipped", "b22_no_carried_acid", "carry_both"):
            for name, initial, t_c, minutes, ph, anchor, anchor_class in CASES:
                initial_amine = sum(
                    float(initial.get(k, 0.0)) for k in AMINE_BEARING)
                # Every caller of `titratable_inventory` is inside ph_state and
                # resolves it as a module global (lines 567, 582, 599, 870), so
                # rebinding the module attribute is sufficient and there is no
                # second reference to keep in step.
                PH.titratable_inventory = _patched_inventory(mode, initial_amine)
                run = S.integrate_sulfur(
                    parameters, t_c + 273.15, initial, np.array([0.0, minutes]),
                    ph=ph, buffer_spec=F.BUFFER_NONE, ph_drift=drift,
                    ph_nodes=41, ph_iterations=4, rtol=1e-7, atol=1e-14,
                )
                meta = run.metadata
                rows.append({
                    "encoding": mode,
                    "case": name,
                    "anchor": anchor,
                    "anchor_class": anchor_class,
                    "label_pH": ph,
                    "initial_amine_mmol_L": initial_amine,
                    "in_situ_start": float(meta["ph_initial_in_situ"]),
                    "in_situ_end": float(meta["ph_final_in_situ"]),
                    "cooled_end": float(meta["ph_final_cooled"]),
                    "miss_pH_units": float(meta["ph_final_cooled"]) - anchor,
                })
                print(f"  {mode:20s} {name:10s} cooled {meta['ph_final_cooled']:7.2f} "
                      f"(anchor {anchor})", flush=True)
    finally:
        PH.titratable_inventory = _ORIGINAL_INVENTORY

    fit_rows = [r for r in rows if r["anchor_class"].startswith("Amendment 7")]
    by_mode = {}
    for mode in ("shipped", "b22_no_carried_acid", "carry_both"):
        block = [r for r in fit_rows if r["encoding"] == mode]
        by_mode[mode] = {
            "cooled_endpoints": {r["case"]: r["cooled_end"] for r in block},
            "worst_miss_pH_units": max(abs(r["miss_pH_units"]) for r in block),
            "all_above_pH_4_5": all(r["cooled_end"] > 4.5 for r in block),
        }

    #: The B2.3 pre-registration sec. 0 printed a three-row table from the probe
    #: that no longer runs. Those digits are pinned here so that the rebuild
    #: either reproduces them or is on the record as not having.
    #: FIREWALL-OK: every number below is a COOLED pH from the B2.3
    #: pre-registration's own published table, computed at B2.2's frozen
    #: parameters against Amendment 7's three DECLARED FIT drift anchors and
    #: Kang's already-disclosed diagnostic. No hold-out row, and no volatile
    #: level of any kind, appears in this block or anywhere in this file.
    B23_PREREG_TABLE = {
        # FIREWALL-OK: cooled pH values from the B2.3 prereg's own table, not levels
        "b22_no_carried_acid": {"zhou_pH6": 3.29, "zhou_pH7": 3.49,  # FIREWALL-OK: a cooled pH, not a volatile level
                                "zhou_pH8": 5.07, "kang_120C": 11.57},
        "shipped": {"zhou_pH6": 2.94, "zhou_pH7": 3.02,
                    "zhou_pH8": 3.47, "kang_120C": 4.85},
        "carry_both": {"zhou_pH6": 5.04, "zhou_pH7": 5.11,
                       "zhou_pH8": 5.89, "kang_120C": 9.11},
    }
    by_case = {(r["encoding"], r["case"]): r["cooled_end"] for r in rows}
    reproduction = {}
    for mode, cells in B23_PREREG_TABLE.items():
        reproduction[mode] = {
            case: {"b2_3_prereg": v,
                   "b2_4_rebuild": round(by_case[(mode, case)], 2),
                   "matches_to_2dp": abs(by_case[(mode, case)] - v) < 0.005}
            for case, v in cells.items()
        }
    all_match = all(cell["matches_to_2dp"]
                    for mode in reproduction.values() for cell in mode.values())

    #: The two anchors that actually ACIDIFY are pH 6 and pH 7 (measured 3.22 /
    #: 3.42). The pH-8 anchor finishes at 5.07 and is a weaker discriminator,
    #: so a worst-miss statistic over all three is dominated by the row that
    #: discriminates least. Both are reported.
    def _acidifying_mean_abs_miss(mode: str) -> float:
        block = [r for r in fit_rows
                 if r["encoding"] == mode and r["case"] in ("zhou_pH6", "zhou_pH7")]
        return sum(abs(r["miss_pH_units"]) for r in block) / len(block)

    verdict = {
        "encodings_are_distinct": len({
            tuple(sorted(by_mode[m]["cooled_endpoints"].items()))
            for m in by_mode
        }) == 3,
        "reproduces_the_b2_3_prereg_table": all_match,
        "reproduction_cell_by_cell": reproduction,
        "mean_abs_miss_on_the_two_acidifying_anchors_pH": {
            mode: _acidifying_mean_abs_miss(mode)
            for mode in ("shipped", "b22_no_carried_acid", "carry_both")
        },
        "prereg_expectation": (
            "the rebuilt 'carry both' encoding leaves Zhou's three cooled "
            "endpoints ABOVE pH 4.5 on all three anchors, reproducing the "
            "DIRECTION of the B2.3 prereg's third row but not its digits"
        ),
        "prereg_expectation_held": by_mode["carry_both"]["all_above_pH_4_5"],
        "shipped_worst_miss_pH_units": by_mode["shipped"]["worst_miss_pH_units"],
        "carry_both_worst_miss_pH_units": by_mode["carry_both"]["worst_miss_pH_units"],
    }

    return {
        "artifact": "b2_4_amine_fate_probe",
        "replaces": "scratch/b23_encoding_probe.py, which raises KeyError: 'AMN' "
                    "on the shipped tree (D1 sec. 7)",
        "parameters": "B2.2's OWN frozen vector, unchanged",
        "amine_bearing_species": list(AMINE_BEARING),
        "reconstruction": (
            "released(t) = (Cys + ARP + TTCA at t=0) - (Cys + ARP + TTCA at t), "
            "added back as an ammonium centre at pKa 9.25. No new species, no "
            "new parameter, no reaction changed."
        ),
        "rows": rows,
        "by_encoding_on_the_three_FIT_anchors": by_mode,
        "verdict": verdict,
    }


def main() -> int:
    payload = run_probe()
    dest = ROOT / "results/validation/kinetic_core_b2_4_amine_fate_probe.json"
    dest.write_text(json.dumps(payload, indent=2, default=str))
    print(f"wrote {dest}")
    print(json.dumps(payload["verdict"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
