"""
WAVE B45 -- THE TWO COMPLETENESS PROBES, MADE REPRODUCIBLE.

`results/validation/kinetic_core_b45_prereg.md` quotes two numbers as evidence:
the fraction of hexanal the engine's declared binding block removes on the pot
Shi et al. 2022 heated, and the fraction of charged cysteine the engine retains
on the pot Baldus et al. 2017 heated. Neither number came from a fit and neither
moved a constant, so neither has a fit report or a ship rule to pin it. This
generator exists so that they cannot drift silently anyway: it recomputes both
from the shipped engine and writes them where a test can read them back.

NEITHER PROBE IS SCORED AND NEITHER IS A BENCHMARK. Both are completeness
checks -- they ask whether the model contains a channel at all, not whether a
constant is the right size. The measured comparators are recorded beside the
model's answers so the comparison in the pre-registration can be audited without
re-reading the papers.

P1 -- the aldehyde binding block against Shi's pot. Shi acidified a soy protein
isolate to pH 4.5 at 30 mg/mL and held it at 95 C for 5 min; headspace hexanal
rose from 23 +/- 0.59 to 66 +/- 1.7 ug/L, a 2.87x RELEASE. The engine's block
(`matrix_sites`) is a one-way covalent adduct to lysine, so it can only remove
hexanal, and the question is how much.

P2 -- the sulfur lane's thiol removal against Baldus's pot. Baldus held ~228-250
uM cysteine in 0.1 M acetate at pH 5.5, air-saturated, and reported that "Cys was
almost completely degraded in 5 minutes during heating from 40-60 C" in the
presence of 18 uM Cu(II)EDTA. The probe charges cysteine alone with no sugar and
no metal -- the engine has no metal-catalysed channel to charge -- and reads what
survives. `hydrogen sulfide` is the target because a sugar-free cysteine pot is
correctly refused for any target it cannot reach.
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.kinetic_core import matrix_sites as MS  # noqa: E402
from src.kinetic_core.engine import (  # noqa: E402
    FormulationSpec, ProcessSpec, ThermalProgram, predict,
)

OUT = ROOT / "results" / "validation" / "kinetic_core_b45_probe.json"

#: (label, segments) -- Shi's own hold first, then three the paper did not run,
#: to show the block stays small even well outside her conditions.
P1_PROGRAMMES = (
    ("shi_hold_95C_5min", ((5.0, 95.0),)),
    ("shi_spi_prep_95C_10min", ((10.0, 95.0),)),
    ("long_cook_100C_45min", ((45.0, 100.0),)),
    ("extreme_140C_60min", ((60.0, 140.0),)),
)
P2_HOLDS_MIN = (5.0, 30.0, 180.0)


class _Proc:
    """The two attributes `matrix_sites.resolve` reads. Not a ProcessSpec: the
    binding block is evaluated on its own, with no integration behind it."""

    def __init__(self, matrix: str, protein_g_per_l: float, ph: float) -> None:
        self.matrix, self.protein_g_per_l, self.ph = matrix, protein_g_per_l, ph
        self.protein_sites = None


def probe_p1() -> dict:
    out = {}
    for matrix in ("soy_isolate", "pea_isolate"):
        charged, note = MS.resolve(_Proc(matrix, 30.0, 4.5))
        if charged is None:
            out[matrix] = {"charged": False, "note": note}
            continue
        rows = {}
        for label, segments in P1_PROGRAMMES:
            bf = MS.bound_fraction("HEXANAL", charged, segments)
            rows[label] = None if bf is None else {
                "bound_fraction": bf["bound_fraction"],
                "bound_fraction_corners": bf["bound_fraction_corners"],
                "classes": bf["classes"],
            }
        out[matrix] = {"charged": True, "amine_mmol_per_l": charged.amine, "programmes": rows}
    out["measured_comparator"] = {
        "source": "shi2022_extraction.md (10.1111/jfpp.16555), Table 2 row 1 and Fig. 4",
        "acidic_spi_30mg_ml_ph4p5_95C_5min_ug_per_l": {"before": 23.0, "before_sd": 0.59,
                                                       "after": 66.0, "after_sd": 1.7},
        "release_factor": 66.0 / 23.0,
        "acidic_11s_ug_per_l_figure_only": {"before": 40.0, "after": 220.0},
        "direction": "RELEASE on heating; the engine's block can only BIND",
    }
    return out


def probe_p2() -> dict:
    rows = {}
    for minutes in P2_HOLDS_MIN:
        spec = FormulationSpec(
            name=f"b45_p2_cysteine_only_{minutes:g}min",
            precursors={"Cysteine": 0.25},
            process=ProcessSpec(thermal=ThermalProgram(((minutes, 95.0),)), ph=5.5),
        )
        prediction = predict(spec, ["hydrogen sulfide"])
        if not prediction.declaration.is_answerable:
            rows[f"{minutes:g}min"] = {"answerable": False}
            continue
        state = prediction.species_mmol_per_l
        cys = float(state["Cys"])
        rows[f"{minutes:g}min"] = {
            "answerable": True,
            "cys_mmol_per_l": cys,
            "cys_fraction_remaining": cys / 0.25,
            "ox": state.get("OX"),
            "oxv": state.get("OXV"),
        }
    return {
        "charged_cysteine_mmol_per_l": 0.25,
        "temperature_c": 95.0,
        "ph": 5.5,
        "holds": rows,
        "measured_comparator": {
            "source": "baldus2017_extraction.md (10.1021/acs.jafc.6b05472)",
            "system": "~228-250 uM cysteine, 0.1 M acetate pH 5.5, air-saturated 8 mg/L O2, 18 uM Cu(II)EDTA",
            "statement": "Cys was almost completely degraded in 5 minutes during heating from 40-60 C",
            "t0_measured_um": 220.0,
            "h2o2_peak_um": 71.0,
            "trace_cu_from_300um_cysteine_um": 0.26,
            "all_metals_in_ultrapure_buffer_um": "< 0.08 (ICP-OES detection limit)",
        },
    }


def main() -> int:
    payload = {
        "artifact": "kinetic_core_b45_probe",
        "prereg": "results/validation/kinetic_core_b45_prereg.md",
        "what_this_is": (
            "Two COMPLETENESS probes from wave B45. Neither is scored, neither is a benchmark, and "
            "neither moved a constant. They ask whether the model contains a channel at all. P1: the "
            "declared aldehyde-binding block against a measured hexanal RELEASE. P2: the sulfur lane's "
            "thiol removal against a measured metal-catalysed loss the model has no channel for."
        ),
        "p1_hexanal_binding": probe_p1(),
        "p2_cysteine_removal": probe_p2(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    p1 = payload["p1_hexanal_binding"]["soy_isolate"]["programmes"]["shi_hold_95C_5min"]["bound_fraction"]
    p2 = payload["p2_cysteine_removal"]["holds"]["5min"]["cys_fraction_remaining"]
    print(f"wrote {OUT.relative_to(ROOT)}: P1 binds {p1:.4%} of hexanal on Shi's hold "
          f"(measured: +187 % release) | P2 keeps {p2:.2%} of cysteine at 5 min "
          f"(measured: almost completely degraded)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
