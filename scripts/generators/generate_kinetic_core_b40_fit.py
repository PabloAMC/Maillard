#!/usr/bin/env python
"""
Build Wave B40 (2026-09-11): B39's fed 3-deoxyglucosone fit re-run with the 3-DG exits' pH term
switched on (trunk_conditions.THREE_DEOXY_EXIT_PH_TERM, declared from Martins 2003 Table 3).
Pre-registration: results/validation/kinetic_core_b40_prereg.md. Same rows, sigmas, bounds and
starts as B39; the Leitzen hold-out is never read.

Run inside docker:  PYTHONPATH=/workspace python scripts/generators/generate_kinetic_core_b40_fit.py
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import trunk_conditions as TC  # noqa: E402

PREREG = data_paths.VALIDATION_DIR / "kinetic_core_b40_prereg.md"
OUT_JSON = data_paths.VALIDATION_DIR / "kinetic_core_b40_fit_report.json"
OUT_MD = data_paths.VALIDATION_DIR / "kinetic_core_b40_fit_report.md"
B39_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b39_fit_report.json"


def _b39():
    spec = importlib.util.spec_from_file_location("b39", ROOT / "scripts" / "generators" / "generate_kinetic_core_b39_fit.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--max-nfev", type=int, default=250)
    ap.add_argument("--starts", type=int, default=2)
    a = ap.parse_args(argv)
    assert PREREG.exists(), "B40 is pre-registered; write the prereg before running"
    TC.THREE_DEOXY_EXIT_PH_TERM = True          # the wave's one structural switch, for this process
    TC.THREE_DEOXY_EXIT_PH = dict(TC.THREE_DEOXY_EXIT_PH, **TC.THREE_DEOXY_EXIT_PH_REJECTED_B40)   # B40: BOTH exits, as pre-registered
    m = _b39()
    m.build._term_set_by_wrapper = True
    t0 = time.time()
    p = m.build(a.quick, a.max_nfev, a.starts)
    p["artifact"] = "kinetic_core_b40_fit_report"; p["wave"] = "B40"; p["prereg"] = data_paths.rel(PREREG)
    p["provenance"] = provenance.provenance_block("kinetic_core_b40_fit_report", generated_by="scripts/generators/generate_kinetic_core_b40_fit.py",
                                                  wave="B40", inputs=[PREREG, B39_REPORT, ROOT / "data/lit/extraction_dossiers/martins2003_extraction.md",
                                                                      ROOT / "data/lit/extraction_dossiers/mittelmaier2010_extraction.md"])
    p["declaration"]["three_deoxy_exit_ph_term"] = {k: {"exponent": v[0], "band": list(v[1]), "source": v[2]} for k, v in TC.THREE_DEOXY_EXIT_PH.items()}
    b39 = json.loads(B39_REPORT.read_text())
    d39 = b39["frozen_parameters"]["fed_3deoxy"]["log10_k_tdg_ddg_100C"] - b39["frozen_parameters"]["start_log10"]["log10_k_tdg_ddg_100C"]
    d40 = p["frozen_parameters"]["fed_3deoxy"]["log10_k_tdg_ddg_100C"] - p["frozen_parameters"]["start_log10"]["log10_k_tdg_ddg_100C"]
    p["predictions"] = {
        "P1_rows_fit_with_tmax": p["predictions"]["P1_rows_fit"],
        "P2_k_tdg_ddg_moves_less_than_b39": {"b39_delta_dex": d39, "b40_delta_dex": d40, "held": bool(abs(d40) < abs(d39))},
        "P5_pinned": p["predictions"]["P5_three_pinned_reverse_pair_collinear"],
    }
    p["wall_seconds"] = round(time.time() - t0, 1)
    OUT_JSON.write_text(json.dumps(p, indent=2, default=str) + "\n")
    OUT_MD.write_text(m.render(p).replace("# Wave B39", "# Wave B40").replace("kinetic_core_b39_prereg.md", "kinetic_core_b40_prereg.md"))
    print(f"wrote {data_paths.rel(OUT_JSON)}: cost {min(mm['cost'] for mm in p['members']):.2f} | chi2_red {p['laplace']['chi2_reduced']:.2f} | "
          f"x = {json.dumps({k: round(v, 3) for k, v in p['frozen_parameters']['fed_3deoxy'].items()})} | tmax {p['fed_pot_before_after']['after']['TDG:tmax_DDG']} | {p['wall_seconds']} s")
    for k, v in p["predictions"].items():
        print(f"  {k}: {'HELD' if v['held'] else 'REFUTED'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
