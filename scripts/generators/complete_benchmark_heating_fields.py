"""
Complete `conditions.heating` on the bundles whose source records a heat-up (wave B64, 2026-10-28).

Write the thermal history each source states into its bundles, with the quote and where it is.

A curation script, like `complete_benchmark_vessel_fields.py`: the blocks are declared here, from
`docs/validation/thermal_history_facts.yml` (every quote read in the PDF), and written into the bundles;
`--check` reports drift and writes nothing. A block is read by the engine only when
`panel.THERMAL_HISTORY = "as_recorded"` (`src/kinetic_core/thermal_history.py`). Bundles whose source states
no come-up carry no block and stay isothermal; bundles whose stated come-up a profile cannot yet encode
(Leitzen 2021's F0, Lin 2021/2022's unstated dry-bath come-up) are listed in the B64 prereg, not guessed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

HOFMANN_1998_RAMP = {
    "time_convention": "ramp_only",
    "start_C": 20.0,
    "ramp_min": 20.0,
    "provenance_class": "primary_source_pdf",
    "provenance_note": (
        "B64 (2026-10-28). Hofmann & Schieberle 1998 (p. 236) writes 'thermally treated for 20 min at 145 °C "
        "in a laboratory autoclave (200 mL; Type II; Roth, Germany)'. The same laboratory runs the same pot "
        "in the same autoclave 'by raising the temperature within 20 min from 20 to 145 °C' (Hofmann & "
        "Schieberle 1995, p. 2; 1997, p. 901, where the Results call that pot 'thermally treated for 20 min "
        "at 145 °C'). No paper states a hold. Read as a ramp with no hold (p ~ 0.8; ramp + 20 min hold "
        "~ 0.12, carried as B64's alternative). Cooling ('after cooling of the reaction mixtures to "
        "20 °C') is unstated in rate and omitted. Record: docs/validation/thermal_history_facts.yml."
    ),
}

#: benchmark_id -> heating block. Every Hofmann 1998 aqueous bundle, scored or not.
BLOCKS = {
    path.stem: HOFMANN_1998_RAMP
    for path in sorted(data_paths.BENCHMARKS_DIR.rglob("*hofmann1998*.json"))
    if "quarantined" not in path.parts
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[3])
    parser.add_argument("--check", action="store_true", help="report drift, write nothing")
    args = parser.parse_args()
    drift = []
    for path in sorted(data_paths.BENCHMARKS_DIR.rglob("*hofmann1998*.json")):
        if "quarantined" in path.parts:
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        conditions = payload.get("conditions") or {}
        assert float(conditions.get("temp_C")) == 145.0 and float(conditions.get("time_min")) == 20.0, path
        block = BLOCKS[path.stem]
        if conditions.get("heating") == block:
            continue
        drift.append(path.stem)
        if not args.check:
            payload.setdefault("conditions", {})["heating"] = block
            path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.check and drift:
        print(f"heating blocks stale on {len(drift)} bundle(s): {drift}", file=sys.stderr)
        return 1
    print(f"{'stale' if args.check else 'wrote'}: {len(drift)} bundle(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
