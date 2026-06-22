"""Diagnostic: measured vs predicted lipid-aldehyde ppb on the in-panel anchors
and the frozen external hold-out.

Read-only. Establishes the baseline bias the lipid-oxidation recalibration
(Workstream A) must close. The hold-out bundles are printed for context only —
they are NEVER used to fit any parameter.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark  # noqa: E402
BENCH = ROOT / "data" / "benchmarks"

LIPID_MARKERS = {"hexanal", "nonanal", "2-pentylfuran", "1-hexanol", "hexanol"}

IN_PANEL = [
    "pea_isolate_40C_PratapSingh2021",
    "soy_isolate_40C_PratapSingh2021",
    "pea_isolate_uht_140C_Trikusuma2019",
    "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
    "soy_isolate_ribose_cysteine_100C_45min_Internal2026",
]
HOLD_OUT = [
    "external_validation/external_validation_bi_2020_raw_pea_hexanal",
    "external_validation/external_validation_bi_2020_roasted_pea_hexanal",
    "external_validation/external_validation_li_2026_spi_wg_hme_control",
    "external_validation/external_validation_liu_2023_ppi_offnote_baseline",
]


def _is_lipid(name: str) -> bool:
    low = str(name).lower()
    return any(tok in low for tok in LIPID_MARKERS)


def _report(group: str, ids: list[str]) -> list[float]:
    print(f"\n=== {group} ===")
    fold_errors: list[float] = []
    for bid in ids:
        path = BENCH / f"{bid}.json"
        if not path.exists():
            print(f"  (missing) {bid}")
            continue
        ev = evaluate_benchmark(path)
        for c in ev.comparisons:
            if c.matched_name is None or not _is_lipid(c.matched_name):
                continue
            meas, pred = float(c.measured_ppb), float(c.predicted_ppb)
            if meas <= 0 or pred <= 0:
                continue
            fold = max(pred / meas, meas / pred)
            fold_errors.append(fold)
            print(f"  {Path(bid).name:52s} {c.matched_name:16s} "
                  f"meas={meas:10.1f}  pred={pred:12.1f}  fold={fold:9.1f}x")
    return fold_errors


def main() -> None:
    in_fold = _report("IN-PANEL anchors (calibration-eligible)", IN_PANEL)
    out_fold = _report("EXTERNAL hold-out (frozen — context only)", HOLD_OUT)
    for label, fe in (("in-panel", in_fold), ("hold-out", out_fold)):
        if fe:
            med = sorted(fe)[len(fe) // 2]
            print(f"\n{label}: median fold error = {med:.1f}x over {len(fe)} lipid points")


if __name__ == "__main__":
    main()
