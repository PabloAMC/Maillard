"""Refresh the synthetic internal reproducibility snapshots from the current model.

The Internal2026 benchmarks and the ProtocolPilot2026 intake payloads are
SYNTHETIC comparators: their volatile values are frozen model output used only
for drift detection (see the 2026-08-26 audit and tasks/audit_remediation.md).
They were previously hand-frozen in commits with no generator; this script makes
the refresh reproducible and forces the provenance labels to stay honest.

Run after any intentional model change that moves the matrix_precursor_augmented
path, then commit the diff so the drift baseline is explicit in review:

    python scripts/generators/refresh_internal_reproducibility_snapshots.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark_payload, load_benchmark
from src.matrix_experiment_intake import materialize_matrix_experiment_benchmark

GENERATOR_TAG = "reproducibility_snapshot_refresh_v5"
NOTES = (
    "SYNTHETIC reproducibility snapshot: values are frozen output of the current "
    "model (refreshed via scripts/generators/refresh_internal_reproducibility_snapshots.py "
    "after the 2026-08-27 G-wave chemistry rebuild). Folded in, newest first: "
    "(G1) the fabricated lipid radical chain was deleted -- 93% of the emitted network, "
    "5500 -> 369 steps -- and hexanal, norfuraneol, glyoxal, CML, DMHF and the correct "
    "disulfide regioisomer became reachable; the accepted 1-deoxyosone -> norfuraneol -> "
    "MFT route (van den Ouweland & Peer 1975) replaced the fabricated one-step shortcut, "
    "which is demoted but retained as a lump; eight families that had been silently "
    "switched off by the 45 kcal/mol DEFAULT_BARRIER fallthrough (acrylamide, CML/CEL, "
    "furanone, thiamine, GSH) were given explicit estimated barriers and are now live. "
    "(Wave H) `Enolisation_2_3_Amadori` stopped collecting the amine-nucleophile "
    "ionisation and Amadori water-activity corrections it was matching by name alone, "
    "which is what had left the new MFT route inert; `thiol_addition_norfuraneol` was "
    "refit 28.60 -> 26.85 against cys_ribose_140C_Hofmann1998 only; the Methional "
    "hydrolysate observability factor was re-derived 0.0045 -> 0.05623 against the two "
    "literature xylose HVP benchmarks. (2026-08-27 earlier) the projection retune: an "
    "apparent-Arrhenius conversion extent replaced the saturating severity sigmoid and "
    "the double-applied Boltzmann selectivity on path_span was de-duplicated. "
    "(2026-08-26) the audit bug fixes -- marker-injection mM basis, 2-pentylfuran SMILES, "
    "full MolWt. This is a drift-detection baseline, not a measurement, and carries zero "
    "calibration evidence."
)


def _refresh_contract_ranks(contract: dict, predicted: dict) -> None:
    targets = contract.get("observable_targets", [])
    ranked = sorted(
        targets,
        key=lambda item: -float(predicted.get(str(item.get("name")), 0.0)),
    )
    for rank, item in enumerate(ranked, start=1):
        item["expected_rank"] = rank
    contract["observable_targets"] = ranked
    notes = str(contract.get("notes", ""))
    marker = "Expected ranks track the current model snapshot"
    if marker not in notes:
        contract["notes"] = (
            notes + " " if notes else ""
        ) + "Expected ranks track the current model snapshot (drift baseline), not measured truth."


def refresh_internal(protein: str) -> dict:
    path = ROOT / "data" / "benchmarks" / f"{protein}_isolate_ribose_cysteine_100C_45min_Internal2026.json"
    bench = load_benchmark(path)
    evaluation = evaluate_benchmark_payload(bench)
    predicted = dict(evaluation.predicted_ppb or {})

    reference = bench.get("reference_volatiles") or {}
    for compound, entry in reference.items():
        if compound in predicted:
            entry["conc_ppb"] = float(predicted[compound])
    _refresh_contract_ranks(bench.get("matrix_ranking_contract") or {}, predicted)

    meta = bench.setdefault("source_metadata", {})
    meta["generator"] = GENERATOR_TAG
    meta["origin"] = "internal_reproducibility_candidate"
    meta["notes"] = NOTES
    meta.pop("measurement_date", None)

    path.write_text(json.dumps(bench, indent=2, sort_keys=True) + "\n")
    return predicted


def refresh_protocol_pilot(protein: str, predicted: dict) -> None:
    yaml_path = ROOT / "data" / "protocols" / f"{protein}_iso_protocol_pilot_intake.yaml"
    payload = yaml.safe_load(yaml_path.read_text())

    measured = payload.get("measured_volatiles") or {}
    for compound, entry in measured.items():
        if compound in predicted:
            entry["conc_ppb"] = float(predicted[compound])
    _refresh_contract_ranks(payload.get("comparison_contract") or {}, predicted)

    yaml_path.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True))

    bench = materialize_matrix_experiment_benchmark(yaml_path)
    out = ROOT / "data" / "benchmarks" / f"{protein}_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json"
    out.write_text(json.dumps(bench, indent=2, sort_keys=True) + "\n")


def main() -> int:
    for protein in ("pea", "soy"):
        predicted = refresh_internal(protein)
        refresh_protocol_pilot(protein, predicted)
        print(f"{protein}: refreshed Internal2026 + ProtocolPilot2026 snapshots")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
