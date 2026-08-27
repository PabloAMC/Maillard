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

GENERATOR_TAG = "reproducibility_snapshot_refresh_v7"
NOTES = (
    "SYNTHETIC reproducibility snapshot: values are frozen output of the current "
    "model (refreshed via scripts/generators/refresh_internal_reproducibility_snapshots.py "
    "after the 2026-08-27 Wave N route correction). Folded in, newest first: "
    "(Wave N) THE MFT ROUTE WAS CORRECTED ON ISOTOPE EVIDENCE. The norfuraneol + H2S -> MFT "
    "step (`Thiol_Addition_Norfuraneol`) was RETIRED: Cerny & Davidek 2003 "
    "(10.1021/jf026123f) spiked authentic norfuraneol into a [13C5]ribose/cysteine system "
    "and the MFT came out mainly 13C5-labelled, i.e. norfuraneol is 'unimportant as an "
    "intermediate'; Cerny & Davidek 2004 (10.1021/jf035265m) positionally confirmed "
    "1,4-dideoxypento-2,3-diulose with [1-13C]ribose. Two families replace it: "
    "`Deoxyosone_Reduction` (1-deoxy-2,3-pentodiulose + 2[H] -> 1,4-dideoxypento-2,3-diulose "
    "+ H2O, barrier 28.0 ESTIMATED) and `Thiol_Addition_Pentodiulose` (1,4-dideoxyosone + "
    "H2S -> MFT + 2 H2O, barrier 28.60 ESTIMATED and UNCONSTRAINED, deliberately NOT the "
    "26.85 that Wave H fitted through the contradicted route, and NOT tuned). The H2S step "
    "no longer consumes a reducing-equivalent token; the reduction moved upstream. Norfuraneol "
    "is retained as a genuine furanone product, demoted out of the MFT lane. Honest effect on "
    "cys_ribose_140C_Hofmann1998: MFT 235.32 -> 151.87 ppb against 342 (1.45x under -> 2.25x "
    "under) and FFT 219.96 -> 243.72 against 200 (1.10x over -> 1.22x over). The numbers got "
    "WORSE and were NOT clawed back: the 1.45x had been bought with a barrier fitted through "
    "a route the isotope evidence contradicts. "
    "(Wave I) the norfuraneol + H2S -> MFT step no longer draws its 2[H] reducing "
    "equivalent from a pool whose only source was pyrazine chemistry -- the reductant is "
    "now sourced from the thiol redox couple (2 cysteine -> cystine + 2[H]), so MFT no "
    "longer goes to zero when aminoketone condensation is disabled. Measured effect on "
    "cys_ribose_140C_Hofmann1998: MFT 61.25 -> 235.32 ppb against a 342 ppb reference "
    "(5.58x under -> 1.45x under) and FFT 61.44 -> 219.96 against 200 (3.26x under -> "
    "1.10x over), with NO barrier changed. The Lipid_Schiff_Base SMARTS finally implements "
    "the exclusions its comment claimed (no alpha-hydroxy carbonyls, C3+ only), taking the "
    "family from 75.2% of the core no-lipid network to 27.2%. The demoted one-step MFT "
    "shortcut no longer fires for pentoses. The Methional hydrolysate observability factor "
    "was REVERTED 0.05623 -> 0.0045: the two benchmarks it had been re-derived against are "
    "fabricated and are now quarantined. 10 of 14 Monte-Carlo barrier priors were exact "
    "no-ops and now route correctly, widening every confidence interval. "
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
    "refit 28.60 -> 26.85 against cys_ribose_140C_Hofmann1998 only (that refit is now "
    "STALE -- it was performed against a network in which the MFT route was throttled by "
    "the H2 coupling Wave I removed, so its 'no barrier value can fix this' conclusion "
    "does not carry over and re-running it is an open owner item). "
    "(2026-08-27 earlier) the projection retune: an "
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
