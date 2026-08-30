import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_observable_closure_audit
from src.presentation import render_matrix_observable_closure_audit_markdown


def test_matrix_observable_closure_audit_labels_transfer_vs_mechanistic_actions():
    payload = build_matrix_observable_closure_audit([
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_uht_140C_Trikusuma2019.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json",
    ])

    by_id = {row["benchmark_id"]: row for row in payload["benchmarks"]}
    pea_candidate = by_id["pea_isolate_ribose_cysteine_100C_45min_Internal2026"]
    compounds = {row["compound"]: row for row in pea_candidate["compounds"]}
    watchlist = {row["benchmark_id"]: row for row in payload["mechanistic_refinement_watchlist"]}

    assert compounds["2-furfurylthiol"]["closure_action"] == "class_level_transfer_acceptable"

    # RE-PINNED 2026-08-27 (Wave S1 fix 2 -- THE MATRIX REGISTRY BECAME REACHABLE).
    # Hexanal: "mechanistic_blocker" -> "class_level_transfer_acceptable", and the
    # mechanistic-refinement watchlist went from one entry to EMPTY. NOTHING ABOUT THE
    # CHEMISTRY CHANGED. Until this wave, species on the `matrix_precursor_augmented` lane
    # reached `describe_matrix_calibration` labelled by canonical SMILES ("CCCCCC=O"), so
    # the name-keyed registry missed and every row came back
    # `calibration_evidence_strength = "heuristic"` / `fallback_mode = "class_level"`.
    # `_matrix_closure_action` routes internal-candidate + heuristic straight to
    # `mechanistic_blocker`, which is how Hexanal earned the flag. With the lookup fixed the
    # row now resolves a real compound-specific record and reads
    # `calibration_evidence_strength = "process_state_mismatch"` (the registry's honest label
    # for a record reached through the intermediate_matrix -> ambient_slurry fallback), and
    # with `evidence_state = "externally_benchmarked"` that lands in the
    # class-level-transfer branch instead.
    #
    # THIS IS PINNED, NOT ENDORSED, AND IT IS A LOSS OF SIGNAL. `_matrix_closure_action` has
    # NO branch for `process_state_mismatch`: a factor reached only by a process-state
    # fallback is scored exactly like a genuine class-anchored transfer, so a flag that was
    # correctly warning about an uncalibrated lane has been switched off by a lookup repair.
    # Carried as an open item in tasks/audit_remediation.md (Wave S1); the repair is to the
    # rule set, not to this assertion, and it was deliberately NOT made in the same pass as
    # the structural fix that exposed it.
    assert compounds["Hexanal"]["closure_action"] == "class_level_transfer_acceptable"
    assert compounds["Hexanal"]["calibration_evidence_strength"] == "process_state_mismatch"
    assert payload["mechanistic_refinement_watchlist"] == []
    assert watchlist == {}

    markdown = render_matrix_observable_closure_audit_markdown(payload)
    assert "Matrix Observable Closure Audit" in markdown
    assert "class_level_transfer_acceptable" in markdown
    # RE-PINNED 2026-08-27 (Wave S1 fix 2): the renderer emits the "Mechanistic Refinement
    # Watchlist" SECTION only when the watchlist is non-empty; with it empty it prints the
    # count line instead. Pinning the count line is strictly stronger than pinning the
    # heading -- it asserts the number, not merely that a section exists.
    assert "Mechanistic watchlist count: 0" in markdown
    assert "Mechanistic Refinement Watchlist" not in markdown