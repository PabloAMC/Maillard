# Legacy lane archive (retirement step B5b, 2026-09-03)

Everything in this directory was produced by the retired SMIRKS/Hammond screening lane and its
validation harness (`src/benchmark_validation.py`, `src/uncertainty_propagation.py`,
`src/external_validation.py`, the matrix calibration layer, the family/matrix promotion reports).
The code that wrote these files was deleted at commit "Retirement step B5b"; the files are kept
as the historical record of what that lane claimed and scored, and NOTHING in the tree reads them.
The kinetic core's own artifacts live in `results/validation/` (`core_panel_scores.*`,
`core_prediction_uncertainty.*`, `cutover_final_exam.*`, `kinetic_core_b*_*`).

To see the code that produced any file here: `git show pre-data-restructure-2026-09-01:<path>` or
any commit before B5b.

## Archived later (2026-09-03, backlog pass)

- `validation/experiment_requests/` — written by `src/experiment_request.py`, deleted at B5b; nothing read them.
- `validation/active_learning_requests.json` — written by `doe_generator.export_doe_requests`, which had no caller (removed).
- `validation/family_coverage.png` — no producer and no reader in the tree.
