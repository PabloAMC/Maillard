# H0 nested-tau: what the implementation reads

These are all read-only inputs. They were read in a laptop checkout of branch `thiol-lolo-2026-11-01`; since PR #35 merged, they are on main.
EX-056's records were written at commit `7b1bc7cd`; HEAD is now `3cbfcdaf`.

## Saved artifacts (the primary analysis needs nothing else)

| file | fields used |
|---|---|
| `results/validation/kinetic_core_thiol_lolo_<lab>.json`, for each of dalian_polytechnic_chen_d, dfa_hofmann_schieberle, jiangnan_zhang_x, nestle_blank, ochanomizu_murata and reading_mottram | see the next table |
| `results/validation/kinetic_core_thiol_lolo.json` | `provenance.inputs[].sha256` (the input-hash check); `verdict`; `sets` (reproduction of the pooled numbers) |
| `results/validation/kinetic_core_thiol_lolo_prereg.md` | sec. 4.5: the definitions of the shipped arm (a) |
| `results/validation/kinetic_core_ph_acid_submodel_fit_{A,B,C}.json` | `leave_one_laboratory_out`: read only to document that it is NOT usable as a second split (it is aggregated only) |

Fields used from each per-fold record:

| field | what it is |
|---|---|
| `rows[]` | `ledger_id`, `species`, `t_min`, `method`, `paired`, `engine_raw.warm` / `.reset`, `null_raw.typical_yield` / `.no_change`: raw held-out residuals for ALL held-out rows, in held-out order |
| `models.{engine_warm, engine_reset, null_typical_yield, null_no_change}.per_row[]` | (`centred`, `var`) for the paired rows, in paired order: the reproduction target for arm (a) |
| `models.*.terms.class_noise` | the fold's pooled class noise per model, which gives s_i |
| `models.*.terms.between_tau` / `between_tau_used` / `fallbacks` | arm (a)'s in-sample tau |
| `models.*.log_density`, `models.*.coverage` | the reproduction targets |
| `paired_rows`, `leak_prone`, `complete`, `provenance.git` | bookkeeping |
| `warm.laplace.{sigma, correlation, at_bound, identified}`, `reset.laplace` | S2 only |

## Curated inputs re-read (no engine run)

- Row sigma, `series` and table/pot are not in the records, so the rows are rebuilt with
  `fit_thiol_lolo.load_fit_table()` and `make_folds()`. This is `fit_rows.load` plus `build_submodel_probe.plan`:
  metadata only, no integration.
- The tables: `data/lit/fit_rows/{deoxyosone_submodel, thiol_submodel, thiol_zhou2023_actz, ph_acid_davidek_acids,
  ph_acid_submodel}.yml` and `data/keys/laboratories.yml`.
- Their sha256 must equal the provenance hashes in `kinetic_core_thiol_lolo.json`. On 2026-11-01 the thiol and
  deoxyosone tables and `src/kinetic_core/fitting.py` matched. If any hash differs, read the files at 7b1bc7cd
  (`git show`).
- The held-out rows are joined to `rows[]` by position, and (`ledger_id`, `species`, `t_min`, `method`) is asserted
  equal row by row.

## Code imported (functions only)

- `fit_thiol_lolo`: `lab_of`, `level_of`, `_model`, `predictive_model`, `row_intervals`, `score`.
- `fit_joint_deoxyosone_thiol`: `series_of`, `method_of`.
- `fit_deoxyosone_submodel.SCALE_TAU`.
- `src.kinetic_core.fitting.MeasurementModel` (`log_density`, `coverage`, `sigma`) and `estimate_between_tau`. The
  latter is used for the identity test only.
- `fit_joint_deoxyosone_thiol.py` changed after 7b1bc7cd (the A3 phase 1 Config refactor; `series_of` and `method_of`
  are untouched). The per-row reproduction check (prereg sec. 5.1) is the guard.

## Needs more than saved artifacts

- **PRIMARY (prereg v2, sec. 3.1): strict nesting.** 10 leave-two-labs-out refits (`fit_thiol_lolo.run_fold` logic
  with two labs held out), ~4 - 8 core-hours, niced, at most 8 cores.
  - Guarded by prereg sec. 4.2: the cost reproduction, two re-predictions of saved folds, and one re-run of the
    Reading fold.
  - If the current engine fails the guard, the refits run pinned at 2e031fdb (7b1bc7cd's engine plus the generator).
  - The strict nulls need only the reloaded rows and pots (`load_fit_table`), with no refit.
- **Optional: S2,** the parameter-covariance term. It needs held-out-pot Jacobians at each fold optimum (engine
  evaluations, no refit).
- The one-deep sensitivity and the controls need the saved artifacts only.

## Outputs (proposed; written only when H0 is run)

- `scripts/generators/score_thiol_lolo_nested_tau.py`
- `results/validation/kinetic_core_thiol_lolo_nested_tau_controls.json` (committed before the real scoring)
- `results/validation/kinetic_core_thiol_lolo_nested_tau.{json,md}`
