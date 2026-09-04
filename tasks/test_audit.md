# Test-suite audit (step 7 / plan item B6) — 2026-09-03

Question per test: **what change would make it fail?** A test that no plausible change fails is not a
test. Inventory: 49 files, 590 test functions after the post-retirement steps; measured line coverage
81 % of `src/` (78 % of `src/kinetic_core/`), 11 % of `scripts/` (the frozen wave generators, on
purpose). Automated pass: 93 functions assert only shape (`in` / `isinstance` / `exists`) or nothing
at all; 12 files read tracked `results/validation/` artifacts.

Verdicts: **keep** (asserts values or behaviour), **keep, label** (contract on a frozen artifact —
say so in the docstring), **merge**, **rewrite** (shape-only tests that should assert a value),
**delete** (tautology, or tests a deleted surface).

## The kinetic core (the engine's own contracts)

| file | tests | verdict | what a maintainer should know |
| --- | ---: | --- | --- |
| `test_kinetic_core_b1.py` | 31 | keep | Trunk network: element balance, Amadori split, fit-report round trip. 2 shape-only (`melanoidin_sink_is_terminal`, `no_dft_module_is_imported`) are structural guards; fine. |
| `test_kinetic_core_b2.py` | 45 | keep, label | Sulfur network. 6 shape-only, of which `the_fit_script_reads_no_benchmark_or_holdout_file` and `the_holdout_scorer_contains_no_optimiser` read the FROZEN generators by text: label them as freeze checks (now redundant with the manifest guard — **merge into `test_wave_generators_frozen.py`** later). |
| `test_kinetic_core_b2_1.py` … `b2_4.py` | 26 / 32 / 29 / 20 | keep, label | Each wave's contract on its own frozen report (prereg exists, amendments named, residuals within sigma). They are regression records: label the docstrings "frozen wave; fails only if the report or the network changes". `b2_3`'s `a_primary_reading_names_the_pdf` and `b2_4`'s prereg test are prose checks — **rewrite** to assert the specific claim or delete. |
| `test_kinetic_core_b3.py` | 45 | keep, label | Acrylamide lane. 7 shape-only; `no_dft`, `no_fructose_specific_parameter_exists` are policy guards (keep); `the_holdout_scorer_declares_the_fructose_rows_low_confidence` and `the_fit_script_reads_no_holdout_file` read generator text — **merge** into the freeze guard. |
| `test_kinetic_core_b4.py` | 53 | keep | Observable layer (OAV, thresholds, intervals). 9 shape-only are "no X is tabulated" policy guards; keep. |
| `test_kinetic_core_b5_cutover.py` | 23 | keep | The front door's refusal contract and the frozen exam schema. `the_prereg_exists…` — prose check, **rewrite** to assert the two families it names are still refused. |
| `test_kinetic_core_b6.py` | 41 | keep | Lipid lane; interval widths; the disclosure. 6 shape-only are policy guards. |
| `test_kinetic_core_b7.py` | 42 | keep | Furanic channel; balance; the B7 fit round trip. |
| `test_kinetic_core_b8.py` | 27 | keep, label | Final sulfur wave; Amendment 16–18 guards. Section 5 removed at B5b (safety.py). |
| `test_kinetic_core_q1_keyspaces.py` | 13 | keep | The key-space contract that hid two report bugs; high value. |
| `test_kinetic_core_scoring.py` | 13 | keep | B3 scorecard contract; the fit-target record check now reads the generated record. |
| `test_kinetic_core_uncertainty.py` | 25 | keep | Envelope contract incl. the B8 covariance partition and headspace-band gating. |
| `test_user_scoring.py` | 8 | keep | Step 6. |

## Front door, reports, guards

| file | tests | verdict | note |
| --- | --- | --- | --- |
| `test_comparative_cli_2026_08.py` | 29 | keep | Reliability wiring reads the directional artifact; CLI end-to-end via subprocess. |
| `test_v1_reports.py` | 20 | keep, rewrite 5 | The 5 shape-only tests (`report_carries_the_model_version…`, `every_absolute…with_an_interval`, …) check for substrings in HTML — **rewrite** to parse the table and assert the interval column is present for every absolute row. |
| `test_cli_scripts.py` | 2 | keep | Smoke: every surviving entry point parses `--help`. |
| `test_core_headline_guards.py` | 11 | keep | Pins every core headline and recomputes the scorecard live; the freshness mechanism. |
| `test_data_headline_guards.py` | 1 | keep | The provenance census. |
| `test_wave_generators_frozen.py` | 2 | keep | Step 3. |
| `test_wave_u_maillard_path_holdout.py` | 4 | keep | Hold-out hygiene; count re-pinned 16. |

## Data integrity (curated inputs)

| file | tests | verdict | note |
| --- | --- | --- | --- |
| `test_compound_keys.py`, `test_paper_keys.py`, `test_chemistry_family_keys.py`, `test_data_readme.py`, `test_data_path_strings_resolve.py` | 6 / 4 / 3 / 1 / 1 | keep | The keys layer and the ghost/path checks; each fails on a concrete data drift. |
| `test_rdkit_logic.py` | 2 | keep | InChIKey resolution. |
| `test_public_api_and_sugar_classification_smoke.py` | 1 | **delete** | Imports the five public names; `test_cli_scripts` and the core tests already import them. Zero information. |

## Literature side (registries and their markdown reports)

| file | tests | verdict | note |
| --- | --- | --- | --- |
| `test_family_ingestion_plan.py`, `test_chemistry_family_scope.py`, `test_family_strategy_policy.py`, `test_matrix_family_coverage.py`, `test_matrix_family_next_action.py`, `test_matrix_family_priority_ranking.py`, `test_mycoprotein_reference.py`, `test_scope_gap_guard.py`, `test_dha_lysinoalanine_external_package.py`, `test_extrusion_external_closure.py`, `test_structural_unlock_triage.py`, `test_deep_research_runtime_queue.py`, `test_literature_backlog.py`, `test_deep_research_tracker.py` | 2–4 each | **merge + rewrite** | Fourteen files with the same two-test shape: "the artifact builds" and "the markdown mentions phrase X". The second kind (10 functions) asserts prose and fails only if someone edits a sentence — **delete those**; keep one builder test per registry, moved into ONE file `tests/unit/test_literature_registries.py` parametrised over the generators. |
| `test_experiment_value.py` | 19 | keep | The `rank` verb's arithmetic. |
| `test_extrusion_process.py`, `test_extrusion_benchmark_protocol.py` | 17 / 7 | keep | Extrusion process model and the DoE generator; value-asserting. |
| `test_gap_heatmap.py` | 4 | keep | One writes-a-file test is fine as a smoke. |
| `test_lipid_oxidation_runtime.py`, `test_lipid_oxidation_saturation.py` | 2 / 6 | keep | The lipid substrate-limited kinetics used by the extrusion side. |

## Actions executed in this step

1. Deleted `test_public_api_and_sugar_classification_smoke.py`.
2. Deleted the ten "markdown mentions phrase" prose tests across the literature-side files (list in the commit).
3. Labelled the wave contract files (`b2`, `b2_1`–`b2_4`, `b3`, `b8`) as frozen-wave regression records in their
   module docstrings.

## Executed 2026-09-03 (second pass, "go ahead with the backlog")

1. `tests/support.py` holds the two helpers three files had copied: `executable_code` / `strip_prose`
   (tokenize-based, replaces two regex `_strip_prose` and one `_executable_code`) and `wave_generator`, a
   context manager that imports a row-splicing wave generator (B8, B9) and restores B2.3's row table on exit.
   `test_kinetic_core_b8.py` uses it; the wave suites now pass in any order (verified with B8 first).
2. The three per-wave "hold-out scorer contains no optimiser" copies are one parametrised test in
   `test_wave_generators_frozen.py`; the B2.3 optimiser-budget and B3 fructose-prose text checks are gone
   (the hash manifest pins those files whole).
3. `tests/conftest.py`: fifteen dead DFT / Cantera / recommendation-engine fixtures and an auto-marker list
   naming five deleted files removed; the scientific-directory marker rule kept.
4. Ten single-test literature-side files merged into `tests/unit/test_literature_registries.py` with their value
   assertions verbatim plus one parametrised build-and-render smoke test (their renderer imports were dead).
   The three `tests/scientific` literature files and `test_deep_research_tracker.py` stay: they pin registry state.
5. `test_v1_reports.py` reviewed and KEPT: its ~60 substring checks are presentation contracts (CSS classes,
   refusal cards, negative checks for "inf"), not prose; every one has a structured counterpart the same file
   already exercises. Not worth a rewrite.
6. Coverage: `[tool.coverage.run]` in `pyproject.toml` (`concurrency = multiprocessing`, `parallel`, wave
   generators omitted), `coverage` in `environment.yml`, and a `coverage` alias in `docker_maillard.sh`.

## Left for later (backlog, in order)

(Empty as of 2026-09-03. New test findings go to the improvement backlog in `tasks/data_restructure_plan.md` §7.)
