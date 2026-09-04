# The fit / hold-out waves are frozen

> Key to the identifiers (B1..B9, the lettered audit waves, Amendment n): `docs/guides/GLOSSARY.md`, Part 3.

`generate_kinetic_core_b*_{fit,holdout,reports,scorers}.py` and `probe_amine_fate_b2_4.py` produced every
frozen parameter of the kinetic core (`results/validation/kinetic_core_b*_fit_report.json`) and every
pre-registered hold-out score beside it. Since 2026-09-03 they are **frozen**:

- they are not re-run by any command, lane or test (`coverage` reports them at 0 % on purpose);
- their SHA-256 hashes are recorded in `results/validation/wave_generators_manifest.json` and
  `tests/scientific/test_wave_generators_frozen.py` fails when a file no longer matches;
- **a change to any of them is a new wave**, not an edit: copy, rename to the next wave id, pre-register
  (`results/validation/kinetic_core_<wave>_prereg.md`), run, freeze, then rebuild the manifest with
  `python scripts/generators/build_wave_manifest.py` and add a line below.

Three non-wave derivatives import them and fit nothing: `generate_kinetic_core_b8_laplace.py` (a covariance at the
frozen B8 optimum), `generate_kinetic_core_b8_profile.py` (slice profiles around it) and `generate_core_fit_targets.py`
(reads the row tables into the fit-target index's vocabulary).

Known dead paths inside frozen files (recorded 2026-09-03, NOT edited because a change is a new wave):
`generate_kinetic_core_b8_reports.py:91` and `generate_kinetic_core_b2_4_scorers.py:79` import
`generate_cutover_final_exam` and `generate_kinetic_core_b6_holdout.py:247` shells out to it; that script was deleted at
B5b, so their `run_exam()` paths raise if called. The exam artifact `cutover_final_exam.{json,md}` is frozen history.

| date | change | why |
| --- | --- | --- |
| 2026-09-03 | manifest created after annotating the sulfur fit rows in `generate_kinetic_core_b2_3_fit.py` with `benchmark_id` / `benchmark_compound` (no numeric change) | step 1 of the post-retirement plan: fit rows declare their bundles |
| 2026-09-03 | **B9** `generate_kinetic_core_b9_fit.py`: B8's objective minus the eight Hofmann 1998 Table 1 LEVEL rows (54 rows), same free set, bands, weighting and protocol; pre-registered in `results/validation/kinetic_core_b9_prereg.md` | the owner's rule: primary evidence fits, end-to-end levels validate |
| 2026-09-03 | **B9 result**: both starts converge (cost 18.74 on 54 rows; B8's vector scores 18.87 there); ships per its prereg (hold-out 4/30 not worse; the four returned Hofmann bundles 1/8 within 3x). Finding: glucose and fructose MFT predict zero without the level rows. Laplace 20/23 identified (chi2_red 1.21); slice profile 4 quadratic / 9 asymmetric / 3 flat / 7 bound-limited. Envelope 10/44 literature rows, 10/43 out of sample. Active bounds: k_dimer_mft, k_dimer_fft, Ea_decay_thiol_sink (ceiling), acid yield (floor). Manifest rebuilt (25 files). | the engine reads `kinetic_core_b9_fit_report.json` |

