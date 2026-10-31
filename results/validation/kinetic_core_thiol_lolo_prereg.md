# Leave-one-laboratory-out adequacy of the frozen thiol lane (pre-registration)

*Plan step B8 "Adequacy" and roadmap rows 1 and 9 (`docs/guides/STATUS_AND_PLAN.md` sec. 5), judged by the
Acceptance rule (`docs/guides/STRATEGY_2026-10-29.md` sec. 9). Written 2026-10-31 BEFORE any fold was fitted; the only
engine-free work done first was listing the fit table by laboratory (counts in sec. 3). Generator:
`scripts/generators/fit_thiol_lolo.py`. Nothing here changes the engine, `src/` or `data/`.*

## 1. Question and decision

**Question.** Does the admitted thiol structure, with its constants refitted without a laboratory, predict that
laboratory's rows within honest intervals, and better than two chemistry-free null models?

**Decision it feeds.** Plan step B8 / roadmap row 1's stop rule: the thiol lane is frozen at the admitted structure
and shipped with the measured out-of-sample spread as its uncertainty. This test supplies that spread (the interval
width at which coverage is honest) and says whether the mechanism earns anything out of sample (the null
comparison). It does not choose between structures and moves no constant of the shipped engine.

## 2. The frozen structure

- **Arm R3En as admitted in EX-024** (`docs/reference/experiment_register.yml`): the joint deoxyosone + thiol
  structure with the pentose + amine entry, the 2-acetylthiazole step, MFT through the intact C5 skeleton (DDP) with
  no pH law on the DDP-forming step; Strecker H2S share 0.042 (CR-11); inline pH coupling; sub-model 1's switches off.
  Coordinates (16), bounds and declared starts exactly as `fit_joint_deoxyosone_thiol.coordinates(False,
  ARMS["R3En"])`: `mgo_share_nf`, the seven multiplier offsets (`k_dpo_ddp`, `k_dpo_c2c3`, `k_dpo_nf`, `k_nf_mft`,
  `fft_formation`, `k_h2s_loss`, `k_mft_decay`), `offset_k_cys_actz`, `offset_k_aa_frag`, the five entry coordinates
  (`log10_k_pent_amine`, `entry_ph_slope`, `log10_k_dpo_aa`, `aa_ph_slope`, `amine_return`), `offset_ddp_route`.
- **The full optimum**: `results/validation/kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en_fixed.json` (status 2,
  cost 1075.7, reduced chi-square 4.48; 521 rows on 72 pots).
- **The measurement model, fixed**: R3En's own (one level, no `within_tau`): `fit_joint_deoxyosone_thiol.series_of`
  / `method_of`, `fit_deoxyosone_submodel.SCALE_TAU` (sida and calibrated absolute; quinoxaline figures `semi_figure`
  tau 0.3; `zhou_free` and other semi-quantitative series free, tau inf), every series' extra noise FIXED at
  `s042_e_conv`'s estimate (no re-estimation inside a fold), so each fold's objective is R3En's objective restricted to
  the training rows. The B7a two-level model is not used: R3En was admitted under the one-level model.

## 3. The fit table and the laboratories

**Rows.** Exactly the 521 fit rows of the full R3En fit (the generator rebuilds them with the joint fit's own
selection and refuses to run unless the count per data set x observable equals the record's): the deoxyosone and
thiol sub-model tables (deduplicated), Zhou 2023's pH 7 2-acetylthiazole row, and Davidek 2008's acetic acid /
xylose rows (the entry table). Hexose pots are dropped as in the fit. NOT in this test: the declared hold-outs
(every Zhai 140 C pot, Yiltirak 2026, Hofmann 1998 T1 / T2 pH 3 / 7, Cao 2017's pH pot) and the score-only rows
(Kumazawa 2003, Hofmann T8 / thiamine pot, Zhou 2023 pH 6 / 8 2-acetylthiazole, Zhou's endpoint pH): they stay the
separate, once-scored hold-out set.

**Laboratory = one research group**, read from each row's `laboratory` field, which `src/kinetic_core/fit_rows.py`
maps from the ledger through the hand-curated `data/keys/laboratories.yml` (group, not institution). The table holds
six laboratories:

| fold (laboratory key) | papers in the fit table | rows | method classes |
|---|---|---|---|
| `jiangnan_zhang_x` | Zhai 2019, 2021, 2023 Food Chem, 2023 JAFC; Zhou 2022, 2023, 2023c, 2024, 2025b | 454 | semi_figure, calibrated, zhou_free |
| `dalian_polytechnic_chen_d` | Wang 2026b | 19 | semi_figure, zhou_free |
| `nestle_blank` | Davidek 2008 | 19 | sida |
| `dfa_hofmann_schieberle` | Hofmann 1998 (T3, T4, T5, T10), Schieberle 2000 | 16 | sida, calibrated |
| `reading_mottram` | Desclaux 2006 | 8 | semi_figure |
| `ochanomizu_murata` | Mikami 2017 | 5 | calibrated |

**Duplicate-data papers are merged into the ORIGINAL's laboratory** (conflicts register, `duplicate_confirmed`):
Kang 2026 (Zhejiang Gongshang, Fu / Zhai Y) reprints Jiangnan's Zhai 2023 Food Chem data (CF-005), so any Kang row
counts as `jiangnan_zhang_x` (the generator's `COPY_LAB`); no Kang row is in this table. The in-table duplicates
(CF-005: Zhai 2023 JAFC Fig. 3 col. E = Zhai 2023 Food Chem Table 1; CF-066: Zhai 2021 / 2023) and same-laboratory
re-publications (CF-019, CF-006) all sit inside one laboratory, so no copy of a held-out row stays in its own fold's
fit. Hofmann 1998 and Schieberle 2000 are one laboratory (DFA Garching); Wang 2026b is a separate group (Chen D,
Dalian Polytechnic; no duplicate recorded). A unit test checks that no laboratory's rows enter its own fold's fit.

**The imbalance is the result, not a defect to fix.** Jiangnan holds 87 % of the rows; its fold is fitted on 67 rows
of five laboratories and most deoxyosone-partition and 2-acetylthiazole coordinates lose their only data. That fold
says how far the lane transfers to the laboratory that built it, and it dominates every pooled number; every metric is
therefore also reported per fold.

## 4. Per fold

1. **Refit** every one of the 16 coordinates on the other laboratories' rows, under the fixed measurement model of
   sec. 2, warm-started from the full optimum: `fitting.fit`, budget 1200 evaluations, `ftol` 1e-4 (EX-024's
   settings). **Converged** = optimiser status != 0 (ftol / xtol / gtol). A budget-exhausted fit is continued from its
   own end point, up to three continuations; a fold still budget-exhausted then has NO verdict, and the test's
   overall verdict is INCOMPLETE (reported, never read as a pass or a fail).
2. **Identification in the fold**: `fitting.laplace` at the fold optimum on the training rows. A coordinate
   identified and off its bounds in the full fit but unidentified or on a bound in the fold is **leak-prone**: the
   training data do not pin it, so a warm start leaves it at a value the held-out laboratory helped set.
3. **Reset arm** (only for folds with a leak-prone coordinate): refit with every leak-prone coordinate moved to the
   arm's DECLARED start (`coordinates()`), the others warm, same budget and continuation rule. Its predictions are the
   honest ones where the warm fit leaks; both are reported.
4. **Predict** the held-out laboratory's rows from each converged fold fit (warm, and reset where it exists).
5. **Predictive distribution** of a held-out row, per model M (the engine fold fit and each null of sec. 5), using
   only the training laboratories:
   - row error: the row's own sigma plus, in quadrature, M's **pooled series noise** for the row's method class: the
     median over the training series of that class of M's extra noise, estimated on M's training residuals with the
     joint fit's estimator (`noise_with_floor`: within-series scatter, 0.10 dex floor; the 2-acetylthiazole floor
     0.15). A held-out series never contributes to its own noise;
   - **between-laboratory error** (the Acceptance rule's between-paper error per method class, D4: pooled across
     laboratories, no laboratory terms): `fitting.estimate_between_tau` on M's training residuals with the level group
     = paper x species for absolute rows (sida, calibrated); a class with fewer than two training groups takes the
     largest tau estimated in that fold (declared, conservative). Quinoxaline figures keep their declared calibration
     width 0.3 (`semi_figure`); free-scale series stay free (scored on shape, one degree of freedom each);
   - scoring: `fitting.MeasurementModel.log_density` and `.coverage` with these terms.
   Sensitivity (engine only, flagged as in-sample): the held-out series keep `s042_e_conv`'s fixed noise.

## 5. The two null baselines (adapted to fit rows; the panel's own baselines, `null_model_baselines.json` and
`contrast_scores.json`, score panel rows, not these)

- **Typical yield** (`build_null_model_baselines.py`'s definition): per species, the median over the training
  laboratories' rows of y = log10(value / limiting precursor), the limiting precursor being the pot's smallest
  positive initial charge (mmol/L); predicted value = 10^median x the held-out pot's limiting precursor.
- **No change** (`build_contrast_scores.py`'s p = 0): per species, the median over the training rows of log10(value);
  the prediction is the same for every pot, time and condition.
- A row of a species no training laboratory measured has no null answer; it is counted and left out of the paired
  comparison. Each null gets its own noise and between-laboratory terms from its own training residuals (sec. 4.5).

## 6. Metrics (all printed from the run)

Per fold and pooled over folds, on the rows all three models answer with a finite residual (counts of rows left out,
by reason, printed):
- **held-out log10 error** (predicted minus measured, dex): n, mean, RMS and median |r|, raw and after the level the
  measurement model integrates out (free series centred on their weighted mean, absolute rows raw); per method class
  as well;
- **coverage of the central 80 %, 90 % and 95 % predictive intervals**, each with its **median full width (dex)**;
- **held-out log predictive density** (natural log, summed over the fold's rows) for the engine and both nulls, and
  the **paired differences** engine - null per fold; across folds the mean, its standard error (sd / sqrt(n folds))
  and the per-row mean;
- per fold: the fit's status, evaluations, continuations, cost, its Laplace identification, the leak-prone list, and
  every coordinate's shift from the full optimum (dex and full-fit sigmas).
- *Addendum after the smoke run (2-iteration Jiangnan fold, no verdict-bearing fit), before the full run:* the
  Jiangnan fold's 1-deoxypentosone (112 rows) and 2-acetylthiazole (1 row) have no training laboratory, so no null
  answers them and the paired set drops them. The engine's coverage, widths and error on EVERY finite held-out row
  are therefore also reported, as a diagnostic; the gate of sec. 7 is unchanged (paired rows).

## 7. Decision rule (Acceptance rule, STRATEGY sec. 9 items 2 - 3; STATUS sec. 5 B8)

Computed on the pooled held-out rows of the **verdict fits**: the warm fits, and separately the set with each
leak-prone fold's warm fit replaced by its reset fit; a criterion passes only if it passes on BOTH sets.
- (a) **Calibration**: pooled 90 % coverage in [0.80, 0.97] (the rule's gate is the 90 % interval; 80 % and 95 % are
  reported beside it with their widths and are not gates).
- (b) **Beats each null**: mean paired fold difference of held-out log density > 2 se; |mean| < 2 se is
  indistinguishable; mean < -2 se the null wins.
- **ADEQUATE** = (a) and (b) against both nulls. **CALIBRATED, NOT BETTER THAN A NULL** = (a) but not (b) for at
  least one null: ship the thiols with the measured spread (sec. 6 widths) as their uncertainty, with no
  out-of-sample claim for the mechanism. **NOT ADEQUATE** = (a) fails (coverage below 0.80: the lane's intervals are
  overconfident across laboratories; above 0.97: wider than needed). **INCOMPLETE** = a fold without a converged fit.
- With six folds, one of them 87 % of the rows, (b) is a weak test; the per-row differences and the per-fold table
  are read beside it, never instead of it.

## 8. Expectations (written before the run)

| expectation | P |
|---|---|
| every fold converges within the continuation allowance | 0.85 |
| the Jiangnan fold has >= 3 leak-prone coordinates (2-acetylthiazole, partition and DDP-route coordinates lose their data) | 0.85 |
| at least one other fold has a leak-prone coordinate (Davidek pins the entry's DPO-cleavage constants; DFA the norfuraneol -> MFT step) | 0.7 |
| the reset fit predicts Jiangnan worse than the warm fit (held-out log density lower) | 0.75 |
| Jiangnan's held-out raw RMS > 2x its in-sample raw RMS at the full optimum | 0.6 |
| pooled 90 % coverage in [0.80, 0.97] (warm set) | 0.35 |
| pooled 90 % coverage below 0.80 (warm set) | 0.55 |
| engine beats the typical-yield null (b) | 0.35 |
| engine beats the no-change null (b) | 0.45 |
| ADEQUATE | 0.15 |
| CALIBRATED, NOT BETTER THAN A NULL | 0.2 |

## 9. Compute

Six folds, run as separate resumable jobs (`--fold <lab>`; `--all` runs every unfinished fold in its own process,
then the summary). Per non-Jiangnan fold ~65 pots per evaluation (~5 s at 6 workers on the full fit: 128 evaluations
in 660 s), ~50 - 150 evaluations from the warm start + 33 for the Laplace; the Jiangnan fold's training pots are few.
Estimated 2 - 5 h wall clock at nice 10, OMP_NUM_THREADS=1, MAILLARD_WORKERS=6, beside another job. Outputs:
`results/validation/kinetic_core_thiol_lolo_<lab>.json` per fold, `kinetic_core_thiol_lolo.{json,md}` the summary.
