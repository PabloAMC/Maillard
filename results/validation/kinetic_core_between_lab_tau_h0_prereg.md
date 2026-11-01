# H0 pre-registration: out-of-sample between-laboratory tau for the thiol lane, from EX-056 (v3)

*Inference track H, stage H0 (`docs/guides/STATUS_AND_PLAN.md`, "Inference track H"; premortem
`scratchpad/inference_track_premortem.md` items 1 - 2).*

- Versions: v1 written 2026-11-01; v2 the same day, after the independent review (`review.md`); v3 the same day,
  after the re-review of v2 (same file). v1 and v2 are kept as `prereg_v{1,2}_superseded.md`.
- The revision log is sec. 13.
- Written BEFORE any H0 number is computed. No repo file is changed by this document.

## 0. What the author had seen when writing this

- EX-056's report and summary JSON (`results/validation/kinetic_core_thiol_lolo.{md,json}`, records at commit
  7b1bc7cd, generator committed in 2e031fdb): the pooled and per-fold coverages, the RMS per method class, and the
  in-sample training taus.
- The SCHEMA of the six per-fold records, and one `per_row` entry.
- No nested tau, no lab-weighted score, and no per-row residual table was computed or read.
- Forking-path risk: the design was chosen knowing the aggregate result (NOT ADEQUATE; held-out error ~1 dex). This
  is declared, and it is one reason the verdict is PROVISIONAL (sec. 9).

## 1. Question and estimand

**Question.** How wide must the thiol lane's (R3En) prediction interval be for a laboratory it has never seen? And,
at that width, has it been shown to beat the two chemistry-free nulls?

**Observable classes** (R3En's measurement model: `fit_deoxyosone_submodel.SCALE_TAU`,
`fit_joint_deoxyosone_thiol.method_of`). EX-056's paired held-out rows fall into these classes:

| class | methods | level term (shipped) | H0 estimand | labs with paired held-out rows |
|---|---|---|---|---|
| ABS (absolute) | `sida`, `calibrated` | per paper x species, N(0, tau^2), one in-sample tau per method | **tau_ABS** (dex): RMS level offset of a paper x species group in a NEW lab; one tau for both methods | DFA (16), Jiangnan (81), Ochanomizu (5) |
| FIG (quinoxaline figures) | `semi_figure` | per series, declared width 0.30 | **tau_FIG**: the same for a figure series, floored at 0.30 | Dalian (4), Jiangnan (123), Reading (8) |
| FREE (relative) | `zhou_free` | per series, free (shape only) | no level term | Dalian (12), Jiangnan (137) |
| all classes | | none | **omega** (dex): a new lab's extra row-level (shape) error inside a group, added in quadrature | all 5 |

- **tau_ABS is the decision quantity** (rule D2). omega and tau_FIG complete arm (b)'s predictive distribution; they
  are reported, and the ship rule does not read them.
- tau is an RMS level error, not a pure variance: it includes any bias common to all labs. A zero-centred interval
  for a new lab must carry that bias.
- **The unit is the laboratory** (EX-056's `lab_of`; Kang 2026 counts as Jiangnan).
  - Five labs are scorable: Dalian, DFA, Jiangnan, Ochanomizu and Reading.
  - Nestle has no paired row. It stays in every training set, is never an outer lab, and enters only diagnostic D-all.
- **Arithmetic limit:** at most 3 labs inform each class, so each outer lab's nested tau rests on exactly 2 inner
  labs. With a lab-level truth (f = 1), tau_hat^2/tau^2 ~ chi^2_2/2: the median is 0.83 tau and the 5 - 95 % range is
  0.23 - 1.73 tau. A Gaussian plug-in interval then covers at most ~0.76, because a_L / tau_hat ~ t_2 (review item 2).
  This is the binding limit of H0.

## 2. The estimator (fixed)

Notation, per model M (engine or null) and per fold F (the fit that produced the residual):
- r_i: the raw held-out residual (predicted minus measured, dex);
- s_i^2 = sigma_row,i^2 + n_{M,F}(class_i)^2: the row's sigma plus F's pooled class noise for M (EX-056 sec. 4.5,
  `training_terms` on F's training residuals);
- groups g: the measurement model's level groups. ABS uses `level_of` (paper x species); FIG and FREE use
  `series_of`.

**Step 1: omega, from within-group shape (the level cancels whatever its prior).**
- For every group with >= 2 finite rows, with w_i = 1/s_i^2 and W = sum w:
  - centred_i = r_i - sum(w r)/W;
  - e_i = centred_i^2 - (s_i^2 - 1/W);
  - c_i = 1 - 2 w_i/W + sum_j w_j^2/W^2.
- Then E[e_i] = c_i omega^2 when Var(error_i) = s_i^2 + omega^2. (s_i^2 - 1/W is the shipped centred variance of
  `row_intervals`.)
- Lab estimate: omega_l^2 = sum_{i in l} e_i / sum_{i in l} c_i.
- Over a lab set S: omega^2(S) = max(0, mean_{l in S} omega_l^2). Labs are weighted equally; the result is truncated
  once, after averaging.

**Step 2: tau_c for c in {ABS, FIG}, by moments on lab-weighted group means.**
- With s~_i^2 = s_i^2 + omega^2(S), per group g of class c:
  - m_g = the precision-weighted mean of r;
  - v_g = 1 / sum(1/s~^2);
  - d_g = m_g^2 - v_g.
- Per lab: dbar_l = mean_g d_g.
- tau_c^2(S) = max(0, mean_{l in S_c} dbar_l), where S_c is the labs in S that hold class c groups.
- tau_FIG,used = max(0.30, tau_FIG).

**Minimum-data rule** (it replaces the shipped ">= 2 groups per class"; the fallbacks are recorded):
- A class needs >= 2 groups in S. A single-row ABS group counts.
- With 2+ labs in S_c the estimate is normal.
- With 1 lab in S_c it is a "single-lab tau", flagged.
- With 0 labs, or < 2 groups, the class takes EX-056's shipped in-sample tau for that fold, flagged.
- By design every outer lab has exactly 2 inner labs for ABS and for FIG. The fallbacks can arise only from strict
  masks (sec. 3.2).

**How this differs from the shipped `fitting.estimate_between_tau`**, which estimates tau^2 = max(0, mean_g(m_g^2 -
v_g)) over groups, per method, on in-sample residuals, with no tau for a method with < 2 groups:
1. The residuals are held-out (sec. 3), not in-sample.
2. Labs are weighted equally, then groups equally within a lab; the shipped estimator weights groups equally.
3. omega^2 enters v_g and the weights; the shipped estimator has no omega term.
4. sida and calibrated are pooled into one tau_ABS. The shipped estimator gives one per method; here sida is held out
   only at DFA, so a separate nested sida tau would have 0 inner labs.
5. The minimum-data rule above replaces the shipped >= 2-groups rule.
6. FIG gets an estimated level width, floored at 0.30; the shipped width is a fixed 0.30.

**Why the method of moments, and not REML.**
- It is closed form and unbiased before truncation, under the model, for any fixed weights. There is no optimiser and
  no boundary pathology.
- REML on a laboratory random effect is degenerate with 2 - 3 labs (P(estimate = 0) is high), and the row-weighted
  version gives Jiangnan (341 of 383 rows) the answer.
- Equal lab weights target the next laboratory, not the next row. The price is variance: Ochanomizu's 5 rows weigh as
  much as Jiangnan's 341.
- Sensitivities, reported only:
  - (i) a lab-weighted Paule-Mandel root, mean_l mean_g m_g^2/(v_g + tau^2) = 1;
  - (ii) the row-weighted estimator;
  - (iii) per-method sida / calibrated taus where defined.

**Identity unit tests** (written before the real run):
1. With omega = 0, one lab, one method and every group of that method, the H0 function returns the shipped
   `estimate_between_tau` for that method to 1e-12.
2. With omega = 0 and one lab holding sida and calibrated groups, the untruncated pooled d-bar equals
   (n_sida d-bar_sida + n_cal d-bar_cal)/(n_sida + n_cal), where n counts groups.
3. E[e_i] = c_i omega^2, checked by simulation (10^5 draws, within 3 Monte-Carlo SE).
4. With every lab the same size, the lab-weighted estimate equals the group-weighted one.

## 3. Nesting

1. **Primary: strict nesting, for the engine and both nulls.** For each outer lab L and each inner lab L' != L (both
   scorable), the inner residuals of L' come from a model trained without L AND L'.
   - **Engine:** 10 leave-two-labs-out refits, one per unordered pair of scorable labs. Each fit serves both nests (L
     outer, L' inner, and the reverse).
   - **Settings:** EX-056's exactly. R3En's 16 coordinates; s042_e_conv's fixed series noise; warm start from
     s042_r3en_fixed; budget 1200; ftol 1e-4; up to 3 continuations; a Laplace on the training rows; the same
     leak-prone rule and a reset arm where it fires. The pair fit's class noise comes from `training_terms`, computed
     at the scoring step (sec. 4.3).
     **Converged = optimiser status > 0** (scipy's ftol / xtol / gtol). Status 0 is a budget exhausted after its
     continuations; status -1 is a failure. Neither counts. (EX-056 wrote "!= 0", and its saved statuses are all 2 or
     3.) A pair fit that does not converge leaves its two nests without strict inner
     residuals; those nests take the minimum-data fallback, flagged, and the verdict is capped at INCONCLUSIVE.
   - **Nulls:** `null_residuals` with train = every lab except L and L' (no refit).
   - L is then scored with EX-056's fold-L residuals (the fit without L), unchanged.
2. **Masks, fixed before scoring.**
   - The inner calibration rows of L' for outer L: L''s rows that are finite under the strict engine fit (L, L') AND
     both strict nulls (L, L'). A strict null can lose a species that only L and L' measured.
   - The outer scored rows of L: EX-056's paired mask for L.
   - All three models use the same masks.
3. **Sets.** EX-056's two sets are 'warm' and 'reset_where_leak_prone'. Within a set, the inner residuals come from
   the SAME set: in the reset set, a pair fit's reset variant is used where it exists, else its warm one. A rule is
   evaluated per set; the two are combined by sec. 6's precedence.
4. **Sensitivity: one-deep nesting.** The inner residuals of L' are taken from EX-056's fold L', which was trained
   WITH L. This needs no refit. It is not a one-sided bias, and that is why it is not primary:
   - fold L''s fitted constants absorb roughly the training labs' row-share-weighted offsets, so
     r_L' ~ delta_L' - sum_l pi_l delta_l, a sum that includes pi_L delta_L;
   - Jiangnan holds ~90 % of every other fold's training rows, so for outer L = Jiangnan, tau(S) ~
     RMS(delta_DFA/Och - delta_J). That grows with Jiangnan's own error, and pushes its coverage toward nominal by
     construction;
   - the reverse holds for DFA and Ochanomizu, and the nulls (medians trained with L) leak the same way.
   - One-deep is reported beside strict. A difference |tau_onedeep - tau_strict| > 0.3 dex for any outer lab is
     printed as the measured leak.
5. **Leak toy (numpy, on the real design; runs at step 7.2).**
   - Lab offsets delta_l ~ N(0, tau^2) per class for the six labs.
   - A "fit" trained on lab set T predicts an offset of kappa sum_{l in T} pi_l delta_l, where pi_l is l's share of
     T's rows (the real counts: Jiangnan 454, Dalian 19, DFA 16, Nestle 19, Ochanomizu 5, Reading 8) and kappa in {0,
     0.5, 1} is how much of the offsets the constants absorb.
   - Row noise as in the design. tau in {0.3, 0.8, 1.2}.
   - A synthetic 7th lab is added: never in any training set, ABS rows only, with the median scorable lab's row and
     group counts, predicted by the all-six-labs fit. It is the baseline for the correlation a genuinely new lab
     legitimately shows (T2).
   - Both nestings are run through the same code path as the real data. Criteria in sec. 8 (T1, T1b, T2).
6. **Ship-rule quantity** (not nested): tau_ABS(all) and omega(all), from EX-056's single-out residuals of all five
   scorable labs (in effect the 3 ABS labs). Each lab's residual is already a genuine new-lab error: the fit that
   predicted it never saw that lab. Nesting is for scoring only. The ship rule reads the larger value of the two
   sets.
7. **Instability flag** (the premortem's early sign). Recompute tau_ABS(all) with each ABS lab dropped. Ratios are
   zero-safe, t~ = max(tau, 0.10 dex). UNSTABLE if max t~ / min t~ > 2, or if a drop-one estimate lands on the other
   side of 0.8 dex. Whether UNSTABLE carries any consequence is settled by the controls (P5, D2).

## 4. Reproduction (the guards)

1. **Post-processing guard**, run before anything else. Rebuild from the reloaded rows and the saved `terms`, to a
   relative 1e-9:
   - for every model x fold x set: the per-row (centred, var), the coverage counts, and the joint `log_density` on
     the paired rows;
   - the `engine_*_all_rows` `log_density` and coverage (the D-all rows).
   The JOINT `log_density` is the check that catches a changed `series_of` / `level_of` (FIG / ABS grouping): the
   per-row values alone do not see grouping. Any mismatch stops the run. The curated inputs' sha256 must equal the
   EX-056 provenance; if they do not, read them at 7b1bc7cd.
2. **Engine guard for the 10 refits.** Try the current engine (HEAD of `thiol-lolo-2026-11-01`; `fit_joint_deoxyosone_thiol.py`
   changed after the records, in the A3 phase 1 Config refactor). It must pass all three checks:
   - (a) EX-056 sec. 4.0: the cost at s042_r3en_fixed on all 521 rows equals 1075.66 to a relative 1e-3;
   - (b) re-predicting the held-out rows of two saved folds (Reading warm, DFA warm) from their saved theta reproduces
     the saved `engine_raw` to 1e-6 dex;
   - (c) re-running one single-out fold (Reading, the cheapest, ~7 min) THROUGH THE PAIR-FIT RUNNER of sec. 4.3,
     configured with one held-out lab, reproduces its saved cost to a relative 1e-3 and its held-out `engine_raw` to
     0.01 dex. The runner computes no score, and Reading's EX-056 scores are already published.
   - If any of the three fails, the refits run in a worktree PINNED at 2e031fdb. That commit is 7b1bc7cd plus the
     generator that wrote the EX-056 records; `git diff 7b1bc7cd 2e031fdb -- src scripts/generators/fit_joint_deoxyosone_thiol.py
     scripts/generators/build_submodel_probe.py scripts/generators/fit_deoxyosone_submodel.py data environment.yml`
     is empty (checked 2026-11-01). (a) - (c) are repeated there.
   - If the pinned engine also fails, strict nesting cannot be done. The verdict is then capped at INCONCLUSIVE, with
     one-deep reported as "leaky". D2 still runs, since it does not need nesting.
3. **The pair-fit runner writes residuals only.** It is a new script, `scripts/generators/fit_thiol_lolo_pairs.py`
   (proposed). `fit_thiol_lolo.py` itself, the frozen EX-056 generator, is NOT modified.
   - **What it imports** from `fit_thiol_lolo`: `load_fit_table`, `lab_of`, `engine_state`, `converge`, `_fit_record`,
     and the constants (ARM, FULL_FIT, NOISE_FROM, BUDGET, FTOL, MAX_CONTINUATIONS, REPRODUCTION_TOL).
   - **What it reuses** from `run_fold`: the problem construction, the reproduction check, the warm and reset variants,
     the Laplace on training rows and the leak-prone rule. It is a copy of `run_fold` up to and excluding its "scoring
     (prereg secs. 4.4 - 6)" block.
   - **What it writes** per pair record (`results/validation/kinetic_core_thiol_lolo_pair_<labA>__<labB>.json`):
     - the optimiser fields (status, cost, evaluations, budget, history, theta), the Laplace and the leak-prone list;
     - the engine raw residuals for every row of BOTH held-out labs, per variant (`held_engine_raw`);
     - the raw TRAINING residuals per variant (`train_engine_raw`), which the scoring step needs for
       `training_terms`;
     - the row keys (ledger_id, species, t_min, method).
   - **What it never computes, writes or prints:** null residuals, `training_terms`, `log_density`, coverage,
     per-row intervals, error statistics or paired differences.
   - **Its stdout** carries only optimiser fields: status, cost, evaluations and seconds per segment.
   - **The convergence check** reads only `status`, `budget_exhausted`, `cost` and `evaluations`.
   - **Unit-tested before launch** (`tests/unit/test_thiol_lolo_pairs.py`, proposed): a pair (L, L') yields training =
     data minus both labs; the single-lab case reproduces `make_folds`; the record (on a stub engine) contains none of
     the keys {models, log_density, coverage, per_row, paired_log_density_difference, error, null_raw}; and its stdout
     contains no residual statistic.

## 5. Arms and scores

**Arms** (on EX-056's paired rows, so every arm scores the same rows):

| arm | centre | predictive distribution |
|---|---|---|
| (a) shipped | the fold optimum (Laplace point estimate) | EX-056 sec. 4.5: s_i^2 + in-sample training tau (ABS, per method), 0.30 (FIG), free (FREE) |
| (b) nested (primary) | the same | Gaussian: s_i^2 + omega^2(S); tau_ABS(S); max(0.30, tau_FIG(S)); FREE free; strict inner residuals |
| (c1) typical-yield null, shipped / nested | the null's training medians | the same two treatments, from the null's OWN (strict) residuals |
| (c2) no-change null, shipped / nested | the same | the same |

- "Shipped Laplace" means the fold optimum plus the measurement-model variance, which is what EX-056 scored.
  EX-056's `predictive_model` adds no parameter covariance.
- **Implementation.** `fitting.MeasurementModel` with `extra_noise[series] = sqrt(class_noise^2 + omega^2)`,
  `between_tau = {sida, calibrated: tau_ABS}`, and a `scale_tau` copy with `semi_figure` set to max(0.30, tau_FIG).
  `fit_thiol_lolo.row_intervals` gives the per-row (centred, var).
- **Sensitivities** (each reported for every arm, none verdict-bearing):
  - **St (t predictive).** Each ABS / FIG group's level prior is a Student t with df = the number of inner labs (2),
    scale tau(S). The group's joint density comes from 1-D numerical integration over the level; the intervals come
    from numerical quantiles of the row marginal. This accounts for a scale estimated from 2 labs (sec. 1).
  - **St4 (t(4) likelihood).** Each row's marginal predictive is a t(4) with the same variance (scale^2 = var
    (4 - 2)/4), on the centred residual. It is a per-row marginal, not a joint density, so it is a robustness check
    against heavy tails only.
  - **S2 (parameter term, optional).** Add J_L Sigma_L J_L' from the fold's saved Laplace `sigma` and `correlation`.
    Coordinates with infinite sigma or at a bound are held fixed. It needs ~33 held-out-pot evaluations per fold
    variant: an engine run, no refit.
  - **One-deep nesting** (sec. 3.4).

**Scores** (natural log; dex for widths), per lab and then averaged with equal lab weights:
1. **Held-out log predictive density per row**: LPD_L = `log_density`(L's rows) / n_L. The free groups' improper
   constant drops out identically across arms on the same rows, so paired differences are exact.
2. **Interval score at 90 %** (Gneiting & Raftery 2007) on the centred residual, with h = 1.645 sqrt(var_i):
   IS_i = 2h + 20 max(0, |centred_i| - h). Lower is better.
3. **Coverage** at 80, 90 and 95 % (90 % is the gate), lab-weighted AND pooled by row.
4. **Mean and median full 90 % width.**
5. **ABS-only paired differences** over the 3 ABS labs, reported. Their sign-flip test has 8 patterns, so the minimum
   p is 1/8: they can never be significant, and they say so.
6. **Relative-branch score** (what D2's relative-only branch would ship). Every group with >= 2 rows is centred on
   its weighted mean, with variance s_i^2 - 1/W + c_i omega^2(S). The shape-only LPD (the restricted density) and the
   coverage are reported for the engine against both nulls, with the D3 machinery.
7. **Concentration check.** The share of the lab-weighted paired LPD difference (per-row marginal Gaussian
   contributions) carried by the top 5 % of rows by |contribution|. If it is > 50 %, the difference is flagged "driven
   by few rows".
8. **Diagnostic D-all** (engine only, no null): (a) vs (b) on every finite held-out row, including Jiangnan's 112
   1-deoxypentosone rows and Nestle's 19 rows (6 labs). Coverage and width only.

**Uncertainty and tests.**
- **Exact cluster bootstrap.** All 126 multisets of the 5 outer labs, each weighted by its multinomial probability 5!
  / (prod k_l!) / 5^5. The per-lab scores stay fixed. The output is a 90 % percentile interval; there is no seed. It
  is reported, and it is crude with 5 labs.
- **Exact sign-flip test** on the 5 per-lab paired differences: one-sided p = #{the 32 sign patterns whose mean >=
  the observed} / 32. The minimum p is 1/32.
- **t-type SE** = sd over labs / sqrt(5), with 4 df. A 2-SE rule is two-sided p ~0.12 (one-sided ~0.06).
- **Coverage SE.** At nominal p = 0.9 with rows fully correlated within a lab: sqrt(p(1 - p)/K) = 0.134 for K = 5,
  0.17 for ABS alone (K = 3), 0.21 for H0' (K = 2, sec. 9). EX-056's per-fold 90 % coverages ran from 0.12 to 0.69,
  so 0.13 is the operative floor.
  - **What this means:** a lab-weighted coverage has a 95 % interval of about +-0.26. It cannot tell 0.90 from 0.65,
    so coverage is a guard against gross failure (D1) only. The paired proper-score differences carry the decision,
    and even they have 4 df.

## 6. Decision rules (declared in advance)

- **D1 Coverage guard.** (b)'s 90 % coverage is >= 0.70 lab-weighted AND >= 0.70 pooled by row. The Acceptance-rule
  band [0.80, 0.97] is reported beside it; at SE >= 0.13 it could not be a power-bearing gate.
- **D2 Ship rule for the absolute class.**
  - **The upper bound** is by test inversion on the grid tau0 = 0.10, 0.15, ..., 1.60 dex. At each tau0, 2000
    synthetic data sets on the real design (f = 1, the widest case; omega at omega_hat(all); both sets) give the null
    distribution of tau_ABS(all). The one-sided 95 % upper bound U is the largest tau0 with P(tau_hat <= tau_obs |
    tau0) >= 0.05. The lower bound is the mirror. f = 0 is reported too.
  - (i) If tau_ABS(all) > 0.8 dex, the thiol lane ships **relative within-lab predictions only**: within-study ratios,
    ranks and time-course shapes, with omega in their intervals. Absolute concentrations for a new laboratory are
    refused. tau is biased low by the declared leaks (sec. 7b), sampling error (x0.23 - 1.73) dominates, and the rule
    reads the point estimate.
  - (ii) If tau_ABS(all) <= 0.8 and U > 0.8, absolute predictions ship only with tau_ABS-widened intervals, labelled
    "provisional; tau biased low".
  - (iii) If U <= 0.8, absolute predictions ship with tau_ABS-widened intervals, labelled provisional (sec. 9).
  - **UNSTABLE** (sec. 3.7) moves branch (ii) or (iii) to branch (i), but ONLY if the controls give P(UNSTABLE |
    homogeneous truth, tau = 0.8) <= 0.3. Otherwise the flag is uninformative: it is reported as "fires under a
    homogeneous truth with P = x" and changes nothing.
  - **The mis-scaled-noise fallback** (N5): if N5 fails, branch (i) also requires the multi-row-group tau_ABS(all)
    (single-row groups excluded) > 0.8.
- **D3 Shown better than a null.** Against null k, under (b) for both: the lab-weighted mean paired LPD difference
  (b) - (c_k) is > 2 SE (t-type, 4 df) AND its exact sign-flip one-sided p <= 2/32. Otherwise:
  - |mean| <= 2 SE: "not shown better";
  - mean < -2 SE with sign-flip p <= 2/32 in the null's favour: "the null is better".
  - **Sign-disagreement trigger:** only when the LPD criterion passes. If the interval-score difference then has the
    opposite sign (its point estimate favours the null), that null's comparison is INCONCLUSIVE. When the LPD
    criterion does not pass, the interval score is reported and triggers nothing.
- **Verdict per set, in this precedence** (the first that applies):
  1. **INVALID**: a hard control criterion fails (P1, N1, T1, T1b), or a reproduction guard or a runner unit test
     fails (secs. 4.1, 4.3). There
     is no verdict; the estimator or the code is repaired by amendment before any real number is read.
  2. **INCONCLUSIVE**: any of
     - strict nesting is impossible or incomplete (secs. 3.1, 4.2);
     - P3 fails, so D1 cannot be interpreted;
     - the sign-disagreement trigger fires for a null whose LPD criterion passed.
  3. **NOT CALIBRATED**: D1 fails. The lab-to-lab error is not one shared Gaussian inflation; that points to
     structure, for H1's masking tests.
  4. **CALIBRATED, NOT SHOWN BETTER THAN A NULL (D3 power = x)**: D1 passes, but D3 fails against at least one null.
     x is the D3 power from control D3b. The lane ships with (b)'s widths as its uncertainty and no out-of-sample
     claim for the mechanism. D2 still applies.
  5. **ADEQUATE on EX-056's 5 labs, nested**: D1 passes, and D3 shows (b) better than both nulls.
- **Combining the warm and reset sets:** the reported verdict is whichever set's verdict comes FIRST in the
  precedence above. Their disagreement is printed. A verdict that changes when any one outer lab is dropped (the
  jackknife) is annotated FRAGILE; with 5 labs, one of them 89 % of the rows, that is expected (P ~0.6).
- **Every verdict carries PROVISIONAL** (sec. 9).
- **Reported, not ruled:**
  - (b) vs (a) on LPD and IS (does calibration help at all);
  - omega(all) and tau_FIG(all), with their drop-one ranges;
  - the relative-branch score. If omega(all) > 0.5 dex, the owner is told that relative predictions also carry a
    large shape error. If D2 fires and the relative branch is not shown better than the nulls, the owner is told that
    the branch D2 ships has no evidence of beating a null either.
- **Multiplicity.** Requiring every criterion to pass is an intersection-union test, so no correction is needed.

## 7. Order of operations (blinding)

1. Reproduction guard 4.1, then the design table (groups, rows and s_i per lab x class; no residual touched).
2. The controls of sec. 8 and the leak toy, on synthetic residuals only. Then the unit tests of the pair-fit runner
   (sec. 4.3) and of the estimator (sec. 2).
3. **FREEZE, before any pair fit is launched:** commit this prereg (v3), the control outputs and the two scripts; the
   commit hash is recorded in the pair records' provenance.
4. Engine guard 4.2 (no score is computed in it), then the 10 pair fits through the residual-only runner. Each
   record is committed as it completes. A pair fit is judged on convergence only, never on its residuals or its
   effect on tau; a non-converged fit is handled by sec. 3.1, not re-run with other settings.
5. Only then the scoring step (`score_thiol_lolo_nested_tau.py`), run once: `training_terms`, the strict nulls, the
   nested quantities, every score. Any later change is an amendment with its reason, scored separately.

## 7b. Direction of the remaining biases in tau

- With strict nesting, these leaks remain. Each makes held-out residuals smaller than a genuinely new lab's would be:
  1. the structure R3En was chosen on all six labs (EX-024);
  2. the series noise (s042_e_conv) and the fit-row selection were set on all labs;
  3. the warm starts are the full optimum, which every lab helped set (the reset arm covers only leak-prone
     coordinates);
  4. the null definitions were chosen with all labs in view.
  So **tau is biased low by these declared leaks**. That is a direction, not a bound: sampling error (x0.23 - 1.73 at
  2 - 3 labs) dominates it.
- **One-deep nesting is NOT in this list.** Its leak goes both ways (sec. 3.4), which is why it is a sensitivity.
- **Row noise.** Too-small s_i would push tau UP through v_g. omega absorbs this only for groups with >= 2 rows, not
  for single-row ABS groups. Control N5 measures it, and D2 has a fallback for it.

## 8. Controls (synthetic residuals on the real design; secs. 7.2 and 3.5)

The real design: the same labs, rows, groups, s_i, masks and both sets, through the same code path. R = 1000
replicates per configuration unless stated.
- Truth: b_g = a_l + u_g, with a_l ~ N(0, f tau^2) and u_g ~ N(0, (1 - f) tau^2), so the lab share of the level
  variance is f in {0, 0.5, 1}.
- e_i ~ N(0, s_i^2 + omega^2). FREE series get an arbitrary level.

**Positive control:** tau in {0.3, 0.8, 1.2}, omega in {0, 0.3}, f in {0, 0.5, 1}.

| id | criterion | pass | if it fails |
|---|---|---|---|
| P1 (hard at tau >= 0.8) | mean UNTRUNCATED tau_hat^2 / tau^2, for tau_ABS(all) AND for each outer lab's nested 2-lab tau_ABS(S), per configuration and set | within +- max(10 %, 3 Monte-Carlo SE) at tau in {0.8, 1.2}. Family-wise: 2 tau x 2 omega x 3 f x 6 estimates x 2 sets = 144 checks. One failure is re-run at R = 4000 for that configuration and must then pass; >= 2 failures, or a failure by > 5 MC SE, fails P1. At tau = 0.3 it is reported only (the omega truncation and the v_g noise bias it, review N3) | INVALID |
| P1b (reported) | median and 5 - 95 % of tau_hat / tau, against the chi^2_2/2 benchmark at f = 1 | none | none |
| P2 | the one-sided 95 % test-inversion bound U (sec. 6) is >= the true tau | in >= 93 % of replicates at tau = 0.8, f = 1 and f = 0 (inversion at f = 1 must be conservative) | D2's branch (ii) / (iii) bounds are reported from f = 1 and labelled "anti-conservative" |
| P3 | (b)'s mean lab-weighted 90 % coverage at tau = 0.8, omega = 0.3, f = 0.5 | >= 0.70 (the D1 guard) | D1 cannot be interpreted, so the verdict is capped at INCONCLUSIVE; St's coverage is reported |
| P4 (reported) | the expected coverage shortfall of (b), and of St, per configuration | none | none |
| P5 (reported; used by D2) | P(UNSTABLE) at tau = 0.8 and 1.2, f in {0, 0.5, 1} | none | D2 reads it |

**Negative control:** tau = 0 (no level offset).

| id | criterion | pass | if it fails |
|---|---|---|---|
| N1 (hard) | tau_ABS(all) > 0.3 dex | in <= 10 % of replicates | INVALID |
| N2 | each outer lab's nested tau_ABS(S) > 0.3 | in <= 15 % | flagged |
| N3 (heavy tails) | e_i a scaled Student t(3) with the same variance: tau_ABS(all) > 0.3 | in <= 15 % | the real tau_ABS is flagged "may be inflated by outliers" |
| N4 (reported) | AR(1) rho = 0.5 within a time course | none | none; a shared error within a course loads into tau, which is correct for prediction but is not a lab effect |
| N5 (mis-scaled row noise) | the true row SD is 1.5 x s_i and, separately, 0.67 x s_i, omega_true = 0. At 1.5 x: tau_ABS(all) > 0.3 | in <= 20 % | D2's fallback (multi-row groups only) |
| N5b (reported) | at 0.67 x: (b)'s coverage and width | none | none |

Unlike P1 and N1, N5 is simulated off the estimator's own model and can fail: single-row ABS groups have no within
term to absorb the noise.

**Leak toy (sec. 3.5)**

| id | criterion | pass | if it fails |
|---|---|---|---|
| T1 (hard; invariance to the scored lab's own offset) | in the toy, every draw is regenerated with ONLY delta_L changed: delta_L -> delta_L + D, for D in {-3, -1, +1, +3} dex and per class, with every other offset, every row error and every fit's absorption kappa held fixed. tau(S), omega(S) and tau_FIG(S) for outer L are recomputed through the real strict code path | identical to a relative 1e-12 (exactly invariant under strict nesting: no inner residual contains delta_L), for every outer lab, kappa, tau and set | INVALID (a leak path in the strict code) |
| T1b (hard; the same test on the real pipeline) | canary: on synthetic residual files with the real record layout, add D = +3 dex to every row of lab L wherever it appears in ANY record (pair records, EX-056 fold records), then rerun the scoring step's calibration | tau(S), omega(S) and tau_FIG(S) for outer L unchanged to a relative 1e-12, for every L | INVALID |
| T2 (reported) | per outer lab: corr(\|outer level error\|, tau(S)) under strict and under one-deep, against the synthetic 7th lab's corr; plus the one-deep invariance test of T1 (expected to fail: its size, the change in tau(S) per dex of D, is the leak); plus the one-deep minus strict differences in mean coverage and LPD, per kappa | none | none. A strict corr is legitimate through Jiangnan, which is an inner lab and dominates every fit predicting L, so it is compared with the 7th lab, not with 0 |

**D3 size and power**

| id | criterion | pass | if it fails |
|---|---|---|---|
| D3a (size) | engine and null residuals simulated identically (tau = 0.8, omega = 0.3, f = 0.5) | P(D3 says "better") <= 0.10 | D3's sign-flip threshold drops to 1/32 |
| D3b (power, reported) | engine tau 0.8 vs null tau 1.2, otherwise identical; and engine tau 0.5 vs null 1.2 | none | none; this is x in "NOT SHOWN BETTER (D3 power = x)" |

## 9. The second split

| candidate | what exists | usable now? |
|---|---|---|
| pH-acid sub-model 1 LOLO (`results/validation/kinetic_core_ph_acid_submodel_fit_{A,B,C}.json`, `leave_one_laboratory_out`; btbu_xie, jiangnan_zhang_x, nestle_blank, wageningen_van_boekel) | per ledger x observable n / mean / RMS only | **No.** No per-row residuals, predictive variances, fold noise or nulls. Fits ran at budget/2 with 2 rounds and no convergence record. The observable is pH units. The structure fails both pH controls (CF-084). 2 of its 4 labs are EX-056 labs. It would need a refit with per-row output, and even then it would test the procedure, not tau_ABS. |
| leave-one-PAPER-out (`--lopo` in `fit_joint_deoxyosone_thiol.py`, `fit_deoxyosone_submodel.py`) | code only; the records' `holdout_*` keys are within-lab pot hold-outs | **No.** A paper is not a lab, and nothing is saved. |
| P1 TEST tables keyed 2026-11-01 (`data/lit/fit_rows/heldout_larralde1995.yml`, `heldout_lee2024_fig2.yml`, `heldout_zhang2021_table3.yml`; branch `ingest-2026-11-01` only) | 3 new labs | **No.** None is an ABS-class thiol-lane observable: Zhang 2021's MFT / FFT are semi-quantitative with one internal standard, at one time point, with lipid arms. `fit_rows.load` refuses every row by design. Each table is another consumer's TEST set (E5, E13, E7). |
| repeated grouped 2/3 : 1/3 lab partitions | needs >= 10 labs; the thiol table has 6 | **No.** |

- **H0's verdict is PROVISIONAL.**
- **Upgrade condition (H0').** The first P1 batch that keys >= 2 NEW labs with ABS-class thiol-lane rows triggers
  H0'. Their split hash is frozen first.
  - They are predicted from s042_r3en_fixed, which never saw them: one engine run of their pots, no refit.
  - The intervals use H0's tau_ABS(all), omega(all) and tau_FIG(all) unchanged, scored by secs. 5 - 6.
  - This is a genuinely external test. Its own power is low: with 2 labs the coverage SE is ~0.21, and the sign-flip
    test has 4 patterns, so its minimum p is 1/4.

## 10. Compute and cost

| item | work (h) | compute | refit? |
|---|---|---|---|
| post-processor (`scripts/generators/score_thiol_lolo_nested_tau.py`, proposed), guards, identity tests | 2.5 - 3.5 | seconds | no |
| controls, leak toy, D3 size / power, D2 test-inversion grid (31 tau0 x 2000 x 2 sets) | 1 - 1.5 | ~0.5 h CPU (numpy) | no |
| engine guard (4.2 a - c: one cost evaluation, two re-predictions, one Reading fold) | 0.5 | ~0.2 h | 1 check refit |
| **10 leave-two-labs-out refits** (pair-fit wrapper + unit test, launch, convergence check) | 1 | **~4 - 8 core-hours; ~2 - 4 h wall**, niced, at most 8 cores (2 - 3 pair fits at a time), started when the machine is free | **yes, 10** |
| scoring, report (`results/validation/kinetic_core_thiol_lolo_nested_tau.{json,md}`), independent review | 2 | minutes | no |
| **H0 total** | **7 - 8.5** | **~5 - 9 core-hours; ~3 - 5 h wall** | **10 (+1 check)** |
| S2 parameter term (optional) | 1 | ~0.5 - 1 h CPU | no (engine evaluations only) |

- **Refit-cost basis.** The six EX-056 folds took 384 - 3800 s of wall time each (2.2 h in all), reset arms
  included.
  - The 6 pairs without Jiangnan train on ~490 - 505 rows (~4 s per evaluation; 40 - 460 evaluations + 33 for the
    Laplace).
  - The 4 pairs with Jiangnan train on 51 - 62 rows. Each evaluation is quick, but leak-prone resets are likely there.
- **Environment.** Everything runs through `./scripts/docker_maillard.sh run "..."`, on a branch from
  `thiol-lolo-2026-11-01` (or the 2e031fdb pin, sec. 4.2), after that checkout's test suite finishes. The work done
  while the fits run is the post-processor and the controls; the fits import nothing that is being edited.

## 11. Expectations (written before any H0 number)

| claim | P |
|---|---|
| guard 4.1 passes first time | 0.85 |
| the current engine passes 4.2 (else the 2e031fdb pin) | 0.7 (the pin passes: 0.95) |
| all 10 pair fits converge | 0.75 |
| P1 / N1 / T1 + T1b pass | 0.9 / 0.85 / 0.95 (T1 fails only on an indexing bug) |
| P3 passes (the Gaussian plug-in reaches 0.70) | 0.55 |
| N5 passes | 0.5 |
| P(UNSTABLE \| homogeneous truth) > 0.3, so the flag is uninformative | 0.8 |
| D3b power (engine 0.8 vs null 1.2) < 0.5 | 0.75 |
| tau_ABS(all) > 0.8 dex (D2 branch (i)) | 0.5 |
| (b) passes D1 | 0.45 |
| (b) shown better than both nulls (D3) | 0.05 |
| (b) better than (a) on lab-weighted LPD | 0.9 |
| (b)'s median 90 % width > 2x (a)'s | 0.85 |
| verdict: ADEQUATE / CALIBRATED-NOT-SHOWN-BETTER / NOT CALIBRATED / INCONCLUSIVE / INVALID | 0.05 / 0.30 / 0.25 / 0.32 / 0.08 |
| H0 scored by the 11-07 checkpoint | 0.65 (the 10 refits are the new risk) |

## 12. Threats to validity

The design has a small-K problem that no estimator removes:
- **The data.** Five scorable labs, three per class, two inner labs per nested calibration, and one lab holding 89 %
  of the rows. Every lab-weighted number has 2 - 4 effective degrees of freedom; the exact bootstrap has only 126
  multisets and the sign-flip test only 32 patterns. tau_ABS(all) also rests on correlated residuals: DFA's and
  Ochanomizu's single-out fits both share Jiangnan's offset, so the effective number of labs is below 3.
- **Optimistic biases.** Strict nesting removes the scored lab from its own calibration. The structure, the noise
  model, the row selection and the warm starts were still chosen on all labs, so tau is biased low by those leaks.
  Sampling error dominates that bias.
- **The Gaussian plug-in.** With a scale estimated from 2 labs, its nominal 90 % interval covers ~0.76 at best
  (sec. 1). P3 and St
  measure the shortfall, but do not remove it.
- **The predictive model.** It treats a lab's groups as independent, while a lab shares one effect across its species
  and figures. Every arm's LPD for a many-group lab (Jiangnan) is therefore overconfident in the same way; paired
  differences are distorted less. The lab-level component cannot be separated with <= 3 labs per class (f in the
  controls only bounds the damage).
- **Absorbed structure.** The known defects of EX-056 sec. 7b (the H2S headspace bias, the pot-independent MFT loss,
  the pH allocation CF-084) load into tau as if they were lab noise. A large tau is not evidence that labs differ; H1's
  masking tests must separate structure from lab effects.
- **Nulls.** Their medians follow the training labs' species mix.
- **Author knowledge.** The prereg was written after the aggregate EX-056 result was known (sec. 0).
- **No second split.** Nothing replicates H0 until P1 delivers new ABS-class labs, and H0' itself will be weak.

## 13. Revision log

### v2 -> v3 (answering `review.md`, "Re-review of v2")

| # | change | re-review item |
|---|---|---|
| 22 | T1 now tests INVARIANCE: perturbing only the scored lab's own offset (D in {-3, -1, +1, +3} dex) must leave tau(S), omega(S) and tau_FIG(S) unchanged to a relative 1e-12. T1b is the same canary on the real pipeline. The outer-error correlation moves to T2, against a synthetic 7th lab | N1 (blocking) |
| 23 | freeze (prereg, controls, scripts) BEFORE any pair fit is launched; the order is rewritten (sec. 7) | N2 (blocking) |
| 24 | a residual-only pair-fit runner, a new script that copies `run_fold` up to its scoring block, with the frozen generator unmodified. It writes optimiser fields and raw residuals only, prints optimiser fields only, and is unit-tested for absent score keys. `training_terms` and the nulls move to the scoring step (secs. 3.1, 4.2c, 4.3) | N2 (blocking) |
| 25 | converged = status > 0 (sec. 3.1) | nit under fix 3 |
| 26 | P1's tolerance is max(10 %, 3 MC SE) and hard at tau >= 0.8 only, with a family-wise allowance (one re-run at R = 4000; >= 2 failures or > 5 SE fails); tau = 0.3 is reported only | N3 |

### v1 -> v2 (answering `review.md`)

| # | change | review item |
|---|---|---|
| 1 | Strict nesting (10 leave-two-labs-out refits; nulls retrained without both labs) is PRIMARY for the engine and both nulls; one-deep is a sensitivity (sec. 3.1, 3.4) | 1, blocking fix 1 |
| 2 | v1's 7b.4 claim (one-deep only biases tau low) withdrawn; the two-way leak explained (sec. 3.4); 7b rewritten | 1 |
| 3 | numpy leak toy on the real row shares, with hard criterion T1 and reported T2 (sec. 3.5, sec. 8) | 1, 5(b) |
| 4 | masks defined before scoring; a strict null may not answer a row (sec. 3.2) | 1 |
| 5 | engine guard for the refits: three reproduction checks, the 2e031fdb pin (= 7b1bc7cd + the generator; diff checked empty), and the cap if both fail (sec. 4.2) | 1, blocking fix 1 |
| 6 | refit cost stated: ~4 - 8 core-hours, niced, at most 8 cores, when the machine is free (sec. 10) | blocking fix 1 |
| 7 | "two changes only" corrected to six named differences (omega in v_g, pooled sida + calibrated, the minimum-data rule replacing >= 2 groups, FIG estimated) and the identity tests adjusted (sec. 2) | 2, blocking fix 2 |
| 8 | chi^2_2/2 sampling distribution and the ~0.76 plug-in ceiling stated (sec. 1, 12) | 2 |
| 9 | verdict precedence written out; warm / reset combined by precedence; inner residuals from the same set (secs. 3.3, 6) | 4, 7, blocking fix 3 |
| 10 | the LPD / IS sign trigger fires only when the LPD criterion passes (sec. 6 D3) | 4, blocking fix 3 |
| 11 | zero-safe ratios, max(tau, 0.10); the S1 "> 2x" rule removed, now that strict is primary (sec. 3.7) | 4, blocking fix 3 |
| 12 | D2's upper bound defined by test inversion on a 0.10 - 1.60 dex grid, at f = 1 (sec. 6), with control P2 | 4, blocking fix 3 |
| 13 | P(UNSTABLE \| homogeneous truth) reported (P5); UNSTABLE acts only if <= 0.3 | 4, blocking fix 3 |
| 14 | controls that can fail: mis-scaled row noise at tau = 0 (N5, N5b, D2 fallback); D3 size and power (D3a, D3b); P1 on the untruncated mean, applied to the nested 2-lab tau as well | 5, blocking fix 4 |
| 15 | exact enumeration of the 126 multisets (no seed) and the exact 32-pattern sign-flip test inside D3; the 2-SE rule's p ~0.12 stated | 3 |
| 16 | t predictive (St, df = inner labs) and t(4)-likelihood (St4) sensitivities; top-5 % concentration check | 2, 7 |
| 17 | ABS-only paired differences (3 labs, 8 patterns, never significant) reported | 7 |
| 18 | relative-branch (centred) score against the nulls, reported with D2 | 7 |
| 19 | the reproduction guard names the joint `log_density` as the grouping check and covers the D-all rows (sec. 4.1) | 7 |
| 20 | wording: "lower bound" -> "biased low by the declared leaks; sampling error dominates; point estimate"; "NOT BETTER" -> "NOT SHOWN BETTER (D3 power = x)"; "ADEQUATE on EX-056's 5 labs, nested"; PROVISIONAL kept; H0''s own power stated | 6 |
| 21 | compute and expectations revised (sec. 10, 11); the 11-07 forecast lowered to 0.65 because of the refits | consequence of 1 |
