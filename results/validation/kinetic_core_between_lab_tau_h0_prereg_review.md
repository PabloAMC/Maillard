# Independent review: H0 nested-tau prereg (2026-11-01)

Read: prereg.md, README.md, the plan's track H, the premortem, `fitting.estimate_between_tau` / `log_density`, and `fit_thiol_lolo.py` (`training_terms`, `predictive_model`, `row_intervals`, `null_residuals`, `run_fold`). I looked at no numbers.

## Verdicts

**1. Out-of-sample: BLOCKER.** One-deep nesting reuses the scored lab's information, and does not only bias tau low (7b.4 is wrong). Fold L' was trained with L, so r_L' ≈ δ_L' − Σ π_l δ_l, which contains δ_L. Jiangnan holds most of every other fold's training rows. Hence tau_ABS(S) for outer L = Jiangnan ≈ RMS(δ_DFA/Och − δ_J): the same pairwise distance as Jiangnan's own error. The same holds the other way for DFA and Ochanomizu. Coverage of the scored lab is pushed toward nominal by construction, and the nulls (medians trained with L) get the same leak. Strict nesting is cheap: 10 refits taking 3–6 CPU h, and the nulls need none. Make S1 primary for every arm, and one-deep the sensitivity. If S1 cannot reproduce on the changed engine, run it at 7b1bc7cd. Also add a numpy toy of "fit = row-share-weighted mean of lab offsets" on the real design to put a size on the leak. Noise weights, warm starts and structure are leaks the prereg declares, and that is fine. Note that a strict null may not answer some rows, so define the paired mask before scoring.

**2. Estimator: fix.** The omega identity checks out: Var(centred_i) = s_i² − 1/W + c_i ω². The FREE log density is the restricted density, so paired differences are exact. With 2 inner labs and f = 1, taû²/tau² ≈ χ²₂/2. The median is 0.83·tau, and the 5–95 % range is 0.23–1.73·tau. A plug-in Gaussian then gives ~0.76 coverage at best, because a_L/taû ~ t₂. Truncating once after averaging is right. But "two changes only" is false. There are also (iii) omega inside v_g, (iv) sida and calibrated pooled into one tau (the shipped code gives one per method), and (v) the shipped ≥ 2-groups rule, which is now unstated. State these, and adjust the identity unit test to match. Suggestion: add a Student-t predictive (df = number of inner labs) as a declared sensitivity.

**3. Scores and bootstrap: OK, with a fix.** The log density (joint per group) and the 90 % interval score are proper, and every arm gets the same treatment. Replace the seeded bootstrap with exact enumeration of the 126 multisets (multinomial weights), and add an exact sign-flip test (32 patterns). D3 at 2 SE with 4 df is two-sided p ≈ 0.12; say so.

**4. Decision rules: fix.**
- No precedence among the verdicts. Write the order out (e.g. INVALID > INCONCLUSIVE > NOT CALIBRATED > …), and say what happens when the warm and reset sets disagree.
- "LPD and IS disagree in sign → INCONCLUSIVE" will fire on noise when both are null. Apply it only when the LPD criterion passes.
- The ratio tests (UNSTABLE max/min > 2; S1 "> 2x") are undefined or infinite when an estimate truncates to 0. Floor them, e.g. max(tau, 0.1), or use an absolute difference.
- D2's "90 % upper bound" is unspecified (3-lab bootstrap? P1 spread at a 3-point grid?). Define it as test inversion on a fine tau grid (0.1–1.6 dex).
- Report P(UNSTABLE | homogeneous truth) from the controls. If it is ~0.6, D2 is decided by the design, not by the data.

**5. Controls: fix.** The data are simulated from the estimator's own model, so P1 and N1 can hardly fail. N1's median ≤ 0.10 holds trivially, because ~50 % of estimates truncate to 0. Add three things:
- (a) mis-scaled row noise (true 1.5× and 0.67× s_i) at tau = 0. This tests 7b.2's "omega absorbs it", which fails for single-row ABS groups.
- (b) the leak toy from item 1.
- (c) D3 size and power: identical arms, and engine tau 0.8 vs null 1.2. Report P(D3 passes) for both.
Also apply P1 to the nested 2-lab tau, not only to tau(all).

**6. Wording: fix.** Keep PROVISIONAL and the H0′ design, which is a genuinely external test. Change:
- "lower bound" to "biased low by the declared leaks; sampling error (×0.3–1.7) dominates; a point estimate".
- "NOT BETTER THAN A NULL" to "NOT SHOWN BETTER (D3 power = x)".
- "ADEQUATE (for a new laboratory)" to "… on EX-056's 5 labs, nested".
State H0′'s own power: 2 labs give a coverage SE of ~0.2.

**7. Missing.**
- Classes and comparisons: requiring every pass is an intersection–union test, so no multiplicity correction is needed. Do report ABS-only paired differences (3 labs).
- Heavy tails: add a t(4)-likelihood LPD sensitivity and a concentration check on the top 5 % of rows.
- Reset vs warm: say explicitly that the inner labs' residuals come from the same set.
- Reproduction guard: it is good, but it covers paired rows only, and it does not check the D-all rows. The row-level `per_row` checks do not see FIG/ABS grouping. Only the joint `log_density` reproduction catches a changed `series_of` or `level_of`, so name that check as the guard.
- The D2 relative-only branch ships shape predictions that nothing scores. Report coverage and LPD of the centred rows against the nulls.

## Required fixes (blocking)

1. Strict nesting (S1) is primary for engine and nulls; one-deep is a sensitivity; add the leak toy; correct 7b.4.
2. Correct the estimator's "two changes" claim (pooling, omega in v_g, minimum groups), and the identity test.
3. Make the rules unambiguous: verdict precedence, warm/reset disagreement, the sign-disagreement trigger, zero-safe ratios, and the definition of D2's upper bound.
4. Add controls that can fail: mis-scaled noise, and D3 size and power. Apply P1 to the nested tau.

## Suggestions

Exact bootstrap enumeration plus a sign-flip test; a t-predictive and a t-likelihood sensitivity; ABS-only differences; scoring of centred (relative) predictions; the wording changes in item 6; P(UNSTABLE) under the null.

## Re-review of v2

**Blocking fixes from v1:**
1. Strict nesting primary, with the leak toy: **yes** (secs. 3.1–3.5, 4.2, 7b).
2. The estimator's differences and the identity tests: **yes** (sec. 2, six named differences, tests 1–4).
3. Unambiguous rules: **yes**. Precedence, sets combined by precedence, the sign trigger gated on LPD, zero-safe ratios, D2's bound by test inversion, P5 gating UNSTABLE. One nit: "converged = status != 0" counts scipy's status -1 (a failure) as converged; use status > 0.
4. Controls that can fail: **yes** (N5, D3a/D3b, P1 on the nested tau). A tolerance problem in P1 follows below.

**New problems in v2:**
- **N1. T1 tests the wrong null: BLOCKING.** Even under correct strict nesting, outer L's error (≈ δ_L − 0.9κ δ_J for L ≠ Jiangnan) and its calibration tau(S) share δ_J. The reason is that Jiangnan is an inner lab whose own residual carries δ_J, and Jiangnan dominates the fit that predicts L. So corr(|outer error|, tau(S)) is legitimately positive (my guess is ~0.2–0.4 at κ = 1). A genuinely new lab would show the same correlation. As written, T1 would declare the run INVALID. Fix: test invariance to the scored lab's own offset. Perturb δ_L with the other offsets held fixed; tau(S) must be unchanged (exactly, in the toy), or corr(δ_L², tau(S)²) ≤ 0.1. Keep the outer-error correlation as T2-style reporting, and compare it with a synthetic 7th lab.
- **N2. Blinding order: BLOCKING (cheap).** `run_fold` writes held-out `engine_raw`, `log_density` and coverage, and step 3 commits those pair-fit records before the step-4 freeze. Those residuals are the inputs to the nested tau. Fix: freeze the prereg and the controls (with their hash) BEFORE the pair fits are launched. The pair-fit wrapper should then write residuals without computing or printing scores, and its convergence check should read only optimiser fields.
- **N3. P1 tolerance: fix (recommended before freeze).** "±10 % at every configuration" now covers tau = 0.3. There, the noise from v_g and the truncation of omega² (biased up at omega = 0, which pulls tau² down) can exceed 10 %. There are also ~100 checks at about a 3-SE threshold. Make the threshold max(10 %, 3 Monte-Carlo SE), and keep the hard criterion at tau ≥ 0.8 (report 0.3), or declare a family-wise allowance.

**Verdict: REVISE.** Fix N1 and N2 (and ideally N3); then the design can be frozen without a further full review.

## Re-review of v3

- **T1 null: resolved.** T1 now tests invariance to delta_L, exactly, in the toy. T1b is the same canary on the real record layout. The outer-error correlation moved to T2 and is compared with a synthetic 7th lab.
- **Blinding order: resolved.**
  - The freeze (prereg, controls, both scripts) comes before any pair fit.
  - A separate residual-only runner is used; the frozen generator is not modified.
  - Unit tests check that score keys are absent and that stdout carries no residual statistic.
  - Guard 4.2c goes through the same runner; Reading's scores are already public.
- **status > 0: resolved.**
- **P1 tolerance: resolved.** The tolerance is max(10 %, 3 MC SE), hard only at tau >= 0.8. With 144 checks at ~3 SE, about 0.4 failures are expected under no bias, so P(>= 2 failures) is about 0.06. The one-re-run allowance makes that acceptable.
- **New problems:** none blocking. One minor point: the runner writes raw training residuals before scoring. That is fine, because they are inputs, not held-out scores. The scoring step should still be the first code to read the held-out residuals of pair records.

**Verdict: FREEZE.**
