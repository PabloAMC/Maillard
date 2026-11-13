# Scoring PREREG amendment A6 on post-freeze sources (STATUS_AND_PLAN sec. 0, step 10)

*Revision 2, 2026-11-13, frozen by the commit that adds it, BEFORE the scorer is run; the scorer runs once per
post-freeze batch. Revision 1 was premortemed by an independent reviewer (ADOPT WITH EDITS; all seven edits are in).
A6 itself is in `PREREG.md`, amendments, frozen at 3b2f80ad. No rule's prediction on any post-freeze source has
been computed. The keyed values and the table's readings of batch 1 (the ten sources keyed in step 9, `new_labs.yml`,
third registration) HAVE been seen, so this tests the rules' predictions, not a blind look at the data.*

## 1. The question, and what this batch can answer

Does banding a pH cell by its keyed control state (A6-R1) predict pH contrasts from sources keyed after the freeze
better than the table's own rule, and better than an arbitrary split of the same sizes?

**Stated in advance (premortem):** batch 1 is expected to be **too few to score**. Counting rows only (no fit), R1
has 17 eligible contrasts from 2 sources, 15 of them from one paper, and all buffered. This protocol is frozen anyway,
so batch 1 gets a formal record and batch 2 is scored under exactly the same rules.

**R2 (final pH) is not scored as a predictive rule.** A6-R2 says final pH is a covariate and its slope an alternative
reading, never the step. Its comparison would be across two targets, slope per final-pH unit against slope per set-pH
unit, whose units and variances differ. So an interval-score difference between them measures units, not the rule.
Its training side is also thin: 59 of 306 training rows carry a final-pH slope. R2 stays an audit column
(`lever_field_readings`). **R3** cannot be scored on batch 1: none of its 49 categorical pairs differs in base mixture.

## 2. Data

- **Training:** every lever-table observation from a source that is not in any post-freeze batch.
- **Test:** every pH observation from the batch's sources, as the regenerated table holds it.
- **The cell** of a test contrast is (odorant, pH, band). `ph_control` comes from `build_lever_field_readings.py`
  (the keyed contrast field first).
- **New-laboratory predictions only:** for each test contrast, the training rows of its own laboratory are removed
  from its cell before any fit, in both arms.

## 3. Predictions

Both arms use the table's numeric model (`build_lever_table.numeric_fit`, prior scale 0.5). The prediction is the
fit's new-laboratory 90 % interval, widened by +/- 1.645 sqrt(v) for the contrast's own measurement variance v.

- **Baseline:** the cell's training rows.
- **R1:** the cell's training rows with the test contrast's control class. The classes are `buffered`, or
  `unbuffered` (`adjusted_unbuffered` + `unadjusted`); `unstated` gets no prediction.
- **A fit needs at least 2 rows from at least 2 laboratories.**
- **Scored set:** the contrasts on which both arms predict AND whose R1 rows differ from the full cell's rows. A cell
  where every training row has one class gives R1 = baseline and is dropped; the count dropped is reported.

## 4. Eligibility, score and placebo

- **Eligible only if** the scored set has at least 20 contrasts, at least 4 sources and at least 3 laboratories, and
  no source holds more than 40 % of it. Otherwise the verdict is **"too few to score"**, and no prediction is
  reported. The scorer computes eligibility from row counts before any fit.
- **Primary score:** the interval score, IS (Gneiting & Raftery 2007, alpha = 0.1):
  IS = (u - l) + 20 (l - y) 1[y < l] + 20 (y - u) 1[y > u], where l and u are the widened interval's ends and y is
  the contrast's slope. It is averaged within each source, then over sources with equal weight; lower is better. The
  contrast-weighted mean, coverage and median width are reported beside it.
- **Uncertainty:** a paired bootstrap over sources (1,000 resamples, seed 20261113) gives a 90 % interval for the
  source-weighted difference. It is the main statement of uncertainty.
- **Placebo:** within each cell, the training rows' control labels are reassigned, keeping the cell's class counts.
  - If the number of distinct joint labellings over the scored cells is at most 1,000, every labelling is enumerated
    exactly; otherwise 1,000 random draws are taken (seed 20261113).
  - Fewer than 20 distinct labellings makes the placebo **uninformative**, and the verdict says so.

## 5. Verdicts (fixed now)

- **too few to score:** sec. 4's eligibility fails.
- **placebo uninformative:** eligible, but fewer than 20 distinct labellings.
- **R1 supported:** the source-weighted mean IS(R1) is below the baseline's, the difference is below the placebo's
  5th percentile, and R1's coverage is at least 0.80. The coverage floor is a sanity check only; with few sources
  its own uncertainty is about +/- 0.15.
- **R1 not supported:** otherwise, marked "worse than baseline" if the difference is positive.
- **What a verdict changes:** nothing in the table. A supported rule needs a second eligible batch under this
  protocol before it is proposed as a main rule in a dated amendment. A rule not supported on an eligible batch
  leaves the candidate list.

## 6. Batch 2

The next post-freeze sources should fill what batch 1 lacks:
- **pH studies from laboratories not yet in the scored cells**, unbuffered ones first (batch 1's eligible contrasts
  are all buffered);
- **studies that print the final pH**, so that R2's audit column grows.

Candidates and the download list are kept with the plan (STATUS_AND_PLAN sec. 0), not here.

## 7. Files

- The scorer is `scripts/generators/score_a6_post_freeze.py --batch <n>`; its output is
  `results/validation/a6_post_freeze_score_batch<n>.{json,md}`. Each batch's record is run once and pinned by its
  commit in the freshness gate.
- The scorer is tested on synthetic data only (`tests/unit/test_a6_score.py`).
