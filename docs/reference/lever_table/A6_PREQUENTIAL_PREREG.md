# Scoring PREREG amendment A6, revision 3: prequential (STATUS_AND_PLAN sec. 0, step 12)

*Revision 3 of [A6_SCORING_PREREG.md](A6_SCORING_PREREG.md). Revision 2 is frozen at 4d7027aa, stays the record for
batches 1 and 2, and is not edited. Drafted 2026-11-13 and premortemed by an independent reviewer: ADOPT WITH EDITS,
and all ten edits are in (sec. 7). Frozen by the commit that adds it, BEFORE its scorer is run and BEFORE any batch-3
source is chosen. A6 itself (`PREREG.md`, amendments, frozen at 3b2f80ad) is unchanged: this revises how R1 is scored,
not R1. C's rules (secs. 2 - 5) are frozen with this commit. A revision written after any record of this file is read
applies only to batches registered after it, and starts a new confirmation set.*

## 0. Why a revision, and what has been seen

Revision 2 trained on pre-freeze sources only. R1 can only act on a cell whose training rows hold both control classes
from at least two laboratories each, and few pre-freeze pH cells do. So both batches were too few to score
([batch 1](../../results/validation/a6_post_freeze_score_batch1.md): 13 contrasts, 2 sources;
[batch 2](../../results/validation/a6_post_freeze_score_batch2.md): 10 contrasts, 3 sources), and more batches of
the same kind would not have changed that.

**What has been seen:** every keyed value of batches 1 and 2, the lever table's readings of them, and the row counts
in sec. 4. **What has not:** any R1, baseline or placebo prediction for a post-freeze contrast. The first draft of this
revision proposed scoring batches 1 and 2 together because a source-only count said that set was eligible. The
premortem's laboratory cap (sec. 3) makes that set ineligible, and the cap was adopted before the count was rerun.
So the first look gets a formal record of "too few to score", and nothing decides before C.

## 1. The question (unchanged)

Does banding a pH cell by its keyed control state (A6-R1: `buffered` against `unbuffered` =
`adjusted_unbuffered` + `unadjusted`) predict pH contrasts from sources keyed after the freeze better than the table's
own rule, and better than an arbitrary split? R2 and R3 stay unscored, for the reasons in revision 2 sec. 1.

## 2. Data: prequential training

- **Batches.** Every post-freeze source belongs to the batch it was keyed in. The batches are one definition:
  `score_a6_prequential.BATCHES`. Batches 1 and 2 are revision 2's (hou2017ph is batch 2, the batch it was keyed in).
  Pre-freeze sources are batch 0. **Every source keyed after this freeze, for any reason, goes into the next batch;
  none is withheld.**
- **Training for a batch-k contrast:** every lever-table observation of batches 0 to k - 1, and nothing from batch k
  or later. Two kinds of rows leave its cell first:
  - its own laboratory's rows, with laboratories as registered in `split v1` and `new_labs.yml` (from author lines
    and affiliations) at the freeze commit;
  - any row of another source whose slope equals the contrast's within 5e-4. That is a duplicated dataset, such as a
    thesis and its paper. The count is reported.
- **Inputs are pinned by the record's provenance** (hashes of the contrast files, the split, `new_labs.yml` and this
  file). A later correction to a scored batch's values gets a dated note in the plan and does not trigger a rescore.
- Unchanged from revision 2 secs. 2 - 3:
  - the cell is (odorant, pH, band), and `ph_control` comes from `build_lever_field_readings.py`;
  - both arms use `build_lever_table.numeric_fit` (prior scale 0.5) and its new-laboratory 90 % interval, widened by
    +/- 1.645 sqrt(v);
  - a fit needs at least 2 rows from 2 laboratories;
  - a contrast whose R1 rows equal its cell's rows is dropped and counted.

## 3. Two sets, each scored once

- **First look, D: batches 1 and 2.** It is recorded once, when this revision is frozen. It commits nothing (sec. 5).
- **Confirmation, C: batches 3, 4, ... pooled.**
  - After each new batch is keyed, C's eligibility is computed from counts only (no fit) and recorded.
  - C is scored once, at the first batch that makes it eligible, and never again. No prediction on C is computed
    before that.
- **Eligibility, all from counts:**
  - at least 20 contrasts, 4 sources and 3 laboratories;
  - no source and no laboratory above 40 % of the set (the laboratory cap is new);
  - **at least 20 distinct usable placebo labellings** (sec. 5; new). A set that would come back "placebo
    uninformative" is never scored.

## 4. Counts behind this revision (row counts only, 2026-11-13)

- **D:** 41 contrasts (batch 1: 13; batch 2: 28, up from 10 under revision 2 because batch 1 trains it).
  - 8 sources and 7 laboratories, largest source share 29 %.
  - **Largest laboratory share 44 % (btbu_beijing: yang2011 12, hou2017ph 6): ineligible.**
  - 6 contrasts dropped as one-class; 0 duplicate rows. Placebo: 762 of 1,000 drawn labellings usable.
- **Cells where C can act.** In the current table, ten pH cells hold both classes from at least 2 laboratories each:
  - 2-methylpyrazine at pH 6 - 7.5 and > 7.5;
  - furfural at 4.5 - 6 and 6 - 7.5;
  - FFT at 4.5 - 6;
  - trimethylpyrazine at 6 - 7.5;
  - 2,5-dimethylpyrazine and pyrazine at > 7.5;
  - benzaldehyde at 6 - 7.5;
  - phenylacetaldehyde at 4.5 - 6.

  In nine of them the scarcer class has 2 or 3 laboratories, usually `unbuffered`. D's sets of cells add the thiols
  at 4.5 - 6 and 6 - 7.5, through batch 2's rows.
- **Batch-3 selection rule (fixed before any search is run):**
  1. **Search.** A fixed search string is run on a fixed date, both recorded in the plan before running.
  2. **Screen.** Papers are screened in the order the search returns them, on their methods only. A paper is taken if
     its design gives pH contrasts that would land in one of these cells, from a laboratory new to that cell.
     Unbuffered designs come first.
  3. **Take every qualifying paper.** Every paper that qualifies is taken, in search order, until the batch is full.
  4. **Keep a screening log.** It records every paper looked at, and why it was kept or dropped.
  5. **Key the control first.** Each reader keys `ph_control` from the methods before any value is keyed.
- **Spreading laboratories.** A source from a laboratory already above 25 % of C's count is taken only once
  laboratories below that share have been exhausted from the search.

## 5. Score, placebo and verdicts

- **Score:** the interval score (alpha 0.1), averaged within each laboratory and then over laboratories with equal
  weight. The source-weighted, contrast-weighted, coverage and median-width figures are reported beside it.
  **Uncertainty:** a paired bootstrap over laboratories (1,000 resamples, seed 20261113).
- **Placebo:** within each scored cell and each batch of origin, the true control class of whole blocks is reassigned.
  - **Block:** one source's training rows of one class in that cell.
  - **Kept fixed:** the number of buffered blocks in each stratum, so every prequential training set keeps its own
    block counts.
  - **Draws:** exact enumeration up to 1,000 distinct labellings, otherwise 1,000 random draws (seed 20261113).
  - **Unfittable draws:** a draw in which any scored contrast's relabelled R1 subset cannot fit is discarded and
    counted. It is never replaced by the baseline.
- **Verdicts (per set):**
  - **too few to score:** sec. 3 fails;
  - **R1 supported:** the laboratory-weighted mean IS(R1) is below the baseline's and below the usable placebo's
    5th percentile;
  - **R1 not supported:** otherwise, marked "worse than baseline" if the difference is positive.
- **Coverage below 0.80** is reported as a flag, not a verdict condition. With about 8 laboratories its own
  uncertainty is about +/- 0.15, and the interval score already penalises misses.
- **Diagnostics, which change no verdict:**
  - the difference per batch, per cell and per laboratory;
  - a leave-one-cell-out sign count.
- **What the verdicts decide. C decides alone.**
  - **R1 supported on C:** R1 is proposed as a main rule in a dated amendment.
  - **R1 not supported on C:** R1 leaves the candidate list.
  - **D, whatever its record says, commits nothing.** It is reported beside C.
- **Low batch yield.** If the batch-3 search yields fewer than 5 qualifying papers, keying more for A6 is
  deprioritised in the plan (STATUS_AND_PLAN sec. 0). C stays open, and every later pH source joins it (sec. 2).

## 6. Files

- **Scorer:** `scripts/generators/score_a6_prequential.py`, with `--set first_look` (D) or `--set confirmation`
  (C).
- **Outputs:**
  - `results/validation/a6_prequential_score_first_look.{json,md}`;
  - one `a6_prequential_confirmation_through_batch<k>.{json,md}` per batch k >= 3, of which at most one holds
    predictions.

  Each record is run once and pinned by its generator in the freshness gate.
- **Revision 2's scorer,** `score_a6_post_freeze.py`, is not edited, because it generates the batch-1 and batch-2
  records. The new scorer reuses its pieces.
- **Tests:** synthetic data only (`tests/unit/test_a6_prequential.py`).

## 7. Premortem edits folded in

1. **Placebo:** stratified by cell and batch of origin, moving source blocks, with unfittable draws discarded.
2. **Usable labellings** are counted in eligibility.
3. **Laboratories:** a laboratory-weighted score and bootstrap, and a 40 % laboratory cap.
4. **Laboratory identity** is frozen as registered, and duplicated datasets are removed.
5. **Inputs and batch membership** are pinned. hou2017ph is batch 2.
6. **C's rules** are frozen against later revisions. Every later source goes into C.
7. **Batch-3 screening:** a fixed search, screening on methods only, a log, and `ph_control` keyed first.
8. **Decisions:** C decides alone, and D commits nothing.
9. **Per-batch figures and D's selection note** are in the record.
10. **Diagnostics:** per cell and leave-one-cell-out. **The coverage floor** is demoted to a flag.
