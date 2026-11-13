# PREREG A6 (R1) scored on post-freeze batch 2 (STATUS_AND_PLAN sec. 0, step 10)

*Run once under [A6_SCORING_PREREG.md](../../docs/reference/lever_table/A6_SCORING_PREREG.md), revision 2. Training: 1643 observations from sources in no post-freeze batch; test: 235 pH observations from batch 2's sources.*

## Verdict: too few to score

Eligibility (>= 20 contrasts, >= 4 sources, >= 3 laboratories, no source above 40%): 10 contrasts, 3 sources, 3 laboratories, largest source share 80%. By source: hou2017ph 1, laemont2023 1, vanba2013 8.
Contrasts dropped because every training row of their cell has one control class (R1 = baseline): 12.

No prediction was computed: the protocol decides eligibility from row counts first.

R2 (final pH) is not scored as a predictive rule, and R3 cannot be scored on this batch (A6_SCORING_PREREG sec. 1). A verdict changes nothing in the table.
