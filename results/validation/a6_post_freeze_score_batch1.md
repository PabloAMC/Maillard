# PREREG A6 (R1) scored on post-freeze batch 1 (STATUS_AND_PLAN sec. 0, step 10)

*Run once under [A6_SCORING_PREREG.md](../../docs/reference/lever_table/A6_SCORING_PREREG.md), revision 2. Training: 1643 observations from sources in no post-freeze batch; test: 82 pH observations from batch 1's sources.*

## Verdict: too few to score

Eligibility (>= 20 contrasts, >= 4 sources, >= 3 laboratories, no source above 40%): 13 contrasts, 2 sources, 2 laboratories, largest source share 92%. By source: wang2026arp 1, yang2011 12.
Contrasts dropped because every training row of their cell has one control class (R1 = baseline): 4.

No prediction was computed: the protocol decides eligibility from row counts first.

R2 (final pH) is not scored as a predictive rule, and R3 cannot be scored on this batch (A6_SCORING_PREREG sec. 1). A verdict changes nothing in the table.
