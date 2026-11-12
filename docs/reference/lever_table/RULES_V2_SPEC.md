# Lever table: sensitivity readings and fragility flags (STATUS_AND_PLAN sec. 0, step 5b-2)

*Revision 2, 2026-11-12, after an independent premortem (ADOPT WITH EDITS) of revision 1. Written AFTER the
all-data table and the disagreement audit were read: nothing here is a blind rule, and nothing here moves a
reading of the table.*

## What changed between revision 1 and this one

Revision 1 proposed a second table with five changed rules, generated beside the first and adopted if its
between-laboratory SD was lower. The premortem showed that this comparison cannot tell "more correct" from "more
split": dropping censored arms removes the extreme values, splitting cells moves conflicting laboratories into
one-laboratory cells that leave the statistic, and the leave-one-laboratory-out check rests on two cells. Every rule
was also chosen after reading the data, TEST laboratories included, so a re-banded DEVELOPMENT-only table would be a
rival informed by TEST contrasts.

So the table keeps its rules and its readings. Each cell gains **sensitivity readings** and a **fragility flag**:
what the same model says under a defensible alternative reading of the same contrasts. A formulator sees which
answers rest on a choice.

## The sensitivity readings (in `build_lever_table.py`)

| reading | what it does | the question it answers |
|---|---|---|
| without censored arms | each contrast is refitted on its measured arms only; a contrast left with fewer than two has no alternative and drops out | does the answer rest on a detection limit used as a value? (PREREG sec. 3 calls these arms bounds; 89 contrasts use one) |
| whole-ladder slopes | a non-monotone ladder contributes its slope over ALL its arms instead of over its longest monotone run | does the answer rest on which stretch of a rise-then-fall ladder was kept? |
| buffered studies only (pH) | the cell is refitted on contrasts whose keyed `buffer` field names a buffer | does the answer rest on unbuffered pots, whose initial pH the reaction did not keep? |
| without censored pairs (categorical) | a level's reading without the pairs that rest on a censored arm | the same as the first row, for categorical levers |

Fitted only where a cell has three laboratories (a reading can only change there). A cell whose reading would
change is marked **[fragile: ...]** with the alternative reading. Also reported per cell, as counts: contrasts with
a censored arm; contrasts with a **supported turning point** (a rising and a falling stretch that each span more
than two standard errors end to end, so that a flat ladder does not turn on noise); for pH, contrasts by buffer
state. A cell with supported turning points in two or more laboratories is marked "depends on the level".

Censoring is informative: a lever that drives an odorant below detection is its strongest effect. The limit used
as a value is a conservative size with the right sign, so the main reading keeps it; the alternative shows what
remains without it. Neither is the censored likelihood, which this Gaussian-grid model does not carry.

## What is NOT done here, and what it needs

| gap seen in the reading | why it is not a rule yet | what it needs |
|---|---|---|
| a "none" arm that is a different pot in each study (9 of the 16 widest size spreads) | the pots can only be told apart from free-text precursor names ("RG-ARP", "TVPH", "soybean peptides", empty), which the two readers did not compare | a per-arm field (`base_amine`, `base_sulfur`) keyed by two readers blind to the values |
| pH bands by buffered against unbuffered, or by final pH | "buffered" is free text with more than 40 spellings and no capacity; 102 contrasts are `ph_basis: label` | a keyed three-way field, and the final pH where the paper prints it |
| the step from zero to the first dose | its size depends on the first dose (0.00125 to 670, mixed units), and many such arms already sit in categorical "none vs X" cells | amounts in one unit; then a list beside the cell, never pooled |
| fed-intermediate and stability experiments in formation cells | not marked in the keyed files | a `system` field at keying |
| a one-sided likelihood for censored arms | not in the model | mu on the grid for numeric cells |

## How a future rule change is tested

No clean test exists on contrasts already read. A rule change is frozen first and scored on sources keyed AFTER the
freeze, one rule at a time, on an identical contrast set and a fixed cell partition, against a placebo (random
bandings with the same cell sizes), per lever class, with interval width beside coverage.

## Step 5a

The rival the engine must beat is the DEVELOPMENT-only table under the rules as they stand (PREREG with A1 - A5),
frozen by hash in step 5a's own pre-registration. The sensitivity readings are additive fields; they change no
estimate of that table. Step 5a opens neither all-data table.
