# Step 8: the fields the disagreement audit found binding (spec, revision 2)

*2026-11-13. STATUS_AND_PLAN sec. 0, row 8. Literature only; no engine. Revision 2 folds in an independent premortem
of revision 1 (ADOPT WITH EDITS, all six edits taken). Written AFTER the all-data table, the disagreement audit and
RULES_V2_SPEC were read: nothing here moves a reading of the table, and nothing here can be promoted to a main rule
from the sources already read.*

## 1. Why, and what this step can and cannot do

Two of the four answered cells are fragile on the pH basis (FFT x pH, furfural x pH), and three of the nine
conflicts between laboratories turn on facts held only as free text: whether a pH was buffered, what pH the pot
ended at, and whether a "none" arm is the same base pot as its reference. Step 7.0 showed that every later layer
inherits the table's intervals.

**What it cannot do (stated up front):** RULES_V2_SPEC allows no rule to be promoted from contrasts already read, and
no reading adds a laboratory. MFT x pH stays at two laboratories, and the perception layer stays parked. **The
lasting value of this step is a frozen candidate rule and a keying template that records these facts**, so that the
rule can be scored on sources keyed after the freeze. The readings on the current sources are an audit: which flags
change, and why.

**Premortem's main failure stories, and what this revision does about each:**
- the flags do not change (about 55 %), so the cost is kept small;
- a final-pH slope divides by about 0.1 pH unit and produces junk (about 35 %), so final pH is never the slope
  denominator;
- reader cost on free-text base mixtures (about 25 %), so the fields are derived and readers only check them;
- post-hoc readings get treated as evidence (about 20 %), so they are labelled as an audit everywhere.

## 2. The fields: derived, with reader checks only where derivation cannot reach

| field | how it is obtained | readers check |
|---|---|---|
| **F1 `ph_control`** per contrast: `buffered` / `adjusted_unbuffered` / `unadjusted` / `unstated` | a written, deterministic map from the keyed `buffer` text (a named buffer with a concentration is `buffered`; "none (pH adjusted ...)" is `adjusted_unbuffered`; "unstated" is `unstated`) | the sources whose text is a bare "none" or ambiguous: shibamoto1977, ames2001b, wei2019flax, hofmann2000 |
| **F2 `final_ph`** per arm, with its location | arms already keyed `ph_basis: final` | the sources expected to print it elsewhere: xu2007 (in the arm labels), mottram2002, zhou2023 (Fig. 2, digitised) |
| **F3 base mixture** per arm: the amine, sulfur and carbonyl sources other than the lever's own level, as compound classes | the keyed `precursors` dict minus the lever level, names resolved through `data/keys/` and classed amine / sulfur / carbonyl / other | the sources behind the conflict cells and the widest-spread cells (about 10, listed in the check file) |

Two readers, blind to each other and to the table's readings, check those sources and write
`docs/reference/lever_table/field_checks/<source>.reader_{a,b}.yml`. A resolver writes `<source>.yml`. Agreement is
counted as exact on F1, within 0.1 pH unit on F2, and on compound classes (not names) on F3. **Kill rule:** F1
agreement below 80 % on the four sources stops the step, and the field definition is rewritten.

**Nothing is written into the resolved contrast files, and `build_lever_table.py` is not edited.** Both are inputs
of the frozen step-5a and 5c records ([lessons](../../../tasks/lessons.md), "Never edit the generator of a frozen
record"). A new generator, `scripts/generators/build_lever_field_readings.py`, derives the fields, reads the check
files and writes `results/validation/lever_field_readings.{json,md}`. **The keying template gains the three fields**,
so every source keyed from now on records them directly.

## 3. The frozen candidate rule (PREREG amendment A6)

Committed before the new generator is run. The rule:

1. **A6-R1, pH banded by control.** A pH cell is split into a `buffered` band and an `unbuffered` band
   (`adjusted_unbuffered` plus `unadjusted`). `unstated` enters neither band.
2. **A6-R2, final pH carried, never the step.** Each pH contrast carries its arms' final pH and their spread
   (Δfinal) as a covariate and a label. A slope on final pH is reported only where Δfinal ≥ 0.5 pH unit, as an
   alternative reading, and never in place of the keyed basis.
3. **A6-R3, matched base mixture.** A categorical pair whose arms differ in base-mixture classes beyond the lever's
   own level leaves the pair's level reading.

**How it is scored (RULES_V2_SPEC sec. "How a future rule change is tested"):**
- only on sources keyed after this commit, one rule at a time, on an identical contrast set and a fixed cell
  partition;
- against a placebo, a random banding with the same cell sizes;
- per lever class, on the coverage of the post-freeze contrasts by the new-laboratory 90 % interval, with interval
  width reported beside it.

No score is computed before at least 20 post-freeze contrasts fall in cells the rule touches.

## 4. Done means

- A6 committed;
- the template extended;
- the check files resolved;
- `lever_field_readings` generated: per answered, fragile and conflict cell, which reading or flag changes under
  each of A6-R1 - R3, and why, labelled as an audit;
- the main table and every frozen record byte-identical;
- a backlog item and a download list for **new thiol x pH sources that print initial and final pH**, the only route
  to a scored test and to a third laboratory on the thiol cells.

## 5. Size (wall clock)

- A6 and the template: 0.5 h;
- the generator: 2 h;
- reader checks and resolution: 1 h of work plus about 1 h of reader time;
- record and plan: 0.5 h.

**About 4 h.**
