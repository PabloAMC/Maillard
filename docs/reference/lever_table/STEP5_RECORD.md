# Steps 5a - 5c: where the engine can be judged, and why it pauses (STATUS_AND_PLAN sec. 0, step 6)

*2026-11-13. One record of the three step-5 diagnostics and what they imply. Nothing here moves a constant, a
reading of the lever table, or a shipped output.*

## The question

Layer 3 of the product (sec. 0) shows the kinetic engine in a lever x odorant cell only where held-out contrasts
show that it beats the lever table there. Steps 5a - 5c asked, on the data on disk, where that can be tested at all.

## What was found

**5a, the display check** ([step5a_display_check.md](../../../results/validation/step5a_display_check.md); rules
frozen first in [STEP5A_PREREG.md](STEP5A_PREREG.md)). 55 cells have TEST contrasts from two laboratories (396
contrasts, 15 sources, 11 laboratories). **The engine runs through the front door on none of them.** 184 measure an
odorant it has no target for (benzaldehyde, phenylacetaldehyde, the methylbutanals, C7 - C8 aldehydes,
ethylpyrazines, ...), 148 are in a hydrolysate or food whose composition is not on file, 64 print no pH. No cell is
displayed, and no cell can be judged: a coverage finding, not an accuracy one. In the secondary reading (the added
pure precursors only), 10 contrasts survive the engine's own refusals, all furfural x sugar dose, where the engine
ties "no change" and the DEVELOPMENT-only table (median error 0.17 against 0.15 dex per doubling).

**5b, the disagreement audit** ([lever_disagreement_audit.md](../../../results/validation/lever_disagreement_audit.md),
[DISAGREEMENT_READING.md](DISAGREEMENT_READING.md), [RULES_V2_SPEC.md](RULES_V2_SPEC.md)). Of 99 cells with two or
more laboratories, 9 are conflicts between laboratories. None needs measurement noise or a laboratory effect: each
comes from what the lever meant in each study (an unbuffered initial pH, a "none" arm that is a different pot,
opposite halves of a rise-then-fall dose curve) or from a rule of the table (a detection limit used as a value). Two
of the four answered cells are fragile to those rules (both pH). The binding gap is a third laboratory, not
disagreement.

**5c, the engine on the conflict cells** ([step5c_discordant_cells.md](../../../results/validation/step5c_discordant_cells.md);
rules frozen in the generator's docstring and committed with its item list before any engine output). The 9
conflict cells hold 32 contrasts or pairs, both partitions. **Through the front door the engine runs none**: 18 are
in a hydrolysate or food with no composition on file, 6 measure an odorant with no target, 2 charge a precursor the
engine does not have, and the 6 that map (MFT and its disulfide x pH: xu2007 against zhou2023) are refused by the
engine itself on the falsified pH balance (S3). In the secondary reading with the pH held at the keyed value, the two
thiol x pH cells run: for MFT the engine's between-laboratory difference has the measured sign (-0.88 against -0.98
dex per pH unit), which the frozen rule scores "explained by conditions"; for the disulfide it is not (-1.79 against
-0.55). **Neither is evidence that conditions explain the conflict.** The engine predicts a fall with pH in both
laboratories, three to ten times steeper than in xu2007's pots, and so gets zhou2023's sign wrong; holding the pH removes the
very cause 5b found (zhou2023's unbuffered pots ended near pH 3.3); and xu2007 is a source the engine was already
scored on. This is the step-1 finding again: the engine's pH response is structurally too steep.

## What it implies

1. **What stands between the engine and any held-out judgement is inputs and coverage, not a rate law.** In order of
   the contrasts each would open on disk: sourced compositions of protein hydrolysates and foods (148 of 5a's 396,
   18 of 5c's 32); engine targets for the Strecker and lipid aldehydes and the ethylpyrazines the literature reports
   (184 of 396); a printed pH (64); lanes that compose (methional and pyrazines with a pentose; a thiol with no
   sulfur source); and, for thiols, an answer to the falsified pH balance, which refuses every thiol x pH ratio.
2. **None of these can be judged on data on disk.** Each would add contrasts the engine can be run on, but the
   contrasts it would be run on are already read (the all-data table is open, 5a and 5c have looked at their
   residuals), and any structure fixed against them would be tuned to them. With no new laboratory in sight, that
   is the tuning loop the 2026-08 audit found.
3. **So the engine pauses.** No engine constant or structure is changed until outside data exist. The paused cases
   (M0a, M4, M5) stay paused: 5a named no cell where any of them could be judged.
4. **What carries the product meanwhile** is layers 1 and 2: the lever table with its sensitivity readings, and the
   experiment map. The "what to record" table of [EXPERIMENTS_BRIEF.md](../../guides/EXPERIMENTS_BRIEF.md) gained
   three rows from 5a - 5b (the buffer and final pH, the reference arm's base mixture, every dose including zero), and
   the P3 design ([FACTORIAL_LEVER_STUDY.md](../../guides/FACTORIAL_LEVER_STUDY.md), "Changes after steps 5a - 5b")
   now asks for the isolate's composition, without which the engine could not be judged on P3 either.

## When the pause ends

Any one of: a laboratory runs P3 (or any designed study) and its data reach the repository; the blind-test partner
shares data (DEVELOPMENT, so it informs but does not judge); papers from laboratories absent from split v1 are
registered before reading and keyed blind, giving a cell two TEST laboratories the engine can run. On the first of
these, the input gaps of point 1 are worth closing BEFORE the data are read, so that the engine is frozen when it
meets them.

## Bookkeeping

5a and 5c are two of the step-5 scoring looks at TEST contrasts that sec. 0 allows (at most three before
2027-01-23). Neither changed the engine, and the pause means no third look is planned.
