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

**Amendment (2026-11-14).** For a designed study, the trigger is the laboratory's commitment (an agreed protocol
and lot), not the arrival of its data. Closing the input gaps takes engine work of the order of 100 h (STATUS_AND_PLAN
W3, A2 + H3); a laboratory's method work and runs take weeks; starting on commitment lets the engine be frozen, and
its freeze recorded, before any value of that study is read, which is what this section asks. The gaps are judged
against the study's own design (the P3 dry run, `results/validation/p3_engine_dry_run.md`), never against contrasts
already read. Literature-only work (scoping, keying blind, the perception layer) was never paused.

**Amendment (2026-11-17, owner).** The owner lifts the pause for one named change: **thiamine's H2S release outside
its measured range is refused, not set to zero.** Today `THIAMINE_H2S_SHARE_FLOOR_PH` = 6.25 sets the share to 0 below
pH 6.25 at every temperature. That zero is measured only near 100 C. At 120 C the two assays of one laboratory
disagree by more than 25x at pH 6 (Dwivedi & Arnold 1971, methylene blue, about 0.36 of the thiamine degraded;
Dwivedi & Arnold 1972, 35S, at most 0.014), and nothing is measured above 121 C. Outputs that depend on that share
there are refused, as the pH-balance refusal (S3) does for its rivals. No constant is fitted and no TEST contrast is
read, so the step-5 looks are untouched.

A premortem (2026-11-17) ruled out the two fits first scoped:

- a fitted share(pH, T): not identifiable (the 120 C profile is not monotone in pH, the share drifts with time, and
  the two assays disagree in both numerator and denominator); its hold-outs were seen at design, and three of the four
  E18 TEST contrasts cannot run through the front door. A declared band (about 40 - 50 h) waits for P3 or the
  partner's thiamine series, which measure this cell directly;
- 3-mercapto-2-pentanone's sugar route (Cerny & Briffod 2007): no shipped output (not a P3 analyte), its concentration
  pair is circular with the fit row, and its constant is confounded with the parked M1/M2 routes. It waits for M1/M2.

Everything else stays paused: water activity for the sulfur lane, methional and the Strecker arm, trimethylpyrazine
and the L1 lipid sink lack the evidence. The forecasts already sent outside the repository are frozen and are scored
as written.

**Record (2026-11-25): partner blind test, scoring and fit plan frozen before any data.** The scoring rules, the
scoring script with its synthetic tests, and the plan for what the partner's data may later fit (at most the
thiamine H2S share above and glutathione's H2S release, each with fixed priors and acceptance tests) were frozen
before any measured value was seen, after an independent premortem. The files describe an unpublished design and
stay private; the SHA-256 of their manifest is
`34df69a9a15a2380a5cddfa90c7a67af3fda844edd8f2ead1c0c04d085145307`. Data and any fitted value stay private until
the partner publishes or agrees in writing; this record will say whether a fit was accepted, without numbers.

**Amendment (2026-11-26, owner).** The owner lifts the pause for two named changes in the sulfur lane: **cysteine
reducing thiol disulfides back to the free thiol (thiol-disulfide exchange), and cysteine raising the yield of
2-methyl-3-furanthiol on thiamine's route.** A partner's development data (private) showed both effects in a
low-water dough; the published record already holds both (CF-089; Zhang 2024 Fig. 1-2; Grosch & Zeiler-Hilgart 1992,
spent). Both are built from the literature, not fitted to the partner's data, and go through a pre-registered,
engine-free feasibility test on mechanism signatures before any engine change
([kinetic_core_b65_cysteine_redox_prereg.md](../../../results/validation/kinetic_core_b65_cysteine_redox_prereg.md)).
An implementation needs a second owner approval after the feasibility verdict. No TEST contrast is read. Everything
else stays paused.

*Outcome (2026-11-26).* The screen found no admissible step. Exchange held at equilibrium is refuted. Zhang 2024's
cysteine boost needs either a structure not on the list or a pentose route outside its band. The dimer share may form
during sampling. No engine change was made, and the pause is back in force for these changes until one of the two
measurements now in `decisive_measurements.yml` exists.

## Bookkeeping

5a and 5c are two of the step-5 scoring looks at TEST contrasts that sec. 0 allows (at most three before
2027-01-23). Neither changed the engine, and the pause means no third look is planned.
