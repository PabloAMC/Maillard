# Step 5a: the step-5 display rule on the engine as shipped (pre-registration)

*Written 2026-11-12, before any engine output on a TEST contrast of the lever table exists. Frozen by the commit that
adds it, together with the generator (`scripts/generators/build_step5a_display_check.py`) and the item set it built
without integrating anything (`results/validation/step5a_display_check_items.md`). A change after the first score
run is an amendment, dated and placed at the end, and never edits the text above it.*

## 1. What this is

STATUS_AND_PLAN sec. 0 step 5 displays the kinetic engine in a lever x odorant cell only where, on that cell's TEST
contrasts, its log-ratio error is below both "no change" and the DEVELOPMENT-only lever table's prediction, the
90 % interval of each difference excluding zero. Step 5a applies that rule, **as a diagnostic**, to every cell that
has TEST contrasts from two or more laboratories. The engine is the shipped one: no constant, structure, switch or
band is changed, nothing is fitted, and nothing ships on the result. This is one of the step-5 scoring looks
(at most 3 before 2027-01-23).

The session that wrote this file and runs the score opened only the DEVELOPMENT-only table
(`results/validation/lever_table_development_only.json`); it did not open the all-data table, the experiment map or
the disagreement audit's artifacts (PREREG.md sec. 2).

## 2. Inputs, frozen

- Rules of the lever table: `docs/reference/lever_table/PREREG.md` with amendments A1 - A5, as on `origin/main` at
  commit `8ea2b58558b0dc154cf83e5771731eebfd6b85f1` (file sha256
  `f2ad0acac2c2918a81f3609c64746aebcecd1ccbb60eb76da1d10620d5bb11ef`).
- The rival: `results/validation/lever_table_development_only.json` at that commit, sha256
  `8e4fceef7f0a3c1533f6f8004713264ad6a5579315d307d727834bd5c5dac5e0`. Only its `fit` estimates are used (a numeric
  cell's `estimate`; a categorical cell's per-level `estimate` and its `reference`), whatever the cell's reading
  ("too few laboratories" included); any `fragility` or sensitivity field is ignored. The estimates used are hashed
  by the generator (`RIVAL_FIT_DIGEST`,
  `e5a8dfe97195ec2d595d8c69e7aec5abef09938716e3c9018d86d4b5b1a13b3f`); the score run stops if that digest changes.
- Membership: `docs/reference/contrast_benchmark/split.yml` (sha256
  `1098b4b420330aa737497f6c855454c732dc394ce2007be6fd94dd9933cb36ac`) and
  `docs/reference/lever_table/new_labs.yml` (sha256
  `cf57620a5ca17e6703490ed3adc3ec7793b8628f51948660032c1be75a0f2842`), through `lab_partition()`.
- The TEST contrasts: the resolved files of `docs/reference/lever_table/contrasts/`, through `observations()` of
  `build_lever_table.py` (so the measured value of every contrast is the one the lever table uses).

**Sources removed from the TEST evidence** (their contrasts stay in the lever table):

| source | laboratory | TEST contrasts | why |
|---|---|---|---|
| `xu2007` | cau_gao_yanxiang | 21 | flagged `engine_seen` in `new_labs.yml`: read in step 0, a row of the directional claims panel (XU07-PH-01, -02), a scored table of the contrast benchmark |
| `kumazawa2003` | ogawa_chiba | 4 | fitted: its pH grid is a fit row of the thiol sub-model (`data/lit/fit_rows/thiol_submodel.yml`; `parameters_sulfur`, `kumazawa2003_ph_grid`) |

The other 19 TEST sources were checked against `docs/reference/FIT_HOLDOUT_DECLARATION.md`, the fit rows
(`data/lit/fit_rows/`), `docs/validation/directional_claims_panel.yml` and the benchmark panel; none is a fit row, a
hold-out row or a panel row. (`haleva1999` is named in the engine as the paper that would close a zero-rate edge; no
value of it is used. Stems such as `zhang2024` and `li2022` match other papers' dossiers, not these sources.)

## 3. The cells and the items

After the exclusions the lever table holds 750 TEST contrasts in 193 cells (odorant x lever class x band, the
table's own cell); **55 cells have TEST contrasts from two or more laboratories.** By lever class: precursor dose 20
(banded by the precursor varied, A5), sugar identity 16, sulfur source 9, amino acid identity 8, pH 1, time 1. The
list, with each cell's laboratories and sources, is `results/validation/step5a_display_check_items.md`.

An **item** is what is predicted: a numeric cell's contrast (one within-study slope per lever step, over the arms
the lever table's `numeric_observation` used: two arms, or the longest monotone run of a ladder), or a categorical
cell's pair (a level against the study's first-listed arm). 396 items. The measured value `m` of an item is the
lever table's own number (log10 ratio per declared lever step). An item with a censored arm keeps the lever table's
value (the limit entered as a value; the table's known departure from PREREG sec. 3), is flagged, enters the sign
score and is left out of the error medians, because its `m` is a bound.

**Which odorants the engine can simulate at all.** 25 odorants appear in the 55 cells. Eleven have an engine target
(`engine.TARGET_ALIASES`): MFT, bis(2-methyl-3-furyl) disulfide, 2-acetylthiazole, furfural, methional,
2-methylpyrazine, 2,5-dimethylpyrazine, unsubstituted pyrazine, hexanal, nonanal, 2-pentylfuran (25 cells). Fourteen
have none and their 30 cells are **not simulable**: benzaldehyde, phenylacetaldehyde, 2- and 3-methylbutanal,
2-methylpropanal, 2-furanmethanol, 2-acetylfuran, 1-octen-3-ol, heptanal, octanal, pentanal, 2-ethylpyrazine,
2-ethyl-3,6-dimethylpyrazine, and "Pyrazines" (a class total).

## 4. From a keyed contrast to an engine run

One rule set, applied by the generator to every arm of every item; an arm that fails a rule makes the item
**not simulable**, with that reason, and the item is listed, never dropped. Nothing is filled in by judgement.

1. **Target.** The odorant's engine target (sec. 3). None: not simulable.
2. **pH.** The keyed `ph`, used as the pot's pH at the start, whatever its basis (initial, label, unstated); a
   `held` pH is clamped; a pH printed only as final, or not printed: not simulable. No default pH.
3. **Temperature and time.** The keyed `temp_C` and `time_min`, as one isothermal hold (no heat-up or cool-down is
   keyed). Either missing: not simulable.
4. **Precursors.** Every keyed precursor other than water must be a front-door precursor name of the engine
   (`engine.PRECURSOR_ALIASES`, after dropping an L- / D- prefix and a parenthetical) or an ingredient on file
   (`data/species/ingredients.yml`, expanded by `src/ingredients.py` at its mass share). Its charge in mM comes from
   the printed amount, first form that matches: mmol in mL; g per 100 mL, g per 100 g or % (1 kg/L); a molarity; a
   plain mass over the pot's mass; a ratio to protein times the keyed protein content. The pot's mass behind a plain
   mass is, first that applies: a printed total; the keyed masses over a printed solids share; a printed water
   volume plus the masses; else the sum of the keyed masses. Molar masses are those of the named free compound
   (thiamine as its hydrochloride). A zero charge is left out. An amount not printed, or a name that is neither a
   precursor nor an ingredient: not simulable ("charge not derivable", "unmapped precursor").
5. **Buffer.** "x M phosphate" is that phosphate buffer; "none ..." is a declared unbuffered pot (the engine
   computes the pH trajectory); any other text leaves the engine's own default.
6. **Matrix: a food or a hydrolysate.** Two readings, both scored:
   - **Front door (decides display).** The pot must be a solution of pure compounds: the keyed matrix is "buffer"
     or "water", and no keyed precursor is a matrix material (a name that is no engine precursor and holds
     hydrolysate, protein, extract, powder, meat, trimmings, carcass, solids or "free amino acids"). Anything else
     is not simulable: "matrix material not on file". This is what the tool does today for a formulator: an
     ingredient with no sourced composition is refused by name.
   - **Added precursors (secondary; never decides display).** The matrix material is left out and the pure
     precursors the study added are charged in water, everything else as above. It asks a narrower question: does
     the engine's response to the added lever point the right way when the background is missing? A cell that
     "wins" here is not displayed; it marks where a matrix front door would be worth building.
     A sensitivity run of this reading holds the pH at the keyed value (clamped) instead of rule 5.
7. **The envelope.** `engine.declare_envelope` must answer every arm for the target (a lipid aldehyde with no lipid
   carrier, a trunk-only target in a pentose pot, a thiol with no sulfur source are refused there). An extrapolated
   answer is run. No matrix sites, vessel, atmosphere or water activity are set (none is keyed in a usable form).
8. **The lever must be expressed.** If all arms of an item map to identical engine inputs (ribose against xylose:
   both are the engine's one pentose), the engine has no term for the lever and the item is not simulable; an
   answer of two identical numbers is not a prediction (the rule of `engine.axis_refusal`).
9. **Categorical swaps** need no further rule: each arm's level is in its keyed precursors, so rules 4 and 8 decide.
10. **The run.** `engine.compare(arm, first arm of the item, [target])`, the shipped path of a ratio, with its own
    refusals (the falsified pH balance, an unidentified hexose entry, an axis the lane has no term for). A refused
    or undefined ratio on any arm makes the item not simulable ("the engine refuses the ratio"), counted apart
    from the mapping failures because it is only known after integration. A thiol ratio the engine marks "not
    resolved" still has a point value, which is used.

**Coverage, known before any integration** (`step5a_display_check_items.md`):

| reading | items that map | cells that keep two laboratories | why the others do not map |
|---|---|---|---|
| front door | **0 of 396** | **0 of 55** | no engine target 184; matrix material not on file 148; pH not printed 64 |
| added precursors | 25 of 396 | 1 of 55 (furfural x sugar dose) | no engine target 184; the engine refuses the arm 95 (no lipid carrier; trunk-only target in a pentose pot; no sulfur source); pH not printed 64; no engine precursor charged 18; charge not derivable 4; no term for the lever 4 (pentose against pentose); unmapped precursor 2 |

So the outcome of the display rule is already fixed by the mapping: **under the front-door reading no cell can be
displayed**, because every TEST contrast of the 55 cells is in a food or hydrolysate matrix, lacks a printed pH, or
measures an odorant the engine does not carry. The score run can add only the secondary reading: one cell judged on two
laboratories (furfural x sugar dose, 12 contrasts), and 13 more items in three cells with one laboratory each
(furfural x sugar identity 9, furfural x cysteine dose 2, bis(2-methyl-3-furyl) disulfide x sugar identity 2), shown
and not judged.

## 5. The three predictions and the errors

For an item with measured `m` (dex per lever step):

- **Engine** `p`: a categorical pair, log10 of the engine's ratio of the two arms; a numeric contrast, the weighted
  least-squares slope of the engine's log10 ratios on the lever's steps, over the same arms and with the same
  weights as `m`. An engine change below 0.001 dex over the whole contrast is numerical and read as zero.
- **No change**: 0.
- **Table** `t`: the DEVELOPMENT-only table's estimate for the item's cell (numeric), or the difference of its two
  levels' estimates (categorical; the reference level is 0). Where that table has no cell, or a level is absent or
  not connected, the table has no estimate and its prediction is 0: it then coincides with "no change", and the
  count of items with a table estimate is printed.

Errors are absolute: |p - m|, |m|, |t - m|. Per cell, per lever class and overall, on the simulated items with no
censored arm: the median of each error; the two differences, median engine error minus median rival error
(negative: the engine is closer); the number of contrasts and laboratories. **Sign accuracy**, on the items whose
measured change over the whole contrast exceeds the measurement noise (log10 1.4, the directional panel's fold):
the share where the prediction has the measured sign, for the engine and the table (a prediction of exactly zero is
a miss; "no change" has no sign).

**Interval.** A paired bootstrap over laboratories: the cell's (or class's) laboratories are resampled with
replacement, 2,000 times (seed 20261112), every item of a drawn laboratory entering with all three of its errors;
the 90 % interval is the 5th to 95th percentile of each difference. With two laboratories the resamples are A+A, A+B
and B+B, so the interval excludes zero only when each laboratory alone and both together agree: the rule is strict
there by construction, which is the reason step 5 asks for two laboratories.

**Readings.** The engine "beats" a rival where the interval of the difference lies below zero, "loses" where it
lies above, "ties" otherwise; with one laboratory the difference is shown and not judged.
**A cell is displayed** when, under the front-door reading, its simulated items come from two or more TEST
laboratories and the engine beats both rivals. The share of the cell's TEST contrasts that were simulated is
printed beside every reading.

## 6. Outputs

`results/validation/step5a_display_check_items.{json,md}` (this commit) and
`results/validation/step5a_display_check.{json,md}` (the score run), both from
`scripts/generators/build_step5a_display_check.py`. The score run asserts that the item set and the rival's fit
digest are the ones recorded here. The rival's file is not a hashed input of the artifact: the digest stands in for
it, so that additive fields in the table's file (PR 130) do not make this record stale.

## 7. What this cannot show

- The secondary reading is not the tool: a formulator's hydrolysate pot is refused at the front door, so nothing
  under it is a display.
- Two laboratories per cell: an interval that excludes zero is weak evidence, one that includes it is no evidence
  of equality.
- The measured values are the lever table's (peak areas mostly; a ratio within one study and method).

## Amendments

(none)
