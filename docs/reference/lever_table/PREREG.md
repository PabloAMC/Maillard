# The lever table from data: pre-registration (STATUS_AND_PLAN sec. 0, step 2a)

*Written 2026-11-10, before any contrast value is keyed for this table; an independent premortem (ADOPT WITH EDITS) was folded in before the freeze. Frozen by the commit that adds it; a change
after the first value is keyed is an amendment, dated and placed at the end, with its reason, and never edits the
text above it.*

## 1. What the table is

The first layer of the product (STATUS_AND_PLAN sec. 0): for each **lever class x odorant** cell, the direction and
size of the lever's effect, with an interval, **estimated from the literature's own within-study contrasts and
nothing else**. No engine output enters it. A cell with too little data stays empty and says so; the experiment
map (step 4) then names what would fill it.

The question a cell answers: *"In a heated Maillard system, if I move this lever by one declared step and keep
everything else fixed, how does this odorant change?"*. Each cell states the scope it was estimated on (model
system, meat, plant protein; temperature range; laboratories).

## 2. Held-out membership, fixed before keying

- TEST membership of every atlas contrast is the frozen split v1:
  `docs/reference/contrast_benchmark/split.yml`, file sha256
  `1098b4b420330aa737497f6c855454c732dc394ce2007be6fd94dd9933cb36ac`, internal hash (`sha256` field)
  `022d0cab5a53ececb27a1fc6837ff275fbab177cf2161ae1f8a4ed0341f273a7`. The split is never redrawn.
- A contrast from a paper that is not in the atlas today takes the partition of its laboratory in `split.yml`
  (`lab_ids.yml` for the mapping); a laboratory absent from the split is TEST.
- The table is generated twice: on **all** keyed contrasts (the product), and on **DEVELOPMENT only** (the rival the
  engine must beat in step 5). A session that runs step 5 opens only the DEVELOPMENT-only table.

## 3. The unit and how it is keyed

**Unit.** One contrast = one odorant measured in two or more arms of one study (one table or figure, one
laboratory, one method) that differ in ONE lever, everything else as reported equal. Arms that differ in two levers
are keyed but enter only a cell defined by that pair ("mixed"), never a single-lever cell.

**Per arm:** the value as printed; its unit; whether it is an absolute concentration (calibrated or SIDA) or a peak
area or relative response; n and the SD if printed; the lever's value; the conditions (pH as labelled and whether it
is the initial, final or held pH; temperature; time; precursors and their concentrations; buffer; matrix;
headspace); the source location (table, figure, page). A figure is digitised only where the paper prints no table;
a digitised arm carries an extra 0.05 dex of error. Each contrast is flagged **focal** (the lever and the odorant are
the paper's stated question) or **incidental**.

**Censored arms.** "Not detected" with a printed LOD or LOQ is left-censored at that limit. "tr", "<x" or "not
detected" with no limit is left-censored at the lowest quantified value of that compound in the same table or
column. With no quantified value in the table it is keyed but enters no estimate, and is counted per cell.

**Dependence.**
- A multi-arm study contributes its arms jointly: the shared arm's variance enters as a covariance term.
- A ladder of three or more levels enters as ONE within-study slope over its monotone range (sec. 4), not as
  adjacent pairs.
- A time course contributes the single time the paper reports as its main result, else the longest time; the rule
  is applied before the values are read.
- Data reprinted in a second paper (a `duplicate` kind in `docs/reference/conflicts_register.yml`, or found while
  keying) is keyed once.
- Odorants measured in the same pots enter different cells; no cell pools across chemical families.

**Two readers, operationally.** Readers A and B are separate subagents with fresh contexts. Each gets only: the
source path, secs. 3 - 4 of this file, and the empty contrast template. Neither opens
`docs/reference/pathway_atlas/`, `results/`, the other reader's file, or any dossier's `quantity` or `contrast`
field (they anchor the reading). Each writes `contrasts/<source>.reader_a.yml` or `.reader_b.yml`. A third subagent
compares the two and re-reads the source only for the rows that differ, writing `contrasts/<source>.yml`. Only the
fields that define the cell, the lever step and the band are compared (lever values, odorant, unit type, arm
values; a difference above 0.05 dex); free-text conditions are not. The sources are the PDFs in the main checkout's
`data/articles/` (worktrees have none). Where no PDF exists, the source is the dossier's verbatim quote, flagged
`source: dossier_only`, keyed once and not double-read.

**Order and budget.** Thiols (MFT, FFT, H2S and their disulfides) x pH, x sulfur source, x sugar identity first.
At most 28 h of keying; at least 8 h are kept for the model and its run. At 12 h the keying rate is recorded and the
cell order is cut to fit. After the 40 h, a new estimate for the remaining cells goes into the backlog before going on.

## 4. The estimand

For each contrast, the log10 ratio of the odorant between arms, normalised to one declared lever step. Every
contrast records the lever's midpoint.

| lever class | step | bands (pooled only within a band) |
|---|---|---|
| pH | +1 pH unit; initial-pH and held-pH studies are separate sub-cells | pH < 4.5, 4.5 - 6.0, 6.0 - 7.5, > 7.5 |
| temperature | +10 C | heating time < 30 min, 30 - 120 min, > 120 min |
| time | a doubling of heating time | temperature < 100 C, 100 - 130 C, > 130 C |
| precursor dose | a doubling of the varied precursor's concentration | -- |
| water | per 0.1 a_w, as reported | -- |
| buffer concentration | a doubling | -- |
| categorical: sugar identity, sulfur source, buffer type, lipid (added vs none; oxidation state), oxygen (air vs inert), matrix (protein or food vs buffer) | a level against another level | -- |

**Ladders.** A study with three or more levels contributes one slope: the least-squares fit of log10 y on the lever
over its monotone range, with its standard error. A ladder whose response changes sign inside its range is flagged
non-monotone; its monotone pieces enter their own bands.

**Categorical levers are a network.** Every within-study pair is keyed against the study's own first-listed arm.
The model has a fixed effect per level and a study intercept, and the cell reports each level against glucose
(sugar), L-cysteine (sulfur source), phosphate (buffer type), none (lipid), inert (oxygen) or buffer (matrix) as a
contrast of the fitted effects. A study with none of those arms still enters through the levels it shares with other
studies; a level connected to the reference by no chain of studies is reported as "not connected".

**Peak areas.** A within-study ratio of one compound by one method cancels the response factor, so peak-area
contrasts are kept; each cell reports their share and an estimate without them.

## 5. The model

Per cell (and band), a hierarchical model on the normalised log ratios y:

    y_c = mu + a_lab(c) + b_study(c) + e_c,
    a ~ N(0, tau_lab^2), b ~ N(0, tau_study^2), e_c ~ N(0, s_c^2 + sigma_resid^2)

For a categorical cell, mu is replaced by the per-level effects of sec. 4.

- s_c: the contrast's measurement variance (with its covariance for shared arms), from printed SDs where they exist;
  otherwise the declared per-replicate CV of step 0 (10 % for isotope dilution, 20 % otherwise; n = 2 where
  unstated).
- tau_study is fitted only where at least two laboratories have at least two studies each; otherwise it is merged
  into tau_lab.
- A Bayesian fit with half-normal priors on the scales (0.5 dex). **Prior sensitivity:** every cell is refitted at
  0.25 and 1.0 dex; a cell whose 90 % new-laboratory prediction interval changes width by more than 30 % is labelled
  "prior-dominated". REML is a check only; where REML gives tau = 0 the Bayesian fit is reported.
- Reported per cell: mu with its 90 % interval; the **90 % prediction interval for a new laboratory** (the number a
  formulator should read); P(a new laboratory sees the same sign); the counts of contrasts, studies and laboratories;
  the matrix and temperature scope; the share of peak-area contrasts; mu on incidental contrasts only (the selection
  check, sec. 6); the prior-dominated label.

**Reading a cell.**
- **"answered by the literature"**: at least three laboratories, the 90 % interval of mu excludes zero, AND
  P(a new laboratory sees the same sign) >= 0.8 from the posterior predictive. Three laboratories, because with two
  the between-laboratory variance has one degree of freedom and its interval is uninformative; 0.8, because a
  formulator acts on a direction they will see in their own pot four times in five. The direction is stated, and the
  size as the new-laboratory interval.
- **"direction unclear"**: at least three laboratories, either condition fails.
- **"too few laboratories"**: fewer than three; the contrasts are listed, nothing is pooled.

**Validation.** Leave-one-laboratory-out on every cell with at least four laboratories: each laboratory's contrasts
predicted from the others. Reported: the coverage of the 90 % new-laboratory interval, and the median absolute error
against "no change". A cell whose coverage falls below 0.7 is flagged: its interval is too narrow to show.

## 6. What this table cannot do, stated in every output

- **Selection.** Published contrasts are a selected sample: levers that did nothing are under-reported. A pooled
  effect is an upper estimate of the typical effect.
- **Scope.** Most contrasts are model systems in buffer. A cell estimated only on model systems says so; it does not
  claim to transfer to a pea-protein patty.
- **No mechanism.** A cell reports what happened, not why; it does not predict a combination of levers that no study
  varied together.

## 7. Files and rules

- `docs/reference/lever_table/contrasts/*.yml`: the keyed contrasts, one file per source, the two readers' files
  beside the resolved one.
- `docs/reference/lever_table/new_labs.yml`: laboratories absent from the split (TEST by sec. 2). `lab_ids.yml` and
  `split.yml` are not edited and `build_contrast_split.py` is not run.
- `scripts/generators/build_lever_table.py`: the model. It writes `results/validation/lever_table.{json,md}` (all
  data) and `results/validation/lever_table_development_only.{json,md}`. Its input paths go into
  `src/data_paths.py`; it resolves odorant names through `src/compound_keys.py` and refuses a name that does not
  resolve; it checks every contrast file against a schema before fitting. Nothing is written under `data/`.
- This file: the frozen rules.

## Amendments

*(none)*
