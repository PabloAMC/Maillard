# Kinetic core, Build Wave B3 -- the acrylamide hold-out scorecard

Frozen parameters from `results/validation/kinetic_core_b3_fit_report.json`. **Nothing was fitted in this script and there is no optimiser in it.**

- gating rows: **0 / 4 pass** (level band 3.0x, ratio band 2.0x)
- diagnostic rows: 3
- unscoreable rows: 3 -- and each one names WHY, because 'unscoreable' and 'pass' are not the same thing

## What the scorecard actually says

**3 of the 4 gating rows were pre-registered as failures BEFORE the fit was run**, with the mechanism named: this module has no glutamine promotion route and was deliberately not given one, because the promotion's temperature shape is the hold-out itself. Those rows failing is the expected outcome and is not new information; what they buy is a localisation -- the missing thing is a MECHANISM, not a mis-set constant.

**The genuinely open row is `knol2010_molar_yield_on_asparagine`.** Result:

- `knol2010_molar_yield_on_asparagine`: predicted 0.007807 against 0.03, a factor of 3.84 LOW, against a 3.0x band -- **FAIL**.

Read the three UNSCOREABLE rows as what they are. Two of them (Knol 2010's step table, Knol 2009's real-food band) are CORPUS gaps: the numbers are not transcribed anywhere in the extraction dossiers, so the declaration's strongest acrylamide hold-out is only partly available. The third (sucrose) is a scope limit this module chose deliberately in order not to read the hold-out. None of the three is a pass.

## Row by row

| row | gating | observed | predicted | fold | result |
|---|---|---:|---:|---:|---|
| `dv1_fructose_over_glucose_acrylamide` | DIAGNOSTIC | 2.073 | 0.2686 | 0.13x | FAIL |
| `dv1_fructose_peak_acrylamide_ppb` | DIAGNOSTIC -- REPORTED, NOT SCORED | -- | 1.273e+06 | -- | n/a |
| `dv1_sucrose_system` | UNSCOREABLE | -- | -- | -- | n/a |
| `dv2_glutamine_promotion_120C` | GATING | 2.67 | 0.003317 | 0.00x | FAIL |
| `dv2_glutamine_promotion_200C` | GATING | 1.2 | 0.02771 | 0.02x | FAIL |
| `dv2_glutamine_kF_ratio_160C` | GATING | 2.255 | 0.9999 | 0.44x | FAIL |
| `knol2010_molar_yield_on_asparagine` | GATING | 0.03 | 0.007807 | 0.26x | FAIL |
| `knol2010_formation_barrier_ceiling` | DIAGNOSTIC (one-sided) | 93 | 190.6 | -- | FAIL |
| `knol2010_table2_seven_steps` | UNSCOREABLE | -- | -- | -- | n/a |
| `knol2009_real_food_band` | UNSCOREABLE | -- | -- | -- | n/a |

## Intervals on the low-confidence rows

- `dv1_fructose_over_glucose_acrylamide`: median **0.2705**, 95 % band **0.2442 - 0.2979** (a factor of 1.2), from 200 draws.
- `dv1_fructose_peak_acrylamide_ppb`: median **1.23e+06**, 95 % band **7.596e+05 - 1.735e+06** (a factor of 2.3), from 200 draws.

Read the two bands together, and note that they are NOT the same kind of statement. The RATIO band is narrow because most of the parameter uncertainty cancels between the two sugars -- the two systems share every constant except the isomerisation that carries fructose into the glucose lane. The LEVEL band is wider. NEITHER contains the structural uncertainty, which here is the dominant term: the model has no fructose-specific chemistry at all, and the isomerisation constants that do all the work were measured 40-80 C below where they are being evaluated. A narrow band around a structurally wrong route is a precise wrong answer, and that is the honest reading of this row.

## What each row means

### `dv1_fructose_over_glucose_acrylamide` -- DIAGNOSTIC

Source: De Vleeschouwer 2009 I Table 3, fructose k_Fref 7.40 +/- 9.48 e-3 against the FIT glucose 3.57 +/- 1.38 e-3 min^-1

LOW CONFIDENCE, and that is the required answer rather than a hedge. The source's own 95 % HPD on the fructose constant is +/- 9.48 against an estimate of 7.40, i.e. IT SPANS ZERO: the target is not a determination and a model that matched it precisely would be matching noise. The declaration says so explicitly ('fructose's HPDs span zero, which the model should reproduce as LOW CONFIDENCE, not as a fitted value'). The model carries NO fructose-specific parameter -- fructose reaches acrylamide only by isomerising to glucose through Martins' measured constants, evaluated 40-80 C above the window they were measured in. The interval below, not the point value, is the prediction.

### `dv1_fructose_peak_acrylamide_ppb` -- DIAGNOSTIC -- REPORTED, NOT SCORED

Source: no absolute acrylamide concentration is printed for the De Vleeschouwer systems; only rate constants

There is nothing to score this against: the source publishes constants, not concentrations. It is reported so that the interval width on an absolute level is visible next to the interval on the ratio.

### `dv1_sucrose_system` -- UNSCOREABLE

Source: De Vleeschouwer 2009 I Table 3, sucrose column

The module has NO SUCROSE SPECIES. Sucrose enters the acrylamide lane only by hydrolysing to glucose and fructose, and the only measurement of that hydrolysis rate in the corpus is in this same held-out column. Adding a sucrose species would therefore have meant either inventing its hydrolysis rate or reading the hold-out. Reported UNSCOREABLE rather than failed, because the two are different findings: this one says the module's scope is narrower than the declaration's, not that its chemistry is wrong.

### `dv2_glutamine_promotion_120C` -- GATING

Source: De Vleeschouwer 2009 II, glutamine promotion of acrylamide at a_w 0.92, 120 C

PRE-REGISTERED FAILURE, direction and cause both known before the fit. Competition in this module can only SUPPRESS: a competitor either eats the shared glucose or scavenges the acrylamide, so the predicted ratio cannot exceed 1. Glutamine is measured to PROMOTE. The model was deliberately not given a promotion mechanism, because the promotion's temperature SHAPE (growing with T in liquid, shrinking with T at a_w 0.92) is the B5.5 sign-crossing whose a_w-0.92 half is THIS ROW. A term fitted to the Claeys half would have been a term built toward this hold-out. The failure localises a missing MECHANISM, which is what the declaration says this row is for.

### `dv2_glutamine_promotion_200C` -- GATING

Source: De Vleeschouwer 2009 II, glutamine promotion of acrylamide at a_w 0.92, 200 C

PRE-REGISTERED FAILURE, direction and cause both known before the fit. Competition in this module can only SUPPRESS: a competitor either eats the shared glucose or scavenges the acrylamide, so the predicted ratio cannot exceed 1. Glutamine is measured to PROMOTE. The model was deliberately not given a promotion mechanism, because the promotion's temperature SHAPE (growing with T in liquid, shrinking with T at a_w 0.92) is the B5.5 sign-crossing whose a_w-0.92 half is THIS ROW. A term fitted to the Claeys half would have been a term built toward this hold-out. The failure localises a missing MECHANISM, which is what the declaration says this row is for.

### `dv2_glutamine_kF_ratio_160C` -- GATING

Source: De Vleeschouwer 2009 II Table 3, glutamine k_Fref 8.05e-3 against the control's 3.57e-3 min^-1 (a 2.25x promotion)

PRE-REGISTERED FAILURE, the same one, on the formation constant rather than on the yield. The model's ratio can approach 1 and cannot exceed it, because suppression is the only thing a competitor in this module can do.

### `knol2010_molar_yield_on_asparagine` -- GATING

Source: Knol 2010: acrylamide yield ~3 % of the initial asparagine in an aqueous 0.1 M asparagine/glucose system

THE ROW TO READ. It is a THIRD lab, and -- more importantly -- it is at a TEN TIMES higher precursor loading than the Claeys system the fitted partition was identified on. Because initiation here is genuinely bimolecular, the predicted molar yield RISES with loading; a model with a fixed yield fraction (which is what the shipped lane had) cannot move on this axis at all. Z1 sec. E diagnosed exactly that saturation as the central defect, and this row is the first out-of-sample test of the fix.

### `knol2010_formation_barrier_ceiling` -- DIAGNOSTIC (one-sided)

Source: Knol 2010's largest published activation energy is 93 +/- 12 kJ/mol

A ONE-SIDED CEILING, not a level. What is known about Knol 2010 from the corpus is that no barrier in it exceeds 93 +/- 12; the seven step constants themselves are not transcribed. The model's apparent formation barrier is compared against that ceiling. Note the tension this row inherits: Claeys measures 168.25 +/- 3.80 for the same transformation and Claeys is a FIT row, so the model is being pulled in two directions by two datasets that disagree by 75 kJ/mol.

### `knol2010_table2_seven_steps` -- UNSCOREABLE

Source: Knol 2010 Table 2, 7 steps x 5 temperatures + Ea + HPD

NOT TRANSCRIBED ANYWHERE IN THE CORPUS. The extraction dossiers carry only Knol 2010's three organic-acid and isomerisation barriers, which orchestrator decision 2 assigns to Module 4 as FIT. This is a CORPUS gap, not a model limitation, and the declaration's strongest acrylamide hold-out is therefore only partly available: what can be scored is the end-to-end yield row above.

### `knol2009_real_food_band` -- UNSCOREABLE

Source: Knol 2009, 9.3e3 - 2.6e4 ug/kg dm in real food

No asparagine or glucose content for Knol's matrix is transcribed in the corpus, so a prediction into this band would begin by inventing the precursor loading it is most sensitive to. The source's own author refused to transfer his model to real food -- 'the parameters are only applicable for specific experimental conditions such as time-temperature profile of frying, potato genotype, slice thickness and initial concentration of precursors' (inventory sec. B8.5). Reporting UNSCOREABLE is the honest outcome of a test the evidence cannot support, and it is NOT a pass.
