# Sub-model 1: pH and acid production (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2, sub-model 1. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`.
Sections marked TO DO are filled in at the step the procedure names; nothing here is changed after the fit except
under "Amendments".*

## 1. Frame (step 1)

**Question.** Does the model's pH machinery reproduce the measured pH trajectories and organic-acid production of
heated sugar + amine systems, across laboratories, buffers and temperatures, with constants the data identify?

**Why first.** Every step the M4 table flags as unclear (`results/validation/sulfur_reaction_evidence.md`) is
pH-dependent: the pentose + cysteine branches, 1-deoxypentosone's branches, norfuraneol and thiol loss, H2S loss.
A wrong pH trajectory biases every later sub-model's constants.

**What the model has now.** A charge balance solved from the state (`ph_state`; the inline coupling or B60's converged
fixed point), titratable species with temperature-dependent pK (B51), and acid made by sugar fragmentation: in the
trunk, measured steps (Martins: 3-DG -> formic acid, 1-DG -> acetic acid, fructose -> acids); in the sulfur lane, a
fitted acid yield per osone-sink event (`ph_acid_yield_per_sink_event`, B60 free) and the ARP's secondary-ammonium pKa
(`ph_arp_secondary_ammonium_pKa`, B60 free).

**Decision it changes.** Whether the shipped pH machinery goes forward unchanged; whether the acid yield needs a
measured basis per pathway rather than one fitted yield; whether the pH a sub-model uses is predicted or must be taken
from the data as a latent condition (the misfit protocol's "conditions" step).

**Done-test** (the operating procedure's mechanism checks):
1. own species across laboratories: pH(t) and the measured acids (formic, acetic, titratable acid, NaOH uptake) within
   the measurement model's errors in at least three of the four laboratories, by leave-one-laboratory-out;
2. no structural residual flag (trend with T, buffer strength, sugar, amine) left unexplained;
3. every fitted constant identified, none on a bound; the trunk's measured acid constants (Martins) not contradicted
   beyond their intervals;
4. carbon closes where the source gives totals (Davidek's acids + DFG loss; Brands' mass balance);
5. one frozen out-of-condition prediction: the pH trajectory of a laboratory and system left out entirely (candidate:
   Cao 2017 cysteine + xylose at 90 C, or Nguyen 2016 caseinate at 130 C), written down before it is scored;
6. rival arm: one fitted acid yield per sink event (current) vs pathway-specific yields from the measured trunk steps.

**Budget.** About 6-10 h of work and 2-4 h of compute (strategy, "Next steps" in hours).

**Expected result (written before the probe).** The charge balance is sound (B51 and B60 checked it), so pH should be
right where the acid supply is right: in pH-stat and buffered systems. Unbuffered drift is the test of the acid supply
and is where the current single yield is most likely to fail, especially across sugar types (pentose vs hexose) and
amines (glycine vs cysteine). Probability the shipped machinery passes the done-test unchanged: about 30 %.

## 2. Data (step 2: look up before building)

Candidates from the measurement ledgers (2026-10-29). Each enters a keyed fit table (`data/lit/fit_rows/`) only after
two-reader verification against the paper (`data/lit/source_verification.yml`).

| ledger id | system | what | quantification | lab |
| --- | --- | --- | --- | --- |
| davidek2002_tab5_acids_naoh_uptake | DFG 0.1 M, 90 C, pH-stat 5-8 | NaOH uptake, formic / acetic acid | calibrated | Nestle |
| martins2003b_fig6_organic_acids | DFG 10 mM, 100 / 120 C | formic, acetic acid | calibrated | Wageningen |
| martins2003b_fig7_ph_drift | DFG 10 mM, 100 / 120 C | pH(t) | calibrated | Wageningen |
| martins2005b_fig4_ph_drift_ladder | glucose + glycine 0.2 M, 100 C | pH(t) at five starting pH | calibrated | Wageningen |
| martins2005b_tab1_phstat_vs_drift_constants | as above | pH-stat vs drift rate constants | calibrated | Wageningen |
| brands2001_fig3_ph_titratable_acid | glucose / fructose 150 mM + caseinate, 120 C | pH(t), titratable acid | calibrated | Wageningen |
| nguyen2016_figs1_ph_mass_balance | glucose / lactose + caseinate, 120 / 130 C | pH(t) (read from the SI's vector figures) | calibrated | Wageningen |
| cao2017_fig4c_ph | cysteine + xylose 0.40 M, 90 C | pH(t) | calibrated pH meter | Xie group |
| zhao2017b_fig2c_ph_diluted | cysteine + xylose 0.40 M, 90 C | pH(t) | calibrated pH meter | Xie group |
| deng2021_fig6_phosphate_arm_ph_trajectory | MG-ARP 0.2 M, 100 C, +/- phosphate | pH(t) | calibrated pH meter | Deng |
| zhang2025_fig3_diarp_120c_dicarbonyls_ph | di-ARP of lysine, 120 C | pH(t), dicarbonyls | calibrated | Jiangnan |
| final-pH endpoints | Kang 2026, Zhou 2023, Chen 2000 / 2002, Xu 2008, Zhai 2023 | final pH only | pH meter | several |

Lookups still TO DO (step 2): the conflicts register for every source above; earlier waves that fitted pH (B2.4's
pH-endpoint weighting, B51, B60) and what they found; Knol 2010 Fig. 2 and Hou 2017 (figure-only pH series:
digitise only if the done-test needs them).

## 3. Cheap probe (step 3) -- done 2026-10-29, nothing refitted

`scripts/generators/build_ph_acid_submodel_probe.py` runs every pot of the READY table (41 pots, 396 rows) once on
the sulfur lane at the shipped constants, on the pot's measured times, and predicts each row by its observation
operator (`src/kinetic_core/observations.py`, tested against hand calculations in `tests/unit/test_observations.py`):
cooled-sample pH, pH after Zhao's 6x dilution, the pH-stat's NaOH uptake at 90 C, titratable acid to pH 8.3 relative to
the unheated sample, formic + acetic acid. 7 pots are skipped (lactose, the di-ARP, heptenal, and Deng's methionine ARP,
which the sulfur lane does not carry); 325 rows are scored; t = 0 pH rows are the charge balance's input and are not.
Three arms, ~40 s each on 8 cores:

| data set (observable) | n | shipped | + trunk pH terms in the lane | + and `k_aa_frag` off |
| --- | --- | --- | --- | --- |
| Davidek 2002, pH-stat 5-8 (FA+AA) | 53 | 0.61 | 0.36 | 0.31 |
| Davidek 2002 (NaOH uptake) | 49 | 0.50 | 0.52 | **0.26** |
| Martins 2003b DFG (acetic) | 32 | 0.80 | 0.51 | 0.49 |
| Martins 2003b DFG (formic) | 23 | 0.57 | 0.51 | 0.53 |
| Martins 2003b DFG (pH) | 16 | 0.10 | 0.10 | 0.06 |
| Martins 2005b Glc + Gly pH ladder 4.8-7.5 (pH) | 52 | 0.75 | 0.44 | **0.38** |
| Brands 2001 sugar + caseinate (pH) | 14 | 0.19 | 0.19 | 0.18 |
| Brands 2001 (titratable acid) | 13 | 0.66 | 0.66 | 0.64 |
| Brands 2001 (FA+AA) | 13 | 0.43 | 0.43 | 0.43 |
| Nguyen 2016 glucose + caseinate (pH) | 24 | 0.06 | 0.06 | 0.06 |
| Cao 2017 Cys + xylose 0.4 M (pH) | 12 | 0.93 | 0.93 | 0.93 |
| Zhao 2017b, same pot, 6x diluted (pH) | 18 | 0.77 | 0.77 | 0.77 |
| Zhou 2023 ARP +/- Cys, final pH | 6 | 0.98 | 0.98 | 0.98 |

RMS residual: pH units for pH, dex (log10 predicted / measured) for amounts; a negative mean is a model below the data. Artifacts:
`kinetic_core_ph_acid_submodel_probe{,_trunk_ph_dynamic,_trunk_ph_dynamic_without_k_aa_frag}.{json,md}`.

**Findings.**

1. **The sulfur lane's acid supply has no pH dependence.** Shipped, Davidek's DFG makes 20.0 mmol/L of formic +
   acetic acid at every pH from 5 to 8; the paper measures 3.3 to 79 (water) and 11 to 79 (phosphate). The trunk lane
   carries Martins 2003 Table 3's measured pH terms on the three Amadori-decay steps and the 3-DG -> formic exit
   (`trunk_conditions`, B12 / B40), but the sulfur lane -- the one lane whose pH moves -- never applied them.
   `sulfur.TRUNK_PH_TERMS` (new, OFF) applies them at the lane's computed pH on the label scale; it halves the Davidek
   and Martins 2005b errors with no constant refitted. This is a code-path omission, not a missing measurement.
2. **The fitted acetic-acid sink `k_aa_frag` is contradicted by two data sets and required by one.** B1 fitted it
   (t1/2 = 59 min at 100 C, Ea at its 20 kJ/mol lower bound) on Martins 2005's glucose + glycine acids. In DFG systems,
   where the precursor runs out, it makes acetic acid peak and fall (Martins 2003b 100 C: predicted 4.6 -> 0.84 mmol/L
   from 45 to 180 min against a measured plateau at 4.8; Davidek pH 8: 60 -> 32 against 44 -> 79). Acetic acid is
   thermally stable at these temperatures; no chemistry the corpus documents removes it. Switched off, Davidek's
   NaOH-uptake error halves (0.52 -> 0.26 dex). But on its own fit data (Martins 2005, Glc + Gly 200 mM, 80-120 C)
   the sink is doing work: without it acetic acid at 100 C / 120 min is 16.8 mmol/L (with the dynamic pH terms;
   19.0 without them) against 11.6 measured. The dynamic pH feedback explains about a third of what the sink absorbs.
   So the sink is a surrogate for something else in the glucose system -- candidates, each with a direction: the
   1-deoxyglucosone supply runs too long (1-DG captured by glycine / melanoidin, which only the glucose system has in
   excess: predicts the excess grows with glycine); the acetic-acid yield of the 1-DG exit is lower than 1 (predicts a
   constant ratio at all times, not a late excess); a reading error at late times in Fig. 5.10 (predicts nothing
   systematic across temperatures). It goes into the fit as a rival arm (section 4), not deleted on this probe.
3. **Pentose + cysteine at 0.4 M is over-acidified by 0.8-0.9 pH units, in two papers of one laboratory, unchanged by
   either arm.** Cao 2017 and Zhao 2017b (same pot; Zhao reads a 6x dilution) measure 4.7 -> 4.15 over 5 h at 90 C;
   the lane goes to 3.6 in 30 min and 3.4 by 5 h. The two papers agree with each other through the dilution operator
   (diluted vs undiluted within 0.2 units in both data and model), so it is the acid supply, not the observation. The
   lumped pentose acid pool (`ACID`, yield per sink event, B2.2) was calibrated on Zhou's 20 mM pots; at 400 mM it
   produces too much acid too fast. A yield that is not proportional to concentration, or acid made by a step with a
   different order, are the candidates.
4. **Zhou 2023's six final-pH endpoints -- B2.2's own calibration anchors -- are missed by up to 2.1 units**: ARP alone
   at pH 8 ends at 7.10 predicted vs 4.95 measured, at pH 6 at 3.19 vs 3.98. Two constants cannot carry the pH and
   cysteine dependence of the pentose acid supply.
5. **Casein systems: pH right, titratable acid low by 3x.** Nguyen's pH is within 0.06; Brands' titratable acid is
   under-predicted 2-6x while its pH is within 0.2. The lane carries neither the caseinate's titratable groups nor its
   lysine (LYSP is a trunk species), and the amine-free sugar loss (`k_glc_frag`) makes no acid. The pH agreement is
   therefore partly cancellation (too little acid, too little buffer); not fitted here (section 4).
6. **What holds:** the charge balance and the observation operators. Nguyen's and Martins 2003b's drifts are within
   0.06-0.10 pH units; the pH-stat's base uptake matches the measured acids where the acids match.

**Against the expected result** (written before the probe): the shipped machinery fails the done-test, as expected
(P = 0.70 that it would). The expectation that "pH should be right in pH-stat and buffered systems" was WRONG: the
pH-stat data fail worst, because the failure is in the acid supply's pH dependence, not in the charge balance.

## 4a. Design (step 4) -- from the probe (2026-10-29)

- **Arm A (base):** `TRUNK_PH_TERMS` on, the trunk's three Amadori pH slopes and the 3-DG -> formic slope fitted with
  Martins 2003 Table 3's values and bands as priors (`trunk_conditions.PH_EXPONENT_BAND`); `k_aa_frag` kept.
- **Arm B (rival):** as A with `k_aa_frag` off. The two are compared by the evidence audit and leave-one-laboratory-
  out, and the glucose-system excess (finding 2) is tested by its candidates' directions, not absorbed.
- **Data to add before the fit:** Martins 2005's glucose + glycine formic and acetic acid at 80 / 100 / 120 C
  (`data/lit/timeseries/martins2005_glucose_glycine_80_100_120C_pH68.yml`, B1's fit data) as keyed rows, two-reader
  verified, so both arms are judged on the data the sink came from.
- **Pentose acid supply (finding 3-4):** the acid yield per sink event and the ARP pKa (B2.2) refitted on Cao, Zhao and
  Zhou together, with a concentration-order rival; if neither closes Cao at 0.4 M and Zhou at 20 mM, the supply's form
  is wrong and the misfit protocol's mechanism step applies.
- **Casein systems** stay diagnostic until the lane carries protein-bound lysine and the caseinate's buffer.
## 4b. Design (step 4) -- decisions taken from the first reading (2026-10-29)

- **Rows without a printed error** (Davidek 2002's NaOH uptake, 49 values; Chen 2000 / 2002 and Xu 2008 final pH):
  no error is invented. Proposed: each such series gets its own noise scale, estimated from the data with the other
  constants (profiled in the likelihood), so it weighs as much as its own scatter allows. To be built in
  `fitting.py` and tested on synthetic data before use; until then these rows stay out.
- **Observation operators** the fit must compute, not approximate: the sum of formic + acetic acid (`FA+AA`); titratable
  acid to pH 8.3 (`TA_to_pH8.3`, a titration of the predicted state); the pH of a 6x dilution read at 25 C
  (`pH_6x_dilution_25C`, Zhao 2017b: the charge balance at the diluted composition). Each is a function of the
  predicted state, written once, tested against a hand calculation.
- **Out of scope here:** lactose (Nguyen 2016's M2 arm) and the di-ARP of lysine (Zhang 2025) are not engine species;
  their rows are kept in the table and reported, not fitted.
- **Laboratories** (leave-one-out groups) from `data/keys/laboratories.yml`: Wageningen (205 rows), Nestle (46),
  BTBU Xie group (18), Jiangnan Zhang X group (24). Wageningen dominates; the done-test's "three of four" therefore
  also requires the fit to hold with Wageningen left out.
- **Duplicates:** Martins 2003 Part II replots the 2003b acids (CF-018), and the trunk's measured acid constants come
  from those data; the rival arm must not count them twice.
## 5. Design review (step 5) -- done 2026-10-29, before any fit

An independent reviewer (read-only, the prereg, the fit script, the probe, the engine changes, the table) raised
seven BLOCKING issues, nine SHOULD-FIX and notes. Every blocking one is addressed below; the fit as first written was
not run (a two-iteration smoke test checked mechanics only and was deleted).

| # | issue | response |
| --- | --- | --- |
| B1 | priors from Martins 2003 Table 3 are fitted to the 2003b rows in this table: data used twice | priors dropped; slopes FLAT in [0, 1.5] |
| B2 | Martins 2005's trunk constants average the drift from pH 6.8 ("concealed in the rate constants"); sloping them along the drift counts it twice | a free log10 offset per sloped step (the constant-pH value at the reference is fitted, not the tabulated one) |
| B3 | a pH-stat run (Davidek, held at 90 C) handed the pH tables its hot set-point, drifting pots their 25 C pH | `integrate_sulfur(clamp_ph_label=...)`: the held state's 25 C pH, computed from the SID that holds the set-point hot |
| B4 | glycine, DFG and the Schiff base not in the charge balance (the module's own gap list: "a defect the day the trunk gets a pH state") | `ph_state.TRUNK_SOLUTES_IN_BALANCE` (OFF by default): glycine at GKL's constants; AMA and SB as the Amadori-compound class (as ARP and CAMA) |
| B5 | Davidek's NaOH uptake counts acids the lane cannot make (NaOH / (FA+AA) = 2.2 at pH 5 to 1.1 at pH 8) | NaOH uptake diagnostic only (reported at the optimum, not fitted) |
| B6 | Zhao 2017b's hexanal / decadienal pots were fitted as blanks (the aldehyde binds cysteine) | excluded, as heptenal |
| B7 | uncentred per-series noise turned bias into noise; leave-one-out kept the held-out lab's noise: done-tests pass by construction | noise per data set x observable from scatter WITHIN each series around its own mean; re-estimated on the kept rows in each fold; done-tests on RAW residuals against the tolerances below |
| S1 | one shared Amadori slope cannot fit the pattern | per-step slopes for k_ama_tdg and k_ama_odg (k_ama_mgo tied to k_ama_odg: no row here sees MGO), the formic exit separate |
| S3 | slopes could not fail identification (log10k threshold 3) | new kinds in `fitting`: ph_slope (0.25 dec / pH), log10_offset (0.5 dex); B8's table untouched |
| S4 | k_aa_frag's Ea pinned at its bound; A/B compared by laboratory | Ea free in arm A (20-160 kJ/mol); A vs B read per data set |
| S5 | done-test 2 needs buffer / sugar / amine attributes | added to the diagnostics; phosphate catalysis (Davidek pH 5: 3.3 vs 11 mmol/L) pre-registered as an expected flag; done-test 4 (carbon) cannot be run on this table and is stated as not run; done-test 6's pathway-yield rival deferred to the pentose block's follow-up |
| S6 | one start | two starts (`fitting.two_starts`, the second on the other side of every band) |
| S7 | inline vs fixed point not checked at the fitted settings | re-checked with both switches on: max 7e-4 over 409 rows; the optimum is re-scored on the fixed point |
| S8 | failures silent | failed pot runs counted and reported |
| S9 | thread oversubscription | BLAS / numba threads pinned to 1, niced, 9 workers (owner: the machine is free overnight) |

**Declared before the fit (done-test 1 and 5, raw residuals):** a data set x observable passes when its RMS raw
residual is within 0.30 pH units (pH, diluted pH) or 0.30 dex (x2) for amounts. Done-test 1 passes when at least
three of four laboratories pass on every data set they hold, in their leave-one-laboratory-out fold.

**Hold-outs (done-test 5), one per block, from data not looked at:**
- hexose block: Martins 2005b Fig. 5 (constant pH, pH-stat) and Fig. 7 at initial pH 4.8 / 5.5 / 6.0 / 7.5 (formic
  and acetic acid, 100 C). Being read now by two readers; the fitted model's curves on a dense time grid are frozen
  to a file before their rows are compared. The pH-stat series tests B2 directly (a constant-pH reference).
  DISCLOSURE: the second reader's report quoted four hold-out values (pH 7.5 and the pH-stat, formic and acetic acid
  at 120 min) to the lead before the predictions were frozen. No design choice depended on them: the fit's code and
  this section were committed before the fits started, and the frozen curves come mechanically from the fitted
  vectors.
- pentose block: none unseen that the lane can simulate (Chen 2000 / 2002's oven heating has no stated internal
  temperature). The BTBU (Cao, Zhao) fold is reported as "seen at the probe", not as a blind test.

**Expected result (written now):** arm B (no acetic-acid sink) wins on the DFG data and loses less than 0.1 dex on
Martins 2005 once the reference offsets are free (P = 0.55); the Amadori slopes land inside Martins' 0.37-0.92 band
(P = 0.6); the pentose acid yield stays poorly identified or on a bound (P = 0.6), which would send the 0.4 M
over-acidification to the mechanism step; done-test 1 passes (P = 0.3).
## Results (2026-10-30, overnight; arm A still running at 04:20)

Fits: `kinetic_core_ph_acid_submodel_fit_{B,C}.md` (arm B's pre-amendment run kept as `..._B_superseded_amendment1`).
Hold-out (frozen per arm, committed before scoring): `kinetic_core_ph_acid_submodel_holdout_score.md`.

| | arm B (sink off) | arm C (+ glucose-step slopes) |
| --- | --- | --- |
| constants on a bound | 2 Amadori slopes (1.5), acid yield (1) | condensation slope (1.5), acid yield (1) |
| Amadori slopes vs Martins' 0.37-0.92 | outside | inside (0.68 / 0.90 / 0.79) |
| residual trend with pH (in-sample) | rho +0.40, p 2e-13 | rho 0.005, p 0.92 |
| done-test 1 in-sample (10 data set x observable) | 4 within | 4 within |
| leave-one-laboratory-out | 0 of 4 laboratories pass | Nestle 1/1, Wageningen 3/6, others 0 |
| HOLD-OUT acetic (5 runs) | 5 within (RMS 0.07-0.27) | 1 within; pH 4.8 -2.34 dex, 5.5 -1.42 |
| HOLD-OUT formic (4 runs) | 0 within (x2 to x10 low, worst at low pH) | 1 within (pH 7.5); pH 5.5 -1.94 dex |
| pentose block (Cao, Zhao, Zhou) | unchanged, outside | unchanged, outside |

**Reading.** Against the expected result (sec. 5): arm C did bring the Amadori slopes inside the measured band and
remove the pH trend (both predicted at P = 0.5), and the pentose block did not move (P = 0.95) -- but only by putting
the condensation slope on its bound, three times Martins' measured 0.50, and the hold-out rejects that: the
low-pH drifting runs lose 10-200x of their acid. The in-sample ladder rows are pH only; a pH can be reproduced by
"less acid made" or "less glucose consumed", and only acid amounts across pH separate them. Arm B transfers acetic
acid to every hold-out run, including the constant-pH mode the fit never saw, and fails formic acid at low pH (its
formic-exit slope sits on the bound). The done-test fails for both arms; no switch is proposed to ship.

**What it changes (next round, owner's call).** (1) The hold-out is spent for B and C: Martins 2005b Fig. 7's acids
(formic and acetic at four pH) become FIT data for the next round, which needs a new hold-out (candidates: Knol
2010; Brands 2002's ladder if read). (2) The condensation's pH dependence must be identified from acid amounts, not
from pH. (3) The 3-DG -> formic exit is not a single pH-sloped step on this evidence (Martins 2005b's own Table 3
gives it pD 0.09 +/- 1.10): literature first. (4) The pentose block waits for sub-model 3's `k_dpo_ddp`.

## Amendments

**Amendment 1 (2026-10-30, 01:05; after arm B's first fit, before arm A's).** Arm B's diagnostics showed two flaws
in my implementation of review items B3 / B4, not in the chemistry: (i) DFG (AMA) and the Schiff base were titrated
as the pentose-ARP class, whose amine pKa is a CALIBRATED constant the fit also moves for Zhou's pentose pots (arm B
put it at 5.8): the pentose and hexose blocks were coupled through a number neither measures, and Davidek's
pH-stat pots then saw a distorted pH (their NaOH-uptake diagnostic fell to -22 dex). Fixed with DFG's MEASURED
constants (Mossine 1994: amine 8.18, carboxyl 2.20; glycine's enthalpies carried as the class value, declared).
(ii) The trunk slopes read the label-scale pH (review B3). Martins' slopes were measured in 0.1 M phosphate, whose
label sits within ~0.1-0.3 of the in-situ pH, so they are in-situ slopes; on the label scale an unbuffered DFG
pH-stat imports its own ~1-unit thermal shift. The trunk slopes now read the in-situ pH (B3's label conversion stays
for the other measured tables). Effect at the start vector: Davidek FA+AA +0.35 -> +0.16 dex mean (RMS 0.58 ->
0.31); NaOH uptake -0.14 -> -0.49 (diagnostic; it counts acids the lane cannot make, B5). Arm B's first result is
kept as `kinetic_core_ph_acid_submodel_fit_B_superseded_amendment1.{json,md}` (slopes on their 1.5 bound, acid
yield on its bound); both arms are re-run. Nothing else in sec. 5 changes; the tolerances stand.

**Amendment 2 (2026-10-30, 02:20; after arm B's re-run, before arm C runs): arm C.** Arm B (re-run) puts two of three
Amadori-decay slopes on their 1.5 bound (Martins measured <= 0.92) and leaves a residual trend with pH (rank
correlation +0.40, p = 2e-13: too low at pH 5-5.5, too high at 7-7.5) that differs between the DFG pots (pH-stat,
Martins 2003b) and the glucose + glycine pots (the ladder). Martins 2005b's Table 3 (fitted on the ladder) finds the
ladder's pH dependence on the GLUCOSE-CONSUMING steps -- the condensation (pD 0.50 +/- 0.14), glucose -> fructose
(0.88 +/- 0.27), the sugar's acids (0.54 +/- 0.36) -- which the lane leaves pH-independent; the DFG pots start after
those steps and cannot see them. Arm C = arm B + a free slope and reference offset on `k_schiff` and one shared free
slope and offset on `k_glc_fru` and `k_fru_acids` (12 constants; Table 3's values as starts only, never priors: it is
fitted on the ladder, review B1). Prediction written now: arm C brings the Amadori slopes inside 0.2-1.0 (P = 0.5),
removes the pH trend (P = 0.5), and does not change the pentose block (P = 0.95). Same tolerances and hold-outs.
