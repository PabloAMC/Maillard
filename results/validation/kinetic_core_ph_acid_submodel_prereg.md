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

## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- decisions taken from the first reading (2026-10-29); the rest TO DO

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
## 5. Design review (step 5) -- TO DO
## Amendments
