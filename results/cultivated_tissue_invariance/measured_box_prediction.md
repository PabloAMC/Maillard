# Does the kinetic layer change what to add to cultivated tissue?

*Generated 2026-10-09 by `src/cultivated_tissue_invariance.py`. Pre-registration: `results/validation/cultivated_tissue_invariance_prereg.md`. A PREDICTION RUN on the current box (9 sourced, 0 secondary, 5 stub ranges). It is not a resolution of the declared run's verdict (pre-registration section 7): the same statistics are reported against the same thresholds so the two runs can be read side by side.*

## Verdict, by programme

| programme | verdict | draws evaluated | engine refused | < 2 restorable | top(E)=top(N) | mean τ | disagreements | dominant reversal | its envelope survival |
|---|---|---|---|---|---|---|---|---|---|
| 100C_20min | **T1** | 175 / 200 | 0 | 25 | 49 % | 0.149 | 51 % | ribose over glucose (64 % of disagreements) | 93 % |
| 140C_5min | **T1** | 175 / 200 | 0 | 25 | 54 % | 0.240 | 46 % | ribose over glucose (74 % of disagreements) | 97 % |

Read against section 3 of the pre-registration: T1 needs ≥ 20 % disagreement, one reversal carrying ≥ 50 % of it, and ≥ 80 % envelope survival; T2 needs ≥ 90 % top agreement and mean τ ≥ 0.75; T3 is ≥ 50 % of draws refused. Anything else is indeterminate and is not a licence to build.

## Where each precursor lands

**100C_20min** (mean position, 1 = restore first; over the draws the engine answered)

| precursor | naive N | engine E |
|---|---|---|
| ribose | 1.580 | 1.185 |
| cysteine | 2.700 | 1.875 |
| thiamine | 2.506 | 2.129 |
| glucose | 1.593 | 2.389 |

Disagreements, engine's winner over naive winner: ribose over glucose ×57; cysteine over glucose ×13; thiamine over ribose ×7; thiamine over glucose ×5; ribose over thiamine ×4; cysteine over ribose ×3.

**140C_5min** (mean position, 1 = restore first; over the draws the engine answered)

| precursor | naive N | engine E |
|---|---|---|
| ribose | 1.580 | 1.105 |
| cysteine | 2.700 | 1.625 |
| thiamine | 2.506 | 2.659 |
| glucose | 1.593 | 2.253 |

Disagreements, engine's winner over naive winner: ribose over glucose ×59; cysteine over glucose ×13; ribose over thiamine ×4; cysteine over ribose ×3; cysteine over thiamine ×1.

## What the engine refused

Structural, before any draw:

- trunk arm (methional, pyrazines, furaneol on glucose + glycine): methional refused, pyrazines ~1e-13 ug/L, no water threshold on any target; no decision metric (A1)
- cultivated fat: no route from tissue lipid to any lane (A1)
- hexanal: needs a lipid carrier; no tissue-fat carrier exists (A1)
- leucine, IMP, ribose-5-phosphate: not species in any core lane (A2)

## Gap map: what has not been measured

Every stub in the box this run swept. A stub is a range with no measurement behind it; the engine's ranking depends on the cultivated-side values of the four rankable precursors, so a stub there is a gap in the answer, not only in the table.

| tissue | precursor | engine can rank it | swept range mM | what would close it |
|---|---|---|---|---|
| cultivated_muscle | ribose | yes | 0.01 to 2.0 | free ribose by GC-MS (oxime-TMS) or enzymatic assay on the washed construct extract; CE-MS panels do not carry it |
| cultivated_muscle | cysteine | yes | 0.02 to 0.5 | free cysteine by CE-TOFMS with thiol protection at extraction, as Muroya 2019 did for beef; report cystine alongside |
| cultivated_muscle | thiamine | yes | 0.0003 to 0.01 | thiamine and its phosphates by HPLC-fluorescence (thiochrome) on the same extract |
| cultivated_muscle | glucose | yes | 0.1 to 10.0 | free glucose by enzymatic assay or GC-MS on the same extract; state the harvest wash |
| cultivated_muscle | ribose-5-phosphate | no | 0.01 to 0.1 | on the CE-TOFMS panel with cysteine; it was quantified in beef by that method |

**One experiment closes the rankable gaps.** 4 of the four precursors the engine can rank have no published measurement in cultured muscle. Muroya 2019 quantified cysteine, ribose 5-phosphate, IMP and leucine in beef on one CE-TOFMS run; the same panel on washed cultured bovine myotubes, with free ribose and glucose by GC-MS or enzymatic assay and thiamine by thiochrome HPLC on the same extract, beside a beef sample handled identically, turns every cultivated stub into a sourced range in one campaign. Three biological replicates, two harvest washes (none; PBS), one ageing arm (24 h at 2 °C) to see whether the IMP-to-ribose route runs in a construct at all.

## The composition box that was swept

| tissue | precursor | lo mM | hi mM | label | dossiers | note |
|---|---|---|---|---|---|---|
| beef | ribose | 0.33 | 2.2 | sourced | koutsidis2008b, koutsidis2008a | Koutsidis 2008b Table 2: 0.25 mmol/kg at day 1 to 1.67 at day 21 (n = 16 steers, GC-MS), 0.33-2.2 mM at 75 % moisture; Koutsidis 2008a Table 1: 0.57-1.08 mmol/kg across 30 steers at 10 d (0.76-1.44 mM) sits inside; the cited 0.26 mg/g point (Aliani 2013 via Hwang 2026, 2.3 mM) sits at the top edge and is no longer a corner |
| beef | cysteine | 0.002 | 0.23 | sourced | muroya2019, koutsidis2008b, koutsidis2008a | Muroya 2019 Table 1: 1.6 nmol/g at D0 to 107 nmol/g at D14, n = 3 steers, CE-TOFMS (D0 at the detection floor); Koutsidis 2008b Table 4: 0.05-0.16 mmol/kg over 21 d, n = 16, GC-MS; Koutsidis 2008a Table 3: 0.05-0.17 mmol/kg across 30 steers at 10 d, whose top is the upper corner (0.23 mM); the span is ageing plus animal spread |
| beef | thiamine | 0.00044 | 0.004 | sourced | lombardiboccia2005 | Lombardi-Boccia 2005 Table 2: total thiamine 0.01-0.08 (+/- 0.01) mg/100 g across five raw beef cuts by HPLC after acid hydrolysis, lowest mean to highest mean + SD; not detected in any cut after cooking; about twofold below the stub |
| beef | glucose | 2.4 | 15.0 | sourced | bischof2023, koutsidis2008b, koutsidis2008a | Bischof 2023 Table 1 with the alpha- and beta-glucose rows SUMMED (the first read took one anomer row as the total): 4.43 +/- 2.61 to 10.01 +/- 1.42 umol/g wet across two breeds and 28 d, mean -/+ SD = 2.4-15 mM; Koutsidis 2008b Table 2 (7.33-10.3 mmol/kg, 9.8-13.7 mM) and 2008a Table 1 (6.94-10.6 mmol/kg across 30 steers) sit inside; the cited 1.48 mg/g (11 mM) too |
| beef | leucine | 0.29 | 3.2 | sourced | muroya2019, bischof2023, koutsidis2008b, koutsidis2008a | Muroya 2019: 263-827 nmol/g (0.35-1.1 mM); Bischof 2023 (corrected rows): 0.29 +/- 0.07 to 1.53 +/- 0.52 umol/g (0.29-2.7 mM); Koutsidis 2008b: 0.43-1.75 mmol/kg; Koutsidis 2008a: 0.78-2.40 mmol/kg across 30 steers (to 3.2 mM); unrankable by the engine |
| beef | IMP | 0.1 | 10.0 | sourced | muroya2019, bischof2023, koutsidis2008b, koutsidis2008a | Muroya 2019: 78 nmol/g pre-rigor to 7574 nmol/g at D1 (0.10-10 mM); Bischof 2023 (corrected rows, 1.1-4.8 mM), Koutsidis 2008b (3.5-8.4 mM) and 2008a (3.3-5.9 mM) inside it; unrankable by the engine |
| beef | ribose-5-phosphate | 0.005 | 0.1 | sourced | muroya2019, koutsidis2008b | Muroya 2019: non-detect at D0, 57-70 nmol/g at D1-D14 (to 0.093 mM); Koutsidis 2008b: 0.04 mmol/kg flat over 21 d (0.053 mM); lower corner set at 0.005 because the draw is log-uniform; unrankable by the engine |
| cultivated_muscle | ribose | 0.01 | 2.0 | stub | — | NO MEASUREMENT FOUND in cultured muscle of any species (read of 2026-09-14); range unchanged from the stub box |
| cultivated_muscle | cysteine | 0.02 | 0.5 | stub | — | NO MEASUREMENT FOUND: Joo 2022 prints cysteine only as a percent of total amino acids with no absolute total and no free/hydrolysed statement; Kim 2024b's free-amino-acid table omits it; range unchanged |
| cultivated_muscle | thiamine | 0.0003 | 0.01 | stub | — | NO MEASUREMENT FOUND in cultured muscle; DMEM carries ~12 uM thiamine HCl, what a washed construct retains is unknown; range unchanged |
| cultivated_muscle | glucose | 0.1 | 10.0 | stub | — | NO MEASUREMENT FOUND in cultured muscle; note Joo 2022 proliferated in glucose-free DMEM; range unchanged |
| cultivated_muscle | leucine | 0.1 | 1.0 | sourced | kim2024b | Kim 2024b Table 4: 35.8 mg/kg free leucine in a pig gelatin-scaffold construct at 90 % moisture (0.30 mM), CONFOUNDED by 49.6 mg/kg in the scaffold-only arm; one point, pig, widened tenfold; unrankable by the engine |
| cultivated_muscle | IMP | 0.0003 | 3.0 | sourced | kim2024b, joo2022 | two primaries four decades apart: Kim 2024b pig construct 0.11 mg/kg (0.00035 mM); Joo 2022 bovine 2D tissue 1.98 mmol/kg (2.6 mM); the box carries both corners; unrankable by the engine |
| cultivated_muscle | ribose-5-phosphate | 0.01 | 0.1 | stub | — | NO MEASUREMENT FOUND in cultured muscle; range unchanged; unrankable by the engine |

Fixed in every spec: pH 6.0, a_w 0.98, matrix water, phosphate 0.03 mol/L. 200 draws per programme, seed 0, 50 envelope draws per reversal.
