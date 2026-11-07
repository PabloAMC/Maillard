# Pre-registration: do the experiments we call decisive shrink the answers a formulator uses?

**Status: DRAFT v1 (2026-11-07), not frozen. DESIGN AND TIMING ONLY.** No simulation has run. The only runs made for
this file are one timed lever prediction and one timed JAX network update (sec. 9.1), and a declaration-only check of the
decision recipes (no integration; sec. 2.3). It freezes when its independent review is answered (sec. 10). After the
freeze, nothing in secs. 2 - 8 may change. A change after any number exists is an amendment, dated, with its reason, and
the run reports both versions.

*Plan: [STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 1 (what formulators need), sec. 5 W5.4 (the lever
table) and W6 (outside data), sec. 7 (items 4, 5, M1.1, S3, S4), sec. 8 (odds). Designs:
[decisive_measurements.yml](../../docs/reference/decisive_measurements.yml); the M4 case
[sec. 4.4 and 5.4](../../docs/reference/cases/M4_glucose_mft_protein_amine.md) and its power check
[kinetic_core_m4_reduced_power.md](kinetic_core_m4_reduced_power.md); the M5 case
[secs. 4.4 - 5.4](../../docs/reference/cases/M5_buffer_fat_water_activity.md) and
[m5_power_check.md](m5_power_check.md); the M3 case and
[kinetic_core_m3_strecker_layer_prereg.md](kinetic_core_m3_strecker_layer_prereg.md). The `rank` verb:
[src/experiment_value.py](../../src/experiment_value.py) and its output
[experiment_value_ranking.md](experiment_value_ranking.md). The judge:
[src/kinetic_core/predictive_densities.py](../../src/kinetic_core/predictive_densities.py).*

## 0. Why, and what this is not

**The worry (owner).** Suppose we ran every experiment the cases call decisive. We might still be unable to answer a
formulator's questions. The cases' power checks show only that each design settles its own structural question
(M4: four rival pairs at 0.70 - 0.87; M5: >= 0.88 - 0.96). No record checks whether settling those questions shrinks
the uncertainty on the lever answers of sec. 1 and W5.4. So neither the experiment list nor `rank` has been
validated against the deliverable.

**What this file fixes, before any number exists:**
- the decisions (sec. 2);
- the designs (sec. 3);
- the truths, the synthetic data and the update (secs. 4 - 5);
- the metrics (sec. 6) and the thresholds that answer the worry (sec. 7);
- the OAV hook (sec. 8);
- the cost (sec. 9).

**What it is not.**
- It is not a comparison of rivals. Nothing here credits or ships anything.
- It spends no data. No measured decisive, reserved or held-out value is read. The truths are drawn from bands,
  never from data.
- It is not a W6 request. The contrast design (sec. 3.3) is a yardstick for the on-record list, not a proposal to send.

**Seen by the author.** The author read the plan, the cases and their power records, `decisive_measurements.yml`,
PR #108's description and its sensitivity record (`ph_balance_s3_sensitivity.md` on branch `s3-ph-balance-refusal`),
`rank`'s current output, and the S1 record's tier counts. The author computed no engine output on any decision, except
the timed D01 lever (sec. 9.1), whose concentrations were not printed or read.

## 1. The questions and hypotheses

- **Q1 (the worry).** For each design on record, and for their union: what fraction of the formulator decisions does it
  make materially more certain?
- **Q2 (targeting).** At the same budget, does a design aimed straight at the decisions do better than the on-record
  list?
- **Q3 (`rank`'s first validation).** Does `rank`'s ordering of designs predict the shrinkage they realise?
- **Q4 (refusals).** Which of today's refusals would each design lift?

The run answers each with the pre-registered verdicts of sec. 7. The author's forecasts are given here so that the run
scores them too:
- P(the worry is confirmed, sec. 7.1) ~0.65;
- P(the contrast design beats the best on-record design by >= 5 decisions, sec. 7.2) ~0.55;
- P(`rank` validated, sec. 7.3) ~0.1.

## 2. The decisions

### 2.1 The base recipes

All recipes use the engine's own vocabulary; every concentration is in the reacting phase, in mM. "Pea" means
`matrix: pea_isolate`, `protein_g_per_l: 50`, `protein_type: pea_iso`: the engine charges the isolate's measured
reactive-site pools (`data/species/protein_matrices.yml`).
- The engine does not yet treat protein amines as a class (M4 part A, structure only). So pea acts through its sites and
  the matrix thresholds, not as an amine source.
- That limit is the engine's, and the decisions inherit it.

| id | recipe |
|---|---|
| BR-S | L-cysteine 10 + D-ribose 10, pea, unbuffered, label pH 5.5, 120 C, 60 min, a_w 0.98, sealed, air headspace |
| BR-G | BR-S with D-glucose 10 in place of ribose |
| BR-T | D-glucose 50 + glycine 50, pea, label pH 6.5, 120 C, 60 min, a_w 0.98 |
| BR-P | D-glucose 100 + L-phenylalanine 100, 0.5 M phosphate pH 7.0, 100 C, 60 min (the M3 layer's design class) |

**Why 120 C.** M5's declared S3 consequence (case sec. 4.5) refuses, above 130 C, the transfer of 0.5 M phosphate-fitted
constants to unbuffered foods. A 140 C base would refuse almost every decision for that reason alone, which would hide the
other questions. 140 C enters as a lever (D10) instead.

### 2.2 The decisions

Each decision's output is the **log10 ratio** of the odorant between arm B and arm A, at the end of the hold. D08 is a
ratio of ratios. No output is absolute, because H0 refuses absolute thiol levels for a new laboratory.

| id | recipe A -> B (the lever step) | output | note's group |
|---|---|---|---|
| D01 | BR-S -> cysteine 20 | MFT | meaty |
| D02 | BR-S -> cysteine 20 | FFT | roasty-meaty |
| D03 | BR-G -> cysteine 20 | MFT | meaty |
| D04 | BR-S -> BR-G (ribose -> glucose) | MFT | meaty |
| D05 | BR-S -> BR-G | FFT | roasty-meaty |
| D06 | BR-S -> label pH 5.0 (-0.5) | FFT | roasty-meaty |
| D07 | BR-S -> label pH 5.0 | MFT | meaty |
| D08 | BR-S -> label pH 6.5 (+1.0) | MFT / FFT (the meaty-roasty balance) | balance |
| D09 | BR-S -> 120 min | MFT | meaty |
| D10 | BR-S -> 140 C | FFT | roasty-meaty |
| D11 | BR-S -> + 0.2 M phosphate, label pH 5.5 | furfural | roasty |
| D12 | BR-S -> + 0.2 M phosphate | MFT | meaty |
| D13 | BR-S -> + linoleate hydroperoxide 1 mM (the oxidised-fat proxy) | MFT | meaty |
| D14 | BR-S + pea isolate as lipid carrier -> + linoleate hydroperoxide 1 mM | hexanal | off-note |
| D15 | BR-S -> + glutathione 2 mM (the yeast-extract proxy) | MFT | meaty |
| D16 | BR-S -> + thiamine 2 mM | MFT | meaty |
| D17 | BR-S -> a_w 0.85 (a drier product) | FFT | roasty-meaty |
| D18 | BR-T -> label pH 7.5 | 2,5-dimethylpyrazine | roasted |
| D19 | BR-P -> 120 min | phenylacetaldehyde | Strecker |

**The proxies, declared and unsourced.**
- **Fat.** The engine refuses "sunflower oil" by name (an unmapped precursor). Its lipid lane makes every product from a
  declared hydroperoxide pool. So the fat lever is a 1 mM linoleate-hydroperoxide charge.
- **Yeast extract.** "Yeast extract" is refused by name (W5.1 is not built). The proxy is its thiol carrier, 2 mM
  glutathione. W5.1's composition data will replace the dose.
- Both doses are declared, not sourced. They size the step, not the claim.

### 2.3 What is refused today

"Engine" is the shipped code, from the declaration-only check of every arm (`declare_envelope` and `axis_refusal`, no
integration; 2026-11-07). "Rule" means a refusal that a frozen case or a draft PR declares under S3 / S4 but that is not
yet in the code.

| id | engine today | H0 | refused or widened by rule (not yet in code) | liftable by data? |
|---|---|---|---|---|
| D01, D02, D15, D16 | answered (ratio, omega in its interval) | ratio: not refused | **PR #108 (draft):** exposed. A balance-solved unbuffered sulfur pot is refused per pot where the carried balance rivals B1 / B13 move it beyond omega or flip its sign. 41 % of pH-lever ratios and 61 % of buffer-lever ratios on balance-solved pots are dependent. Dose ratios: 263 / 378. **M5-B:** widened (0.5 M phosphate constants in an unbuffered food, <= 130 C) | yes: the union's sign and width (sec. 6 item 3) |
| D03, D04 | answered (the incumbent's C2 + C3 route) | not refused | **M4 (S3):** the incumbent's hexose MFT route is falsified with no rival credited, so glucose-MFT outputs are refused (M4 sec. 4.3). Also PR #108 and M5-B as above | yes, only if a hexose rival is credited (M1.1) |
| D05 | answered (FFT through `r_glc_fur`, not falsified) | not refused | PR #108; M5-B widened | yes |
| D06, D07, D08 | answered | not refused | **PR #108:** pH-lever ratios, the most exposed class (1085 / 2666 dependent, 323 sign flips). M5-B widened | yes |
| D09 | answered | not refused | PR #108; M5-B widened | yes |
| D10 | answered | not refused | **M5-B:** refused above 130 C (no kappa(T) band there). PR #108 | only by a phosphate-factor measurement above 130 C, which no design here makes |
| D11 | answered | not refused | **M5-B:** R0 falsified at pH 5 - 7. The union of R1 / R2c / R4' / R7 spans zero where R4' excludes the furfural steps, so the **direction is refused**. PR #108 | yes (M5-B separates R1 from R4') |
| D12 | answered | not refused | M5-B: widened, direction shown only where the rivals agree. PR #108 | yes |
| D13 | answered (sulfur + lipid co-integration) | not refused | **M5-F:** the direction is shown (every carried rival lowers MFT); the **size is refused** (capture has no band) | yes: a two-sided bound from the data (sec. 6 item 3) |
| D14 | answered (lipid lane; the null model beats it on held-out rows, sec. 2 of the plan) | -- | none | -- |
| D17 | **refused** (axis: the sulfur lane has no a_w term) | -- | M5-W: sulfur a_w comparisons stay refused; no design written | **no, structural**: needs an a_w term on the sulfur lane |
| D18 | answered (trunk; E4 pH band) | -- | M5-B: widened (pH 5 - 8) | yes |
| D19 | **refused** (no phenylalanine or phenylacetaldehyde species) | -- | M3 phase 1 is a forced-input layer: it cannot take a new recipe | **no, structural**: needs M3 phase 2 (owner first) |

- **H0 refuses none of these decisions**, since each is a ratio. H0 matters for the OAV hook (sec. 8), because an odour
  activity needs an absolute level.
- **Which refusals stand when the run starts is set at step 0** (sec. 5.6), by the rules in force that day. The table is
  the expectation, not the input.
- **Structural refusals (D17, D19)** cannot lift on data. They are kept as witnesses. Their metrics are "not computable"
  and they are excluded from the denominators of sec. 7.

## 3. The designs

### 3.1 The budget unit

The unit is the **heating**: one independently prepared and heated vessel, as each record counts it.
- The records do not all mean the same by it. M4's arms sample time courses within a heating; M1b's B1 uses one vial per
  time point.
- The contrast design (sec. 3.3) counts one sealed vial per cell and time, the more expensive reading. So budget equality
  is conservative for it.
- Analytical cost (labels, SIDA, trapping) is not equalised. It is reported per design.

### 3.2 The designs evaluated

| code | design (record) | heatings | informs (structural questions, sec. 4.1) | completed here? |
|---|---|---|---|---|
| G-M4 | **M4 minimum deciding design**: arm 2 (fed 1-DG vs fed glucosone, N2) + arm 3 (fed Cys-hexose Amadori vs the matched glucose + Cys pot), n = 3 (M4 sec. 5.4) | 96 | Q-hex | no |
| G-M4+ | M4 all three arms (variant; adds arm 1's positional labels) | 204 | Q-hex | no |
| G-M5B | M5 part B: 10 conditions (phosphate inert zero / 5 / 50 mM / 0.2 / 0.5 M x pH 5.5, 7.0) x 4 readouts, n = 3 (M5 sec. 5.4) | 30 | Q-buf | **temperature: 120 C** (the records state none; inside Guan's bounded 110 - 130 C) |
| G-M5F | M5 part F: 6 arms, both phases read, n = 3 | 18 | Q-fat | no |
| G-M5W | M5 part W: >= 3 moistures (10 - 70 % wb) x n = 3 (acrylamide) | 9 | Q-aw | no |
| G-M5 | M5 B + F + W | 57 | Q-buf, Q-fat, Q-aw | -- |
| G-M3 | M3 W6 pots: (a) fed ARP-Phe, ARP measured, argon vs air; (b) trapped fed ARP-Phe under argon, MGO / GO quantified, with the post-incubation control; (c) fed 1-DG + Phe, argon vs air | 18 | Q-str | **n = 3 per cell and the time grid 0, 10, 30, 60, 120 min** (records: "not stated"); (c)'s conditions: Hofmann 2000 Table 4's fed-donor class (0.1 M + Phe, pH 7, 98 C, 60 min) |
| G-M1c | `thiol_supply_structure`: arms 1 - 3, n = 3 | 150 | Q-sup | no (22 + 4 + 24 cells) |
| G-SPL | `pentose_split_vs_removal`: B1 (3 pH x in-situ / ex-situ OPD; centre point 4 replicate vials per mode) + B3 | 13 | Q-split | no |
| G-PHL | `ph_law_outside_trunk`: 4 recipes x 3 pH x n = 3 | 36 | Q-phlaw | no |
| G-LOSS | `thiol_loss_cysteine_dependence`: [2H3]-MFT at 0 and 30 mM Cys, 100 and 120 C, plus the mid-cook spike | 15 | Q-loss | **label pH 5.5, n = 3, spike at 120 C** |
| G-FUR | `furfural_fft_cysteine_contrast`: arms (ii), (iii) | 9 | Q-fft | **n = 3** |
| G-REC | **the on-record union**: G-M4 + G-M5 + G-M3 + G-M1c + G-SPL + G-PHL + G-LOSS + G-FUR | 394 | all of the above | -- |
| G-C | **the contrast design** (sec. 3.3) | 96 | every question the decisions touch, at once, in the matrix | defined here |
| G-C+REC | G-C + G-REC | 490 | -- | -- |

**Not simulated, with the reason.**
- `norfuraneol_decay_ph_law`, `methylglyoxal_self_loss`: the records fix no concentrations, temperatures or pH.
- `pyrazine_ring_aminoketone`, `thiamine_route_thiol_role`, `glycation_barrier_low_moisture`: "design not yet
  written".
- `thiol_absolute_second_lab`: absolute levels only, with no pots designed. It cannot move a ratio decision. It moves
  only the OAV hook's weights (sec. 8).

"Completed here" fills only what the run needs to generate data. It is not a design request.

### 3.3 The contrast design G-C

One laboratory. Pea isolate 50 g/L in water, unbuffered, label pH set with HCl / NaOH, 20 mL sealed vials with 10 mL
fill, air headspace. One vial per cell and time, n = 3 independent vials.

**Core factorial (24 cells), at 120 C:**
- cysteine 10 / 20 mM;
- sugar: ribose 10 or glucose 10 mM;
- label pH 5.0 / 5.5 / 6.5;
- time 60 / 120 min.

**Augmentation (8 cells),** at the BR-S and BR-G points (cysteine 10, pH 5.5, 60 min), one cell per sugar each:
- + 0.2 M phosphate at label pH 5.5;
- 140 C;
- + linoleate hydroperoxide 1 mM;
- + glutathione 2 mM.

**Size.** 32 cells x 3 = **96 heatings**, equal to G-M4.

**Readouts.**
- SIDA GC-MS: MFT, FFT, bis(2-methyl-3-furyl) disulfide (also TCEP-reduced total), furfural, 2-acetylthiazole.
- HS-SPME-GC-MS against an internal standard, ratios only: hexanal, nonanal, 2-pentylfuran, (E,E)-2,4-decadienal,
  2,5-dimethylpyrazine, methylpyrazine, phenylacetaldehyde, methional.
- Free cysteine (NEM at the quench) and the sugar (HPAEC).
- pH read cooled on every vial, and in situ where an electrode allows.
- LOD / LOQ per assay.

**What it covers.** It covers D01 - D15 directly. D16 (thiamine), D18 (pyrazine) and D19 it reaches only through shared
coordinates. That limit is part of the comparison.

### 3.4 The measurement model of every design (sec. 7 item 5)

One model, shared by the generator (sec. 5.1) and the update (sec. 5.3), per readout:

| term | value | applies to |
|---|---|---|
| run-to-run noise | the design's printed or requested CV per independent heating: SIDA 10 % (M4, M1c, G-C); M5-B 8 % (its sizing CV; 15 % sensitivity); M5-F 15 % (25 % sensitivity); M5-W 20 %; HS-SPME ratios 20 % (30 % sensitivity); a design that states none takes the same-class reproducibility its case cites, else 20 % | every level |
| model discrepancy | omega = 0.419 dex per value, shared by a cell's replicates, independent between cells (H0; 0.390 as a sensitivity). The same for every rival | every level and ratio |
| label shares | +-1 % absolute, plus a 1 % mixed-channel floor; the d3 / M+3 overlap as declared in M4 sec. 5.1 | labelled arms |
| trapping | the OPD artefact as a blank-subtracted term with its own variance (M4 sec. 5.1; M3b sec. 3.2) | trapped arms |
| floors | the request's LOD where stated; else a declared transfer of 1 % of the lowest predicted level in that design's table, in required agreement over x0.1 - x10. Censored densities (`censored_gaussian_log_density`) | every level |
| between laboratories | none: every design is one laboratory. H0's tau is never applied to a ratio | -- |
| pH | cooled and in-situ scales kept apart (`require_ph_scale`); a measured trajectory enters as data, never as a clamp in the scoring model | pH readouts |

**The misspecified generator.** This variant adds the same-observable reproducibility term, else the lane's within-lab
term, a second time beyond what the update assumes (sec. 7, false-credit check). Every metric is reported under it as a
sensitivity.

## 4. Truths

### 4.1 The structural questions and their carried rivals

A truth is one **world**: one rival per structural question, plus every continuous coordinate. Rivals are those carried
under S3: named in a case before results, and not falsified. The incumbent is a rival only where it is not falsified.

| question | carried rivals (equal weights unless sourced) | falsified, never a truth | decisions it moves (declared from structure, before any run) | runnable today? |
|---|---|---|---|---|
| Q-ph (M1b part (i), the pot's balance) | B1 (phi_cys = 1), B13 (best draws as PR #108; their S3 bands where they exist) | B0 + the shipped centre ledger (CF-084, CF-134) | D01 - D16 (unbuffered sulfur pots) | yes (NumPy, inline coupling) |
| Q-sup (M1c, thiol supply) | R0, C1, X, X0, R-rel, R-cal, R-sink, R-K | -- | D01 - D10, D12, D15 | R0, C1 yes; the rest with M1c's build |
| Q-loss (in-pot MFT loss) | the frozen-form sink banded; a temperature-only law; M1's C2 (free cysteine slows the loss) | -- | D01, D03, D07, D09, D10, D12, D13 | C2 yes (EX-061 switch); T-only with the joint-refit amendment |
| Q-route (M2, MFT route shares) | M2 v2's rivals, carried as a factor (M4 sec. 5.2, row "NF and DDP shares") | -- | D01, D04, D07, D09 | as a factor |
| Q-split (M1b part (ii)) | the enolisation split carried to pentoses; base-favoured 3-DP removal; base-catalysed NF decay (B61's switch); "deferred" (pH-flat) | -- | D05, D06, D08, D11 | B61 switch yes; the others with M1b(ii)'s build |
| Q-hex (M4 part G) | R1, R2, R4', R5, and the temperature-dependent R1 / R5 mixture (f at 145 C on [0, 1], dEa on [-60, +60] kJ/mol) | R0's and R3's C2 + C3 route; R4 as written | D03, D04 | with A2b (W3, M4 structure at parity) |
| Q-buf (M5 part B) | R1, R2c, R4', R7 | R0 at pH 5 - 7; R2 with K <= 0.3 M; R4 as written; R8 pure | D11, D12 (and every unbuffered pot, through the 0.5 M-fitted constants) | with M5's buffer switch (~0.5 day, M5 sec. 4.5) |
| Q-fat (M5 part F) | R2 (capture), R3 (partition + capture), R6 (oil-phase oxidation) | R0; partition alone; R4 as written | D13 (and D14 as a sensitivity) | **no**: the two-phase build is deferred until a decider exists. A reduced lever operator (sec. 4.3) stands in, built for this run only |
| Q-aw (M5 part W) | R0 inside its window, R1, R3, R4 | R0 and R2 above a_w ~0.9 (acrylamide lane) | none answered (D17 is structural) | -- |
| Q-str (M3) | R1 - R4 | -- | D19 only (structural) | the phase-1 layer only |
| Q-phlaw (M0a / M4 part A) | R2, R3, R5 of M0a | -- | D18 (trunk) | parked case; as M4 part A |
| Q-fft (CF-065) | the furfural-lifetime and the FFT-survival readings, carried as two bands on `k_fur_decay` and an FFT survival term | -- | D02, D05, D06, D10 | bands only |

**A rival not runnable on the run date is not dropped silently.** The decision is computed on the runnable subset and
flagged "envelope incomplete". Sec. 7's verdicts are computed twice: on all decisions (primary), and on the
complete-envelope decisions only (sensitivity).

### 4.2 The continuous coordinates

Every coordinate comes from the S1 record (`results/validation/incumbent_evidence_tiers.json`, 306 rows) or from a
case's declared bands:
- **measured, analogue:** the measurement's own interval (log-normal where the source prints an SD; uniform on the
  analogue band otherwise);
- **bounded:** uniform within the stated side(s); an open side is closed by the reaction-class distribution;
- **fitted-identified:** its S1 profile interval;
- **fitted-unidentified, with a band that is not a search bound** (93 rows): uniform on the band;
- **fitted-unidentified on a search bound** (26 rows), or a rival constant with no band (M5-F capture, R6's oxidant
  load): the reaction-class distribution where a class exists (B55's `ea_class_prior` for Ea; `reaction_families.yml`
  for log k). Otherwise a **simulation-only range** (the search bound), labelled "truth range not S4-sourced". Decisions
  whose prior-ensemble total-effect index on such coordinates exceeds 0.1 are flagged and reported separately;
- **conflicting:** the band spanning the measurements' own intervals;
- **structural, physical constant:** fixed. **Withdrawn:** excluded;
- **case bands** (M1c, M4 sec. 5.2, M5 sec. 5.2, M3 prereg): as declared there, both ends of every required-agreement
  range included as forms.

**The incumbent's point is never used.** B60's fitted value is not a truth, a band centre or a start. It enters only
where the S1 band's own source is a measurement that happens to agree with it.

**Joint draws.**
- Scrambled Sobol over the continuous coordinates.
- Rival indices stratified, so each rival of each question holds >= N / |R_q| worlds.
- Common random numbers: one world's theta serves both arms of every decision and every pot of every design.

### 4.3 The reduced lever operator for Q-fat (this run only)

Until M5's two-phase build exists, D13's rivals act as a first-order loss of water-phase MFT and FFT. It is
proportional to the hydroperoxide charge and is rival-specific:
- R2: MFT-selective;
- R6: O2-dependent;
- R3: partition (log K 1.0 - 2.5) plus R2.

The constants are drawn on a simulation-only range spanning the lipid falls of Elmore 2002 and Yiltirak 2026 (x2 - x17;
seen, spent data, used here only to scale the range). The operator is declared before any run. D13's results are
labelled "reduced operator".

## 5. Synthetic data, the update and its checks

### 5.1 Generating the data

- **Ensemble.** N = 2048 worlds (sec. 4).
- **Truths.** Every world serves as a truth once. Its data are scored against the other N - 1 worlds (leave-one-out),
  so truths are drawn outside the scoring ensemble, as the judge requires.
- **Per truth i and design g:** y_g(i) = the world's prediction on g's cells + omega discrepancy per cell + run noise
  per replicate, censored at the floors (sec. 3.4).
  - A union design (G-M5, G-REC, G-C+REC) takes its parts' data from the same truth world, so its likelihood is the
    product of theirs.
  - K = 1 noise draw for the metrics.
  - K = 4 for the nulls (sec. 5.4), so each rival as truth gives >= 1000 null replicates: >= 256 worlds x 4 at
    |R_q| <= 8.
- **Simulators.**
  - **Decisions and `rank`'s rows:** the one network (JAX), or the NumPy engine where a module is not yet ported, with
    each rival behind its switch.
  - **Design cells:** the case mini-models where they exist (M4: `src/m4_reduced_power.py`; M5:
    `scripts/diagnostics/m5_power_check.py`; M3: the phase-1 layer, `scripts/generators/fit_m3_strecker_layer.py`), the
    network otherwise.
- **Parity, before any metric.** On 64 worlds, the mini-model's design outputs must match the network's on the same cells
  within 0.1 dex (a quarter of omega). A design that fails is generated by the network instead, at the cost stated in
  sec. 9.

### 5.2 Choice of update: importance reweighting of the ensemble

**Primary: importance reweighting.** Each world j != i gets weight w_j(i, g) ∝ p(y_g(i) | theta_j), under sec. 3.4's
measurement model. The weights are Pareto-smoothed through the judge's own estimator (`importance_lppd`, `psis`), with
its acceptance rule (k-hat < 0.7, ESS >= 100). Why this and not Laplace:
1. **The rivals are discrete.** A world is a structure index plus coordinates. Laplace approximates one mode of one
   structure. Reweighting carries the mixture over structures, and it is the only form in which S3's union (sec. 6.1) and
   M1.1's T can be computed.
2. **The bands are bounded.** Uniform bands with mass at their edges, and a posterior pressed against a band wall (S4's
   band-edge criterion), are where a Gaussian at the mode fails.
3. **Lever outputs flip sign.** A Gaussian on the coordinates, pushed through a sign-flipping output, misstates P(sign).
4. **Cost.** The forward ensemble is paid once and reused for every design, truth and noise draw. Laplace pays a fit per
   (truth, design, rival). The timing run gives 49 s for a 5-coordinate, 1-pot fit (sec. 9.1). Scaled to a design's ~20
   pots and ~30 relevant coordinates, that is ~0.9 h per update, or ~10^5 core-h over 2048 truths x 15 designs x ~4
   rivals.

**Fallback for weight collapse.** An informative design collapses the weights. Where more than 20 % of a design's truths
fail the acceptance rule, an enrichment ensemble is drawn for that design:
- M = 512 worlds from a defensive mixture: 0.8 x a Gaussian mixture fitted to the pooled stage-1 posteriors of that
  design's truths, plus 0.2 x the prior;
- each forward-run on the decision, `rank` and design cells;
- combined with the prior ensemble by multiple importance sampling (balance heuristic, exact proposal densities).

A truth still failing after that is **unresolved**. It is never read as estimated:
- in sec. 7.1's test it counts as **shrinkage >= 25 %**, against confirming the worry;
- in sec. 7.3 it is excluded, with a sensitivity that counts it at the design's maximum.

### 5.3 The posterior quantities

For decision d, design g and truth i:
- **L_d(theta):** the lever log-ratio.
- **I_post:** the weighted 5 - 95 % quantiles of L_d (`weighted_quantile`): what the data could teach.
- **I_ship:** what the shipping rules would show:
  - the union, over the surviving rivals of each question d depends on, of each rival's own weighted 5 - 95 % interval;
  - the other questions' surviving rivals are held at equal mass inside each rival's interval (the per-question union
    rule). This understates the full cartesian union, by the same rule in prior and posterior. The direction of its bias
    on the shrinkage ratio is not known; it is reported, not corrected;
  - a rival survives unless it is falsified (posterior-predictive p < 0.05 after Holm across the design's units, within
    that rival's weighted draws);
  - if one rival is credited under M1.1 (sec. 5.4) against every survivor, I_ship is that rival's interval alone.
- **Both intervals have two parts:**
  - the **model part**, without omega: the only part a design can shrink;
  - the **total**, which convolves each world's lever with the shipped ratio-interval omega term (`thiol_relative_only`'s
    rule): what the formulator reads.

### 5.4 Credit inside the simulation (M1.1)

- **T.** For rivals a and b of one question, T is the pooled log predictive density difference on the design's decisive
  units, from within-rival weights.
- **Null.** T's null under b is the set of T values from b-truths, K = 4 noise draws each, with its 95th percentile's
  Monte Carlo error reported.
- **Credit.** a is credited if T exceeds that percentile and a is adequate on every decisive unit (p >= 0.05 after Holm).
  b by the mirror rule. Both, or neither: undecided.
- **Units.** A unit's role (credit / falsification-only / not decisive) follows sec. 7 item 5, through `UnitNoiseModel`.
  A design's absolute levels are falsification-only, since no design here prints a run-to-run SD of the same
  observable that the update assumes.

**Declared departure from the case judge.** The prior is the S4 bands, not conditioned on each case's spent rows (the
judge uses p(y_dec | y_spent)). The prior is therefore wider than the judge's.
- The designs look more informative here than they would in the case: **the departure favours the designs, so it biases
  against confirming the worry.**
- It also lets one simulation serve every design on equal terms.

### 5.5 How the update's accuracy is checked (all before any metric is read)

| check | what | pass |
|---|---|---|
| A1 per update | k-hat and ESS of every (truth, design) weight set, after enrichment | k-hat < 0.7 and ESS >= 100; the unresolved fraction is reported per design |
| A2 calibration (simulation-based) | over the truths, the coverage of I_post (model part) for the true L_d, and the uniformity of the truth's rank in its posterior (10 bins, chi-square, Holm across d x g cells) | coverage in [0.85, 0.95] and Holm p >= 0.01, in >= 80 % of cells |
| A3 Laplace cross-check (the JAX fit path) | G-C and G-M4; 10 truths each; one rival; a sub-problem of <= 5 coordinates and <= 4 pots, the rest held at the truth. `fit.fit` with the exact Jacobian, then `fit.laplace`, against reweighting a dedicated 1024-draw sub-ensemble on the same sub-problem | posterior mean within 0.25 posterior sd and sd ratio within [0.8, 1.25] in >= 8 of 10 truths per design. A Laplace fit pressed against a band wall is reported, not counted |
| A4 exact quadrature | a 2-coordinate sub-problem (`quadrature_lppd`) against reweighting | lppd within 0.05 nats |
| A5 mini-model parity | sec. 5.1 | 0.1 dex |

**If A2 fails** in more than 20 % of cells, or A3 in either design, the run reports "update not validated" and no verdict
of sec. 7.

### 5.6 Run order

0. **Freeze.** This file's commit is recorded. Then **step 0**: the refusal state of every decision in the prior
   ensemble, under the rules in force that day. PR #108's rule is used if merged, else its draft.
1. **Builds and parity.** List what is runnable (sec. 4.1); the parity checks (A5).
2. **Pilot.** 32 worlds, timed. The cost is re-projected (sec. 9.3).
3. The ensemble.
4. The synthetic data.
5. The updates and enrichment; A1.
6. Checks A2 - A4.
7. Metrics and verdicts.
8. Write-up, independent adversarial review (sec. 7 item 8), EX entry.

No metric is computed before step 6 passes.

## 6. Metrics, per decision and design

All are means over the truths, with a 90 % bootstrap interval over truths (1000 resamples). "Prior" is the same
quantity under the ensemble with no data.

1. **Expected shrinkage:** s(d, g) = E_i[1 - width(I_ship,model | y_g(i)) / width(I_ship,model)]. Reported with the
   median and with:
   - the I_post version;
   - the total-interval version;
   - the omega-dominance flag: the model-part width is below the omega-part width, so even full shrinkage barely moves
     what the formulator reads.
2. **P(sign correct).**
   - The true sign is sign(L_d(theta_i)), read only beyond 0.01 dex (PR #108's resolution). Truths inside that are
     excluded from the denominator.
   - Reported as three probabilities: P(shown and correct), P(shown and wrong) and P(refused).
   - "Shown" means S3's rule: every surviving rival agrees in sign.
3. **P(the refusal lifts),** for each data-liftable refusal of sec. 2.3. The probability that, after the update, the
   rule that refused d no longer does:
   - **PR #108:** the surviving balance rivals agree in sign and their union is computable (its lift clause assumed met
     on CF-084's controls; the design moves only the union);
   - **M4:** a hexose rival credited against every survivor (sec. 5.4);
   - **M5-B direction (D11):** the survivors agree in sign;
   - **M5-F size (D13):** the capture constant's 90 % posterior interval lies strictly inside its simulation-only range
     (a two-sided bound, which would give S4 a band);
   - **M5-B above 130 C (D10):** only a phosphate-factor measurement above 130 C lifts it. None of the designs makes one,
     so P = 0 by construction, reported as such.
4. **`rank`'s ordering and its change.**
   - (a) **V_rank(g):** the sum of `rank`'s VoI scores over the rows design g targets (frozen mapping, sec. 6.1).
   - (b) The Kendall distance between `rank`'s row ordering recomputed under the prior ensemble's envelope and under
     each posterior. Rows refused under step 0's rules drop, as `rank` drops them.
   - (c) V_real(g): the mean of s(d, g) over the decisions.
   - (d) The Spearman correlation between V_rank and V_real across designs (sec. 7.3).

### 6.1 The frozen `rank` mapping (today's 24 rows, `experiment_value_ranking.md`)

**Strict rule.** A design targets a row when it reads the row's compound in the row's lane. **System rule**
(sensitivity): the design also targets a row when it re-runs the row's own pot or precursor system, whatever it reads.

| design | rows targeted (strict) | V_rank strict | added under the system rule |
|---|---|---|---|
| G-M4, G-M4+ | none (no 2-acetylthiazole read; its trunk-lane dicarbonyls sit in sulfur-lane pots) | 0 | -- |
| G-M5B | 23 (3-DG, trunk) | 0.29 | 24 (MGO) if read |
| G-M5F | 9, 12 (soy HME hexanal, nonanal), 15 - 17 (pea UHT decadienal, hexanal, nonanal) | 6.13 | 1, 18 (2-pentylfuran) |
| G-M5W | 4, 8, 13, 14, 22 (acrylamide) | 7.81 | -- |
| G-M3 | none | 0 | -- |
| G-M1c | none | 0 | 3 (arm 2 re-runs Zhang 2025b's pot): 4.00 |
| G-SPL, G-LOSS, G-FUR | none | 0 | -- |
| G-PHL | 2, 5 (DMHF), 19, 20, 7, 21 (HMF), 4, 8, 13, 14, 22 (acrylamide), 23 (3-DG), 24 (MGO) | 20.11 | -- |
| G-C | 3 (2-acetylthiazole), 15 - 18 (pea UHT), 1, 9, 12 (soy HME) | 20.88 | -- |
| G-REC | the union of its parts (the G-M5B and G-M5W rows lie inside G-PHL's) | 26.24 | -- |

**Already true, before any run.** The single designs are the twelve of sec. 3.2 that are not unions: G-M4, G-M4+,
G-M5B, G-M5F, G-M5W, G-M3, G-M1c, G-SPL, G-PHL, G-LOSS, G-FUR and G-C. Under the strict rule, `rank` gives zero
value to 7 of the 12 (G-M4, G-M4+, G-M3, G-M1c, G-SPL, G-LOSS, G-FUR). This follows from the mapping alone: `rank` reads absolute-level
misses on its 24 rows, and only one of them (2-acetylthiazole) is a meaty thiol-lane row. It is recorded here as a
coverage finding, not as a result of the simulation.

## 7. Pre-registered thresholds and verdicts

**The denominator** is the 17 decisions whose lever interval is computable: D01 - D16 and D18. D17 and D19 are
structural.

**Why these thresholds.**
- A 25 % narrower interval is a variance ratio of 0.56. That is worth ~1.8 times the information already held: the
  smallest change that would alter a formulator's next experiment, rather than polish an interval.
- "At least 50 % of decisions" means a majority of the lever table.
- Both numbers come from the owner's statement of the worry, and are fixed here before any number exists.

### 7.1 The worry (Q1): primary test

- **Decision-relevant (primary):** s(d, g) >= 0.25 (sec. 6 item 1, I_ship model part).
- **WORRY CONFIRMED** if every single on-record design (G-M4, G-M4+, G-M5B, G-M5F, G-M5W, G-M5, G-M3, G-M1c, G-SPL,
  G-PHL, G-LOSS, G-FUR) and the on-record union G-REC are each decision-relevant for **fewer than 9 of the 17**
  decisions.
- **WORRY REFUTED** if any of them reaches >= 9.
- **Robustness.** A cell is borderline when its bootstrap interval spans 0.25. A verdict stands only if it survives
  flipping every borderline cell against itself. Otherwise the verdict is **UNDECIDED**, with the borderline count.
- **Secondary (the composite).** d is decision-relevant for g if any of these holds:
  - s >= 0.25;
  - P(the refusal lifts) >= 0.5 for a data-liftable refusal;
  - P(shown and correct) rises by >= 0.25 over the prior.

  In every case P(shown and wrong) must stay <= 0.10. Above that the design is **harmful** for d: confident and wrong.
  The worry test is repeated on the composite and reported beside the primary. The primary decides.
- **Also reported:** the same verdicts on I_post, on the total interval, on the complete-envelope decisions, under the
  misspecified generator, and at omega 0.390.

### 7.2 Targeting (Q2)

- **LIST MIS-TARGETED** if G-C is decision-relevant for >= 5 more decisions than the best single on-record design with
  a budget <= 96 heatings (G-M4, G-M5 and its parts, G-M3, G-SPL, G-PHL, G-LOSS, G-FUR), and the bootstrap 90 %
  interval of that difference excludes 0.
- **Not mis-targeted** if the difference is <= 0 with its interval excluding +5. Otherwise inconclusive.
- **Reported beside it:**
  - G-C against G-REC (394 heatings, four times the budget);
  - what G-REC adds to G-C (G-C+REC vs G-C): the value of the mechanism experiments once the contrasts are known.

### 7.3 `rank`'s first validation (Q3)

- **Coverage.** If V_rank(g) = 0 for at least half of the single designs, `rank` **CANNOT RANK THE DESIGNS THAT DECIDE
  ITS OWN QUESTIONS**. By the mapping (sec. 6.1) this outcome is expected; the run confirms it under the system rule
  too.
- **Ordering.** The Spearman rho between V_rank and V_real over the 12 single designs (ties at mid-rank).
  - **VALIDATED** if rho >= the one-sided 5 % permutation critical value (~0.50 at n = 12).
  - **CONTRADICTED** if rho <= 0 and the design with the largest V_real has V_rank = 0 or the lowest V_rank.
  - Otherwise **INCONCLUSIVE**.
  - Repeated on the five designs with V_rank > 0 (G-M5B, G-M5F, G-M5W, G-PHL, G-C). At n = 5 only rho >= 0.9 reaches
    one-sided p < 0.05; this is stated so the small-n result is not over-read.
- **The ordering change** (sec. 6 item 4b) is reported per design, with no threshold.

### 7.4 Refusals (Q4)

P(the refusal lifts) for each pair of a data-liftable refusal and a design, as a table. The structural refusals (D17,
D19) and D10's above-130 C refusal are listed with the build or measurement that would lift them. No threshold: this
section is the record W6 and the lever table read.

## 8. The OAV hook (perception layer; no thresholds added now)

The owner has put a perception layer in scope: odour activity values from measured thresholds. This section fixes how
the decisions will be re-weighted once those thresholds exist, so that adding them changes the inputs, not the rules.

**The weight of decision d:**

    w_d = sum_j w_j * DR( max(C_A,d(theta_j), C_B,d(theta_j)) * 10^(tau_j) / ODT_m(d) )

where:
- **DR** is `rank`'s own `_decision_relevance` (log10 of the ratio to threshold, clipped to [0.5, 5]). `rank` and the
  decisions then share one definition of relevance.
- **C_A, C_B** are the world's absolute levels in the two arms.
- **tau_j** is a draw from H0's between-laboratory term (0.964 dex) for thiol-lane species. That keeps the weight honest
  for a new laboratory.
  - The absolute level is **integrated, never shown**: H0's refusal of absolute thiol levels stands.
  - A laboratory with a `calibrate` overlay uses its calibrated level instead, with no tau.
- **ODT_m(d)** is the **measured** threshold of the output odorant in the decision's matrix (pea isolate, aqueous where
  the perception layer has no matrix value). It is log-normal with the source's reported spread, else the spread between
  sources for that compound and matrix class. Never borrowed from another matrix, as the engine's rule. A decision with no
  measured threshold gets no weight: it is excluded from the weighted metrics and counted.
- **Ratio of ratios** (D08): w = the larger of the two odorants' weights, since the more potent one dominates the balance.
- w_j are the posterior weights of sec. 5.2. The prior weight uses the ensemble.

**What changes when thresholds arrive.**
- The weighted fraction of decision-relevant decisions, sum_d w_d 1[s(d, g) >= 0.25] / sum_d w_d, re-runs sec. 7.1's
  test (>= 0.5 refutes) and sec. 7.2's (on weighted counts).
- Each decision is also flagged sub-threshold when P(OAV < 1 in both arms) > 0.9, and threshold-crossing when
  P(the lever moves OAV across 1) > 0.5.
- The note groups of sec. 2.2 (meaty, roasty, off-note, roasted, Strecker, balance) report their own fractions.
- No re-run of the ensemble is needed. The weights read the stored absolute levels of every world.

**What is not done now.** No threshold is added, read or chosen. Until the perception layer lands, w_d = 1 for every d.
The hook is a function signature and the stored per-world absolute levels. Its threshold file is named by the perception
layer's own PR and read through `src/data_access`.

## 9. Cost

### 9.1 The timing runs (2026-11-07, this machine, niced, threads 1)

| run | env | measured | load (1 min) |
|---|---|---|---|
| one lever prediction: D01, two arms through `engine.predict`, shipped constants (B0 fixed-point coupling) | `maillard` (native) | 5.78 s (arm A 3.16 s incl. warm-up, arm B 2.62 s) | 3.7 |
| one JAX network update: the composed sulfur network, pot `thiamine_cys_glucose_120C_Bolton1994`, 5 coordinates, 16 log10 observables (`build_jax_network_join._fit_problem`) | `maillard-jax` (jax 0.12, float64, CPU) | residual 1.13 s warm (5.73 s first); exact Jacobian 3.78 s warm (12.7 s first); **fit from a displaced start 49.0 s** (7 residual + 7 Jacobian evaluations, converged, cost 28.1 -> 0.89); Laplace 0.002 s; batch predictor **2.71 s per draw and pot** (8 draws 21.7 s; first call 11.2 s) | 3.0 rising to 11.9 (other jobs) |

**Reference, not re-run:** PR #108's record gives 0.41 - 0.43 s per pot for the inline coupling (B1, B13), against
3.14 s for B0's fixed point. The truths never use B0 (falsified), so the sulfur pots' cost lies between the inline
NumPy figure and the JAX predictor's.

**The timings are PROVISIONAL.** One shared machine whose load rose during the JAX run; they are upper-biased.

### 9.2 The projection

Per world:
- **136 network integrations:**
  - 16 distinct decision arms (D17, D19 are not run);
  - 11 `rank` benchmark pots (24 rows);
  - 109 network design cells: G-M1c 50, G-SPL 7, G-LOSS 5, G-FUR 3, G-PHL 12, G-C 32;
- **plus ~60 mini-model cells** (G-M4 32, G-M5 19, G-M3 6) at <= 0.1 s.

Per-integration cost: central 1.5 s, low 0.42 s (inline NumPy), high 2.9 s (the B0-like engine cost; the JAX predictor
measured 2.7 s).

| item | central | low | high |
|---|---|---|---|
| pilot, 32 worlds | 1.9 | 0.6 | 3.6 |
| ensemble, 2048 worlds | 119.5 | 35.9 | 227.8 |
| enrichment (central: 3 designs; high: 6) x 512 worlds x ~47 integrations | 30.1 | 8.4 | 116.2 |
| updates, nulls, metrics (array work) | 2.0 | 2.0 | 2.0 |
| A3 Laplace check: 20 updates of ~4 pots (~3.3 min each) plus a 2 x 1024-draw sub-ensemble | 4.5 | 2.1 | 7.7 |
| A4 quadrature | 0.2 | 0.2 | 0.2 |
| contingency (10 %) | 15.8 | 4.9 | 35.8 |
| **total, core-hours** | **174** | **54** | **393** |
| **wall-clock at 2 processes** (this task's limit) | **87 h** | **27 h** | **197 h** |
| **wall-clock at 8 cores** (the owner's standing cap) | **22 h** | **7 h** | **49 h** |

**Attended work: 28 - 40 h.**
- The harness: world sampler from the S1 record and the case bands, rival switch registry, generators per design, PSIS
  update, metrics, `rank` mapping: 14 - 20 h.
- Parity checks: 3 - 4 h.
- Q-fat's reduced operator: 4 - 6 h.
- Pilot and checks A1 - A4: 2 - 3 h.
- Write-up: 3 - 4 h.
- Answers to the review: 2 - 3 h.

**Not charged here.** Prerequisite builds already on the plan:
- A2b with M4's structure at parity (sec. 6, 11-21 to 12-05);
- M5's buffer switch;
- M1c's rival switches;
- the OPD trap module for G-SPL's in-situ arm. If it is not built, that arm is dropped and G-SPL is marked reduced.

The earliest start is therefore the week of 12-08.

### 9.3 The pre-declared reductions

- If the pilot projects more than 48 h of wall-clock at the cores granted, N falls to 1024 and M to 256. The central
  total falls to ~92 core-h: ~11.5 h at 8 cores, ~46 h at 2.
- Monte Carlo error then grows by ~sqrt 2. It is reported, and the borderline rule of sec. 7.1 absorbs it.
- Paid Colab (CPU, many cores) is an option to propose to the owner, not to assume: draws only, scored locally. JAX batch
  draws gave no CPU gain at equal tolerance (PR #74), so a GPU is out until a CPU-vs-GPU benchmark says otherwise.

## 10. Before the freeze

1. **An independent review** whose job is to break this design (sec. 7 item 8). In particular, it should check:
   - whether the decision set or the base recipes favour one design family;
   - whether the shrinkage threshold is the right reading of the worry;
   - the per-question union rule;
   - the departure in sec. 5.4;
   - the `rank` mapping.
2. **Answers to every objection** in this file's revision log.
3. **The freeze commit recorded here.** Then sec. 5.6's step 0.

## Revision log

- 2026-11-07, v1 DRAFT: written with the two timing runs and the declaration-only check. No simulation.
