# Pre-registration: do the experiments we call decisive shrink the answers a formulator uses?

**Status: FROZEN v3.1 (2026-11-07).** The freeze commit is the commit that set this header; its sha is recorded in
PR #110's body, since a commit cannot hold its own sha.
- **Run so far:** step 0 only (sec. 5.0), in separate commits after the freeze, written to
  `decision_relevance_step0.{json,md}`.
- **Not run:** the early variant and the full run.
- **Before the freeze,** the only runs were one timed lever prediction and one timed JAX network update (sec. 9.1),
  and a declaration-only check of the decision recipes (sec. 2.3).
- **History of the text:** v2 applied the premortem's edits 1 - 10. v3 answered the confirmation read of v2. v3.1
  answers the re-confirmation of v3 ("confirmed with edits", items 1 - 7). The changelog is at the end.
- **From here on,** nothing in secs. 2 - 8 and 12 may change. A change after any number exists is an amendment,
  dated, with its reason, and the run reports both versions.

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
- the cost (sec. 9);
- a partial early variant (sec. 12).

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
scores them too (revised in v3).
- **They are informed by seen summaries**, so they are not blind. The plan's sec. 2 gives the median within-study
  log-ratio error as 0.31 dex in sample and 0.66 dex held out (W1.1). That suggests a ratio discrepancy between ~0.45
  and ~1 dex, against 0.59 dex under iid omega.
- **step 0's ratio error sigma_r** (sec. 5.0.1), thiol lane: ~0.55 dex (90 % range 0.35 - 0.9). That forecast was
  made before stage 0c ran;
- P(fewer than 2 decisions are eligible, so the verdict is "CONFIRMED (omega-bound)", sec. 7.1) ~0.7;
- P(the worry is confirmed, either way) ~0.85;
- P(Q2 is informative and the contrast design wins, sec. 7.2) ~0.15;
- P(`rank` validated, sec. 7.3) ~0.05;
- the effective number of independent decisions (sec. 6 item 5): ~4 (90 % range 2.5 - 7).

**The verdicts assume a closed world:** the truth is one of the carried rivals, inside its bands (sec. 4.1). The
leave-one-rival-out run of sec. 5.7 reports how much a verdict leans on that.

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
| D13 | BR-S -> + linoleate hydroperoxide 1 mM (the ROOH lever) | MFT | meaty |
| D14 | BR-S + pea isolate as lipid carrier -> + linoleate hydroperoxide 1 mM | hexanal | off-note |
| D15 | BR-S -> + glutathione 2 mM (the GSH lever) | MFT | meaty |
| D16 | BR-S -> + thiamine 2 mM | MFT | meaty |
| D17 | BR-S -> a_w 0.85 (a drier product) | FFT | roasty-meaty |
| D18 | BR-T -> label pH 7.5 | 2,5-dimethylpyrazine | roasted |
| D19 | BR-P -> 120 min | phenylacetaldehyde | Strecker |

**Two levers are proxies (v2: relabelled, with their limits).** Both doses are declared, not sourced. They size the
step, not the claim.
- **D13 / D14 are the "ROOH lever", not a fat lever.** The engine refuses "sunflower oil" by name (an unmapped
  precursor), and its lipid lane makes every product from a declared hydroperoxide pool. The lever is a 1 mM
  linoleate-hydroperoxide charge.
  - At 10 % oil, 1 mM ROOH is ~10 mmol/kg oil, i.e. a peroxide value of ~20 meq O2/kg. That is twice Codex's
    <= 10 meq/kg for refined oils.
  - So it is a **rancid-fat lever**. It says nothing about adding fresh oil, oil type or oil level.
- **D15 is the "GSH lever", not a yeast-extract lever.** "Yeast extract" is refused by name (W5.1 is not built). The
  lever is 2 mM glutathione, yeast extract's thiol carrier.
  - It misses yeast extract's IMP -> ribose route, its free amino acids and thiamine. W5.1's composition data will
    replace it.

**What the base recipes leave out.** The 120 C sealed aqueous base ignores high-moisture extrusion, searing and crust
(W5.2), and heme (W5.3). The decisions inherit that scope.

**Decisions the engine cannot represent at all (a scope finding, v2).** These formulator levers have no decision here,
because the engine refuses them or has no term:
- fresh oil, oil type (saturated vs PUFA) and oil level;
- heme / leghemoglobin catalysis;
- yeast extract as an ingredient;
- protein-bound amines as a class (M4 part A);
- cooking method (pan, grill, extrusion);
- a_w in the sulfur lane (D17);
- Strecker aldehydes and methional (D19).

They are reported with the results as the part of the lever table that no design in this file can reach.

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

### 2.4 The off-grid decision set (v2; rebuilt in v3; primary for Q2)

**Why it exists.** The contrast design G-C (sec. 3.3) is built on the grid of sec. 2.2: the same cysteine, pH, time,
phosphate, ROOH and GSH steps. Scoring G-C on D01 - D16 would partly score its own cells, which teaches to the test.

**What v3 fixed.** v2's set overstated its claim: arm A of 9 of its 10 decisions was BR-S, itself a G-C cell, and
its B arms interpolated G-C's grid. v3 therefore uses a base that is no cell of any design.

**BR-X:** L-cysteine 12 + D-ribose 12, pea, unbuffered, label pH 5.8, 125 C, 45 min, a_w 0.98, sealed, air. Neither
its arms nor any B arm below is a cell of any design in sec. 3.

| id | recipe A -> B | output | relation to G-C's grid |
|---|---|---|---|
| O01 | BR-X -> cysteine 18 | MFT | inside the grid's span, off its cells (interpolation) |
| O02 | BR-X -> cysteine 18 | FFT | interpolation |
| O03 | BR-X -> label pH 6.3 | FFT | interpolation |
| O04 | BR-X -> 135 C | MFT | interpolation |
| O05 | BR-X -> 90 min | MFT | interpolation |
| O06 | BR-X -> + 0.1 M phosphate, label pH 5.8 | furfural | interpolation (G-C holds 0 and 0.2 M at one point) |
| O07 | BR-X -> + linoleate hydroperoxide 0.5 mM (ROOH lever) | MFT | interpolation (G-C holds 0 and 1 mM at one point) |
| O08 | BR-X with cysteine 5 -> cysteine 8 | MFT | outside the grid (below cysteine 10: extrapolation) |
| O09 | BR-X -> label pH 7.0 | FFT | outside the grid (above pH 6.5) |
| O10 | BR-X at 100 C, 90 min -> + glutathione 1 mM (GSH lever) | MFT | outside the grid (below 120 C) |

**What this still allows.** The interpolating seven still favour a design whose grid spans the lever space. That is a
legitimate property of a targeted design, so they are kept and reported apart from the extrapolating three (O08 - O10).

- No off-grid decision rests on the falsified hexose MFT route.
- Their refusals follow sec. 2.3's rows for the same lever. 135 C lies above M5-B's 130 C edge, so O04 carries D10's
  refusal.
- They add 12 distinct arms per world (sec. 9.2).
- Q1 is also reported on them, as a secondary set.

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
- + linoleate hydroperoxide 1 mM (the ROOH lever);
- + glutathione 2 mM (the GSH lever).

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
| model discrepancy delta (v3.1) | For the arms i of one study: delta_i = u + e_i. **e_i ~ N(0, sigma_r^2 / 2), independent per arm, with sigma_r the ratio error measured at stage 0c (sec. 5.0.1) for the arm's lane: the primary.** u is shared within the study, with variance max(omega^2 - sigma_r^2 / 2, 0) (omega = 0.419, H0); it cancels in every ratio. Sensitivities: sigma_r at stage 0c's leave-one-study-out minimum and maximum; and sigma_r in {0.593 (as shipped: iid omega), 0.42, 0.19}. Between designs, and between a design's cells and a decision's arms (another laboratory), everything is independent. The same for every rival | every level and ratio |
| laboratory offset (v3.1) | none in the primary. **Sensitivity (pre-registered check for any REFUTED verdict, sec. 7.1):** every design is a new laboratory, so each design's absolute-level readouts get one offset per observable, N(0, tau^2) with tau = 0.964 dex (H0's tau_ABS), shared by all its cells. It cancels in within-design ratios. A REFUTED verdict must hold with it | absolute levels |
| label shares | +-1 % absolute, plus a 1 % mixed-channel floor; the d3 / M+3 overlap as declared in M4 sec. 5.1 | labelled arms |
| trapping | the OPD artefact as a blank-subtracted term with its own variance (M4 sec. 5.1; M3b sec. 3.2) | trapped arms |
| floors | the request's LOD where stated; else a declared transfer of 1 % of the lowest predicted level in that design's table, in required agreement over x0.1 - x10. Censored densities (`censored_gaussian_log_density`) | every level |
| between laboratories | none: every design is one laboratory. H0's tau is never applied to a ratio | -- |
| pH | cooled and in-situ scales kept apart (`require_ph_scale`); a measured trajectory enters as data, never as a clamp in the scoring model | pH readouts |

**Why sigma_r decides the answer (v3.1).** A lever is a ratio of two arms of one study, so its discrepancy SD is
sigma_r. Under iid omega it is sqrt 2 x 0.419 = 0.593 dex, and a perfect design still shows only levers beyond ~0.98
dex (sec. 5.0.2). So sigma_r, not the designs, can decide whether any experiment helps.
- v2 assumed iid without a source. v3 measured a correlation rho against H0's omega.
- **v3.1 measures sigma_r directly** (sec. 5.0.1). omega cancels from it.
- rho = 1 - sigma_r^2 / (2 omega^2) is reported only as a derived diagnostic. It is not separately identified, since
  H0's omega is itself a within-group excess estimated on overlapping thiol rows.

**The shipped rule** (`thiol_relative_only`, iid omega on a ratio, 0.593 dex) is reported as the sensitivity "as
shipped today". Shipping the measured sigma_r would be its own sec. 7 change, not this file's.

**The misspecified generator.** This variant adds the same-observable reproducibility term, else the lane's within-lab
term, a second time beyond what the update assumes (sec. 7, false-credit check). Every metric is reported under it as a
sensitivity.

## 4. Truths

### 4.1 The structural questions and their carried rivals

A truth is one **world**: one rival per structural question, plus every continuous coordinate. Rivals are those carried
under S3: named in a case before results, and not falsified. The incumbent is a rival only where it is not falsified.

| question | carried rivals (equal weights, fixed in v3: no sourced rival weight exists on disk) | falsified, never a truth | decisions it moves (declared from structure, before any run) | runnable today? |
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
  whose prior-ensemble total-effect index on such coordinates exceeds 0.1 are flagged. **v3: flagged decisions are
  excluded from the primary denominators** and reported as a sensitivity. A search bound is not a band (S4), and it
  widens the lever spread that the oracle ceiling (sec. 5.0) rewards;
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

### 5.0 Step 0: the ratio error and the oracle ceiling (v3; rewritten in v3.1; runs first)

Step 0 has two stages, run in this order, each fixed here before it runs. It comes before the early variant (sec. 12)
and before any design is simulated. Its generator writes `results/validation/decision_relevance_step0.{json,md}`.

#### 5.0.1 Stage 0c: the ratio error sigma_r of the engine, between arms of one study

**What it is.** A noise-model estimate on spent and seen data. It credits nothing, fits no constant, scores no rival
and moves no shipped output. It sets only the discrepancy model of sec. 3.4.

**The units.** Groups g = (study, compound, lane) with k_g >= 2 arms, all of one laboratory (`lab_ids.yml`; a group
whose arms carry no common id is excluded). Sources:
1. `results/validation/contrast_scores.json`, its within-study pairs (each pair gives e_ab = engine_dex - measured_dex,
   the ratio residual). Within a group, the arm residuals r_i - mean(r) are recovered from the pairs by least squares
   on the group's pair graph (exact when the graph is consistent, which the pair definition guarantees). A group with a
   disconnected graph is split into its connected components.
2. `results/validation/core_panel_scores.json`: groups of answered rows of one study (benchmarks citing one source),
   same compound. These are mostly non-thiol: the thiol lane's absolute rows are refused there by H0.
3. `results/validation/kinetic_core_contrast_benchmark_score_v1.json`, its answered within-study cells, read as ratio
   residuals as in source 1.

An arm seen in two sources is counted once, preferring source 1.

**Inclusion.** The primary takes out-of-sample units only (`out_of_sample` true, or not a fit row of the shipped
calibration). In-sample units are a sensitivity: their residuals were shrunk by the fit, so they understate sigma_r, in
favour of the designs. Excluded: censored arms and arms whose engine value is missing.

**The estimator, per lane (thiol = sulfur lane; trunk; acrylamide; lipid).**
- **The pooled within-group variance:** sigma_w^2 = sum_g sum_i (r_gi - mean_g r)^2 / sum_g (k_g - 1).
- **The measurement part,** sigma_m^2 = the mean over arms of each arm's measurement variance (item 1's rule).
  - An arm's value is a mean of n replicates, so its variance is SD_rep^2 / n, in dex (log10(1 + CV)).
  - **A printed SD is read as one replicate's SD (SD_rep)**, unless the source says it is of a mean (an SE). An SE is
    multiplied by sqrt n to give SD_rep.
  - An unstated n is 2, as sec. 7 item 5 rules.
  - With no printed SD, the declared assay CV is taken per replicate: SIDA 10 %, HS-SPME / external standard 20 %.
  - The step-0 record states, per arm, which of these cases applied.
- **THE PRIMARY QUANTITY:** the ratio error sigma_r = sqrt(2 (sigma_w^2 - sigma_m^2)), the discrepancy SD of a
  within-study log ratio. omega does not enter it.
- **The floor (item 2).** sigma_r is never below sigma_floor = sqrt(2) x sigma_m, the ratio SD that measurement noise
  alone gives. The floor is sourced: it is the arms' own measurement model, never 0. When sigma_w^2 < 2 sigma_m^2 the
  floor is used, and the record says so.
- **Derived diagnostic only:** rho = 1 - sigma_r^2 / (2 omega^2). It is not identified separately from omega (item 1).
- **Robust sensitivity:** sigma_w from the median absolute deviation of the within-group residuals (x 1.4826).

**The minimum number of studies (item 3): G_min = 5 per lane.**
- With G studies, any interval over study-level estimates has G - 1 degrees of freedom. At G = 5, t_0.95,4 = 2.13;
  below it the interval is about twice the normal one or wider.
- Cluster-robust inference is known to over-reject badly with fewer than ~5 clusters, even with small-sample
  corrections (Cameron, Gelbach & Miller 2008; Cameron & Miller 2015).
- Five is also the number of scorable laboratories H0 itself rests on.
- **The author knew, before setting it,** that `contrast_scores.json` holds 5 studies in all (76 sulfur pairs, 5
  acrylamide), and states that here.
- **Below G_min,** the lane's sigma_r is **UNDECIDED** and is reported descriptively only.

**The interval (item 3: a small-G method).** A t-interval over study-level estimates.
- For each study s: its own within-group variance v_s, with df_s = sum over its groups of (k_g - 1).
- **The pooled mean:** v-bar = sum_s df_s v_s / sum_s df_s. Its standard error comes from the df-weighted spread of the
  v_s about v-bar, and the 90 % interval uses t with G - 1 df.
- The interval is mapped to sigma_r through sigma_r = sqrt(2 (v - sigma_m^2)), with the floor applied at the lower end.

**The leave-one-study-out check (item 4: symmetric).** sigma_r is recomputed dropping each study in turn, giving
LOSO_min and LOSO_max.
- The layer-1 verdict (sec. 7.1) must be the same at sigma_r = LOSO_min, at the point estimate and at LOSO_max.
  Otherwise it is **UNDECIDED**.
- **One choice, no bootstrap draws (item 4):** every later use takes the point estimate sigma_r as fixed; the LOSO
  bracket is the robustness check. v3's "draw rho from the bootstrap" is withdrawn.

**Which sigma_r a decision uses (item 7).**
- A decision uses its output's lane's sigma_r when that lane has >= G_min studies.
- **The pooled rule:** otherwise it uses the sigma_r pooled over every lane (all groups in one sigma_w^2, with every
  lane's studies counted toward G_min).
- If the pooled set also has fewer than G_min studies, the decision's eligibility is UNDECIDED.

**Cost.** It reads JSON already on disk, with no engine run: < 0.1 core-h of compute and ~2 - 3 h attended.

**Seen by the author before freezing it:** one pair of `contrast_scores.json` (Schibilsky 2019 HMF, pH 5 vs 8), and
that file's study and lane counts, printed while inspecting its structure. The plan's summary error statistics are in
sec. 1.

#### 5.0.2 Stage 0a: the oracle ceiling, per decision

**The oracle** is a perfect design: it reveals each truth's L_d(theta_i) exactly. What remains is the discrepancy of
the formulator's own lever:
- **The bench truth:** T_d(i) = L_d(theta_i) + delta_d(i), with delta_d ~ N(0, sigma_r^2) and sigma_r from stage 0c
  (the decision's lane, or the pooled rule);
- **The oracle's interval:** L_d +- 1.645 sigma_r;
- **Shown** when |L_d| > 1.645 sigma_r (the rival is known, so S3 agrees trivially);
- **Correct** when sign(L_d) = sign(T_d).

**The ceiling.** G*_d is the oracle's net gain of sec. 7: P_right - P_wrong|worth, given |T_d| >= 0.3 dex. It is
computed with 1000 delta draws per world, from a fixed seed, together with the oracle's P_wrong.

**Eligible decisions.** A decision is **eligible** when its oracle meets sec. 7's bar: G*_d >= 0.6 and P_wrong <= 0.1.
- **Only eligible decisions judge the designs** (secs. 7.1 - 7.2).
- A decision that is not eligible is **omega-bound**: no design can show it, under the current discrepancy.

**Reported per decision:**
- G*_d and the oracle's P_wrong, at sigma_r's point estimate and at LOSO_min and LOSO_max;
- **the reduction needed:** sigma*_d, the largest sigma_r at which the oracle meets the bar (bisection on [0.01, 3]
  dex), and the ratio sigma_r / sigma*_d;
- for the table, the sigma_r at which >= 2 decisions, and a majority of them, would become eligible.

**The step-0 ensemble (v3.1: declared exactly).**
- **512 worlds.** Each is one draw of the shipped envelope sampler (`uncertainty.sample_draws(512, seed = 20261107)`).
  That sampler draws the bands the S1 record cites as its own band sources, so it is the S1-band sampler.
- **One pH-balance rival per world:**
  - B1 and B13 alternate across the worlds;
  - within a rival, one of the passing points of PR #108's points record, chosen uniformly;
  - that point's PhDrift replaces the draw's own drift, and its Config changes (inline coupling, amine survival, and
    for B13 the acid supply) are applied.
  - The record is read at the pinned commit 5e61fed0 (`git show`). B0 (falsified) is never used.
- **Declared departures from sec. 4.2 at step 0:**
  - the 25 sulfur coordinates whose S1 band is "fixed" stay at the envelope's treatment (held);
  - 19 bands are normals centred on the shipped fit, so their centre is the incumbent's value; no other band exists on
    disk for them;
  - Q-sup, Q-loss, Q-split and Q-fft rivals beyond their bands are not switched in.
  - All of these narrow the prior spread of L_d. That tends to lower the ceilings, so it biases step 0 toward
    "omega-bound".
- **Hence the binding layer-1 verdict is the full run's recomputation**, on its own complete ensemble. Step 0's
  verdict is reported as "step 0 reading (envelope incomplete)", in the words of sec. 7.1.
- **Computable at step 0** are the decisions whose rivals run today:
  - D01, D02, D05, D06, D07, D08, D09, D10, D14, D15, D16, D18 (12 of 17);
  - O01 - O05 and O08 - O10 (8 of 10).
- **Not computable at step 0:**
  - D03 and D04 (the falsified hexose route; Q-hex not built);
  - D11, D12 and O06 (Q-buf: R0 falsified, the switch not built);
  - D13 and O07 (Q-fat: R0 falsified, no rivals).

  Each is listed with its reason.
- **The arms** are the 13 distinct arms of the computable D-decisions and the 10 of the computable O-decisions, 23 in
  all. They are run through `engine.predict` with the world's draw. Each arm's value is read from the run's internal
  molar state (`species_mmol_per_l`, through the declaration's `mapped_targets`), so only ratios are used and no
  absolute level is displayed (H0).
- **Execution:** nice -n 10, 2 processes, threads pinned to 1, a 300 s limit per arm. A failed arm leaves that world's
  decision missing; the failed fraction is reported.

**Cost.** 512 worlds x 23 arms x 0.42 - 1.5 s = **1.4 - 4.9 core-h**: 0.7 - 2.5 h wall at 2 processes. Attended:
~4 - 6 h for the generator and its toy-input unit test.

### 5.1 Generating the data

- **Ensemble.** N = 2048 worlds (sec. 4).
- **Truths.** Every world serves as a truth once. Its data are scored against the other N - 1 worlds (leave-one-out),
  so truths are drawn outside the scoring ensemble, as the judge requires.
- **Per truth i and design g:** y_g(i) = the world's prediction on g's cells + the discrepancy delta per cell (u + e_i,
  sec. 3.4) + run noise per replicate, censored at the floors.
- **The two targets (v3, declared).** The generator's reality is the model plus delta, so two truths exist:
  - **L_d(theta_i)**, the model's lever: the target of every **model-part** quantity (shrinkage, the model-part
    interval, A2's model-part coverage);
  - **T_d(i) = L_d(theta_i) + delta_d(i)**, the formulator's lever at the bench: the target of every **total**
    quantity (P_right, P_wrong, the total interval's coverage). delta_d ~ N(0, sigma_r^2) as in sec. 5.0.2, independent of
    every design's delta. Its arms are another laboratory's pots; that choice is conservative for designs that
    re-run the decision's own recipe.
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

### 5.2 The update: two co-primary estimators (v2)

**Estimator 1: the Bayes-linear (ensemble-Kalman) preposterior.** It needs no weights, so it cannot collapse.
- **Inputs.** For each world j, simulate y_g(j) under sec. 3.4's measurement model. Levels are taken in log10, with
  censored values set to their floor.
- **The estimate.** The best linear predictor of L_d from y_g:
  - m_BL(y) = E L_d + Cov(L_d, y_g) Cov(y_g)^-1 (y - E y_g);
  - Var_BL(L_d | y_g) = Var(L_d) - Cov(L_d, y_g) Cov(y_g)^-1 Cov(y_g, L_d).
- **Covariances.** Ledoit-Wolf shrinkage of Cov(y_g), with its analytic intensity and no tuning (v3: ridge-CV dropped,
  so there is one choice and no analyst freedom). Cell means of the replicates are used, so p_g <= ~500 for G-REC.
- **Cross-fitting.** 5-fold: the coefficients are fitted on four folds and applied to the fifth. Each truth's m_BL and
  the reported Var_BL are therefore out-of-fold.
- **Shrinkage:** s_BL(d, g) = 1 - Var_BL(L_d | y_g) / Var(L_d), with Var_BL the out-of-fold mean squared error.
  - A linear predictor's mean squared error is at least the Bayes posterior's expected variance. So **s_BL is a
    conservative bound on the Bayes shrinkage**, up to Monte Carlo error.
  - Width ratios use sqrt(1 - s_BL).
- **Its decision rule.** The formulator reads N(m_BL(y), Var_BL + sigma_r^2), the total interval for T_d. The rule is scored empirically over
  the truths (sec. 6), so its P(wrong) is a frequency of that rule, not a model claim.
- **Limits.** It approximates the mixture posterior (I_post), not the S3 union (I_ship). Credit (sec. 5.4) and the S3
  union are PSIS-only.

**Estimator 2: importance reweighting (PSIS), kept from v1.** Each world j != i gets weight w_j(i, g) ∝ p(y_g(i) |
theta_j), under sec. 3.4's measurement model. The weights are Pareto-smoothed through the judge's own estimator
(`importance_lppd`, `psis`), with its acceptance rule (k-hat < 0.7, ESS >= 100). It is the only estimator that carries
the S3 union and M1.1's T.

**Why not Laplace (reasons 1 - 3 decide; the cost is corrected in v2):**
1. **The rivals are discrete.** A world is a structure index plus coordinates. Laplace approximates one mode of one
   structure.
2. **The bands are bounded.** Uniform bands with mass at their edges, and a posterior pressed against a band wall (S4's
   band-edge criterion), are where a Gaussian at the mode fails.
3. **Lever outputs flip sign.** A Gaussian on the coordinates, pushed through a sign-flipping output, misstates P(sign).
4. **Cost (v2, corrected).** The timing run gives 49 s for a 5-coordinate, 1-pot fit (sec. 9.1). Scaled to a design's
   ~20 pots and ~30 relevant coordinates, that is ~0.9 h per update.
   - About 100 truths per design are enough for the metrics, so the total is ~100 x 15 designs x ~4 rivals x 0.9 h
     ~= 5 x 10^3 core-h. v1's 10^5 assumed every one of the 2048 truths.
   - Still ~30 times the whole ensemble; but the rejection rests on reasons 1 - 3.

**No enrichment (v2).** v1's enrichment mixture was fitted to the pooled posteriors of all truths, so it sat close to
the prior and could not rescue a collapsed truth. It is dropped. Estimator 1 covers the collapse instead.

**Unresolved truths are UNDECIDED (v2).** A truth whose PSIS weights fail the acceptance rule is unresolved for
estimator 2.
- It is never read as estimated, and never counted as relevant or as not relevant.
- v1's rule (unresolved counts as shrinkage >= 25 %) let the fallback refute the worry. It is removed.

**How the estimators combine, per (d, g) cell (v3: BL alone can never make a cell relevant).**
- **Relevant** only when PSIS is validated (sec. 5.5), resolved on >= 80 % of the truths, and the cell meets sec. 7's
  bar with S3's every-rival sign agreement.
- **Not relevant** when:
  - PSIS is validated and resolved (>= 80 %) and the cell misses the bar; or
  - PSIS is not validated or not resolved, and the validated BL rule misses the bar. Such a cell is flagged "BL-only".
- **UNDECIDED** otherwise. That includes a cell where BL meets the bar but PSIS cannot confirm it under S3.
- A verdict of sec. 7 that UNDECIDED cells could flip is itself UNDECIDED.
- **Sensitivity:** BL-only cells re-read as UNDECIDED. BL is a conservative bound, so it can miss a bar the Bayes rule
  meets. This sensitivity shows how much the "not relevant" reading leans on it.

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
  - the **total**, which convolves each world's lever with the ratio discrepancy N(0, sigma_r^2), sigma_r from
    stage 0c (sec. 5.0.1): what the formulator reads, and what covers T_d. "As
    shipped today" (iid, 0.59 dex, about +-0.98 dex at 90 %) is a sensitivity.
- **"Shown" (v2).** A direction is shown when:
  - the total 90 % interval excludes zero;
  - and, for I_ship, every surviving rival agrees in sign (S3).

  The BL rule (sec. 5.2) uses its own total interval, N(m_BL, Var_BL + sigma_r^2).

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
| A1 per update | k-hat and ESS of every (truth, design) PSIS weight set | k-hat < 0.7 and ESS >= 100; the unresolved fraction is reported per design |
| A2 calibration (simulation-based) | over the truths: the coverage of I_post's **model part for L_d**, and of its **total for T_d** (v3: each interval against its own target); and the uniformity of the truth's rank in its posterior (10 bins, chi-square, Holm across d x g cells) | coverage in [0.85, 0.95] for both, and Holm p >= 0.01, in >= 80 % of cells |
| A3 Laplace cross-check (the JAX fit path) | G-C and G-M4; 10 truths each; one rival; a sub-problem of <= 5 coordinates and <= 4 pots, the rest held at the truth. `fit.fit` with the exact Jacobian, then `fit.laplace`, against reweighting a dedicated 1024-draw sub-ensemble on the same sub-problem | posterior mean within 0.25 posterior sd and sd ratio within [0.8, 1.25] in >= 8 of 10 truths per design. A Laplace fit pressed against a band wall is reported, not counted |
| A4 exact quadrature | a 2-coordinate sub-problem (`quadrature_lppd`) against reweighting | lppd within 0.05 nats |
| A5 mini-model parity | sec. 5.1 | 0.1 dex |
| A6 Bayes-linear (v2; targets fixed in v3) | the out-of-fold MSE against the in-fold MSE (overfitting); the coverage of N(m_BL, Var_BL) for **L_d** and of N(m_BL, Var_BL + sigma_r^2) for **T_d** (never the model part against T_d, nor the total against L_d: v2's window failed exactly when a design was informative); s_BL against PSIS's variance shrinkage on the cells where PSIS is resolved (s_BL must not exceed it beyond Monte Carlo error, since BL is the conservative bound) | out-of-fold / in-fold MSE <= 1.25 and coverage in [0.85, 0.95], in >= 80 % of cells; the bound holds in >= 90 % of resolved cells |

**Validation per estimator (v2).**
- **PSIS is validated** if A1 (resolved >= 80 % of truths per design), A2, A3 and A4 pass.
- **BL is validated** if A6 passes. A2's coverage test also applies to BL's own interval.
- If **neither** is validated, the run reports "update not validated" and no verdict of sec. 7.
- If **only one** is validated, that estimator alone decides, and the run says so.

### 5.6 Run order

0. **Freeze.** This file's commit is recorded. Then **step 0**, in this order:
   - stage 0c (sec. 5.0.1): sigma_r per lane;
   - stage 0a (sec. 5.0.2): the oracle ceilings and the eligible set;
   - the refusal state of every decision in the prior ensemble, under the rules in force that day. PR #108's rule
     is **pinned to commit 5e61fed0** (its head on 2026-11-07, open). If #108's refusal set moves before the run
     (a new commit, or its merge with changes), step 0's refusal stage reruns at the new commit. **The new
     commit's refusal set is the primary**, as it is the rule the engine will apply; 5e61fed0's is reported beside
     it.
   Then the early variant (sec. 12) may run; then the full run's steps 1 - 8.
1. **Builds and parity.** List what is runnable (sec. 4.1); the parity checks (A5).
2. **Pilot.** 32 worlds, timed. The cost is re-projected (sec. 9.3).
3. The ensemble.
4. The synthetic data.
5. The updates (BL and PSIS) and the leave-one-rival-out run (sec. 5.7); A1.
6. Checks A2 - A4 and A6.
7. Metrics and verdicts.
8. Write-up, independent adversarial review (sec. 7 item 8), EX entry.

No metric is computed before step 6 passes.

### 5.7 The open world: leave one rival out (v2)

Every verdict above assumes the truth is among the carried rivals, inside their bands. That is the closed-world
assumption, and it is stated in sec. 1.

**The test.** For each question q and each of its rivals r:
- draw the truths from the worlds with rival r;
- score them, under both estimators, against the ensemble with every r-world removed.

This is the case where the true structure is one nobody carried. It reuses the same forward ensemble at no extra
forward cost.

**Reported per (d, g):** P_wrong and P_right (against T_d, sec. 6 item 1) under leave-one-rival-out, as the
worst rival r and the mean over r.

**What it can change.** A cell that is relevant under sec. 7's bar but has P(shown and wrong) > 0.1 under
leave-one-rival-out is marked **"relevant only in a closed world"**.
- The worry verdict is reported both with and without those cells.
- The primary verdict counts them as relevant. The marking is reported beside it.

## 6. Metrics, per decision and design

All are means over the truths, with a 90 % bootstrap interval over truths (1000 resamples). "Prior" is the same
quantity under the ensemble with no data. Each is computed under both estimators where it is defined.

1. **The answer a formulator acts on (the primary metric), on the TOTAL interval, against T_d (v3).**
   - **The bench lever** is T_d(i) (sec. 5.1). It is **worth acting on** when |T_d| >= 0.3 dex (a factor of 2).
   - **P_right(d, g)** = P(shown, and the shown sign is correct | worth acting on).
   - **P_wrong|worth(d, g)** = P(shown, and the shown sign is wrong | worth acting on).
   - **P_wrong(d, g)** = P(shown and wrong) over every truth with |T_d| > 0.01 dex (PR #108's sign resolution).
   - **The net gain** G(d, g) = P_right - P_wrong|worth (sec. 7 derives the bar from it).
   - P(refused) is reported beside them.
   - "Shown" is defined in sec. 5.3: with I_ship, S3's rule plus a total interval that excludes zero; with BL, the BL
     total interval excludes zero.
   - **The oracle ceiling** G*_d and the eligibility of each decision come from step 0 (sec. 5.0.2).
2. **Expected shrinkage of the model part (secondary in v2).**
   - PSIS: s(d, g) = E_i[1 - width(I_ship,model | y_g(i)) / width(I_ship,model)], with the median, the I_post version
     and the total-interval version.
   - BL: s_BL(d, g) (sec. 5.2).
   - The omega-dominance flag: the model-part width is below the omega-part width.
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
   - (c) **V_real(g) = the mean s_BL(d, g) over the 17 computable decisions** (v3: one choice). It is continuous and
     defined whether or not any decision is eligible, so the Spearman of sec. 7.3 cannot be undefined by ties at zero.
     The count of decisions meeting the bar is reported beside it, not used.
   - (d) The Spearman correlation between V_rank and V_real across designs (sec. 7.3).
5. **The effective number of independent decisions (v2).** Fifteen of the 17 share BR-S and Q-ph, so they are not
   17 independent tests. Two estimates, both from the prior ensemble:
   - n_eff,L = (sum lambda)^2 / sum lambda^2, the participation ratio of the eigenvalues lambda of the correlation
     matrix of the 17 L_d;
   - n_eff,rel: the same ratio on the correlation, over truths and designs, of the per-decision "shown and correct"
     indicators.

   Both are reported beside every count. The reviewer's estimate is 3 - 5; the author's ~4.
6. **The per-question count (v2).** Each decision is assigned the structural question carrying the largest share of
   its prior variance (a first-order variance decomposition over the rival index and the coordinates of each question).
   Every "k of 17" is also given as k per question, so that 9 decisions all resting on Q-ph read as one question.
7. **Information per budget unit (v2).** For each design:
   - sum_d s_BL(d, g) per heating;
   - the same per (cell x readout): each readout of each cell counts once;
   - the number of eligible decisions meeting the bar, per 100 heatings.

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

**Already true, before any run (v2 wording): `rank` cannot value designs at all.**
- **What `rank` does.** It ranks benchmark ROWS, not designs. Each row (benchmark, compound) is scored
  (miss + 0.3 x width) x DR on the absolute 5 - 95 % quantiles of the envelope (`_voi_score`, `src/experiment_value.py`
  lines 273 - 285). Here miss is the log distance of the measured value outside the interval (zero inside), width the
  interval's width in dex, and DR the threshold-anchored relevance.
- **What it lacks.** It has no notion of a design, of a lever, or of what a measurement would change. V_rank above is
  this file's construction on top of it, so Q3 tests that construction, not a design ranking `rank` makes.
- **The design sets (v3: one definition, used everywhere).**
  - **The on-record singles (11):** G-M4, G-M4+, G-M5B, G-M5F, G-M5W, G-M3, G-M1c, G-SPL, G-PHL, G-LOSS, G-FUR.
  - **The unions:** G-M5, G-REC, G-C+REC.
  - **The `rank` set (12):** the 11 on-record singles plus G-C.
- **The coverage.** Even so, under the strict rule the construction gives zero value to 7 of the 12 designs in the
  `rank` set (G-M4, G-M4+, G-M3, G-M1c, G-SPL, G-LOSS, G-FUR). Only one of `rank`'s 24 rows (2-acetylthiazole) is
  a meaty thiol-lane row.
- **The minimal fix** is a follow-up (sec. 11). It is not implemented here.

## 7. Pre-registered thresholds and verdicts

**The denominator.**
- The 17 decisions whose lever interval is computable: D01 - D16 and D18. D17 and D19 are structural.
- **v3: the designs are judged only on the ELIGIBLE decisions** of step 0 (sec. 5.0.2), E <= 17, with the flagged
  decisions of sec. 4.2 excluded.
- Every count is reported with n_eff and the per-question count (sec. 6 items 5 - 6).

**The bar (v3: re-derived; the threshold is a chosen value, and says why).**
- **The formulator's choice.** A formulator either runs a lever's first bench test in the direction the tool shows, or,
  when the tool shows none, runs it blind and is right half the time (v2 treated "not shown" as a failure; it is a
  coin).
- **The success rate.** For a lever worth acting on, the first test goes the right way with probability
  S = P_right + 0.5 x P(not shown | worth) = 0.5 + 0.5 x (P_right - P_wrong|worth) = 0.5 + 0.5 x G.
- **The bar.**
  - **G >= 0.6, that is S >= 0.8.** The tool must remove at least 60 % of the failed first tests a blind formulator
    makes (0.5 -> 0.2). **0.8 is a chosen value**, not a measured one: it is the smallest gain judged worth building a
    lever table for. Read on P_right alone, it is v2's 0.7 when P_wrong|worth is at its cap of 0.1.
  - **P_wrong <= 0.1, over all truths.** At most one shown direction in ten may send a formulator the wrong way. A
    calibrated 90 % interval that excludes zero errs at most 5 % one-sided; 0.1 allows for miscalibration and no
    more.
  - **|T_d| >= 0.3 dex** defines "worth acting on": a factor of 2.
- **A design is decision-relevant for an eligible d** when G and P_wrong meet the bar after the update, and G rises by
  >= 0.1 over the prior.
- **"A majority":** >= ceil(E / 2) of the eligible decisions (9 of 17 when all are eligible). Every verdict is worded
  "on E of 17" (v3.1).
- **Secondary:** model-part shrinkage s >= 0.25 (a variance ratio of 0.56, ~1.8 times the information held) on either
  estimator, over all 17 decisions.

### 7.1 The worry (Q1)

The verdict has two layers, read in order.

**Layer 1, the ceiling (step 0, before any design is simulated).** It is computed at sigma_r's point estimate, LOSO_min
and LOSO_max, and must agree at all three (sec. 5.0.1). Otherwise it is UNDECIDED. If any decision's sigma_r is
UNDECIDED (fewer than G_min studies), layer 1 is UNDECIDED for the decisions it governs. At step 0 it is the "step 0
reading (envelope incomplete)", on the computable decisions; the full run's recomputation binds. If **fewer than 2
decisions are eligible**, the verdict is, verbatim:

> **"CONFIRMED (omega-bound): the lever table cannot show directions under the current discrepancy, whatever the
> experiments."**

It is reported with:
- the reduction needed: sigma*, the ratio discrepancy at which >= 2 decisions, and a majority, would become eligible;
- the factor by which the measured sigma_r would have to fall to reach sigma*.

In that case no design verdict is given. Q2 is "not informative (omega-bound)". The design metrics are reported
descriptively, and Q3 still runs on s_BL.

**Layer 2, the designs (only if E >= 2; worded "on E of 17").** The primary is sec. 7's bar on the eligible decisions, under the
combination rule of sec. 5.2.
- **WORRY CONFIRMED (designs)** if every on-record single (the 11 of sec. 6.1) and both on-record unions (G-M5,
  G-REC) are decision-relevant for **fewer than ceil(E / 2)** eligible decisions.
- **WORRY REFUTED** if any of them reaches ceil(E / 2), **with E >= 3**. One decision must never refute: at E = 2 a
  "majority" would be one decision. **At E = 2, REFUTED is not available** and the outcome is CONFIRMED (designs) or
  UNDECIDED. **A REFUTED verdict must also hold under the laboratory-offset sensitivity** (sec. 3.4, tau = 0.964 on
  design absolute levels). If it does not, the verdict is UNDECIDED.
- **UNDECIDED** otherwise. That covers two cases:
  - UNDECIDED cells could flip the verdict (sec. 5.2);
  - the borderline rule fails: a cell is borderline when the bootstrap interval of G spans 0.6 or that of P_wrong
    spans 0.1, and a verdict stands only if it survives flipping every borderline cell against itself.
- **The omega-bound share** (17 - E decisions) is reported with every layer-2 verdict, as the part of the lever table
  no design can reach.

**Reported beside both layers, not deciding:**
- the per-question count and n_eff;
- the closed-world-only cells (sec. 5.7);
- the BL-only sensitivity (sec. 5.2);
- the secondary (model-part shrinkage) verdict on both estimators;
- the composite: P(the refusal lifts) >= 0.5 for a data-liftable refusal also counts;
- the off-grid set;
- the complete-envelope decisions;
- the flagged decisions of sec. 4.2;
- the misspecified generator;
- omega 0.390;
- sigma_r in {0.593, 0.42, 0.19}, and at the LOSO bracket;
- "as shipped today" (iid).

**Harmful.** A design for which P_wrong rises above 0.1 on any decision is listed as harmful there, whatever else it
does.

### 7.2 Targeting (Q2), judged on the eligible off-grid decisions

- **The primary set** is the eligible members of O01 - O10 (sec. 2.4), E_O of them.
- If **E_O < 2**, Q2 is **"not informative (omega-bound)"**. v2 would have scored 0 - 0 as a tie.
- **LIST MIS-TARGETED** if G-C is decision-relevant for >= m more eligible off-grid decisions than the best on-record
  single with a budget <= 96 heatings (G-M4, G-M5B, G-M5F, G-M5W, G-M3, G-SPL, G-PHL, G-LOSS, G-FUR), or the union
  G-M5, and the bootstrap 90 % interval of that difference excludes 0.
  - m = max(1, ceil(0.3 x E_O)).
  - The 0.3 is **the author's choice** (v2's "3 of 10", rescaled; not asked for by the premortem; see the revision
    log).
- **Not mis-targeted** if the difference is <= 0 with its interval excluding +m. Otherwise inconclusive.
- **Reported:**
  - the interpolating (O01 - O07) and extrapolating (O08 - O10) subsets apart;
  - the secondary on-grid set (D01 - D16, D18), flagged "on-grid: partly G-C's own cells";
  - G-C against G-REC (394 heatings);
  - G-C+REC against G-C;
  - the information per heating and per cell x readout (sec. 6 item 7).

### 7.3 `rank`'s first validation (Q3)

- **Coverage.** `rank` scores rows, not designs (sec. 6.1). If V_rank(g) = 0 for at least half of the `rank` set, the
  verdict is **`rank` CANNOT VALUE DESIGNS**. By the mapping this is expected; the run confirms it under the system rule
  too.
- **Ordering.** The Spearman rho between V_rank and V_real = mean s_BL (sec. 6 item 4c) over the 12 designs of the
  `rank` set, ties at mid-rank. It does not depend on eligibility.
  - **VALIDATED** if rho >= the one-sided 5 % permutation critical value (~0.50 at n = 12).
  - **CONTRADICTED** if rho <= 0 and the design with the largest V_real has V_rank = 0 or the lowest V_rank.
  - Otherwise **INCONCLUSIVE**. That includes the case where V_rank or V_real is constant, so that the correlation is
    undefined.
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
- w_j are the PSIS posterior weights of sec. 5.2. For the Bayes-linear estimator, and for the prior weight, w_j is the
  ensemble's equal weight.

**What changes when thresholds arrive.**
- The weighted fraction of decision-relevant decisions, sum_d w_d 1[d meets sec. 7's bar under g] / sum_d w_d,
  re-runs sec. 7.1's test (>= 0.5 refutes) and sec. 7.2's (on weighted counts).
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

### 9.2 The projection (v2)

Per world:
- **148 network integrations:**
  - 28 distinct decision arms: the 16 of sec. 2.2 (D17 and D19 are not run) and 12 off-grid arms (sec. 2.4);
  - 11 `rank` benchmark pots (24 rows);
  - 109 network design cells: G-M1c 50, G-SPL 7, G-LOSS 5, G-FUR 3, G-PHL 12, G-C 32;
- **plus ~60 mini-model cells** (G-M4 32, G-M5 19, G-M3 6) at <= 0.1 s.

Per-integration cost: central 1.5 s, low 0.42 s (inline NumPy), high 2.9 s (the B0-like engine cost; the JAX predictor
measured 2.7 s). v2 drops the enrichment ensemble. The Bayes-linear estimator, the leave-one-rival-out run and the
omega-correlation forms are array work on the same forward runs.

| item | central | low | high |
|---|---|---|---|
| step 0: stage 0a, 512 worlds x 23 arms (stage 0c: < 0.1) | 4.9 | 1.4 | 9.5 |
| pilot, 32 worlds | 2.0 | 0.6 | 3.9 |
| ensemble, 2048 worlds | 129.7 | 38.8 | 247.6 |
| array work: BL (5-fold, Ledoit-Wolf), PSIS, nulls, leave-one-rival-out, the sigma_r and tau sensitivities, metrics | 4.0 | 4.0 | 4.0 |
| A3 Laplace check: 20 updates of ~4 pots (~3.3 min each) plus a 2 x 1024-draw sub-ensemble | 4.5 | 2.1 | 7.7 |
| A4 quadrature | 0.2 | 0.2 | 0.2 |
| contingency (10 %) | 14.6 | 4.7 | 27.5 |
| **total, core-hours** | **161** | **52** | **302** |
| **wall-clock at 2 processes** (this task's limit) | **80 h** | **26 h** | **151 h** |
| **wall-clock at 8 cores** (the owner's standing cap) | **20 h** | **7 h** | **38 h** |

**Attended work: 34 - 48 h** (step 0c adds 2 - 3 h; stage 0a's sampler is counted in the harness).
- The harness: world sampler from the S1 record and the case bands, rival switch registry, generators per design,
  metrics, `rank` mapping: 14 - 20 h.
- The two estimators with checks A1 - A6, and the leave-one-rival-out run: 5 - 7 h.
- Parity checks: 3 - 4 h.
- Q-fat's reduced operator: 4 - 6 h.
- Pilot: 1 - 2 h.
- Write-up: 3 - 4 h.
- Answers to the review: 2 - 3 h.

**Not charged here.** Prerequisite builds already on the plan:
- A2b with M4's structure at parity (sec. 6, 11-21 to 12-05);
- M5's buffer switch;
- M1c's rival switches;
- the OPD trap module for G-SPL's in-situ arm. If it is not built, that arm is dropped and G-SPL is marked reduced.

The earliest start of the full run is therefore the week of 12-08. The early variant (sec. 12) needs none of these
builds.

### 9.3 The pre-declared reductions

- If the pilot projects more than 48 h of wall-clock at the cores granted, N falls to 1024. The central total falls to
  ~90 core-h: ~11 h at 8 cores, ~45 h at 2.
- Monte Carlo error then grows by ~sqrt 2. It is reported, and the borderline rule of sec. 7.1 absorbs it.
- Paid Colab (CPU, many cores) is an option to propose to the owner, not to assume: draws only, scored locally. JAX batch
  draws gave no CPU gain at equal tolerance (PR #74), so a GPU is out until a CPU-vs-GPU benchmark says otherwise.

## 10. Before the freeze

1. **An independent review** whose job is to break this design (sec. 7 item 8). Done for v1 (premortem: "freeze with
   edits"); v2 applied edits 1 - 10. The confirmation of v2 returned "not confirmed", and v3 answers its items 1 - 7.
   **A confirmation read of v3** checks that each item is applied as asked. In particular it checks:
   - stage 0c's estimator, unit selection and leave-one-study-out rule;
   - the oracle ceiling and the eligibility rule;
   - the two targets (L_d and T_d) in A2 and A6;
   - the BL-cannot-credit rule;
   - the rebuilt off-grid set.
2. **Answers to every objection** in this file's revision log.
3. **The freeze commit recorded here.** Then sec. 5.6's step 0.

## 11. Follow-ups (not done here)

1. **`rank` as a design ranker: the minimal fix** (from sec. 6.1; a code change with its own PR, not implemented here):
   - **Lever rows.** Add rows for levers, not only for pots. Their width is the paired-draw ratio width: both arms from
     the same draw, so shared errors cancel as in `compare`. Add a sign-disagreement term: the share of draws whose sign
     disagrees with the median's.
   - **Design value.** V(g) = sum_d DR_d x Cov(L_d, y_g) (Cov y_g + R)^-1 Cov(y_g, L_d) / Var L_d. This is sec. 5.2's
     Bayes-linear shrinkage weighted by DR, with R the design's measurement covariance (sec. 3.4).
   - **Per-draw outputs.** The envelope must keep per-draw outputs, not only quantiles, so that the covariances can be
     formed. `core_prediction_uncertainty.json` stores quantiles only today.
   - Validated against this file's V_real once the run exists.
2. **The scope list of sec. 2.2** (oil type and level, heme, yeast extract, protein amines, cooking method, sulfur a_w,
   Strecker) feeds W5's lever-table scope note.

## 12. The early variant (PARTIAL; pre-registered, v2; v3: after step 0)

**What it is.** A subset that runs before A2b, on rivals that already run in the NumPy engine. **It runs after step 0**
(sec. 5.0), with step 0's sigma_r and eligible set. It **cannot decide sec.
7.1**: it lacks every on-record design but three, and eight of the 17 decisions. A2b is needed only for Q-hex (D03 /
D04); the rest of what it lacks waits on M5's switch, M1c's build and the reduced operator.

- **Decisions:** D01, D02, D05, D06, D07, D08, D09, D15, D16 (9; 8 distinct arms).
- **Rivals:**
  - Q-ph: B1 and B13, inline NumPy coupling. **pH is NOT clamped**: a clamped pH removes Q-ph, the largest exposure,
    and would bias the result toward refuting the worry;
  - Q-sup: R0 and C1 only (envelope incomplete: X, X0, R-rel, R-cal, R-sink, R-K are not built);
  - Q-loss: C2 (with the banded frozen-form sink);
  - Q-split: B61's NF hydroxide decay against pH-flat;
  - Q-fft: its two bands;
  - every other coordinate from sec. 4.2.
- **Designs:**
  - G-C without its two phosphate cells (30 cells, 90 heatings);
  - G-LOSS;
  - G-FUR;
  - G-SPL with its ex-situ OPD arm only, plus B3 (4 cells).
- **Exclusions inside G-C: the author's own change (v2), not asked for by the premortem; see the revision log.** Two kinds of readout would
  be generated by a structure already falsified, so they are dropped from the synthetic data:
  - the MFT readouts of glucose cells (the hexose MFT route, Q-hex);
  - the thiol readouts of the ROOH cells (Q-fat, whose R0 is falsified and whose rivals are not built).

  Their other readouts are kept.
- **Analysis.**
  - 1024 worlds, about 53 integrations each (the 8 arms and ~45 design cells; no `rank` pots, since Q3 needs every
    design).
  - Bayes-linear and PSIS as in sec. 5.2.
  - Checks A1, A2, A4 and A6. A3 needs JAX and is not run.
  - The leave-one-rival-out run of sec. 5.7.
  - The metrics of sec. 6 items 1, 2, 5, 6, 7.
- **Cost.** 1024 x 53 x 0.42 - 1.5 s = **6 - 23 core-h**, plus ~0.5 h of array work: **~1 - 3 h wall at 8 cores**,
  ~3 - 12 h at 2 processes. Attended: ~8 - 12 h, most of it the harness the full run reuses.
- **What it reports, and the only reading allowed.**
  - Per-decision G, P_right, P_wrong, eligibility and s_BL for its four designs.
  - n_eff on its nine decisions.
  - Whether the estimators pass their checks on real ensemble shapes. That makes it the full run's pilot for A1, A2
    and A6.
  - **An indicator, never a verdict (v3: numeric).** "omega-bound" means not eligible at step 0 (oracle G* < 0.6 or
    oracle P_wrong > 0.1, sec. 5.0.2).
    - If **>= 7 of its 9 decisions are omega-bound** and G-C meets the bar for none of the eligible rest, it records:
      "the full run's Q1 is likely settled by omega".
    - If **fewer than 2 are eligible**, it records step 0's layer-1 sentence (sec. 7.1) as an early reading.
    - Whether to run the full simulation stays the owner's call.
- **Its results cannot amend this file.** The full run's rules stay as frozen.

## Revision log

- 2026-11-07, v1 DRAFT: written with the two timing runs and the declaration-only check. No simulation.
- 2026-11-07, **v2 DRAFT**: the independent premortem returned "freeze with edits". It estimated P ~0.55 that v1 would
  end "update not validated" or be decided by its fallback rules. All ten edits are applied. No simulation; the timings
  are unchanged.
  1. **Estimator** (sec. 5.2):
     - Bayes-linear (ensemble-Kalman) preposterior added as co-primary: cross-fitted, Ledoit-Wolf, a conservative
       bound; check A6;
     - the enrichment mixture is dropped (fitted to pooled posteriors, it sat near the prior);
     - unresolved truths are UNDECIDED, never relevant; v1's rule that let the fallback refute the worry is removed;
     - the combination rule is stated per cell.
  2. **Primary metric** (secs. 5.3, 6, 7): on the TOTAL interval, P_right (shown and correct given |L| >= 0.3 dex)
     >= 0.7 with P_wrong <= 0.1. Model-part shrinkage is secondary. An omega-bound flag is added, and an omega
     correlation sensitivity, rho in {0.5, 0.9} between neighbouring cells (sec. 3.4).
  3. **No teaching to the test:**
     - an off-grid decision set O01 - O10 (sec. 2.4), primary for Q2, with the margin restated as 3 of 10;
     - the effective number of independent decisions and a per-question count reported beside every "k of 17" (sec. 6
       items 5 - 6).
  4. **Open world:** the closed-world assumption is stated (sec. 1); a leave-one-rival-out run (sec. 5.7) reports
     P_wrong, and cells relevant only in a closed world are marked.
  5. **Thresholds:** the bar is derived from the primary metric, as what a formulator acting on the shown direction
     gets right or wrong (sec. 7).
  6. **Proxies** relabelled the "ROOH lever" (1 mM at 10 % oil ~20 meq/kg: a rancid-fat lever) and the "GSH lever" (no
     IMP / ribose route). The base's scope is stated (no extrusion, searing or heme), and the levers the engine cannot
     represent are listed as a scope finding (sec. 2.2).
  7. **Laplace cost** corrected to ~5 x 10^3 core-h (about 100 truths per design). The rejection rests on reasons
     1 - 3.
  8. **`rank`:** reworded as "cannot value designs at all", since it ranks rows by (miss + 0.3 x width) x DR on
     absolute quantiles. Its minimal fix is a follow-up (sec. 11), not implemented.
  9. **Budget:** information per heating and per cell x readout (sec. 6 item 7).
  10. **The early variant** (sec. 12), PARTIAL: 9 decisions, NumPy rivals, pH not clamped, 6 - 23 core-h.

  **The author's own changes in v2, not asked for by the premortem** (marked in v3, at the confirmation's request):
  - **(i) the early variant's dropped readouts**: G-C's glucose-cell MFT and ROOH-cell thiol readouts. They would be
    generated by falsified structures;
  - **(ii) Q2's margin of 3 of 10**: v1's 5-of-17 rescaled by the author. v3 rescales it again, to
    max(1, ceil(0.3 x E_O)).

  Also in v2: the cost is re-projected for the off-grid arms and the dropped enrichment (152 core-h central; sec. 9.2),
  and the forecasts of sec. 1 are revised.
- 2026-11-07, **v3 DRAFT**: the independent confirmation of v2 (d299af43) returned "NOT CONFIRMED". All ten edits were
  present, but two new defects decided the result by construction. Every item is applied. No simulation; the timings
  are unchanged.
  1. **The omega ceiling.**
     - **(a)** Step 0's oracle ceiling G*_d per decision (sec. 5.0.2), and a separate layer-1 verdict, "CONFIRMED
       (omega-bound)" (sec. 7.1).
     - **(b)** The designs are judged only on eligible decisions. If fewer than 2 qualify, layer 1 decides and Q2 is
       "not informative".
     - **(c)** NEW stage 0c: rho of the engine's discrepancy between arms of one study, from residuals on disk
       (sec. 5.0.1). It has its own estimator, a cluster-bootstrap interval and a leave-one-study-out rule, fixed before
       it runs. It is a noise-model estimate on spent data, not a credit. The measured rho is the primary noise model
       (sec. 3.4); rho = 0, 0.5 and 0.9 are sensitivities.
     - **(d)** The verbatim omega-bound sentence and the reduction needed are pre-registered.
  2. **Target:** declared two targets, L_d for model-part quantities and T_d = L_d + delta_d for total quantities
     (sec. 5.1). A2 and A6 test each interval against its own target.
  3. **Estimator rule:** BL alone can no longer make a cell relevant, only "not relevant" (flagged BL-only) or
     UNDECIDED (sec. 5.2).
  4. **The bar re-derived:** "not shown" is a 50 % blind test, so S = 0.5 + 0.5 x G. The bar is G >= 0.6 (S >= 0.8, a
     chosen value with its reason) and P_wrong <= 0.1 (sec. 7).
  5. **The author's own changes** in v2 are marked in its entry above: the dropped readouts and the 3-of-10 margin.
  6. **The early variant:** a numeric omega-bound threshold (not eligible at step 0, >= 7 of 9), and step 0 runs first
     (secs. 5.6, 12).
  7. **Analyst freedom removed:**
     - PR #108 is pinned to 5e61fed0, with a rerun if its refusal set moves;
     - V_real = mean s_BL only;
     - one set of designs is defined and used everywhere (sec. 6.1);
     - flagged simulation-only decisions are excluded from the primary;
     - rival weights are equal, fixed;
     - the off-grid set is rebuilt on BR-X, a base that is no design's cell, and split into interpolating and
       extrapolating decisions (sec. 2.4);
     - Ledoit-Wolf only.

  Also in v3:
  - the cost is re-projected for 12 off-grid arms and step 0 (161 core-h central; sec. 9.2);
  - the forecasts of sec. 1 are revised and marked as informed by seen summaries.
- 2026-11-07, **v3.1, FROZEN**: the re-confirmation of v3 (ccd02490) returned "confirmed with edits". Every edit is
  applied, and the file is frozen in the commit that set the header (its sha is in PR #110's body).
  1. **The estimator.** rho was not identified separately from H0's omega, so the ratio error sigma_r =
     sqrt(2 (sigma_w^2 - sigma_m^2)) is now the primary quantity, and omega cancels from it. rho is a derived
     diagnostic. Each printed SD is read as one replicate's (an SE is converted), with an unstated n of 2.
  2. **The floor.** sigma_r >= sqrt(2) x sigma_m, the measurement-noise ratio SD; never 0.
  3. **Few clusters.** G_min = 5 studies per lane, justified; below it, UNDECIDED. A t-interval over study-level
     estimates replaces the bootstrap.
  4. **LOSO symmetric.** Layer 1 must agree at LOSO_min, the point estimate and LOSO_max. Only the point estimate is
     used downstream; the bootstrap draws are withdrawn.
  5. **One decision must not refute.** REFUTED needs E >= 3, and every verdict is worded "on E of 17".
  6. **The laboratory offset.** A tau = 0.964 sensitivity on design absolute levels; a REFUTED verdict must hold under
     it.
  7. **The last choices.**
     - #108's newer commit is the primary for a rerun.
     - Non-thiol decisions use their lane's sigma_r at >= G_min studies, else the pooled sigma_r, else UNDECIDED.

  **Also in v3.1:** stage 0a's ensemble is declared exactly (sec. 5.0.2): the envelope sampler, B1 / B13 at #108's
  pinned points, its declared departures, and its 12 + 8 computable decisions. The sigma_r sensitivities replace the
  rho ones. The step-0 cost is re-projected (23 arms).
