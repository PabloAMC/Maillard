# Build Wave B7 — the furanic channel — HOLD-OUT REPORT

Generated 2026-08-29 at `49df685` on `audit-remediation`. Pass band 3×, taken unchanged from the other module scorecards.

Pre-registered in `results/validation/kinetic_core_b7_prereg.md`. The seventh test -- the cutover exam's 5 HMF + 2 DMHF bundles -- is scored by scripts/generators/generate_cutover_final_exam.py and reported beside this panel. Under the Amendment 9 clause 1 / Amendment 10 clause 1 precedent those seven rows are SEEN_DIAGNOSTIC (prereg sec. 0.1) and may never gate.

## Summary

| # | hold-out | pre-registered | outcome |
|---|---|---|---|
| H1 | Kocadagli NaCl arm | 0/3 rate ratios, 2/3 conversions in band | **0/3 and 2/3** |
| H2 | Gursul 27 °C zero accumulation | PASS (< 100 µg/L) | **FAIL** — 1.17e+03 µg/L (prereg MISS, disclosed) |
| H3 | Hamzalioglu matrix-vs-water | FAIL by construction | **FAIL**, 9.33× on the collapse; the cysteine sub-test passes at 1.00× |
| H4 | shu1988 × wang2008 paired | sink REFUSED, formation right sign wrong reason | **sink REFUSED, formation FLAT** — a pre-registration correction, disclosed |
| H5 | apriyantono held-vs-drifting | 2 of 4 channels scored, direction pass | **2 of 4 scored, BOTH RATIOS 1.000 — the drifting arm does not drift** |
| H6 | norfuraneol ≫ DMHF | PASS | **PASS** |
| H7 | cutover exam, 5 HMF + 2 DMHF | ≥5 of 7 answered; HMF over, DMHF 1–3 decades under | **7 of 7 answered, 1 in band, 4 inside the declared interval — HMF direction MISSED** |

## H1_kocadagli_nacl_arm

**Role.** HOLD-OUT (Amendment 12: 'its NaCl arm ... = sharpest hold-outs')

**Pre-registered.** 0 of 3 rate-ratio cells inside the 3x band; 2 of 3 mole-conversion cells inside it (the 180 and 200 C ones), by arithmetic and not by chemistry

THE MODEL PREDICTS EXACTLY 1.000, BY CONSTRUCTION AND WITH NO FREE PARAMETER. It carries no ionic-strength term, no salt species and no activity coefficient. Nothing was run to produce this number, which is why it is a clean test: there is no way to have tuned it.

| quantity | T (°C) | measured | predicted | fold | in band |
|---|---:|---:|---:|---:|---|
| k(Fru->Int) NaCl / glucose | 160 | 3.91 | 1.00 | 3.91× | **no** |
| k(Fru->Int) NaCl / glucose | 180 | 3.88 | 1.00 | 3.88× | **no** |
| k(Fru->Int) NaCl / glucose | 200 | 4.06 | 1.00 | 4.06× | **no** |
| HMF mole-conversion NaCl / glucose | 160 | 3.50 | 1.00 | 3.50× | **no** |
| HMF mole-conversion NaCl / glucose | 180 | 1.90 | 1.00 | 1.90× | yes |
| HMF mole-conversion NaCl / glucose | 200 | 1.06 | 1.00 | 1.06× | yes |

> The size of a declared gap, not the quality of a fit. Kocadagli's own finding is that NaCl SWITCHES the flux: it multiplies k(Fru->Int) by 3.9-4.1x, flat across 40 C -- the paper's best-behaved number -- while simultaneously HALVING k(Glc->3-DG). A model with no ionic term cannot express a switch.

## H2_gursul_aktag_2020_27C_zero_accumulation

**Role.** HOLD-OUT (Amendment 12: "Gursul's 27 C zero-accumulation row")

**Pre-registered.** PASS -- predicted HMF below the declared 100 ug/L floor

| T (°C) | predicted HMF (µg/L) | fructose-limb share | 3-DG-limb share |
|---:|---:|---:|---:|
| 27 | 1171 | 0.0007 | 0.9993 |
| 37 | 6265 | 0.0015 | 0.9985 |

> ★ THE PRE-REGISTRATION PREDICTED A PASS AND THIS IS A FAIL. It argued from the ingested activation energies (100.4 on Fru->Int, 151.4 on Int->HMF) that 27 C / 24 weeks would land below 100 ug/L. It lands above. The pre-registration's arithmetic was done on the fructose limb alone; at 27 C the model's HMF comes almost entirely from the OTHER limb, whose terminal step (3,4-DG -> HMF) carries Ea = 0 BY THE AUTHORS' OWN CHOICE -- Kocadagli fixed it to zero with the footnote 'does not follow Arrhenius equation'. A zero-barrier terminal step cannot switch off as temperature falls, and that is exactly what this row was designed to catch: Amendment 12 calls it 'the cheapest and most informative single test in the module', and it earned the description. WHAT IT DOES NOT MEAN: it is not evidence that the ingested Ea are wrong. It is evidence that carrying an author-fixed Ea = 0 on a terminal step, which K5a sec. 7.3 shows is the ONLY defensible value for that edge because no usable one exists in any paper of the cluster, leaves the node with no low-temperature shut-off. The gap is G1 (no usable HMF activation energy in any real matrix), measured for the first time.

Declared floor **100 µg/L** — Gursul Aktag's printed LOD/LOQ pair is INTERNALLY IMPOSSIBLE -- an LOD of 10 mg/L against an LOQ of 30 ug/L, i.e. three orders the wrong way round and above its own calibration range -- so no detection limit can be taken from the paper and this floor is declared here.

> The paper's own sugar table is not transcribed in the dossier, so this is a DECLARED representative juice charge. The row is a TEMPERATURE-axis ceiling test and is robust to the charge within a factor of two; it is not a level test and is not scored as one.

> The 37 C row is reported ALONGSIDE and is NOT scored: Gursul Aktag's 37 C maxima (16.2 / 3.8 / 12.2 mg/L across three juices) span 4.3x between juices on composition the model does not carry, and the 27 C row is the one Amendment 12 names. Printing the 37 C prediction anyway is what makes the 27 C result interpretable: it shows the model has NOT simply switched HMF off at low temperature.

## H3_hamzalioglu_matrix_vs_water_selectivity

**Role.** HOLD-OUT (Amendment 12: the same-method matrix-vs-water pair)

**Pre-registered.** FAIL by construction, with the missing ingredient named: the model carries ONE amine on the sulfur lane and NO moisture term, so it predicts a collapse factor of exactly 1.000

Measured Cys/Lys selectivity: **11.4×** in water, **1.23×** in roasted coffee — a **9.3× collapse**. Predicted collapse **1.000**, fold error **9.33×**.

> ★ AN HONEST NUANCE THE PAIRED FORM WOULD HIDE. The CYSTEINE constant itself does NOT move between water and coffee at 25 C -- both are 0.103 day^-1 -- so the model's no-moisture-term prediction of 1.000 is RIGHT for this sub-test. What collapses is the SELECTIVITY, and it collapses because arginine and lysine go UP in the low-moisture matrix by 6.7x and 9.3x, not because cysteine goes down. The model cannot express that because it carries neither amine. Reporting only the paired collapse would have scored a coincidental pass as part of a fail; reporting only the sub-test would have scored a coincidence as a pass.

**Why the model cannot form the ratio.** The sulfur lane tracks exactly one amine, cysteine. Arginine and lysine are not species in it, so a Cys : Arg : Lys ratio has no representation at all. And the sink edge carries no moisture, matrix or water-activity term, because the only source that measures the axis is this hold-out.

## H4_shu1988_x_wang2008_paired_sink_and_net_pH

**Role.** HOLD-OUT, PAIRED, cross-paper (Amendment 12)

**Pre-registered.** sink half REFUSED (Edge C rate is exactly zero); formation half predicted to RISE with pH from the sulfur lane's base-favoured 2,3-enolisation, i.e. the right sign for the wrong reason

> ★ THE PRE-REGISTRATION'S REASONING ON THE FORMATION HALF APPLIED TO THE WRONG SYSTEM, AND THIS IS DISCLOSED RATHER THAN QUIETLY DROPPED. It argued that the model would predict a base-favoured rise through the sulfur lane's 2,3-enolisation pH factor. That factor sits on the PENTOSE Amadori route, and Wang & Ho's system contains NO SUGAR AT ALL -- it is 1.4 M methylglyoxal plus 1 M cysteine. In that system Edge A is silent, Edge B carries no pH term (single-temperature, single-pH-fit-forbidden), and Edge C is zero, so the model's actual prediction is FLAT. The measured curve rises 3.5x. Recorded as prereg correction C1.

### Sink half — REFUSED

Edge C ships at EXACTLY ZERO, so the model makes no claim about the pH shape of a sink it does not run. This is a SHARPER refusal than the pre-B7 one: the species (DMHFS), the edge and the balanced stoichiometry all exist and what is missing is named -- a magnitude, and specifically Haleva-Toledo et al. 1999 (JAFC 47:4140-4145), the only identified source that quantifies DMHF inhibition by cysteine against pH IN A BUFFER.

> Shu & Ho's pH 7.1 zeros are ambiguous ON THE AUTHORS' OWN READING: they argue 'secondary reactions occurred readily at this pH' and that the primary products were CONSUMED rather than never formed. The paper cannot separate 'the sink is off' from 'the sink ran and its products were eaten', so even a model with a sink would be scored against net survival rather than coupling rate.

### Formation half — scored

| pH | measured (mg/mol MG) | predicted (mg/mol MG) |
|---:|---:|---:|
| 3 | 20 | 0.6203 |
| 5 | 34 | 0.6203 |
| 8 | 70 | 0.6208 |

Measured pH 3→8 trend **3.50×**, predicted **1.001×**.

> FLAT against a measured 3.5x rise. The model has no pH term anywhere on the DMHF node, and it has none for a stated reason: no paper in either cluster varies pH within one system on a furanone edge with a measurable rate, so a pH term would be fitted across labs and matrices at once.

**Paired verdict.** The pair was designed to say WHICH of the sink and the formation edge is mis-signed if a model reproduces one and not the other. This model reproduces NEITHER, and the pair says why: it has no sink at all and no pH term on the formation edge. That is a cleaner answer than a coincidental half-pass would have been.

## apriyantono1993_xylose_lysine_pH_trajectory_pair

**Role.** HOLD-OUT, pH-TRAJECTORY, scored as ONE PAIRED LOG-RATIO TEST (Amendment 12, the named role)

**Pre-registered.** 2 of 4 chemical channels scored, 2 refused with the missing species named; furfural UP on drift with a large magnitude under-shoot; norfuraneol DOWN on drift; the pH state itself drifting in the right direction but far less than 2.3 units

**Declared limitation.** THE PAPER'S AMINE IS LYSINE, which lives only on the acrylamide lane while the pentose lives only on the sulfur lane -- the two do not compose. The pair is therefore run SUGAR-ONLY and the amine's contribution is an uncontrolled, declared difference. Only direction and order of magnitude are scored, which is also what K5b sec. 8.1 recommends: the ratio form is immune to the single-alkane internal standard and to the SDE recovery bias, because both arms went through the identical isolation.

**And a second.** The HELD arm received unreported 3 M NaOH throughout the hour, so it is not sodium- or volume-matched to the drifting arm. K5b sec. 8.1's own caveat, carried.

| channel | measured drift/held | predicted | fold | scored |
|---|---|---:|---:|---|
| total volatiles | 143× | — | — | REFUSED |
| 2-furfural | 274× | 1 | 274× | yes |
| N-containing volatiles (16 compounds) | 0.01333× | — | — | REFUSED |
| norfuraneol | DOWN (trace -> not detected) | 1 | — | yes |

* **total volatiles** — refused: the core has no total-volatile observable. Apriyantono's total is 58 identified compounds across nine classes, most of which are not species in any lane.
* **N-containing volatiles (16 compounds)** — refused: no pyrazine, pyrrole, pyridine or pyrrolizine species exists in any lane. Four whole compound classes go from present to not-detected in this experiment and the core can see none of them.

### The pH-state channel — this is what exams B2.2/B2.3

Measured: held at 5.0, drifting **4.9 → 2.6** (2.3 pH units). Predicted drifting endpoint: **4.899999992456287**.

> THE B2.2/B2.3 pH STATE ITSELF. A buffered pH ladder asks 'does the model get the pH-5 rate right'. This pair asks 'does the model's internal pH EVOLVE correctly, and does the chemistry follow it' -- and a model that treats pH as a fixed input can pass every point of a ladder and still fail this.

> Apriyantono's norfuraneol cell must NOT be scored against the norfuraneol >> DMHF ordering (K5b sec. 8.5): both terms are at the detection floor, the amine is lysine rather than glycine or alanine, the isolation is Likens-Nickerson SDE (close to a worst case for a water-miscible labile furanone), and the authors themselves read the trace as consumption into coloured products. The paper is SILENT on the ordering, not contradictory.

> ★ THE PAIR COLLAPSES ONTO THE pH-STATE CHANNEL, AND THAT CHANNEL FAILS OUTRIGHT. Both chemical channels come back at a ratio of 1.000 to ten significant figures -- not because the model has no pH response (it has: `r_pent_tdp` is acid-catalysed and `r_pent_dpo` is base-catalysed, and they are the two halves of the very enolisation fork Apriyantono invokes) but because THE MODEL'S DRIFTING ARM DOES NOT DRIFT. Its predicted pH moves by about 1e-8 units against a measured 4.9 -> 2.6. THE CAUSE IS THE DECLARED LIMITATION, NOT A pH-STATE DEFECT: the acid that drives Apriyantono's 2.3-unit fall is generated by the xylose/LYSINE Maillard chemistry, and lysine cannot be put in the same lane as a pentose, so this pair had to be run SUGAR-ONLY. A sugar-only pot at 100 C for one hour makes almost no titratable acid, so B2.2's pH state is asked to integrate a source term that is nearly zero and correctly returns nearly nothing. WHAT THE TEST THEREFORE ESTABLISHES: not that the pH state is wrong, but that THE CORE CANNOT REACH THIS EXPERIMENT AT ALL. That is a lane-algebra gap, it is now measured rather than suspected, and closing it means putting a second amine on the sulfur lane -- which is a modelling addition, not a conservation fix, and therefore not this wave's licence.

**Pre-registered outcome.** PARTIAL MISS, disclosed. The pre-registration predicted 'DIRECTION PASS, MAGNITUDE UNDER-SHOOT' on furfural and on the pH state. The magnitude under-shoot is right and is total; the direction is not merely wrong, it is UNRESOLVED, because a ratio of 1.000000000 has no direction. Scoring it as a direction pass would have been the dishonest reading available here.

**Record correction.** research_round3_channels.md sec. C.2 records this paper as 'RATIO-ONLY' on the strength of its abstract's g/kg shares. Its Table 1 reports all 58 compounds in nmol per mol of xylose -- absolute molar yields on the limiting sugar, against an internal standard, +/-16 %. The furfural pH effect is 274x, not the ~1.9x the mass shares imply, because the denominator collapsed 143x at the same time. Amendment 12 already overturns the verdict; this panel is the first artefact to score against it.

## H6_norfuraneol_much_greater_than_DMHF

**Role.** HOLD-OUT, ORDINAL ONLY (Amendment 8, kept by Amendment 12)

**Pre-registered.** PASS by a large margin

Norfuraneol **1.175** mmol/L vs DMHF **0.02831** mmol/L at Blank's own conditions — **PASS**.

> Blank 1997 p. 2646: 'in all samples analyzed, 4-hydroxy-5-methyl-3(2H)-furanone was the main reaction product (data not shown)'. Corroborated across two papers and six systems and QUANTIFIED IN NEITHER. Wang & Ho, Poisson and Shu & Ho do not measure norfuraneol at all. The ordering is supported twice and quantified zero times, so the hold-out can only ever be 'norfuraneol > DMHF' and never a ratio -- the predicted ratio above is printed for information and is NOT scored.

> The DMHF magnitude WAS fitted to this same paper's Table 1. What was not fitted, and what this tests, is the RELATION between DMHF and a species whose constants come entirely from B2's sulfur fit on a different corpus. A calibration that had put DMHF above norfuraneol would have failed here while still fitting Table 1 perfectly.

## H7_cutover_exam_delta

**Role.** SEEN_DIAGNOSTIC under prereg sec. 0.1 -- the exam artefact prints every measured value and this builder saw it while locating the seven refused bundles. These rows are reported and MAY NEVER GATE.

**Pre-registered.** at least 5 of the 7 rows become ANSWERED; HMF OVER-predicted; DMHF under-predicted by 1-3 decades; the two Schibilsky pH arms predict the SAME HMF because the module carries no pH term; no previously-answered row moves by more than 1.5x

| bundle | compound | T (°C) | measured (ppb) | predicted (ppb) | fold | dir | band | in interval |
|---|---|---:|---:|---:|---:|---|---|---|
| `mp_holdout_fructose_asparagine_180C_Lin2022` | 5-Hydroxymethylfurfural (HMF) | 180.0 | 1.228e+04 | 1942 | 6.32× | UNDER | fail | yes |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019` | 5-Hydroxymethylfurfural (HMF) | 130.0 | 5.725e+04 | 2.796e+04 | 2.05× | UNDER | **PASS** | yes |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019` | DMHF | 130.0 | 1153 | 21.84 | 52.80× | UNDER | fail | no |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019` | 5-Hydroxymethylfurfural (HMF) | 130.0 | 1.006e+05 | 2.796e+04 | 3.60× | UNDER | fail | yes |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019` | DMHF | 130.0 | 5894 | 21.84 | 269.86× | UNDER | fail | no |
| `mp_holdout_glucose_asparagine_180C_30min_water_Chang2021` | 5-Hydroxymethylfurfural (HMF) | 180.0 | 7000 | 2108 | 3.32× | UNDER | fail | yes |
| `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021` | 5-Hydroxymethylfurfural (HMF) | 121.0 | 1.74e+04 | 1459 | 11.93× | UNDER | fail | no |

Every one of the seven was `out_of_envelope` before this wave.

> ★ THE PRE-REGISTRATION PREDICTED HMF **OVER**-PREDICTION AND EVERY ROW UNDER-PREDICTS. The reasoning was that the module has no validated HMF sink at cooking temperature (K5a G2: the 50-150 C window is empty), so a source-only node must over-shoot. It does not, and the missing sink is therefore NOT the binding constraint -- THE SOURCE TERMS ARE. The two formation limbs are ingested from an AMINE-FREE, FREEZE-DRIED, essentially anhydrous glucose melt and are being asked to run in aqueous and food matrices; the melt -> matrix transfer loses more flux than the absent sink adds back. That is a quantified, directional result about a transfer nobody had tested, and it is the opposite of what this wave expected.

### The pH pair — a pre-registered structural miss

The two Schibilsky bundles differ ONLY in pH (5.0 vs 8.0) and the module carries NO pH term on any furanic edge -- K5a declared gap G8: six distinct pH values appear across the seven papers of the cluster and NO SINGLE PAPER VARIES pH, so a pH term would have to be fitted across labs and matrices at once, which k3 sec. B.2 forbids at family level.

Predicted HMF and DMHF are **identical** at pH 5.0 and pH 8.0. Measured pH effect: **1.76×** on HMF and **5.11×** on DMHF. The size of the measured pH effect IS the size of the gap, now measured for the first time on both compounds at once.

### What the trunk perturbation cost

The furanic block hangs on the TRUNK, so four of its steps put a new drain on B1-fitted species and every previously-answered row moves a little. B1 is NOT refit -- this wave has no licence to move its constants -- so the movement is reported rather than absorbed.

**18** previously-answered rows moved; largest **1.2621×** against a pre-registered ceiling of 1.5× — **within**.

| bundle | compound | before | after | fold |
|---|---|---:|---:|---:|
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7` | 2-Furfurylthiol (FFT) | 0.01886 | 0.01495 | 1.2621× |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3` | 2-Furfurylthiol (FFT) | 94.23 | 86.44 | 1.0902× |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7` | 2-Methyl-3-furanthiol (MFT) | 10.57 | 9.918 | 1.0661× |
| `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3` | 2-Methyl-3-furanthiol (MFT) | 32.9 | 31.82 | 1.0340× |
| `mp_holdout_fructose_asparagine_180C_Lin2022` | Acrylamide | 228.3 | 225.4 | 1.0125× |
| `mp_holdout_glucose_asparagine_180C_30min_Chang2021` | Acrylamide | 4041 | 4031 | 1.0025× |

