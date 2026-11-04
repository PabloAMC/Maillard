# Case M1: thiol loss and its oxidant; pH drift

**Rescoped 2026-11-04.** This case drafted the plan's original M1 claim. An audit of the six refused sink waves
(CF-107) then showed four of them are one degenerate fit optimum, so their refusals are not evidence against the
structures; M1 now compares combinations C1 - C3. Its evidence led the plan to rescope M1 to "thiol loss or slow supply" (no-fit checks M1.0 first) and to move part D (the pH balance) into M1b; see STATUS_AND_PLAN sec. 5, W2. The claim below is kept as drafted.

*Plan step M1 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E12 (F28-a, -b, -c) and CF-084).
Items 1 - 5 of sec. 7, written 2026-11-04 before any result of the change exists. No engine was run for this case,
with or without a changed switch. Items 6 - 9 are appended later and never edited into this part.*

**Seen by the author.** Held-out or keyed held-out results read while writing: the W1.1 score
(`results/validation/kinetic_core_contrast_benchmark_score_v1.md`), its Grosch 1992 Table 2 items (measured and
engine signs) and its Meynier 1995 items; the engine-vs-measured Yiltirak 2026 numbers quoted in CF-068, in WAVES.md
(B10, B11, B16, B51), in `SULFUR_FLUX_MAP.md` sec. 5 and in the INTRODUCTION figure 04 caption; the B27 gate's T1
ratios at B9 (`kinetic_core_b27_ship_rule.md`); EX-056's per-fold results. Measured values only (no engine number)
for the v2 E12 TEST contrasts F28#8, F28#20, F28#25 and F28#26, read in the atlas file.

## 1. The claim

The meaty thiols (MFT, FFT) are lost in a heated cysteine + pentose pot through an oxidation channel whose oxidant the
pot itself makes, plus the air in the vessel, and the disulfide this makes gives the thiol back. The shipped engine
instead removes thiol through one first-order, pot-independent step and starves the disulfide channel of oxidant.
The pH the pot drifts to is part of the same claim: it sets the thiolate share that every loss step uses.

Network steps that change (names from [sulfur.py](../../../src/kinetic_core/sulfur.py)):

| part | step(s) | today | M1 |
|---|---|---|---|
| A. oxidant from the pot's sugar intermediates | `ch_redox_mp3p` (B27; phi = `k_redox_mp3p`), or a new step from DPO / TDP / DDP | inert (phi 0) | on, with a declared source and band (open question 1) |
| B. headspace O2 | `OXR`, `OXV`, `ox_supply`, `ch_cys_ox`, `ch_red_ox_*` (B11) | inert; `OX` seeded at the ambient 1.0 mM, zero in all 14 fed mol % rows | on, `OXR` charged from each pot's vessel volumes |
| C. reversible disulfide | `ch_dimer_release_mft` / `_fft` | `k_dimer_release` = 0 | declared from Guth 1995, not fitted |
| D. pH drift | the sulfur lane's charge balance (`ph_state`, B51), `ph.trunk_solutes_in_balance` | balance on in the sulfur lane; trunk solutes off | on by default, with the trunk's acids as titratable pools (open question 4) |
| not changed by the claim | `ch_mft_decay`, `ch_thiolate_loss_*`, `ch_dimer_*` constants | shipped B60 | whether a refit moves them is open question 2 |

## 2. All the evidence, by kind

Laboratory ids from [lab_ids.yml](../contrast_benchmark/lab_ids.yml). "v1" = [split.yml](../contrast_benchmark/split.yml)
(contrast_split_v1); "v2 E12" = [split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml), error E12,
whose DEV labs are DFA Garching, Nestle / Cerny and Reading, and whose TEST contrasts are F28#8, F28#20, F28#25,
F28#26. Atlas cells are 0-based indices of [F28_thiol_losses.yml](../pathway_atlas/F28_thiol_losses.yml) unless named.
"Cal. input" = listed in v2's `calibration_input_papers` (never a v2 TEST contrast).

### 2.1 Isotope labelling

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Hofmann 2001 / 2002 (F28#14; hofmann2001, hofmann2002 dossiers) | [2H2]-FFT ends on coffee melanoidin as covalent thioether; disulfide < 1.5 % of the loss; DTE gives no thiol back (30 C) | against "the disulfide is the main sink"; storage temperature, coffee matrix | DEV (v1, v2 E12); cal. input |
| Davidek 2006a / 2008 (F07#7, F07#8) | acetic acid ~70 % from sugar C1 / C2 (the 1-deoxyosone cleavage) | for part D: names the main acid of the drift | DEV (v1) |
| none | no labelling of thiol loss in a sugar + cysteine pot at process temperature (F28 `unknown_regimes`) | gap | -- |

### 2.2 Fed intermediates

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Whitfield & Mottram 1999 (whitfield1999 dossier; B27 prereg sec. 1) | fed norfuraneol + Cys or + H2S, 140 C, pH 4.5: 26 disulfides; ~7.5 % of all thiol and ~35 % of MFT as disulfide; the authors argue aerial oxidation is unfavourable at high H2S and name the alpha-dicarbonyls as the oxidant | for parts A and C; against part B in H2S-rich pots | DEV (Reading); cal. input |
| Kumazawa 2003 (F28#3; fit rows `kumazawa2003_*` in [thiol_submodel.yml](../../../data/lit/fit_rows/thiol_submodel.yml)) | FFT alone in buffer, 121 C, pH 3 - 7, canned without deoxidation: FFT 99.5 -> 0.1 %, disulfide rising with pH, a stated mass-balance shortfall | for the thiolate law; the must-not-break pot for any global oxidant (B27 sec. 9.3: ambient x10 moves it -0.44 -> -3.39 sigma) | v1 TEST; cal. input (a fit row since B16). Used here only as a constraint already in every objective, never to choose a rival |
| Hofmann 1998 T3 / T4 fed precursor + H2S pots (CF-065; thiol_submodel.yml) | the 14 fed mol % rows; OX charged zero | defensible zero for H2S pots (B27 sec. 2); carry the H2S headspace bias (CF-065) | DEV (DFA); cal. input |
| Mottram 1996, Mottram 1996b, Adams 2001 (F28#21, #22, #23) | fed furyl disulfides come back as free thiol with albumin, beef or denatured ovalbumin; in beef 83 - 90 % of thiol units accounted | for part C in protein matrices; no sugar + cysteine pot | DEV (Reading) |
| Mikami 2017, Nakamura 2020 (F15#7, F15#8) | fed norfuraneol alone, pH 8, 100 C: decay slowed by degassing and by DTPA | O2 and metals act on the thiol precursor too: a side effect of part B | DEV (Ochanomizu); Mikami cal. input |

### 2.3 Trapping

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Guth 1995 (F28#1; guth1995 dossier sec. 4) | 100 C, 2 h, water, 4-vinylpyridine trap: MFT-MFT 86 % gone, 62 % of two thiols recovered; apparent loss 0.0165 /min (MFT-MFT), >= 0.034 /min (FFT-FFT); D2O puts one D in trapped MFT (hydrolysis). Without the trap the MFT-MFT pool is conserved (95.5 %) | for part C (the release rate); the release shows only when thiol is drawn off | DEV (DFA) |
| Esterbauer 1975 (F28#17) | GSH Michael adducts of enals revert (DTNB-driven), 20 C | for reversible capture as a rival form | v1 TEST (Graz); not used to choose |
| Anantharamkrishnan 2020b (F28#24) | FFT 1:1 and 2:1 adducts on beta-lactoglobulin (MS) | protein capture; no sugar pot | v1 DEV (Minnesota); cal. input |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Zhou 2022 (fit rows `zhou2022_axarp_cys_ph7_120C`, MFTD) | MFT disulfide 34.7 -> 43.7 ug/L over 30 - 120 min, SD printed | a disulfide time course: the reservoir's shape | DEV (ho_network) |
| Wang 2026b (fit rows `wang2026b_xyl_cys_120c_ph55`) | MFTD, FFTD and the mixed disulfide all rise 60 - 150 min, SD printed | as above, a second laboratory | Dalian Polytechnic (EX-056 fold); no atlas cell and no `lab_ids.yml` entry, so in no frozen split: a fit row, development by use |
| Hou 2026 (hou2026 dossier; F06#8; CF-084) | unbuffered GX-ARP 5 mM, 120 C: pH 7.0 -> 5.05 (20 min) -> 3.55 (120 min); acid keeps forming after the ARP is gone | for part D's shape test (acid lag) | DEV (ho_network) |
| Cao 2017, Zhang 2025b, Zhai 2019 ([SULFUR_FLUX_MAP sec. 4](../../guides/SULFUR_FLUX_MAP.md)) | free cysteine and sugar last hours at 90 - 120 C where the engine destroys them in seconds | for rival R3 (supply timing), not for M1 | DEV (BTBU, ho_network) |

### 2.5 Mass balance

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Grosch 1992 (F28#2) | MFT disulfide <= 0.2 ug (thiamine), n.d. with Cys, 1.6 +- 1.2 ug full system, against 20 - 100 ug MFT missing: <= 1 - 5 % of the loss (100 C reflux, pH 5.7) | against "the disulfide carries the loss" | DEV (DFA). Its Table 2 is also a keyed held-out table for E18 / CF-089, scored in W1.1 |
| van Seeventer 2001 (F28#9, #10, #11) | 50 C storage: MB and FFT losses close as thiol + disulfide; MFT's does not; air ~ argon for MFT | storage channel only (atlas convention) | DEV (Quest); cal. input |
| Mottram 1996b (F28#22) | beef: disulfide a reversible reservoir, 83 - 90 % accounted | for part C, matrix only | DEV (Reading) |
| Cao 2017 (cao2017 dossier) | 0.4 M Cys + Xyl, 90 C, 5 h: no measurable cysteine destroyed | for R3 | DEV (BTBU) |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | for / against | split status |
|---|---|---|---|
| Zhou 2023 (F28#4; [thiol_zhou2023_ph_series.yml](../../../data/lit/fit_rows/thiol_zhou2023_ph_series.yml)) | initial pH 6 / 7 / 8, ARP + Cys, 120 C: MFT disulfide share 8.6 / 6.5 / 9.6 % (thiol equivalents), nearly pH-flat while MFT swings 3x | for "a disulfide share of ~7 - 10 %"; against a thiolate-only dimer law for MFT | DEV (ho_network); cal. input |
| Zhang 2024 (F28#5) | cystine vs cysteine arm: MFT disulfide share x6.3 (54.2 vs 8.6 %), free MFT x1.5 higher | for "oxidant -> disulfide", against "oxidant -> less free thiol" | DEV (BTBU); cal. input; figure read-off only |
| Charles-Bernard 2005 (F28#0, F28#16; CF-081) | coffee solids, 25 C: air vs N2 slows loss 15 - 45x; hydroxylamine 4 - 13x; jointly 65 - 236x; DTPA no effect | for an O2-dependent loss, at 25 C in coffee | DEV (Nestle) |
| Blank 2002 (F28#7) | Fenton factorial, 37 C: FFT loss needs Fe + H2O2 + reductant | metal route exists; no Maillard pot | DEV (Nestle) |
| Ehrenberg 1989 (F28#8) | cysteine stable under air without metal; Cu / Fe autoxidise it (Michaelis-Menten, 37 C) | the form CF-068 proposes for `ch_cys_ox` | v1 DEV; **v2 E12 TEST**: see open question 3 |
| Grosch 1992 (F28#12) | SDE vs reflux: in-pot loss 0.044 /min (thiamine), 0.034 /min (+ Cys), 100 C | a measured processing loss at the engine's size, possibly thiamine-specific | DEV (DFA) |
| Elmore 2002 (F28#6, #19) | 140 C, 80 mL air headspace: MFT 807 +- 508, disulfide 202 +- 213 ng; PUFA lowers MFT 14 - > 80x | qualitative for the disulfide; lipid is E7 | DEV (Reading) |
| Xu 2011 (F28#20) | lard lowers MFT x21 / x143 | lipid, E7 | v1 DEV; **v2 E12 TEST**: read only in item 7 |
| Wang 2022 (F28#25; CF-013) | 100 - 110 C rising, 120 - 140 C peak and fall | the temperature switch of the net balance | v1 DEV; **v2 E12 TEST**: read only in item 7 |
| Liu 2023b (F28#26) | 168 C: MFT halves 20 -> 60 min, FFT flat | as above | v1 DEV; **v2 E12 TEST**: read only in item 7 |
| Kumazawa 2003 (F28#3) | pH 3 -> 7 (see 2.2) | thiolate law | v1 TEST; constraint only |
| Stack 2018 (F28#18) | thiol + quinone K falls with T | reversible capture lets go when hot (B17a) | v1 TEST; not used to choose |
| Yiltirak 2026 (fit rows `yiltirak2026_rib_cys_ph55_buffer_*`; CF-068) | equal-cook ladder 100 / 110 / 120 / 130 C, 25 + 25 mM, pH 5.5, 3 mL in 20 mL air: MFT 6.88 / 3.29 / 2.4 / 1.71 ug/L, FFT 1.28 - 1.68 | between-lab level gap x120 vs Schieberle; within-lab T shape | Reading = DEV in v1 and v2 E12; a declared hold-out of every fit (never fitted), but its engine scores have been read many times. Not "unseen" |
| Meynier 1995 (CF-084, meynier1995 dossier) | 0.2 M pyrophosphate, 140 C: measured drift <= 0.2 pH units; the engine gives -0.79 / -0.65 / -0.39 | part D's negative control | DEV (Reading) in v1; the same statement prints in the keyed held-out table `heldout_meynier1995_tab1_2_4.yml` |
| Zhou 2023 Fig. 2 (CF-084) | ARP alone, labels 6 / 7 / 8 -> 3.99 / 4.13 / 4.99 | part D, label sensitivity (engine ~4x too large) | DEV (ho_network) |
| Chen 2000 / Chen & Ho 2002 (CF-054) | same recipe, end pH 6.05 vs 4.96 | ~0.6 pH units run-to-run spread from an alkaline start | DEV (ho_network); not a fit row |

### 2.7 Fitted time courses

| source (record) | what it shows | for / against | split status |
|---|---|---|---|
| Schieberle 2000 Table IV (fit rows `schieberle2000_rib_cys_ph5_100C`; CF-013, CF-055) | 100 C, 0.5 M phosphate pH 5: MFT 4.5 / 13.8 / 156 / 179 ug at 30 / 60 / 360 / 720 min, still rising; FFT 2.0 / 3.1 / 110 / 132 | the reference pot of the plan row | DEV (DFA); cal. input; a fit row of R3En and B16 - B27 |
| US 10,076,126 (F28#13; CF-013, CF-015) | co-heated [2H] standards, 110 C: MFT still rising at 4 h; bound lambda(110 C) <= 0.02 - 0.05 /min | for a slower loss in Cys-rich pots | DEV (Nestle); not measurement-grade (P 0.15, EX-055); unscored until a tracer operator exists |

### 2.8 Earlier tests of the same parts (engine records, not literature)

| record | part | result |
|---|---|---|
| B11 (WAVES.md; `kinetic_core_b11_ship_rule.md`) | B: headspace O2 reservoir, two consumers | DO NOT SHIP: consumers to their floors, unidentified; Yiltirak and Bolton unchanged |
| B16 (`kinetic_core_b16_ship_rule.md`) | one sink barrier | DO NOT SHIP: peaks at 6 h, breaks the 145 C fed pots |
| B17 variant b (`kinetic_core_b17_ship_rule.md`) | C: reversible disulfide, fitted | DO NOT SHIP: release to its floor; the model holds 0.04 - 0.9 % as disulfide vs 6.5 - 9.6 % measured |
| B17 variant a (`kinetic_core_b17a_ship_rule.md`) | saturable thioether | DO NOT SHIP: site source dead; Stack's K(T) lets go at 100 C |
| B25 (`kinetic_core_b25_ship_rule.md`) | irreversible thiol + deoxypentosone | DO NOT SHIP: switched off at every temperature |
| B27 (`kinetic_core_b27_prereg.md`, `kinetic_core_b27_ship_rule.md`) | A: dicarbonyl-sourced oxidant (`ch_redox_mp3p`) | GATED, NOT FITTED: in ambient pots the consumers use < 1 % of the oxidant pool, so no oxidant SOURCE reaches the measured disulfide share (G1); Whitfield's 35 % needs phi = 1, the declared disqualifying ceiling (G2). Its own sec. 9 "budget" reading is withdrawn in its outcome section |
| B49 (WAVES.md) | thiol branch | "every sulfur structure tried since B16 ends at the same constant" |
| EX-040, EX-050, EX-055 | the loss law | a slower loss helps the lag but not enough; P(k_mft_decay too fast >= 3x at <= 110 C, Cys-rich) 0.7; a pot-dependent law is flexible, not identifiable (P 0.1) |
| EX-056 sec. 7b defects 2 and 3 | the loss law; pH acid allocation | NOT ADEQUATE; both defects named as levers before any rerun |
| CF-084 | D | the shipped balance fails both controls; remedy: do not score at engine in-situ pH until both pass |

**Sources left out, with reasons.** F28 storage cells at <= 60 C enter only as context (atlas rule: a storage cell
never sets a processing law). Lipid cells (Elmore 2002 PUFA, Xu 2011, Wang 2022 lipid arm, F30) belong to E7 / M5.
Thiamine pots (Grosch 1992 thiamine arms, jhoo2002) are a separate electrophile (F27). Dossiers outside the atlas
families F28, F15, F06 / F07 and CF-054 / 068 / 084 were not screened one by one for this case: a stated gap, not a
judgement that they do not bear.

## 3. The rivals

| id | explanation | what it says about the thiols | about the pH |
|---|---|---|---|
| R0 | **current model** (B60) | one first-order pot-independent loss (`k_mft_decay` 0.041 /min at 100 C), a thiolate loss, a dimer on ambient OX 1.0 mM (zero in fed pots), no release | B51 charge balance, carboxyl carried out at ARP loss, amine declared destroyed |
| R1 | **M1 as claimed** | parts A - D of sec. 1 | as R0 plus trunk acids in the balance |
| R2 | cysteine-suppressed capture (EX-050's k_E E / (1 + Cys / K_c)) | loss slower in Cys-rich pots; no O2 term; disulfide untouched | -- |
| R3 | supply timing (SULFUR_FLUX_MAP sec. 4, EX-027, B59b `CYS_AMADORI_ROUTE`) | the loss is roughly right; formation is spread over hours through the TTCA / cysteine-Amadori reservoir, so the 100 C pot keeps rising | -- |
| R4 | laboratory and work-up offsets (CF-068 rivals) | DFA reads high or Yiltirak's post-heating standard misses bound thiol; no chemistry | -- |
| R5 (pH) | amine fate (CF-084 explanation 1) | -- | released NH3 survives in buffered cysteine pots: less acid than R0 |
| R6 (pH) | acid lag (CF-084 explanation 2) | -- | acid forms downstream of the ARP, after it is gone |

R1 - R3 can hold together. The case asks which of them the decisive data need, not which one is true alone.

## 4. The decisive data

Directions are read from each mechanism, not computed.

| # | measurement (record) | R0 | R1 | R2 | R3 | separates |
|---|---|---|---|---|---|---|
| D1 | Schieberle 2000 Table IV, ratios MFT(360) / MFT(30) and MFT(720) / MFT(360) in one pot (scale cancels) | falls after ~1 h (INTRODUCTION fig. 01; B27 T1 at B9: MFT 360 / 30 = -1.80 dex) | adding oxidant and O2 makes free thiol loss FASTER; only the release slows it, and that is bounded by the disulfide share (D3). Net direction: not a rise unless the release dominates | rise continues (Cys 33 mM protects) | rise continues (slow supply) | R0 from R2 / R3; does NOT favour R1 |
| D2 | free cysteine and sugar time courses: Cao 2017, Zhang 2025b, Zhai 2019 (measured intermediates) | consumed in seconds | as R0 | as R0 | consumed over hours | R3 from all others |
| D3 | MFT disulfide share vs initial pH inside Zhou 2023 (pH 6 / 7 / 8; one species pair, so the response factor cancels across pH) | share rises with pH (thiolate law on `ch_dimer_*`; F28 engine note 6) | depends on the pH law of the oxidant-forming step; must be declared before the build (open question 1) | as R0 | as R0 | the sign of the slope tests R0's thiolate dimer |
| D4 | the share's level: Zhou 2023, Zhang 2024 (cysteine arm), Whitfield 1999, the Zhou 2022 / Wang 2026b MFTD time courses | about a decade short at B9 (11x and 8.8x; B27 prereg); at B60 residuals -2.71 sigma (Zhou) and -1.99 sigma (Zhang cysteine arm). (0.04 - 0.9 % is B17's moved vector, not R0) | higher, rising with the oxidant | as R0 | as R0 | R1 from R0 / R2 / R3, but carries the cross-species response factor (HS-SPME, RF 1) |
| D5 | Zhang 2024 cystine vs cysteine (an oxidant added, within study) | free MFT down with more oxidant (irreversible dimer) | share up, free MFT not down (reservoir) | no oxidant term | no oxidant term | R1's part C from R0; figure-only |
| D6 | Whitfield 1999, fed NF + Cys vs + H2S (same paper, duplicates) | no oxidant in fed pots: no disulfide | dicarbonyl oxidant: disulfide in both; an aerial-O2 source predicts less at high H2S | -- | -- | part A from part B |
| D7 | Yiltirak 2026 vs Schieberle 2000 at 100 C / 4 h (between laboratories; vessel air 17 mL over 3 mL vs a closed autoclave) | no vessel term, so no O2 gap (the joint thiol fit s042, not B60, predicts Yiltirak HIGHER, x3.7: CF-068) | O2 lowers Yiltirak | -- | -- | R1 vs R4: **not separable on disk** (no within-study O2 arm at process temperature; CF-068, F28 `unknown_regimes`), unless Xu 2010 (140 C, air vs 10 - 40 MPa N2: neither thiol nor disulfide changes) is such an arm: whether the N2 displaced the air is UNVERIFIED |
| D8 | pH: Meynier 1995 drift <= 0.2 (buffered); Hou 2026 trajectory; Zhou 2023 label sensitivity | too acidic (Meynier -0.4 to -0.8); prompt acid (Hou) | adding trunk acids pushes MORE acidic: Meynier worse unless R5 is also built | -- | -- | R5 is tested by Meynier, R6 by Hou's late acid |

**What the data cannot decide.** R1 against R4 on absolute levels (D7). R1's part B (O2) against part A
(dicarbonyl oxidant) except through D6, one paper at one pH. R2 against R3 on the Schieberle shape alone (D1): only D2
separates them. And the disulfide share bounds part C: if 7 - 10 % of MFT is held as disulfide (Zhou, Zhang), a
release that returned all of it would raise free MFT by a factor of about 1.07 - 1.1, far below the decade-scale
misses of D1 and D7. Part C alone cannot be the fix; it can only be judged on D3 - D5.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1 Schieberle 2000 Table IV | no SD or replicate count printed; only the read error (half the last digit) | thiol_submodel.yml line ~425 comment; schieberle2000 dossier. The case uses the measurement model's per-series noise and says so |
| D2 Cao 2017, Zhang 2025b | Zhang 2025b: n = 3, mean +- SD; its cysteine and xylose calibration curves print identical. Cao 2017: [UNVERIFIED, not read for this case] | zhang2025b dossier secs. 2 - 3; CF-063 |
| D3 / D4 Zhou 2023 | mean +- SD, n = 3 per cell; Duncan letters unreliable | thiol_zhou2023_ph_series.yml header |
| D4 Zhou 2022, Wang 2026b | printed SD per point | thiol_submodel.yml rows |
| D4 / D6 Whitfield 1999 | duplicate analyses, both values printed | whitfield1999 dossier, "Replication" |
| D4 / D5 Zhang 2024 | triplicate, error bars drawn; figure read-off only | Zhang2024 dossier sec. 2a |
| part C input: Guth 1995 | no replicates, no SD; pH not stated; one end point, first order assumed | guth1995 dossier (summary table, secs. 4 and caveats) |
| D7 Yiltirak 2026 | mean +- SD, three heatings | thiol_submodel.yml line ~634; yiltirak2026 dossier |
| D7 across laboratories | between-lab spread from H0 (EX-059, pending) | -- |
| D8 Meynier 1995 | "no change of more than 0.2 pH unit occurred after heating" | meynier1995 dossier; heldout_meynier1995_tab1_2_4.yml |
| D8 Hou 2026 | triplicate; per-point pH SD [UNVERIFIED: not found in the dossier] | hou2026 dossier |
| D8 Zhou 2023 Fig. 2 | sd = half the drawn error bar (three repeats), read error 0.02 pH | thiol_zhou2023_ph_series.yml header |
| Kumazawa 2003 (constraint) | no error bars, SD or n | kumazawa2003 dossier |

## Open questions for the build

1. **Which oxidant step?** B27's co-product of `r_nf_mp3p` was refused by its own gate. A new source from the
   deoxypentosones (the dicarbonyls themselves as the oxidant pool, at mM levels) is a different structure and needs its
   own declared pH law (D3) and stoichiometry before it runs. B27 was gated at B9; whether G1 still holds at B60 is
   UNVERIFIED.
2. **Does `k_mft_decay` move?** If the refit is free to lower it, R1 and R2 share their effect on D1 and the
   comparison loses its separation. Declare before the run whether it is frozen, banded (EX-050: Cys-conditional at
   most) or free.
3. **Ehrenberg 1989 is a v2 E12 TEST contrast (F28#8).** CF-068's remedy takes `ch_cys_ox` in its form. Using it as a
   constant source spends a TEST contrast; either take the O2 consumer from a DEV source or register the choice.
4. **What "on by default" means.** The sulfur lane already solves the charge balance; CF-084 says not to score pots at
   engine in-situ pH until the Meynier and Hou controls pass. Is part D the trunk lane switching from declared to
   charge-balance pH, `ph.trunk_solutes_in_balance` on, or both? And is R5 (amine fate) built with it?
5. **Release stoichiometry.** `ch_dimer_release_*` returns two thiols; Guth recovers 1.43 per MFT-MFT lost (hydrolysis
   gives one thiol and a sulfenic acid). Which is declared?
6. **Vessel volumes.** Only Hofmann 1998 and Yiltirak state fills; Schieberle 2000's autoclave fill is an author request
   (CF-068). B11's "declared default" for the rest decides most of part B.
7. **The INTRODUCTION figures.** "480x -> 18x" (fig. 04) and B16's recorded "198x -> 14x" differ; the figure's source
   table was not found [UNVERIFIED].
