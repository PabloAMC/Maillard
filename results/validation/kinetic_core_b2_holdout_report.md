# Kinetic core, Build Wave B2 -- the pre-registered hold-out scorecard

Predicted from the FROZEN fit in `results/validation/kinetic_core_b2_fit_report.json`. **No parameter was changed after the fit**; there is no optimiser in the scoring script.

**GATING: 7 / 27 passed.** Diagnostic: 1 / 7. Unscoreable or explicitly deferred: 3.

## Pre-registered failures (written down before the numbers)

- van Seeventer 50 C (channel carries no rate on any FIT row)
- Hofmann dry-180 C (no water-activity term in the module)
- Hofmann pH-3/pH-7 MFT (the B2.5 sign conflict; at most one of Hofmann and Zhou can be satisfied)
- Zhang Fig. 2 (methionine is not in the state vector)

## Qualifications the pass/fail column hides

- **`vanseeventer_50C_MFT_per_day` (PASS)** -- THIS PASS SHOULD NOT BE COUNTED AS ONE. It was pre-registered as a failure and it passed the 3x band from the WRONG SIDE, by over-predicting: the model destroys ~99% of the MFT per day against a measured 59%, and it does so through the unassigned-decay and dimerisation channels, NOT through the C-5 oligomerisation channel that van Seeventer actually measured -- which still carries a rate of exactly zero. The FUNCTIONAL FORM is also still wrong: the measured loss is ZERO ORDER and every channel here is first or second order. The row tests the functional form, and on that test the model fails while landing inside the magnitude band by coincidence.

- **`cerny2007_branch_responsiveness` (FAIL)** -- A FAIL THAT CARRIES THE WAVE'S MAIN ARCHITECTURAL RESULT. This is the row a model with fixed branch fractions fails BY CONSTRUCTION, because such a model predicts exactly 1.00x. This model predicts 1.49x against a measured 3.07x. So the branch fraction DOES respond to concentration -- the architectural requirement is met and the no-fixed-fraction design is vindicated -- but it responds about half as strongly as measured. The residual is a magnitude error in the relative reaction ORDERS of the two routes, not the structural error the row was designed to catch. Reported as a fail, because it is one.

- **`zhou_pH8_FFT` (FAIL)** -- THE DOMINANT FAILURE MODE OF THE WHOLE WAVE, and it is one thing rather than many. The structural pH term uses acid-catalysis proportional to [H+] and base-catalysis proportional to the free-amine fraction, i.e. ONE DECADE PER pH UNIT. That gives the right SHAPES -- FFT monotone down, MFT peaked, both confirmed by the shape rows that pass -- but far too much SLOPE, so every absolute level three pH units from the fitted pH collapses or explodes. Most of the failed rows in this table are that single defect seen from different angles: Zhou pH-6 and pH-8 absolutes, the 582-vs-665 pair, the Hofmann pH-3/pH-7 rows, and the Meynier folds. A sub-decade effective slope -- which is what a partly rate-limiting catalysis step gives -- is the obvious next structural change, and it is NOT a free pH parameter.

- **`hofmann2002_brew_80C_FFT` (FAIL)** -- THE FAILURE THE DECLARATION PREDICTED AND ASKED FOR. Its own words: 'it is the one the model will get WRONG in the informative direction'. The model loses FFT 18x too fast in the real brew. The declaration says a failure here localises the electrophile-pool depletion term, and it does: the brew's pool was partly consumed during extraction, so its effective [E] is far below the nominal site density this run was given. That is a STATE error (the wrong initial electrophile pool), not a rate error, and the reversible release step added late in this wave makes it more tractable rather than less.

## Row by row -- no averaging, no dropping

| row | group | gating | observed | predicted | fold | pass |
|---|---|---|---:|---:|---:|---|
| `zhou_pH8_MFT` | Zhou pH-8 column | yes | 0.004604 | 0.001089 | 0.24x | **FAIL** |
| `zhou_pH8_FFT` | Zhou pH-8 column | yes | 0.002849 | 1.23e-07 | 0.00x | **FAIL** |
| `zhou_pH8_MFTD` | Zhou pH-8 column | yes | 0.0002212 | 4.126e-06 | 0.02x | **FAIL** |
| `zhou_pH8_ACTZ` | Zhou pH-8 column | yes | 0.00458 | 2.245e-06 | 0.00x | **FAIL** |
| `zhou_pH8_FUR_arp_alone` | Zhou pH-8 column | yes | 0.004544 | 1.868e-05 | 0.00x | **FAIL** |
| `zhou_MFT_shape_pH8_over_pH7` | Zhou pH-8 column | yes | 0.3309 | 0.5991 | 1.81x | PASS |
| `zhou_FFT_shape_pH8_over_pH7` | Zhou pH-8 column | yes | 0.4291 | 0.0001603 | 0.00x | **FAIL** |
| `zhou_pH6_MFT` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.006105 | 8.679e-05 | 0.01x | **FAIL** |
| `zhou_pH6_FFT` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.007127 | 0.0009885 | 0.14x | **FAIL** |
| `zhou_pH6_MFTD` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.0002638 | 2.268e-08 | 0.00x | **FAIL** |
| `zhou_fig3_grid` | Zhou Fig. 3 grid | no | -- | -- | -- | _not scored_ |
| `zhou_582_vs_665_consistency` | Zhou cross-system pair | yes | 0.8757 | 3.747e-05 | 0.00x | **FAIL** |
| `whitfield_2001_pH65_collapse` | Whitfield 2001 pH collapse | yes | 150 | 0.2845 | 0.00x | **FAIL** |
| `cerny2003_NF_share_ceiling_95C` | Cerny 2003, 95 C | yes | 0.07 | 0.383 | 5.47x | **FAIL** |
| `cerny2003_intact_skeleton_share_95C` | Cerny 2003, 95 C | yes | 0.93 | 0.617 | 0.66x | PASS |
| `route_mix_moves_with_temperature` | Cerny 2003, 95 C | no | -- | -- | -- | **FAIL** |
| `hofmann_dry_180C` | Hofmann dry-180 | yes | -- | -- | -- | _not scored_ |
| `hofmann_ribose_pH3_MFT` | Hofmann pH axis | yes | 0.004844 | 0.0003912 | 0.08x | **FAIL** |
| `hofmann_ribose_pH3_FFT` | Hofmann pH axis | yes | 0.002006 | 0.03273 | 16.32x | **FAIL** |
| `hofmann_ribose_pH7_MFT` | Hofmann pH axis | yes | 0.000219 | 0.01319 | 60.22x | **FAIL** |
| `hofmann_ribose_pH7_FFT` | Hofmann pH axis | yes | 0.0001051 | 7.369e-06 | 0.07x | **FAIL** |
| `meynier_MFT_pH_fold` | Meynier directional | no | 152 | 0.08669 | 0.00x | **FAIL** |
| `meynier_FFT_pH_fold` | Meynier directional | no | 6.1 | 80.48 | 13.19x | **FAIL** |
| `meynier_FUR_pH_fold` | Meynier directional | no | 25 | 125.5 | 5.02x | PASS |
| `hofmann2002_brew_80C_FFT` | thiol sinks | yes | 0.023 | 0.4131 | 17.96x | **FAIL** |
| `vanseeventer_50C_MFT_per_day` | thiol sinks | yes | 0.59 | 0.9942 | 1.69x | PASS |
| `zhang_115C_MFT_consumed_share` | thiol sinks | yes | 0.11 | 0.1208 | 1.10x | PASS |
| `zhou_120C_dimer_share_pH6` | thiol sinks | yes | 0.086 | 0.0005226 | 0.01x | **FAIL** |
| `zhou_120C_dimer_share_pH7` | thiol sinks | yes | 0.065 | 0.01378 | 0.21x | **FAIL** |
| `zhou_120C_dimer_share_pH8` | thiol sinks | yes | 0.096 | 0.007579 | 0.08x | **FAIL** |
| `zhou_dimer_share_pH_invariance` | thiol sinks | yes | 1.48 | 26.37 | 17.82x | **FAIL** |
| `cerny2007b_control_no_cysteine` | Cerny single-route controls | yes | 0.99 | 1 | 1.01x | PASS |
| `cerny2007b_control_no_thiamine` | Cerny single-route controls | yes | 0.95 | 1 | 1.05x | PASS |
| `cerny2007_1x_xylose_share` | Cerny concentration pair | yes | 0.15 | 0.4561 | 3.04x | **FAIL** |
| `cerny2007_2x_xylose_share` | Cerny concentration pair | yes | 0.46 | 0.6781 | 1.47x | PASS |
| `cerny2007_branch_responsiveness` | Cerny concentration pair | yes | 3.067 | 1.487 | 0.48x | **FAIL** |
| `amrani_hemaimi_onoff_switch` | pyrazines | no | -- | -- | -- | _not scored_ |

## Each row's source anchor and qualification

### `zhou_pH8_MFT`

- source: Zhou 2023 Table 1, pH 8 column: MFT 525.62 ug/L
- HS-SPME, absolute_concentration: false

### `zhou_pH8_FFT`

- source: Zhou 2023 Table 1, pH 8 column: FFT 325.22 ug/L

### `zhou_pH8_MFTD`

- source: Zhou 2023 Table 1, pH 8 column: MFTD 50.07 ug/L
- only 1.5x above a pseudo-LOD: its calibration curve zeroes at x = 34.3 ug/L (inventory sec. F row 14)

### `zhou_pH8_ACTZ`

- source: Zhou 2023 Table 1, pH 8 column: ACTZ 582.34 ug/L

### `zhou_pH8_FUR_arp_alone`

- source: Zhou 2023 Table 1, pH 8, ARP alone: 2-furfural 436.63 ug/L

### `zhou_MFT_shape_pH8_over_pH7`

- source: Zhou 2023 T1: MFT 525.62 at pH 8 against 1588.57 at pH 7 => 0.331. THE MAXIMUM AT pH 7 IS THE POINT: a model fitted at the maximum must PREDICT the fall on the alkaline side.
- the single most informative row in the whole hold-out set

### `zhou_FFT_shape_pH8_over_pH7`

- source: Zhou 2023 T1: FFT 325.22 at pH 8 against 757.965 at pH 7 => 0.429, monotone DOWN -- the acid-lane prediction

### `zhou_pH6_MFT`

- source: Zhou 2023 Table 1, pH 6 column: MFT 696.99 ug/L
- DIAGNOSTIC per FIT_HOLDOUT_DECLARATION.md sec.5 decision 1: the pH labels are INITIAL pH of an unbuffered system whose pH-6 and pH-7 runs converge to within 0.2 units by the end of heating

### `zhou_pH6_FFT`

- source: Zhou 2023 Table 1, pH 6 column: FFT 813.65 ug/L
- DIAGNOSTIC per FIT_HOLDOUT_DECLARATION.md sec.5 decision 1: the pH labels are INITIAL pH of an unbuffered system whose pH-6 and pH-7 runs converge to within 0.2 units by the end of heating

### `zhou_pH6_MFTD`

- source: Zhou 2023 Table 1, pH 6 column: MFTD 59.7 ug/L
- DIAGNOSTIC per FIT_HOLDOUT_DECLARATION.md sec.5 decision 1: the pH labels are INITIAL pH of an unbuffered system whose pH-6 and pH-7 runs converge to within 0.2 units by the end of heating

### `zhou_fig3_grid`

- source: Zhou 2023 Fig. 3, Cys + MGO 20 mM each, unbuffered, 120 C, 0/30/60/90 min, pH 6/7/8
- NOT SCORED, and the reason is a SCOPE limit rather than a model failure: seven of the grid's eight panels are pyrazines, thiophenes and methylthiazoles, all of which this wave declares OUT OF SCOPE (sulfur.OUT_OF_SCOPE). The one in-scope quantity the grid supplies is the 2-acetylthiazole cross-system pair, which IS scored below. Declaring this unscored is the honest option; scoring one panel of eight and calling the grid passed would not be.

### `zhou_582_vs_665_consistency`

- source: Zhou 2023: 2-acetylthiazole at 60 min, pH 8 -- ARP-fed 582 ug/L against MGO-fed 665 ug/L (14% apart)
- scores the ARP -> MGO flux with no free downstream parameter

### `whitfield_2001_pH65_collapse`

- source: Whitfield 2001 Table 1: fed NF + cysteine, MFT 0.150 mol% at pH 4.5 -> nd (<0.0010 mol%) at pH 6.5 = a >=150x collapse within one lab, while total volatiles fall only 1.7x
- ONE-SIDED: the measurement is a LOWER bound (the pH-6.5 value is a non-detect), so any predicted collapse of 150x or more passes. NOTE the inventory's own qualification: Whitfield 2001's H2S column should NOT carry a standalone row (its Methods omits HMF -- likely a typo).

### `cerny2003_NF_share_ceiling_95C`

- source: Cerny 2003 Table 3: the norfuraneol share of the MFT flux is <=7% at 95 C / 4 h / pH 5.00, and that is an UPPER bound (NF was spiked in at 1.5x the cysteine)
- ONE-SIDED CEILING, evaluated at CERNY'S conditions (95 C, 4 h) and not at the fit panel's 145 C, exactly as the declaration requires

### `cerny2003_intact_skeleton_share_95C`

- source: Cerny 2003 Table 2: MFT is 49% unlabelled / 46% 13C5 with NO fragment pattern => ~93% of MFT comes from the INTACT-C5 (1,4-dideoxyosone) route at 95 C; 'pathways via ribose fragmentation were not relevant'
- and its partner: Cerny 2004 finds the C2+C3 route 'not relevant' at 95 C while Hofmann 1998 T10 makes it the STRONGEST route at 145 C. The model must therefore get a TEMPERATURE-DEPENDENT route mix right, with one lumped formation Ea.

### `route_mix_moves_with_temperature`

- source: Hofmann 1998 T10 vs Cerny 2004: the C2+C3 route is the STRONGEST MFT route at 145 C (0.24 mol%) and 'not relevant' at 95 C
- DIRECTIONAL: passes if the C2+C3 share RISES with temperature

### `hofmann_dry_180C`

- source: Hofmann 1998 T2 dry rows: ribose/glucose/rhamnose on silica, 180 C / 6 min, FFT 97.2 / 1.4 / 0.4 and MFT 25.1 / 4.2 / 3.1 ug; and T10's C2+C3 dry row at 1553.9 ug (1.39 mol%)
- UNSCOREABLE, PRE-REGISTERED. The module has NO WATER-ACTIVITY TERM, because nothing in its fit corpus varies a_w -- B1's policy 3, inherited. A dry silica system at 180 C is outside the declared domain, and producing a number for it would be inventing an a_w dependence at scoring time. Also note inventory sec. F row 17: the dry protocol is described TWICE with different times (5 vs 6 min) and a 4x different buffer molarity while reporting the same number to four significant figures.

### `hofmann_ribose_pH3_MFT`

- source: Hofmann 1998 T2, ribose pH 3: MFT 553.0 ppb
- PRE-REGISTERED FAILURE. Inventory sec. B2.5: Hofmann's buffered free-sugar system has MFT FALLING as pH rises while Zhou's unbuffered Amadori-fed system has it RISING to a pH-7 maximum, with the absolute levels 64x apart at pH 7. The dossier's own recommendation is NOT to merge them into one pH response curve. A single structural pH shape can satisfy at most one, and this module carries Zhou's.

### `hofmann_ribose_pH3_FFT`

- source: Hofmann 1998 T2, ribose pH 3: FFT 229.0 ppb
- FFT is the row the acid lane SHOULD get right: it falls monotonically with pH in all three papers that measure it.

### `hofmann_ribose_pH7_MFT`

- source: Hofmann 1998 T2, ribose pH 7: MFT 25.0 ppb
- PRE-REGISTERED FAILURE. Inventory sec. B2.5: Hofmann's buffered free-sugar system has MFT FALLING as pH rises while Zhou's unbuffered Amadori-fed system has it RISING to a pH-7 maximum, with the absolute levels 64x apart at pH 7. The dossier's own recommendation is NOT to merge them into one pH response curve. A single structural pH shape can satisfy at most one, and this module carries Zhou's.

### `hofmann_ribose_pH7_FFT`

- source: Hofmann 1998 T2, ribose pH 7: FFT 12.0 ppb
- FFT is the row the acid lane SHOULD get right: it falls monotonically with pH in all three papers that measure it.

### `meynier_MFT_pH_fold`

- source: Meynier 1995: MFT falls >152x over pH 4.5 -> 6.5 (constant-pH, 4-amino-acid, single-sugar factorial, ~59 series; DIRECTIONAL ONLY -- 'response factors not determined')
- DIRECTIONAL, order-of-magnitude band. Meynier has no absolute quantification of any kind (inventory sec. C.8).

### `meynier_FFT_pH_fold`

- source: Meynier 1995: FFT falls 6.1x (constant-pH, 4-amino-acid, single-sugar factorial, ~59 series; DIRECTIONAL ONLY -- 'response factors not determined')
- DIRECTIONAL, order-of-magnitude band. Meynier has no absolute quantification of any kind (inventory sec. C.8).

### `meynier_FUR_pH_fold`

- source: Meynier 1995: furfural falls 15-49x (midpoint 25x) (constant-pH, 4-amino-acid, single-sugar factorial, ~59 series; DIRECTIONAL ONLY -- 'response factors not determined')
- DIRECTIONAL, order-of-magnitude band. Meynier has no absolute quantification of any kind (inventory sec. C.8).

### `hofmann2002_brew_80C_FFT`

- source: Hofmann 2002 Fig. 1: FFT loss in a real coffee brew (50 g powder/L, thermos, 80 C) at 0.023 /min, t1/2 ~30 min
- THE ONLY TEMPERATURE EXTRAPOLATION THE MODULE HAS, and the model has no activation energy for this channel by policy, so the constant is HELD at its 30 C value. A pass would mean the depletable-electrophile structure is doing the work; a failure localises the electrophile-pool depletion term, which is exactly what the declaration says this row is for.

### `vanseeventer_50C_MFT_per_day`

- source: van Seeventer 2001 Table 1: MFT 59% of initial per DAY, ZERO ORDER, 50 C, pH 5.0, air ~ argon (FFT 28%/day)
- PRE-REGISTERED FAILURE. The C-5 oligomerisation channel that carries this loss has NO rate on any FIT row -- its only measurement IS this hold-out -- so the model carries it at ZERO and predicts essentially no loss. The failure localises a MISSING CHANNEL, not a wrong barrier, and 'do not bolt van Seeventer's 59%/day onto the network as the MFT sink' is a verbatim dossier instruction (sec. C.17).

### `zhang_115C_MFT_consumed_share`

- source: Zhang 2024 Figs. 2d/e/f: 3-19% of free MFT already consumed into MEASURED products in the Cys arm (21-55% in the GSH arm); LOWER BOUNDS, since a third sink (melanoidin/quinone) is named and never measured
- SCORED AGAINST THE BAND, not the midpoint, because the measurement is a lower bound and a range. SYSTEM MISMATCH DISCLOSED: Zhang's Fig. 2 system is 1:3:3 Met:VB1:Xyl and METHIONINE IS NOT IN THIS MODULE'S STATE VECTOR, so the MeSH supply is structurally wrong; the Fig. 1 geometry is substituted and that substitution is part of what is being scored.

### `zhou_120C_dimer_share_pH6`

- source: Zhou 2023 Table 1: the dimer carries 8.6% of MFT-equivalents at pH6, and the fraction is NEAR-INVARIANT in pH while [MFT] swings 3.0x
- THE SCORE TO BEAT: Zhou (120 C, unbuffered, Amadori-fed, pH 6-8) and Zhang (115 C, buffered pH 4.9, thiamine/xylose) agree to 1.3x on this channel across two labs, two feedstocks and two buffers. Both are held out, so that agreement stays an out-of-sample fact.

### `zhou_120C_dimer_share_pH7`

- source: Zhou 2023 Table 1: the dimer carries 6.5% of MFT-equivalents at pH7, and the fraction is NEAR-INVARIANT in pH while [MFT] swings 3.0x
- THE SCORE TO BEAT: Zhou (120 C, unbuffered, Amadori-fed, pH 6-8) and Zhang (115 C, buffered pH 4.9, thiamine/xylose) agree to 1.3x on this channel across two labs, two feedstocks and two buffers. Both are held out, so that agreement stays an out-of-sample fact.

### `zhou_120C_dimer_share_pH8`

- source: Zhou 2023 Table 1: the dimer carries 9.6% of MFT-equivalents at pH8, and the fraction is NEAR-INVARIANT in pH while [MFT] swings 3.0x
- THE SCORE TO BEAT: Zhou (120 C, unbuffered, Amadori-fed, pH 6-8) and Zhang (115 C, buffered pH 4.9, thiamine/xylose) agree to 1.3x on this channel across two labs, two feedstocks and two buffers. Both are held out, so that agreement stays an out-of-sample fact.

### `zhou_dimer_share_pH_invariance`

- source: Zhou 2023: the dimer share is pH-INVARIANT (8.6 / 6.5 / 9.6% over pH 6-8) while [MFT] swings 3.0x -- a 1.48x spread
- a SHAPE test that a fitted pH term would have been able to fake and a structural one cannot

### `cerny2007b_control_no_cysteine`

- source: Cerny 2007b Table 5: with NO CYSTEINE, MFT is >99% thiamine-derived and <1% xylose-derived
- THE SHARPEST STRUCTURAL TEST IN THE WHOLE SPLIT. A model fitted on the ternary must PREDICT both limiting cases; getting 54:46 right while missing either control means the routes are wrong and the ratio was fitted.

### `cerny2007b_control_no_thiamine`

- source: Cerny 2007b Table 6: with NO THIAMINE, MFT is >95% xylose-derived and <5% thiamine-derived
- the other limiting case

### `cerny2007_1x_xylose_share`

- source: Cerny 2007 Table 5 arm B, 1x (xylose 0.15 / cys 0.05 / thiamine 0.05 M): 85 : 15 thiamine : xylose, i.e. a 15% xylose share

### `cerny2007_2x_xylose_share`

- source: Cerny 2007 Table 5 arm A, 2x (xylose 0.3 / cys 0.1 / thiamine 0.1 M): 54 : 46, i.e. a 46% xylose share

### `cerny2007_branch_responsiveness`

- source: Cerny 2007 Table 5: a 2x change in precursor loading moves the xylose share of MFT from 15% to 46% -- a 3.1x change in the BRANCH FRACTION, one lab, one method, one pH, one temperature
- THE SINGLE HIGHEST-VALUE HOLD-OUT ROW IN THE SET. It scores whether the model's branch fractions RESPOND TO CONCENTRATION AT ALL. A model with fixed branch fractions predicts exactly 1.00x here and fails by construction, which is the point of the row.

### `amrani_hemaimi_onoff_switch`

- source: Amrani-Hemaimi 1995: 3-ethyl-2,5-dimethylpyrazine present with alanine (20 / 19%) and ABSENT with glycine (0 / 0%), 100% 13C-labelled => 'one single reaction route exists'
- DEFERRED EXPLICITLY. Pyrazines are OUT OF SCOPE for this wave (sulfur.OUT_OF_SCOPE): no pyrazine is in the state vector, so neither this hold-out nor the Amrani-Hemaimi Table 2 FIT rows are touched. A later wave that adds the pyrazine lane inherits both, unscored and unfitted.

## The aroma qualification

Zhou 2023 SI Table S2 prints, at pH 7, MFT OAV 3.18e5 against disulfide 3.21e5 -- the dimer's OAV MATCHES the monomer's, because the dimer is 15.6x more potent and carries 6.5% of the thiol equivalents.

Any objective that scores the dimerisation channel as a pure aroma LOSS is wrong by roughly the threshold ratio. This module tracks the dimer as a species with its own potency, so the hold-out scores above are mass scores and the OAV comparison is reported alongside them rather than folded into them.

## Firewall disclosure

The K3 inventory sec. A.3.3 prints Zhou's pH-6 and pH-8 columns next to the pH-7 FIT column, and Zhang2024_extraction.md sec. 3 prints Fig. 2 panels (d) and (e). The build brief directed this wave to read both files, so those hold-out values were SEEN before the fit ran. They were not used: they appear in no fit row, bound or initialisation, and a firewall test greps the fit script and src/kinetic_core/ for the literals. The frozen mp_holdout_* bundles were never opened.
