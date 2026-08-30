# Kinetic core, Build Wave B2.2 -- the pre-registered hold-out scorecard

Predicted from the FROZEN fit in `results/validation/kinetic_core_b2_2_fit_report.json`. **No parameter was changed after the fit**; there is no optimiser in the scoring script.

**GATING: 12 / 33 passed.** Diagnostic: 6 / 9. Unscoreable or explicitly deferred: 3.

**Against B2.1: 15 / 33 gating.** B2.2 adds no row and removes none, so the denominators match and the totals ARE directly comparable. **6 FIXED, 9 REGRESSED.**

## B2.1 versus B2.2, row by row

| row | gating | B2.1 fold | B2.1 | B2.2 fold | B2.2 | |
|---|---|---:|---|---:|---|---|
| `zhou_pH8_MFT` | yes | 0.628x | PASS | 0.343x | PASS | unchanged |
| `zhou_pH8_FFT` | yes | 0.00312x | **FAIL** | 0.381x | PASS | FIXED |
| `zhou_pH8_MFTD` | yes | 0.701x | PASS | 4.33x | **FAIL** | REGRESSED |
| `zhou_pH8_ACTZ` | yes | 0.000719x | **FAIL** | 0.00028x | **FAIL** | unchanged |
| `zhou_pH8_FUR_arp_alone` | yes | 0.000744x | **FAIL** | 7.46x | **FAIL** | unchanged |
| `zhou_MFT_shape_pH8_over_pH7` | yes | 1.28x | PASS | 0.509x | PASS | unchanged |
| `zhou_FFT_shape_pH8_over_pH7` | yes | 0.00524x | **FAIL** | 0.779x | PASS | FIXED |
| `zhou_pH6_MFT` | no | 0.0576x | **FAIL** | 1.82x | PASS | FIXED |
| `zhou_pH6_FFT` | no | 2.15x | PASS | 0.402x | PASS | unchanged |
| `zhou_pH6_MFTD` | no | 0.000259x | **FAIL** | 0.238x | **FAIL** | unchanged |
| `zhou_fig3_grid` | no | -- | -- | -- | -- | unchanged |
| `zhou_582_vs_665_consistency` | yes | 0.00032x | **FAIL** | 9.84e-07x | **FAIL** | unchanged |
| `whitfield_2001_pH65_collapse` | yes | 0.00408x | **FAIL** | 0.00437x | **FAIL** | unchanged |
| `cerny2003_NF_share_ceiling_95C` | yes | 0.736x | PASS | 14.1x | **FAIL** | REGRESSED |
| `cerny2003_intact_skeleton_share_95C` | yes | 1.02x | PASS | 0.0147x | **FAIL** | REGRESSED |
| `route_mix_moves_with_temperature` | no | -- | PASS | -- | **FAIL** | REGRESSED |
| `hofmann_dry_180C` | yes | -- | -- | -- | -- | unchanged |
| `hofmann_ribose_pH3_MFT` | yes | 0.0589x | **FAIL** | 0.0468x | **FAIL** | unchanged |
| `hofmann_ribose_pH3_FFT` | yes | 24.8x | **FAIL** | 1.47x | PASS | FIXED |
| `hofmann_ribose_pH7_MFT` | yes | 1.13x | PASS | 0.00539x | **FAIL** | REGRESSED |
| `hofmann_ribose_pH7_FFT` | yes | 0.379x | PASS | 0.208x | **FAIL** | REGRESSED |
| `meynier_MFT_pH_fold` | no | 0.0159x | **FAIL** | 3.87x | PASS | FIXED |
| `meynier_FFT_pH_fold` | no | 5.31x | PASS | 12.2x | **FAIL** | REGRESSED |
| `meynier_FUR_pH_fold` | no | 0.21x | PASS | 0.416x | PASS | unchanged |
| `hofmann2002_brew_80C_FFT` | yes | 2.88x | PASS | 2.96x | PASS | unchanged |
| `vanseeventer_50C_MFT_per_day` | yes | 1.38x | PASS | 1.66x | PASS | unchanged |
| `zhang_115C_MFT_consumed_share` | yes | 0.473x | PASS | 0.981x | PASS | unchanged |
| `zhou_120C_dimer_share_pH6` | yes | 0.00451x | **FAIL** | 0.131x | **FAIL** | unchanged |
| `zhou_120C_dimer_share_pH7` | yes | 0.747x | PASS | 0.233x | **FAIL** | REGRESSED |
| `zhou_120C_dimer_share_pH8` | yes | 1.12x | PASS | 12.6x | **FAIL** | REGRESSED |
| `zhou_dimer_share_pH_invariance` | yes | 187x | **FAIL** | 72.4x | **FAIL** | unchanged |
| `cerny2007b_control_no_cysteine` | yes | 1.01x | PASS | 1.01x | PASS | unchanged |
| `cerny2007b_control_no_thiamine` | yes | 1.05x | PASS | 1.05x | PASS | unchanged |
| `cerny2007_1x_xylose_share` | yes | 2.74x | **FAIL** | 3.69x | **FAIL** | unchanged |
| `cerny2007_2x_xylose_share` | yes | 1.01x | PASS | 1.19x | PASS | unchanged |
| `cerny2007_branch_responsiveness` | yes | 0.369x | **FAIL** | 0.323x | **FAIL** | unchanged |
| `kang_140C_MFT` | yes | 0.000656x | **FAIL** | 2.15x | PASS | FIXED |
| `kang_switch_on_MFT` | yes | 0.000342x | **FAIL** | 2.45x | **FAIL** | unchanged |
| `kang_140C_FFT` | yes | 1.52e-07x | **FAIL** | 0.0706x | **FAIL** | unchanged |
| `kang_switch_on_FFT` | yes | 3.21e-07x | **FAIL** | 0.0792x | **FAIL** | unchanged |
| `kang_thiols_diverge_from_their_class` | no | -- | PASS | -- | PASS | unchanged |
| `kang_140C_cys_conversion` | no | 1.38x | PASS | 1.24x | PASS | unchanged |
| `sun2019_pH9_over_pH3_free_FFT` | yes | 0.0735x | **FAIL** | 0.0569x | **FAIL** | unchanged |
| `sun2019_temperature_ordering_inversion` | yes | -- | **FAIL** | -- | **FAIL** | unchanged |
| `amrani_hemaimi_onoff_switch` | no | -- | -- | -- | -- | unchanged |

A row marked REGRESSED is not automatically bad and a row marked FIXED is not automatically good: see the qualifications block. B2.2 adds no hold-out row and removes none, so NEW ROW should never appear -- if it does, the two scorers have drifted apart and the comparison is unsound.

## Pre-registered failures (written down before the numbers)

- van Seeventer 50 C -- STILL, and now expected to fail from the OTHER SIDE. The C-5 oligomerisation channel that actually carries this loss still has a rate of exactly zero, because its only measurement IS this hold-out. B2 landed inside the 3x band by over-destroying the thiol through unassigned decay lumps pinned to their 145 C values, and its own report refused to count that pass. B2.1 gives those lumps an activation energy, so they no longer run at 145 C speed in a 50 C store -- which is correct, and which should turn an accidental pass into an honest under-prediction. A scorecard that goes DOWN on this row is the model getting more truthful, not less.
- Hofmann dry-180 C -- unscoreable. Still no water-activity term, because nothing in the fit corpus varies a_w.
- Hofmann pH-3/pH-7 MFT -- the B2.5 sign conflict is a property of the CORPUS, not of the model: Hofmann's buffered free-sugar system has MFT falling with pH while Zhou's unbuffered Amadori-fed system has it peaking at pH 7, with the levels 64x apart at pH 7, and the dossier's own instruction is not to merge them. One pH shape satisfies at most one. Nothing in B2.1 changes that.
- Zhang Fig. 2 -- methionine is still not in the state vector.
- Cerny 2003's norfuraneol-share ceiling at 95 C -- B2 was 5.47x over an upper bound and nothing in B2.1 targets the norfuraneol partition. Expected to stay a fail.
- Sun 2019's pH-9 TEMPERATURE-ORDERING INVERSION -- expected to FAIL, and the reason is worth stating in advance. At pH 9 the thiolate channel dominates and, like every thermal channel in the module, it accelerates with temperature, so free thiol should FALL with T. The measurement has it RISING (0.4 -> 0.5 -> 0.5 at 20 -> 55 -> 90 C). The only term in the module that moves the other way is Stack's measured release, which grows with temperature; whether it can outrun the sinks at pH 9 is exactly what the row asks, and the honest prior is that it cannot.
- Kang 140 C -- DIRECTION UNKNOWN, WHICH IS WHY IT IS THE EXAM. Amendment 5 pre-registers a 3.8x/2.5x UNDER-prediction for a single-Ea branch. This module is not single-Ea at that node: the sulfide supply is gated behind Zheng's measured 133 kJ/mol while the carbon skeleton carries a much smaller lumped barrier, so the thiol lane accelerates late by construction. That mechanism was in the network before Kang was read and it is not tunable -- but it can equally OVERSHOOT. This wave declares that it does not know which side it will miss on, and that only a landing inside 3x counts as the switch-on having emerged.

## Qualifications the pass/fail column hides

- **`vanseeventer_50C_MFT_per_day` (PASS)** -- THIS PASS SHOULD NOT BE COUNTED AS ONE. It was pre-registered as a failure and it passed the 3x band from the WRONG SIDE, by over-predicting: the model destroys ~99% of the MFT per day against a measured 59%, and it does so through the unassigned-decay and dimerisation channels, NOT through the C-5 oligomerisation channel that van Seeventer actually measured -- which still carries a rate of exactly zero. The FUNCTIONAL FORM is also still wrong: the measured loss is ZERO ORDER and every channel here is first or second order. The row tests the functional form, and on that test the model fails while landing inside the magnitude band by coincidence.

- **`cerny2007_branch_responsiveness` (FAIL)** -- A FAIL THAT CARRIES THE WAVE'S MAIN ARCHITECTURAL RESULT. This is the row a model with fixed branch fractions fails BY CONSTRUCTION, because such a model predicts exactly 1.00x. This model predicts 0.99x against a measured 3.07x. So the branch fraction DOES respond to concentration -- the architectural requirement is met and the no-fixed-fraction design is vindicated -- but it responds about half as strongly as measured. The residual is a magnitude error in the relative reaction ORDERS of the two routes, not the structural error the row was designed to catch. Reported as a fail, because it is one.

- **`zhou_pH8_FFT` (PASS)** -- THE ROW B2.1 EXISTS TO FIX. B2 predicted 1.23e-07 against a measured 2.85e-03 -- 0.00x -- and its report named the cause: a structural pH slope of one decade per pH unit, right in shape and far too steep in level. B2.1 changed three things and NONE of them is a fitted pH parameter: the pKa are now evaluated at reaction temperature; no pH-sensitive product has a single fully-catalysed route any more; and the thiol's pH-dependent CONSUMPTION lane, which Kumazawa measures directly at 121 C and which B2 did not have at all, now carries the part of the pH response that belongs to it. B2.1 predicts 0.00109 against the same measured 2.849e-03, a fold of 0.381x. Whether that is a pass or a fail, it is the single number by which this wave should be judged.

- **`hofmann2002_brew_80C_FFT` (PASS)** -- THE FAILURE THE DECLARATION PREDICTED AND ASKED FOR, and B2.1 attacks it from two measured directions rather than by tuning. FIRST: B2 pinned every fitted decay lump to its 145 C value and evaluated it unchanged at 80 C, so k_fft_decay alone predicted several times the measured brew loss; those lumps now carry the lumped Ea, which is a WEAKER claim than temperature-independence, not a stronger one. SECOND: Stack 2018's measured van 't Hoff enthalpy is NEGATIVE, so the thiol-electrophile equilibrium constant FALLS 5x between 30 and 80 C -- heating releases bound thiol. B2 could not express that, because its K was a temperature-independent number that the source does not print. B2 lost FFT 18x too fast; B2.1 is at 3x.

- **`kang_140C_MFT` (PASS)** -- THE NEW EXAM, AND THE ONE WHOSE DIRECTION OF FAILURE MATTERS MOST. Amendment 5 pre-registers that a single-Ea branch UNDER-predicts here by 3.8x, because the measured behaviour is non-Arrhenius: free-thiol formation switches on between 120 and 140 C. This module is not single-Ea at that node -- its thiol-forming flux is [carbon skeleton] x [sulfide] and the sulfide is gated behind Zheng & Ho's MEASURED 133 kJ/mol, the largest barrier in the module by a factor of two, against a lumped formation barrier for everything else. So it can also OVERSHOOT, and an overshoot is a failure of the same seriousness as an undershoot: it would mean the switch-on mechanism is real but mis-sized. MEASURED 5.907 ug/L; PREDICTED 0.0001113 mmol/L, fold 2.15x.

## Row by row -- no averaging, no dropping

| row | group | gating | observed | predicted | fold | pass |
|---|---|---|---:|---:|---:|---|
| `zhou_pH8_MFT` | Zhou pH-8 column | yes | 0.004604 | 0.001581 | 0.34x | PASS |
| `zhou_pH8_FFT` | Zhou pH-8 column | yes | 0.002849 | 0.001086 | 0.38x | PASS |
| `zhou_pH8_MFTD` | Zhou pH-8 column | yes | 0.0002212 | 0.0009579 | 4.33x | **FAIL** |
| `zhou_pH8_ACTZ` | Zhou pH-8 column | yes | 0.00458 | 1.284e-06 | 0.00x | **FAIL** |
| `zhou_pH8_FUR_arp_alone` | Zhou pH-8 column | yes | 0.004544 | 0.0339 | 7.46x | **FAIL** |
| `zhou_MFT_shape_pH8_over_pH7` | Zhou pH-8 column | yes | 0.3309 | 0.1684 | 0.51x | PASS |
| `zhou_FFT_shape_pH8_over_pH7` | Zhou pH-8 column | yes | 0.4291 | 0.3344 | 0.78x | PASS |
| `zhou_pH6_MFT` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.006105 | 0.01109 | 1.82x | PASS |
| `zhou_pH6_FFT` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.007127 | 0.002867 | 0.40x | PASS |
| `zhou_pH6_MFTD` | Zhou pH-6 column (DIAGNOSTIC) | no | 0.0002638 | 6.267e-05 | 0.24x | **FAIL** |
| `zhou_fig3_grid` | Zhou Fig. 3 grid | no | -- | -- | -- | _not scored_ |
| `zhou_582_vs_665_consistency` | Zhou cross-system pair | yes | 0.8757 | 8.614e-07 | 0.00x | **FAIL** |
| `whitfield_2001_pH65_collapse` | Whitfield 2001 pH collapse | yes | 150 | 0.6558 | 0.00x | **FAIL** |
| `cerny2003_NF_share_ceiling_95C` | Cerny 2003, 95 C | yes | 0.07 | 0.9863 | 14.09x | **FAIL** |
| `cerny2003_intact_skeleton_share_95C` | Cerny 2003, 95 C | yes | 0.93 | 0.0137 | 0.01x | **FAIL** |
| `route_mix_moves_with_temperature` | Cerny 2003, 95 C | no | -- | -- | -- | **FAIL** |
| `hofmann_dry_180C` | Hofmann dry-180 | yes | -- | -- | -- | _not scored_ |
| `hofmann_ribose_pH3_MFT` | Hofmann pH axis | yes | 0.004844 | 0.0002268 | 0.05x | **FAIL** |
| `hofmann_ribose_pH3_FFT` | Hofmann pH axis | yes | 0.002006 | 0.002942 | 1.47x | PASS |
| `hofmann_ribose_pH7_MFT` | Hofmann pH axis | yes | 0.000219 | 1.18e-06 | 0.01x | **FAIL** |
| `hofmann_ribose_pH7_FFT` | Hofmann pH axis | yes | 0.0001051 | 2.182e-05 | 0.21x | **FAIL** |
| `meynier_MFT_pH_fold` | Meynier directional | no | 152 | 588.4 | 3.87x | PASS |
| `meynier_FFT_pH_fold` | Meynier directional | no | 6.1 | 74.61 | 12.23x | **FAIL** |
| `meynier_FUR_pH_fold` | Meynier directional | no | 25 | 10.41 | 0.42x | PASS |
| `hofmann2002_brew_80C_FFT` | thiol sinks | yes | 0.023 | 0.06809 | 2.96x | PASS |
| `vanseeventer_50C_MFT_per_day` | thiol sinks | yes | 0.59 | 0.9793 | 1.66x | PASS |
| `zhang_115C_MFT_consumed_share` | thiol sinks | yes | 0.11 | 0.1079 | 0.98x | PASS |
| `zhou_120C_dimer_share_pH6` | thiol sinks | yes | 0.086 | 0.01131 | 0.13x | **FAIL** |
| `zhou_120C_dimer_share_pH7` | thiol sinks | yes | 0.065 | 0.01513 | 0.23x | **FAIL** |
| `zhou_120C_dimer_share_pH8` | thiol sinks | yes | 0.096 | 1.212 | 12.62x | **FAIL** |
| `zhou_dimer_share_pH_invariance` | thiol sinks | yes | 1.48 | 107.2 | 72.43x | **FAIL** |
| `cerny2007b_control_no_cysteine` | Cerny single-route controls | yes | 0.99 | 1 | 1.01x | PASS |
| `cerny2007b_control_no_thiamine` | Cerny single-route controls | yes | 0.95 | 1 | 1.05x | PASS |
| `cerny2007_1x_xylose_share` | Cerny concentration pair | yes | 0.15 | 0.5529 | 3.69x | **FAIL** |
| `cerny2007_2x_xylose_share` | Cerny concentration pair | yes | 0.46 | 0.5479 | 1.19x | PASS |
| `cerny2007_branch_responsiveness` | Cerny concentration pair | yes | 3.067 | 0.991 | 0.32x | **FAIL** |
| `kang_140C_MFT` | Kang 140 C (NEW gating hold-out) | yes | 5.174e-05 | 0.0001113 | 2.15x | PASS |
| `kang_switch_on_MFT` | Kang 140 C (NEW gating hold-out) | yes | 4.256 | 10.43 | 2.45x | **FAIL** |
| `kang_140C_FFT` | Kang 140 C (NEW gating hold-out) | yes | 0.0001002 | 7.078e-06 | 0.07x | **FAIL** |
| `kang_switch_on_FFT` | Kang 140 C (NEW gating hold-out) | yes | 2.785 | 0.2205 | 0.08x | **FAIL** |
| `kang_thiols_diverge_from_their_class` | Kang 140 C (NEW gating hold-out) | no | -- | -- | -- | PASS |
| `kang_140C_cys_conversion` | Kang 140 C (NEW gating hold-out) | no | 0.626 | 0.775 | 1.24x | PASS |
| `sun2019_pH9_over_pH3_free_FFT` | Sun 2019 pH-9 column | yes | 0.1481 | 0.008423 | 0.06x | **FAIL** |
| `sun2019_temperature_ordering_inversion` | Sun 2019 pH-9 column | yes | -- | -- | -- | **FAIL** |
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

### `kang_140C_MFT`

- source: Kang 2026 SI Table S4, 140 C column: MFT 5.907 ug/L (Tier A, external calibration curve printed, raster-verified mu-g/L, nine-subtotal arithmetic closure to +/-0.003)
- A TRUE EXTRAPOLATION: the fit saw 100 and 120 C in this pot and nothing above. The declaration pre-registers a 3.8x (MFT) / 2.5x (FFT) UNDER-prediction for any single-Ea branch. Uncertainty on the observed value is the dossier's replacement +/-15% (sec. 7d), NOT the printed SD, which is demonstrably unreproducible between Tables S4 and S5 for the identical experiment.

### `kang_switch_on_MFT`

- source: Kang 2026 SI sec. 7a: MFT rises 4.26x from 120 to 140 C after rising only 1.12x from 100 to 120 C -- apparent Ea climbing from ~6-7 to 70-98 kJ/mol
- THE ROW THE WHOLE REVISION IS ABOUT. It scores the SWITCH-ON itself rather than the level, so it is immune to the calibration axis and to any offset the fit absorbed on the two lower rungs. A single-Ea branch cannot produce a rising apparent Ea at all.

### `kang_140C_FFT`

- source: Kang 2026 SI Table S4, 140 C column: FFT 11.439 ug/L (Tier A, external calibration curve printed, raster-verified mu-g/L, nine-subtotal arithmetic closure to +/-0.003)
- A TRUE EXTRAPOLATION: the fit saw 100 and 120 C in this pot and nothing above. The declaration pre-registers a 3.8x (MFT) / 2.5x (FFT) UNDER-prediction for any single-Ea branch. Uncertainty on the observed value is the dossier's replacement +/-15% (sec. 7d), NOT the printed SD, which is demonstrably unreproducible between Tables S4 and S5 for the identical experiment.

### `kang_switch_on_FFT`

- source: Kang 2026 SI sec. 7a: FFT rises 2.79x from 120 to 140 C after rising only 1.10x from 100 to 120 C -- apparent Ea climbing from ~6-7 to 70-98 kJ/mol
- THE ROW THE WHOLE REVISION IS ABOUT. It scores the SWITCH-ON itself rather than the level, so it is immune to the calibration axis and to any offset the fit absorbed on the two lower rungs. A single-Ea branch cannot produce a rising apparent Ea at all.

### `kang_thiols_diverge_from_their_class`

- source: Kang 2026 SI sec. 7a: over 120 -> 140 C the sulfur CLASS decelerates (2.57x then 1.68x, apparent Ea falling 57.5 -> 35.2 kJ/mol) while MFT and FFT accelerate. Precursor depletion DEPRESSES apparent Ea, so a saturation artefact cannot produce the thiols' behaviour.
- DIRECTIONAL, and DIAGNOSTIC rather than gating because this module carries only part of Kang's sulfur class -- no thiophenes, no thiazoles beyond 2-acetylthiazole, all declared out of scope. The in-scope proxy for 'the thiols outrun the rest of the pot' is that MFT's 120->140 fold EXCEEDS furfural's, which is a comparison the module can honestly make.

### `kang_140C_cys_conversion`

- source: Kang 2026 SI Fig. S4 (digitised): free-cysteine conversion 62.6% at 140 C / 120 min, against 16.2% and 38.7% at the two fitted rungs
- DIAGNOSTIC, not gating: it is figure-digitised, its n is unstated, and the identification of the system as FREE cysteine rather than the TTCA-bound moiety is 85% confident on the dossier's own reading. The row it really tests is Kang's measured Ea of 55.1 kJ/mol, which this module carries as a fixed barrier on r_cys_thermal.

### `sun2019_pH9_over_pH3_free_FFT`

- source: Sun 2019 Fig. 4A, values printed in the text (p. 453): free 2-FFT in a coffee brew after 1 h at 20 C is 2.7 ug/L at pH 3 and 0.4 ug/L at pH 9, a 6.75x collapse
- RATIO, so the brew's free-thiol loading divides out. B2 had no pH-dependent thiol sink at all and could only express this through formation, which does not run in a spiked brew; B2.1 has the thiolate-mediated channel Kumazawa's FIT grid pays for.

### `sun2019_temperature_ordering_inversion`

- source: Sun 2019 p. 453: free 2-FFT falls with temperature at pH 3 (2.7 > 1.9 > 1.6 at 20/55/90 C) and RISES with temperature at pH 9 (0.4 < 0.5 = 0.5). The ordering INVERTS across the pH axis.
- THE DOSSIER'S OWN VERDICT: 'Any model that gets that inversion right has learned something real; getting it wrong is a clean falsification.' Note the pH-9 differences (0.4 vs 0.5) are one significant figure apart and the paper prints no error bars, so this is scored as an ORDERING and not as a magnitude.

### `amrani_hemaimi_onoff_switch`

- source: Amrani-Hemaimi 1995: 3-ethyl-2,5-dimethylpyrazine present with alanine (20 / 19%) and ABSENT with glycine (0 / 0%), 100% 13C-labelled => 'one single reaction route exists'
- DEFERRED EXPLICITLY. Pyrazines are OUT OF SCOPE for this wave (sulfur.OUT_OF_SCOPE): no pyrazine is in the state vector, so neither this hold-out nor the Amrani-Hemaimi Table 2 FIT rows are touched. A later wave that adds the pyrazine lane inherits both, unscored and unfitted.

## The aroma qualification

Zhou 2023 SI Table S2 prints, at pH 7, MFT OAV 3.18e5 against disulfide 3.21e5 -- the dimer's OAV MATCHES the monomer's, because the dimer is 15.6x more potent and carries 6.5% of the thiol equivalents.

Any objective that scores the dimerisation channel as a pure aroma LOSS is wrong by roughly the threshold ratio. This module tracks the dimer as a species with its own potency, so the hold-out scores above are mass scores and the OAV comparison is reported alongside them rather than folded into them.

## Firewall disclosure

THREE EXPOSURES, ALL DISCLOSED. (1) The K3 inventory sec. A.3.3 prints Zhou's pH-6 and pH-8 columns next to the pH-7 FIT column, and Zhang2024_extraction.md sec. 3 prints Fig. 2 panels (d) and (e) -- B2's exposure, inherited. (2) NEW AND UNAVOIDABLE: kang2026_SI_extraction.md PRINTS THE 140 C COLUMN, and the build brief directed this wave to read that dossier, which is also where the 100 and 120 C FIT rows live. The 140 C values were therefore SEEN before the fit ran. (3) sun2019_extraction.md prints the pH-9 column alongside the pH-3 column. WHAT WAS DONE ABOUT IT: tests/unit/test_kinetic_core_b2_1.py greps every file under src/kinetic_core/ and the fit script for the literals 5.907, 11.439, 60.400, 62.6, 8.195 and Sun's pH-9 column, and fails if any appears outside a line that explicitly marks it as held out; a second test walks the fit script's SYSTEMS dictionary and asserts that no integrated condition is a hold-out condition (no Kang system at 140 C, no Zhou system off pH 7, no system at pH 9 or pH 6.5). The frozen mp_holdout_* bundles under data/benchmarks/external_validation/ were never opened, and a third test asserts that no B2.1 file performs file I/O on that path.
