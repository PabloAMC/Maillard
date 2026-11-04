# M1.0: the shipped engine (R0 = B60) on M1's decisive data, no fit

*Plan step M1.0 (docs/guides/STATUS_AND_PLAN.md sec. 5, W2, row M1); case `docs/reference/cases/M1_thiol_loss_oxidant_ph_drift.md` secs. 3 - 5; CF-107. Configuration: the shipped one, asserted before and after (digest `40b746476905786c`). B60 vector from `results/validation/kinetic_core_b60_fit_report.json`. No constant, switch or row was changed. Wall time 234 s.*

These are R0's numbers only. Which rival (C1 - C3) the data favour is not stated here; that is M1's comparison.

## G1: B27's gate at B60

In each ambient pot (the fit systems B60 was fitted on, `B23.SYSTEMS`), the fraction of the charged oxidant OX (1.0 mM ambient; 62.4 mM cystine in Zhang's GCys arm) that the OX-consuming steps use over the cook, and the relative rise in the pool if the whole mercaptoketone flux delivered one OX each (phi = 1; B27's sum r_nf_mp3p + r_mgo_mp + ch_redox_mp3p, and r_nf_mp3p alone). The rise is an upper bound: it books the whole cook's supply as if present from t = 0.

| pot | OX used, B9 (recorded) | OX used, B60 | rise at phi = 1, B9 (recorded) | rise at phi = 1, B60 (sum / r_nf_mp3p only) | B9 vector, today's numerics: used / rise |
|---|---|---|---|---|---|
| zhou_arp_cys_pH7 | 0.001 % | 0.000 % | 0.35 % | 0.08 % / 0.08 % | 0.008 % / 0.59 % |
| zhang_fig1_cys | 0.487 % | 0.009 % | 0.41 % | 0.07 % / 0.07 % | 1.303 % / 0.49 % |
| zhang_fig1_gcys | 0.100 % | 0.008 % | 0.01 % | 0.00 % / 0.00 % | 0.143 % / 0.01 % |

**G1 HOLDS at B60.** Largest use of the pool 0.009 %; largest rise at phi = 1 0.08 %, against the ~9x (900 %) a decade on the disulfide share would need. G1 holds when the consumers use < 1 % of the charged OX in every ambient pot and the whole mercaptoketone flux at phi = 1 would raise the pool by < 1 % (B27: < 1 % and < 0.5 %; a decade on the share needs ~9x).

OX consumers carrying flux at B60 (mmol/L over the cook): zhou_arp_cys_pH7: ch_dimer_mft 3.95e-06, ch_dimer_fft 1.71e-08; zhang_fig1_cys: ch_dimer_mft 5.11e-05, ch_dimer_fft 4.02e-05; zhang_fig1_gcys: ch_dimer_mft 0.00243, ch_dimer_fft 0.00237. OX budget closure (|change - net flux| / scale): zhou_arp_cys_pH7 2.8e-11; zhang_fig1_cys 1.2e-09; zhang_fig1_gcys 8.4e-10.

**Reading.** As at B9, the ambient pots are not oxidant-limited: the dimer steps use a negligible share of a pool that is a constant multiplier on their rate, and no oxidant SOURCE of the mercaptoketone's size can move that multiplier. The B9 numbers were computed before B51's pH state and B60's closed flux budget; the 'B9 vector, today's numerics' column separates the vector's change from the numerics': under today's numerics the B9 vector uses up to 1.303 % of the pool, so the fall to B60's numbers comes from the B60 vector, not from the numerics (which constant carries it was not resolved here).

## D1: Schieberle 2000 Table IV, the reference pot's own ratios

100 C, 0.5 M phosphate pH 5, ribose 100 + cysteine 33 mM, one pot read at 30 / 60 / 360 / 720 min. Residual = log10(R0 ratio / measured ratio).

| species | ratio | measured | R0 | residual (dex) | read error (dex) | B9 residual (recorded) |
|---|---|---|---|---|---|---|
| MFT | 360 / 30 min | 34.7 (+1.54 dex) | 2.79 (+0.45 dex) | -1.09 | 0.0050 | -1.80 |
| MFT | 720 / 360 min | 1.15 (+0.06 dex) | 0.541 (-0.27 dex) | -0.33 | 0.0018 | -- |
| MFT | 720 / 30 min | 39.8 (+1.60 dex) | 1.51 (+0.18 dex) | -1.42 | 0.0049 | -2.21 |
| FFT | 360 / 30 min | 55 (+1.74 dex) | 127 (+2.10 dex) | +0.36 | 0.0111 | -0.80 |
| FFT | 720 / 360 min | 1.2 (+0.08 dex) | 0.981 (-0.01 dex) | -0.09 | 0.0026 | -- |
| FFT | 720 / 30 min | 66 (+1.82 dex) | 124 (+2.09 dex) | +0.27 | 0.0110 | -1.21 |

- MFT: measured 45, 138, 1560, 1790 ug/L; R0 12.6, 30.8, 35.2, 19.1 ug/L at 30 / 60 / 360 / 720 min. R0 peaks at 158 min (2.92x its 12 h value); still rising 6 -> 12 h: R0 False, measured True.
- FFT: measured 20, 31, 1100, 1320 ug/L; R0 18, 114, 2.28e+03, 2.24e+03 ug/L at 30 / 60 / 360 / 720 min. R0 peaks at 499 min (1.1x its 12 h value); still rising 6 -> 12 h: R0 False, measured True.

**Noise.** Table IV prints no SD and no n; the only stated error is the transcription error (half the last printed digit), 0.001 - 0.011 dex per point (thiol_submodel.yml); a ratio's error is the two in quadrature. Measurement-model series noise (s042): MFT 0.80 dex, FFT 1.16 dex -- the per-series extra noise the joint thiol fit s042 estimated from its own residual scatter (fitting.estimate_series_noise): a property of that fit's misfit, not a measured replicate error; quoted as the case's sec. 5 says, not used to score.

**Reading.** MFT: R0 rises 2.8x from 30 min to 6 h where the pot rises 35x, peaks at 2.6 h and falls from 6 to 12 h (x0.54; measured x1.15, still rising). Residuals -1.09 / -0.33 / -1.42 dex (360/30, 720/360, 720/30), against B9's -1.80 on 360/30; the s042 series noise is 0.80 dex. FFT: R0 rises FASTER than the pot from 30 min to 6 h (+0.36 dex), then is flat (x0.98 6 -> 12 h; measured x1.20), peak at 8.3 h; both FFT residuals lie inside the s042 FFT series noise (1.16 dex). Every residual is far outside the printed read error, which is the only error the source states. Net: R0 at B60 still does not reproduce the reference pot's late MFT rise; its FFT shape is closer than at B9 but plateaus where the pot keeps rising.

## D2: cysteine and pentose over time

Fractions of the charge. 'Free' is the free amino acid or sugar; '+ reservoir' adds TTCA and the cysteine Amadori compound (each holds one cysteine and one pentose), i.e. the material not yet destroyed. R0 forms neither (TTCA formation and the cysteine-Amadori route are off in the shipped configuration).

| pot | pool | measured: time to half | R0: time to half | measured at last time | R0 at last time |
|---|---|---|---|---|---|
| cao2017_90C | cysteine free | < 30 min (already 0.19 at the first measured time) | > 300 min (0.67 left at the end) | 0.141 at 300 min | 6.66e-01 |
| cao2017_90C | pentose free | < 30 min (already 0.20 at the first measured time) | 1.82 s | 0.172 at 300 min | 0.00e+00 |
| cao2017_90C | cysteine free plus reservoir | > 300 min (1.06 left at the last measured time) | > 300 min (0.67 left at the end) | 1.062 at 300 min | 6.66e-01 |
| cao2017_90C | pentose free plus reservoir | > 300 min (1.09 left at the last measured time) | 1.82 s | 1.094 at 300 min | 0.00e+00 |
| zhang2025b_120C | cysteine free | ~5.13 min (log-time interpolation between measured points) | 0.524 s | 0.020 at 50 min | 9.27e-03 |
| zhang2025b_120C | pentose free | < 5 min (already 0.46 at the first measured time) | 0.306 s | 0.030 at 50 min | 0.00e+00 |
| zhang2025b_120C | cysteine free plus reservoir | > 50 min (0.52 left at the last measured time) | 0.524 s | 0.515 at 50 min | 9.27e-03 |
| zhang2025b_120C | pentose free plus reservoir | > 50 min (0.53 left at the last measured time) | 0.306 s | 0.526 at 50 min | 0.00e+00 |
| zhai2019_70C | cysteine free | not measured | 37.3 s | -- | -- |
| zhai2019_70C | pentose free | not measured | 38.1 s | -- | -- |
| zhai2019_70C | cysteine free plus reservoir | not measured | 4.65 min | -- | -- |
| zhai2019_70C | pentose free plus reservoir | not measured | 57 s | -- | -- |
| zhai2019_90C | cysteine free | not measured | 15 s | -- | -- |
| zhai2019_90C | pentose free | not measured | 11.6 s | -- | -- |
| zhai2019_90C | cysteine free plus reservoir | not measured | 1.62 min | -- | -- |
| zhai2019_90C | pentose free plus reservoir | not measured | 16.8 s | -- | -- |
| zhai2019_110C | cysteine free | not measured | 6.49 s | -- | -- |
| zhai2019_110C | pentose free | not measured | 3.94 s | -- | -- |
| zhai2019_110C | cysteine free plus reservoir | not measured | 39.8 s | -- | -- |
| zhai2019_110C | pentose free plus reservoir | not measured | 5.58 s | -- | -- |
| zhai2019_130C | cysteine free | not measured | 2.96 s | -- | -- |
| zhai2019_130C | pentose free | not measured | 1.49 s | -- | -- |
| zhai2019_130C | cysteine free plus reservoir | not measured | 18.6 s | -- | -- |
| zhai2019_130C | pentose free plus reservoir | not measured | 2.07 s | -- | -- |

Zhai 2019 measures TTCA and the cysteine Amadori compound only (no free cysteine or xylose). Their sum is a LOWER bound on the cysteine not destroyed (fraction of the 82.6 mM charge), against R0's reservoir:

| pot | measured reservoir / charge (lower bound on cysteine left; first, last time) | R0 reservoir (mmol/L; first, last time) | R0 cysteine free + reservoir / charge (first, last time) |
|---|---|---|---|
| zhai2019_70C | 0.50 at 20 min, 0.69 at 160 min | 13.2, 6.72 | 0.417, 0.391 |
| zhai2019_90C | 0.60 at 20 min, 0.81 at 160 min | 10.4, 1.03 | 0.414, 0.375 |
| zhai2019_110C | 0.74 at 20 min, 0.52 at 160 min | 5.32, 0.00476 | 0.421, 0.324 |
| zhai2019_130C | 0.61 at 20 min, 0.08 at 160 min | 0.955, 5.15e-09 | 0.411, 0.211 |

Endpoint check, dense grid vs the engine's own one-segment path (max |dex| over Cys, PENT): cao2017_90C 0.0e+00; zhang2025b_120C 0.0e+00; zhai2019_70C 0.0e+00; zhai2019_90C 0.0e+00; zhai2019_110C 0.0e+00; zhai2019_130C 0.0e+00; Schieberle pot (Cys, PENT, MFT, FFT) 0.0e+00.

**Noise.** Cao 2017 and Zhang 2025b: printed SD (n = 3) where an error bar is drawn, plus the digitisation error (route_submodel.yml); Zhang 2025b's cysteine and xylose calibrations print identical (CF-063) and its free series are read as ln(C/C0) with C0 taken as the charge. Zhai 2019: raster read error only (no SD printed).

**Reading.** The pentose: R0 consumes it in seconds in every pot (Cao 90 C 1.82 s, Zhang 120 C 0.306 s), where Cao keeps 1.09 of it (free + reservoir) at 5 h and Zhang 0.53 at 50 min. The cysteine: in Zhang's pot R0 consumes it in 0.524 s (measured free half-life ~5.13 min (log-time interpolation between measured points); 0.52 still accounted at 50 min). In Cao's equimolar 0.4 M pot R0 does NOT consume it in seconds: the pentose runs out first and 0.67 of the cysteine is left free at 5 h, i.e. a third destroyed where the pot destroys none measurably (1.06 accounted); the pot's FREE cysteine is lower than R0's (0.14) because it sits in TTCA, which R0 does not form. Zhai 2019 (70C: free Cys halves in 37.3 s, pentose in 38.1 s; at 160 min R0 holds 0.39 of the cysteine (free + reservoir) where the reservoir alone shows >= 0.69; 90C: free Cys halves in 15 s, pentose in 11.6 s; at 160 min R0 holds 0.38 of the cysteine (free + reservoir) where the reservoir alone shows >= 0.81; 110C: free Cys halves in 6.49 s, pentose in 3.94 s; at 160 min R0 holds 0.32 of the cysteine (free + reservoir) where the reservoir alone shows >= 0.52; 130C: free Cys halves in 2.96 s, pentose in 1.49 s; at 160 min R0 holds 0.21 of the cysteine (free + reservoir) where the reservoir alone shows >= 0.08). R0's reservoir there is only the measured t = 0 TTCA (14.5 mM) dissociating; it never forms more. At 130 C the measured reservoir also empties by 160 min, so the bound says little late in that pot.

## Rows excluded as held out (not read)

- mp_holdout_xylose_cysteine_phosphate_120C_50min_pH65_Zhang2025b (Table 2 blank: MFT, FFT, MFT disulfide, 2-acetylthiazole at 50 min): frozen Maillard-path hold-out bundle (FIT_HOLDOUT_DECLARATION.md Amendment 45); D2 integrates that pot but reads and reports only cysteine, xylose, TTCA and cysteine-Amadori.
- Kumazawa 2003 (F28#3), Esterbauer 1975 (F28#17), Stack 2018 (F28#18): contrast_split_v1 TEST; not used.
- Ehrenberg 1989 (F28#8), Xu 2011 (F28#20), Wang 2022 (F28#25), Liu 2023b (F28#26): split_v2_per_change E12 TEST contrasts; not used.

Every row read is DEVELOPMENT in `contrast_split_v1` and has no TEST contrast in `split_v2_per_change` (laboratories btbu_beijing, dfa_garching, ho_network; checked by this generator).

