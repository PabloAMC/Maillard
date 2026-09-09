# Wave B17 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (reference pot ratios within 0.5 dex and still rising 6 -> 12 h), T2 (no B9 row +0.3 dex) and T6 (release identified, off its bound, slice not flat) hold; T3-T5 reported. Pre-registration `results/validation/kinetic_core_b17_prereg.md`.*

log10 k_dimer_release at 145 C: -8.262 (band [-10.0, 0.5]); cost 931.12; active bounds ['k_dimer_mft', 'k_dimer_fft', 'k_pent_caramel', 'Ea_decay_thiol_sink']

| test | result | pass |
|---|---|---|
| T1 reference pot 100 C | MFT [22.8, 55.6, 110.6, 97.4] ug/L at 30/60/360/720 min; ratios (dex) {'MFT_360_over_30': -0.8547879185589875, 'MFT_720_over_30': -0.9696544513148303, 'FFT_360_over_30': -0.2951333159094513, 'FFT_720_over_30': -0.3467642684335807}; rising 6->12 h {'MFT': False, 'FFT': True} | False |
| T2 B9 rows | worst fed_ribose_h2s_MFT +1.90 dex; 16 rows over 0.3 | False |
| T3 dimer shares | Zhou 2023: { pH 6.0: model 0.04 % vs 8.6 %, pH 7.0: model 0.35 % vs 6.5 %, pH 8.0: model 0.90 % vs 9.6 % }; Zhang 2024 Cys arm: model 0.39 % vs 8.7 %; worst 2.38 dex | False |
| T4 Yiltirak | median fold 13.6 (B9 115; B9 under this rule: 115.2) | True |
| T5 Wang 2022 140 C shape | MFT peak 30 min, decline 1.32 dex; FFT peak 60 min, decline 0.11 dex | False |
| T6 identification | sigma 47665.31623985493 dex, identified False, on bound False, slice flat {'-1.0': 931.1241397488964, '-0.5': 931.124129759625, '+0.0': 931.1242747147551, '+0.5': 931.124253177657, '+1.0': 931.1241825817193} | False |

## The reference pot under B9, for comparison

- MFT [65.9, 130.3, 36.4, 16.2] ug/L; ratios {'MFT_360_over_30': -1.7979941506348263, 'MFT_720_over_30': -2.2078951029456637, 'FFT_360_over_30': -0.8021645937164703, 'FFT_720_over_30': -1.2074159427622306}; rising {'MFT': False, 'FFT': False}
- Wang 140 C under B9: MFT decline 2.96 dex, FFT 2.45 dex
- dimer shares under B9: Zhou { pH 6.0: 0.13 %, pH 7.0: 0.59 %, pH 8.0: 1.97 % }, Zhang 0.43 %
