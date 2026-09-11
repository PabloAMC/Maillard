# Wave B25 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (reference pot ratios within 0.5 dex and still rising 6 -> 12 h), T2 (no B9 row +0.3 dex) and T6 (the addition constant and its barrier identified, off their bounds, slices not flat) hold; T3-T5 reported. Pre-registration `results/validation/kinetic_core_b25_prereg.md`.*

log10 k_add at 145 C: -5.999 (band [-6.0, 0.0]); Ea_add 120.0 kJ/mol (band [10.0, 120.0]); cost 931.27; active bounds ['k_dimer_mft', 'k_dimer_fft', 'k_pent_caramel', 'Ea_decay_thiol_sink', 'log10_k_add', 'Ea_add']

| test | result | pass |
|---|---|---|
| T1 reference pot 100 C | MFT [22.7, 55.3, 109.6, 96.0] ug/L at 30/60/360/720 min; ratios (dex) {'MFT_360_over_30': -0.856791776384977, 'MFT_720_over_30': -0.9741039461725843, 'FFT_360_over_30': -0.27034601682556386, 'FFT_720_over_30': -0.3255665433906534}; rising 6->12 h {'MFT': False, 'FFT': True} | False |
| T2 B9 rows | worst fed_ribose_h2s_MFT +1.90 dex; 16 rows over 0.3 | False |
| T3 dimer shares | Zhou 2023: { pH 6.0: model 0.03 % vs 8.6 %, pH 7.0: model 0.31 % vs 6.5 %, pH 8.0: model 0.72 % vs 9.6 % }; Zhang 2024 Cys arm (figure-derived, reported only): model 0.38 % vs 8.7 %; worst (Zhou) 2.46 dex | False |
| T4 Yiltirak | median fold 13.7 (B9 115; B9 under this rule: 115.2) | True |
| T5 Wang 2022 140 C shape | MFT peak 30 min, decline 1.32 dex; FFT peak 60 min, decline 0.11 dex | False |
| T6 identification | sigma 24389.45452250456 dex / 2732809.3525724374 kJ/mol, identified {'log10_k_add': False, 'Ea_add': False}, on bound {'log10_k_add': True, 'Ea_add': True}, slices {"log10_k_add": "bound_limited", "Ea_add": "bound_limited"} | False |

## The reference pot under B9, for comparison

- MFT [65.9, 130.3, 36.4, 16.2] ug/L; ratios {'MFT_360_over_30': -1.7979941506348263, 'MFT_720_over_30': -2.2078951029456637, 'FFT_360_over_30': -0.8021645937164703, 'FFT_720_over_30': -1.2074159427622306}; rising {'MFT': False, 'FFT': False}
- Wang 140 C under B9: MFT decline 2.96 dex, FFT 2.45 dex
- dimer shares under B9: Zhou { pH 6.0: 0.13 %, pH 7.0: 0.59 %, pH 8.0: 1.97 % }, Zhang 0.43 %
