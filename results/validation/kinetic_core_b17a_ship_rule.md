# Wave B17 variant (a) ship rule: DO NOT SHIP

*Rule: SHIP if T1 (reference pot ratios within 0.5 dex and still rising 6 -> 12 h), T2 (no B9 row +0.3 dex) and T6 (site yield identified, off its bound, slice not flat) hold; T3-T5 reported. Pre-registration `results/validation/kinetic_core_b17_prereg.md` sec. 2 (a).*

log10 electrophile-site yield per osone decayed: 0.048 (band [-4.0, 0.2]; k_mele_site 2.85e-09 /min at 145 C); cost 930.98; active bounds ['k_dimer_mft', 'k_dimer_fft', 'k_pent_caramel', 'Ea_decay_thiol_sink']

| test | result | pass |
|---|---|---|
| T1 reference pot 100 C | MFT [22.0, 53.6, 107.1, 94.4] ug/L at 30/60/360/720 min; ratios (dex) {'MFT_360_over_30': -0.8525210141402069, 'MFT_720_over_30': -0.9672730466160062, 'FFT_360_over_30': -0.30481419394856163, 'FFT_720_over_30': -0.35406301459513984}; rising 6->12 h {'MFT': False, 'FFT': True} | False |
| T2 B9 rows | worst fed_ribose_h2s_MFT +1.91 dex; 17 rows over 0.3 | False |
| T3 dimer shares | Zhou 2023: { pH 6.0: model 0.03 % vs 8.6 %, pH 7.0: model 0.32 % vs 6.5 %, pH 8.0: model 0.93 % vs 9.6 % }; Zhang 2024 Cys arm (figure-derived, reported only): model 0.38 % vs 8.7 %; worst (Zhou) 2.43 dex | False |
| T4 Yiltirak | median fold 13.4 (B9 115; B9 under this rule: 115.2) | True |
| T5 Wang 2022 140 C shape | MFT peak 30 min, decline 1.29 dex; FFT peak 60 min, decline 0.11 dex | False |
| T6 identification | sigma 85765.98474895378 dex, identified False, on bound False, slice flat {'-1.0': 930.9782792472745, '-0.5': 930.9785425677112, '+0.0': 930.9777880612544, '+0.5': 930.9785128014848, '+1.0': 930.9785128014848} | False |

## The reference pot under B9, for comparison

- MFT [65.9, 130.3, 36.4, 16.2] ug/L; ratios {'MFT_360_over_30': -1.7979941506348263, 'MFT_720_over_30': -2.2078951029456637, 'FFT_360_over_30': -0.8021645937164703, 'FFT_720_over_30': -1.2074159427622306}; rising {'MFT': False, 'FFT': False}
- Wang 140 C under B9: MFT decline 2.96 dex, FFT 2.45 dex
- dimer shares under B9: Zhou { pH 6.0: 0.13 %, pH 7.0: 0.59 %, pH 8.0: 1.97 % }, Zhang 0.43 %
