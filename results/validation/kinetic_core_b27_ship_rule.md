# Wave B27 ship rule: NOT FITTED -- DO NOT SHIP

*Rule: SHIP if T1, T2, T3 (Zhou 2023's three dimer shares within 0.3 dex) and T6 hold; evaluated here as a GATE before the fit, at the shipped B9 vector with phi swept over its band. Pre-registration `results/validation/kinetic_core_b27_prereg.md`.*

**The fit was not run.** The decisive test cannot be reached by this structure or by any oxidant source, and the one pot the structure does fix needs its single coordinate at the physical ceiling the pre-registration declared disqualifying.

| gate | result | pass |
|---|---|---|
| G1 T3 reachability | zhou_arp_cys_pH7: consumers use 0.001 % of the pool, whole mercaptoketone flux at phi = 1 adds 0.35 %; zhang_fig1_cys: consumers use 0.487 % of the pool, whole mercaptoketone flux at phi = 1 adds 0.41 %; zhang_fig1_gcys: consumers use 0.100 % of the pool, whole mercaptoketone flux at phi = 1 adds 0.01 %; a decade needs about 9x the pool | False |
| G2 Whitfield share | phi 0.0001: 0.0 %, phi 0.01: 0.6 %, phi 0.1: 5.5 %, phi 0.2: 10.3 %, phi 0.316: 15.4 %, phi 0.5: 22.3 %, phi 0.708: 28.9 %, phi 1: 36.4 %; the 35 % floor is first reached at phi = 1.0 | False |
| G3 fit-free T2 | cost 1539.0 (phi -> 0) -> 1459.3 (phi = 1); worst growth whitfield_nf_h2s_MFT +0.00 dex; 0 rows over 0.3; Kumazawa max |growth| 0.0000 dex | True |
| T1 reference pot at phi = 1 | ratios (dex) {'MFT_360_over_30': -1.7979941793558196, 'MFT_720_over_30': -2.207895111496129, 'FFT_360_over_30': -0.80216811470117, 'FFT_720_over_30': -1.2074193247020377}; rising 6->12 h {'MFT': False, 'FFT': False} | False |
| T5 Wang shape at phi = 1 | MFT decline 2.96 dex, FFT 2.44 dex | False |
| the same two under B9 itself (phi absent) | T1 ratios {'MFT_360_over_30': -1.7979941506348263, 'MFT_720_over_30': -2.2078951029456637, 'FFT_360_over_30': -0.8021645937164703, 'FFT_720_over_30': -1.2074159427622306}, rising {'MFT': False, 'FFT': False}, pass False; T5 MFT decline 2.96, FFT 2.45, pass False | reference |

## G1, read plainly

In the ambient pots the oxidant is NOT a budget: the consumers use well under one per cent of it, so it acts as a constant multiplier on the dimer rate. The whole mercaptoketone flux at phi = 1 would raise that multiplier by under one per cent. What those rows need is the dimer RATE CONSTANT, which sits at its band ceiling and is opposed by Kumazawa's retention rows. No oxidant SOURCE reaches T3.

## G2, read plainly

The pot whose oxidant is genuinely zero. The share climbs with phi and reaches the 35 % floor only at phi = 1.0 exactly -- every mercaptoketone-forming event oxidising a thiol. The pre-registration (sec. 10) declared that a phi pinned at its ceiling means the objective is asking the mercaptoketone flux for more than it can supply, and is evidence against the structure rather than a fitted value.

## What phi = 1 does to the rest of the objective

Rows that improve by more than 0.1 dex: whitfield_mft_disulfide_share_floor (-3.78), whitfield_nf_cys_MFT (-0.19).
The four Kumazawa rows move by at most 0.0000 dex, as section 9 predicted (they carry no norfuraneol).

## What is kept

The step `ch_redox_mp3p`, the parameter `k_redox_mp3p` and the engine hook stay in the code, INERT at phi = 0, the way B17's and B25's refused structures do. The three Whitfield charge corrections are installed only by this wave's generator and are restored on exit, so the shipped objective is unchanged -- and that is a named debt: the next sulfur refit must carry the printed charges (norfuraneol 50, cysteine 50, H2S 97 mmol/L; 0.5 M phosphate pH 4.5; total MFT 0.230 mol %).
