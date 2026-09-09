# Wave B10 ship rule -- **RE-MERGE: barriers not identified**

Prereg: `results/validation/kinetic_core_b10_prereg.md` sec. 5. Rule: SHIP if T1 and T2 and T4 (Kang direction) hold and T3 improved at all (prereg sec. 5).

* route barriers: sugar trunk **49.3**, thiol assembly **55.0** kJ/mol (B9's single lumped barrier: 64.1)
* active bounds at the optimum: ['k_dimer_mft', 'k_dimer_fft', 'Ea_decay_thiol_sink', 'Ea_thiol_assembly']

## T1 -- Laplace identification of the two barriers

* sugar_trunk: identified = False, sigma = 74.25325168326071
* thiol_assembly: identified = False, sigma = 52.37755702672499
* **FAIL**

## T2 -- no B9 row moved more than 0.3 dex

* 51 shared rows; worst growth +0.40 dex on `hofmann_ribose_FUR_insitu`; 1 rows over 0.3 dex, 0 over 0.5 dex -> **FAIL**

## T3 -- the eight Yiltirak levels (FIT-ADJACENT, not out of sample)

* median fold error 115.2x (B9) -> **89.9x** (B10); within 3x 0 -> 0 of 8 -> below 10x: **FAIL**; improved: True

| bundle | compound | measured | B10 predicted | fold |
|---|---|---:|---:|---:|
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 6.88 | 69.5 | 10.1 |
| mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.28 | 132 | 102.8 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 3.29 | 93.4 | 28.4 |
| mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.46 | 146 | 100.0 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 2.4 | 120 | 50.0 |
| mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.68 | 157 | 93.4 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Methyl-3-furanthiol (MFT) | 1.71 | 148 | 86.3 |
| mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026 | 2-Furfurylthiol (FFT) | 1.62 | 164 | 101.2 |

## T4 -- hold-outs

* Kang 140 C over 120 C: MFT predicted x0.77 (B9 x0.78; observed x4.26), FFT x0.09 (B9 x0.01; observed x2.79) -> direction **FAIL**
* Hofmann 2002 brew, FFT loss at 80 C: observed 0.023 /min; B9 0.0094 (2.4x), B10 0.0094 (2.4x) -> not worse: True

## T5 -- leave-Yiltirak-out

* with the folds: {'sugar_trunk': 49.31931339833706, 'thiol_assembly': 55.000278224682674}; without: {'sugar_trunk': 50.375002770083256, 'thiol_assembly': 55.000230204042666}; difference: {'sugar_trunk': -1.1, 'thiol_assembly': 0.0} kJ/mol; in Laplace sigmas: {'sugar_trunk': 0.014217416043265171, 'thiol_assembly': 9.168171013381776e-07}
* shared-row cost: with folds 35.65, without 32.31
