# Wave B11 ship rule -- **DO NOT SHIP: T1 and T3 both fail -- the consumers ship as declared-inert**

Prereg: `results/validation/kinetic_core_b11_prereg.md` sec. 5. Rule: SHIP if T1 (Bolton < 6x), T3 (Yiltirak 130 C both >= 3x better) and T5 (no row +0.3 dex) hold (prereg sec. 5).

* fitted consumers: k_cys_ox **0.00037**, k_red_ox **0.000126** per unit per min (bands 1e-5..1e-1); supply 10 /min; saturation 0.3 mmol/L
* active bounds at the optimum: ['k_dimer_mft', 'k_dimer_fft', 'Ea_decay_thiol_sink']

## T1 -- Bolton 1994 (O2 : thiol 2.07)

* MFT fold error 20.2x (B9) -> **20.2x** (B11); reservoir 81 units (92 mL headspace over 33.3 mL: 24.2 mmol O2 per litre) -> below 6x: **FAIL** (FALSIFIED: above 12x)

## T2 -- Hofmann 1998 Table-1 pH-5 rows (O2 : thiol 0.27; the eight B9 validation rows)

| row | B9 fold | B11 fold |
|---|---:|---:|
| hofmann_ribose_FFT | 1.6 | 1.5 |
| hofmann_ribose_MFT | 2.8 | 2.8 |
| hofmann_xylose_FFT | 2.0 | 1.8 |
| hofmann_xylose_MFT | 3.8 | 3.9 |
| hofmann_glucose_FFT | 9.8 | 10.3 |
| hofmann_glucose_MFT | predicts ZERO | predicts ZERO |
| hofmann_fructose_FFT | 6.4 | 6.6 |
| hofmann_fructose_MFT | predicts ZERO | predicts ZERO |
* worst answered row 9.8x -> **10.3x** (2 hexose rows predict zero under both, the B9 finding); none beyond 8x: **FAIL**

## T3 -- Yiltirak 2026, 130 C / 0.5 h (oxygen in excess)

* 2-Methyl-3-furanthiol (MFT): 99.7x -> **102.5x** (improvement x0.97)
* 2-Furfurylthiol (FFT): 130.7x -> **129.3x** (improvement x1.01)
* both at least 3x better: **FAIL** (FALSIFIED: neither 2x)

## T4 -- Yiltirak 100 C / 4 h (recorded, not a target)

* 2-Methyl-3-furanthiol (MFT): 9.3x -> 9.5x
* 2-Furfurylthiol (FFT): 476.8x -> 408.2x

## T5 -- in-sample discipline

* 51 shared rows; worst growth +0.20 dex on `zhang_fig1_gcys_dimer_over_MFT`; 0 rows over 0.3 dex, 0 over 0.5 dex -> **PASS**

## T6 -- fed-intermediate rows

* 12 rows; worst shift 0.05 dex; all within 2x of their B9 residual: **PASS**

## Laplace

* 15 of 25 coordinates identified (at least 20: False); k_cys_ox identified = False, k_red_ox identified = False (both EXPECTED unidentified: one laboratory's vessel, no oxygen contrast)
