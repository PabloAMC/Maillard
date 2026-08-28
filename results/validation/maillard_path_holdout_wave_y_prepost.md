# Wave Y vs. the frozen out-of-sample baselines (2026-08-28)

Baseline: `results/validation/maillard_path_holdout_frozen_predictions.json` at commit `12f43dd`. READ, never regenerated. Comparing a wave against a baseline the same wave produced proves nothing, which is why this file is committed and append-only.

## maillard_path hold-out — **8 of 22 targets moved**

The 22 targets the frozen pre-registration itself carries. Bundles added after 12f43dd (Waves W and X) are listed under targets_added_since_freeze and are NOT counted as movement -- they have no frozen counterpart to move from.

Targets added since the freeze (10, not scored here): `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 / 2-Furfurylthiol (FFT)`, `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3 / 2-Methyl-3-furanthiol (MFT)`, `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 / 2-Furfurylthiol (FFT)`, `mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7 / 2-Methyl-3-furanthiol (MFT)`, `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 / 2-Furfurylthiol (FFT)`, `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3 / 2-Methyl-3-furanthiol (MFT)`, `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 / 2-Furfurylthiol (FFT)`, `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7 / 2-Methyl-3-furanthiol (MFT)`, `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 / 2-Furfurylthiol (FFT)`, `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5 / 2-Methyl-3-furanthiol (MFT)`

The summary block below is the WHOLE current hold-out against the whole frozen one, so the count rows differ by the post-freeze additions; the movement result above is the like-for-like comparison.

| metric | frozen (12f43dd) | Wave Y (as shipped) |
|---|---:|---:|
| `bundle_count` | 12 | 17 |
| `target_count` | 22 | 32 |
| `median_fold_error` | 6.038752949192657 | 10.86381837964712 |
| `median_abs_log10_error` | 0.7809472625943955 | 1.0359824965004563 |
| `worst_fold_error` | 52.58641429978577 | 506.424391151322 |
| `best_fold_error` | 1.5167561508937306 | 1.2435497830531987 |
| `within_10x` | 12 | 15 |
| `ordinal_pairs_correct` | 8 | 10 |
| `series_directions_correct` | 3 | 3 |

### Targets that moved

| benchmark | compound | frozen | Wave Y |
|---|---|---:|---:|
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026` | 2-Furfurylthiol (FFT) | 15.559603988705325 | 15.43934132907918 |
| `mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026` | 2-Methyl-3-furanthiol (MFT) | 4.535996109820317 | 8.555622507406007 |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026` | 2-Furfurylthiol (FFT) | 33.517542910781415 | 33.174756204391144 |
| `mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026` | 2-Methyl-3-furanthiol (MFT) | 8.490127310987353 | 15.70902971000651 |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026` | 2-Furfurylthiol (FFT) | 62.075063043156405 | 61.28098212094799 |
| `mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026` | 2-Methyl-3-furanthiol (MFT) | 14.37462009088175 | 26.073164111153087 |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026` | 2-Furfurylthiol (FFT) | 85.18999116565296 | 83.9926871854059 |
| `mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026` | 2-Methyl-3-furanthiol (MFT) | 20.645018442692226 | 36.734806615858105 |

## Attribution — who moved them

| ablation | frozen targets still differing |
|---|---:|
| none (as shipped) | **8** |
| Wave X norfuraneol channel replaced with a no-op | **0** |
| Wave Y's five constants reverted in memory | **0** (vs. the shipped run) |

EVERY target that has moved since the freeze is attributable to Wave X's norfuraneol parallel channel and to nothing else: with that one template ablated, all 22 frozen targets are bit-identical to the pre-registration. Wave Y moves none of them, which is what its pre-registration Y-P1 claimed mechanistically -- but Y-P1's literal wording ('0 of 22 targets move') is FALSIFIED against the frozen file, because it assumed the frozen file was still current. It was not, and nobody had scored it.

> Wave X priced its norfuraneol channel against the Wave W panel (0.9241 -> 0.9518 dex) and did not re-score the frozen maillard_path hold-out. The unmeasured price is 8 of 22 frozen targets and a median fold error of 6.0388x -> 10.8638x on the repository's only free-precursor out-of-sample surface. All eight are PENTOSE (ribose/cysteine) rows, which is exactly what a norfuraneol channel should touch: MFT rises 1.78-1.89x on all four bundles and FFT falls 0.992x. That is the channel's fingerprint, and it is why the ablation closes the gap to zero.

## Wave X step-level rows — mean \|log10\| over the 5 scored rows: **0.5177 dex** (Wave X reported 0.518)

| benchmark | analyte | measured | predicted | fold | \|dex\| |
|---|---|---:|---:|---:|---:|
| `hofmann1998_furan2aldehyde_h2s_145C_20min_pH5` | 2-Furfurylthiol (FFT) | 11016 | 2413.65 | 0.2191x | 0.6593 |
| `hofmann1998_norfuraneol_h2s_145C_20min_pH5` *(FIT TARGET — excluded)* | 2-Methyl-3-furanthiol (MFT) | 4224 | 1533.06 | 0.3629x | 0.4402 |
| `hofmann1998_norfuraneol_cysteine_145C_20min_pH5` | 2-Methyl-3-furanthiol (MFT) | 1016 | 1745.31 | 1.7178x | 0.2350 |
| `hofmann1998_c2c3_recombination_145C_20min_pH3` | 2-Methyl-3-furanthiol (MFT) | 310 | 2180.34 | 7.0334x | 0.8472 |
| `hofmann1998_c2c3_recombination_145C_20min_pH5` | 2-Methyl-3-furanthiol (MFT) | 5362 | 2180.34 | 0.4066x | 0.3908 |
| `hofmann1998_c2c3_recombination_145C_20min_pH7` | 2-Methyl-3-furanthiol (MFT) | 6230 | 2180.34 | 0.3500x | 0.4560 |
