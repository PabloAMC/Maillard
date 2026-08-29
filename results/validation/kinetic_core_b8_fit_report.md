# Build Wave B8 -- the fit

Pre-registered in `results/validation/kinetic_core_b8_prereg.md`. The final parameter wave of the kinetic-core rebuild: it implements Amendments 16 and 17 and refits the sulfur lane on the updated FIT set at the single declared weighting B2.4's own exam-blind and panel-blind criterion shipped.

## The objective

- **62 declared FIT rows** (58 of them B2.4's, 4 new: `feng_ARP_conversion_100C`, `feng_ARP_conversion_120C`, `zhai_MFT_fold_120_over_100`, `zhai_FFT_fold_120_over_100`)
- **23 free of 48**
- pH exchange rate **E = 0.7** (sigma_ph = 0.5 pH units)
- final cost **20.0824**
- `sum_r2_level` on the rows B2.4 also scored: **38.718**; on all B8 rows: 38.960

TOTAL COST IS NOT COMPARABLE TO B2.4's: B8 scores four rows B2.4 did not. `sum_r2_level_shared_with_b2_4` is the like-for-like comparator and is summed over exactly the non-pH rows both waves scored.

## The members

| start | arm | start cost | final cost | evals | status | converged | s |
|---|---|---:|---:|---:|---:|---|---:|
| 0 | incumbent | 57.747 | 20.0882 | 425 | 3 | yes | 1526 |
| 1 | local | 512.786 | 20.0824 | 348 | 3 | yes | 1390 |

Best: start **1**.

## The T-structure

| barrier | band | prior centre | fitted | state |
|---|---|---:|---:|---|
| `Ea_decay_thiol_sink` (covalent capture) | [7.0, 102.0] | 60.2 | **102.00** | **ON THE CEILING** |
| `Ea_decay_carbonyl_sink` | [20.0, 250.0] | -- | **166.67** | interior |

Measured barriers the fit cannot move:

| step | Ea (kJ/mol) |
|---|---:|
| `k_dimer_mft` | 122.2 |
| `k_dimer_fft` | 122.2 |
| `k_arp_dpo` | 85.7 |
| `k_arp_tdp` | 85.7 |
| `k_cys_thermal` | 55.1 |

**REFUTED:** Ea_decay_thiol_sink = 248.0 (B2.2), 216.1 (B2.3), 212.9-218.1 (B2.4). Gigl 2021 measures the covalent-capture channel at k(333)/k(279) = 67.2; 248 kJ/mol predicts 3.4e7 -- wrong by 5.1e5x.

## The frozen vector

- lumped formation Ea: **64.08 kJ/mol**
- pH drift: acid yield 0.0000, ARP ammonium pKa 7.092

| log10 k_ref at 145 C | value |
|---|---:|
| `k_arp_dpo` | -2.325 |
| `k_arp_dpo_th` | -1.290 |
| `k_arp_tdp` | -1.790 |
| `k_arp_tdp_th` | -4.354 |
| `k_cys_actz` | -2.974 |
| `k_cys_thermal` | -2.052 |
| `k_ddp_mft` | -6.541 |
| `k_ddp_mft_hs` | -7.693 |
| `k_dimer_decay` | -9.563 |
| `k_dimer_fft` | 0.500 |
| `k_dimer_mft` | 0.500 |
| `k_dpo_c2c3` | -3.016 |
| `k_dpo_ddp` | 0.042 |
| `k_dpo_nf` | 0.434 |
| `k_dpo_ptr` | -4.347 |
| `k_fft_decay` | -0.657 |
| `k_fur_decay` | 0.436 |
| `k_fur_fft` | -0.550 |
| `k_fur_fft_hs` | -0.010 |
| `k_glc_fur` | -4.923 |
| `k_glc_ha` | -0.757 |
| `k_h2s_loss` | -1.351 |
| `k_ha_mp_mft` | -3.870 |
| `k_hmp_mft` | -2.609 |
| `k_hmp_mp2p` | -3.437 |
| `k_mft_decay` | 0.134 |
| `k_mgo_mp` | -3.891 |
| `k_mmft` | -1.984 |
| `k_nf_decay` | -2.182 |
| `k_nf_mft` | -2.042 |
| `k_nf_mp3p` | -2.859 |
| `k_osone_decay` | -1.213 |
| `k_pent_caramel` | -2.597 |
| `k_pent_dpo` | -0.545 |
| `k_pent_tdp` | -0.502 |
| `k_pent_thermal` | -3.437 |
| `k_tdp_fur` | -3.005 |
| `k_thi_hmp` | -2.588 |
| `k_thi_mesh` | -3.273 |
| `k_thiol_decay` | -9.483 |
| `k_thiolate_loss` | -1.834 |
| `k_ttca_cys` | -0.579 |
| `k_ttca_deg` | -1.694 |

## Bound adjacency -- what is estimated and what is only bounded

A parameter pressed against a bound is not an estimate, and a report that prints it like one is misleading. Checked on every coordinate the fit was free to move:

| coordinate | bounds | fitted | state |
|---|---|---:|---|
| `Ea_decay_thiol_sink` | [7.0, 102.0] | 102 | **ON THE CEILING -- bound-limited, not an estimate** |
| `Ea_decay_carbonyl_sink` | [20.0, 250.0] | 166.7 | interior |
| `ph_acid_yield_per_sink_event` | [0.0, 1.0] | 3.764e-06 | **ON THE FLOOR -- bound-limited, not an estimate** |
| `ph_arp_secondary_ammonium_pKa` | [5.0, 11.0] | 7.092 | interior |

- `Ea_decay_thiol_sink` is on its ceiling for the FOURTH consecutive wave (248 against 250, 216, 213-218, now 102 against 102). The band changed; the objective's appetite did not. It still wants more temperature dependence on the thiol sink than the measurement allows, and B8 does not widen the band to give it: Gigl's range is a measurement, the conflict is real, and it is reported rather than dissolved.
- `ph_acid_yield_per_sink_event` sits at ~4e-6 against a floor of 0. This is NOT new in B8: B2.4-half was already at ~1e-4, having come down from B2.3's 0.968. The pH trajectory is now driven almost entirely by the ARP ammonium pKa rather than by acid production, and that constant is a fitted nuisance parameter living in a pKa slot -- `kinetic_core_b2_2_diagnosis.md` sec. 5 item 2 says so and it is still true. **Neither may be quoted as an acid-base constant.**
- No parameter interval is computed by this wave, so `identified` is unknown for all 48. B2.3's report is the last artifact that evaluated them and found 3 of 48. Nothing in B8 adds identifiability: it removes freedom (four barriers become immovable) and adds four rows.

*The blind re-sits are in `results/validation/kinetic_core_b8_holdout_report.md`. No score in that file entered any choice made in this one.*
