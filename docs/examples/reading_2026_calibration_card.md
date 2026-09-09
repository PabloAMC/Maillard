================================================================================================
  CALIBRATION CARD   laboratory: Reading 2026   base: kinetic_core_b9   matrix: water
================================================================================================

  records: fit yiltirak_100C, yiltirak_120C
           validate yiltirak_110C, yiltirak_130C

  RESPONSE FACTORS  (from the fit records' levels; a property of the measurement, not the chemistry)
    2-methyl-3-furanthiol            x0.0402   (log10 -1.40 +/- 0.26, 2 rows)
    2-furfurylthiol                  x0.00209   (log10 -2.68 +/- 0.05, 2 rows)

  KINETIC OVERRIDES  (from contrasts only; pulled toward the shipped value by its shipped sigma)
    b8.k_fur_decay.log10_k_ref_145C          +0.470 -> -0.304   (shift -0.77; prior sigma 0.44, posterior 0.24)
    b8.k_fur_fft.log10_k_ref_145C            -0.631 -> -1.678   (shift -1.05; prior sigma 1.03, posterior 1.02)
    contrasts 2, candidates 16, identified 2
    not identified by the contrasts (left at the shipped value): b8.k_arp_dpo.log10_k_ref_145C, b8.k_arp_tdp.log10_k_ref_145C, b8.k_tdp_fur.log10_k_ref_145C, b8.k_nf_mft.log10_k_ref_145C, b8.k_nf_mp3p.log10_k_ref_145C, b8.k_mgo_mp.log10_k_ref_145C, b8.k_dimer_mft.log10_k_ref_145C, b8.k_dimer_fft.log10_k_ref_145C, b8.k_fft_decay.log10_k_ref_145C, b8.k_osone_decay.log10_k_ref_145C, b8.k_pent_caramel.log10_k_ref_145C, b8.k_ttca_deg.log10_k_ref_145C ...

  HOLD-OUT  (validate records, never fitted)
    median fold error   shipped 115   calibrated 1.16   (within 3x: 0/4 -> 3/4)
      yiltirak_110C                2-methyl-3-furanthiol      measured 3.29   shipped     29.5x   calibrated     1.11x
      yiltirak_110C                2-furfurylthiol            measured 1.46   shipped      367x   calibrated     1.11x
      yiltirak_130C                2-methyl-3-furanthiol      measured 1.71   shipped     99.7x   calibrated     3.79x
      yiltirak_130C                2-furfurylthiol            measured 1.62   shipped      131x   calibrated     1.22x

  note: roles as tagged in the document

  The shipped parameters are untouched. Apply this file with --calibration on compare, predict or score.