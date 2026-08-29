# Build Wave B8 -- blind re-sits

Pre-registered in `results/validation/kinetic_core_b8_prereg.md`, which was written and saved to disk BEFORE any B8 fit member ran and BEFORE any scorer was re-run. Nothing in this file fits anything.

## Three columns, not two

d1_exam_panel_reconciliation.md found 96% of the B2.2 -> B2.3 exam regression was an undeclared objective change rather than the physics that wave shipped. B8 changes BOTH, so reporting only the ends would let the refit take credit for the measurements or the measurements take the blame for the refit.

- **pre_b8** -- The artifacts already on disk at HEAD. PANEL: B2.4-half's own vector under the pre-B8 code, so `pre_b8 -> physics_only` is a clean A/B on ONE vector in which only the code changed. EXAM: B2.3's vector under the pre-B8 code, which is what shipped before this wave -- B2.4's exam artifacts predate B7's trunk change and are not like-for-like with anything B8 produces, so the exam's pre-B8 column carries a vector change as well as a code change and is read as a SHIPPING comparison rather than as a controlled one.
- **physics_only** -- B2.4-half's OWN vector, UNREFITTED, under B8's T-structure. The FOUR MEASURED BARRIERS apply (k_dimer_mft and k_dimer_fft at 122.2 where they had none, k_arp_dpo and k_arp_tdp at 85.7 instead of the lumped 64.08); not one rate constant moves. NOTE, PRECISELY: this column keeps B2.4-half's OWN `Ea_decay_thiol_sink` of 215.7 -- Gigl's (7, 102) band is a constraint on the SEARCH and B2.4 was not searched under it, so clipping it here would mix the measurements in with the band and blunt the decomposition. pre-B8 -> physics is therefore the four measured barriers ALONE; physics -> B8 is the refit, which is also where the thiol sink enters its band.
- **b8** -- B8's refitted vector under B8's T-structure -- the ship

## Pre-registration scorecard -- every claim, held or not

Scored MECHANICALLY, in code, from the artifacts. A pre-registration graded by hand after the numbers are seen is graded by someone who has seen the numbers.

| # | claim | outcome | verdict |
|---|---|---|---|
| **F-1** | Ea_decay_thiol_sink lands ON its new 102 ceiling (within 2 kJ/mol) | 102.00 kJ/mol -- ON THE CEILING | **HELD** |
| **F-2** | sum_r2_level on the shared rows WORSENS against B2.4-half's 15.123 AND lands in 15.5-25.0 | 38.718 (B2.4-half: 15.123); worsened: yes; in the 15.5-25.0 band: NO | **HALF-FALSIFIED** |
| **F-3** | the lumped formation Ea stays inside 40-100 kJ/mol | 64.08 kJ/mol | **HELD** |
| **F-6** | still <= 5 of 48 parameters carry a finite Gauss-Newton interval | NOT COMPUTED -- B8's fit report does not evaluate intervals | **UNSCORED** |
| **H-1** | on the OLD basis B8's gating total is within +/-2 of B2.4-half's | old basis 8 / 32 against B2.4-half 12 / 32 | **FALSIFIED** |
| **H-2** | had they been scored, BOTH retired rows would still have FAILED | kang_switch_on_MFT: FAIL; kang_switch_on_FFT: FAIL | **HELD** |
| **H-3** | the two LOW-TEMPERATURE rows are at risk from the barrier drop; at least one regresses | regressed: neither | **FALSIFIED** |
| **H-4** | kang_140C_FFT stays failing | FAIL at 0.00163x | **HELD** |
| **H-5** | the three zhou_120C_dimer_share_pH* rows move by >= 1.5x in either direction (they are the only rows that directly score the dimer branch the T-structure re-barriers) | zhou_120C_dimer_share_pH6 x0.29; zhou_120C_dimer_share_pH7 x0.36; zhou_120C_dimer_share_pH8 x0.23 | **HELD** |
| **X-1** | THE HEADLINE -- the Yiltirak family median improves by >= 3x, to <= 65x | 193x -> 122x | **FALSIFIED** |
| **X-2** | Yiltirak enters the band on at least 1 of 8 | 0 / 8 | **FALSIFIED** |
| **X-3** | the exam's core paired median improves against 42.23x | 42.2x -> 24.8x | **HELD** |
| **X-4** | the sulfur_hofmann1998_145C family median moves by < 1.5x either way (its rows sit AT the k_ref reference temperature, where a barrier change is nearly inert by construction) | 11.5x -> 12.2x (x1.07) | **HELD** |
| **X-5** | the acrylamide, furanic and lipid families are BIT-IDENTICAL -- no B8 parameter touches them and the covalent retirement is documentation over a term that already contributed exactly 0.0 | 16 points checked, 0 moved | **HELD** |
| **X-6** | the exam's declined points stay declined and the answered count stays 34 | declined 6 -> 6; answered 34 | **HELD** |
| **N-1** | `wang2026_MFT_peak_and_fall_125_over_115` FAILS, and in a named direction (the model predicts a RISE where the measurement falls) | PASS at 1.21x (predicted 0.362 against observed 0.3) | **FALSIFIED** |
| **N-2** | `wang2026_FFT_peak_and_fall_115_over_105` FAILS, and in a named direction (the model predicts a RISE where the measurement falls) | FAIL at 5.43x (predicted 0.815 against observed 0.15) | **HELD** |
| **N-3** | the 13C threshold is reported as a SCOPE GAP rather than as a pass or a fail | scored FAIL on the in-scope proxy; comment states the scope gap: yes | **HELD** |
| **N-4** | the excess-Ea SIGN SPLIT reproduces (MFT and FFT positive, 2-acetylthiazole negative) in the module's aqueous pot | NOT reproduced: {'MFT': -130.80381387888602, 'FFT': -143.80284264851113, 'ACTZ': 59.084879622697095} | **FALSIFIED** |

- **F-1** (HELD): Pre-registered at 60% BECAUSE it predicts the wave's own headline parameter coming back bound-limited again. B2.2/B2.3/B2.4 all wanted 213-248.
- **F-2** (HALF-FALSIFIED): The direction was right and the SIZE was badly underestimated: the shared-row residual is 2.6x B2.4-half's, not the ~1.0-1.7x the band allowed. That gap IS the disagreement between the corpus's measured barriers and the fitted ones they replaced, and 23 free parameters could not absorb it.
- **F-3** (HELD): it was 64.08 in B2.3
- **F-6** (UNSCORED): B2.3's report is the last artifact that computed them (3 of 48). Reported as unscored rather than assumed.
- **H-2** (HELD): reported scored-but-not-counted so the retirement is visibly not a rescue
- **H-3** (FALSIFIED): The reasoning was that a thiol sink barrier falling from 216 to <= 102 with k_ref anchored at 145 C makes the sink relatively FASTER when cold -- the mechanism that destroyed both rows in B2. It did not happen, and the reason is visible in the fitted vector: the refit paid for the lower barrier by moving k_ref, so the sink's ABSOLUTE rate at 50-80 C did not rise. The price the wave expected to pay at low temperature was paid on the pH-8 column instead.
- **X-1** (FALSIFIED): The decay-turnover mechanism is the named candidate for its 0/8: a core whose sink barely runs cannot produce a peak, so it must over-predict every long hold.
- **X-5** (HELD): Named in the prereg sec. 4 as the first thing that would falsify the whole wave.
- **N-1** (FALSIFIED): ★ A PRE-REGISTERED FAIL THAT PASSED, and it is the wave's best single result. The claim was that the model predicts a RISE here; it predicts 0.362 against a measured 0.3 -- A FALL, inside the 2x band. THE TURNOVER EMERGED. This row is the reason the T-structure was built and it is reported as a falsified PREDICTION, not as a win: the wave was wrong about its own mechanism working.
- **N-2** (HELD): FAILED as pre-registered, but note the DIRECTION IS RIGHT: predicted 0.815 is below 1.0, i.e. a fall, against a measured 0.15. The model produces the turnover and undersizes it; the claim said it would produce a RISE, so the prediction held on the pass/fail column and was wrong about the mechanism.
- **N-4** (FALSIFIED): Pre-registered at 45% -- genuinely uncertain, and DIAGNOSTIC either way because Ames measures a low-moisture extrudate and this module has no matrix-transfer term.

## The hold-out panel

| statistic | pre-B8 | physics only | **B8** |
|---|---:|---:|---:|
| gating | 12 / 32 | 7 / 30 | **8 / 30** |
| median \|log10 fold\| (gating) | 0.619 | 1.01 | **0.862** |
| geometric-mean fold (gating) | 13.3x | 21.4x | **17.1x** |
| diagnostic | 7 / 10 | 7 / 16 | **8 / 16** |

⚠️ **THE DENOMINATORS ARE NOT THE SAME ACROSS COLUMNS.** The pre-B8 column is the panel as it stood at HEAD; the two B8 columns run the post-Amendment-17 panel, which RETIRES two gating rows and adds four rows that cannot gate. The retirement accounting immediately below gives both bases on the same scorer so the change is checkable.

### Retirement accounting -- BOTH BASES

| basis | gating |
|---|---:|
| OLD (the two `kang_switch_on_*` rows counted, exactly as they scored) | **8 / 32** |
| NEW (post-retirement) | **8 / 30** |

- RETIRED: `kang_switch_on_MFT` -- would have FAILED at 0.444x.
- RETIRED: `kang_switch_on_FFT` -- would have FAILED at 0.00691x.
- **Rows added back into the gating total: 0.** The three replacement hold-outs are seen-by-extraction and may never gate, so the denominator genuinely falls by 2.

## The both-ways cutover exam

| statistic | pre-B8 | physics only | **B8** |
|---|---:|---:|---:|
| core paired median fold | 42.2x | 41.7x | **24.8x** |
| core geometric-mean fold | 43.3x | 44.2x | **32.1x** |
| core median fold (all) | 40.2x | 40x | **19.1x** |
| worst fold | 3.34e+04x | 3.34e+04x | **3.34e+04x** |
| in band | 4 / 34 | 5 / 34 | **3 / 34** |
| answered / declined | 34 / 6 | 34 / 6 | **34 / 6** |

### By family

| family | pre-B8 median | physics only | **B8 median** | in band (B8) |
|---|---:|---:|---:|---:|
| `acrylamide_180C` | 6.32x | 6.32x | **6.32x** | 1 / 7 |
| `furan_browning_glc_alanine` | 11.9x | 11.9x | **11.9x** | 1 / 5 |
| `matrix_path_lipid` | 1.86e+03x | 1.86e+03x | **1.86e+03x** | 0 / 4 |
| `sulfur_hofmann1998_145C` | 11.5x | 11.7x | **12.2x** | 1 / 10 |
| `sulfur_yiltirak2026_T_ladder` | 193x | 206x | **122x** | 0 / 8 |

### X-5 -- the three families no B8 parameter touches

No acrylamide, furanic or lipid parameter is touched by B8, and the covalent-sink retirement is documentation over a term that already contributed exactly 0.0. Pre-registered as X-5: if any of these moves, the retirement was not the no-op this wave claims and that is a finding, not a rounding error.

**16 points checked across `acrylamide_180C`, `furan_browning_glc_alanine`, `matrix_path_lipid`. VERDICT: BIT-IDENTICAL.**

## The Yiltirak family -- the verdict, point by point

| compound | T (C) | t (min) | pre-B8 | physics only | **B8** | improvement |
|---|---:|---:|---:|---:|---:|---:|
| 2-Furfurylthiol (FFT) | 100 | 240 | 599x | 655x | **496x** | 1.21x |
| 2-Furfurylthiol (FFT) | 110 | 120 | 648x | 713x | **382x** | 1.70x |
| 2-Furfurylthiol (FFT) | 120 | 60 | 469x | 519x | **234x** | 2.01x |
| 2-Furfurylthiol (FFT) | 130 | 30 | 226x | 249x | **149x** | 1.51x |
| 2-Methyl-3-furanthiol (MFT) | 100 | 240 | 38.2x | 38.2x | **8.83x** | 4.33x |
| 2-Methyl-3-furanthiol (MFT) | 110 | 120 | 107x | 108x | **28.1x** | 3.82x |
| 2-Methyl-3-furanthiol (MFT) | 120 | 60 | 154x | 155x | **53.6x** | 2.87x |
| 2-Methyl-3-furanthiol (MFT) | 130 | 30 | 161x | 163x | **95.5x** | 1.69x |

**8 of 8 points improved; 0 entered the band.** Median improvement by compound: 2-Furfurylthiol (FFT) 1.70x, 2-Methyl-3-furanthiol (MFT) 3.82x.

The mechanism appears and it is NOT ENOUGH. Every one of the eight points improves, none enters the band, and the improvement is strongly ASYMMETRIC between the two thiols -- which is the informative part, because the same asymmetry shows up independently on the Wang replacement hold-outs: MFT's turnover reproduces (predicted 0.36 against a measured 0.30, inside the band) and FFT's does not (predicted 0.81 against a measured 0.15). Two scorers, two laboratories, one conclusion: B8's T-structure gives MFT the peak-and-decline the corpus measures and leaves FFT's sink far too weak. The remaining FFT error is therefore localised rather than diffuse, and it is localised on the FURFURAL LIMB -- the leg where Wang measures the precursor SATURATING while the thiol collapses.

## The three replacement hold-outs

All are `seen_diagnostic` and **none gates.** Each carries a written, quantitative prediction made before this scorer existed (`kinetic_core_b8_prereg.md` sec. 3, rows N-1 to N-4).

| row | gates | result | fold |
|---|---|---|---:|
| `wang2026_MFT_peak_and_fall_125_over_115` | **no** | PASS | 1.21x |
| `wang2026_FFT_peak_and_fall_115_over_105` | **no** | **FAIL** | 5.43x |
| `zhai_13C_exogenous_carbon_threshold` | **no** | **FAIL** | --x |
| `ames2001_excess_Ea_class_split` | **no** | **FAIL** | --x |

### `wang2026_MFT_peak_and_fall_125_over_115`

- observed: `0.3`
- predicted: `0.3620415949630796`
- cannot gate because: SEEN-BY-EXTRACTION and disclosed (Amendment 17 clause 3, B8). Printed in k6a sec. 4.3, which this wave was instructed to read. Pre-registered as N-1 with a predicted FAIL in a named direction before the scorer existed. Scored and reported, never counted.
- MFT rises x2.6 over 105->115 C and then FALLS x0.30 over 115->125 C. The model's thiol formation is monotone in temperature once its sink is anchored at 145 C, so a fold BELOW 1.0 here requires a turnover the pre-B8 core could not express at all. THIS IS THE ROW B8's T-STRUCTURE EXISTS FOR.

### `wang2026_FFT_peak_and_fall_115_over_105`

- observed: `0.15`
- predicted: `0.8147067458705215`
- cannot gate because: SEEN-BY-EXTRACTION and disclosed (Amendment 17 clause 3, B8). Pre-registered as N-2. Scored and reported, never counted.
- FFT collapses x0.15 over 105->115 C WHILE ITS OWN PRECURSOR SATURATES -- furfural goes x2.03, x1.10, x1.03 across the same legs. That makes it a SINK observation rather than a supply observation, and it is the second independent one in the corpus (Zhai's 140 C FFT turnover is the first).

### `zhai_13C_exogenous_carbon_threshold`

- observed: `{'100C': 0.0, '120C': 0.2}`
- predicted: `{'100C_ACTZ_proxy': 0.9313120482217044, '120C_ACTZ_proxy': 0.939475689992803}`
- cannot gate because: SEEN-BY-EXTRACTION and disclosed (Amendment 17 clause 3, B8), and additionally OUT OF SCOPE in kind: the three compounds Zhai isotope-traced are thiophenes, which this module declares out of scope, and MFT/FFT were never traced. Pre-registered as N-3 with the scope gap stated in advance. Scored on an in-scope proxy, reported as a SCOPE GAP, never counted.
- ★ SCOPE GAP, AND IT IS THE HONEST ANSWER RATHER THAN A SCORE. The three isotope-traced compounds are thiophenes and this module carries none of them; MFT and FFT -- which it does carry -- were never traced. The proxy scored here is the ADDED-XYLOSE share of 2-acetylthiazole, the only traced compound in scope, and 2-acetylthiazole's own measured isotope pattern is non-monotone (44 / 24 / 46 %) and is called noise-dominated by its own dossier. A pass would therefore be weak evidence and a fail would be weak evidence; what the row genuinely records is that the corpus's best temperature-structure measurement sits on chemistry this module declares out of scope.

### `ames2001_excess_Ea_class_split`

- observed: `{'MFT_excess_sign': '+', 'FFT_excess_sign': '+', 'thiazole_excess_sign': '-'}`
- predicted: `{'MFT': -130.80381387888602, 'FFT': -143.80284264851113, 'ACTZ': 59.084879622697095}`
- cannot gate because: SEEN-BY-EXTRACTION and disclosed (Amendment 17 clause 3, B8). Diagnostic for a second, independent reason: Ames measures a LOW-MOISTURE EXTRUDATE and k6a sec. 5.1 measures a ~2x aqueous-to-matrix gap on exactly these two barriers, with no transfer term in this model. Pre-registered as N-4. Never counted.
- DIRECTIONAL AND DIAGNOSTIC. Scored in this module's own AQUEOUS pot because it has no matrix-transfer term, against a measurement made in a low-moisture starch extrudate where k6a measures the same two barriers ~2x higher than in water. A pass would mean the class split is a property of the network's topology rather than of the medium, which is the interesting reading; a fail localises it to the medium. THREE FURTHER MISMATCHES, all stated in advance and none corrected: (1) the 150 C rung is ABOVE this module's declared 100-145 C window, so the upper leg is an extrapolation; (2) pH 7.5 here is Ames' FEED pH -- the paper MEASURES the extrudate running 1.4-2.6 units below it, and the module is given the feed value because the offset is a property of the extruder and not of this pot; (3) Ames' own nominal ladder is not the real one -- measured die temperatures run -7 to +14 C off target and the xylose pH-7.5 upper leg is 153->174 C, so the true axis is 21 K where the nominal is 30 K and excess energies on that leg are ~40% lower than the nominal ladder implies.

## The T-structure, as fitted

- `Ea_decay_thiol_sink` = **102.00 kJ/mol**, in the declared band [7.0, 102.0] around a prior centre of 60.2 -- **ON THE CEILING**.
- `Ea_decay_carbonyl_sink` = **166.67 kJ/mol** in [20.0, 250.0] -- **interior** (UNCHANGED band; K6a measured no carbonyl-sink barrier).
- Measured barriers the fit cannot move: `k_dimer_mft` = 122.2, `k_dimer_fft` = 122.2, `k_arp_dpo` = 85.7, `k_arp_tdp` = 85.7, `k_cys_thermal` = 55.1.

**REFUTED:** Ea_decay_thiol_sink = 248.0 (B2.2), 216.1 (B2.3), 212.9-218.1 (B2.4). Gigl 2021 measures the covalent-capture channel at k(333)/k(279) = 67.2; 248 kJ/mol predicts 3.4e7 -- wrong by 5.1e5x.

*Objective: 62 declared FIT rows (4 new in B8), 23 free of 48 parameters. TOTAL COST IS NOT COMPARABLE TO B2.4's: B8 scores four rows B2.4 did not. `sum_r2_level_shared_with_b2_4` is the like-for-like comparator and is summed over exactly the non-pH rows both waves scored.*
