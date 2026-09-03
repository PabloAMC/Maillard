# Core directional scores (the kinetic core on the directional claims panel)

Panel `docs/validation/directional_claims_panel.yml` (69 claims, flat tolerance 5 %). Nothing in the core was tuned to this panel.

* **headline (strictly independent, evaluable): 17/26 (65%)**; 26 independent claims not evaluable
* independent, excluding pH and water activity: 13/21 (62%); pH and water activity alone: 4/5 (80%)
* all claims (independent + fit-adjacent + fit-system overlap): 19/37 (51%); 32 not evaluable
* misses where the lane carries no term for the moved axis (identical predictions): 4
* not evaluable, by reason: the claim carries no runnable conditions (prose-only) (16); refused by the engine (4); a predicted concentration is zero; no direction is defined (2); arm 'processed, 140 C' refused (2); arm 'pH 4.5' refused (1); arm 'D-ribose' refused (1); arm 'pH 4' refused (1); arm 'T4 profile, 160 C max' refused (1); arm 'hydrolysate + xylose (no cysteine), 120 C' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone, 1 mmol each in 50 mL, 145 C 20  (1); arm 'norfuraneol + H2S, 1 mmol each in 50 mL, 145 C 20 min pH 5.0' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone at pH 7.0' refused (1)

## Per category (strictly independent claims)

| category | agree | evaluable | rate | not evaluable | misses |
|---|---|---|---|---|---|
| additive_cysteine | 2 | 3 | 0.67 | 2 | CYS-02 |
| lipid_lane | 0 | 0 | - | 1 | - |
| matrix_identity | 0 | 0 | - | 1 | - |
| moisture_aw | 0 | 0 | - | 4 | - |
| ph | 4 | 5 | 0.80 | 5 | MOT-01 |
| process_heating | 0 | 0 | - | 1 | - |
| ranking | 0 | 1 | 0.00 | 2 | MOT-03 |
| scope | 0 | 0 | - | 3 | - |
| sugar_identity | 4 | 8 | 0.50 | 6 | SUG-03, SUG-12, HOF-02, HOF-03 |
| temperature | 5 | 7 | 0.71 | 1 | TEMP-01, TEMP-05 |
| time | 2 | 2 | 1.00 | 0 | - |

## Per category (all claims)

| category | agree | evaluable | rate | not evaluable |
|---|---|---|---|---|
| additive_cysteine | 2 | 3 | 0.67 | 2 |
| lipid_lane | 0 | 3 | 0.00 | 1 |
| matrix_identity | 0 | 2 | 0.00 | 1 |
| moisture_aw | 0 | 0 | - | 4 |
| ph | 4 | 6 | 0.67 | 6 |
| process_heating | 0 | 0 | - | 4 |
| ranking | 0 | 2 | 0.00 | 4 |
| scope | 0 | 0 | - | 3 |
| sugar_identity | 6 | 12 | 0.50 | 6 |
| temperature | 5 | 7 | 0.71 | 1 |
| time | 2 | 2 | 1.00 | 0 |

## Claims

| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |
|---|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | fit_adjacent | MFT | A>B | **agree** | sulfur | 76.3, 0 |  |
| SUG-02 | sugar_identity | fit_adjacent | FFT | A>B | **agree** | sulfur | 57.7, 52.3 |  |
| SUG-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 52.3, 52.3 | identical predictions |
| SUG-04 | sugar_identity | fit_adjacent | MFT | A>B | **disagree** | sulfur | 76.3, 76.3 | identical predictions |
| SUG-05 | sugar_identity | fit_system_overlap | Furfural | A>B | **disagree** | sulfur | 136, 147 |  |
| SUG-06 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-01 | ph | fit_system_overlap | MFT | A>B | **disagree** | sulfur | 46.8, 83.1 |  |
| PH-02 | ph | fit_system_overlap | bis(2-methyl-3-furyl) disulfide | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| PH-03 | ph | independent | FFT | decreasing | **agree** | sulfur | 65.9, 53.8, 10.9 |  |
| PH-04 | ph | independent | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | arm 'pH 4.5' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either: th |
| PH-05 | ph | independent | Furfural | decreasing | **agree** | sulfur | 146, 131, 58.1 |  |
| TEMP-01 | temperature | independent | Acrylamide | decreasing | **disagree** | acrylamide | 21.1, 454, 814 |  |
| TEMP-02 | temperature | independent | Acrylamide | A>B | **agree** | acrylamide | 21.1, 0 |  |
| TEMP-03 | temperature | independent | HMF | increasing | **agree** | acrylamide | 77.7, 399, 577 |  |
| AW-01 | moisture_aw | independent | HMF | decreasing | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in WATER ACTIVITY and no core lane carries an a_w term; the model would return identical arms and call it  |
| AW-02 | moisture_aw | independent | Acrylamide | peak | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in WATER ACTIVITY and no core lane carries an a_w term; the model would return identical arms and call it  |
| CYS-01 | additive_cysteine | independent | FFT | A>B | **agree** | sulfur | 47, 0 |  |
| CYS-02 | additive_cysteine | independent | HMF | A>B | **disagree** | trunk | 676, 674 |  |
| SCOPE-01 | scope | independent | CEL | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SCOPE-02 | scope | independent | Methional | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-01 | matrix_identity | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-02 | matrix_identity | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 0.267, 0.339 |  |
| MAT-03 | matrix_identity | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 0.339, 0.339 | identical predictions; the lane has no term for the moved axis |
| PROC-01 | process_heating | fit_adjacent | Hexanal | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-02 | process_heating | fit_adjacent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | arm 'processed, 140 C' refused: UNREPRESENTED TARGETS: 2-Pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate a |
| PROC-03 | process_heating | fit_adjacent | Nonanal | A>B | **not_evaluable** | - | - | arm 'processed, 140 C' refused: UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydrop |
| PROC-04 | process_heating | independent | 2,5-Dimethylpyrazine | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PROC-05 | ranking | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| TIME-01 | time | independent | HMF | increasing | **agree** | trunk | 14.9, 676, 2.39e+03 |  |
| LIP-01 | lipid_lane | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| LIP-02 | lipid_lane | fit_system_overlap | Hexanal | increasing | **disagree** | lipid | 22.8, 88.6, 88.6 |  |
| LIP-03 | lipid_lane | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-07 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 1.42e+03, 402 |  |
| SUG-08 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 10.1, 0 |  |
| SUG-09 | sugar_identity | independent | FFT | A>B | **agree** | sulfur | 36.9, 33.5 |  |
| SUG-10 | sugar_identity | independent | 2,5-Dimethylpyrazine | decreasing | **not_evaluable** | - | - | arm 'D-ribose' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either:  |
| SUG-11 | sugar_identity | independent | FFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-12 | sugar_identity | independent | HMF | A>B | **disagree** | trunk | 697, 937 |  |
| SUG-13 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 2.09, 0.405 |  |
| PH-06 | ph | independent | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | arm 'pH 4' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either: the  |
| PH-07 | ph | independent | Furfural | flat | **agree** | sulfur | 2.17, 2.09, 2.08 |  |
| TEMP-04 | temperature | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'T4 profile, 160 C max' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound  |
| TEMP-05 | temperature | independent | HMF | increasing | **disagree** | trunk | 684, 644, 654 |  |
| TEMP-06 | temperature | independent | Furfural | A>B | **agree** | sulfur | 2.09, 0.000551 |  |
| AW-03 | moisture_aw | independent | HMF | peak | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in WATER ACTIVITY and no core lane carries an a_w term; the model would return identical arms and call it  |
| CYS-03 | additive_cysteine | independent | MFT | A>B | **agree** | sulfur | 12.5, 0 |  |
| CYS-04 | additive_cysteine | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'hydrolysate + xylose (no cysteine), 120 C' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unr |
| LIP-04 | lipid_lane | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 63.9, 63.9 | identical predictions; the lane has no term for the moved axis |
| ACR-01 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 1.21, 2.83, 2.42 |  |
| ACR-02 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 12.1, 842, 3.43 |  |
| TIME-02 | time | independent | HMF | increasing | **agree** | trunk | 40.5, 684, 937 |  |
| SCOPE-03 | scope | independent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOF-01 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HOF-02 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 1.67e+03, 2.3e+03 |  |
| HOF-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 0.000911, 0.297 |  |
| HOF-04 | moisture_aw | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-01 | ph | independent | MFT | A>B | **disagree** | sulfur | 190, 505 |  |
| MOT-02 | ph | independent | FFT | A>B | **agree** | sulfur | 1.37e+03, 1.08e+03 |  |
| MOT-03 | ranking | independent | MFT | A>B | **disagree** | sulfur | 505, 505 | identical predictions |
| MOT-04 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-05 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-06 | ph | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-01 | ranking | fit_adjacent | MFT | A>B | **not_evaluable** | - | - | arm 'hydroxyacetaldehyde + mercapto-2-propanone, 1 mmol each in 50 mL, 145 C 20 min pH 5.0' refused: UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-prop |
| HOX-02 | ranking | fit_adjacent | MFT | A>B | **not_evaluable** | - | - | arm 'norfuraneol + H2S, 1 mmol each in 50 mL, 145 C 20 min pH 5.0' refused: UNMAPPED PRECURSORS 'Hydrogen sulfide': not a species in any core lane. The core is  |
| HOX-03 | ph | independent | MFT | A>B | **not_evaluable** | - | - | arm 'hydroxyacetaldehyde + mercapto-2-propanone at pH 7.0' refused: UNMAPPED PRECURSORS 'Hydroxyacetaldehyde', 'Mercapto-2-propanone': not a species in any core |
| HOX-04 | ranking | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-05 | ph | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-06 | additive_cysteine | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-07 | ranking | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
