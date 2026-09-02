# Core directional scores (the kinetic core on the directional claims panel)

Panel `docs/validation/directional_claims_panel.yml` (69 claims, flat tolerance 5 %). Nothing in the core was tuned to this panel.

* **headline (strictly independent, evaluable): 18/30 (60%)**; 22 independent claims not evaluable
* independent, excluding pH and water activity: 15/22 (68%); pH and water activity alone: 3/8 (38%)
* all claims (independent + fit-adjacent + fit-system overlap): 22/42 (52%); 27 not evaluable
* misses where the lane carries no term for the moved axis (identical predictions): 7
* not evaluable, by reason: the claim carries no runnable conditions (prose-only) (16); arm 'processed, 140 C' refused (2); a predicted concentration is zero; no direction is defined (1); arm 'pH 4.5' refused (1); arm 'D-ribose' refused (1); arm 'pH 4' refused (1); arm 'T4 profile, 160 C max' refused (1); arm 'hydrolysate + xylose (no cysteine), 120 C' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone, 1 mmol each in 50 mL, 145 C 20  (1); arm 'norfuraneol + H2S, 1 mmol each in 50 mL, 145 C 20 min pH 5.0' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone at pH 7.0' refused (1)

## Per category (strictly independent claims)

| category | agree | evaluable | rate | not evaluable | misses |
|---|---|---|---|---|---|
| additive_cysteine | 3 | 3 | 1.00 | 2 | - |
| lipid_lane | 0 | 0 | - | 1 | - |
| matrix_identity | 0 | 0 | - | 1 | - |
| moisture_aw | 0 | 3 | 0.00 | 1 | AW-01, AW-02, AW-03 |
| ph | 3 | 5 | 0.60 | 5 | PH-07, MOT-01 |
| process_heating | 0 | 0 | - | 1 | - |
| ranking | 0 | 1 | 0.00 | 2 | MOT-03 |
| scope | 0 | 0 | - | 3 | - |
| sugar_identity | 5 | 9 | 0.56 | 5 | SUG-03, SUG-12, HOF-01, HOF-03 |
| temperature | 5 | 7 | 0.71 | 1 | TEMP-01, TEMP-05 |
| time | 2 | 2 | 1.00 | 0 | - |

## Per category (all claims)

| category | agree | evaluable | rate | not evaluable |
|---|---|---|---|---|
| additive_cysteine | 3 | 3 | 1.00 | 2 |
| lipid_lane | 0 | 3 | 0.00 | 1 |
| matrix_identity | 0 | 2 | 0.00 | 1 |
| moisture_aw | 0 | 3 | 0.00 | 1 |
| ph | 3 | 6 | 0.50 | 6 |
| process_heating | 1 | 1 | 1.00 | 3 |
| ranking | 0 | 2 | 0.00 | 4 |
| scope | 0 | 0 | - | 3 |
| sugar_identity | 8 | 13 | 0.62 | 5 |
| temperature | 5 | 7 | 0.71 | 1 |
| time | 2 | 2 | 1.00 | 0 |

## Claims

| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |
|---|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | fit_adjacent | MFT | A>B | **agree** | sulfur | 72.8, 0.185 |  |
| SUG-02 | sugar_identity | fit_adjacent | FFT | A>B | **agree** | sulfur | 78.9, 3.81 |  |
| SUG-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 3.81, 3.81 | identical predictions |
| SUG-04 | sugar_identity | fit_adjacent | MFT | A>B | **disagree** | sulfur | 72.8, 72.8 | identical predictions |
| SUG-05 | sugar_identity | fit_system_overlap | Furfural | A>B | **agree** | sulfur | 157, 4.41 |  |
| SUG-06 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-01 | ph | fit_system_overlap | MFT | A>B | **disagree** | sulfur | 44.6, 79.4 |  |
| PH-02 | ph | fit_system_overlap | bis(2-methyl-3-furyl) disulfide | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| PH-03 | ph | independent | FFT | decreasing | **agree** | sulfur | 90.2, 73.6, 15 |  |
| PH-04 | ph | independent | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | arm 'pH 4.5' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either: th |
| PH-05 | ph | independent | Furfural | decreasing | **agree** | sulfur | 169, 152, 67.5 |  |
| TEMP-01 | temperature | independent | Acrylamide | decreasing | **disagree** | acrylamide | 21.1, 454, 814 |  |
| TEMP-02 | temperature | independent | Acrylamide | A>B | **agree** | acrylamide | 21.1, 2.63e-16 |  |
| TEMP-03 | temperature | independent | HMF | increasing | **agree** | acrylamide | 77.7, 399, 577 |  |
| AW-01 | moisture_aw | independent | HMF | decreasing | **disagree** | trunk | 676, 676, 676 | identical predictions; the lane has no term for the moved axis |
| AW-02 | moisture_aw | independent | Acrylamide | peak | **disagree** | acrylamide | 454, 454, 454 | identical predictions; the lane has no term for the moved axis |
| CYS-01 | additive_cysteine | independent | FFT | A>B | **agree** | sulfur | 3.45, 0 |  |
| CYS-02 | additive_cysteine | independent | HMF | A>B | **agree** | trunk | 676, 286 |  |
| SCOPE-01 | scope | independent | CEL | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SCOPE-02 | scope | independent | Methional | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-01 | matrix_identity | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-02 | matrix_identity | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 0.267, 0.339 |  |
| MAT-03 | matrix_identity | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 0.339, 0.339 | identical predictions; the lane has no term for the moved axis |
| PROC-01 | process_heating | fit_adjacent | Hexanal | A>B | **agree** | lipid | 22.8, 0.000886 |  |
| PROC-02 | process_heating | fit_adjacent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | arm 'processed, 140 C' refused: UNREPRESENTED TARGETS: 2-Pentylfuran -- The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's six-product slate a |
| PROC-03 | process_heating | fit_adjacent | Nonanal | A>B | **not_evaluable** | - | - | arm 'processed, 140 C' refused: UNREPRESENTED TARGETS: nonanal -- the lipid lane exists and nonanal is a species in it, but its ONLY parent is the OLEATE hydrop |
| PROC-04 | process_heating | independent | 2,5-Dimethylpyrazine | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PROC-05 | ranking | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| TIME-01 | time | independent | HMF | increasing | **agree** | trunk | 14.9, 676, 2.39e+03 |  |
| LIP-01 | lipid_lane | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| LIP-02 | lipid_lane | fit_system_overlap | Hexanal | increasing | **disagree** | lipid | 22.8, 88.6, 88.6 |  |
| LIP-03 | lipid_lane | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-07 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 1.34e+03, 59.4 |  |
| SUG-08 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 9.65, 0.0508 |  |
| SUG-09 | sugar_identity | independent | FFT | A>B | **agree** | sulfur | 42.5, 8.44 |  |
| SUG-10 | sugar_identity | independent | 2,5-Dimethylpyrazine | decreasing | **not_evaluable** | - | - | arm 'D-ribose' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either:  |
| SUG-11 | sugar_identity | independent | FFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-12 | sugar_identity | independent | HMF | A>B | **disagree** | trunk | 697, 937 |  |
| SUG-13 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 0.0805, 0.0225 |  |
| PH-06 | ph | independent | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | arm 'pH 4' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound list either: the  |
| PH-07 | ph | independent | Furfural | flat | **disagree** | sulfur | 0.0852, 0.0804, 0.0799 |  |
| TEMP-04 | temperature | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'T4 profile, 160 C max' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unrepresented-compound  |
| TEMP-05 | temperature | independent | HMF | increasing | **disagree** | trunk | 684, 644, 654 |  |
| TEMP-06 | temperature | independent | Furfural | A>B | **agree** | sulfur | 0.0804, 0.000551 |  |
| AW-03 | moisture_aw | independent | HMF | peak | **disagree** | trunk | 676, 676, 676 | identical predictions; the lane has no term for the moved axis |
| CYS-03 | additive_cysteine | independent | MFT | A>B | **agree** | sulfur | 11.9, 1.43e-30 |  |
| CYS-04 | additive_cysteine | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'hydrolysate + xylose (no cysteine), 120 C' refused: UNREPRESENTED TARGETS: 2,5-Dimethylpyrazine -- not a species in any core lane, and not on the named unr |
| LIP-04 | lipid_lane | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 63.9, 63.9 | identical predictions; the lane has no term for the moved axis |
| ACR-01 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 1.21, 2.83, 2.42 |  |
| ACR-02 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 12.1, 842, 3.43 |  |
| TIME-02 | time | independent | HMF | increasing | **agree** | trunk | 40.5, 684, 937 |  |
| SCOPE-03 | scope | independent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOF-01 | sugar_identity | independent | MFT | A>B | **disagree** | sulfur | 10.5, 33.1 |  |
| HOF-02 | sugar_identity | independent | FFT | A>B | **agree** | sulfur | 338, 237 |  |
| HOF-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 0.00744, 0.023 |  |
| HOF-04 | moisture_aw | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-01 | ph | independent | MFT | A>B | **disagree** | sulfur | 181, 481 |  |
| MOT-02 | ph | independent | FFT | A>B | **agree** | sulfur | 1.83e+03, 1.45e+03 |  |
| MOT-03 | ranking | independent | MFT | A>B | **disagree** | sulfur | 481, 481 | identical predictions |
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
