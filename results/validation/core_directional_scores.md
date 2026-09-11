# Core directional scores (the kinetic core on the directional claims panel)

Panel `docs/validation/directional_claims_panel.yml` (92 claims, flat tolerance 5 %). Nothing in the core was tuned to this panel.

* **headline (strictly independent, evaluable): 25/43 (58%)**; 28 independent claims not evaluable
* independent, excluding pH and water activity: 16/31 (52%); pH and water activity alone: 9/12 (75%)
* all claims (independent + fit-adjacent + fit-system overlap): 29/59 (49%); 33 not evaluable
* misses where the lane carries no term for the moved axis (identical predictions): 4
* not evaluable, by reason: the claim carries no runnable conditions (prose-only) (18); refused by the engine (4); a predicted concentration is zero; no direction is defined (4); arm 'pH 4.5' refused (1); arm 'D-ribose' refused (1); arm 'hydrolysate + xylose (no cysteine), 120 C' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone, 1 mmol each in 50 mL, 145 C 20  (1); arm 'norfuraneol + H2S, 1 mmol each in 50 mL, 145 C 20 min pH 5.0' refused (1); arm 'hydroxyacetaldehyde + mercapto-2-propanone at pH 7.0' refused (1); the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, d (1)

## Per category (strictly independent claims)

| category | agree | evaluable | rate | not evaluable | misses |
|---|---|---|---|---|---|
| additive_cysteine | 2 | 3 | 0.67 | 2 | CYS-02 |
| lipid_lane | 0 | 0 | - | 1 | - |
| matrix_identity | 0 | 0 | - | 1 | - |
| moisture_aw | 1 | 2 | 0.50 | 2 | AW-01 |
| ph | 8 | 10 | 0.80 | 4 | MOT-01, CER07-PH-01 |
| process_heating | 0 | 0 | - | 1 | - |
| ranking | 0 | 1 | 0.00 | 2 | MOT-03 |
| scope | 0 | 0 | - | 3 | - |
| sugar_identity | 4 | 10 | 0.40 | 7 | SUG-03, SUG-12, HOF-02, HOF-03, DIC-01, DIC-03 |
| temperature | 7 | 11 | 0.64 | 2 | TEMP-01, TEMP-05, YIL-01, YIL-02 |
| time | 3 | 6 | 0.50 | 3 | WANG22-T-01, WANG22-T-03, WANG22-T-04 |

## Per category (all claims)

| category | agree | evaluable | rate | not evaluable |
|---|---|---|---|---|
| additive_cysteine | 2 | 3 | 0.67 | 2 |
| lipid_lane | 0 | 3 | 0.00 | 1 |
| matrix_identity | 0 | 2 | 0.00 | 1 |
| moisture_aw | 2 | 3 | 0.67 | 2 |
| ph | 9 | 14 | 0.64 | 4 |
| process_heating | 0 | 0 | - | 4 |
| ranking | 0 | 2 | 0.00 | 4 |
| scope | 0 | 0 | - | 3 |
| sugar_identity | 6 | 14 | 0.43 | 7 |
| temperature | 7 | 11 | 0.64 | 2 |
| time | 3 | 7 | 0.43 | 3 |

## Claims

| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |
|---|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | fit_adjacent | MFT | A>B | **agree** | sulfur | 76.3, 0 |  |
| SUG-02 | sugar_identity | fit_adjacent | FFT | A>B | **agree** | sulfur | 57.7, 52.2 |  |
| SUG-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 52.2, 52.2 | identical predictions |
| SUG-04 | sugar_identity | fit_adjacent | MFT | A>B | **disagree** | sulfur | 76.3, 76.3 | identical predictions |
| SUG-05 | sugar_identity | fit_system_overlap | Furfural | A>B | **disagree** | sulfur | 136, 147 |  |
| SUG-06 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-01 | ph | fit_system_overlap | MFT | A>B | **disagree** | sulfur | 46.8, 83.1 |  |
| PH-02 | ph | fit_system_overlap | bis(2-methyl-3-furyl) disulfide | A>B | **disagree** | sulfur | 0.0452, 0.111 |  |
| PH-03 | ph | independent | FFT | decreasing | **agree** | sulfur | 65.9, 53.7, 10.9 |  |
| PH-04 | ph | independent | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | arm 'pH 4.5' refused: PYRAZINE TARGETS '2,5-Dimethylpyrazine' (wave B18) run on the trunk lane only: the sulfur lane's network keeps the topology its fit was ru |
| PH-05 | ph | independent | Furfural | decreasing | **agree** | sulfur | 146, 131, 58.1 |  |
| TEMP-01 | temperature | independent | Acrylamide | decreasing | **disagree** | acrylamide | 13.6, 357, 1.04e+03 |  |
| TEMP-02 | temperature | independent | Acrylamide | A>B | **agree** | acrylamide | 13.6, 0 |  |
| TEMP-03 | temperature | independent | HMF | increasing | **agree** | acrylamide | 78, 403, 585 |  |
| AW-01 | moisture_aw | independent | HMF | decreasing | **disagree** | trunk | 1.1e+03, 1.22e+03, 544 |  |
| AW-02 | moisture_aw | independent | Acrylamide | peak | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in WATER ACTIVITY (0.30 vs 0.60) and the acrylamide lane's a_w term is measured only inside 0.34-0.99 (De  |
| CYS-01 | additive_cysteine | independent | FFT | A>B | **agree** | sulfur | 47, 0 |  |
| CYS-02 | additive_cysteine | independent | HMF | A>B | **disagree** | trunk | 544, 674 |  |
| SCOPE-01 | scope | independent | CEL | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SCOPE-02 | scope | independent | Methional | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-01 | matrix_identity | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-02 | matrix_identity | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 0.267, 0.339 |  |
| MAT-03 | matrix_identity | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 0.339, 0.339 | identical predictions; the lane has no term for the moved axis |
| PROC-01 | process_heating | fit_adjacent | Hexanal | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-02 | process_heating | fit_adjacent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-03 | process_heating | fit_adjacent | Nonanal | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-04 | process_heating | independent | 2,5-Dimethylpyrazine | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PROC-05 | ranking | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| TIME-01 | time | independent | HMF | increasing | **agree** | trunk | 10.3, 544, 2.1e+03 |  |
| LIP-01 | lipid_lane | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.8, 22.8 | identical predictions; the lane has no term for the moved axis |
| LIP-02 | lipid_lane | fit_system_overlap | Hexanal | increasing | **disagree** | lipid | 22.8, 88.6, 88.6 |  |
| LIP-03 | lipid_lane | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-07 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 1.42e+03, 402 |  |
| SUG-08 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 10.1, 0 |  |
| SUG-09 | sugar_identity | independent | FFT | A>B | **agree** | sulfur | 36.9, 33.5 |  |
| SUG-10 | sugar_identity | independent | 2,5-Dimethylpyrazine | decreasing | **not_evaluable** | - | - | arm 'D-ribose' refused: PYRAZINE TARGETS '2,5-Dimethylpyrazine' (wave B18) run on the trunk lane only: the sulfur lane's network keeps the topology its fit was  |
| SUG-11 | sugar_identity | independent | FFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-12 | sugar_identity | independent | HMF | A>B | **disagree** | trunk | 881, 1.37e+03 |  |
| SUG-13 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 2.09, 0.405 |  |
| PH-06 | ph | independent | 2,5-Dimethylpyrazine | increasing | **agree** | trunk | 0.000306, 6.94, 22.7 |  |
| PH-07 | ph | independent | Furfural | flat | **agree** | sulfur | 2.17, 2.09, 2.08 |  |
| TEMP-04 | temperature | independent | 2,5-Dimethylpyrazine | A>B | **agree** | trunk | 0.00159, 4.37e-06 |  |
| TEMP-05 | temperature | independent | HMF | increasing | **disagree** | trunk | 1.04e+03, 757, 680 |  |
| TEMP-06 | temperature | independent | Furfural | A>B | **agree** | sulfur | 2.09, 0.000551 |  |
| AW-03 | moisture_aw | independent | HMF | peak | **agree** | trunk | 1.1e+03, 1.24e+03, 647 |  |
| CYS-03 | additive_cysteine | independent | MFT | A>B | **agree** | sulfur | 12.5, 0 |  |
| CYS-04 | additive_cysteine | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'hydrolysate + xylose (no cysteine), 120 C' refused: PYRAZINE TARGETS '2,5-Dimethylpyrazine' (wave B18) run on the trunk lane only: the sulfur lane's networ |
| LIP-04 | lipid_lane | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 63.9, 63.9 | identical predictions; the lane has no term for the moved axis |
| ACR-01 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 0.716, 1.67, 1.43 |  |
| ACR-02 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 7.64, 935, 4.37 |  |
| TIME-02 | time | independent | HMF | increasing | **agree** | trunk | 55.3, 1.04e+03, 1.37e+03 |  |
| SCOPE-03 | scope | independent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOF-01 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HOF-02 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 1.67e+03, 2.29e+03 |  |
| HOF-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 0.000911, 0.297 |  |
| HOF-04 | moisture_aw | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-01 | ph | independent | MFT | A>B | **disagree** | sulfur | 190, 505 |  |
| MOT-02 | ph | independent | FFT | A>B | **agree** | sulfur | 1.36e+03, 1.07e+03 |  |
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
| YIL-01 | temperature | independent | MFT | decreasing | **disagree** | sulfur | 79.5, 121, 163, 197 |  |
| YIL-02 | temperature | independent | FFT | A>B | **disagree** | sulfur | 1.21e+03, 2.3e+03 |  |
| WANG-01 | temperature | independent | MFT | peak | **not_evaluable** | sulfur | 67.4, 574, 1.47e+03, 2.71e+03, 4.09e+03 | the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, disagree at ph 9); the source does not state it |
| WANG-02 | temperature | independent | FFT | peak | **agree** | sulfur | 149, 3.63e+03, 2.32e+04, 3.55e+04, 2.29e+04 |  |
| MENG-01 | temperature | independent | MFT, FFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-05 | moisture_aw | fit_adjacent | Acrylamide | flat | **agree** | acrylamide | 640, 640 | identical predictions |
| PH-ACR-01 | ph | fit_adjacent | Acrylamide | A>B | **agree** | acrylamide | 5.31e+04, 2.73e+04 |  |
| DIC-01 | sugar_identity | independent | 3-deoxyglucosone, glucosone, methylglyoxal, glyoxal, diacetyl | ranking | **disagree** | trunk | 1.19e+04, 12, 759, 38.1, 0.0285 |  |
| DIC-02 | sugar_identity | independent | glucosone, glyoxal, methylglyoxal, diacetyl | ranking | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| DIC-03 | sugar_identity | independent | 3-deoxyglucosone, glucosone, glyoxal, methylglyoxal | ranking | **disagree** | trunk | 4.7e+04, 119, 161, 2.03e+03 |  |
| RIB-T-01 | time | independent | MFT | decreasing | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| RIB-T-02 | time | independent | FFT | flat | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HEX-T-01 | time | independent | FFT | increasing | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| SCH-T-01 | time | fit_system_overlap | MFT | increasing | **disagree** | sulfur | 65.9, 130, 36.4, 16.2 |  |
| WANG22-T-01 | time | independent | MFT | increasing | **disagree** | sulfur | 185, 498, 391 |  |
| WANG22-T-02 | time | independent | FFT | increasing | **agree** | sulfur | 2.74e+03, 2.01e+04, 3.57e+04 |  |
| WANG22-T-03 | time | independent | MFT | peak | **disagree** | sulfur | 2.51e+03, 693, 32.7, 2.77 |  |
| WANG22-T-04 | time | independent | FFT | peak | **disagree** | sulfur | 1.72e+04, 7.05e+03, 732, 61.5 |  |
| WHI-PH-01 | ph | fit_adjacent | MFT | A>B | **disagree** | sulfur | 4.29e+03, 9.11e+03 |  |
| CER07-PH-01 | ph | independent | FFT | decreasing | **disagree** | sulfur | 2.46e+03, 310, 0.0189, 0.104 |  |
| CER07-PH-02 | ph | independent | furfural | decreasing | **agree** | sulfur | 1.03e+03, 830, 89.9 |  |
| MOT02-PH-01 | ph | independent | MFT | A>B | **agree** | sulfur | 404, 374 |  |
| MOT02-PH-02 | ph | independent | FFT | A>B | **agree** | sulfur | 1.47e+03, 183 |  |
