# Core directional scores (the kinetic core on the directional claims panel)

Panel `docs/validation/directional_claims_panel.yml` (92 claims, flat tolerance 5 %). Nothing in the core was tuned to this panel.

* **headline (strictly independent, evaluable): 17/31 (55%)**; 22 independent claims not evaluable
* independent, excluding pH and water activity: 12/23 (52%); pH and water activity alone: 5/8 (62%)
* all claims (independent + fit-adjacent + fit-system overlap): 26/50 (52%); 42 not evaluable
* misses where the lane carries no term for the moved axis (identical predictions): 3
* not evaluable, by reason: the claim carries no runnable conditions (prose-only) (34); refused by the engine (3); a predicted concentration is zero; no direction is defined (3); arm 'hydrolysate + xylose (no cysteine), 120 C' refused (1); the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, d (1)

## Per category (strictly independent claims)

| category | agree | evaluable | rate | not evaluable | misses |
|---|---|---|---|---|---|
| additive_cysteine | 1 | 1 | 1.00 | 2 | - |
| lipid_lane | 0 | 0 | - | 1 | - |
| matrix_identity | 0 | 0 | - | 1 | - |
| moisture_aw | 1 | 1 | 1.00 | 1 | - |
| ph | 4 | 7 | 0.57 | 2 | PH-07, HOX-03, CER07-PH-01 |
| process_heating | 0 | 0 | - | 1 | - |
| ranking | 0 | 1 | 0.00 | 2 | MOT-03 |
| scope | 0 | 0 | - | 3 | - |
| sugar_identity | 3 | 10 | 0.30 | 5 | SUG-03, SUG-09, SUG-12, HOF-02, HOF-03, DIC-01, DIC-03 |
| temperature | 6 | 9 | 0.67 | 2 | TEMP-01, TEMP-05, YIL-01 |
| time | 2 | 2 | 1.00 | 2 | - |

## Per category (all claims)

| category | agree | evaluable | rate | not evaluable |
|---|---|---|---|---|
| additive_cysteine | 2 | 2 | 1.00 | 3 |
| lipid_lane | 0 | 1 | 0.00 | 3 |
| matrix_identity | 0 | 2 | 0.00 | 1 |
| moisture_aw | 1 | 1 | 1.00 | 4 |
| ph | 9 | 14 | 0.64 | 4 |
| process_heating | 0 | 0 | - | 4 |
| ranking | 1 | 4 | 0.25 | 2 |
| scope | 0 | 0 | - | 3 |
| sugar_identity | 5 | 14 | 0.36 | 7 |
| temperature | 6 | 9 | 0.67 | 4 |
| time | 2 | 3 | 0.67 | 7 |

## Claims

| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |
|---|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | fit_adjacent | MFT | A>B | **agree** | sulfur | 1.68, 0 |  |
| SUG-02 | sugar_identity | fit_adjacent | FFT | A>B | **disagree** | sulfur | 2.56, 888 |  |
| SUG-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 82.3, 82.3 | identical predictions |
| SUG-04 | sugar_identity | fit_adjacent | MFT | A>B | **disagree** | sulfur | 1.68, 1.68 | identical predictions |
| SUG-05 | sugar_identity | fit_system_overlap | Furfural | A>B | **agree** | sulfur | 1.02e+03, 696 |  |
| SUG-06 | sugar_identity | retired | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-01 | ph | fit_system_overlap | MFT | A>B | **agree** | sulfur | 158, 9.97 |  |
| PH-02 | ph | fit_system_overlap | bis(2-methyl-3-furyl) disulfide | A>B | **agree** | sulfur | 10.5, 3.47 |  |
| PH-03 | ph | not_independent | FFT | decreasing | **disagree** | sulfur | 4.73e+03, 0.613, 159 |  |
| PH-04 | ph | retired | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-05 | ph | retired | Furfural | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| TEMP-01 | temperature | independent | Acrylamide | decreasing | **disagree** | acrylamide | 13.6, 357, 1.04e+03 |  |
| TEMP-02 | temperature | independent | Acrylamide | A>B | **agree** | acrylamide | 13.6, 0 |  |
| TEMP-03 | temperature | retired | HMF | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-01 | moisture_aw | retired | HMF | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-02 | moisture_aw | retired | Acrylamide | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| CYS-01 | additive_cysteine | not_independent | FFT | A>B | **agree** | sulfur | 76.8, 0 |  |
| CYS-02 | additive_cysteine | retired | HMF | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SCOPE-01 | scope | independent | CEL | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SCOPE-02 | scope | independent | Methional | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-01 | matrix_identity | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAT-02 | matrix_identity | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 0.267, 0.339 |  |
| MAT-03 | matrix_identity | fit_system_overlap | Hexanal | A>B | **disagree** | lipid | 0.339, 0.339 | identical predictions; the lane has no term for the moved axis |
| PROC-01 | process_heating | fit_adjacent | Hexanal | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-02 | process_heating | fit_adjacent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-03 | process_heating | fit_adjacent | Nonanal | A>B | **not_evaluable** | - | - | refused by the engine: REFUSED -- the two arms differ in pH and the resolved lane(s) (lipid) carry NO pH term by declaration; the model would return identical a |
| PROC-04 | process_heating | independent | 2,5-Dimethylpyrazine | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PROC-05 | ranking | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.9, 22.9 | identical predictions; the lane has no term for the moved axis |
| TIME-01 | time | retired | HMF | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| LIP-01 | lipid_lane | fit_adjacent | Hexanal | A>B | **disagree** | lipid | 22.9, 22.9 | identical predictions; the lane has no term for the moved axis |
| LIP-02 | lipid_lane | retired | Hexanal | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| LIP-03 | lipid_lane | independent | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-07 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 2.3e+03, 263 |  |
| SUG-08 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 8.76, 0 |  |
| SUG-09 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 23.7, 36.5 |  |
| SUG-10 | sugar_identity | retired | 2,5-Dimethylpyrazine | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-11 | sugar_identity | independent | FFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-12 | sugar_identity | independent | HMF | A>B | **disagree** | trunk | 948, 1.54e+03 |  |
| SUG-13 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 20.7, 4.04 |  |
| PH-06 | ph | independent | 2,5-Dimethylpyrazine | increasing | **agree** | trunk | 0.000311, 6.92, 21.8 |  |
| PH-07 | ph | independent | Furfural | flat | **disagree** | sulfur | 22.2, 20.7, 20.6 |  |
| TEMP-04 | temperature | independent | 2,5-Dimethylpyrazine | A>B | **agree** | trunk | 0.0016, 4.37e-06 |  |
| TEMP-05 | temperature | independent | HMF | increasing | **disagree** | trunk | 1.1e+03, 726, 651 |  |
| TEMP-06 | temperature | independent | Furfural | A>B | **agree** | sulfur | 20.7, 0.000551 |  |
| AW-03 | moisture_aw | independent | HMF | peak | **agree** | trunk | 1.3e+03, 1.47e+03, 760 |  |
| CYS-03 | additive_cysteine | independent | MFT | A>B | **agree** | sulfur | 34.8, 0 |  |
| CYS-04 | additive_cysteine | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'hydrolysate + xylose (no cysteine), 120 C' refused: PYRAZINE TARGETS '2,5-Dimethylpyrazine' (wave B18) run on the trunk lane only: the sulfur lane's networ |
| LIP-04 | lipid_lane | retired | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| ACR-01 | temperature | retired | Acrylamide | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| ACR-02 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 7.64, 935, 4.37 |  |
| TIME-02 | time | independent | HMF | increasing | **agree** | trunk | 53.2, 1.1e+03, 1.54e+03 |  |
| SCOPE-03 | scope | independent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOF-01 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HOF-02 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 2.12e+03, 2.83e+03 |  |
| HOF-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 0.246, 1.57 |  |
| HOF-04 | moisture_aw | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-01 | ph | independent | MFT | A>B | **agree** | sulfur | 154, 27.7 |  |
| MOT-02 | ph | independent | FFT | A>B | **agree** | sulfur | 3.34e+03, 85.9 |  |
| MOT-03 | ranking | independent | MFT | A>B | **disagree** | sulfur | 27.7, 27.7 | identical predictions |
| MOT-04 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-05 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-06 | ph | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-01 | ranking | fit_adjacent | MFT | A>B | **disagree** | sulfur | 3.54e+03, 9e+03 |  |
| HOX-02 | ranking | fit_adjacent | MFT | A>B | **agree** | sulfur | 9e+03, 1.81e+03 |  |
| HOX-03 | ph | independent | MFT | A>B | **disagree** | sulfur | 2.47e+03, 4.08e+03 |  |
| HOX-04 | ranking | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-05 | ph | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-06 | additive_cysteine | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-07 | ranking | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| YIL-01 | temperature | independent | MFT | decreasing | **disagree** | sulfur | 47.1, 57, 65.6, 71.3 |  |
| YIL-02 | temperature | independent | FFT | A>B | **agree** | sulfur | 947, 866 |  |
| WANG-01 | temperature | independent | MFT | peak | **not_evaluable** | sulfur | 31.2, 297, 972, 1.56e+03, 2.87e+03 | the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, disagree at ph 9); the source does not state it |
| WANG-02 | temperature | independent | FFT | peak | **agree** | sulfur | 47.9, 1.49e+03, 1.77e+04, 5.26e+04, 3.42e+04 |  |
| MENG-01 | temperature | independent | MFT, FFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-05 | moisture_aw | retired | Acrylamide | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-ACR-01 | ph | fit_adjacent | Acrylamide | A>B | **agree** | acrylamide | 5.31e+04, 2.73e+04 |  |
| DIC-01 | sugar_identity | independent | 3-deoxyglucosone, glucosone, methylglyoxal, glyoxal, diacetyl | ranking | **disagree** | trunk | 1.54e+04, 12, 978, 38.1, 0.0285 |  |
| DIC-02 | sugar_identity | independent | glucosone, glyoxal, methylglyoxal, diacetyl | ranking | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| DIC-03 | sugar_identity | independent | 3-deoxyglucosone, glucosone, glyoxal, methylglyoxal | ranking | **disagree** | trunk | 6.55e+04, 119, 161, 2.55e+03 |  |
| RIB-T-01 | time | independent | MFT | decreasing | **agree** | sulfur | 0.294, 0.00048, 4.55e-06 |  |
| RIB-T-02 | time | independent | FFT | flat | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HEX-T-01 | time | independent | FFT | increasing | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| SCH-T-01 | time | fit_system_overlap | MFT | increasing | **disagree** | sulfur | 17.3, 26, 0.593, 0.0742 |  |
| WANG22-T-01 | time | retired | MFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-02 | time | retired | FFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-03 | time | retired | MFT | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-04 | time | retired | FFT | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WHI-PH-01 | ph | fit_adjacent | MFT | A>B | **disagree** | sulfur | 4.38e+03, 6.36e+03 |  |
| CER07-PH-01 | ph | independent | FFT | decreasing | **disagree** | sulfur | 4.73e+03, 2.38, 0.899, 159 |  |
| CER07-PH-02 | ph | independent | furfural | decreasing | **agree** | sulfur | 2.64e+03, 1.34e+03, 570 |  |
| MOT02-PH-01 | ph | not_independent | MFT | A>B | **agree** | sulfur | 154, 27.7 |  |
| MOT02-PH-02 | ph | not_independent | FFT | A>B | **agree** | sulfur | 3.34e+03, 85.9 |  |
