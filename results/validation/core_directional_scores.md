# Core directional scores (the kinetic core on the directional claims panel)

Panel `docs/validation/directional_claims_panel.yml` (119 claims, flat tolerance 5 %). Nothing in the core was tuned to this panel.

* **headline (strictly independent, evaluable): 21/38 (55%)**; 37 independent claims not evaluable
* independent, excluding pH and water activity: 15/27 (56%); pH and water activity alone: 6/11 (55%)
* all claims (independent + fit-adjacent + fit-system overlap): 31/61 (51%); 58 not evaluable
* misses where the lane carries no term for the moved axis (identical predictions): 3
* not evaluable, by reason: the claim carries no runnable conditions (prose-only) (39); a predicted concentration is zero; no direction is defined (7); refused by the engine (3); arm 'hydrolysate + xylose (no cysteine), 120 C' refused (1); the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, a (1); arm 'lysine 1x (40 mM)' refused (1); a predicted concentration is zero; no ranking is defined (1); arm 'glucose + glutathione' refused (1); arm 'cysteine alone' refused (1); arm 'beef triglyceride alone' refused (1); arm 'egg PC alone' refused (1); arm 'methionine alone' refused (1)

## Per category (strictly independent claims)

| category | agree | evaluable | rate | not evaluable | misses |
|---|---|---|---|---|---|
| additive_amine | 0 | 0 | - | 2 | - |
| additive_cysteine | 1 | 1 | 1.00 | 4 | - |
| lipid_lane | 0 | 0 | - | 3 | - |
| matrix_identity | 0 | 0 | - | 1 | - |
| moisture_aw | 1 | 1 | 1.00 | 3 | - |
| ph | 5 | 10 | 0.50 | 5 | PH-07, HOX-03, CER07-PH-01, CER07-PH-02, MOD08-PH-01 |
| process_heating | 0 | 0 | - | 1 | - |
| ranking | 0 | 1 | 0.00 | 3 | MOT-03 |
| scope | 0 | 0 | - | 3 | - |
| sugar_identity | 5 | 11 | 0.45 | 5 | SUG-03, SUG-12, HOF-02, HOF-03, DIC-01, DIC-03 |
| sulfur_donor | 2 | 2 | 1.00 | 3 | - |
| temperature | 5 | 9 | 0.56 | 2 | TEMP-01, TEMP-05, YIL-01, YIL-02 |
| time | 2 | 3 | 0.67 | 2 | Z23C-T-01 |

## Per category (all claims)

| category | agree | evaluable | rate | not evaluable |
|---|---|---|---|---|
| additive_amine | 0 | 0 | - | 2 |
| additive_cysteine | 2 | 3 | 0.67 | 5 |
| lipid_lane | 0 | 1 | 0.00 | 5 |
| matrix_identity | 0 | 2 | 0.00 | 1 |
| moisture_aw | 1 | 1 | 1.00 | 6 |
| ph | 10 | 17 | 0.59 | 8 |
| process_heating | 0 | 0 | - | 4 |
| ranking | 1 | 4 | 0.25 | 3 |
| scope | 0 | 0 | - | 3 |
| sugar_identity | 7 | 15 | 0.47 | 7 |
| sulfur_donor | 2 | 3 | 0.67 | 3 |
| temperature | 5 | 9 | 0.56 | 4 |
| time | 3 | 6 | 0.50 | 7 |

## Claims

| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |
|---|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | fit_adjacent | MFT | A>B | **agree** | sulfur | 93.6, 0 |  |
| SUG-02 | sugar_identity | fit_adjacent | FFT | A>B | **disagree** | sulfur | 347, 746 |  |
| SUG-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 55.8, 55.8 | identical predictions |
| SUG-04 | sugar_identity | fit_adjacent | MFT | A>B | **disagree** | sulfur | 93.6, 93.6 | identical predictions |
| SUG-05 | sugar_identity | fit_system_overlap | Furfural | A>B | **agree** | sulfur | 1.38e+03, 625 |  |
| SUG-06 | sugar_identity | retired | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-01 | ph | fit_system_overlap | MFT | A>B | **agree** | sulfur | 100, 7.21 |  |
| PH-02 | ph | fit_system_overlap | bis(2-methyl-3-furyl) disulfide | A>B | **agree** | sulfur | 3.9, 1.85 |  |
| PH-03 | ph | not_independent | FFT | decreasing | **disagree** | sulfur | 3.99e+03, 12.1, 113 |  |
| PH-04 | ph | retired | 2,5-Dimethylpyrazine | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-05 | ph | retired | Furfural | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| TEMP-01 | temperature | independent | Acrylamide | decreasing | **disagree** | acrylamide | 13.6, 357, 1.04e+03 |  |
| TEMP-02 | temperature | independent | Acrylamide | A>B | **agree** | acrylamide | 13.6, 0 |  |
| TEMP-03 | temperature | retired | HMF | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-01 | moisture_aw | retired | HMF | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-02 | moisture_aw | retired | Acrylamide | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| CYS-01 | additive_cysteine | not_independent | FFT | A>B | **agree** | sulfur | 52.3, 0 |  |
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
| SUG-07 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 2.4e+03, 420 |  |
| SUG-08 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 17.7, 0 |  |
| SUG-09 | sugar_identity | independent | FFT | A>B | **agree** | sulfur | 60, 32 |  |
| SUG-10 | sugar_identity | retired | 2,5-Dimethylpyrazine | decreasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-11 | sugar_identity | independent | FFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| SUG-12 | sugar_identity | independent | HMF | A>B | **disagree** | trunk | 948, 1.54e+03 |  |
| SUG-13 | sugar_identity | independent | Furfural | A>B | **agree** | sulfur | 17.5, 3.42 |  |
| PH-06 | ph | independent | 2,5-Dimethylpyrazine | increasing | **agree** | trunk | 0.000311, 6.92, 21.8 |  |
| PH-07 | ph | independent | Furfural | flat | **disagree** | sulfur | 18.7, 17.5, 17.4 |  |
| TEMP-04 | temperature | independent | 2,5-Dimethylpyrazine | A>B | **agree** | trunk | 0.0016, 4.37e-06 |  |
| TEMP-05 | temperature | independent | HMF | increasing | **disagree** | trunk | 1.1e+03, 726, 651 |  |
| TEMP-06 | temperature | independent | Furfural | A>B | **agree** | sulfur | 17.5, 0.000551 |  |
| AW-03 | moisture_aw | independent | HMF | peak | **agree** | trunk | 1.3e+03, 1.47e+03, 760 |  |
| CYS-03 | additive_cysteine | independent | MFT | A>B | **agree** | sulfur | 104, 0 |  |
| CYS-04 | additive_cysteine | independent | 2,5-Dimethylpyrazine | A>B | **not_evaluable** | - | - | arm 'hydrolysate + xylose (no cysteine), 120 C' refused: PYRAZINE TARGETS '2,5-Dimethylpyrazine' (wave B18) run on the trunk lane only: the sulfur lane's networ |
| LIP-04 | lipid_lane | retired | Hexanal | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| ACR-01 | temperature | retired | Acrylamide | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| ACR-02 | temperature | independent | Acrylamide | peak | **agree** | acrylamide | 7.64, 935, 4.37 |  |
| TIME-02 | time | independent | HMF | increasing | **agree** | trunk | 53.2, 1.1e+03, 1.54e+03 |  |
| SCOPE-03 | scope | independent | 2-Pentylfuran | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOF-01 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HOF-02 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 543, 746 |  |
| HOF-03 | sugar_identity | independent | FFT | A>B | **disagree** | sulfur | 0.74, 1.2 |  |
| HOF-04 | moisture_aw | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-01 | ph | independent | MFT | A>B | **agree** | sulfur | 97, 45.5 |  |
| MOT-02 | ph | independent | FFT | A>B | **agree** | sulfur | 2.12e+03, 215 |  |
| MOT-03 | ranking | independent | MFT | A>B | **disagree** | sulfur | 45.5, 45.5 | identical predictions |
| MOT-04 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-05 | sugar_identity | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MOT-06 | ph | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-01 | ranking | fit_adjacent | MFT | A>B | **disagree** | sulfur | 3.54e+03, 8.46e+03 |  |
| HOX-02 | ranking | fit_adjacent | MFT | A>B | **agree** | sulfur | 8.46e+03, 1.28e+03 |  |
| HOX-03 | ph | independent | MFT | A>B | **disagree** | sulfur | 2.46e+03, 4.08e+03 |  |
| HOX-04 | ranking | independent | MFT | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-05 | ph | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-06 | additive_cysteine | independent | MFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HOX-07 | ranking | independent | FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| YIL-01 | temperature | independent | MFT | decreasing | **disagree** | sulfur | 10.3, 12, 12.9, 12.7 |  |
| YIL-02 | temperature | independent | FFT | A>B | **disagree** | sulfur | 167, 299 |  |
| WANG-01 | temperature | independent | MFT | peak | **not_evaluable** | sulfur | 7.11, 62.7, 198, 404, 675 | the verdict depends on the unstated input ph (disagree at ph 5, agree at ph 7, agree at ph 9); the source does not state it |
| WANG-02 | temperature | independent | FFT | peak | **agree** | sulfur | 49.5, 1.82e+03, 2.06e+04, 4.56e+04, 2.76e+04 |  |
| MENG-01 | temperature | independent | MFT, FFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| AW-05 | moisture_aw | retired | Acrylamide | flat | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| PH-ACR-01 | ph | fit_adjacent | Acrylamide | A>B | **agree** | acrylamide | 5.31e+04, 2.73e+04 |  |
| DIC-01 | sugar_identity | independent | 3-deoxyglucosone, glucosone, methylglyoxal, glyoxal, diacetyl | ranking | **disagree** | trunk | 1.54e+04, 12, 978, 38.1, 0.0285 |  |
| DIC-02 | sugar_identity | independent | glucosone, glyoxal, methylglyoxal, diacetyl | ranking | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| DIC-03 | sugar_identity | independent | 3-deoxyglucosone, glucosone, glyoxal, methylglyoxal | ranking | **disagree** | trunk | 6.55e+04, 119, 161, 2.55e+03 |  |
| RIB-T-01 | time | independent | MFT | decreasing | **agree** | sulfur | 2.06, 0.031, 0.00132 |  |
| RIB-T-02 | time | independent | FFT | flat | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| HEX-T-01 | time | independent | FFT | increasing | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| SCH-T-01 | time | fit_system_overlap | MFT | increasing | **disagree** | sulfur | 12.6, 30.8, 35.2, 19.1 |  |
| WANG22-T-01 | time | retired | MFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-02 | time | retired | FFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-03 | time | retired | MFT | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WANG22-T-04 | time | retired | FFT | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| WHI-PH-01 | ph | fit_adjacent | MFT | A>B | **disagree** | sulfur | 3.5e+03, 3.64e+03 |  |
| CER07-PH-01 | ph | independent | FFT | decreasing | **disagree** | sulfur | 3.99e+03, 1.18e+03, 1.1, 113 |  |
| CER07-PH-02 | ph | independent | furfural | decreasing | **disagree** | sulfur | 1.88e+03, 2.07e+03, 340 |  |
| MOT02-PH-01 | ph | not_independent | MFT | A>B | **agree** | sulfur | 97, 45.5 |  |
| MOT02-PH-02 | ph | not_independent | FFT | A>B | **agree** | sulfur | 2.12e+03, 215 |  |
| ZHA25-LYS-01 | additive_amine | independent | MFT | peak | **not_evaluable** | - | - | arm 'lysine 1x (40 mM)' refused: LANE CONFLICT: this request needs both the acrylamide and sulfur lanes at once. They do not compose -- the acrylamide network d |
| ZHOU22-T-01 | time | fit_system_overlap | MFT | increasing | **disagree** | sulfur | 75.8, 170, 192, 178 |  |
| ZHOU22-T-02 | time | fit_system_overlap | FFT | increasing | **agree** | sulfur | 38.2, 131, 215, 260 |  |
| ZHOU22-CYS-01 | additive_cysteine | fit_system_overlap | furfural | A>B | **disagree** | sulfur | 534, 798 |  |
| Z23C-GSH-01 | sulfur_donor | fit_system_overlap | MFT | A>B | **disagree** | sulfur | 134, 1.75e+03 |  |
| Z23C-T-01 | time | independent | FFT | peak | **disagree** | sulfur | 25.7, 73.5, 111, 130 |  |
| CHEN02-RANK-01 | ranking | independent | FFT, MFT | ranking | **not_evaluable** | - | - | a predicted concentration is zero; no ranking is defined |
| CHEN02-PH-01 | ph | independent | FFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| CHEN02-PH-02 | ph | independent | MFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| CHEN00-PH-01 | ph | not_independent | FFT | A>B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| XU07-PH-01 | ph | independent | MFT | decreasing | **agree** | sulfur | 45.5, 0.204, 0.00326 |  |
| XU07-PH-02 | ph | independent | FFT | A>B | **agree** | sulfur | 72.5, 0.0404 |  |
| WANG12-01 | sulfur_donor | independent | MFT | A>B | **agree** | sulfur | 72.4, 42.3 |  |
| WANG12-02 | sulfur_donor | independent | FFT | A>B | **agree** | sulfur | 170, 131 |  |
| HIN92-01 | sugar_identity | independent | MFT | A>B | **agree** | sulfur | 1.21e+03, 682 |  |
| ZHG91-01 | sulfur_donor | independent | 2-acetylthiazole | A>B | **not_evaluable** | - | - | arm 'glucose + glutathione' refused: THIS POT CHARGES NO PRECURSOR THAT COULD MAKE '2-acetylthiazole': in the sulfur lane's network no chain of reactions leads  |
| ZHG91-02 | sulfur_donor | independent | furfural | A<B | **not_evaluable** | - | - | a predicted concentration is zero; no direction is defined |
| ZHG88-01 | sulfur_donor | independent | 2-acetylthiazole | A>B | **not_evaluable** | - | - | arm 'cysteine alone' refused: THIS POT CHARGES NO PRECURSOR THAT COULD MAKE '2-acetylthiazole': in the sulfur lane's network no chain of reactions leads from wh |
| FAR92-01 | lipid_lane | independent | hexanal | A>B | **not_evaluable** | - | - | arm 'beef triglyceride alone' refused: UNMAPPED PRECURSORS 'Beef triglyceride': not a species in any core lane. The core is a named small-molecule network; an i |
| FAR92-02 | lipid_lane | independent | 2,4-decadienal | A>B | **not_evaluable** | - | - | arm 'egg PC alone' refused: UNMAPPED PRECURSORS 'Egg phosphatidylcholine': not a species in any core lane. The core is a named small-molecule network; an intact |
| ZHG23B-01 | additive_cysteine | independent | dimethyl disulfide | A>B | **not_evaluable** | - | - | arm 'methionine alone' refused: METHIONINE CHAIN TARGETS 'dimethyl disulfide': the methionine chain (wave B22) did not ship: fitted on Pan 2025's rates, the ide |
| MOD08-PH-01 | ph | independent | MFT | A>B | **disagree** | sulfur | 209, 345 |  |
| MOD08-AW-01 | moisture_aw | independent | MFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| HAR84B-AW-01 | moisture_aw | independent | MFT | peak | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAD95-01 | ph | independent | MFT | increasing | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAR10-CYS-01 | additive_cysteine | independent | HMF, furfural | A<B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
| MAR10-GLY-01 | additive_amine | independent | MFT, FFT | A>B | **not_evaluable** | - | - | the claim carries no runnable conditions (prose-only) |
