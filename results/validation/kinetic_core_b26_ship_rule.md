# Wave B26 ship rule: SHIP

*Rule: SHIP if T1 (the arithmetic), T2 (every active hold-out row improves and no sign inverts) and T5 (nothing else moves) hold; T3 and T4 reported. Pre-registration `results/validation/kinetic_core_b26_prereg.md`.*

Rows installed: kg_hexanal_pea, kg_z_2_penten_1_ol_pea (pooled) and kg_t_2_octenal_pea (quarantined as a binding constant and excluded from the unsaturation fit).

| test | result | pass |
|---|---|---|
| T1 arithmetic | n_alkanal 0.054038 L/g from hexanal@skim_milk, hexanal@pea_protein_1pct; branched 0.019231; reference loading 33.9 -> 33.9 g/L; classes moved ['alkenol', 'branched_alkanal', 'n_alkanal']; unmoved ['diketone', 'ester', 'furanone', 'lactone', 'methyl_ketone'] | True |
| T2 flagship hold-out | rows made worse: none; within 5x 0 -> 0 of 10; signs correct 3 -> 3 | True |
| T3 evidence ceiling | rows now over the 25% cap: ['hexanal'] | reported |
| T4 the loading | hexanal in the Hong paste from 6.43x to 11.85x across a factor of four in the assumed loading; live 8.67x | reported |
| T5 nothing else moves | 6 sealed keys still valueless; penalty 3.7295x on 2 rows; kinetic panel identical True | True |

## The three rows that carry a binding term

| compound | measured | predicted before | predicted after | fold before | fold after | explained share before | after |
|---|---:|---:|---:|---:|---:|---:|---:|
| hexanal | 132.5x | 2.634x | 8.673x | 50.3x | 15.28x | 19.8% | 44.2% |
| 3_methylbutanal | 263.2x | 1.582x | 3.731x | 166.4x | 70.55x | 8.2% | 23.6% |
| 2_methylbutanal | 261.4x | 1.582x | 3.731x | 165.3x | 70.07x | 8.2% | 23.7% |

## The seven rows that do not

Unchanged, and unchanged for the same reason as before: the corpus supplies no binding constant for their classes, so the layer emits exactly 1.0 and reports the whole shift as unexplained residual. This wave adds a protein, not a class.

## What the ceiling says now

> The cap was computed from ONE compound in beef and one dairy protein. A plant isolate that binds an alkanal 22x harder than cow's milk does is a reason to doubt that the cap transfers, not a reason to shrink a measured constant. The layer's flag stays and starts firing; that is the flag doing its job.
