# Wave B28 ship rule: SHIP

*Rule: SHIP if T1 (the arithmetic), T2 (refusals only ever lift, and nothing answered moves) and T5 (no other lane moves) hold; T3 and T4 reported. Pre-registration `results/validation/kinetic_core_b28_prereg.md`.*

| test | result | pass |
|---|---|---|
| T1 arithmetic | oleate columns sum to {'oleate_autoxidised': 100.0, 'oleate_photosensitized': 99.5}; 2-pentylfuran / hexanal {'linoleate_autoxidised': 0.16, 'linoleate_photosensitized': 0.0353}; the B6 six-product slate untouched True; nonanal still a structural zero from linoleate True | True |
| T2 the refusals | refused rows 25 -> 25; lifted 0; newly refused 0; answered rows that moved 0; lifted rows answered degenerately 0 | True |
| T3 1981 against 1989 | worst ME_9_OXONONANOATE 1.60x over five shared products | reported |
| T4 the new rows | 4 scored | reported |
| T5 nothing else moves | non-lipid benchmarks changed: none; rows ADDED by later waves (not a move): {'mp_holdout_glucose_only_autoclave_121C_Steinhagen2021': [('3,4-dideoxyglucosone', 'ppb'), ('3-deoxyglucosone', 'ppb'), ('glucosone', 'ppb'), ('glyoxal', 'ppb'), ('methylglyoxal', 'ppb')]} | True |

## The cross-laboratory check

Frankel 1981 against Frankel 1989, both renormalised onto the five products they both quantify. The C14 oxo-ester is dropped from both: 1981 could not identify it for want of an authentic reference, which is an analytical absence and not a chemical one.

| product | 1981 (%) | 1989 (%) | fold | higher in |
|---|---:|---:|---:|---|
| PENTANE | 13.6 | 20.0 | 1.47x | 1989_higher |
| HEXANAL | 20.6 | 13.8 | 1.50x | 1981_higher |
| ME_OCTANOATE | 20.6 | 21.2 | 1.03x | 1989_higher |
| DECADIENAL | 19.2 | 28.8 | 1.50x | 1989_higher |
| ME_9_OXONONANOATE | 26.1 | 16.2 | 1.60x | 1981_higher |

> Same laboratory, same first author, eight years apart. 210 C neat against 180 C in hexane, a 25 C column start against a -65 C cryotrap, packed against capillary. This is the first external check the lane's own fit source has had, and it is not a two-point Arrhenius: temperature and light-end loss are confounded and these two papers cannot separate them.

## What was lifted

- nothing

> A branch fraction from a NEAT hydroperoxide pyrolysed in an injector port at 210 C, predicting a food. Lifting a refusal is not the same as being right and this wave claims only the first.
