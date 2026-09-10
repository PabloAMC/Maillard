# Wave B29 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (both air/argon ratios on one multiplier), T2, T3 (air exactly 1) and T4 hold.*

| test | result | pass |
|---|---|---|
| T1 the two air/argon ratios | model separates the pots by 0.97x against a printed 2.63x | False |
| T2 the copper arm | worst +0.16 dex | True |
| T3 air is exactly 1 | parameters identical True, observable identical True | True |
| T4 identification | sigma {"log10_f_argon": 0.25, "log10_f_air_cu": 1.62} | False |

| ratio | printed | model | dex |
|---|---:|---:|---:|
| arp_air_over_argon | 9.2 | 5.63 | -0.213 |
| glc_air_over_argon | 3.5 | 5.804 | +0.220 |
| arp_aircu_over_air | 2.5 | 1.729 | -0.160 |
| glc_aircu_over_air | 1.9 | 2.291 | +0.081 |

## What the axis exposed

| pot | share of the Strecker aldehyde made through the oxidative entries |
|---|---:|
| fed_amadori | 100.0% |
| glucose_glycine | 100.0% |

> The trunk's ONLY route to the Strecker aldehyde runs through glucosone, in BOTH pots. So the model's aldehyde is 100 % oxygen-dependent by construction, the two pots cannot have different air/argon ratios whatever the multiplier is, and under argon the model goes to ZERO. Hofmann measures 0.06 mol % from the Amadori compound and 0.04 from the sugar pot UNDER ARGON -- small, and not zero. A NON-OXIDATIVE route to the Strecker aldehyde exists and this model does not have one. That is the finding, and it is a structural gap the oxygen axis exposed rather than a bad multiplier.
