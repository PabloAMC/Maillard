# Wave B29 ship rule: DO NOT SHIP

*Rule: SHIP if T1 (both air/argon ratios on one multiplier), T2, T3 (air exactly 1) and T4 hold.*

| test | result | pass |
|---|---|---|
| T1 the two air/argon ratios | model separates the pots by 0.79x against a printed 2.63x | False |
| T2 the copper arm | worst +0.26 dex | False |
| T3 air is exactly 1 | parameters identical True, observable identical True | True |
| T4 identification | sigma {"log10_f_argon": 42.09, "log10_f_air_cu": 18.97} | False |

| ratio | printed | model | dex |
|---|---:|---:|---:|
| arp_air_over_argon | 9.2 | 2.188 | -0.624 |
| glc_air_over_argon | 3.5 | 2.756 | -0.104 |
| arp_aircu_over_air | 2.5 | 1.387 | -0.256 |
| glc_aircu_over_air | 1.9 | 2.142 | +0.052 |

## What the axis exposed

| pot | share of the Strecker aldehyde made through the oxidative entries |
|---|---:|
| fed_amadori | 54.6% |
| glucose_glycine | 64.1% |

> CORRECTED ON REVIEW, 2026-09-10. The first run measured AKG alone -- glyoxal's Strecker product -- found it 100 % oxidative in both pots, and concluded that a non-oxidative route to the Strecker aldehyde was missing. That was an artefact of the observable: the trunk HAS a non-oxidative route (Amadori -> 1-deoxyosone -> methylglyoxal -> AKM), and summed over both Strecker products the non-oxidative share is 45 % in the Amadori pot and 36 % in the sugar pot. THE REAL FINDING IS THE ORDER. Hofmann's Amadori pot is the MORE oxygen-sensitive (9.2x against 3.5x), so its oxidative share must be the larger. The model's is the SMALLER. That is a statement about the trunk's branching between the oxidative route to glucosone and the non-oxidative routes to the deoxyosones, and it is the reason one multiplier cannot serve both pots.
