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


---

## Two corrections from B32's premise check (2026-09-10)

**1. "Share of the Strecker aldehyde made through the oxidative entries" is not a share.** The
quantity is `(air − argon) / air`, and under argon the NON-oxidative product goes UP — AKM by 81 % in
the fed pot and 93 % in the sugar pot — because the Amadori compound no longer drained to glucosone
flows down the deoxyosone route instead. It is a net difference between two competing channels. The
share that is oxidative by construction is `AKG / (AKG + AKM)`: **74.9 %** (fed) and **81.3 %**
(sugar). The order is wrong on that measure too, so nothing in this ship rule's verdict moves.

**2. The table above is a difference in EXTENT, not in branching.** After 120 min at 100 °C the fed
pot has 0.124 mM of Amadori compound left and no glucose; the sugar pot has 2.13 mM and 73.9 mM of
glucose still unreacted. Run the sugar pot to matched extent and its share falls to **55.4 %**
against the fed pot's 54.6 %. The two pots branch identically. See
`kinetic_core_b32_premise_check.md`; the refit this ship rule's closing sentence pointed at was
proposed, probed and refused without a constant being moved.
