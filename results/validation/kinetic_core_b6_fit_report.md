# Kinetic core B6 -- FIT report (the lipid-oxidation module)

Generated 2026-08-29 at `242bd1e9eae5` (branch `audit-remediation`).

Pre-registration: `results/validation/kinetic_core_b6_prereg.md` -- written before this file existed.

## The one thing to read first

**The branch DISTRIBUTION is measured. The absolute RATE is not.** The rate enters as a bounded input anchored at Schroen & Berton-Carabin 2022's 25 C, hand-fitted, no-standard-error, lumped `k4 = 6e-3 /h`, with a Q10 that is the repo's own assumption carried with a [2, 3] band and a mandatory warning on every prediction. No Q10 number is baked into any stored constant.

## Hold-out exposure disclosure

THE BUILDER OF THIS WAVE SAW THE ALPHA-TOCOPHEROL HOLD-OUT COLUMNS. Frankel 1989 prints the zero-additive column (FIT) and the tocopherol columns (HOLD-OUT) in the SAME TABLE ROWS, and the abstract states the hold-out result in prose. Reading the fit column without seeing the hold-out column is impossible from that PDF. Consequence, adopted in the pre-registration sec. 0 and binding: the tocopherol hold-out is scored seen_diagnostic, never gating (the Amendment 9 clause 1 precedent), and the mitigation is structural -- the donor term has NO FITTED PARAMETER, so there is nothing to tune toward the seen numbers. The nonanal absence stays GATING because it is a structural zero fixed by molecular topology. The two frozen Bi 2020 bundles were never opened.

## The fit

| statistic | value | pre-registered bound |
|---|---:|---:|
| median abs residual (percentage points) | 0.55 | 3.0 |
| worst abs residual (percentage points) | 3.99 | 8.0 |
| fit values / free parameters / df | 18 / 12 / 3 | -- |
| **F-1 met?** | **YES** | -- |


### Measured vs predicted, all eighteen FIT numbers

| system | product | measured % | predicted % | residual pp |
|---|---|---:|---:|---:|
| mixed_ct_tt_9_13 | PENTANE | 16.0 | 12.0 | -3.99 |
| mixed_ct_tt_9_13 | HEXANAL | 11.0 | 10.5 | -0.45 |
| mixed_ct_tt_9_13 | ME_OCTANOATE | 17.0 | 17.6 | +0.64 |
| mixed_ct_tt_9_13 | DECADIENAL | 23.0 | 23.1 | +0.06 |
| mixed_ct_tt_9_13 | ME_9_OXONONANOATE | 13.0 | 13.7 | +0.67 |
| mixed_ct_tt_9_13 | ME_13_OXO_TRIDECADIENOATE | 20.0 | 23.1 | +3.07 |
| pure_ct_13 | PENTANE | 21.0 | 23.2 | +2.18 |
| pure_ct_13 | HEXANAL | 20.0 | 20.3 | +0.35 |
| pure_ct_13 | ME_OCTANOATE | 5.0 | 3.8 | -1.19 |
| pure_ct_13 | DECADIENAL | 3.5 | 5.0 | +1.48 |
| pure_ct_13 | ME_9_OXONONANOATE | 4.3 | 3.0 | -1.34 |
| pure_ct_13 | ME_13_OXO_TRIDECADIENOATE | 46.0 | 44.5 | -1.48 |
| tt_9_13 | PENTANE | 1.5 | 1.5 | +0.00 |
| tt_9_13 | HEXANAL | 13.0 | 13.0 | -0.00 |
| tt_9_13 | ME_OCTANOATE | 13.0 | 13.0 | -0.00 |
| tt_9_13 | DECADIENAL | 30.0 | 30.0 | -0.00 |
| tt_9_13 | ME_9_OXONONANOATE | 26.0 | 26.0 | -0.00 |
| tt_9_13 | ME_13_OXO_TRIDECADIENOATE | 16.0 | 16.0 | -0.00 |

### Frozen branch simplexes

| hydroperoxide | product | share |
|---|---|---:|
| 13_ct | PENTANE | 0.2633 |
| 13_ct | HEXANAL | 0.2311 |
| 13_ct | ME_13_OXO_TRIDECADIENOATE | 0.5056 |
| 13_tt | PENTANE | 0.0492 |
| 13_tt | HEXANAL | 0.4262 |
| 13_tt | ME_13_OXO_TRIDECADIENOATE | 0.5246 |
| 9_ct | ME_OCTANOATE | 0.3244 |
| 9_ct | DECADIENAL | 0.4241 |
| 9_ct | ME_9_OXONONANOATE | 0.2515 |
| 9_tt | ME_OCTANOATE | 0.1884 |
| 9_tt | DECADIENAL | 0.4348 |
| 9_tt | ME_9_OXONONANOATE | 0.3768 |

### Fitted hydroperoxide-pool compositions

| system | 13-ct | 13-tt | 9-ct | 9-tt |
|---|---:|---:|---:|---:|
| mixed_ct_tt_9_13 | 0.456 | 0.000 | 0.544 | 0.000 |
| pure_ct_13 | 0.882 | 0.000 | 0.118 | 0.000 |
| tt_9_13 | 0.000 | 0.307 | 0.000 | 0.693 |

## F-2 -- the shipped `hexanal 0.37`

Fitted hexanal share of the six-product slate at the autoxidation reference composition: **0.1055**. Frankel's own zero-additive range: 11-20 %. Shipped FAST-lane value: 0.37. **REFUTED.**

Compared LIKE FOR LIKE: both are a hexanal share of the six-product slate. The per-isomer simplex values are NOT comparable to 0.37 -- they are shares of a three-product sub-slate -- and quoting one of those instead would be the kind of denominator swap this repository has already been caught by.

## F-3 -- the shipped `nonanal 0.15`

STRUCTURAL ZERO from every linoleate pool. Its only parent is LOOH_OL, whose branch fraction is unmeasured, so it is answered as exactly 0.0 for a linoleate feed and REFUSED for an oleate-bearing matrix.

## The falsified 1:1 pairing (pre-registered as F-1's reason)

| system | pentane : methyl 13-oxo-tridecadienoate |
|---|---:|
| mixed_ct_tt_9_13 | 0.800 |
| pure_ct_13 | 0.457 |
| tt_9_13 | 0.094 |

Spread: **8.5x** against an expected 1.0. Pentane and methyl 13-oxo-9,11-tridecadienoate are the two halves of ONE beta-scission and should be 1:1. Their measured ratio across the three zero-additive FIT columns spans 8.5x (0.80 / 0.46 / 0.094). Either the GC recovery of a C5 alkane and a C14 dienoate differ enormously, or the pairing is wrong. This module does NOT impose the pairing and does NOT absorb the discrepancy into a response factor; it fits free shares inside each isomer's simplex and reports the falsified pairing.

## What was refused

* **Ea for LOOH decomposition from Frankel 1989** -- ONE TEMPERATURE (180 C). k3 sec. C.9, verbatim: 'NOT a yield source ... no absolute yield and no Ea exist in it'. Re-affirmed by research_round3_channels.md sec. F.3 after a further search round.
* **Q10 baked into a rate constant** -- The Q10 is the REPO'S assumption, not the authors'. It is exposed as Q10_ASSUMPTION with a band and a mandatory warning, and no stored constant is ever pre-multiplied by it.
* **oleate -> nonanal branch fraction** -- Measured nowhere in the fit corpus. Frankel 1989 fed linoleate only. The FAST lane's shipped 0.15 has no source. Requests for absolute nonanal in an oleate-bearing matrix are REFUSED.
* **linoleate -> 2-pentylfuran branch fraction** -- Not in Frankel's six-product slate and measured nowhere else in the corpus. The FAST lane's shipped 0.08 has no source.
* **aldehyde -> alcohol reduction (hexanal -> 1-hexanol)** -- No reduction step is measured anywhere in the corpus, and in a thermally processed extrudate the reductant pool is not even identified. Refused.
* **propanal branch from the linolenate channel** -- Schroen's 7 % propanal is a property of RAPESEED OIL's 92.2 mg/g alpha-linolenate, not a branch fraction of a linoleate hydroperoxide. Frankel's slate contains no propanal. Importing it would import an oil's composition as if it were chemistry.
* **a mass yield of hexanal per hydroperoxide** -- Frankel's shares are fractions of SIX MEASURED PEAKS, and two of the scissions' partners (2-nonenal, methyl 12-oxo-10-dodecenoate) are named in his introduction and quantified in none of his tables. The share is therefore a WITHIN-SLATE share, and this module carries the difference in LIPID_FRAG_C rather than pretending the slate closes.

## Contradictions reported, not resolved

* **schroen_1.2pct_vs_frankel_11_20pct** -- Schroen's hexanal is 1.2 % of the secondary-product flux at 25 C in whole rapeseed oil; Frankel's is 11-20 % of a six-peak slate at 180 C from a PURE linoleate hydroperoxide. These differ by ~10-17x and they are NOT the same quantity: different denominator (all secondary products vs six peaks), different feed (whole oil, 62 % oleate, vs pure linoleate hydroperoxide), different temperature (25 vs 180 C) and different mechanism mix (Frankel's heterolytic Hock route is heat-promoted and cannot be active at 25 C). Reported side by side; never averaged, never reconciled by a fitted factor.
* **pentane_vs_me13oxo_pairing** -- Pentane and methyl 13-oxo-9,11-tridecadienoate are the two halves of ONE beta-scission and should be 1:1. Their measured ratio across the three zero-additive FIT columns spans 8.5x (0.80 / 0.46 / 0.094). Either the GC recovery of a C5 alkane and a C14 dienoate differ enormously, or the pairing is wrong. This module does NOT impose the pairing and does NOT absorb the discrepancy into a response factor; it fits free shares inside each isomer's simplex and reports the falsified pairing.
