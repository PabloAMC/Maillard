# Kinetic core B6 -- HOLD-OUT report (the lipid-oxidation module)

Generated 2026-08-29 at `242bd1e9eae5`.

Scored against `results/validation/kinetic_core_b6_prereg.md`, written before the fit.

## Headline

| hold-out | role | outcome |
|---|---|---|
| Frankel alpha-tocopherol two-sided signature | seen_diagnostic | machinery checks PASS |
| Frankel nonanal ABSENCE | **gating** | PASS |
| exam refusals | -- | 17 -> 13 (pre-registered 13) |
| all eight registered row outcomes met? | -- | YES |
| G-1 regression guard (23 answered points unmoved) | -- | PASS |


### Pre-registration scorecard -- INCLUDING what was missed

| registered claim | outcome |
|---|---|
| E-1 fold: at least 3 of 4 hexanal points worse than 3x | **MET** (4/4) |
| E-1 direction: at least 3 of 4 OVER-predicted | **MISSED** (0 over, 4 under) |
| E-1 interval: at least 3 of 4 measured values inside the interval | **MISSED** (2/4) |


**The direction claim was MISSED, and the miss is the wave's single most important result.** The pre-registration expected OVER-prediction, on the strength of this repository's own standing diagnosis for this lane (the FAST lane's documented 36x and 1304x hexanal over-predictions). The B6 lipid lane UNDER-predicts all four points. The standing diagnosis does not transfer: it was a property of the FAST lane's unbounded first-order extrapolation from an invented 0.37 branch fraction, not of lipid chemistry.


**And the four points are not four of a kind.** Two of the four rows are 160 C process points; two are 40 C for 10 min, which measures an as-received isolate's ACCUMULATED storage oxidation, not a process output. A formation model with a declared no-formation-during-heating gap cannot make hexanal in 10 min at 40 C and should not be expected to. Reported as a benchmark-design finding, not corrected by fiat.


| class | rows | fold errors | inside interval |
|---|---:|---|---:|
| 160 C PROCESS points | 2 | 3.7x, 8.7x | 2/2 |
| 40 C / 10 min AMBIENT points | 2 | 3717x, 33392x | 0/2 |


On the two rows that are actually a thermal process, the lane lands at 3.7x and 8.7x with BOTH measured values inside the reported interval. On the two 40 C / 10 min rows it is 3 717x and 33 392x low, which is exactly what `parameters_lipid.LOOH_FORMATION_GAP` says will happen: the module cannot MAKE hydroperoxide during the hold, so on a near-ambient point it predicts that almost nothing forms -- and is then scored against the isolate's ACCUMULATED storage oxidation. Reported, not corrected by fiat.


bi_2020_raw_pea and liu_2023_ppi_offnote_baseline declare IDENTICAL conditions (40 C, 10 min, pH 6.0, a_w 0.95) and IDENTICAL precursors ('Pea Protein Isolate'), yet their measured hexanal differs 9.0x. No model that reads only the recorded fields can separate them, and this one predicts the same number for both. That is a defect in the BENCHMARK RECORDS, not in the model: whatever distinguishes the two samples is not written down.


## Hold-out 1 -- the alpha-tocopherol two-sided signature

**Role: seen_diagnostic (prereg sec. 0 -- the builder saw these columns).** Frankel 1989 prints the FIT column and the tocopherol columns in the same table rows, and states the result in its abstract. Under the Amendment 9 clause 1 precedent a seen hold-out cannot gate.

The donor parameter is NEVER FITTED and has no stored value. Every claim here holds for EVERY d in (0, 1), so there was nothing to tune toward the seen numbers even in principle.

| check | result |
|---|---|
| H1-a total volatile flux DECREASES | HOLDS |
| H1-b hexanal SHARE INCREASES | HOLDS |
| H1-d me-9-oxononanoate share INCREASES | HOLDS |
| H1-d pentane share DECREASES | HOLDS |
| H1-d methyl octanoate share DECREASES | HOLDS |
| H1-e me-13-oxo-tridecadienoate share DECREASES (registered as EXPECTED WRONG) | HOLDS |
| H1-e 2,4-decadienal share DECREASES (registered as EXPECTED WRONG) | HOLDS |
| H1-c the two move in OPPOSITE directions | HOLDS |

### The donor response curve

| d | total flux (rel.) | hexanal share | me-9-oxononanoate share | pentane share |
|---:|---:|---:|---:|---:|
| 0.00 | 1.0000 | 0.1055 | 0.1367 | 0.1201 |
| 0.05 | 0.9621 | 0.1096 | 0.1421 | 0.1186 |
| 0.10 | 0.9242 | 0.1141 | 0.1479 | 0.1170 |
| 0.20 | 0.8484 | 0.1243 | 0.1611 | 0.1133 |
| 0.35 | 0.7348 | 0.1435 | 0.1861 | 0.1063 |
| 0.50 | 0.6211 | 0.1698 | 0.2201 | 0.0967 |
| 0.65 | 0.5074 | 0.2078 | 0.2694 | 0.0829 |
| 0.80 | 0.3937 | 0.2679 | 0.3472 | 0.0610 |
| 0.90 | 0.3180 | 0.3317 | 0.4300 | 0.0378 |
| 0.95 | 0.2801 | 0.3766 | 0.4881 | 0.0214 |
| 0.99 | 0.2497 | 0.4223 | 0.5474 | 0.0048 |

The architecture treats methyl 13-oxo-tridecadienoate and 2,4-decadienal as homolytic CO-PRODUCTS, so their shares fall with the donor. Frankel's Discussion states that no trend was observed for those two, and his Fig. 4 excludes them from both numerator and denominator. This was registered as a PREDICTED FAILURE before the fit, so it cannot be presented as a discovery now. It is the sharpest open question the module leaves: either those two products are not pure homolytic co-products, or their measured flatness is a stability artefact of the kind Frankel himself warns about.

## Hold-out 2 -- the nonanal ABSENCE (gating)

Pure linoleate hydroperoxide feed, 180 C: nonanal = **0 mmol/L** (exact zero: True). **PASS.**

In a REAL matrix the oleate pool is not zero and the oleate -> nonanal branch fraction is measured NOWHERE, so the engine REFUSES absolute nonanal there rather than answering it. Both halves are the same hold-out: honouring an absence means refusing where the parent exists and answering zero where it does not.

In a real (pea) matrix the engine refuses nonanal: True.

## Hold-out 3 -- the exam

| bundle | compound | was | now | pre-registered | met? | measured | core | fold | old lane |
|---|---|---|---|---|---|---:|---:|---:|---:|
| `external_validation_bi_2020_raw_pea_hexanal` | hexanal | refused | answered | ANSWERED | YES | 1.26e+03 | 0.339 | 3.72e+03 | 1.12 |
| `external_validation_bi_2020_roasted_pea_hexanal` | hexanal | refused | answered | ANSWERED | YES | 324 | 88.6 | 3.66 | 2.75e+03 |
| `external_validation_li_2026_spi_wg_hme_control` | hexanal | refused | answered | ANSWERED | YES | 606 | 69.7 | 8.69 | 103 |
| `external_validation_li_2026_spi_wg_hme_control` | nonanal | refused | refused | STILL REFUSED | YES | 72.7 | -- | -- | 131 |
| `external_validation_li_2026_spi_wg_hme_control` | 2-pentylfuran | refused | refused | STILL REFUSED | YES | 5.63e+03 | -- | -- | 2.18 |
| `external_validation_li_2026_spi_wg_hme_control` | 1-hexanol | refused | refused | STILL REFUSED | YES | 20 | -- | -- | 1.24e+03 |
| `external_validation_liu_2023_ppi_offnote_baseline` | hexanal | refused | answered | ANSWERED | YES | 1.13e+04 | 0.339 | 3.34e+04 | 10.1 |
| `external_validation_liu_2023_ppi_offnote_baseline` | nonanal | refused | refused | STILL REFUSED | YES | 0.802 | -- | -- | 105 |

**E-1**: 4 hexanal points answered; 4 worse than 3x (pre-registered at least 3); 0 over-predicted and 4 UNDER-predicted (pre-registered at least 3 OVER-predicted -- **MISSED**); 2/4 measured values inside the reported interval (pre-registered at least 3 -- **MISSED**).


| bundle | measured | point | interval | inside? |
|---|---:|---:|---|---|
| `external_validation_bi_2020_raw_pea_hexanal` | 1.26e+03 | 0.339 | 0.0122 - 9.43 | False |
| `external_validation_bi_2020_roasted_pea_hexanal` | 324 | 88.6 | 4.05 - 1.94e+03 | True |
| `external_validation_li_2026_spi_wg_hme_control` | 606 | 69.7 | 2.13 - 2.28e+03 | True |
| `external_validation_liu_2023_ppi_offnote_baseline` | 1.13e+04 | 0.339 | 0.0122 - 9.43 | False |
**G-1**: 23 points were answered before B6; regressions: 0.

## Lane coupling

```json
{
 "may_cointegrate": true,
 "rule": "direct sum (independent parallel integration, same thermal program)",
 "reason": "Disjoint species sets, and the only candidate coupling (the aldehyde-lysine covalent channel) contributes exactly zero: it is INERT BY RULING, FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2. Revisit the moment the aldehyde-lysine Ea on food proteins is measured -- that is a NAMED WET-LAB GAP.",
 "shared_species": []
}
```

