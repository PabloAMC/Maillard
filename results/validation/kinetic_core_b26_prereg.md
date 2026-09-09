# Pre-registration: wave B26, the plant-protein binding row the matrix layer never had (written 2026-09-09, before anything ran)

## 1. Why

The matrix layer's reversible-binding table carries thirteen per-gram constants and every one of
them is an animal protein: skim milk, calcium caseinate, soy (via dialysis and gel filtration) and
beta-lactoglobulin. The class the panel leans on hardest, `n_alkanal`, is pooled from **exactly one
row**, Meynier's skim-milk hexanal constant. Meanwhile `data/species/protein_matrices.yml` charges a
full `pea_isolate` site table, the flagship hold-out is a plant paste, and the tool is meant for
people formulating with plant protein.

The reading audit found a paper that closes that gap on the very protein: Bi, Pan, Zhang and
colleagues, *Food Chem* 389:133044 (2022), "Non-covalent interactions of selected flavors with pea
protein". It measures a pea isolate against three volatiles at pH 7.6 and 37 °C, and among several
determinations it prints a **phase-ratio-variation partition pair** — a matrix leg and a water leg
in the same run on the same instrument — which is the exact construction the Meynier and
Leksrisompong rows were built from, and for the same stated reason: the absolute static-headspace
scale is suspect by 6 to 17× and the offset cancels in a within-run ratio.

## 2. What enters, and what does not

Three rows, computed from the printed pair by the registry's own formula
`K_g = (K_water/K_matrix − 1) / protein_g_per_L` at the paper's 10 g/L loading. Nothing here is
fitted; every value is arithmetic on printed cells, so seeing any downstream result cannot move it.

| key | compound | printed pair (matrix/gas over buffer/gas) | ratio | K_g, L/g | role |
|---|---|---|---:|---:|---|
| `kg_hexanal_pea` | hexanal | 116.37 / 32.90 | 3.537 | **2.537e-1** | FIT |
| `kg_z_2_penten_1_ol_pea` | (Z)-2-penten-1-ol | 2065.44 / 1460.49 | 1.414 | **4.14e-2** | FIT, new class, no consumer |
| `kg_t_2_octenal_pea` | (E)-2-octenal | 2203.85 / 455.93 | 4.834 | **3.834e-1** | QUARANTINED as a binding constant |

The octenal row is quarantined on the **Meynier precedent**, which is exact: a 2-alkenal held two
hours at 37 °C against a protein carrying about 0.47 mmol lysine amine per gram, measured by
disappearance, is partly irreversible Michael chemistry and not partition. The paper's only
reversibility evidence is a 20 % headspace recovery under guanidine, and a 20 % recovery is not a
demonstration that 80 % was reversible.

**Bi's alkenal/alkanal contrast is EXCLUDED from the unsaturation penalty, and that exclusion is
the wave's second finding.** Constructed the way the two carried observations are — a ratio of the
raw partition ratios — it is 4.834/3.537 = **1.367×**, which would pull the fitted penalty from
3.73× down to about 2.67× and into the 2–3× band the corpus states independently. It is excluded
anyway, because Meynier's and Vega's observations are **same-carbon** pairs (C6 alkenal against C6
alkanal) and Bi's is C8 against C6. Divide out the registry's own measured chain-length slope,
2.81×/CH₂, and two carbons alone would predict about 7.9× on the per-gram constant where Bi
measures 1.51×: the alkenal is **less** bound than chain length by itself would give. A contrast
that inverts once a measured confound is removed is not evidence for a penalty, and it is not
allowed to set one.

## 3. What this moves

One number, by arithmetic: the pooled headspace-family `n_alkanal` constant becomes the geometric
mean of the dairy and pea rows, and the branched-alkanal surrogate follows it down the 2.81×/CH₂
slope. Every other class is untouched. The predicted shift for hexanal in the Hong soy paste is
`1 + K_g × 142 g/L`, so it rises with the pooled constant.

## 4. What counts as success, declared before the run

- **T1 arithmetic.** The pooled `n_alkanal` constant equals the geometric mean of the two rows;
  `branched_alkanal` equals it divided by 2.81; the class's reference loading is unchanged at
  33.9 g/L; no other class constant moves by more than 1e-12.
- **T2 the flagship hold-out, decisive.** On the three Hong 2020 rows where the binding term is
  active (hexanal and the two methylbutanals) every fold error must **decrease**, no sign may
  invert, and the ten-row verdict must not get worse.
- **T3 the evidence ceiling, reported.** Amendment 6 ruling 2 caps the reversible term at ~25 % of
  an observed log-shift. Report the explained share on those three rows before and after.
- **T4 the loading, reported.** K_g is inversely proportional to the 10 g/L the paper states for its
  headspace assay and does not restate for the partition run, and the paper's own two headspace
  routes disagree by 2.5× on one compound. Report both as a band on the hexanal prediction.
- **T5 nothing else moves.** The sealed lupin and mucin keys stay sealed and valueless; the
  unsaturation penalty stays at exactly sqrt(2.81 × 4.95); the kinetic panel scorecard does not
  change at all.

Ship rule: **SHIP if T1, T2 and T5 hold.** T3 and T4 are reported.

## 5. Predictions, before the run

1. T2 holds and the hexanal fold error roughly halves, from about 50× to about 15×. **85 %.**
2. T3 **breaks the ceiling**: the reversible term crosses 25 % of the hexanal log-shift, at roughly
   40–45 %, for the first time on a real hold-out row. **75 %.** If it does, the wave must say which
   reading it takes, and it takes this one in advance: the cap was computed from **one compound in
   beef and one dairy protein**, and a plant isolate binding an alkanal 22× harder than cow's milk
   does is a reason to doubt the cap's transfer, not a reason to shrink a measured constant. The
   layer's flag stays and starts firing; that is the flag doing its job.
3. The residual stays enormous. Even at 15× fold the named terms explain under half of one row and
   nothing at all of six. **95 %.** This wave does not rescue the matrix layer; it removes one
   excuse, that the layer had never been given a plant protein to work with.

## 6. Outcome (2026-09-09, run the same day) — SHIP

**Verdict SHIP.** T1, T2 and T5 all held; T3 and T4 are reported. The full table is in
`kinetic_core_b26_ship_rule.md`.

| | before | after |
|---|---:|---:|
| pooled `n_alkanal` constant | 0.01151 L/g, one dairy row | **0.05404 L/g, two rows, one of them a plant protein** |
| hexanal in the Hong soy paste | 2.634× predicted against 132.5× measured | **8.673×** |
| that row's fold error | 50.3× | **15.3×** |
| the two methylbutanals' fold error | 166× and 165× | **70.6× and 70.1×** |
| reversible term's share of the hexanal log-shift | 19.8 % | **44.2 %** |

**All three pre-registered predictions held.** The fold error roughly halved as prediction 1 said
at 85 %; the evidence ceiling broke, at 44.2 % against a predicted 40–45 % at 75 %; and the residual
stayed enormous, as prediction 3 said at 95 % — six of the ten hold-out rows still have no term at
all and the layer still emits exactly 1.0 for them.

**What it cost, and where that is visible.** The n-alkanal class no longer reproduces either of its
own members: the Meynier dairy row it used to hit at 1.390× is now predicted at 2.832×, over by
2.04×. That is the price of pooling two proteins into one class constant, it is pinned in
`test_pinned_fit_row_reproduction_meynier_hexanal` rather than buried, and the layer takes it
because a class constant assembled from a single cow's-milk measurement was never a class constant.

**The ceiling.** Amendment 6 ruling 2's ~25 % cap was corroborated out of sample by wave B4 and is
now exceeded, on one row, by a term built from a measurement rather than from a fit. The reading
taken here was declared before the run: the cap came from one compound in beef and one dairy
protein, and a plant isolate that binds an alkanal 22× harder than cow's milk is a reason to doubt
that the cap transfers. The constant is not shrunk to fit the cap. The layer's flag stays and now
fires on that row, which is the flag doing its job, and the cap itself is not rewritten by this wave
— re-deriving it is its own question and it is in the backlog.

**A defect this wave found and fixed.** `kinetic_core_b4_frozen_predictions.json` records a blind
prediction, made before its wave read the paired thresholds. Its generator overwrote it
unconditionally on every run, date and all, so any later wave that changed the registry and re-ran
B4 would have silently replaced a pre-registration with a prediction made by somebody who had seen
the answer — and nothing afterwards could have told. The generator now refuses without an explicit
`--refreeze`. The blind record still carries the one-row pre-B26 constant, which is what makes it
the blind record.

**Not done here, and named.** The alkenal contrast is excluded and the unsaturation penalty is
untouched, so the layer's second term is still unvalidated out of sample and still sits above the
band the corpus states. And 37 °C is a mouth temperature: nothing in this wave licenses a pea
binding constant at 90 or 140 °C.
