# Pre-registration: wave B31, what the isolate arrived with (written 2026-09-10, before anything ran)

## 1. Why

Four panel rows ask the model to reproduce a level it did not make. `pea_isolate_40C_PratapSingh2021`
and `soy_isolate_40C_PratapSingh2021` are a **ten-minute hold at 40 °C**. Essentially nothing forms
in ten minutes at 40 °C, so what was measured is what the isolate arrived with, and the model
predicts formation from zero. It misses by **3357× and 6078× on hexanal** in those two pots and has
done since long before this week; nothing had named the cause.

The same thing shows from the other side in the pea beverage. Trikusuma prints three columns and the
panel scores only the heated one: hexanal **331 → 782 µg/L**, 2-pentylfuran **59.4 → 163**, nonanal
**8.24 → 24.0**. Between 42 % and 36 % of every measured level was in the beverage before any heat.
Scoring the total against a formation model charges the model for the raw material.

## 2. What is built

A declared input, not a fitted one. `ProcessSpec.carried_volatiles`: a mapping of compound name to
the concentration in µg/L the pot **starts** with, read from a benchmark's own
`conditions.carried_volatiles`. The engine adds it to the integrated concentration **before** the
matrix-binding factor, because the protein does not know which molecule was carried and which was
made, and both are equally available to bind.

**It is exactly zero where nothing is declared,** so every other row is bit-for-bit what it was.

## 3. What it may and may not be used for

**MAY.** A carried level the source PRINTS, in the same table, same units, same pot, as an unheated
control. Trikusuma prints one. It enters as a declared measurement with the source's own anchor.

**MAY NOT.** A carried level inferred, assumed, or taken from another paper's isolate. The two
PratapSingh pots have **no dossier on this disk and no control column**, so nothing can be declared
for them and nothing will be invented. What happens to those rows is section 4's second test.

## 4. What counts as success, declared before the run

- **T1 zero by default.** Every row that declares no carried volatile is bit-for-bit unchanged:
  identical predicted value, to the last bit, for all 46 scored rows but the three Trikusuma ones.
- **T2 the three declared rows improve.** All three fold errors must fall. A carried level that makes
  a prediction worse would mean the model was right for the wrong reason and this wave would be
  hiding that.
- **T3 the two undeclarable pots.** Report what a rule refusing carried-volatile rows without a
  declared starting state would refuse across the whole panel, and adopt it only if it refuses
  exactly the rows that cannot be declared. A rule that catches answerable rows is the wrong rule.
- **T4 nothing else moves.** The sulfur, trunk and acrylamide lanes bit-for-bit.

Ship rule: **SHIP if T1, T2 and T4 hold.** T3 is decided on its own evidence and reported.

## 5. Predictions, before the run

1. T2 holds on all three. **85 %.** The declared levels are 36–42 % of the measured totals and the
   model under-predicts all three, so subtracting a carried floor can only close the gap.
2. **Nonanal enters the 3× band.** **60 %.** It sits at 3.3× now and its carried share is 34 %;
   7.37 predicted against a formed 15.8 is 2.1×. That would be the first new row inside the band in
   this whole sequence of waves, and the first accuracy gain rather than an honesty gain.
3. T3's rule refuses exactly the two PratapSingh pots and nothing else. **40 %.** I expect it to
   catch the two external isolate rows as well, and if it does the rule is too broad and is not
   adopted.

---

# Outcome (run 2026-09-10, after the pre-registration above)

## The scoreboard

```
before   panel 37 | scored 27 | rows 46 (refused 18) | within 3x  4/46 | out-of-sample  3/45
declared panel 37 | scored 27 | rows 46 (refused 18) | within 3x  7/46 | out-of-sample  6/45
+ T3     panel 37 | scored 23 | rows 39 (refused 25) | within 3x  7/39 | out-of-sample  6/38
```

## T1 zero by default — **HOLDS**

Every one of the 43 rows that declares nothing is bit-for-bit identical. Exactly three
predictions moved and they are the three declared ones.

## T2 the three declared rows improve — **HOLDS, on all three**

```
2-pentylfuran  5.042 ->  64.44   fold 32.3 -> 2.53   in band False -> True
hexanal        22.84 ->  353.8   fold 34.2 -> 2.21   in band False -> True
nonanal        7.366 ->  15.61   fold  3.3 -> 1.54   in band False -> True
```

Prediction 1 (85 %) was right. Prediction 2 (60 %, nonanal enters the band) was right, and
so was the reason it was offered: these are **the first three new rows inside the 3× band in
this whole sequence of waves**, and the first accuracy gain rather than an honesty gain.

## T3 the undeclarable pots — **rule adopted, and my forecast for it was wrong**

Prediction 3 said 40 % that the rule would refuse exactly the two PratapSingh pots, and said
that if it also caught the two external isolate rows it would be **too broad and not adopted**.
It does catch them, and that reasoning was wrong. I had assumed the external rows were
answerable. They are not, and each bundle says so in provenance written long before this wave:

| pot | the bundle's own words |
|---|---|
| `external_validation_bi_2020_raw_pea_hexanal` | "this bundle is RAW pea flour, **never heated** (the 40 C / 10 min block is the HS-SPME incubation)" |
| `external_validation_liu_2023_ppi_offnote_baseline` | "rehydrated to 10 % solids in deionized water and **NEVER HEATED**; the bundle's 40 C / 10 min is the headspace equilibration" |
| `pea_isolate_40C_PratapSingh2021` | "an **UNHEATED** protein powder; the bundle's 40 C / 10 min is the HS-SPME incubation, not a thermal process" |
| `soy_isolate_40C_PratapSingh2021` | as above; "SOURCE NOT ON DISK" |

None of the four is a cook. All four measure what the raw material arrived with, and for none
of them does any source on this disk print something prior to declare. So the criterion **as
written** — "refuses exactly the rows that cannot be declared" — is met. What failed was my
guess about which rows those were, not the criterion.

**The rule needs two independent declarations to agree, because either alone is unsafe.**

1. The bundle's vessel says `closure: "no cook"`. This is a datum the bundles recorded and it
   owes nothing to any prediction — but the string is **overloaded**. Three *hot* bundles carry
   it meaning "there is no vessel to record": a synthetic snapshot, and two commercial products
   whose conditions block is a proxy operating point rather than a cook anyone ran. Keyed on
   the string alone the rule would have reached into those. It cannot be the whole rule.
2. The thermal load cannot form what was measured: the fraction of the hydroperoxide pool that
   decomposes over the program is below 1 %.

| pot | | extent |
|---|---|---|
| 40 C, 10 min | the four headspace incubations | 3.826e-3 |
| 140 C, 6 s | Trikusuma UHT, the mildest real cook | 0.2578 |
| 160 C, 25 s | Li 2026 extrusion | 0.9994 |
| 160 C, 30 min | Bi 2020 roasted pea | 1.000 |

**The threshold is not a knob.** The gap between the first row and the second is a factor of 67
and nothing in the panel lies inside it; 0.01, 0.05 and 0.10 give the identical verdict on every
row. It also survives the Q10 band end to end — at q10 = 2.0, the corner that slows the hot pots
most, the two sides are 2.824e-3 and 2.855e-2, still on opposite sides of the line.

Seven rows are refused, and they are exactly the seven the rule was described to refuse:

```
hexanal        pea_isolate_40C_PratapSingh2021                       was  3357x
2-pentylfuran  pea_isolate_40C_PratapSingh2021                       was  8524x
hexanal        soy_isolate_40C_PratapSingh2021                       was  6078x
2-pentylfuran  soy_isolate_40C_PratapSingh2021                       was 42301x
hexanal        external_validation_bi_2020_raw_pea_hexanal           was  3717x
hexanal        external_validation_liu_2023_ppi_offnote_baseline     was 33392x
nonanal        external_validation_liu_2023_ppi_offnote_baseline     was    7.3x
```

**The refusal is conditional and names its own cure.** Declare what the pot started with and
the row is answered — Trikusuma does exactly that and is scored. What is refused is the pot for
which nothing can be declared, and the honest report of that is "cannot be asked", not a
four-decade miss.

### The reason to be suspicious of this rule, stated rather than buried

Refusing rows raises the headline by arithmetic alone: within-3× goes 7/46 → 7/39 and
out-of-sample 6/45 → 6/38 **without a single prediction improving**. That is the exact shape of
a self-serving rule. Three things are offered against it, and a reader who is unconvinced should
read the scoreboard's `declared` line, which is the wave's accuracy claim with T3 switched off.

* The criterion is **condition-side** and was fixed before any error was looked at. Neither
  clause can see a measurement.
* It refuses on **two independent declarations**, one of which the bundles wrote months ago
  for an unrelated purpose.
* It **leaves every lipid miss in a pot that was cooked standing**, including the panel's
  largest: 2-pentylfuran at 366× in Li 2026, plus hexanal at 8.7×, roasted-pea hexanal at 3.7×
  and nonanal at 3.0×. A rule reaching for misses would have taken the 366×.

The Bi 2020 pair is the control that makes the case: same lab, same compound, same matrix, one
raw and one roasted. Raw missed by 3717× and roasted by 3.7×. The rule refuses the first and
scores the second.

## T4 nothing else moves — **HOLDS**

The sulfur, trunk and acrylamide lanes are bit-for-bit. No refused row was previously
un-refused for any other reason, and no row that was refused became answerable.

## Ship

**SHIPPED.** T1, T2 and T4 hold, which is the ship rule. T3 is adopted on the evidence above,
with its own suspicion recorded.

Frozen by `tests/unit/test_kinetic_core_b31.py` (13 tests), including the empty-gap and
Q10-band checks on the threshold, the overloaded-string check on the three hot `no cook`
bundles, and the guard that the largest lipid miss survives the rule.
