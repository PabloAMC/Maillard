# Pre-registration: wave B29, the first oxygen axis on the trunk (written 2026-09-10, before anything ran)

## 1. Why, and why now rather than earlier

The trunk has no oxygen axis. Every prediction it makes is at whatever oxygen the fits happened to
have, and nothing in the output says so. Two of its constants are already flagged as oxidative —
`k_glc_g` carries the flag `oxidative_entry_in_air`, and `r_ama_g` is the Amadori route to the same
glucosone — so the branch exists and simply has no lever.

Hofmann & Schieberle 2000b measures the lever. One pot, one laboratory, three atmospheres:

| | ARP-Phe | Glc + Phe |
|---|---:|---:|
| Strecker aldehyde, air / argon | **9.2×** | **3.5×** |
| Strecker acid, air / argon | 10× | 9× |
| air + Cu²⁺ / air | 2.5× | 1.9× |

**The companion paper supplies the control that makes this structural rather than empirical.** From
FED dicarbonyls the Strecker ALDEHYDE is oxygen-INDEPENDENT, while the Strecker ACID is not (air over
argon 4.3 / 5.5 / 2.2 / 1.6 from methylglyoxal / glyoxal / 3-deoxyglucosone / glucosone). So the
oxygen dependence of the aldehyde is entirely UPSTREAM of the dicarbonyl — in getting from the sugar
or the Amadori compound to it — which is exactly where this model's two oxidative entries sit. The
source's own control tells us which step to put the lever on.

## 2. The structure

A DECLARED PROCESS INPUT, not a fitted species: `atmosphere`, one of `argon`, `air` (the default,
because every fit row in this model was run in a closed vial in air) or `air_cu`. It multiplies the
two oxidative entries and nothing else:

    k_ama_g  *= f(atmosphere)      the Amadori route to glucosone
    k_glc_g  *= f(atmosphere)      the sugar route to glucosone

    f(air) = 1 by definition. f(argon) and f(air_cu) are the two fitted coordinates.

**Everything else is untouched, and that is the claim being tested.** If the oxygen dependence of a
Strecker aldehyde really is upstream of the dicarbonyl, then multiplying only these two constants
should reproduce BOTH ratios — 9.2 on the Amadori pot and 3.5 on the sugar pot — with a single
`f(argon)`. Two ratios, one coordinate: the wave is over-determined and can fail.

## 3. Fit rows, and the transfer that is declared not hidden

**IN.** Two within-study ratios: the Strecker aldehyde air / argon from the ARP pot (9.2) and from
the sugar pot (3.5). And two more for the copper arm: 2.5 and 1.9. Four ratios, two coordinates.

**THE TRANSFER, STATED.** Hofmann's amino acid is PHENYLALANINE and this model's is glycine. The
aldehyde identity differs and so does its own chemistry. What is transferred is not a yield but the
RATIO OF ONE POT TO ITSELF UNDER TWO ATMOSPHERES, and the argument that it transfers is that the
oxygen sensitivity sits in the sugar chemistry, upstream of the amino acid — which is the companion
paper's fed-dicarbonyl control, not an assumption. **If that argument is wrong the wave fails on the
second ratio**, because the two pots share `f(argon)` and differ only in which oxidative entry
dominates.

**OUT.** The absolute mol % yields. B22b was refused for exactly the mistake of taking a level from
one amino acid into a model built on another, and this wave does not repeat it.

**OUT.** The Strecker ACID ratios (10 and 9). The model has no Strecker acid, and the companion
paper shows the acid IS oxygen-dependent from fed dicarbonyls, so it is a different step with its
own oxygen term. Recorded, not fitted.

## 4. What counts as success, declared before the run

- **T1 the two air/argon ratios, decisive.** Both within **0.2 dex** on a single `f(argon)`. The
  point of the wave is that ONE multiplier explains a 9.2 and a 3.5 because the two pots weight the
  two oxidative entries differently. Fitting them with two multipliers would prove nothing.
- **T2 the copper arm, decisive.** Both within 0.2 dex on a single `f(air_cu)`.
- **T3 air is exactly 1.** With `atmosphere="air"` every prediction in the model is bit-for-bit what
  it is today. No scored panel row may move at all — not by 0.05 dex, by nothing.
- **T4 identification.** Laplace sigma below one decade on both coordinates, neither on a bound.

Ship rule: **SHIP if T1, T2, T3 and T4 hold.**

## 5. Predictions, before the run

1. T3 holds exactly. **99 %.** It is a multiplier defined as 1; if this fails something is wired
   wrong, not discovered.
2. T1 holds. **35 %, and this is the wave's real question.** For one multiplier to give 9.2 and 3.5
   the two pots must already differ in how much of their glucosone comes through the Amadori route
   against the direct sugar route, by roughly the right amount. Nothing was tuned to make that so —
   those constants come from a hazelnut paper and a milk paper — so this is close to an out-of-sample
   test of the trunk's own branching.
3. T2 holds given T1. **60 %.** Copper is a catalyst on the same step; if the step is right the
   ratio should follow.
4. The wave SHIPS. **30 %.** Two waves this week were refused for structures their sources named,
   and this one asks more of the trunk than either.

## 6. Outcome (2026-09-10, run the same day) — DO NOT SHIP, and the axis found a missing route

**Verdict DO NOT SHIP** (`kinetic_core_b29_ship_rule.md`). T2 and T3 held; T1 and T4 failed.

| test | result |
|---|---|
| T1 the two air/argon ratios | **−0.21 and +0.22 dex**, both just outside the 0.2 window. And the shape is worse than the residuals: the model separates the two pots by **1.03×** against a printed **2.63×** |
| T2 the copper arm | held, −0.16 and +0.08 dex |
| T3 air is exactly 1 | **held exactly.** A pot that declares no atmosphere is bit-for-bit what it was before this axis existed, in the parameters and in the observable |
| T4 identification | the argon multiplier identified at σ 0.25; **the copper one is not, at σ 1.6** |

**Prediction 1 held at 99 %, prediction 2 failed at 35 %,** which is roughly what 35 % should feel
like. Prediction 4 (30 % that it ships) was right to be low.

### What the axis exposed, which is worth more than the multiplier

The two ratios came out equal because they are **forced** equal. Measured directly:

| pot | share of the Strecker aldehyde made through the oxidative entries |
|---|---:|
| fed Amadori compound | **100.0 %** |
| glucose + glycine | **100.0 %** |

**The trunk's only route to the Strecker aldehyde runs through glucosone, in both pots.** So the
model's aldehyde is entirely oxygen-dependent by construction: no multiplier can give two pots
different air/argon ratios, and under argon the model goes to **zero**.

Hofmann measures **0.06 mol % from the Amadori compound and 0.04 from the sugar pot under argon** —
small, and not zero. **A non-oxidative route to the Strecker aldehyde exists and this model does not
have one.** That is a structural gap the axis exposed, not a bad multiplier, and it is the reason the
two ratios differ in the source at all: the pot with more non-oxidative background is the less
oxygen-sensitive one.

It also explains the shape of the miss. Fitting a single multiplier to two ratios that the model
forces equal lands it between them, which is exactly the −0.21 / +0.22 pair.

### What stays

The `atmosphere` input stays on `ProcessSpec` and the lever stays on the two oxidative entries,
because T3 held exactly and the axis costs nothing at air. `ATMOSPHERE_FACTORS` carries only
`air: 1.0`, so **any pot declaring argon or air + copper RAISES** rather than quietly returning the
air answer — which is the same discipline every other refusal in this model follows. A successor
needs the non-oxidative route to the dicarbonyls first; until then the axis has one setting and says
so.

## 7. Correction on review, 2026-09-10: the "missing route" in section 6 was an artefact of the observable

Section 6 says the model's Strecker aldehyde is "100 % oxygen-dependent by construction" and
concludes that "a non-oxidative route to the Strecker aldehyde exists and this model does not have
one". **That conclusion is retracted.** The run behind it measured `AKG` alone, which is the Strecker
product of GLYOXAL — one dicarbonyl's Strecker, not the Strecker aldehyde. Glyoxal comes only through
glucosone, so of course it read as fully oxidative. Hofmann's phenylacetaldehyde is the aldehyde
whichever dicarbonyl did the Strecker, and the model's equivalent is the total over both Strecker
products, `AKG + AKM`. The trunk **has** the non-oxidative route: Amadori → 1-deoxyosone →
methylglyoxal → `AKM`. Hofmann's own Figures 2 and 3 show exactly this split — the 1-deoxyosone
dominates under argon and glucosone under air with copper — and the model's branching agrees with
them in kind.

Re-run with the correct observable:

| | fed Amadori | glucose + glycine |
|---|---:|---:|
| oxidative share of the Strecker flux (model) | **54.6 %** | **64.1 %** |
| air / argon, printed | **9.2** | **3.5** |
| air / argon, model at the argon multiplier's FLOOR | 2.19 | 2.76 |

Both multipliers ran to their bounds and neither is identified. **The corrected finding is about
ORDER and about a FLOOR, and it is a sharper refutation than the wrong one.** Hofmann's Amadori pot
is the more oxygen-sensitive, so its oxidative share must be the larger; the model's is the smaller.
And the model's non-oxidative share puts a ceiling of roughly 1 / (non-oxidative share) — about 1.8
to 2.8 — on any air/argon ratio it can produce, whatever the multiplier, which is why 9.2 is out of
reach at the floor of the band. So the trunk's branching between the oxidative route to glucosone
and the non-oxidative routes to the deoxyosones is wrong in the direction the two pots differ, and
too heavily non-oxidative on the Amadori side. That is a statement about `k_ama_g` against
`k_ama_odg` and `k_ama_mgo`, and a successor should look there rather than for a route that is not
missing.

**Prediction 2's 35 % was, in the end, about the right thing** — whether the trunk's own branching
already differs between the two pots by the right amount — and the answer is no, in the wrong
direction. The first run got the verdict right for the wrong reason, and the wrong reason was
written into four places before a review caught it. The axis, T3 and the raise-on-argon behaviour
are unchanged by this correction.
