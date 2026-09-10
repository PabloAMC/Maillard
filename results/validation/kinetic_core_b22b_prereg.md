# Pre-registration: wave B22b, the route Deng's own experiment names (written 2026-09-10, before anything ran)

## 1. Why, in the words of the wave that failed

B22 built methional as free dicarbonyl times methionine, with an identity ratio to glycine's Strecker
constants, and the data refuted it. Its outcome says where the route actually is:

> No single ratio serves the two pots: methional does not form as free dicarbonyl times methionine
> with the Strecker constants measured on fed dicarbonyls. Deng's own experiment says where it does
> form: the methionine Amadori compound alone gives 1.4 to 2.6 times more methional than methionine
> plus glucose, so the route is the Amadori compound's own decomposition (the sugar moiety supplies
> the dicarbonyl in the same molecule), which is a first-order step in a methionine Amadori compound
> the core does not carry.

That is a strong refutation and a specific successor, and the experiment that names it is a **fed**
one — Deng charged the isolated Amadori compound at 200 mmol/L and followed methional over five
times. A fed intermediate at one temperature is the cleanest kind of row this repository fits.

## 2. The structure

One trunk species, two steps:

    MARP   N-(1-deoxy-D-fructos-1-yl)-methionine, the methionine-glucose Amadori compound

    r_marp_mtal   MARP -> MTAL + fragments      k_marp_mtal   first order
    r_marp_loss   MARP -> melanoidin pools      k_marp_loss   first order

**Why two steps and not one.** Deng's series rises to 120 minutes and then FALLS: 0.244, 0.520,
1.887, 1.477, 0.822 µmol/L at 30, 60, 120, 180, 240 minutes. A single first-order decomposition of a
fed pool cannot fall — it saturates. Something removes the product or the precursor. The methional's
own removal already exists in the network (`r_mtal_msh`, B22's retro-Michael, fitted and inert), so
the second step here is the Amadori compound's own competing loss, which is what every other Amadori
compound in this model has.

**What is NOT built.** No route from methionine and glucose TO the Amadori compound. That would need
a methionine-specific glycosylamine and Amadori rearrangement, and B22 already charges methionine as
glycine for the Amadori chemistry by declaration. The Met + Glc arm of the same table is therefore a
**check**, not a fit row: does the trunk's own Amadori route, charged the declared way, land near
the 1.4 to 2.6 fold ratio Deng measures between the two arms?

## 3. Fit rows

**IN.** Deng 2022 Table 1's five MG-ARP methional points at 120 °C, 200 mmol/L, initial pH 7.5,
unbuffered: 0.244 / 0.520 / 1.887 / 1.477 / 0.822 µmol/L at 30 / 60 / 120 / 180 / 240 min. Five rows,
two free coordinates.

**Sigma 0.30 dex, and the reason is printed.** The values are semi-quantitative — headspace SPME
against dichlorobenzene with a response factor taken as 1. That is a real number with a soft scale,
and 0.3 dex is what this repository gives such a scale elsewhere.

**OUT.** The 100 and 130 °C methional series, because Deng prints them in a figure only and Table 2
has no methional row at all. So this wave has ONE temperature and **declares** its barrier rather
than fitting one; the declared value is B22's own methional-release barrier, and the declaration is
recorded on the parameter.

**OUT.** Pan 2025's rates, which refuted B22. They are a different pot with no Amadori compound in
it, and re-using them here would be asking this structure to explain the pot that killed the last one.

## 4. What counts as success, declared before the run

- **T1 the shape, decisive.** All five rows within **0.4 dex**, AND the model must PEAK between 60
  and 180 minutes. A monotone rise fails whatever the residuals: the peak is the reason there are
  two steps.
- **T2 the two-arm ratio, decisive.** The Met + Glc arm, charged the declared way and NOT fitted,
  must give less methional than the fed Amadori arm at every one of the five times — the direction
  Deng measures. The magnitude is reported.
- **T3 identification.** Laplace sigma below one decade on both coordinates, neither on a bound.
- **T4 nothing else moves.** No scored panel row changes by more than 0.05 dex.

Ship rule: **SHIP if T1, T2, T3 and T4 hold.**

## 5. Predictions, before the run

1. T1's peak requirement holds. **80 %.** Two first-order steps on a fed pool produce a peak almost
   by construction; the question is whether it lands in the window.
2. T1's residuals all within 0.4 dex. **55 %.** Five points, two constants, and the last two points
   fall by 2.3x while the first three rise by 7.7x — a shape with a sharp turn.
3. T2 holds on direction. **40 %, and this is the one I expect to fail.** The trunk's Amadori route
   was fitted on glycine, and B22 already showed that charging methionine as glycine overshoots
   Deng's pot by 1.5 to 2.9 decades. If the declared substitution is wrong, the check says so, and
   that is worth more than a pass.
4. The wave SHIPS. **30 %.**

## 6. Outcome (2026-09-10, run the same day) — DO NOT SHIP, and a fed first-order pool cannot make this shape

**Verdict DO NOT SHIP** (`kinetic_core_b22b_ship_rule.md`). T3 and T4 held; T1 failed; **T2 could not
be evaluated at all**.

| test | result |
|---|---|
| T1 the shape | worst residual +0.51 dex, and the model **peaks at 30 minutes** against a printed 120. The model is essentially FLAT at 0.78 µmol/L |
| T2 the two-arm check | **VACUOUS** — see below |
| T3 identification | both coordinates identified, σ 0.09, neither on a bound |
| T4 nothing else moves | panel bit for bit |

### T2 was ill-posed, and that is my error rather than a result

Section 4 said the binary Met + Glc arm "must give less methional than the fed Amadori arm", and
section 2 said it "rests on charging methionine as glycine for the Amadori chemistry". **It rests on
nothing.** Wave B22 did not ship, so its four steps are inert at zero and the model has no route
from methionine and glucose to methional at all. The binary arm is structurally zero and the ratio
divides by nothing — the fit report duly prints 7.8e29. That is a division by zero dressed as a
finding, and the pre-registration should have noticed that its comparison arm was switched off
before making it a decisive test.

### The structural result, and it was tested rather than assumed

The model's flatness is not a fitting failure. Methional has no sink here — the only one in the
model is B22's retro-Michael, inert because B22 did not ship — so the first attempt was asked to
make a rise-and-fall with nothing that can fall. Two probes:

**A methional sink alone never produces the peak.** Sweeping it from inert to 0.1 per minute at the
fitted constants:

| methional sink | 30 | 60 | 120 | 180 | 240 min | peak |
|---|---:|---:|---:|---:|---:|---|
| inert | 0.78 | 0.78 | 0.78 | 0.78 | 0.78 | 30 min |
| 1e-3 /min | 0.69 | 0.61 | 0.48 | 0.38 | 0.30 | 30 min |
| 1e-2 /min | 0.23 | 0.07 | 0.006 | 0.001 | 0.000 | 30 min |

It only makes the curve fall from the start.

**Nor can three constants.** A grid over the Amadori decomposition, its competing loss AND the
methional sink — 4 decades × 4.5 decades × 3.5 decades — gives a best of 0.625, 0.618, 0.603, 0.588,
0.574, peaking at 30 minutes, at a cost WORSE than the two-coordinate fit. Against a printed
0.244 → 1.887 → 0.822.

**Why, and it is a statement about first-order kinetics rather than about this chemistry.** The
source's series RISES 7.7× between 30 and 120 minutes. A first-order decomposition of a FED pool is
fastest at t = 0 and decelerates from there — it is concave from the first instant and cannot show an
induction period. Producing one needs a SEQUENTIAL intermediate between the Amadori compound and the
aldehyde, or an acceleration the model does not have. **So the refutation is of the one-step reading
of Deng's own proposal, not of the proposal.** The Amadori route may still be right; a single
first-order step from it is not.

### What stays, and what a successor needs

One species and two steps stay in the network at zero; `FROZEN_B22B` is empty and `B22B_SHIPPED` is
False. A B22c needs (i) an intermediate between the Amadori compound and methional, to buy the
induction period, and (ii) a comparison arm that is not switched off — which means it cannot be run
until something supplies a methionine-to-methional route at all, and B22's is refused.
