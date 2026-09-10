# Pre-registration: wave B24b, the branch that refused B24 (written 2026-09-10, before anything ran)

## 1. Why, in the words of the wave that failed

B24's outcome names its own successor:

> A next pre-registration (B24b) needs hydroxyacetone as a trunk species (the Strecker of an amino
> acid on methylglyoxal makes it; Hofmann's Table 4 gives its reaction with 1-pyrroline at four pH
> values, Table 5 the intermediate's own conversion), the tetrahydropyridine as the competing
> product, and a loss of 1-pyrroline sized on experiment 3. With those, Table 9's switch is the
> within-study shape to fit.

B24 fitted the two fed-pyrroline rows (+0.18, +0.30 dex) and failed the proline chain: −1.29, +0.13,
+1.16 dex across a hundredfold methylglyoxal ladder. Its whole-chain yield rose almost linearly with
the methylglyoxal charge where the source's rises threefold, because the arm had no competing
branch and no loss of 1-pyrroline. Fed 1-pyrroline in fivefold excess, it made 41 mol % of the
methylglyoxal into the product against a printed 0.33 — 2.1 decades out.

## 2. The structure, and the one thing that makes it fittable

Two trunk species, three steps, all appended last so the existing state vectors are unchanged:

    HA     hydroxyacetone (acetol)
    ATHP   2-acetyl-1,4,5,6-tetrahydropyridine

    r_mgo_ha        MGO + amine -> HA            k_mgo_ha    (the Strecker-type reduction)
    r_pyrl_ha_athp  PYRL + HA   -> ATHP          k_ha_athp   (pH-gated, see below)
    r_pyrl_loss     PYRL        -> FRAG_C        k_pyrl_loss (first order)

**The branch is EXCLUSIVE and the source says so in words, not by inference.** Schieberle & Hofmann
2005 state that hydroxyacetone gives ONLY the tetrahydropyridine and methylglyoxal gives ONLY the
pyrroline product. So the switch is not a fitted competition between two rates on one substrate; it
is a competition for the METHYLGLYOXAL between the amino acid (which turns it into hydroxyacetone,
committing it to the tetrahydropyridine) and the pyrroline (which acylates, committing it to the
2-acetyl-1-pyrroline). That is why an excess of proline drives the tetrahydropyridine and an excess
of methylglyoxal drives the pyrroline product, and it is why the shape is a ratio and not a level.

## 3. Fit rows, and what stays out

**IN (within-study, one pot, one laboratory).** Hofmann & Schieberle 1998b Table 9's **three
AP : ATHP ratios** — 0.16, 0.51 and 12.8 at 4, 40 and 400 mmol/L methylglyoxal against 400 mmol/L
proline, pH 7, 100 °C, 30 min. A ratio inside one analysis cancels the response factor and the
extraction, which is exactly why the ratio and not the level is the target. Plus experiment 3's
**0.33 mol % of methylglyoxal** with 1-pyrroline in fivefold excess, which is what sizes the
pyrroline loss.

**A CHECK, not a fit row — see section 3b.** Schieberle & Hofmann 2005 Table 2's tetrahydropyridine
ladder: <0.1, 0.9, 10.8, 38.4 µg at pH 3, 5, 7, 9. Scored as ratios to the pH-7 rung against B18's
transferred pH term. The pH-3 rung is censored and enters as a one-sided bound, not a point.

**OUT.** The three AP levels of Table 9 as absolute mol % — B24 already showed the arm can be made
to fit any one of them and not the shape. The Schieberle 2005 time course of the final oxidation
(26 %, 72 %, >99 % over 5, 30, 120 min at 25 °C) — it is a different step at a different
temperature and this wave does not touch it. Chan 1994's apparent barrier, which B24 showed is the
trunk's dicarbonyl supply and not the step.

## 3b. Two design decisions taken while building, before anything ran

**B24's two constants are HELD at their frozen optimum and not refitted.** The question this wave
asks is whether ADDING the branch fixes the shape, and refitting the two constants that already
worked would let the fit buy the shape with them instead. Holding them makes T3 a real test rather
than a tautology: the new sinks change the 1-pyrroline pool, so the two fed rows' PREDICTIONS still
move even though their constants do not.

**The pH ladder moves from a fit target to a CHECK, and that is a demotion on purpose.** Section 3
had it as fit rows. Fitting it needs a piecewise pH slope of its own — the ladder is not log-linear,
falling 0.54 decades per unit below pH 7 and rising 0.28 above — which is two more coordinates on
two independent ratios, and a coordinate per data point is not a fit. Instead the step takes B18's
existing pH term, which was fitted on a different chemistry, and the ladder is scored as an
OUT-OF-SAMPLE check of that transfer. If B18's slopes do not describe this step, the check says so
and the number is worth more than a fitted slope would have been.

So: **two free coordinates on four rows** (three ratios plus the pyrroline-excess row), with the pH
ladder and B24's two fed rows as checks.

## 4. What counts as success, declared before the run

- **T1 the switch, decisive.** All three Table 9 ratios within **0.3 dex**, and the ORDERING strictly
  increasing with the methylglyoxal charge. Getting the ordering wrong fails outright whatever the
  residuals: the ordering is the whole finding.
- **T2 the pyrroline excess.** Experiment 3 within 0.5 dex — that is, no worse than 3x, against
  B24's 2.1 decades.
- **T3 B24's two fed rows keep their fit**, within 0.3 dex of where B24 left them.
- **T4 identification.** Laplace sigma below one decade on both new coordinates, neither on a bound.
- **T4b the pH transfer, reported.** How well B18's pH slopes describe a ladder they were not fitted
  on. A miss here is a result about the transfer, not about this wave's structure.
- **T5 nothing else moves.** No scored panel row changes by more than 0.05 dex; the sulfur, lipid and
  acrylamide lanes bit for bit.

Ship rule: **SHIP if T1, T2, T3 and T5 hold.** T4 is reported, because a coordinate that turns out
unidentified on a within-study ratio is worth knowing about but does not by itself make the shape
wrong.

## 5. Predictions, before the run

1. T1 holds on the ORDERING. **85 %.** The mechanism is a competition for one substrate and the
   ladder spans a hundredfold, so the sign of the trend is close to forced once the branch exists.
2. T1 holds on all three residuals within 0.3 dex. **45 %.** Three ratios, three free constants, but
   the ratios span 80x and the middle rung has to land as well as the ends.
3. T2 holds. **60 %.** A first-order loss is the crudest possible sink and experiment 3 is a single
   point; if it needs a second-order loss on the pyrroline this fails and says so.
4. The wave SHIPS. **40 %.** Three refused sink structures on the sulfur lane are a standing
   reminder that a named mechanism and a fittable one are different things.

## 6. Outcome (2026-09-10, run the same day) — DO NOT SHIP, and the loss is refuted, not merely unfitted

**Verdict DO NOT SHIP** (`kinetic_core_b24b_ship_rule.md`). T3 and T5 held; T1 and T2 failed.

| test | result |
|---|---|
| T1 the switch | **ORDERING HELD** — 0.596, 0.916, 1.248, strictly increasing. Magnitude did not: the model spans **2.1×** across the ladder against a printed **80×**. Worst residual −1.01 dex |
| T2 the pyrroline excess | 32.5 mol % against 0.33, **+1.99 dex** — B24 was +2.1, so this barely moved |
| T3 B24's fed rows | **+0.05 and +0.16 dex.** Adding the branch did not break what B24 already fitted |
| T4 identification | the branch constant identified (σ 0.94); **the loss unidentified at σ 8.0 and sitting on its bound** |
| T4b the pH transfer | **FAILS**: B18's slopes give +0.99 dex at pH 5 and −0.55 at pH 9 |

**Prediction 1 held (85 %): the ordering is right.** The competition for the methylglyoxal is the
correct mechanism and produces the correct sign. Predictions 2 (45 %) and 3 (60 %) failed, and
prediction 4 (40 % that it ships) was right to be low.

### The loss is refuted, and this is the wave's real result

The loss coordinate ran to its bound, which normally means "widen the band and re-run". It was tested
instead of widened, and it cannot work at any value:

| log10 k_loss | experiment 3, AP mol % of methylglyoxal | experiment 1, AP mol % of pyrroline | ratio |
|---:|---:|---:|---:|
| −1.56 (the bound) | 32.5 | 32.3 | 0.99 |
| 0.0 | 2.63 | 2.61 | 0.99 |
| +2.0 | 0.027 | 0.027 | 1.00 |
| +3.0 | 0.0027 | 0.0027 | 1.00 |

**A first-order loss on 1-pyrroline moves the two experiments together at every rate.** The source
demands they differ by **428×** per pyrroline — experiment 1 gives 28.7 mol % of its pyrroline,
experiment 3 gives 0.067 — and no value of a first-order loss produces a ratio other than 1. So
widening the band would have been chasing a fit that cannot exist, and the honest statement is that
**the loss STRUCTURE is refuted**, not that its band was too tight. The dossier had already flagged
that Table 7's rate law is not bilinear; this is that flag, sized.

### And a second finding, which is a cost of this wave's own design

The switch's magnitude is limited by a constant this wave **held**. The ratio is set by free
methylglyoxal against the hydroxyacetone the amino acid made from it, and how much gets made is
`k_mgo_pro` — held at B24's optimum by the design decision in section 3b. That decision was right as
a test and it is why the shape cannot reach 80×. A B24c should refit all four coordinates jointly,
and should expect the two B24 constants to move.

### What stays

The two species and three steps stay in the network at zero. `parameters_proline.FROZEN_B24B` is
empty and `B24B_SHIPPED` is False, so nothing is installed and 2-acetyl-1-pyrroline stays refused
with B24's verdict. What B24b adds to the record is that the branch is the right mechanism, that a
first-order pyrroline loss is the wrong sink, and that B18's pH term does not transfer to this step.
