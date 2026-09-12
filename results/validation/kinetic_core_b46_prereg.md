# Pre-registration: wave B46, the per-lane offset diagnostic (written 2026-09-11, BEFORE the diagnostic ran)

## 1. Why

`core_prediction_uncertainty.json` states, in its own words, that widening every parameter prior all
the way to uncapped moves coverage from 19 % to 21 % against a nominal 90 %, and concludes:

> "What is left is model-structure error -- systematic per-lane offsets a draw around a wrong centre
> cannot reproduce -- not parameter uncertainty."

That sentence names a signature and then nobody measures it. This wave measures it.

**The reasoning behind the diagnostic.** A missing *elementary step* usually perturbs a handful of
rows, because it sits on one branch. A whole lane sitting consistently on one side is the signature
of a missing *process class*, and that list is short enough to work through by hand: a sink, a source,
a catalysis, or a measurement-channel mismatch. Both channels this repository has found so far fit
that description. Wave B45's thiol finding was a catalysis; its hexanal finding was a measurement
channel. Neither was an elementary step, and neither would have been produced by enumerating
chemistry. So the question worth automating is not "what other reactions are possible" -- the
repository already enumerates 270 such steps and leaves 244 of them unused for want of a rate -- but
**"where is a lane wrong in a way a parameter cannot explain, and what does that error track?"**

## 2. What is computed

For each of the 44 scored rows that carry both a prediction and a measurement, the **signed** offset

    dex = log10(predicted / measured)

Positive means the model reads high. Every published score on this branch so far is an *unsigned*
fold error, which cannot distinguish a lane that is randomly wrong from a lane that is consistently
wrong, and only the second kind points at a missing process. Then, per lane:

- the median signed offset, and the **sign consistency**: the share of rows on the majority side;
- the spread of the signed offset;
- the Spearman correlation of the signed offset against each of **temperature, time, pH and water
  activity**, the four covariates the bundles state.

The rows split 19 sulfur, 11 acrylamide, 8 lipid, 6 trunk.

**Nothing here is a fit and nothing is scored.** No constant, tolerance or benchmark value may move
in this wave. The artifact is a diagnostic, in the same class as B38's identifiability audit.

## 3. Predictions

- **P1, the offsets are systematic.** At least two of the four lanes show sign consistency **≥ 75 %**.
  If every lane sits near 50 %, the residual is scatter rather than structure and the premise of this
  wave is wrong.
- **P2, the sulfur lane reads LOW.** Its median signed offset is **negative**. The reason is on
  record: the guide has said since before this wave that the model "removes the meaty thiols far
  faster than any real pot does, at 100 °C and at 140 °C alike", and B38 found that sink unreachable
  by refit with its barrier and both dimerisation rates on their ceilings.
- **P3, the fat lane's offset tracks TEMPERATURE.** |Spearman rho| ≥ **0.6** against temperature.
  This is the sharpest prediction here and the one most worth being wrong about. The fat lane's whole
  rate is one constant read off a graph at room temperature and carried to cooking temperature by a
  rule of thumb its source licenses only for 15–40 °C. If that extrapolation is wrong, the error must
  grow with temperature. If the offset does **not** track temperature, the Q10 is not the dominant
  fat-lane error and the section this repository wrote about it last week is over-weighted.
- **P4, the fat lane is the worst by median absolute offset**, worse than any other lane.
- **P5, at least one lane tracks at least one covariate at |rho| ≥ 0.6.** This is the wave's own
  success condition. If no lane tracks any covariate, the diagnostic found nothing, and this
  pre-registration commits to saying exactly that rather than reaching for a weaker threshold
  afterwards.
- **P6, nothing moves.** No constant, no tolerance, no benchmark value, no engine behaviour.

## 4. What is deliberately NOT claimed

A correlation here does **not** identify a mechanism. It localises where to look and rules out the
parametric explanation, which the envelope has already ruled out globally. In particular a lane whose
offset tracks temperature is consistent with a wrong barrier, a wrong Q10, a missing
temperature-dependent channel, **or** a measurement method whose efficiency changes with temperature.
Wave B45 met that last case: a purge strips a hot cell better than a cool one. The diagnostic cannot
separate those and will not pretend to.

Sample sizes are small. The fat lane has 8 rows and the trunk 6. A Spearman rho on 6 points is a hint
and is reported as one; the artifact prints n beside every correlation and the outcome section is
required to quote it.

## 5. Outcome (written 2026-09-11, after the diagnostic ran)

**Two of six predictions were refuted, and a third survived its own threshold only to fail a stricter
one I had to add after seeing the data.** The tally first, then what it means.

| lane | n | pots | median signed dex | reads | sign consistency | systematic |
|---|---:|---:|---:|---|---:|---|
| sulfur | 19 | 11 | **+0.84** | high | 79 % | YES |
| trunk | 6 | **1** | −0.90 | low | 83 % | YES |
| lipid | 8 | 3 | −0.45 | low | 88 % | YES |
| acrylamide | 11 | 7 | −0.44 | low | 73 % | no |

- **P1, the offsets are systematic — HELD.** Three of four lanes clear 75 % sign consistency.
  Acrylamide misses at 73 %. The premise of the wave survives.
- **P2, the sulfur lane reads LOW — REFUTED, and this is the wave's main finding.** It reads **HIGH**,
  by a median of +0.84 dex, and it is the thiols themselves: **2-furfurylthiol +1.48 dex (about
  thirtyfold high, 8 rows) and 2-methyl-3-furanthiol +0.78 dex (about sixfold high, 10 rows).**
- **P3, the fat lane's offset tracks temperature — HELD on the pre-registered criterion, WITHDRAWN on
  the stricter one.** rho = −0.87, which clears the threshold I wrote. But those 8 rows are **3 pots
  at 2 temperatures**. A rank correlation across two levels is a two-group comparison, not a trend,
  and dropping the known-broken alkylfuran does not help (rho = −0.88) because the problem is the
  design, not an outlier. See amendment 1.
- **P4, the fat lane is the worst by median absolute offset — REFUTED, and badly.** It is the **best**
  lane, at 0.45 dex. Sulfur is worst at 1.47, then trunk at 0.90 and acrylamide at 0.80.
- **P5, at least one lane tracks a covariate — HELD**, and it survives the stricter criterion:
  sulfur against temperature, rho = −0.82 over **6 distinct temperatures from 11 pots**.
- **P6, nothing moves — HELD.** No constant, tolerance, benchmark value or engine behaviour changed.

### The one claim with a design behind it

**The sulfur lane over-predicts the meaty thiols, and the over-prediction shrinks as temperature
rises.** rho = −0.82 across 6 temperatures, robust to dropping either thiol on its own (−0.76 without
2-furfurylthiol on 11 rows, −0.91 without 2-methyl-3-furanthiol on 9). This is the only lane in the
panel with enough distinct pots and levels to support a trend statement at all.

**It fits wave B45's finding exactly, from the other side.** B45 measured that the engine keeps
99.49 % of charged cysteine after five minutes at 95 °C where a real pot with trace copper has lost
nearly all of it. Cysteine the model fails to lose is cysteine available to make thiols. So the model
makes too much of them. The two observations are one story: **a removal channel missing upstream shows
up downstream as over-production**, and it matters most where thermal chemistry is slowest, which is
the low-temperature end, which is where the offset is largest.

**It also sharpens a claim this repository has been making for several waves.** The experiments guide
says the model "removes the meaty thiols far faster than any real pot does", and B38 found that sink
unreachable by refit with its barrier and both dimerisation rates on their ceilings. Both remain true
of the **fed** thiol. This diagnostic measures the **reacting** pots, and there the model reads six to
thirty times high. Those are not in conflict, but the branch has been framing the whole sulfur problem
as a sink that is too weak to pull the level down, and a sink already at its ceiling cannot be the
only answer. The formation side is now equally implicated, and the diagnostic says so.

### What is NOT claimed

- **The trunk lane's 6 rows are ONE pot.** Its 83 % sign consistency is a statement about that pot and
  nothing more. It is reported and it is not a lane offset.
- **The fat lane's temperature correlation is withdrawn**, per amendment 1. Its correlations against
  time and pH clear the stricter bar on 3 levels, which is the minimum and no more than that.
- **No correlation here identifies a mechanism.** Sulfur's temperature dependence is equally
  consistent with a missing low-temperature removal channel, a wrong barrier somewhere in formation,
  or a measurement whose efficiency changes with temperature. The diagnostic localises; it does not
  identify.

## 6. Amendment 1 — made AFTER seeing the first run, and labelled so it cannot hide a result

The pre-registered rule in sec. 2 called a lane/covariate pair a "track" on |rho| ≥ 0.6 alone. **That
was insufficient and the first run showed why**: the fat lane returned rho = −0.87 against temperature
from rows spanning only **two** temperatures, and the trunk returned 83 % sign consistency from a
single pot. No threshold on |rho| can distinguish a trend from a two-group comparison.

The artifact now also records, for every lane, how many **distinct pots** it rests on, and for every
correlation how many **distinct levels** of that covariate exist, and carries a second verdict
requiring at least three levels. **The pre-registered verdict is kept beside the stricter one,
unchanged**, so this amendment cannot be used to quietly delete a result: both are printed, and P3 is
recorded above as having passed the first and failed the second.

This is a change to a criterion made after seeing data, which is exactly the move this repository
distrusts. It is recorded here rather than folded into sec. 2 for that reason. What it does **not** do
is change any lane's sign, median or sign consistency, or promote any pair that did not already pass
the pre-registered test.

## 7. What follows

Nothing ships. Two things are now worth doing and neither is a refit:

- **The sulfur lane's formation side has never been examined on its own.** Every sulfur wave since B8
  has gone after the sink. The diagnostic says the reacting pots read high and the sink is already on
  its ceiling, so the next sulfur wave should test formation, not removal.
- **The panel cannot answer a design question it was never built for.** Three of four lanes cannot
  support a trend statement: the trunk has one pot, the fat lane two temperatures, acrylamide two.
  That is a benchmark-portfolio gap, not a model gap, and it belongs in the data wishlist: the corpus
  needs rows that vary one condition at a time within a lane.
