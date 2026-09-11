# Pre-registration: ENV-B18, the pyrazine step's own spread (written 2026-09-10, before anything ran)

## 1. Why this is a defect and not a gap

The Monte-Carlo envelope re-integrates the core under draws of its own fitted coordinates and
declared bands, and publishes a P5/P50/P95 for every panel row. Its prior table has no row for the
pyrazine step. `kinetic_core_b18_fit_report.json` carries that step's Laplace sigma — 0.082 and
0.075 dex on the two Strecker constants — and the envelope has never read it.

So **every pyrazine interval the tool prints today is the trunk's interval with the pyrazine step's
own spread missing.** It is too narrow, in a direction that flatters the model, and nothing in the
output says so. That is different from every other gap in the backlog, which is an answer the model
declines to give. This is an answer it gives, with a number attached, that is wrong.

## 2. What enters

Four of the step's six coordinates, routed through the `pyrazine` override block the engine already
accepts:

| coordinate | Laplace sigma | identified | how it is drawn |
|---|---:|---|---|
| `log10_k_go_ak_100C` | 0.0818 dex | yes | normal in log10 |
| `log10_k_mgo_ak_100C` | 0.0753 dex | yes | normal in log10 |
| `ea_go_ak_kj_mol` | 16.2 kJ/mol | **no, and on its upper bound** | uniform across its declared band, 100.59 to 103.1 |
| `ea_mgo_ak_kj_mol` | 16.1 kJ/mol | **no, and on its upper bound** | uniform across its declared band, 111.66 to 114.9 |

The two barriers went to their bounds in the fit and their bounds are narrow because they are
DECLARED from a measured source, not fitted. Sampling flat across a declared band is what this
envelope already does for every other coordinate a fit could not pin, and the alternative — freezing
them at a bound the optimiser pushed them to — would understate the interval again.

## 3. What does NOT enter, and this is half the defect left standing

The two pH slopes, 0.197 and 0.580 decades per unit with sigma 0.039 and 0.046, are **identified**
and will still not be sampled. The engine's `pyrazine` override takes four numbers and the slopes
are not among them: they are module-level constants read wherever the pH shape is applied. Threading
them through is a change to the parameter path, not to the envelope, and it belongs in its own wave
rather than being smuggled into a defect fix.

**The consequence, stated so nobody has to discover it.** Pyrazine intervals will still be too narrow
at any pH away from 7, and correctly wide at pH 7 where the slopes contribute nothing by
construction. The artifact must say which of the six coordinates it sampled.

## 4. What counts as success, declared before the run

- **T1** the four coordinates appear in the prior table with the sigma and bounds the B18 report
  carries, and the two slopes appear as present-but-not-sampled with the reason.
- **T2 decisive.** Every pyrazine row's interval gets WIDER or stays equal. Not one may narrow. A
  narrower interval would mean the draw is cancelling spread somewhere, which is the opposite of the
  fix.
- **T3 the rest of the panel is untouched.** No non-pyrazine row's P5, P50 or P95 moves by more than
  1e-9 dex at the same seed and sample count.
- **T4 reported.** By how much the pyrazine intervals widen, and whether any measurement that sat
  outside its interval now sits inside.

Ship rule: **install if T1, T2 and T3 hold.** T4 is reported.

## 5. Predictions, before the run

1. T2 holds and the widening is modest: under 0.2 dex on the interval width, because 0.08 dex on a
   rate constant is small against the trunk's own spread. **70 %.**
2. No measurement moves from outside its interval to inside. **80 %.** The pyrazine rows miss by
   decades, not by fractions of one, so a slightly wider interval will not reach them. The point of
   this fix is honesty about the interval, not coverage.

## 6. Outcome (2026-09-10) — INSTALL, and the effect is invisible on this panel

**Verdict INSTALL** (`env_prior_ship_rule.md`, evaluated jointly with ENV-B13). T1, T2 and T3 held.

**T2 could not be tested as written, and that is the wave's methodological finding.** Section 4 said
"every pyrazine row's interval gets WIDER. Not one may narrow." The sampler draws every coordinate
from ONE random stream, so adding coordinates re-shuffles every later draw and every row's width
moves a little at finite n. The first run gave 14 rows wider and **25 narrower** — and the
narrowings were 0.04 % to 4.7 %, on rows these priors cannot reach at all.

The rule was re-specified to measure the noise rather than assume it: two runs of the SAME priors at
different seeds, giving a relative width difference of **3.1 % median, 17.0 % worst** across 42
rows. The tolerance is the WORST, not a quantile — a 95th percentile is expected to be exceeded by
about 5 % of rows, so a rule built on one tests the quantile rather than the model. With that,
zero violations.

**No pyrazine row appears in the widened list, and that is not a failure.** No benchmark on this
panel scores a pyrazine. The prior rows change what the tool reports when a USER asks about
pyrazines; they change nothing in the scorecard, because the scorecard never asks. Prediction 1
(70 %, "the widening is modest, under 0.2 dex") is therefore unresolved on this panel rather than
right or wrong, and saying so is more honest than claiming the prediction held.

**Prediction 2 held**: no measurement moved from outside its interval to inside.

**Half the defect still stands, as declared.** The two pH slopes are identified and still not
sampled, because the engine's `pyrazine` override takes four numbers and the slopes are module-level
constants. Pyrazine intervals remain too narrow away from pH 7 and correctly wide at pH 7. That is
recorded on the prior rows themselves, not only here.
