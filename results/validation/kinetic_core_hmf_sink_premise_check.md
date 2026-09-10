# The HMF sink's centre: proposed, probed, and refused (2026-09-10)

## What was proposed

ENV-B13 gave `k_hmf_self` a five-decade band because two laboratories disagree about it, and moved
no centre. The five hydroxymethylfurfural rows on the panel are under-predicted by 2.05× to 11.9×.
The proposal was to move the sink's centre with a stated matrix — decide which laboratory's value
applies to a plant-protein cook, rather than carrying the disagreement as a band.

**Nothing was changed. One probe refutes it.**

## The probe

`k_hmf_self` ships at 8.969e-7 /min, which is **log10 = −6.047 against a band floor of −6.05**. It is
already at the bottom of its own band, three thousandths of a decade off the floor. So the only
direction a centre-move has room to go is FASTER — and a faster sink destroys more HMF, when every
one of the five rows wants MORE.

The limiting case settles it. Switch the sink off entirely — slower than any value the band could
ever license:

| row | shipped | sink OFF | change | fold error |
|---|---:|---:|---:|---:|
| `mp_holdout_fructose_asparagine_180C_Lin2022` | 1947 | 1947 | **1.000×** | 6.31 |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019` | 2.796e4 | 2.796e4 | **1.000×** | 2.05 |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019` | 2.796e4 | 2.796e4 | **1.000×** | 3.60 |
| `mp_holdout_glucose_asparagine_180C_30min_water_Chang2021` | 2192 | 2192 | **1.000×** | 3.19 |
| `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021` | 1459 | 1459 | **1.000×** | 11.93 |

**Identical to four significant figures on all five.** At the shipped value the sink carries no flux
at all, so it is not removing the missing HMF and no choice of its centre can put any back. ENV-B13's
own record said a faster sink moves these rows the wrong way; this adds the other half — a slower one
moves them by nothing, because there is no slower.

## Where the gap actually is, and it is not a sink

The last step into HMF, `k_ddg_hmf`, is one of the four constants ENV-B13 declined to band **because
two laboratories agree on it to 1.13×** — the best cross-laboratory agreement this trunk has. So
neither the making of HMF from 3-deoxyglucosone nor its removal is a plausible home for a 2–12×
deficit. What is left is the SUPPLY of deoxyosone, upstream of both.

One signal in the five rows points at something sharper, and it is visible without any fit:

* Four of the rows are sugar **+ amine** pots — fructose + asparagine, glucose + alanine twice,
  glucose + asparagine. They miss by 2.05× to 6.31×.
* The fifth, `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021`, is **555 mM glucose and
  nothing else — no amine at all** — and it is the worst of the five at **11.9×**.

The trunk reaches HMF through the Amadori compound, which needs an amine. An amine-free pot makes HMF
by acid-catalysed dehydration of the sugar, and this model has no such route: what it predicts for
that pot is whatever leaks through a network that should be near-silent. That the amine-free row is
the worst by a factor of two over the next is consistent with a missing route rather than a mis-set
constant, and it is a testable claim: a pot with no amine should be REFUSED for HMF, or given the
dehydration route with a source.

## Verdict

**REFUSED, and redirected.** The proposal cannot work: the sink is already inert. Two items for the
backlog in its place, both better-targeted than the original —

1. **The amine-free pot.** `Steinhagen2021` asks a route-less network for a number and gets 11.9× too
   little. Either refuse it (with the cure named, as B31's uncooked pots are) or supply the
   acid-catalysed dehydration route from a source. The refusal is the honest short move; the route is
   a wave with a real chance of shipping, since sugar dehydration to HMF is among the best-measured
   reactions in this whole area.
2. **The deoxyosone supply.** The four amine pots miss by 2–6×, and neither end of the HMF step is
   in doubt. That points at how much 3-deoxyglucosone the trunk makes, which is a different wave.
