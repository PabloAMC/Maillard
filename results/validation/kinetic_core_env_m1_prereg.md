# Pre-registration: ENV-M1, one random stream per coordinate (written 2026-09-11, before the change was made)

## 1. Why

The ENV method note (2026-09-10) said the pre-registrations' rule "not one row may narrow" could not
be tested as written, because the sampler draws every coordinate from ONE stream: adding a prior row
re-shuffles every later draw, so rows a new prior cannot reach still move at finite n. The fix then
was to MEASURE the shuffle's size from two seeds of identical priors and use the observed maximum as
the floor. ENV-B34's run shows that estimator is not stable enough to carry a verdict:

| measurement | width floor, worst | median floor, worst |
|---|---:|---:|
| seed pair at ENV-B13's verdict (45 rows) | 16.99 % | 0.121 dex |
| seed pair today, same code (45 rows) | 11.28 % | 0.139 dex |
| ENV-B34's unreached rows, actual shuffle | up to 17.8 % | up to **0.178 dex** |

Seven rows the ENV-B34 priors **cannot touch** — three furfurylthiol rows on the sulfur lane and four
lipid rows, three of them in one bundle moving by an identical 0.1783 dex, the signature of a stream
shift and nothing else — exceeded a floor that another seed pair would have set above them. The rule
fired on the sampler, not on the wave, and a rule that can go either way depending on which seed pair
measured the floor is not a rule.

## 2. What is changed

`draw_from_rng` consumes exactly ONE integer from the parent generator per draw (its entropy for that
draw) and then derives an independent generator **per coordinate key** from (entropy, sha256(key)).
The two shared lipid quantiles get keys of their own; the sulfur lane's joint Laplace block gets one
key for its multivariate draw. A coordinate's stream then depends only on the draw index and its own
name: **adding, removing or reordering a prior row cannot move any other coordinate's value.**

Consequence for every envelope ship rule from now on: rows a wave's priors cannot reach are
**bit-identical** between BEFORE and AFTER, so "not one row may narrow" is testable as written, and
the measured floor becomes a sanity check rather than the verdict.

## 3. What counts as success, declared before the run

- **M1** with ENV-B34's rows filtered out of the prior list, the envelope regenerated under the new
  streams (BEFORE′) and with them in (AFTER′) agree **to the last bit** on every row none of the four
  constants can reach; the reached rows differ.
- **M2** the new streams are deterministic: two runs at the same seed are identical.
- **M3** every envelope number the guards pin moves — this is a re-seeding of the whole table — and
  is re-pinned in the same commit with this document as the reason. No coverage claim is made about
  the direction.

Ship rule: **ADOPT if M1 and M2 hold.** ENV-B34 is then re-judged on BEFORE′/AFTER′ with T2 and T3
as written in its own pre-registration; the floor is reported beside them, no longer decisive.

## 4. Predictions

1. M1 holds. **90 %.** It is a property of construction; the 10 % is an unkeyed draw I have missed.
2. ENV-B34 then INSTALLS: 3-DG inside its interval, HMF rows wider, unreached rows unmoved. **85 %.**

---

# Outcome (2026-09-11): ADOPTED

- **M1 holds at the draw level**: removing any of the `b13.`, `b18.` or `b34.` blocks leaves every
  shared coordinate bit-identical, a reversed prior table gives the same draw, and the joint sulfur
  block and the two observable multipliers are untouched (`tests/unit/test_kinetic_core_env_m1.py`).
- **M1 holds on the artifacts**: with ENV-B34's rows filtered out (`--exclude-prior-prefix b34.`) and
  in, the eight lipid rows — the only ones those constants cannot reach — are bit-identical in width
  and median. The eleven Maillard-lane rows that moved are rows the constants DO reach through shared
  glucose, which the shared stream had hidden behind noise.
- **M2 holds**: two draws at one seed are identical.
- **M3**: the whole envelope table is re-seeded once; every pinned number moved and is re-pinned in
  this commit with this document as the reason. Coverage read 14/44 under the old streams and 16/44
  under the new for the same priors; neither number is a claim about the model.

Prediction 1 (M1 holds, 90 %) — right. Prediction 2 (ENV-B34 then INSTALLS, 85 %) — right.
