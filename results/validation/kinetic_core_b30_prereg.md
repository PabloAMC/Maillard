# Pre-registration: wave B30 (W8), pH on the thiol FORMATION steps (written 2026-09-10, before anything was fitted)

## 1. Why, and the premise checked first

The backlog's W8 says the corpus holds three pH contrasts on the pentose-cysteine path and "the lane
has no pH term on it". **The second half needed checking before anything was built, and it is
wrong.** The lane has a pH mechanism on formation, and a principled one: B2.1 gave the sulfide two
protonation states and let both add, on the argument that "HS⁻ is by far the better nucleophile and
H₂S dominates only by abundance below the pKa". The note on those steps says the crossover "is a
rate ratio, not a pH parameter", and that this was "half of why B2's thiols collapsed at high pH".

So the question is not whether a pH term exists. It is whether the one that exists is right.

## 2. The decisive test, declared before it was run

Whitfield's fed norfuraneol pot is the only place in the corpus where the SAME pot is measured at two
pH values by the same laboratory: pH 4.5 in Whitfield 1999, which is already a fit row, and pH 6.5 in
Whitfield 2001. Free 2-methyl-3-furanthiol falls from **0.150 to below 0.001 mol %**, at least
**150× down**; the mercaptoketones fall from 74.5 to 0.03, about **2500× down**.

**T0, and it gates everything else: the model must get the SIGN right.** Run the shipped lane on that
pot at pH 4.5 and 6.5. If the ratio is not a fall, no pH slope may be fitted on top, because fitting
a slope over a wrong-signed mechanism buys agreement with the mechanism still wrong.

Only if T0 passes do the rest follow: T1, Cerny 2007's five-point ladder as ratios; T2, Mottram
2002's pyrophosphate pair; T3, identification; T4, nothing else moves.

## 3. Outcome (2026-09-10) — T0 FAILS, and nothing was fitted

| arm | pH 4.5 | pH 6.5 | ratio 4.5 / 6.5 |
|---|---:|---:|---:|
| **as shipped** | 0.0383 | 0.0777 | **0.49** |
| thiolate loss switched off | 0.0384 | 0.0990 | 0.39 |
| hydrosulfide branch switched off | 0.0383 | 0.0777 | **0.49, unchanged** |
| **measured** | 0.150 mol % | < 0.001 | **≥ 150** |

**The model gives twice as much thiol at the higher pH where the measurement collapses by at least
150-fold.** The sign is wrong and the ratio is out by about 300×. No slope was fitted, by the rule
above.

### Two attributions, and the second is the useful one

**The thiolate loss is not the answer.** Switching it off moves the ratio from 0.49 to 0.39 — it
pushes the right way and carries almost none of the effect.

**The hydrosulfide branch is not even active here.** Switching it off changes the answer by nothing
at all, to four figures. The reason is structural: the two-branch mechanism was built for the
deoxypentosone route (`r_ddp_mft_hs`) and the furfural route (`r_fur_fft_hs`), and **the norfuraneol
route has no hydrosulfide partner** — `r_nf_mft` and `r_nf_mp3p` are single steps in neutral H₂S.
So on the one pot in the corpus where a pH pair is actually measured, the lane's pH mechanism is
absent, and what is left of the pH response comes only from how much sulfide the cysteine releases.

### And adding the missing partner would make it worse, which is why this is not a two-line fix

Giving the norfuraneol steps a hydrosulfide branch would push the SAME way as the branch that is
already there: more hydrosulfide at higher pH means faster addition means MORE thiol at pH 6.5. The
measurement wants less. So the collapse is not in the nucleophile at all. It is in the substrate or
in the sulfide budget — norfuraneol's own enolisation and stability, or the sulfide going somewhere
else as it deprotonates — and neither is in this model.

## 4. What this wave leaves

Nothing is installed and nothing is fitted; the lane is untouched. What it produces is a sharper
statement of the gap than the backlog had:

* the lane's pH mechanism on formation is a **nucleophile-speciation** mechanism;
* it is **structurally absent** on the norfuraneol route, which is the route the only measured pH
  pair runs on;
* and its **sign is wrong** for that pair, so supplying the missing branch would move the model
  further from the measurement, not closer.

A successor needs a pH term on the SUBSTRATE side — norfuraneol's enolisation, or the sulfide
budget — and it needs it before any slope is fitted to Cerny's ladder or Mottram's pair, because
those two would otherwise be fitted on top of a mechanism this test has already refuted.
