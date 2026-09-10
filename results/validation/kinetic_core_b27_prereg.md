# Pre-registration: wave B27, the oxidant the three refused sink waves were never given (written 2026-09-09, before anything ran)

## 1. Why, and what changed

Three thiol-sink structures have been refused on the same 64 rows: B17 variant (b), the reversible
disulfide; B17 variant (a), the saturable thioether; B25, irreversible addition to the pot's own
deoxypentosones. B25's record concluded that the objective's weighting was the problem and named two
remaining routes, a laboratory experiment or an owner decision to re-weight.

The reading audit changes that conclusion, in two steps.

**First, a fact about the objective, verified by enumerating it rather than by reading a comment.**
The oxidative dimerisation channel `ch_dimer_*` is first order in an explicit oxidant pool `OX`.
That pool is charged at the ambient convention of 1.0 in 24 of the 41 systems that carry a fit row,
at a cystine equivalent of 62.4 in one, and at **exactly zero in all fourteen fed mol % rows** — the
twelve at 145 °C and the two Whitfield pots at 140 °C. The headspace reservoir `OXR` is zero in
**every** system in the repository. So the channel that makes disulfide carries structurally zero
flux in precisely the rows that dominate the objective, while both dimerisation constants sit pinned
at the top of their bands and the model still makes ten to two hundred times too little disulfide.
Three sink structures were asked to reproduce a disulfide share in pots where disulfide formation
was arithmetically impossible.

**Second, a named oxidant.** Whitfield & Mottram 1999, read for the first time on 2026-09-09, argues
the same thing from the bench and proposes what supplies the equivalents: **the pot's own
α-dicarbonyls, reduced to hydroxyalkanones on the way to the mercaptoketones**. In that pot the
mercaptoketone flux is 244 µg per 10 mg of norfuraneol against 15 µg of free 2-methyl-3-furanthiol,
so if each reduction delivers one oxidising equivalent the supply scales with a flux **sixteen times**
the thiol's own.

## 2. What is defensible about the zero, and what is not

This must be said before anything is changed, because half of the zero is right.

**Defensible.** Whitfield argues explicitly that aerial oxidation is unfavourable in these pots
*because of the high hydrogen sulfide concentrations*, and that is exactly why the paper needs a
different redox system at all. Every fed pot charged with H₂S therefore has a chemical reason to
carry no ambient oxidant, and this wave does not give it one.

**Not defensible.** `fed_nf_cys_MFT` and `whitfield_nf_cys_MFT` are **cysteine** pots with no H₂S,
and they carry zero while the Hofmann pentose + cysteine pot in the same table, from the same
laboratory, carries the ambient 1.0. Nothing in the generator states a reason for the difference —
the `MELE` charge two lines above it names exactly which pots carry it and why, and the oxidant
charge does not. That is an inconsistency between two comparable pots, not a declaration.

## 3. The structure

Two changes, both minimal, both inert at their defaults.

**(a) A source step.** No reaction anywhere makes an α-dicarbonyl from norfuraneol, so a
dicarbonyl-sourced oxidant would be identically zero in the very pot that has to show the disulfide
share. The paper **measures** three diketones in that pot (Table 1 rows 1–3: 2,3-pentanedione,
2,4-pentanedione and 3,4-hexanedione), so the substrate is real and the step is missing. One step,
`r_nf_dicarbonyl`, on one fitted constant with a declared barrier.

**(b) The redox couple.** One step, `ch_redox_dicarbonyl`, second order, consuming an α-dicarbonyl
and delivering oxidant equivalents to `OX`, on one fitted constant and one declared barrier. Zero by
default, so every earlier wave reproduces bit for bit.

**What this wave does NOT do.** It does not charge ambient oxygen into any H₂S pot. It does not
re-weight the objective; the row weights are untouched, which is what makes the result readable
against the three refusals. It does not touch the `k_thioether` pool, whose own measured capacity is
about 320× below what the `MELE` pool is sized on and which is its own backlog item.

## 4. Fit rows

**One new row, and it is the first of its kind on this lane.** Whitfield 1999's MFT disulfide share,
about 35 % of total MFT at 140 °C and pH 4.5, computed by summing the paper's own printed disulfide
rows as MFT equivalents. It is a within-study ratio in a pot the objective already simulates, and it
is the only disulfide constraint anywhere outside 115–120 °C. It travels with its caveat: the paper
assumes unit response factors, and disulfides are far less volatile than thiols under a 60 °C
dynamic headspace, so **35 % is a lower bound.**

**Three charge corrections that must ride in this wave**, because the generators are frozen by
manifest and an edit outside a wave is not allowed:

| what | carried | printed |
|---|---:|---:|
| Whitfield norfuraneol | 20 mmol/L | **50** |
| Whitfield hydrogen sulfide | 40 mmol/L | **~97** |
| `whitfield_nf_cys_MFT` target | 0.150 mol % (**free** MFT) | free 0.150, **total 0.230** |

The concentrations are 2.5× low and the H₂S steps are second order, so the error propagates. Sized
on the shipped B9 vector the predicted yield moves 1.12×, which is 0.05 dex against that row's own
sigma of 0.5 — a tenth of one sigma, so B9 stands and this is a correction, not a rescue. The target
correction is 0.19 dex, a third of a sigma. The pot's buffer is also upgraded from ASSUMED to the
printed 0.5 M phosphate at pH 4.5.

## 5. What counts as success, declared before the run

- **T1** the three refused waves' own T1: the reference pot at 100 °C within 0.5 dex of Table IV as
  ratios to the 30-minute point, and both thiols still rising between 6 and 12 hours.
- **T2** every B9 fit row within 0.3 dex of its B9 residual.
- **T3 the decisive one, and it is new.** Zhou 2023's three dimer shares within 0.3 dex. All three
  refused waves failed this by 2.4 decades and every one of them failed it the same way, which is
  the signature of a channel that cannot run rather than of a constant that is wrong.
- **T4** the new Whitfield share row within 0.3 dex, read as a **lower bound**: over-prediction is
  not a failure of this row and is reported.
- **T5** Yiltirak's median fold below 20 and Wang's 140 °C decline under one decade, as before.
- **T6** both new coordinates identified: Laplace sigma below one decade, off their bounds, slices
  not flat. All three refusals failed T6 with sigma between 2e4 and 9e4 decades on a flat slice.

Ship rule: **SHIP if T1, T2, T3 and T6 hold.** T4 and T5 are reported.

## 6. Predictions, before the run

1. T3 improves by more than one decade on at least two of Zhou's three shares. **60 %.** The
   argument for it is that the channel currently cannot run at all in the dominant rows; the
   argument against is that Zhou's own pots already carry the ambient 1.0, so their shares are
   short for a different reason and this wave will not touch them.
2. T6 holds, unlike all three refusals. **65 %.** The earlier coordinates were unidentifiable
   because the objective had no gradient on them; a coordinate that feeds a channel with a fit row
   attached should have one.
3. The wave does **not** ship. **50 %.** T2 is the historical killer: every structure that changed
   the 100 °C behaviour moved a fed row by 1.9 dex. If this one ships it will be because it changed
   what the dimer channel does without changing what the thiols do.
4. Whatever the verdict, the record gains something the three refusals did not have: a statement of
   whether the sink question was ever a sink question. **95 %.**

## 7. Status

**PRE-REGISTERED, NOT RUN.** Written 2026-09-09 with the fit deliberately not started, so the
declaration above cannot be revised after seeing an outcome. The two structural facts it rests on
were verified by enumeration, not quoted: the oxidant charge per system and the objective's row
composition.

## 8. Amendment, 2026-09-09 evening: four papers read, STILL NOT RUN

Four papers were fetched and read specifically to anchor this wave's one unanchored constant. The
amendment is written before any fit has been started, and the fit is still not started.

**The direct answer: not one of them supplies a rate, an order or a barrier for the oxidation of a
thiol to a disulfide by an alpha-dicarbonyl.** Every constant in all four belongs to the ADDUCT
branch, and in three cases the adduct was structurally characterised, so that is not an inference.

| paper | what it actually supplies | branch |
|---|---|---|
| Wondrak 2002 | four MEASURED second-order constants at 37 C, pH 7.4 (L-cysteine + phenylglyoxal **0.63 +/- 0.04 M^-1 s^-1**, D-penicillamine 24.8, and two more); one temperature, so no barrier | adduct: 2-acylthiazolidines, isolated, NMR and MALDI-TOF. The word "disulfide" never appears |
| Zheng 2022 | no constant; extent only. Glutathione falls about 18.1 / 8.6 / 1.7 % against methylglyoxal / glyoxal / 3-deoxyglucosone, instantly at 1:1, then flat for six hours | adduct, four conjugates by LC-TOF-MS -- plus the batch's only direct disulfide measurement, and it is unquantified |
| Zheng 2023 | prints 4.1e4 M^-1 s^-1 and a reverse 7.5e-3 s^-1, but AS CITATIONS to another laboratory, conditions not restated; measures no constant itself | adduct: hemithioacetal formation and hydrolysis |
| Coukos 2021 | nothing quantitative. An ordering, and a turnover at 2-6 h at 25 C. Every axis is a peak area normalised within its own series, so no constant could be refitted from it | adduct, proven by mass arithmetic |

**Two cautions that now belong to the corpus.** The two thiol-times-dicarbonyl second-order constants
it holds -- Wondrak's 0.63 and the constant Zheng 2023 cites from Lo 1994 at 4.1e4 -- differ by
**6.5e4** and are not the same object. And that cited forward/reverse pair implies an equilibrium
constant of 5.5e6 M^-1, which predicts essentially all of the glutathione bound at Zheng 2023's own
1:10 ratio against a measured 34.7 %, an inconsistency of at least two decades that the paper does
not remark on.

**How this moves the odds, and why less than it first appears.** Section 6's prediction 3 ("the wave
does not ship", 50 %) rises a few points. It does not rise twenty, for four reasons.

1. Zheng 2022's "limited" is never given a number. The disulfide channel was instrumented with four
   transitions against an authentic standard and the result is one adjective and "data not shown".
   That adjective spans at least an order of magnitude in the quantity that decides this.
2. **The same group partly reverses it.** Over 48 hours, with the disulfide on its own calibration
   curve, Zheng 2023 finds glutathione going to its disulfide as the DOMINANT sink. Both papers
   attribute that to oxygen autoxidation, and neither ran the dicarbonyl-free blank that would test
   the attribution. The control the whole question turns on was not run, twice.
3. A null at 37 C is the weakest kind of null for a route expected to carry the higher barrier, and
   no paper here has a second temperature.
4. **A structural mismatch that cuts FOR this wave, not against it.** All three dicarbonyls Zheng
   tested are alpha-OXOALDEHYDES: they carry an aldehyde carbon, and that is the carbon that makes
   the hemithioacetal which won. The species section 3(a) would source from norfuraneol --
   2,3-pentanedione, 2,4-pentanedione, 3,4-hexanedione, the three this pot's own paper measures --
   are **alkyl diketones with no aldehyde carbon at all**. Remove the aldehyde and you remove the
   branch that beat the disulfide.

**A BOUND DECLARED BEFORE THE FIT, which is the useful thing this reading produced.** The redox
constant has no anchor, so the fit could put it anywhere its band allows and call the result an
estimate. It cannot. The prior these papers establish is that the redox branch is MINOR against the
adduct branch at 37 C on an oxoaldehyde. Accordingly: **if the fitted `ch_redox_dicarbonyl` constant
lands within two decades of Wondrak's measured adduct constant of 0.63 M^-1 s^-1, transported to the
lane's reference temperature, the wave must report that the objective has put a redox step at or
above the speed of the adduct step that outcompetes it -- and treat that as evidence against the
structure rather than as a fitted value.** T6's identification test is not sufficient on its own for
a coordinate with no independent measurement anywhere.

**One more route, noted and structurally unavailable here.** Coukos's mercaptomethylimidazole
crosslink is a fourth, irreversible carbon sink for a thiol. It needs a guanidine, and the pot that
must show the 35 % disulfide share has none, so it cannot run there.

## 9. Amendment, 2026-09-10: four probes, STILL NOT RUN, and the wave's structure is vindicated by them

Four probes on the shipped B9 vector. **No fit was started and no constant moved.** They were run to
decide whether this wave is worth building, and they sharpen its case in a way section 1's
enumeration could not.

### 9.1 The enumerated premise re-verified, independently

35 systems, 58 active rows. **OX is charged in 18 systems and is exactly zero in 17**, including all
**fourteen** fed mol % rows. `OXR` is zero in every system in the repository. Section 1's claim holds
to the row.

### 9.2 The dimer channel is not broken — it is starved, and almost exactly first order in oxidant

| dimer-share row | its pot's OX | target | shipped | ×10 | ×100 |
|---|---:|---:|---:|---:|---:|
| `zhou_pH7_dimer_over_MFT` | 1.0 (ambient) | 0.0323 | 0.00294 | 0.0292 | 0.275 |
| `zhang_fig1_cys_dimer_over_MFT` | 1.0 (ambient) | 0.0429 | 0.00487 | 0.0482 | 0.441 |
| `zhang_fig1_gcys_dimer_over_MFT` | **62.4 (a measured cystine charge)** | 0.2711 | **0.285** | 2.03 | 9.26 |

Two things follow, and the second is new.

**The shipped misses are 11×, 8.8× and 0.95× — not the 2.4 decades the refused waves recorded.** Those
failures were measured at each refused wave's OWN fitted vector, which had moved the thiol chemistry;
at the shipped vector the deficit is one decade, not two and a half.

**The one pot whose oxidant is not a convention gets its disulfide share right.** `zhang_fig1_gcys` is
charged from a MEASURED cystine loading of 62.4 and lands at 0.95× of target. The two pots that are
short are the two on the ambient convention of 1.0. That is the strongest evidence this question has
that the deficit is an oxidant-BUDGET problem and not a broken channel — and it comes from the
corpus, not from an argument.

### 9.3 The obvious cheap fix is refused, and refusing it is what makes the case for this wave

The one-line change would be to raise the ambient convention. Measured over the whole objective:

| ambient × | total cost | rows worse by > 0.3 dex |
|---:|---:|---:|
| 1 | 587.3 | — |
| 3 | **583.4** | 0 |
| 10 | 591.1 | 1 |
| 30 | 616.3 | 4 |
| 100 | 670.5 | 4 |

At ×10 the two short dimer rows go to **−0.11 and +0.13 sigmas** — essentially exact, a 1.0 and 0.89
decade gain, which is section 6's prediction 1 met twice over. And the bill arrives immediately:

```
zhou_pH7_dimer_over_MFT         -0.997 dex   residual -2.60 -> -0.11
zhang_fig1_cys_dimer_over_MFT   -0.894 dex   residual -2.36 -> +0.13
kumazawa_FFT_retention_pH5_4    +0.296 dex   residual -0.69 -> -2.16
kumazawa_FFT_retention_pH6_0    +0.590 dex   residual -0.44 -> -3.39
```

Kumazawa's rows are the purest measurement in the corpus: **1 ppm of 2-furfurylthiol alone in buffer,
heated 121 °C for 10 min, scored as the heated half over the unheated half of the same solution**, with
no formation chemistry anywhere. Raising a global oxidant makes more disulfide there too, and the
thiol stops surviving. So a global oxidant cannot be the answer, and **the ambient convention is not
quietly changed by this wave or any other**.

### 9.4 Why that refusal is an argument FOR this wave's structure, not against it

The two measurements are not in conflict; they are in different pots, and the difference between the
pots is exactly the one this wave's oxidant is selective on.

* Zhou's and Zhang's pots are **Maillard pots full of α-dicarbonyls**.
* Kumazawa's pot is **a pure thiol in buffer with no dicarbonyl in it at all**.

A dicarbonyl-sourced oxidant (`ch_redox_dicarbonyl`) supplies equivalents in the first and **exactly
none in the second, by construction**. That is the whole difference between the change refused in 9.3
and the change proposed in section 3, and it was not part of the argument when this wave was written.

**A prediction this adds, declared now and before the fit.** T2 (every B9 fit row within 0.3 dex of
its B9 residual) is the historical killer, and 9.3 shows precisely which rows would kill it —
Kumazawa's four. Because those pots carry no dicarbonyl, **this wave should move them by nothing at
all**, and if a fitted `ch_redox_dicarbonyl` does move them it will mean the step is being used as a
disguised global oxidant. **T2 is therefore re-read as a targeted test rather than a broad one: the
four Kumazawa rows must move by less than 0.05 dex.** Probability that they do: **80 %.**

Section 6's prediction 1 (T3 improves by more than one decade on at least two of Zhou's three shares)
rises from 60 % to **80 %**, because 9.2 measures the channel's oxidant response directly and 9.3
shows a decade is available on exactly two of the three. Prediction 3 (the wave does not ship, 50 %)
falls to **40 %**.

**STILL NOT RUN.** Written before the fit was started.

## 10. Amendment, 2026-09-11, written BEFORE the fit is started: the structure is simplified to ONE coordinate

Section 3 proposed two steps and two fitted constants: a source step `r_nf_dicarbonyl` making an
α-dicarbonyl from norfuraneol, and a second-order redox couple consuming it to deliver oxidant
equivalents. Reading the lane before building shows that is one step too many, and the reason is
in section 1's own sentence: **the oxidant supply "scales with a flux sixteen times the thiol's" —
the mercaptoketone flux — and the lane already carries that flux as `r_nf_mp3p`.** Whitfield's
Figure 6 says exactly this: the reduction of the dicarbonyl on the way to the mercaptoketone IS the
redox system. A separate dicarbonyl species with its own source rate would put a coordinate in the
vector that nothing measures (only the product of the two constants would be observable), and the
identification test T6 would then be failing on a coordinate the structure never needed.

**What is built instead.** One step, `ch_redox_mp3p`: NF + H₂S → MP3P + OX, running in parallel to
`r_nf_mp3p` with rate `φ · k_nf_mp3p`, while `r_nf_mp3p` runs at `(1 − φ) · k_nf_mp3p`. Same
barrier, same pH factor (`neutral_h2s`), so **the total mercaptoketone flux is unchanged for every
φ** — the 16.3 : 1 ratio row cannot move — and each mercaptoketone formed through the redox branch
delivers **one oxidant equivalent** to `OX`, in the pool's own units, where one equivalent makes one
disulfide. **φ is the one fitted coordinate**, the fraction of mercaptoketone-forming events that
oxidise a thiol on the way. It is a branching ratio and has no barrier of its own. At φ = 0 every
earlier wave reproduces bit for bit. Band: log10 φ ∈ [−4, 0].

**The declared ceiling replaces section 8's Wondrak bound.** φ cannot exceed 1: one dicarbonyl
reduction is one disulfide. **If the fit pins φ at its ceiling, the objective is asking for more
oxidant than the mercaptoketone flux can supply, and that is evidence AGAINST the structure, to be
reported as such and not as a fitted value.** Section 8's bound on a second-order redox constant no
longer applies because no such constant exists in the vector.

**T6 is re-read for one coordinate**: Laplace sigma on log10 φ below one decade, off both bounds,
slice not flat.

**Two row kinds are added to the objective's vocabulary**, both inert on every existing row:
`molpct_total` (a mol % over a list of species with multiplicities, so total MFT = MFT + 2 × dimer)
and `floor` (one-sided: no penalty ABOVE the target), for the Whitfield share row, which section 4
declared a lower bound. The three charge corrections of section 4 are installed by this wave's own
generator at import, with the original values restored on exit, so the frozen B2.3 generator is not
edited: norfuraneol 50, cysteine 50, H₂S 97 mmol/L; the buffer declared as the printed 0.5 M
phosphate at pH 4.5; the `whitfield_nf_cys_MFT` target as the printed total, 0.230 mol %, over the
printed basis of 50.

Predictions stand as revised in section 9. **Still not run at the time of writing.**

---

# Outcome (2026-09-11): the gate fired before the fit. NOT FITTED, DO NOT SHIP.

`kinetic_core_b27_ship_rule.md` is the record. The structure of section 10 was built (`sulfur.py
ch_redox_mp3p`, `parameters_sulfur.apply_dicarbonyl_redox`, the engine hook, the two row kinds, the
charge corrections installed by the generator) and the fit was **not started**, because two probes at
the shipped B9 vector settle the ship rule before any constant could move:

**G1 — the decisive test T3 is unreachable, by this structure or by any oxidant source.** In Zhou's
and Zhang's pots the dimer step consumes **0.001 %, 0.49 % and 0.10 %** of the charged oxidant. The
pool there is not a budget; it is a constant multiplier on the dimer rate. The entire mercaptoketone
flux at φ = 1 would raise it by **0.35 %, 0.41 % and 0.01 %**, against the roughly ninefold a decade
on the share needs. What those rows want is the dimer RATE CONSTANT, which is at its band ceiling and
is opposed by Kumazawa's retention rows. Since the rule is SHIP only if T3 holds, no fit result could
have shipped.

**G2 — in the one pot the structure does fix, φ must sit at its physical ceiling.** Whitfield's
cysteine pot, whose oxidant is genuinely zero, climbs from 0 % disulfide-bound MFT to 5.5 % at φ = 0.1,
15.4 % at φ = 0.316, 22.3 % at φ = 0.5 and reaches the 35 % floor **only at φ = 1.0 exactly** — every
mercaptoketone-forming event oxidising a thiol. Section 10 declared that outcome disqualifying before
the fit, and it is honoured.

**G3 — the fit-free T2 held, and section 9's targeted prediction held.** At φ = 1 no row of the 65
worsens by more than 0.3 dex (worst +0.00), and the four Kumazawa rows move by **0.0000 dex**, as
predicted at 80 %, because they carry no norfuraneol. T1 and T5 fail at φ = 1 and fail identically
under B9 itself; φ changes them in the fourth decimal.

## Section 9 was wrong about what it had found, and this is the correction

Section 9 said the disulfide deficit was "an oxidant BUDGET problem", on two observations: the pot
charged with a measured cystine loading of 62.4 landing at 0.95× of target, and a ×10 ambient oxidant
fixing the two short rows exactly. Both observations are right and the reading was wrong. **The
consumers use under one per cent of the pool in every ambient pot**, so 62.4 is a 62-fold multiplier
on the rate and ×10 is a tenfold one — neither is a budget being met. The deficit in the ambient pots
is a RATE deficit that the oxidant term happens to multiply. Only in the fed pots, where OX is exactly
zero and the dimer flux is exactly zero, was the budget reading correct — and there the structure works,
but only at its ceiling. Prediction 1's revision from 60 % to 80 % rested on the wrong reading and
should not have been made.

## What the sink question was, then (prediction 4, the 95 % one)

It was two questions wearing one name. In the **fed** pots it was an oxidant question: three sink
structures were scored on rows where disulfide formation was arithmetically impossible, exactly as
section 1 said. In the **ambient** pots it never was: the channel has all the oxidant it can use and
makes ten times too little disulfide anyway, and the only lever is a rate constant already pinned by
the purest measurement in the corpus. That second half is the finding this wave leaves behind.

## What is kept, and one named debt

The step, the parameter and the hook stay, inert at φ = 0, as B17's and B25's refused structures do.
The three Whitfield charge corrections are installed only by this wave's generator and restored on
exit, so the shipped objective is unchanged. **That is a debt**: the printed charges (norfuraneol 50,
cysteine 50, H₂S 97 mmol/L; 0.5 M phosphate at pH 4.5; total MFT 0.230 mol %) are facts, and the next
sulfur refit must carry them.

Predictions, scored: 1 (T3 improves by a decade on two shares) — **wrong**; 2 (T6 holds) — not tested,
no fit; 3 (does not ship, 50 % then 40 %) — **right, for a reason the wave did not anticipate**; 4
(a statement of whether the sink question was a sink question) — **delivered**, and the answer is
"half of it".
