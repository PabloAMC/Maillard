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
