# Where the sulfur lane's carbon and sulfur go: a flux map, and what it found

*For a reader who wants to know which reactions actually carry the chemistry in the sulfur lane — not
which ones exist, but which ones the cysteine, the pentose and the sulfide flow through, how fast, and
whether that agrees with what was measured. Written 2026-10-28. The numbers come from
`scripts/generators/build_sulfur_flux_attribution.py`
([results/validation/sulfur_flux_attribution.md](../../results/validation/sulfur_flux_attribution.md)),
which integrates every reaction's rate over six reference pots and checks that the budget closes.
Companion documents: [REACTION_TREES.md](REACTION_TREES.md) (which steps exist),
[PH_BUFFERS_AND_WATER.md](PH_BUFFERS_AND_WATER.md) (how the pH is computed).*

1. [The question and the method](#1-the-question-and-the-method)
2. [The paths, as the model runs them](#2-the-paths-as-the-model-runs-them)
3. [Finding 1: the pH trajectory was not converged (B60, fixed)](#3-finding-1-the-ph-trajectory-was-not-converged-b60-fixed)
4. [Finding 2: the trunk runs in seconds; the pots say hours (B59, open)](#4-finding-2-the-trunk-runs-in-seconds-the-pots-say-hours-b59-open)
5. [Smaller findings](#5-smaller-findings)
6. [Why nobody saw it, and the general lesson](#6-why-nobody-saw-it-and-the-general-lesson)
7. [What happens next](#7-what-happens-next)

## 1. The question and the method

A fit can reproduce its target yields through the wrong paths: a mol % of 2-methyl-3-furanthiol (MFT) at
20 min says nothing about whether the cysteine that made it lived for 20 minutes or 2 seconds. So the
question is not "does the model hit the number" but "which reactions carry the flux, and is that what the
pot did".

For a pot, the flux through reaction *j* is the integral of its rate over the run,
F<sub>j</sub> = ∫ r<sub>j</sub>(t) dt (mmol/L). Each species' change is then a stoichiometric sum,
c<sub>i</sub>(end) − c<sub>i</sub>(0) = Σ<sub>j</sub> ν<sub>ij</sub> F<sub>j</sub>, which is a check the
budget must pass: if it does not close, the fluxes were computed with rates the trajectory did not
see. Ranking the F<sub>j</sub> that make or consume a species gives its sources and sinks.

Two technical requirements turned out to matter, and the first version of the tool violated both:

- **The output grid must resolve the fastest transient.** Rates are integrated by the trapezoid rule on
  the output grid; a transient lasting seconds needs points at milliseconds. The tool uses the union of
  a linear grid and a log-spaced one starting at 10⁻⁷ × the run time.
- **The rates must be the ones the solver used.** The rate constants depend on the in-situ pH, which the
  integrator holds at a set of time nodes (sec. 3). The run now keeps those node constants
  (`SulfurRun.rate_nodes`) and one function interpolates them for both the ODE and the audit
  (`sulfur.rate_vector_at`), so the two cannot disagree.

With both, every budget closes to better than 0.03 %.

## 2. The paths, as the model runs them

The lane's carbon trunk, simplified to what carries flux in a cysteine + pentose pot:

```
pentose ──(cysteine-catalysed 2,3-enolisation, r_pent_dpo)──▶ 1-deoxypentosone (DPO)
pentose ──(acid-catalysed 1,2-enolisation, r_pent_tdp)──────▶ 3-deoxypentosone (TDP) ──▶ furfural ──+H2S──▶ FFT
DPO ──(Strecker with cysteine, r_dpo_ddp; consumes Cys)──▶ 1,4-dideoxypentosone (DDP) + NH3 + CO2 + [S → inert]
DPO ──(r_dpo_nf)──▶ norfuraneol (NF) ──+H2S──▶ MFT
cysteine ──(thermolysis, r_cys_h2s, Zheng & Ho measured)──▶ H2S
```

What the flux map shows on the B60 calibration (`results/validation/sulfur_flux_attribution.md`; on B51's vector the picture is the same to a few per cent):

| pot | cysteine half-life | pentose half-life | main cysteine sink | main pentose sink |
| --- | --- | --- | --- | --- |
| Hofmann ribose pH 5, 145 °C, 20 min | 2.4 s | 2.2 s | Strecker `r_dpo_ddp` 99 % | `r_pent_tdp` 60 %, `r_pent_dpo` 40 % |
| Zhang 2025b xylose pH 6.5, 120 °C, 50 min | 0.5 s | 0.3 s | Strecker 100 % | `r_pent_dpo` 100 % |
| Cao 2017 xylose 0.4 M, 90 °C, 5 h | > 5 h (33 % lost) | 1.8 s | Strecker 54 %, thermal 45 % | `r_pent_tdp` 81 % |
| Zhai 2019 xylose 2:1, unbuffered, 90 °C | > 160 min (49 % lost) | 9 s | Strecker 88 % | `r_pent_tdp` 74 % |

Three things stand out. All of the sugar chemistry is over within seconds; afterwards the pot sits
nearly static while slow cysteine thermolysis trickles H2S into a pool of norfuraneol and furfural. MFT
comes 98–100 % from norfuraneol + H2S; the DDP route (the "intact-skeleton" precursor in Cerny 2003's
labelling) carries none, because its fitted constant sits at 3 × 10⁻⁷. And the in-situ pH collapses on
the same seconds timescale: Cao's pot goes 5.28 → 3.91 within 18 s (0.1 % of the run), against a measured 5.49 → 4.67 at
30 min and 4.15 at 5 h.

## 3. Finding 1: the pH trajectory was not converged (B60, fixed)

The pH is not an input; it is solved from a charge balance against the state (a conserved strong-ion
difference against every titratable group). Because the rates depend on pH and the acid that moves the pH
comes from the rates, it is a fixed point, solved by sweeps: integrate with a guessed pH trajectory,
recompute the pH from the states, repeat. Two numerical defects sat in this loop.

**The nodes.** The pH was held at 41 times equally spaced over the run. If the pH collapses within 0.1 %
of the run, the first node after the start sits at 2.5 % of it, so through the entire transient every
pH-dependent constant was evaluated at the starting pH. Doubling the resolution moved Hofmann's MFT 39x;
the fix puts 41 more nodes on a geometric grid from 10⁻⁵ × the run time (verification: 41 vs 81 nodes
agree within 0.001 dex on all 52 fit rows; the old grid was off by up to 0.83 dex on 17 of them).

**The sweeps.** The iteration was plain Picard, stopped after a fixed number of sweeps (2 in the fitting
loop, 4 elsewhere) without checking the residual. In unbuffered pots it oscillates: the Kang TTCA pot
alternates final pH 5.52 / 6.55 / 5.56 / 6.37 on sweeps 2–5 — high pH speeds the base-catalysed
enolisation, which makes acid, which lowers the pH, which slows it — and settles only after ~20–40 sweeps.
The fix is under-relaxation (move the guess 70 % of the way to the new value) plus iterating to a 10⁻³ pH
residual; it converges in ≤ 10 sweeps on every fit system, and the converged answer is the same for every
relaxation factor tried, so the damping picks no branch.

Every fit from B2.2 to B58 optimised its constants against these errors, so they absorbed them. B60
refits B51's objective on the converged numerics (kind 1 under the acceptance rule: shipped on
verification, `results/validation/kinetic_core_b60_verification.md`; verdict in
`kinetic_core_b60_prereg.md` section 6). The verdict carries the diagnosis. B51's reported cost of 20.4
was flattered by the numerics: the same vector integrated correctly costs 33.2, and the refit brings that
down only to 27.1. The rows that got worse are the cysteine + pentose thiol rows at 100–120 °C (Kang,
Zhai's temperature fold; up to +2.5 σ over-predicted), which the oscillating pH had held down. No
constant in the free set can slow the trunk, so the fit cannot fix them: they are finding 2.

## 4. Finding 2: the trunk runs in seconds; the pots say hours (B59, open)

Three papers in the corpus measure what happens to the cysteine and the sugar in exactly these pots:

| pot | measured | model |
| --- | --- | --- |
| Cao 2017, Cys + Xyl 0.4 M, 90 °C, pH 5.5 → 4.15 | free Cys + TTCA + Cys-Amadori = 0.425 M at 5 h: **no measurable cysteine destroyed**; xylose likewise; ~60 % of both held in TTCA | 33 % of the cysteine and **all** of the xylose destroyed, the xylose in ~2 s |
| Zhang 2025b, 40 + 40 mM, 120 °C, pH 6.5 | 20.8 mM of cysteine still accounted for at 50 min: **~48 % destroyed in 50 min**, xylose 1:1 with it | 99 % destroyed, in ~1 s |
| Zhai 2019, Xyl 165 + Cys 83 mM, 90 °C | TTCA holds 65–77 % of the cysteine at 160 min | TTCA ≈ 0 (B58 showed that even with TTCA formation switched on, the sugar is gone before it can bind) |

So at 90–120 °C the model destroys cysteine and sugar roughly 10³–10⁴ times faster than the pots do.

**Why the fit allowed it.** Every sulfur fit row is an endpoint: a yield of a thiol at the end of a run,
almost all at 145 °C. An endpoint fixes the *ratios* in which material splits between branches, not how
fast it moves: a trunk running in 2 seconds and one running in 2 hours give the same yields if the
branch ratios are the same. The optimiser therefore had no reason to keep the trunk slow, and the
temperature dependence that would have slowed it at 90 °C is one shared activation energy (64 kJ/mol)
fitted to that same endpoint data. The two measured points above imply an apparent Ea of at least
~150 kJ/mol (confounded with pH), in line with the 100–150 kJ/mol typical of sugar degradation in
Maillard systems.

**What the chemistry says the path is.** de Roos (1992, 2005) and Cao 2017 describe it: cysteine and
a pentose condense quickly and reversibly to the thiazolidine (TTCA; Cao's equilibrium constant
≈ 60 L/mol), which is a reservoir, not a dead end; the open Schiff base rearranges slowly to the
cysteine Amadori compound; and the Amadori compound breaks down to the deoxyosones, **releasing the
cysteine** (an Amadori 2,3-enolisation frees the amino acid). The cysteine is then lost only through
slower channels: Strecker degradation with the dicarbonyls, thermolysis. The model instead has the
cysteine act as an instant catalyst on the free sugar (`r_pent_dpo`, `r_pent_tdp`) and then be destroyed
by the dicarbonyl it made — a shortcut that only works because nothing constrains its speed.

**Consequences.** Every prediction that depends on how long cysteine survives inherits the error: H2S
release (thermolysis needs free cysteine), the thiols, the pH drift, and the temperature extrapolation
(a too-low Ea under-predicts how much slower everything gets below 145 °C). The B58 result fits this:
switching on TTCA formation, which shelters even a trace of cysteine for the first seconds, raised the
Hofmann pH 5 thiols 10–25x.

## 5. Finding 3: the pH sweep and the thiol branch ratio (B60 calibration, open)

On converged numerics the scorecard moved a lot. The Hofmann ribose and xylose pH 5 hold-outs went from
37–118x off to 1.5–3.6x, so B51's "~100x low, the thiolate sink" was a numerical artefact. The panel went
from 16/50 to 19/50 within 3x. What remains in the cysteine + pentose rows is three path-level errors that
the flux map names directly. The Hofmann pH 3 and pH 7 arms are now pots in the attribution generator.

| Hofmann ribose + Cys, 145 °C, 0.5 M phosphate | pH 3 | pH 5 | pH 7 |
| --- | --- | --- | --- |
| measured MFT / FFT, ppb | 553 / 229 | 198 / 121 | 25 / 12 |
| model MFT / FFT, ppb | 25 / 1 523 | 94 / 347 | 0.005 / 0.17 |
| pentose to DPO (1,2- vs 2,3-enolisation) | 0.4 % | 40 % | 89 % |
| H2S made over the run, mM | 1.43 | 0.13 | 0.0002 |
| cysteine destroyed by Strecker | 3 % of its loss | 98 % | 100 % |

- **pH 7: H2S starvation (4 900x low).** The base-catalysed `r_pent_dpo` makes 49 mM of DPO in
  seconds. Its Strecker step `r_dpo_ddp` consumes all 33 mM of cysteine, leaving none for thermolysis,
  and 13.5 mM of norfuraneol sits unused for want of H2S. The Strecker step sends cysteine's sulfur
  (mercaptoacetaldehyde) to an inert pool. In the literature that sulfur is the classic Maillard H2S
  source: mercaptoacetaldehyde → H2S + acetaldehyde, and Du 2023's labelling supports it. Hofmann's
  smooth ~20x fall in MFT from pH 3 to pH 7 needs H2S at pH 7, and the only cysteine left to supply it
  went through Strecker. The split is the "smaller finding" of section 6, but it is NOT small: it is the
  pH 7 arm.
- **pH 3: the MFT branch is shut (MFT 22x low, FFT 6.6x high).** The acid factor sends 100 % of the
  pentose down `r_pent_tdp` to the 3-deoxyosone. DPO, and with it norfuraneol and MFT, gets 0.4 %.
  Hofmann's MFT is at its HIGHEST at pH 3. So the model's 1,2- vs 2,3-enolisation pH law is too steep, or
  MFT has a low-pH source the model lacks. Two candidates, both unmeasured in the corpus as rate laws:
  phosphate catalysis of the enolisation (0.5 M phosphate in every Hofmann arm), and the 2,3-enolisation
  of the Cys-Amadori compound, which frees cysteine and is the de Roos path B59 adds.
- **FFT over MFT is 6–60x too high (Hofmann pH 5, Yiltirak).** Measured FFT/MFT is 0.6 at Hofmann pH 5
  and 0.2–0.9 across Yiltirak's 100–130 °C pots. The model gives 3.7 and 13–20, and the excess grows as
  the temperature falls (Yiltirak FFT 234x high at 100 °C / 4 h, 59x at 130 °C / 30 min; MFT within
  1.5–7x). FFT is 95–99 % furfural + H2S, MFT 100 % norfuraneol + H2S, and both draw on the same H2S. So
  the mis-weighting sits in the furfural branch's temperature dependence (`r_tdp_fur`, `r_fur_fft_hs`)
  against the norfuraneol branch, or in the FFT-specific sink. It is not a pH effect.

All three are statements about paths, not constants, and none can be fixed by refitting B60's free set.
The first two are B59's subject: the Strecker sulfur split, and the enolisation pH law through the
Cys-Amadori route. Hofmann's pH sweep is the natural pre-registered hold-out for both, since its shape
(MFT falling ~20x and FFT ~19x from pH 3 to 7) is what a correct mechanism must reproduce. The third
needs a measured furfural + H2S → FFT barrier, or FFT time courses at two temperatures; it is on the
download list.

## 6. Smaller findings

- **The Strecker sulfur is discarded.** `r_dpo_ddp` sends the cysteine's sulfur (as mercaptoacetaldehyde)
  to an inert pool, deliberately ("rather than being silently promoted to H2S"). The textbook route is
  mercaptoacetaldehyde → H2S + acetaldehyde, and Du 2023's labelling supports it (84 % of a thiophene's
  C2 unit is cysteine-derived acetaldehyde), but the split between H2S and ring products (thiophenes,
  trithiolanes) is unmeasured. It decides the Hofmann pH 7 arm (section 5).
- **One flux tool did not close.** `sulfur_flux_budget`, which scores the Cerny thiamine-share fit row,
  integrated on a coarse linear grid with constants that ignored B51's ionic strength; in Cao's pot it
  booked 2 575 mM of 3-deoxypentosone from 400 mM of xylose. Fixed with B60.
- **The DDP route is off.** All MFT comes via norfuraneol. Cerny's labelling (intact C5 skeleton) cannot
  tell the two apart, since norfuraneol is also intact C5, so this is not a contradiction — but the
  route REACTION_TREES.md calls "the intact-skeleton route" carries nothing.

## 7. Why nobody saw it, and the general lesson

The scores looked like ordinary model error: thiol predictions 10–100x off at some conditions, with a
plausible mechanism to blame (B51 blamed the thiolate-acting thiol sink). The flux map shows the
thiolate sink takes 0–4 % of MFT's loss in the pH ≤ 5.5 pots (28 % of MFT's and 78 % of FFT's in Zhang's pH 6.5 pot); what governs them is how long free cysteine
survives the first seconds. The lesson for the project: **check paths and timescales, not only yields**.
A model fitted to endpoints should be audited for the speed of its intermediates against any time course
in the corpus, and its numerics should be shown converged at the speed of its fastest transient — which
is why the flux attribution is now a generator with a closure check, and why B60 records its convergence
artifacts.

## 8. What happens next

1. **B60** (shipped 2026-10-28): refit on the converged numerics; scorecards and envelope regenerated.
2. **B59** (planned, pre-registered before it runs): the de Roos path — PENT + Cys ⇌ TTCA, Schiff base →
   Cys-Amadori → deoxyosones with cysteine released — in place of the direct cysteine-catalysed
   enolisation; Cao 2017's time courses as fit rows (free Cys, free Xyl, TTCA, Cys-Amadori, pH at six
   times); Zhang 2025b, Zhai 2019 and de Roos's purified Cys-Amadori decay as hold-outs; the trunk's Ea
   released from the shared value; the Strecker sulfur split (section 5); Hofmann's pH 3 / 5 / 7 sweep
   as a shape hold-out. Backlog: `tasks/data_restructure_plan.md`, "B59".
3. Papers that would sharpen B59: Edwards & Wedzicha 1992 (cysteine + dicarbonyl kinetics), de Roos 1992
   and Mizutani 1960 (thiazolidine pH dependence, H2S release), Hofmann & Schieberle 1995 (the cysteamine
   route) — see `LITERATURE_SEARCH_2026-10_ROUND3.md`.
