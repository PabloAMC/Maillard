# Wave B11 pre-registration — oxygen as an input (programme step R2(a))

> **Renumbered B10 -> B11 on 2026-09-06 (owner: "take the most sensible decision long term").** The
> temperature structure (R2(c)) runs first as wave B10 (`kinetic_core_b10_prereg.md`): the probe below
> showed the oxidant channels inert at trace thiol, the corpus holds four unread thiol temperature
> ladders, and every prediction away from 145 C depends on the slopes. B10 also charges the
> ambient oxidant consistently (finding 2.1), so B11 starts from a consistent baseline. Nothing else
> in this draft changed.

*Written 2026-09-06 before any B10 number existed. STATUS: DRAFT FOR OWNER SIGN-OFF. No B10
generator exists yet; nothing below has been fitted. Programme: `tasks/data_restructure_plan.md`,
"Reaction-modelling programme". Licence for the data it reads: FIT_HOLDOUT_DECLARATION.md
Amendment 19 (the vessel block).*

## 1. What this wave is for

Step R1 put the physical state of every panel pot on record: fill, vessel, atmosphere, water. The
scorecard now prints the oxygen each closed vessel held per mol of thiol charged: Yiltirak 2026
1.98, Bolton 1994 2.07, Hofmann 1998 Table-1 systems 0.27, Hofmann fed-intermediate systems 1.32.
The core is near-unbiased on Hofmann's Table-1 pots and over-predicts the thiols 20x (Bolton) and
9-480x (Yiltirak). B11 asks whether making oxygen an INPUT the network consumes explains the
oxygen-rich misses without breaking the oxygen-poor fits.

## 2. Findings that shape the design (measured 2026-09-06 on the shipped B9 lane, no fit)

**2.1 A fit/deploy inconsistency.** Every fit system in the B2.3-B9 objective was integrated with
the oxidant pool `OX` charged at `OX_AMBIENT_MMOL_L = 1.0` (Zhang 2024's cystine arm at 62.4).
The engine's `predict` passes only the mapped precursors as the initial state, so every panel
prediction and every user prediction runs at `OX = 0`. The two oxidant-consuming channels
(`ch_dimer_mft`, `ch_dimer_fft`) therefore carry zero flux in deployment and non-zero flux in the
fit, and both `k_dimer_*` sat on active bounds in B9. B10 removes the inconsistency (it is a one-line consistency fix); B11 builds on it.

**2.2 The probe.** The shipped lane was integrated with `OX` = 0, 1, the pot's own O2 in mmol per
litre of liquid, and twice that (`scratch/ox_probe_2026-09-06.py`, gitignored; the numbers are
reproduced here). MFT and FFT in ug/L; the dimers are reported as the fraction of the thiol they
remove.

| pot | O2 in pot (mmol/L) | measured MFT / FFT | OX = 0 (engine today) | OX = 1 (the fit's ambient) | OX = O2 | OX = 2 x O2 |
|---|---:|---|---|---|---|---|
| Hofmann ribose + cys, 145 C, pH 5 | 9.0 | 198 / 121 | 731 / 834 | 731 / 831 | 727 / 803 | 724 / 775 |
| Yiltirak buffer, 100 C / 4 h | 49.6 | 6.88 / 1.28 | 64 / 616 | 64 / 610 | 63 / 447 | 61 / 367 |
| Yiltirak buffer, 130 C / 0.5 h | 49.6 | 1.71 / 1.62 | 171 / 214 | 171 / 212 | 155 / 153 | 143 / 126 |
| Bolton thiamine + cys + glc, 120 C | 24.2 | 11.7 / - | 236 / 534 | 236 / 533 | 235 / 523 | 234 / 514 |
| Hofmann norfuraneol + cys (fed, mM) | 26.4 | 1020 / - | 2360 / 0 | 2353 / 0 | 2201 / 0 | 2077 / 0 |

Reading, in order of consequence:

1. **At trace thiol levels the shipped oxidant channels cannot matter.** Both are SECOND order in
   thiol. At 1-700 ug/L (1e-8 to 6e-6 M) they remove under 10 % of MFT and at most 1.7x of FFT even
   with fifty times the ambient charge, and the pool is not consumed (`OX` left = `OX` charged).
   Charging the pot's oxygen into the network as it stands changes nothing the panel scores.
2. **The only first-order thiol sink, `ch_thiolate_loss_*`, carries no oxygen dependence.** It was
   fitted on Kumazawa 2003's pH grid, which the dossier confirms was canned "without the
   deoxidization process" (air in the can, volume unstated) at 1 ppm FFT, i.e. oxygen in large
   excess. The same constant runs unchanged in an oxygen-poor pot. If oxygen is the hidden input,
   THIS channel is where it has to enter: first order in thiol, first order in dissolved oxygen.
3. **The Yiltirak 100 C / 4 h FFT miss (480x) is not an oxygen miss.** At the same pot MFT is 9x over.
   The split is the frozen 64 kJ/mol formation barrier applied to the furfural + H2S route at 100 C
   over four hours; that is step R2(c)'s problem (temperature structure) and this wave does not
   claim it. The 130 C / 0.5 h rung (100x / 132x, both thiols alike) is the rung this wave can
   speak to.
4. **The fed-intermediate rows are where the dimer channel is real.** At mM thiol the second-order
   channel removes a material share (norfuraneol + cys: ~5 % of MFT-equivalents at OX = 1, ~50 % at
   the pot's own O2). Those rows are in the objective; they will re-fit `k_dimer_*` once `OX` is
   charged consistently, and the wave must show they do not degrade.

## 3. The structural change

Three additions to the sulfur lane, each declared, none fitted to a hold-out.

**(a) Oxygen as a two-pool state.** `OXR` = headspace reservoir in mmol per litre of liquid,
charged from the vessel block (`vessel.py`: headspace air at 1 atm / 20 C, plus dissolved O2);
`OX` = dissolved O2, held at its saturation value at temperature while `OXR > 0` and drained by
every O2-consuming flux; when the reservoir is exhausted `OX` falls with the consumption. The
gas-liquid equilibration is taken as fast on the cook's timescale (stirred 3 mL tubes; a 100 mL
autoclave charge is the weakest case and is said so). Saturation: Henry's law for O2 at the segment
temperature under the sealed tube's air partial pressure (the fixed air charge at the higher
temperature, ~0.29 atm O2 at 130 C), giving ~0.2-0.4 mmol/L; carried as a DECLARED constant with a
band (0.1, 1.0) mmol/L sampled by the envelope, not fitted.

**(b) Oxygen consumers.** (i) `ch_thiolate_loss_mft/fft` become first order in `OX`, normalised to
`[O2]_ref` = the saturation value at Kumazawa's 121 C, so B9's `k_thiolate_loss` keeps its meaning
at Kumazawa's condition and the fit only has to move it if the other rows ask. `ch_dimer_*` read the
same dissolved `OX`. (ii) NEW `ch_cys_ox`: cysteine + 1/2 O2 -> 1/2 cystine-equivalent (to
`FRAG_S`), constant `k_cys_ox`, declared band (1e-4, 1e-1) /min at 145 C (thiol autoxidation is
metal-catalysed; Bagiyan 2004 gives initial rates, not constants, so the band is wide and the
coordinate is expected to be unidentified). (iii) NEW `ch_red_ox`: the reductone pool (Amadori +
deoxyosones, the species the network already carries) + O2 -> oxidised fragments, constant
`k_red_ox`, declared band (1e-4, 1e-1) /min at 145 C, likewise expected unidentified. Both sinks
exist so that an oxygen-POOR pot can run out of oxygen; without them a 9 mmol/L charge would never
deplete and the wave could not distinguish the regimes. A trace-metal multiplier for tap water is
NOT added: no measurement supports a number; the water source stays a recorded caveat.

**(c) The vessel on every fit system.** The fit's systems live in the frozen generators as
literals; B10's generator carries its own vessel table with provenance per system. Known today:
Hofmann 1998 (200 mL autoclave, 100 mL Table-1 pots / 50 mL fed pots); Kumazawa 2003 (air in the
can, volume unstated: declared `oxygen in excess`, which at 1 ppm FFT no can could contradict).
To be read before the wave runs, all PDFs on disk: Kang 2026 (sealed pressure vessels, volume not
in the dossier), Zhou 2023, Zhang 2024, Feng 2022, Zhai 2023 (pressure bottles), Whitfield 1999,
Cerny 2007, van Seeventer 2001, Yaghmur 2005. Rule for a system whose source leaves the vessel
unstated: `OXR` set to Hofmann's 9 mmol/L with a (1, 100) mmol/L band sampled by the envelope, and
the system is listed as NOT identifying `k_cys_ox` / `k_red_ox`.

## 4. What is fitted, what is held

| | B10 (temperature wave, the baseline) | B11 |
| --- | --- | --- |
| objective rows | 60 | **60, unchanged** (no level row enters; the vessel is an input, not a target) |
| network | B9 | + `OXR` state, + `ch_cys_ox`, + `ch_red_ox`; `ch_thiolate_loss_*` and `ch_dimer_*` read dissolved O2 |
| free set | 23 | **25** (+ `k_cys_ox`, `k_red_ox`); `k_thiolate_loss`, `k_dimer_mft`, `k_dimer_fft` re-fit in the new structure |
| declared bands | B9 | unchanged + the two new bands above; `[O2]_sat` and unstated `OXR` are envelope bands, not fit coordinates |
| pH weighting, optimiser, budget, starts | B9 | unchanged; start 0 = B9's optimum with the new constants at their band centres, start 1 = the perturbation protocol |
| every other coordinate | B9 | unchanged, including the single formation Ea (R2(c) is a separate wave) |

## 5. Pre-registered tests and falsifiers

The wave's own claim is the OXYGEN-REGIME CONTRAST. Scored on the panel after the fit is frozen,
with today's numbers as the baseline.

- **T1 Bolton 1994** (O2 : thiol 2.07; today 20.2x over): fold error falls below **6x**.
  Falsifier: stays above 12x.
- **T2 Hofmann Table-1 pH-5 rows** (O2 : thiol 0.27; today MFT 3.7x, FFT 6.9x over): neither
  worsens beyond **8x**. Falsifier: either exceeds 8x.
- **T3 Yiltirak 130 C / 0.5 h** (O2 in excess; today MFT 100x, FFT 132x over): both fold errors fall
  by at least **3x**. Falsifier: neither improves by 2x.
- **T4 Yiltirak 100 C / 4 h MFT-vs-FFT split** (9x vs 480x): NOT a target of this wave; recorded so
  the wave is not credited or blamed for it. Expected: the split persists until R2(c).
- **T5 in-sample discipline:** no B9 objective row's |residual| grows by more than 0.3 dex; the
  Kumazawa grid stays inside its sigma. Falsifier: any row moves more than 0.5 dex.
- **T6 fed-intermediate rows** (Hofmann T3/T4/T10 mol% rows, mM thiol, dimer channel live): every one
  stays within 2x of its B9 residual.
- **Secondary, not gating:** directional panel not below 17/26; envelope coverage not below 5/33;
  Laplace identifies at least 20 of the 25 coordinates; `k_cys_ox` and `k_red_ox` are EXPECTED to be
  unidentified and are then drawn across their bands by the envelope (the 2026-09-04 rule).

**Ship rule.** B11 ships as the engine's sulfur report if T1, T3 and T5 hold. If T1 AND T3 both fail,
the oxygen hypothesis as structured here is refuted for the sulfur lane: the vessel block stays (it
is a fact about the pots), the two new constants are removed, the fit/deploy `OX` inconsistency is
fixed by charging the ambient value consistently, and R2 proceeds to (c). Partial outcomes (one of
T1/T3) ship with the finding recorded and the failing test named in the model card.

## 6. Forecasts, revised after the probe

P(T1 passes) 0.35. P(T3 passes) 0.25. P(B11 ships under the rule) 0.30. P(any Yiltirak rung lands
within 3x after B11 alone) 0.05, down from the 0.30 written in the programme before the probe: the
probe showed the second-order channels are inert at trace thiol and the 100 C FFT miss is a
temperature-structure miss.

## 7. Cost and order

1. Vessel reads for the nine fit-system papers above (about half a day; each becomes a dossier line
   with a verbatim quote).
2. Engine and network: `VesselSpec` on `ProcessSpec`; the two-pool O2 state; three reaction edits;
   bands; charge-closure ledger entries; unit tests that the shipped B9 numbers are reproduced
   exactly when the vessel is absent and `OX` is charged at the ambient value (one day).
3. `generate_kinetic_core_b11_fit.py` (B9's shape: import, extend, freeze), two starts (about an
   hour of compute), consolidation.
4. Laplace, profile, fit targets, scorecard, directional, envelope, model card, README re-pin
   (half a day).

## 8. What is read and what is not

As B9: the fit reads its row table and its own vessel table (literals with provenance) and nothing
under `data/benchmarks/`. The panel bundles' vessel blocks are read only by the scorer, after the
fit is frozen. The hold-out guard runs on the B11 generator like any other.
