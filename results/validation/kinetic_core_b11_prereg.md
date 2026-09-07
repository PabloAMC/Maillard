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

## 9. Amendments before the run (2026-09-07, written after B10 re-merged and before any B11 number)

1. **Base wave.** B10 did not ship its route split, so B11 builds on B9: 54 objective rows, B9's
   23 free coordinates plus the two consumers (25 free). The six Yiltirak folds are NOT in the
   objective (they were B10's and stay with its record).
2. **Units.** Dissolved oxygen keeps the ambient unit every fit system was integrated at
   (`OX = 1.0` = the air-saturated liquid, `OX_SAT_MMOL_L` = 0.3 mmol/L declared, band 0.1-1.0).
   The reservoir `OXR` is in the same units per litre of liquid (Yiltirak 165, Hofmann's 100 mL
   pot 30, Bolton 81). The two-pool state is mass-action: a dissolved-oxygen VACANCY `OXV` is
   created by every consumer and refilled from the reservoir at a declared fast rate, so `OX`
   never exceeds saturation and falls only when the reservoir is spent. With the consumers at
   zero every wave before B11 reproduces bit for bit; B9's `k_thiolate_loss`, `k_dimer_*`
   therefore keep their meaning and are re-fitted only through the depletion they now see.
3. **The vessel table of the fit systems.** Read from the PDFs on disk on 2026-09-07: only
   Hofmann 1998 states its volumes (200 mL autoclave; 100 mL Table-1 pots, 50 mL fed pots ->
   30 and 88 units). Cerny 2007: 1.00 mL in Teflon vials of unstated volume. Whitfield 1999:
   flame-sealed 5 mL ampoules, fill unstated. Kang 2026, Zhai 2023, Feng 2022: pressure-rated
   glass vessels or bottles, volumes unstated. Zhou 2023, Zhang 2024: glass vials, unstated.
   Kumazawa 2003: a can without deoxidisation, volume unstated. van Seeventer 2001: a 2 L
   autoclave then closed bottles under air. Yaghmur 2005: a 100 mL vial. Every unstated system
   is charged with the declared default (30 units) and marked as identifying neither consumer.
   **Consequence, stated in advance:** the objective contains ONE laboratory's vessel and no
   oxygen contrast; T1 (identification of the consumers) is EXPECTED to fail, and the wave's
   value is the declared structure with its bands in the envelope plus the out-of-sample tests.
4. **Ship rule unchanged** (sec. 5). If T1 fails and T3 does not improve, the structure ships
   as declared-inert (consumers zero) with the priors marked "not sampled: B11 not shipped",
   which is what the engine carries today; the reservoir arithmetic and the vessel plumbing
   stay because they are facts about the pots.

## 10. OUTCOME (2026-09-07, after the run) — DO NOT SHIP: the consumers ship as declared-inert

Two starts, 600 evaluations each (budget-exhausted, status −9), `kinetic_core_b11_fit_report.json`,
`kinetic_core_b11_laplace_covariance.json`, `kinetic_core_b11_ship_rule.md`.

- **The fit walked both consumers toward their floors:** k_cys_ox 3.7e-4, k_red_ox 1.3e-4 per unit per
  minute (log10 −3.43 / −3.90 from a −3.0 start); cost 19.71 against B9's 18.74 on the same 54 rows.
  The Laplace gives sigma 8.1 and 19.9 dex on the two consumers — **unidentified, as sec. 9.3 said**; only
  15 of 25 coordinates identified (B9: 20 of 23), the budget-exhausted optimum being less curved.
- **T1 Bolton 1994: 20.2x → 20.2x. FAIL.** With the consumers this small, no reservoir (Bolton's 81
  units included) is drained on Bolton's hour at 120 C.
- **T3 Yiltirak 130 C: MFT 99.7x → 102.5x, FFT 130.7x → 129.3x. FAIL** (neither improves 2x: falsified).
- **T5 in-sample: worst +0.20 dex** (`zhang_fig1_gcys_dimer_over_MFT`) — passes; the dimer channels now
  vacate dissolved oxygen, and the fit paid for it with k_dimer_decay +1.85 dex and k_thiol_decay +0.6.
- **T6 fed rows: PASS** (every one within 2x of its B9 residual).
- **T2 Hofmann Table-1 (not gating):** the two hexose MFT rows predict zero under B9 and B11 alike (the B9
  finding); the answered rows do not move.

**Ruling under sec. 5:** T1 AND T3 both fail, so the oxygen hypothesis AS STRUCTURED HERE — first-order
consumers draining a headspace reservoir — is refuted for the sulfur lane on this corpus. Per sec. 9.4 the
two consumers ship at ZERO (declared-inert), the two-pool state, the vessel plumbing and the reservoir
arithmetic stay (they are facts about the pots), the engine keeps reading B9, and the envelope priors for
`sulfur.oxygen.*` stay marked "not a free coordinate of the shipped report". The B11 report, members,
Laplace and ship rule are kept as the record.

**What the negative result says.** One laboratory's vessel in the objective cannot identify an oxygen
consumer; and the between-lab gaps (Bolton 20x, Yiltirak 100x) are not closed by any first-order
consumption within the declared bands — the bands would have to reach 10^-1 and above, where the pots'
cysteine is destroyed (probe: at 50 per unit per minute cysteine falls 33 → 0.04 mM) and every fed row
breaks. R2 proceeds to (c): the paper-level response factor (R3) is now the leading explanation for the
between-lab level gap, and the wet-lab oxygen axis (R7) is the only measurement that can settle oxygen.
