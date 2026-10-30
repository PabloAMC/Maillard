# Sub-model 3: the deoxyosones and their fates (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Do the pentose-derived dicarbonyls (1-deoxypentosone, 3-deoxypentosone, 1,4-dideoxypentosone),
norfuraneol and furfural form and disappear at the measured rates, given the route's supply?

**Decision.** Whether the unpinned branch constants (`k_dpo_ddp`, `k_dpo_nf`, `k_nf_decay`; B2.3) can be pinned on
literature data, and whether 1,4-dideoxypentosone, which no laboratory has measured, can stay a modelled species or
must be declared an assumption with a band.

**Rival arms.** Norfuraneol from 1-deoxypentosone (current) vs also from the Cys-Amadori compound directly; furfural
from 3-deoxypentosone only vs also from the pentose's acid-catalysed dehydration.

**Expected result.** Norfuraneol and furfural time courses pin their formation and loss (identified, ~65 %);
1,4-dideoxypentosone stays unpinned (~80 %), and its MFT route is reported as an assumption.

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 2 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2) -- first reading (2026-10-29); adjudicated table below

First reader: agent-sm3-reader1-2026-10-29, against the PDFs. Per-check verdicts and quotes:
`data/lit/source_verification/deoxyosone_submodel.yml`. Rows: `data/lit/fit_rows/deoxyosone_submodel.yml` (status:
draft; 24 pots, 315 rows). Quantification is the class the paper supports; the ledger's class, which the loader reads,
is given in brackets where they differ.

| ledger id | system | what (rows) | quantification | laboratory | verdict |
|---|---|---|---|---|---|
| zhai2019_fig3b_dicarbonyls_ttca | TTCA 5 mM +- Cys 5 mM, water, pH 7.4 initial, 120 C, 10-60 min | 3-DX, 1-DX, GO (36) | semi-quantitative: hexanedione IS, no dicarbonyl standard [ledger: calibrated] | Jiangnan (Zhang X) | corrected (quantification); transcribed |
| zhai2021_fig5_dicarbonyls_ph_ladder | TTCA 10 mM +- xylose 10 mM, water, pH 5.5 / 6 / 7 / 8 initial, 100 C, 20-120 min | 3-DX, 1-DX, MGO (144); GO left out | 3-DX / 1-DX semi-quantitative (diacetyl IS) [calibrated]; MGO calibrated (10 rows below LOQ, flagged) | Jiangnan | corrected; GO fails (every value below the paper's LOD) |
| zhai2023foodchem_fig2_deoxyosones | TTCA 10 mM +- xylose 10 mM, water, pH 7 initial, 120 / 140 C, 20-120 min | 3-DX, 1-DX (44) | semi-quantitative [calibrated]; the paper says "external standards", Zhai 2021 says diacetyl IS | Jiangnan | corrected; 100 C arms are Zhai 2021 Fig. 5 pH 7 re-plotted (not transcribed) |
| zhai2023foodchem_fig3_cys_da_mgo_go | as above | MGO 120 / 140 C (24); GO, Cys, diacetyl left out | MGO calibrated | Jiangnan | corrected; 100 C re-plot; GO below LOD |
| hofmann1998_t5_fa_nf_pairs | ribose 100 + Cys 33 mM, 0.5 M phosphate pH 5, 145 C, 20 min (endpoint) | furfural, NF (2) | calibrated: internal standard / labelled homologue with measured response factors [unknown] | DFA Garching | corrected (quantification); furfural row flagged CF-016 |
| mikami2017_fig6b_norfuraneol_decay_ph8 | NF 1.7-1.8 mM alone, 0.2 M phosphate pH 8, 100 C, 0.5-6 h | NF decay, MGO (5) | NF calibrated (pure compound held); MGO quinoxaline vs authentic quinoxaline, OPD equimolar (flagged) | Ochanomizu (Murata) | pass; NF at 3 / 6 h plotted zero (left out) |
| zhou2025b_figs1_3_axarp_xylose_120c | Ala-xylose Amadori 20 mM +- xylose 20 mM, water, pH 8 initial (measured after), 120 C, 30-120 min | 3-DX, 1-DX, GO, MGO, furfural (40) | furfural calibrated (curve printed); 3-DX / 1-DX semi-quantitative (hexanedione IS) [calibrated]; GO / MGO calibrated | Jiangnan (Zhang X / Ho) | corrected; NOT a re-plot of Zhou 2024 but a 7-24x conflict with it (rows flagged) |
| zhou2025b_fig2a_xylose_alone_dicarbonyls | xylose 20 mM alone (inferred, checked against the printed 1.1 %), pH 8, 120 C | 3-DX, GO, MGO (12) | as above | Jiangnan | corrected (conditions inferred) |
| desclaux2006_fig2_xylose_glycine_100c | xylose + glycine 240 mmol/kg each, 0.5 M phosphate pH 6, 100 C, 30-120 min | NF, 3-DX (8), as ug/kg from "umol per mmol initial sugar" | semi-quantitative: GC / GC-MS, no calibration stated [unknown] | Reading (Mottram) | corrected (units, class); transcribed |
| cao2017 (fig4ab, fig4c, table1) | Cys + xylose 0.4 M, 0.2 M phosphate pH 5.5, 90 C | no sub-model 3 species over time (furfural only at 3 h, IS area ratio) | -- | BTBU (Xie) | no data for this sub-model |
| hou2017 (uv/pH, table2, fig3/4) | Gly-Amadori +- Cys, 0.2 M phosphate pH 5.5, 140 C | A294 / A420 / pH; furfural as peak area at 30 min only | -- | BTBU (Xie) | no data for this sub-model |
| deng2021_fig5_dicarbonyls_unbuffered | Met-glucose Amadori +- Met | hexose dicarbonyls (3-DG, 1-DG, GO, MGO) | -- | Jiangnan | out of scope (hexose) |
| schieberle1998_t4_ribose_aqueous | the Hofmann 1998 ribose pot | furfural 5.2 ug per pot vs Hofmann 67.5 (CF-016) | SIDA claimed | DFA Garching | read for CF-016 only; not transcribed (same pot) |

Not yet read (pentose systems with these species in the ledgers; candidates for a second pass): zhou2022 Fig. 1
(AX-ARP dicarbonyls; a third reading of the Zhou 2024 / 2025b pot, which could settle the 10x conflict), zhou2023c
Fig. 4 (ARP + GSH), tang2019 Figs. 5 / 7 (GSH-xylose), lu2019 Fig. 6 (Gly-Gly-arabinose), nakamura2020 Fig. 9
(xylose + lysine, NF and DX), zhang2025 Figs. 3 / S5-S7 (di-ARP), zhai2023c Fig. 5, mikami2017 Fig. 3B (xylose +
lysine NF), xu2008 Table 1 (NF endpoint vs pH).

What the table holds: 300 of 315 rows from one laboratory (Jiangnan), 8 Reading, 5 Ochanomizu, 2 Garching, so a
leave-one-laboratory-out test has three small hold-out groups. No laboratory in the corpus measures
1,4-dideoxypentosone or 2,3,4-pentanetrione. 2,3-Pentanedione over time in a pentose system exists in one unread
source: Zhang 2025 Tables 1-3 / S2-S4 (di-xylulosyl-lysine ARP 20 mM +- xylose 20 mM, 100-140 C, HS-SPME with an
authentic standard, ug/L +- SD), a lysine system; the engine has no 2,3-pentanedione species. Desclaux 2006 names
pentanedione among its measured dicarbonyls but does not show it. 3-DX / 1-DX absolute levels rest on an assumed
internal-standard response in every Jiangnan paper.
## Data after adjudication (2026-10-29: two readers, reconciled, adjudicated; table READY)

The paragraphs above describe the first reading. The adjudicated table holds **56 pots, 483 rows**: Jiangnan 52 pots
(284 semi-quantitative, 180 calibrated), Reading 1 (8 semi-quantitative), Ochanomizu 1 (5 calibrated), Dalian 1
(4 semi-quantitative), Garching 1 (2 calibrated). By species: 3-DX (TDP) 146, MGO 135, 1-DX (DPO) 130, GO 44,
furfural 21, norfuraneol 7. Decisions that change the fit:

- **CF-062 fixed:** every 3-DX / 1-DX series from the eight Jiangnan papers is semi-quantitative (internal standard,
  no authentic deoxyosone standard) and enters shape-only (`fitting.SCALE_TAU`); MGO and GO there are calibrated
  (purchased standards). Seven mixed figures were split into calibrated and semi-quantitative ledger ids.
- **GO from Zhai 2021 / 2023 excluded:** every plotted value is below the paper's own LOQ (CF-061).
- **Counted once:** the Zhai 2021 -> 2023 Food Chem -> 2023 JAFC triple print, by transcription (CF-066); the JAFC
  3-DX 140 C / 120 min point is a copy slip and is dropped.
- **CF-060 open:** Zhou 2024 vs Zhou 2025b, same nominal pot, every dicarbonyl ~10x apart; 3-DX spans ~10x across
  four Zhou papers. Until the laboratory answers, the Zhou calibrated rows get a per-paper scale in the fit.
- **Double use to avoid:** Hofmann 1998's rows are also B2.3 fit rows and Mikami 2017's norfuraneol decay is B61's
  measured constant; the sub-model 3 fit frees neither alongside those.
- 1,4-dideoxypentosone and pentanetrione: nothing measures either (only presence in Nedvidek 1992).

## Targets from the conflict analyses (2026-10-29, docs/reference/conflict_analyses/)

- **Norfuraneol at low pH (CF-012):** the shipped engine's norfuraneol at pH 3 is 55x below its pH 5 value, while
  three buffered laboratories see MFT (made from norfuraneol) fall only gently with pH. A no-refit probe must check
  norfuraneol's pH dependence first; the measurement that settles it is a norfuraneol pH series at 140 C (the designed
  protocol's pH 4 cell).
- **Furfural (CF-016):** Hofmann 1998's furfural (internal standard) is 13x Schieberle 1998's for the same run;
  Meynier's and Cerny's furfural / MFT ratios side with Schieberle. Until settled, Hofmann's furfural row enters with
  the conflict's width, not its printed 10 %.

## 3. Cheap probe (step 3) -- done 2026-10-29, nothing refitted

`scripts/generators/build_submodel_probe.py --table deoxyosone_submodel` runs every pot of the READY table on the sulfur lane at the shipped constants and predicts each row in its printed unit (molar masses for ug/L; semi-quantitative series compared in shape only, their scale removed). Two arms: shipped, and with the four pre-registered steps that ship OFF switched on, nothing refitted (B58 TTCA equilibrium, B59b cysteine Amadori route, B59a Strecker sulfur, B61 norfuraneol hydroxide decay). Residual = predicted minus measured (dex). Artifacts: `kinetic_core_deoxyosone_submodel_probe{,_with_b58_with_b59b_with_b59a_with_b61}.{json,md}`.

**Findings.**
1. **Deoxyosone SHAPES are right; their levels are unknowable here.** 3-DX and 1-DX are semi-quantitative (CF-062)
   and compared in shape only: 3-DX within 0.02-0.40 dex RMS across eight data sets; 1-DX within 0.10-0.41 with the
   route switches on (it was 0.7-2.5 shipped: B59b's route through the cysteine Amadori compound is what gives 1-DX
   its time course).
2. **Methylglyoxal is 10^2-10^4 too low in every pentose system, on calibrated data:** Zhai 2021 Fig. 5 (-3.9 dex
   mean with the switches on), Zhai 2023 Food Chem (-4.1), Zhai 2023 JAFC (-4.3), Zhou 2022 (-2.1), Zhou 2023c
   (-1.9), Zhou 2024 (-2.8), Zhou 2025b (-1.8; xylose alone -3.1). **Corrected 2026-10-29 (literature check):** the
   pentose C2 + C3 step EXISTS (`r_dpo_c2c3`, 1-deoxypentosone -> hydroxyacetaldehyde + MGO) but ships at log10 k
   -3.0 at 145 C against k_dpo_nf +0.43, i.e. ~3e-4 of 1-DX's exits; B2.3 / B60 fitted it on MFT rows only and no
   MGO row ever touched it. So this is a branching constant (and possibly its parent), not a missing node. The pH
   ladder (Zhai 2021 Fig. 5: MGO 13 -> 57 uM with 1-DX 25 -> 74 while 3-DX falls 520 -> 86) says MGO tracks the
   1-deoxy branch. 1-deoxypentos-2,3-diulose has no direct retro-aldol to C3 + C2 (its retro-aldol gives C4 + HCHO),
   so the parent is open. Design: three pre-registered single-step topologies -- (a) refit `r_dpo_c2c3`, (b) a
   base-favoured PENT / ARP -> MGO + HA step with an amine-free twin, (c) norfuraneol -> MGO (Mikami 2017) -- fitted
   on the 135 calibrated MGO rows with the MGO sinks held, judged against the 1-DX pH trend, the norfuraneol rows
   and Cerny's >= 93 % intact-C5 MFT labelling. Glycolaldehyde (HA) has no rows anywhere: Hofmann 1999 (doi
   10.1007/s002170050468) would pin it.
3. Glyoxal (44 rows) is not an engine species; furfural from the xylose-alanine ARP at 120 C is ~8x low (Zhou
   2025b, calibrated); norfuraneol (Hofmann, Mikami, Desclaux) is within 0.02-0.12 with B61 on.
## 4. Design (step 4) -- first draft 2026-10-30, SUPERSEDED by sec. 5 (kept as written)

**The joint hypothesis (lead, checkpoint 2).** 1-deoxypentosone (DPO) leaves the lane mainly by
`r_osone_decay_dpo` (DPO -> the lumped acid pool `ACID` + C4) and only ~3e-4 by `r_dpo_c2c3` (-> MGO + HA). Two
independent misfits point the same way: pentose systems are OVER-acidified (sub-model 1: Cao / Zhao 0.8-0.9 pH units
too acid at 0.4 M; the smoke fit pushed the acid yield per sink event toward zero) and MGO is 10^2-10^4 too LOW here.
If DPO's branching is wrong -- fragmentation to C3 + C2 underweighted, acid formation overweighted -- both close
together. Prediction: one shift of branching moves pH up in Cao / Zhao AND MGO up in the Jiangnan pots, with the
1-DX shape unchanged (its total loss rate fixed). It is falsified if the MGO rows need a shift that leaves the
pentose pH where it is, or that drives the pH above the data.

**Arms (each fitted on the 135 calibrated MGO rows; the MGO sinks held at shipped values; 1-DX / 3-DX in shape only):**
- (a) refit `k_dpo_c2c3` alone (one rate + Ea);
- (a') the joint form: `k_dpo_c2c3` and `k_osone_decay` (DPO branch) together, DPO's total exit held -- scored on
  the MGO rows AND the sub-model 1 pentose pH rows (Cao, Zhao, Zhou), the cross-check above;
- (b) a base-favoured PENT / ARP -> MGO + HA step with an amine-free twin (the xylose-alone pot), parallel to (a);
- (c) norfuraneol -> MGO (Mikami 2017's measured decay, its MGO share a lower bound) -- expected insufficient alone.
Discriminating constraints written now: MGO must track 1-DX across Zhai 2021's pH ladder (5.5-8); the norfuraneol
rows (Hofmann T5, Desclaux, Mikami) must not drop when DPO is diverted; Cerny's >= 93 % intact-C5 MFT at 95 C bounds
the HA made 1:1 with MGO; Zhai 2023 JAFC's three temperatures identify the Ea. CF-060 (Zhou 2024 vs 2025b, ~10x)
gets a per-paper scale on the Zhou calibrated rows. Glyoxal is not an engine species and stays out.

**Expected result (written now):** (a') beats (a) and lifts the pentose pH by >= 0.3 units (P = 0.45); (b) is
needed for the xylose-alone and the pH-8 pots (P = 0.35); no single arm brings every MGO data set within x3
(P = 0.6).
## 5. Design review (step 5) -- done 2026-10-30, before any fit; the design above is superseded

An independent reviewer (read-only; two small diagnostics at the shipped B60 constants with B58 / B59a / B59b / B61
on) found the drafted design unworkable. Recorded as found, then the revised design.

**BLOCKING**
1. **The joint "DPO -> acid" hypothesis is false in the shipped engine.** The lumped acid pool enters the balance at
   the yield 7.5e-5 per sink event; in Cao's pot it adds ~1e-4 mM of acid, while the carried carboxyls (CBX) hold
   ~70 mM. `k_osone_decay` = 0 leaves Cao's 5 h pH at 3.96; `k_dpo_c2c3` x 1000 moves it +0.01. DPO leaves mainly by
   the cysteine Strecker step `r_dpo_ddp` (11.4 L/mmol/min x [Cys] x a base factor at 100 C, pH 7), not by the acid
   sink (4.7e-4 /min) -- the draft's premise was wrong. The coupling that EXISTS is `r_dpo_ddp`: it destroys
   cysteine's amine and carries its carboxyl into CBX; switched off, Cao's 5 h pH goes 3.96 -> 4.20 (data 4.15). This
   is also why sub-model 1's acid yield cannot close the pentose over-acidification (its fit's pentose block is
   expected to fail for a reason outside its parameters).
2. **Per-pot scales hid the largest misfit.** Shape-only per series (ledger|pot|species) gave every pot its own
   scale and discarded all cross-pot contrasts. The hidden 1-DX offsets on Zhai 2021's ladder run -0.5 dex (pH 5.5)
   to -4.0 dex (pH 8): no response factor moves 3.5 decades with pH; the model's 1-DX collapses at high pH. With
   `k_dpo_ddp` x 1e-3, 1-DX at pH 7 / 8 goes -3.4 / -4.2 -> -0.8 / -1.4 dex and MGO -4.9 / -5.8 -> -2.7 / -3.1:
   about half the MGO gap is a PARENT deficit, not a branch. Probe finding 1 above ("shapes are right") holds only
   within each pot.
3. **The textbook MGO route is missing from the rivals:** 3-deoxypentosone's C3-C4 retro-aldol (MGO + glycolaldehyde,
   base-catalysed; the hexose analogue `r_tdg_mgo` is in the engine). DPO -> MGO + HA has no documented direct
   mechanism (reviewer ~70 %); the note citing Hofmann T10 for it cites the downstream recombination, not the split.
4. **Glycolaldehyde has no sink but MFT formation.** Raising C3 + C2 flux 10^2-10^4 would over-produce fragment-route
   MFT or force B2.3's `k_ha_mp_mft` / `k_mgo_mp` down: an HA sink with a declared band (Hofmann 1999 would pin it)
   is needed, or MFT stays out of the arm comparison.
5. **Cerny 2003 is a declared hold-out** (FIT_HOLDOUT_DECLARATION.md) with no rows here: Cerny 2004's in-situ
   branching (declared FIT) is the constraint; Cerny 2003 is scored once, after freezing.

**SHOULD-FIX (adopted):** production, not sinks, is short (~10^4; `r_mgo_mel`, B1-fitted, dominates the sink and
agrees with the decline phases; its band is propagated); the pH-tracking constraint written as a statistic (predicted
d log MGO / dpH at 120 min vs measured 0.13 / 0.16 dex per pH unit); the ladder numbers in sec. 3 to be recomputed
from the table; furfural added as a falsifier; leave-one-laboratory-out CANNOT be run for MGO (132 / 135 rows
Jiangnan) -- declared not run, leave-one-paper-out instead, and a search for a second laboratory's pentose MGO time
course; CF-060: each Zhou paper's MGO its own scale, Zhai absolute, plus a sensitivity with Zhai free; possible OPD
derivatisation artefact (4-12 h at room temperature with mM ARP present) to check against any unheated control;
the engine state named and frozen (route switches, sub-model 1's switches) and the winner re-scored in the other;
SDs of 1-3 % get the per-data-set noise; information counted by pots (~30), not rows.

**Revised design (supersedes sec. 4):**
- Step 0 (before any fit): scales shared per (ledger id, species) across a figure's pots with a finite tau of 0.3
  dex for the semi-quantitative quinoxalines (choices register: tau 0.15 / 0.3 / inf); an HA sink with a declared
  band; Mikami's measured NF -> MGO share (>= 0.6 mol %) inside B61's decay in every arm; the engine state named;
  then re-probe.
- Arms: A0 null; A1 `k_dpo_c2c3` alone (expected to fail or to starve NF and DDP); A2 `k_dpo_ddp` + `k_dpo_c2c3`
  (the parent arm; cross-check, unfitted, that Cao / Zhao's pH rises >= 0.2 at 5 h); A3 3-DX -> MGO + HA (rate, Ea,
  hydroxide factor); A4 PENT / ARP / CAMA -> MGO + HA, base-catalysed, with an amine-free twin.
- Each scored on: the MGO data sets; the 1-DX / 3-DX contrasts across pots; the pH-slope statistic; the +/- Cys and
  +/- Ala ratios; furfural; the norfuraneol rows; every re-scored B60 fit row and hold-out (the constants they move
  were frozen there).
- Done-test tolerances (declared now): Zhai MGO absolute |mean| <= 0.3 dex and RMS <= 0.4 per data set with no
  one-sign flag in any pH or temperature bin; Zhou MGO in shape RMS <= 0.3; 1-DX / 3-DX per figure RMS <= 0.35;
  pH slope within +/- 0.1 dex per pH unit; +/- Cys ratio within x2; B60 fit rows no more than 10 % worse (whitened
  RMS); Cerny 2003 fragment-route share <= 7 % (scored once); constants off their bounds with profile intervals
  under 1 dex. Hold-outs: Zhai 2023 JAFC's 140 C pots and the Zhou 2023 2-acetylthiazole pair.
- Expected (revised): an arm with parent freedom (A2 or A3) beats A1 (P = 0.75); the `k_dpo_ddp` arm lifts Cao by
  >= 0.2 pH (P = 0.6); some arm brings every Zhai MGO set within x3 (P = 0.4).

**Cross-sub-model consequence.** `k_dpo_ddp` sits on the intact-C5 MFT route (DDP) and on cysteine's loss: lowering
it moves sub-model 2 (cysteine) and sub-model 4 (MFT) too. Sub-models 2-4 are coupled through this one constant;
it is fitted here, and the route and thiol tables are re-scored with it before either is fitted.
## Results (2026-10-30, 05:40; arms A2 and A3 -- A4 / A1 not run tonight)

Fits: `kinetic_core_deoxyosone_submodel_fit_{A2,A3}.md` (engine state: B60 + B58 / B59b / B59a / B61 on; sub-model 1
switches off; per-figure scales, tau 0.3; Zhou MGO per-paper scale).

| | A2: k_dpo_ddp + k_dpo_c2c3 | A3: 3-DX -> MGO + HA |
| --- | --- | --- |
| cost (start -> optimum); reduced chi-square | 1114 -> 640; 3.2 | 1245 -> 845; 4.7 |
| constants | k_dpo_ddp x 10^-4 (ON ITS BOUND); k_dpo_c2c3 x 10^2.47; NF -> MGO share 1.4 % | k_tdp_mgo 10^-5.1 /min at 100 C, pH 5; Ea 40 (on bound); NF -> MGO 1.4 % |
| Zhai MGO, absolute (tolerance mean <= 0.3, RMS <= 0.4) | 2021 ladder -0.04 / 0.13; Food Chem +0.18 / 0.34; JAFC +0.21 / 0.35 -- all within | +0.35 / 0.58 to -0.70 / 1.59 -- outside |
| Zhou MGO, shape (RMS <= 0.3) | 0.50-0.82 -- outside | 0.53-2.16 -- outside |
| 1-DX per figure (RMS <= 0.35) | Zhai 2021 0.47, JAFC 0.44, Zhai 2019 0.81, Food Chem 1.21, Zhou 0.5-1.6 | 0.6-3.2 |
| 3-DX per figure | 0.4-1.6 (Zhai over-predicted +0.5 to +1.2) | similar |
| hold-out (Zhai 2023 JAFC, 140 C) | MGO -0.45 / 0.62; 1-DX -1.67; 3-DX +0.16 / 0.47 | MGO -0.68 / 0.99; 1-DX -2.58 |
| norfuraneol (Mikami, Hofmann) | within 0.01-0.27 | within 0.01-0.02 |

**Reading.** A2 -- the parent arm the review recommended -- closes the calibrated MGO gap on Zhai's three data sets
(from -2 dex to within tolerance) and halves to quarters the 1-DX misfit across the pH ladder, by making 1-DX's
cysteine Strecker step (`k_dpo_ddp`) 10^4 slower and its C3 + C2 split ~300x faster. A3 (the textbook 3-DX route)
does worse on every MGO set and leaves 1-DX untouched: the data favour the 1-deoxy parent, as the pH ladder said.
Mikami's NF -> MGO share fits at 1.4 % (his printed 0.6 % is a lower bound), consistent. NOT YET ACCEPTED: (1)
`k_dpo_ddp` sits ON its lower bound -- the data want the step even slower or gone, and it carries cysteine loss,
the intact-C5 MFT route (DDP) and the Strecker H2S; the review required every B60 MFT / FFT row and the route and
thiol tables re-scored before acceptance (running, unfitted). (2) Leave-one-paper-out not run (time). (3) 3-DX
remains over-predicted in Zhai's pots and the hold-out's 1-DX at 140 C is 47x low: a temperature dependence the
held Ea of the carbonyl family does not carry. Expected result (sec. 5): A2 or A3 beats A1 -- A1 not run; A2 lifts
Cao's pH -- being checked; some arm brings every Zhai MGO set within x3 -- yes (A2, P was 0.4).

**Cross-checks of A2, unfitted (2026-10-30, 06:00; `kinetic_core_{thiol,route,ph_acid}_submodel_probe_sm3{A2,null}_*`).**
Route table: unchanged (cysteine, TTCA, CAMA within 0.0-0.1 dex of the switches-on probe). Sub-model 1's pentose
pots: Cao's pH -0.20 (null, route switches on) -> -0.03 / RMS 0.07; Zhao +0.11 / 0.15 -- both within tolerance (the
review's prediction held); Zhou's ARP endpoints worse (+0.39 -> +0.81). Thiol table: A2 harms the MFT data it feeds
through DDP -- Schieberle 2000 MFT -1.24 -> -1.55 dex and FFT -0.16 -> -0.67; Yiltirak MFT +0.53 -> +1.51; Zhou
2022 MFT +0.15 -> +0.64; Hofmann's norfuraneol +0.02 -> +0.27 (1-DX no longer drained by cysteine goes to NF). The
glucose rows' zero MFT is untouched (the missing hexose -> C5 entry). **Verdict: A2 is not accepted alone.**
`k_dpo_ddp` must be fitted JOINTLY with the steps it feeds (k_ddp_mft, the DDP and NF exits; sub-models 3 + 4 on
the deoxyosone and thiol tables together), with the MGO rows and the MFT rows both in the cost. That is the next fit.

## Amendments
