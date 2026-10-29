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

## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
