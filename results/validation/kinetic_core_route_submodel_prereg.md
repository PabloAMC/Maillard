# Sub-model 2: the cysteine-pentose route (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Does de Roos's measured route (TTCA <-> Cys + pentose; TTCA -> cysteine Amadori compound -> 1- and
3-deoxypentosone + Cys) reproduce the time courses of TTCA, the Cys-Amadori compound, free cysteine and the pentose
across the three laboratories that measured them, and does it replace the shipped step that treats cysteine as an
instant catalyst on the pentose (`r_pent_tdp`, `r_pent_dpo`: 59 % and 45 % of a reference pot's flux, unpinned since
B2.3; `docs/reference/unpinned_constants.yml`)?

**Decision.** Whether the route ships (switch on) and the two catalyst steps retire; which reading of phosphate
catalysis vs pH the data support (stage 1's 1.87x Cao-vs-Zhai tension, CF-011).

**Rival arms.** (a) Phosphate catalyses the Amadori steps (Davidek 2002's 1.3-3.5x) vs (b) the difference is pH
alone vs (c) laboratory. Above pH 6 the decay is extended or flat (CR-14).

**Expected result (before any probe).** The route reproduces Cao 2017 and Zhai 2019's intermediates within the
measurement errors (stage 1 got there with the procedure); the phosphate question stays open on literature data
(~60 %), and the designed experiment's buffer arm is what settles it.

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 1 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2)

First reader done 2026-10-29 (agent-sm2-reader1), against the papers; second reader pending. Checks per source:
`data/lit/source_verification/route_submodel.yml`. Rows: `data/lit/fit_rows/route_submodel.yml` (status: draft;
45 pots, 423 rows; the loader refuses it until both readers sign). Laboratory = the group that produced the data,
not the paper's affiliation.

| ledger id | system | what | quantification | laboratory | verdict | rows |
|---|---|---|---|---|---|---|
| cao2017_fig5_ttca_cama_glyarp | Cys + Xyl 0.40 M, 0.2 M Na phosphate pH 5.5, 90 °C, 0.5-5 h | TTCA, CAMA (blank arm) | calibrated (LC-MS SIM, synthesized standards) | BTBU (Xie) | pass | 12 |
| cao2017_fig4ab_cys_xyl | same pot | free Cys, Xyl (blank arm) | calibrated (HPLC-ELSD) | BTBU (Xie) | pass | 12 |
| zhao2017b_fig2abde_xyl_cys_ttca_cama | same recipe, a sibling run | Xyl, Cys, TTCA, CAMA (blank) | calibrated (HPLC-ELSD) | BTBU (Xie) | pass | 23 |
| zhang2025b_fig2abc_ttca_cama_lysarp | Cys + Xyl 40 mM, 0.2 M phosphate pH 6.5, 120 °C, 5-50 min | TTCA, CAMA (blank) | calibrated (LC-MS/MS MRM) | BTBU (Xie) | pass | 15 |
| zhai2019_fig3a_ttca_cama | Xyl 165 + Cys 82.6 mM, water, pH 7.4 start, 70-130 °C, 20-160 min | TTCA, CAMA | calibrated (HPLC-ELSD) | Jiangnan (Zhang X) | pass | 62 |
| zhai2019_fig3b_dicarbonyls_ttca | TTCA 5 mM ± Cys 5 mM, water, pH 7.4 start, 120 °C | 3-DX (TDP), 1-DX (DPO) | **semi-quantitative** (internal standard, no authentic standard; ledger says calibrated) | Jiangnan | corrected | 24 |
| zhai2021_fig3a_ttca_tladder | TTCA 10 mM, water, natural pH 5.5, 120 / 140 °C | TTCA | calibrated (HPLC-ELSD) | Jiangnan | pass (100 °C = Fig. 4b) | 12 |
| zhai2021_fig4b_ttca_ph_ladder | TTCA 10 mM ± Xyl 10 mM, water, pH 5.5 / 6 / 7 / 8 start, 100 °C | TTCA | calibrated | Jiangnan | pass | 41 |
| zhai2021_fig5_dicarbonyls_ph_ladder | same systems | 3-DX, 1-DX | **semi-quantitative** (diacetyl internal standard; ledger says calibrated) | Jiangnan | corrected | 96 |
| zhai2021_fig6_cys_xyl | same systems + Xyl-alone controls | released Cys, Xyl | calibrated (HPLC-ELSD) | Jiangnan | corrected (conditions missing in the ledger); vector-read here | 96 |
| kang2026_si_fig_s3_cys_release_from_ttca | TTCA 10 mM, water, pH 7 start, 120 / 140 °C | released Cys | calibrated (the original's) | **Jiangnan** (= Zhai 2023 Food Chem Fig. 3a; ledger says Zhejiang Gongshang) | corrected | 12 |
| kang2026_si_fig_s4_free_cys_depletion | Cys 10 mM alone, water, pH 7 start, 100-140 °C | free Cys | calibrated (the original's) | **Jiangnan** (= Zhai 2023 Food Chem SI Fig. S2) | corrected | 18 |
| zhai2023foodchem_fig2_deoxyosones / _fig3_cys_da_mgo_go | TTCA 10 mM ± Xyl, water, pH 7, 100-140 °C | 3-DX, 1-DX; Cys | semi-quant. (3-DX / 1-DX); calibrated (Cys) | Jiangnan | pass / corrected; **not transcribed** (no digitisation; 100 °C series = Zhai 2021, TTCA-alone Cys = Kang S3) | 0 |
| deroos2005_fig4_xylulose_cysteine_ph_stability | purified pentose CAMA, boiling 0.5 M phosphate, pH 3-6, 2 h | % CAMA left | not stated (a ratio) | Givaudan (de Roos) | pass; **not transcribed** (no unit for a fraction, no concentration, no points) | 0 |
| du2024_fig2_cama_aldehyde_thiazolidines | Cys + glucose (+ aldehydes) | glucose Cys-Amadori | calibrated | CAAS (Wang F) | **fail** (conditions not stated; glucose, not pentose) | 0 |
| zhou2025b_figs1_3_axarp_xylose_120c | alanine-xylose ARP 20 mM ± Xyl, pH 8, 120 °C | ARP, 3-DX, 1-DX | calibrated (ARP), semi-quant. (DX) | Jiangnan | pass / corrected; **not transcribed** (figure only; alanine route, outside this sub-model) | 0 |

**What the table holds.** 45 pots, 423 rows: TTCA 105, Cys 90, PENT 60, TDP 60, DPO 60, CAMA 48. Two
laboratories, not three: Jiangnan 361 rows (Kang's SI figures are Jiangnan data re-published), BTBU 62. The phosphate
question (CF-011) is therefore a two-laboratory contrast with buffer, pH and concentration all confounded with
laboratory; leave-one-laboratory-out has two folds. 154 rows carry sigma < 0.01 dex (small printed SDs or vector read
errors with no SD), so the fit's measurement model must carry its per-series scale or a floor.

**Duplicates found (count once).** Zhai 2021 Fig. 3a 100 °C = Fig. 4b pH 5.5; Zhai 2021 Fig. 6a pH 7 (TTCA and
TTCA-Xyl) = Zhai 2023 Food Chem Fig. 3a 100 °C; Zhai 2021 Fig. 5 pH 7 = Zhai 2023 Food Chem Fig. 2 100 °C; Zhai 2023
Food Chem Fig. 3a TTCA arms = Kang SI Fig. S3; Zhai 2023 Food Chem SI Fig. S2 = Kang SI Fig. S4; Zhai 2019 Fig. 1B =
Fig. 3A 90 °C. Each is transcribed once.

**Before the table can be ready:** a second reader; ledger quant_class for zhai2019_fig3b, zhai2021_fig5,
zhai2023foodchem_fig2 (3-DX / 1-DX semi-quantitative); Kang S3 / S4 re-labelled to Jiangnan with same_data_as flags
(the loader also refuses the string "Zhejiang Gongshang (Fu Y / Zhai Y)", absent from `data/keys/laboratories.yml`);
zhai2021_fig6's conditions in the ledger. Left out pending a digitisation: Zhai 2023 Food Chem Figs. 2-3 at 120 / 140 °C
(TTCA+Xyl Cys, 3-DX, 1-DX: 60 points with SD bars), Zhai 2021 Fig. 7a (13C5-TTCA, no ledger entry).
## 3. Cheap probe (step 3) -- done 2026-10-29, nothing refitted

`scripts/generators/build_submodel_probe.py --table route_submodel` runs every pot of the READY table on the sulfur lane at the shipped constants and predicts each row in its printed unit (molar masses for ug/L; semi-quantitative series compared in shape only, their scale removed). Two arms: shipped, and with the four pre-registered steps that ship OFF switched on, nothing refitted (B58 TTCA equilibrium, B59b cysteine Amadori route, B59a Strecker sulfur, B61 norfuraneol hydroxide decay). Residual = predicted minus measured (dex). Artifacts: `kinetic_core_route_submodel_probe{,_with_b58_with_b59b_with_b59a_with_b61}.{json,md}`.

| data set (species) | n | shipped RMS | switches on: mean / RMS |
| --- | --- | --- | --- |
| Cao 2017 (TTCA, CAMA, Cys, PENT) | 24 | CAMA never formed | -0.16..+0.13 / 0.03-0.31 |
| Zhao 2017b, same pot (TTCA, CAMA, Cys, PENT) | 23 | Cys +0.61 | -0.20..+0.19 / 0.07-0.22 |
| Zhai 2019 T-ladder (TTCA / CAMA) | 32 / 29 | 3.37 / never formed | -0.29 / 0.53 ; +0.71 / 0.78 |
| Zhai 2021 Fig. 4b pH ladder (TTCA) | 41 | 0.71 | -0.35 / 0.41 |
| Zhai 2021 Fig. 3a T-ladder pH 5.5 (TTCA) | 12 | 4.93 | -1.12 / 1.16 |
| Zhai 2021 Fig. 6 (Cys / PENT) | 48 / 48 | 1.20 / 1.13 | +1.07 / 1.15 ; +0.10 / 0.17 |
| Zhang 2025b (TTCA, CAMA, PENT, Cys) | 31 | PENT 10.0 | -0.10..+0.66 / 0.25-0.79 |
| Zhai 2023 Food Chem / JAFC (Cys / CAMA) | 24 / 45 | 0.40 / never formed | +0.43 / 0.65 ; -0.32 / 0.83 |

**Findings.**
1. **The shipped route is wrong because two measured steps ship off.** Without B58 (TTCA <-> Cys + pentose at its
   measured equilibrium) a charged TTCA can only decay (Zhai 2021 Fig. 3a: 4.9 dex), and without B59b the cysteine
   Amadori compound is never formed (92 rows predicted zero). With both on, and nothing refitted, the route lands
   within ~x3 on most rows of five papers from two laboratories. This is independent evidence for B58 / B59b; the
   evidence audit for switching them on belongs to Phase 3's single regeneration (strategy sec. 7.6), and this probe
   is its input.
2. **Open misfits with the switches on:** TTCA at 120-140 C, pH 5.5 decays ~13x too fast (Zhai 2021 Fig. 3a); free
   cysteine in Zhai 2021 Fig. 6 is lost ~12x too slowly (+1.07: the model keeps too much); the cysteine Amadori
   compound is ~5x too high in Zhai 2019 and Zhang 2025b and ~2x low in Zhai 2023 JAFC.
3. **Decision:** sub-model 2's fit is deprioritised. Its remaining misses are factor 2-13, while sub-model 3's
   methylglyoxal miss is 10^2-10^4 on calibrated data from seven papers (a missing route) and feeds the MFT path;
   that comes first (tasks/overnight_log_2026-10-29.md, checkpoint 1).
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
