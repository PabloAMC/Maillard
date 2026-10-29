# Sub-model 4: thiol formation and loss (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Do MFT and FFT form from their precursors (norfuraneol / dideoxypentosone and furfural with H2S) and
disappear (decay, dimerisation, binding) at the rates the new time series show, and is H2S supply and loss right?

**Decision.** Whether the unpinned thiol constants (`k_fur_fft_hs`, `k_mft_decay`, `k_ddp_mft`, `k_ddp_mft_hs`,
`k_h2s_loss`; B2.3) can be pinned; how the laboratories' rising-vs-falling thiol curves (CF-013) are explained.

**Rival arms.** Thiol loss as first-order decay vs dimerisation (second order in thiol) vs binding to melanoidins
or protein (Gigl 2021); semi-quantitative series (the Zhai / Kang ladder, CF-004) read by shape only vs the absolute
ones (Feng 2022, isotope-dilution data) alone.

**Expected result.** Formation and loss of FFT and MFT become identified on the time series (~55 %); the H2S budget
stays the weakest link (~60 % unpinned), and the binding question needs the matrix layer (Phase 4).

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 3 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2) -- first reader done 2026-10-29 (agent-sm4-reader1), second reader pending

Verification: `data/lit/source_verification/thiol_submodel.yml` (24 entries, every check with its quote; read against
the PDFs / SI, not the dossiers). Keyed table: `data/lit/fit_rows/thiol_submodel.yml`, `status: draft`, **40 pots,
119 rows, 5 laboratories**: 63 absolute (Garching SIDA / IS-calibrated 33, Jiangnan external curves 30), 50
semi-quantitative with response factor 1 (within-compound shape only), 6 heated/unheated survival ratios (response
factor cancels). Verdict: P = all five checks pass, C = passed after a correction, F = a check failed.

| ledger id | system | what | quantification | laboratory | verdict | rows |
|---|---|---|---|---|---|---|
| hofmann1998_t1_cys_carbohydrates | Cys 33 + ribose / xylose / glucose / fructose 100 mM, 0.5 M phosphate pH 5, 145 C, 20 min | FFT, MFT | SIDA | Garching (Hofmann / Schieberle) | P | 8 |
| hofmann1998_t2_ph_ladder | same, ribose / glucose at pH 3 and 7 | FFT, MFT | SIDA | Garching | P | 8 |
| hofmann1998_t3_fft_c5_h2s | H2S 20 + ribose / 3-deoxyribosulose / furfural 20 mM, pH 5, 145 C | FFT | SIDA | Garching | P | 3 |
| hofmann1998_t4_mft_ribose_nf_h2s | H2S 20 + ribose / NF 20 mM, pH 5, 145 C | MFT | SIDA | Garching | P | 2 |
| hofmann1998_t5_fa_nf_pairs | Table 1 pots (ribose, glucose) | furfural, NF | IS with measured RF (ledger "unknown") | Garching | C | 4 |
| hofmann1998_t8_c2c3_mft_fft_ph | hydroxyacetaldehyde + mercapto-2-propanone 20 mM, pH 3 / 5 / 7, 145 C | FFT, MFT | SIDA | Garching | P | 6 |
| hofmann1998_t10_mft_routes | NF + Cys, thiamin 20 mM, pH 5, 145 C | MFT | SIDA | Garching | P | 2 |
| schieberle2000_t4_100c_time_course | Rib 100 + Cys 33 mM, pH 5, 100 C, 30-720 min | MFT, FFT time course | SIDA inferred (not stated) | Garching | C, no error printed | 0 |
| feng2022_table1_rgarp_volatiles | ribose-GSH Amadori 20 mM, water pH 7, 120 C 30-120 min + 100 C 120 min | MFT, FFT, furfural | external curves (SI S1) | Jiangnan (Zhang X / Cui / Ho) | P | 15 |
| feng2022_table2_rgarp_gsh_volatiles | same + GSH 20 mM | MFT, FFT, furfural | external curves | Jiangnan | C | 15 |
| zhai2023foodchem_table1_ttca_tladder | TTCA 10 mM, water pH 7, 100 / 120 / 140 C, 120 min | MFT, FFT, furfural | IS, f' = 1 (CF-004) | Jiangnan | P | 9 |
| zhai2023foodchem_table1_ttcaxyl_tladder | TTCA 10 + xylose 10 mM, same | MFT, FFT, furfural | IS, f' = 1 | Jiangnan | P | 9 |
| kang2026_si_table_s4_ttca_tladder | reprint of Zhai Table 1 (CF-005) | -- | f' = 1 | Jiangnan / Zhejiang Gongshang | F (duplicate) | 0 |
| zhai2023foodchem_fig1a_100c_timecourse_grid | TTCA +/- Xyl, 100 C, 20-120 min | heat map | f' = 1 | Jiangnan | F (cells not verifiable, no error) | 0 |
| wang2026b_tables1_22_volatiles | Xyl 200 + Cys 200 mM, 0.2 M NaH2PO4 pH 5.5, 120 C, 60-150 min | MFT, FFT, MFTD, FFTD, mixed disulfide, furfural | IS, RF 1 | Dalian Polytechnic (Lin / Chen) | P | 23 |
| wang2026b_table1_oav_thiols_dimers | OAVs of the same data | -- | -- | Dalian | F (1e6 unit error, copy) | 0 |
| zhang2026_tables6_s8_volatiles | GSH 60 + glucose 60 mM, 0.2 M NaH2PO4 pH 6.5, 130 / 140 / 150 C, 10-120 min | MFT (90, 120 min), furfural (120 min) | IS, RF 1 | BTBU (Xie J) | P | 9 |
| kumazawa2003_fig3_fft_survival_ph_grid | FFT 8.8 uM, McIlvaine pH 3-6.4, canned (air), 121 C, 10 min | FFT residual % | heated / unheated ratio (RF cancels; ledger semi-quant) | Ogawa & Co. | C | 6 |
| kumazawa2003_t2_fft_survival_t_time_ph | same, pH 6.0 / 6.8, 121 C 10 min vs 123 C 20 min | FFT residual % | ratio | Ogawa | C, no error printed | 0 |
| kumazawa2003_t3_fft_products_ph_grid | Fig. 3 pots | FFT, difurfuryl disulfide area ratios | area ratio, RF 1 | Ogawa | F (no unit, no error) | 0 |
| gigl2021_fig2_fig4_coffee_melanoidin_fft_binding | FFT 5 mM + coffee HMW melanoidins, phosphate pD 5.5, RT, 1 min-96 h | % free FFT | qHNMR | TU Munich (Hofmann / Frank) | P, no error printed | 0 |
| zhang2024_fig1_thiamine_xylose_redox_series | thiamine + xylose +/- Cys / GSH / cystine, pH 4.9, 60 min | MFT, MFT-MFT | external curves | BTBU (Cao Y) | F (temperature not stated) | 0 |
| zheng1994_table1_cysteine_h2s_release_k | Cys 0.1 M, pH 3-9, 80-110 C, N2 headspace | H2S release k | sulfide ISE, Na2S curve | Rutgers (Ho) | P, constants not rows | 0 |
| zheng1994_figs1_4_cysteine_ln_ct_c0 | same | ln(Ct/C0) | ISE | Rutgers | F (not digitised, no unit) | 0 |

Not on disk: Li, Schieberle & Steinhaus 2016 (34S-sulfide SIDA for H2S); `data/articles/li2016.pdf` is a different
paper (Li, Lund et al., quinone-thiol kinetics) and the ledger has no entry. The table holds **no H2S measurement**:
the H2S arm rests on Zheng & Ho's constants and Hofmann's fed-H2S pots.
## Notes from the first reading (2026-10-29)

- **No H2S measurement on disk.** The table has 63 absolute rows (Hofmann 1998 SIDA / internal standard, Feng 2022
  external curves), 50 semi-quantitative (shape only: Wang 2026b, Zhai 2023, Zhang 2026) and 6 Kumazawa survival
  ratios, and nothing that measures H2S. So `k_h2s_loss` cannot be pinned from this table: its expected-result line
  above ("H2S stays the weakest link") becomes ~85 % unless Li, Schieberle & Steinhaus 2016 (Eur. Food Res. Technol.
  243:69, doi 10.1007/s00217-016-2723-z, SIDA for H2S) is obtained; `li2016.pdf` on disk is a different paper.
  Otherwise the designed protocol's H2S analyte is the route.
- Feng 2022's RG-ARP carries sulfur, the engine's ARP does not (CF-058): map it before its rows enter.
- Schieberle 2000's "factor of 13" is likely a column slip for 7.9 (CF-055), which the B16 row used.

## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
