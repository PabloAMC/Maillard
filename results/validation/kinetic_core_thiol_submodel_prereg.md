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

## Targets from the conflict analyses (2026-10-29, docs/reference/conflict_analyses/)

- **Temperature ordering (CF-017):** measured thiols order 100 < 120 < 140 C in 16 of 16 comparisons; the shipped
  engine gets 140 >= 120 C in only 1 of 8. A structural misfit (loss outrunning formation at high temperature): the
  first thing the no-refit probe reports, and a done-test item.
- **The MFT sink (CF-013):** US 10,076,126 added its labelled standards before heating and they survived, which
  bounds the real MFT loss: the engine's sink (k_mft_decay, unpinned since B2.3) is at least 5x too fast. A rival arm
  with the sink capped by that bound.
- **FFT across sources (CF-015):** MFT agrees across sources, FFT does not (the 3:1 xylose recipe favours FFT; one
  patent's standards were heated). FFT rows from the two patents enter with that width.

## H2S literature read (2026-10-29: four dossiers, 12 ledger entries)

The heated-system H2S search (2026-10-29) found no H2S time series in a Maillard pot; the four papers it recommended
are now read (`data/lit/extraction_dossiers/{li2022,obata1960,cuevasanta2015,suleimenov1994}_extraction.md`):

- **Li 2022** (free cysteine 5.7 mM, sealed 1 mL vial, headspace GC-FPD against Na2S standards, read cold): H2S
  3.8 +/- 0.3 uM at 90 C / 10 min, 10.7 uM at 120 min (single measurement). With no sink the shipped `k_cys_h2s`
  predicts 2.6 and 31 uM: the 10-min point agrees (1.5x), the 120-min point is 0.26-0.34x. A cold, sealed read against
  external standards rules out headspace partitioning as the cause, and storage loss would hit both times alike, so
  the data imply a time-dependent loss: a first-order sink of ~0.025 /min at 90 C, about 10x the shipped
  `k_h2s_loss`. One unreplicated point carries this; it is a lead for the probe, not a constant.
- **Obata 1960** agrees with `k_cys_h2s` within 0.9-1.3x (liquid near 100 C), and shows cystine itself releasing H2S
  (0.4x cysteine at 1 h, 1.3-1.5x at 3-6 h).
- **Cuevasanta 2015:** HS- + cystine at 0.6 M-1 s-1 (25 C, pH 7.4; pH dependence only through HS-). A real sink only
  where cysteine is oxidised (air headspace); negligible under N2 (Zheng & Ho's pots). No temperature dependence and
  no reverse rate measured.
- **Suleimenov & Krupp 1994:** H2S's Henry constant (gas/liquid ~0.9-1.1 by concentration at 90-145 C). The engine
  has no H2S headspace: in a 15 mL vial with 5 mL liquid at 120 C, ~31 % of total sulfide at pH 7 would sit in the gas.
  This is a mechanism the lane lacks, with a measured constant, for the design step.

Two rows enter the table once double-read: Li 2022's two free-cysteine H2S points (and the Cys + Met arm's
replicated pair, if Met is carried). Expected result for the H2S arm revised from ~85 % "stays the weakest link" to
~70 %: the sink's existence has a lead, its rate does not.

## 3. Cheap probe (step 3) -- done 2026-10-29, nothing refitted

`scripts/generators/build_submodel_probe.py --table thiol_submodel` runs every pot of the READY table on the sulfur lane at the shipped constants and predicts each row in its printed unit (molar masses for ug/L; semi-quantitative series compared in shape only, their scale removed). Two arms: shipped, and with the four pre-registered steps that ship OFF switched on, nothing refitted (B58 TTCA equilibrium, B59b cysteine Amadori route, B59a Strecker sulfur, B61 norfuraneol hydroxide decay). Residual = predicted minus measured (dex). Artifacts: `kinetic_core_thiol_submodel_probe{,_with_b58_with_b59b_with_b59a_with_b61}.{json,md}`.

**Findings.**
1. **Hexose + cysteine makes no MFT and no norfuraneol** (Hofmann 1998 T1 / T2 / T5: glucose, fructose; predicted
   exactly 0 against 19-25 ug/L MFT and 135 ug/L norfuraneol). **Corrected 2026-10-29 (literature check):** the lane
   DOES run the trunk's amine-free hexose steps (r_glc_tdg, r_fru_odg, r_tdg_mgo), so glucose + Cys makes MGO and
   1-mercapto-2-propanone; MFT is zero because hexose hydroxyacetaldehyde comes only from r_glc_c2c3 (at its band
   floor) and there is no hexose -> C5 entry. Labelling on disk (Du 2023: Cys + [13C6]glucose, MFT 0 : 5 = 52 : 48,
   mixed isotopomers < 1 %; Zhao 2019b with GSH: 50 : 47) says hexose MFT keeps an intact glucose C5: the missing
   step is a hexose -> C5 entry (C-6 or C-1 loss), not more C2 + C3. Candidate: ODG or Glc -> DPO (+ C1); falsified
   if it predicts the ribose NF / MFT ratio (~10^3) against Hofmann's 7 for glucose, in which case the entry is at
   1,4-dideoxypentosone (bypassing NF). Hofmann T9's mercaptopropanone from glucose + Cys (pH 3 / 5 / 7) is not yet
   a fit row and tests the existing hexose -> MGO flux.
2. **The pH direction of the thiols is wrong at pH 7** (ribose + Cys, Hofmann T2): MFT 46x low, FFT 46x high; at
   pH 3 both within 2x. CF-012 / CF-017 (conflict analyses) already point here.
3. **MFT at 100 C is 10-60x low early** (Schieberle 2000, isotope dilution), and FFT's shape is wrong (too low
   before 1 h, 3-5x high at 6-12 h).
4. The fed-precursor pots (H2S + furfural / norfuraneol / C2 + C3, Hofmann T3 / T4 / T8 / T10) are within 0.03-0.8
   dex: they were B2.3's fit data, so this is in-sample. Semi-quantitative Jiangnan series are within 0.2-0.8 dex
   in shape. Yiltirak 2026's buffered arm: FFT ~120x high.
5. **Order:** the thiol constants are not fitted before the two missing routes (finding 1 here, finding 2 of the
   deoxyosone probe) are defined, because a thiol fit would absorb them (the k_aa_frag lesson).
## 4. Design (step 4) -- the JOINT sub-model 3 + 4 fit (2026-10-30, owner: "run the joint submodel fits")

**Why joint.** Sub-model 3's best arm (A2: `k_dpo_ddp` x 10^-4 on its bound, `k_dpo_c2c3` x 300) closes Zhai's
calibrated methylglyoxal but harms the MFT / FFT rows (Schieberle 2000 MFT -1.24 -> -1.55 dex, FFT -0.16 -> -0.67;
Yiltirak MFT +0.53 -> +1.51; Zhou 2022 MFT +0.15 -> +0.64). The link is NOT the intact-C5 DDP route (`k_ddp_mft`
ships at 0) but the H2S budget: with B59a on, `r_dpo_ddp` releases a share of the cysteine sulfur as H2S, and
both MFT (from norfuraneol) and FFT (from furfural) need H2S. One constant cannot be set on one table without
breaking the other.

**Engine state (named):** B60 + the route switches B58 / B59b / B59a / B61 on; sub-model 1's switches off; inline
pH coupling (checked against the fixed point at the optimum).

**Coordinates (7, log10 multipliers on the shipped k_ref unless stated; Ea held at the carbonyl family's):**
`k_dpo_ddp`, `k_dpo_c2c3`, the NF -> MGO share (0-1), `k_dpo_nf` (DPO -> norfuraneol), `k_nf_mft` (norfuraneol +
H2S -> MFT), `k_fur_fft` (furfural + H2S -> FFT), `k_h2s_loss`. `k_ddp_mft` stays at its shipped 0 (B60 fitted it
there; the hexose / DDP question is a separate arm once the hexose -> C5 entry exists). Bounds: multipliers 10^-4 to
10^3.

**Rows:** the deoxyosone and thiol tables together (MGO, 1-DX, 3-DX, furfural, norfuraneol, MFT, FFT, mercapto-
propanone). OUT: glyoxal (no engine species); the hexose + cysteine pots (glucose / fructose: MFT structurally zero
until the hexose -> C5 entry exists, declared); the hold-out -- Zhai 2023 JAFC's 140 C pots in BOTH tables.
**Measurement model:** semi-quantitative series share one scale per figure (tau 0.3 dex, the Jiangnan thiols as
well as the quinoxalines); each Zhou paper's MGO its own scale (CF-060); calibrated and SIDA rows absolute; extra
noise per data set x observable from within-series scatter.

**Done-tests (declared now, raw residuals):** Zhai MGO |mean| <= 0.3 and RMS <= 0.4 per data set; SIDA / calibrated
MFT and FFT (Hofmann fed-precursor pots, Schieberle 2000, Zhou 2022) RMS <= 0.3 dex (x2) per data set; semi-
quantitative thiol and deoxyosone figures RMS <= 0.35 in shape; the pentose pH (Cao, Zhao; sub-model 1's table,
scored unfitted) within 0.3; every constant off its bound with a Laplace sigma under 1 dex; the hold-out (JAFC 140 C)
scored once, after the fit. Leave-one-paper-out on the winner.

**Expected (written now):** the joint fit keeps Zhai's MGO within tolerance AND brings Schieberle's MFT within x3
(P = 0.35); `k_dpo_ddp` comes off its bound once the H2S budget can move (P = 0.5); `k_h2s_loss` rises (Li 2022's
single point suggested ~10x faster at 90 C; P = 0.6 that it moves up by > x3).
## 5. Design review of the joint fit (step 5) -- done 2026-10-30, before any fit; sec. 4 is superseded where noted

An independent reviewer (read-only; flux budgets of the null engine and of sub-model 3's A2, closing within 1e-4)
found 7 BLOCKING, 6 SHOULD-FIX, 4 NOTE.

**The H2S explanation is half right.** In Schieberle 2000 the B59a Strecker share supplies 82 % of the H2S; A2 removes
it (supply / 4) and FFT / MFT fall -- as sec. 4 says. In Zhou 2022, Yiltirak and Zhai, A2 RAISES MFT: 1-DX no longer
drained by cysteine goes to norfuraneol (NF x3 to x30, MFT x2.6 to x47). The DDP route carries < 0.1 % of MFT
everywhere (`k_ddp_mft` 10^-6.5, frozen in B60).

**BLOCKING, and the response:**
1. FFT forms 89-99.5 % through `r_fur_fft_hs`, not `k_fur_fft` -> one common multiplier on both FFT branches.
2. `ch_mft_decay` destroys 80-99 % of the MFT formed (MFT ~ formation / k_mft_decay); CF-013: the shipped loss is
   >= 5x too fast -> `k_mft_decay` free (rate; Ea held), plus the CF-013-capped rival arm.
3. Declared hold-outs were fit rows (Hofmann T2 pH 3 / 7 and T1 xylose; Yiltirak; every Zhai 140 C pot, Food Chem
   AND JAFC, in both tables; the Hofmann T1 level rows are validation-only since B9) -> all held out. This also
   weakens sub-model 3's A2 / A3 "140 C hold-out": the Food Chem 140 C pots of the same vessels were fit rows there
   (recorded in that prereg).
4. A2's MGO destroys 2-acetylthiazole's credibility: Zhou 2022 ACTZ predicted ~20,600 ug/L against 11.7 measured
   (Zhou 2023 pH 7, B60 fit row) -- `k_cys_actz` was fitted at an MGO 10^3 too low -> ACTZ scored as a FALSIFIER of
   every arm (the row is not yet in a two-reader keyed table; freeing `k_cys_actz` waits for it).
5. Strecker H2S is attached only to 1-DX + Cys (B59a); moving 1-DX flux to MGO deletes H2S that should follow ->
   pre-registered as a topology caveat; an MGO + Cys Strecker-H2S arm is owed (needs a network step).
6. tau 0.3 for the Jiangnan thiols contradicts CR-15 / CF-065 (semi-quantitative = tau inf) -> thiols tau inf per
   figure (keeps the cross-pot contrasts); the deoxyosone quinoxalines' tau 0.3 gets its choices-register entry.
7. Ten rows sit in both tables (Hofmann T5 FUR / NF, Zhou 2022 FUR, Wang 2026b FUR) -> deduplicated.
SHOULD-FIX adopted: a declared between-run noise floor of 0.10 dex where a series has too few rows to estimate its
own (Garching sd 10 % would otherwise dominate; CF-016 / CF-055); Zhou thiols get the per-paper scale Zhou MGO has
(CF-014: 7.8x within-lab reproducibility); rows the coordinates cannot move are score-only (Kumazawa survival,
Hofmann T8 C2 + C3, T10 thiamine, T9 glucose, Zhang 2026); `k_dpo_ddp`'s bound widened to 10^-8; the B60 rows and
declared hold-outs re-scored before any acceptance. The B59a share is run at BOTH CR-11 alternatives (0.042, 0.080).
NOT adopted tonight: returning Hofmann T1 ribose pH 5 to the fit (an owner decision); the glycolaldehyde sink (no
data until Hofmann 1999 is read).

**Expected results, revised by the reviewer's flux arguments:** `k_h2s_loss` up by > x3: P ~ 0.25 (the fit must
replace lost Strecker H2S, which pushes it DOWN; a fall against Li 2022's "x10 faster" would flag a missing H2S
source); `k_dpo_ddp` off its bound: ~0.3; Schieberle MFT within x3 with `k_mft_decay` free: ~0.45. Expected failures
written now: 3-DX (no coordinate touches it), Hofmann T8, Hofmann pH 7 MFT (no HS- route for norfuraneol), CF-017's
temperature ordering (every Ea held).
## 6. Round 1 read, and the round-2 design (2026-10-30, before any round-2 fit)

**How round 1 must be read (a correction).** The first report's residual tables are raw, BEFORE each semi-
quantitative series' free scale, although its header said "after". For a scale-free series only the scatter around
its fitted offset counts. Read that way (s042, the first run), the Jiangnan thiol figures fit in shape (sd 0.1 - 0.3
dex: Zhai JAFC Fig. 3 MFT 0.23 / FFT 0.17 - 0.20; Food Chem Fig. 1a MFT 0.27, FFT 0.73; Wang 2026b 0.1 - 0.2; Zhou
2022 MFT 0.47, FFT 0.20), and the misses are in the ABSOLUTE rows: Garching's fed-precursor pots fit (Hofmann T3 /
T4 / T5 / T10 within 0.24), Schieberle 2000 (100 C time course) is under (MFT -0.85 at every time; FFT -1.2 to -3.1,
worst early), Yiltirak 2026 (hold-out) is over (MFT +2.35, FFT +1.9), and the Hofmann pH ladder (hold-out) has the
WRONG pH DIRECTION: measured MFT / FFT fall x22 / x19 from pH 3 to 7; predicted MFT is flat and FFT rises x180.
2-Acetylthiazole (falsifier) is x1000 high. The report generator now prints, per series, its scale class and the
residual after the fitted offset; hold-out rows of a scaled series take the offset fitted on that series' fit rows.
The four Feng 2022 pots (RG-ARP, precursor mapping unresolved, CF-058) were scored as constant 25-dex failures; they
are now skipped by the probe planner (the optimum is unaffected: no coordinate moves them). The s042 and capped runs
are repeated with both fixes so the three round-1 runs compare on one footing.

**Flux diagnosis at the s042 optimum** (`kinetic_core_thiol_flux_diagnosis_s042.md`; every budget closes to 1e-4):
1. H2S comes 100 % from cysteine thermolysis (`r_cys_h2s`, MEASURED, Zheng & Ho 1994, rising with pH). The B59a
   Strecker H2S is gone: `k_dpo_ddp` went to its 10^-8 bound (1-DX carbon routed to MGO), and `k_h2s_loss` fell to
   its lower bound -- the review's predicted signature of a missing H2S source (sec. 5, B5).
2. MFT is all `r_nf_mft` and sits at formation / decay (`ch_mft_decay` takes 92 - 99 %); FFT is all `r_fur_fft_hs`,
   lost to `ch_fft_decay` at pH 5 and to the thiolate loss at pH 7.
3. Furfural comes only from 3-DX (`r_tdp_fur`), and the model makes LESS of it at pH 3 than at pH 7 (0.060 vs 0.130
   mmol/L, Hofmann 145 C), while H2S rises x15 from pH 5 to 7: the inverted pH direction is the carbonyl supply and
   the thiol loss, not the measured H2S release.
4. 2-Acetylthiazole has one source, `r_cys_actz` (Cys + MGO): 0.089 mmol/L (~11 mg/L) in Zhou 2023's pot against
   11.7 ug/L. Its constant was fitted (B60) when MGO was 10^3 too low.
5. Schieberle 2000 vs Yiltirak 2026 (both SIDA, ribose + cysteine, 0.5 M phosphate, 100 C): measured MFT ~800 ug/L
   (Schieberle, 4 h, interpolated) vs 6.9 ug/L (Yiltirak, 4 h) -- x120, where 4x less ribose and pH 5.5 vs 5 explain
   about x10 and the model predicts Yiltirak HIGHER. Yiltirak's units and SIDA are verified by two readers (SI Table
   S3, ug/L; labelled standards added after heating; 3 mL under 17 mL air). A between-laboratory conflict, recorded
   as CF-068; the headspace oxygen (about 50 mM-equivalent O2 against 25 mM cysteine) is a candidate, not a finding:
   the shipped oxygen consumers are off (`k_cys_ox` = `k_red_ox` = 0) and no measured cysteine autoxidation rate is
   in the corpus.
6. Zhou 2023's final pH (ledger `zhou2023_fig2_final_ph`): ARP + Cys started at pH 6 / 7 / 8 ends at 3.22 / 3.42 /
   5.07; the lane ends the pH 7 pot at 4.54. The ARP + cysteine pots acidify ~1 pH unit more than the lane says.

**Round-2 arms (each on the s042 engine state, the same rows plus what is named; independent design review before
any of them runs):**
- **R2a, the H2S source follows the carbon:** `MGO_CYSTEINE_STRECKER` on, its rate a coordinate (log10 of k at
  100 C / pH 5, flat prior, bounds 10^-6 to 10^0 L/(mmol min); Ea held at 64.08), its H2S share CR-11's (0.042 /
  0.080 with the run). Expected: `k_h2s_loss` comes off its lower bound (P 0.5); Schieberle MFT within x3 (P 0.4);
  the pH direction NOT fixed (the step carries the "base" pH factor; P 0.8 it stays inverted).
- **R2b, 2-acetylthiazole fitted:** `k_cys_actz` a coordinate (log10 multiplier, -6 to +1), fitted on the Zhou 2023
  pH 7 row (absolute, the paper's per-paper scale not applied: one row cannot carry a scale; CF-014's x7.8 within-lab
  reproducibility as its noise floor, 0.9 dex); pH 6 and 8 SCORED, not fitted -- a test of the step's pH structure
  (measured x240 from pH 6 to 8). A leaf (ACTZ takes < 2 % of the cysteine anywhere), so R2b rides with R2a without
  confounding it. Expected: the pH 7 row within x3 after the fit (P 0.9); pH 6 and 8 within x3 unfitted (P 0.2).
- **R2c, the pH of the thiol layer (after the two-reader Zhou 2023 pH table exists):** fit rows Zhou 2023 MFT / FFT /
  MFTD / furfural at initial pH 6 / 7 / 8 (per-paper scale, tau inf: the pH contrast counts) and its final pH (pH
  rows: absolute); coordinates two `ph_slope`s, on furfural formation (`r_tdp_fur`) and on the thiol assembly steps
  (`r_nf_mft`, `r_fur_fft_hs` together), sign free, applied as in sub-model 1 (10^(slope (pH_in_situ - 5)));
  the Hofmann pH ladder stays a hold-out. Expected: both slopes negative (P 0.65); the hold-out ladder's direction
  right for FFT (P 0.45) and MFT (P 0.35).
- **Not this round:** the headspace-oxygen arm (no measured rate); the hexose -> C5 entry (Hofmann T1 / T2 / T5 /
  T9 glucose rows are out until it is designed); the glycolaldehyde sink (Hofmann 1999's alanine pots are not
  simulable, CF-067).

**Stop rule (owner, 2026-10-30: "you decide on the priorities").** If R2a - R2c leave the absolute MFT / FFT rows
outside x3 in more than one laboratory, the thiol layer ships as a declared limitation with its between-lab spread as
the error, and effort moves to the lanes not yet built from data (pyrazines, lipid aldehydes, Strecker aldehydes).

## 7. Design review of round 2 (step 5) -- done 2026-10-30, before any round-2 fit; sec. 6 is superseded where noted

An independent reviewer (read-only; it re-derived the sec. 6 numbers from the s042 JSON -- all confirmed -- and ran
the Zhou 2023 pH-series pots at the s042 optimum) found 3 BLOCKING, 3 SHOULD-FIX, 1 NOTE.

**BLOCKING, and the response:**
1. **R2c would fit on the wrong pH.** On the cooled scale the lane's final pH in Zhou 2023's pots is right at
   initial pH 6 and wrong above it: ARP + Cys 3.22 / 4.38 / 7.24 (lane) vs 3.23 / 3.43 / 5.05 (measured) at pH
   6 / 7 / 8; ARP alone 3.19 / 4.55 / 7.10 vs 3.98 / 4.10 / 4.96. The lane's furfural falls 2.66 dex from pH 6 to
   8 where the measurement falls 0.43, so a slope fitted on these rows would come out POSITIVE for the wrong reason,
   and the disulfide residual swings by 3.7 dex through thiolate dimerisation at the wrong pH. -> R2c is DEFERRED
   until the in-situ pH of unbuffered near-neutral ARP / cysteine pots is right (sub-model 1's acid supply), and the
   disulfide rows are out of it (they also sit near their calibration's zero point; table adjudication note 2).
2. **The scale key threw the contrasts away.** Round 1 keyed each semi-quantitative scale by ledger id, not by
   figure as sec. 5 (B6) declared: Zhai JAFC Fig. 3's TTCA and TTCA + Cys arms got separate scales (fitted offsets
   MFT +1.44 vs +2.67 -- a 1.2 dex error in the cysteine-addition contrast, hidden), single-row series carried
   nothing, and Zhou 2023's pH 7 column would have been scaled apart from pH 6 / 8. -> Fixed: paper x figure x
   species (Zhou's thiols and MGO: paper x species). All three round-1 runs are repeated with it (an implementation
   deviation from sec. 5 corrected, not a new choice); the first s042 report stays as the superseded record.
3. **R2b is not a leaf.** `r_cys_actz` takes 12 - 24 % of the methylglyoxal consumed (Zhou 2023 15.6 %, Zhou 2022
   14.8 %, Zhai TTCA + Cys 12.3 %, Hofmann pH 7 24.2 %) and shares R2a's MGO x Cys rate law. CF-014's x7.8 is MFT's
   reproducibility; 2-acetylthiazole itself reproduces to x1.4 between the two Zhou papers. -> floor 0.15 dex;
   R2b runs ALONE first, then R2a + R2b, and the report gives their correlation. "pH 7 within x3, P 0.9" is
   withdrawn as a test (one constant on one row); the test is pH 6 / 8, and it reads the pH model as much as the step.

**SHOULD-FIX, adopted:**
4. Sec. 6 flux claim 3 is wrong for FFT: from pH 3 to 7 (Hofmann, 145 C) FFT formation rises x1600 = furfural
   x2.2 x measured H2S x15 x ~50 from the HS- branch (its constant 43x the neutral one, tied to it by one multiplier
   in round 1); the furfural pH direction comes from the cysteine-Amadori decay (de Roos's curve extended above pH 6,
   CR-14), not `r_tdp_fur`. Every Garching ribose / xylose pot is under on FFT at acid pH while fed furfural + H2S
   (T3) fits. -> the replacement for R2c, when the pH is right: the HS- / neutral FFT ratio free (two multipliers)
   and the CR-14 alternatives, not one shared slope (Hofmann needs ~-0.34 dex / pH for MFT, ~-0.9 for FFT).
5. R2a's ceiling: its extra H2S is at most share / (1 - share) x the MGO consumed -- in Schieberle x1.7 (s042) or
   x2.4 (s080) the cysteine H2S, where MFT within x3 needs x2.3. -> expected "Schieberle MFT within x3" lowered to
   P 0.2; success is read on the Strecker share of Schieberle's H2S budget with `k_mft_decay` within 1 sigma of its
   round-1 value; `k_h2s_loss` is a flat direction (Laplace sigma inf), so "off its bound" is dropped as a test. The
   capped round-1 run is read before R2a runs.
6. The stop rule degenerates to "Garching fails" (the only absolute fit laboratory; Yiltirak is a hold-out in CF-068).
   -> rows that count: absolute MFT / FFT, fit rows and declared hold-outs separately; US 10,076,126 (Nestle, SIDA,
   110 C) scored at bracketing buffer strengths as a third laboratory; shape criteria as sec. 4.

**NOTE:** the Hofmann pH ladder is no longer blind for a pH structure chosen after seeing it; Meynier & Mottram 1995
(CF-012) is the candidate blind pH hold-out if it can be keyed. The MGO Strecker step's Ea (64.08) is the carbonyl
family's, not a measured value for this step. The diagnosis mislabelled Hofmann xylose pH 5 (a hold-out, off by
MFT +0.72 / FFT -0.68) as "the control that fits"; corrected. `ch_mft_decay` takes 92 - 99 % of MFT only at
pH <= 5.5 (47 % at Hofmann pH 7, 71 % in the Zhai pot).

**Order now:** (1) round 1 repeated with figure scales (s042, capped, s080); (2) R2b alone, then R2a + R2b, on the
round-1 winner; (3) in parallel, the in-situ pH of unbuffered near-neutral pots (sub-model 1), which gates R2c.

## Amendments
