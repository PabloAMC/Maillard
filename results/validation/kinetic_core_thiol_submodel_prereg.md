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

## 8. Round 1 results, scales per figure (2026-10-30; s042, s042_capped, s080)

The three runs agree: every constant within its Laplace sigma, reduced chi-square 5.38 / 5.40 / 5.32 (21.8 with the
ledger-keyed scales). `k_dpo_ddp` goes to its 10^-8 bound in all three, so the B59a Strecker share (0.042 vs 0.080)
is not identified -- the H2S it carries rides on a step the data switch off. Base for round 2: s042.

| constant (log10 x shipped) | s042 | capped | s080 |
|---|---|---|---|
| NF -> MGO share | 0.016 +- 0.005 | 0.015 | 0.016 |
| k_dpo_c2c3 | +2.20 +- 0.18 | +2.21 | +2.26 |
| k_dpo_nf | -0.19 +- 0.19 | -0.15 | -0.15 |
| k_nf_mft | -0.68 +- 0.15 | -0.61 | -0.70 |
| FFT formation (both branches) | -1.28 +- 0.26 | -1.25 | -1.28 |
| k_h2s_loss | -0.17 +- 0.43 | -0.11 | -0.21 |
| k_mft_decay | -0.92 +- 0.42 | -0.70 (cap) | -0.91 |

**Done-tests (sec. 4), s042:** Zhai MGO |mean| 0.09 / 0.16 / 0.26 and RMS <= 0.31 -- PASS. Absolute MFT / FFT within
0.3: Garching fed-precursor pots PASS (T3 0.19, T4 0.29, T10 0.09); Schieberle 2000 FAIL (MFT -1.58, FFT -1.72);
Zhou 2022 is per-paper scaled (shape MFT 0.03, FFT 0.06 -- PASS). Semi-quantitative shapes <= 0.35: PASS except Zhai
JAFC Fig. 3 MFT (TTCA arm -1.00, TTCA + Cys +0.62 after their now-shared scale: the lane's MFT answer to added
cysteine is ~1.6 dex too strong) and Food Chem Fig. 1a FFT (0.71). Constants off their bounds: all but `k_dpo_ddp`.
`k_mft_decay` 10^-0.9 (x8 slower than shipped) is CF-013's >= x5 from an independent route; `k_h2s_loss` is now
identified at the shipped value. **Hold-outs, scored once:** Hofmann pH ladder MFT +0.09 / 0.31 (PASS), FFT
+0.10 / 1.66 (the direction still inverted); Hofmann T1 MFT +0.24, FFT -0.52; Yiltirak MFT +1.41, FFT +1.90
(CF-068); Zhai 140 C, on the scale fitted at 100 - 120 C: FFT -1.44 / -1.47 (the lane's FFT rises too little with
temperature), MFT TTCA + Cys +1.10. 2-Acetylthiazole x1000 (+2.97), unchanged.

**Expected vs found (sec. 4-5):** `k_h2s_loss` up by > x3: no (P 0.25 stated); `k_dpo_ddp` off its bound: no (0.3);
Schieberle MFT within x3: no (0.45). Expected failures written in sec. 5 that occurred: 3-DX (no coordinate; now
within its tau 0.3 scale), Hofmann pH 7 FFT, the temperature ordering (every Ea held).

**Flux at the corrected optimum** (`kinetic_core_thiol_flux_diagnosis_s042.md`, regenerated): in Schieberle's pot
73 % of the H2S goes to `r_h2s_loss`, and MFT is late, not just low -- -2.4 dex at 30 min, -2.0 at 60, -1.2 at 6 h,
-0.7 at 12 h. A slope over time cannot come from a multiplier on one step: an early MFT source is missing or too
slow in a ribose + cysteine pot at 100 C (norfuraneol is the only MFT parent the lane uses; Hofmann T10's
route split and the intact-C5 DDP route, frozen at 10^-6.5, are the candidates). This is a new structural lead for
round 3; it does not change round 2's order (R2b, then R2a + R2b).

## 9. Round 2 read, and the round-3 design: MFT through the intact C5 skeleton (2026-10-30, before any round-3 fit)

**Round 2.** R2b (`k_cys_actz` free on Zhou 2023 pH 7): x10^-3.06 +- 0.39, nothing else moved (every round-1
constant within 0.03 dex) -- accepted as a record; its pH 6 / 8 test fails (about +0.6 / -1.5 dex: the step's pH
dependence is too weak). R2a + R2b: the MGO + cysteine Strecker rate goes to its lower bound (10^-6): the data reject
the step, as the review's ceiling argument (sec. 7, item 5) predicted. Reduced chi-square 5.48 / 5.33. Owner,
2026-10-30: "take the wisest decision to make the model sound and as accurate as possible (mechanism wise)".

**Why round 3 is this.** The lane makes 100 % of its MFT from free norfuraneol + H2S (`r_nf_mft`); the intact-C5 route
(1-deoxypentosone + cysteine -> 1,4-dideoxypentosone, DDP -> + H2S -> MFT) is frozen at k_ddp_mft 2.9e-7, ~7500x
below k_nf_mft, because B60 fitted it there and every round since held it fixed (sec. 4). The labelling evidence says
the opposite. Cerny & Davidek 2003 (ribose + cysteine, 0.5 M K phosphate, pH 5.00, 95 C, 4 h): MFT keeps the intact
ribose skeleton (95 % unlabelled + 13C5, the 1-4 label channels at the isotopic-impurity level), and in pot D, where
unlabelled norfuraneol and [13C5]ribose compete at equal charge for the same cysteine, only 7 % of the MFT is
unlabelled: free norfuraneol supplies ~7 %, the ribose ~93 % by another route ("only to a small extent from
4-hydroxy-5-methyl-3(2H)-furanone"). Cerny 2004 adds that the C2 + C3 route (hydroxyacetaldehyde + mercapto-
propanone) "was not relevant" at 95 C. A route through accumulating norfuraneol must start slowly, which is the
shape of the round-1 miss in Schieberle 2000 (-2.4 dex at 30 min, -0.7 at 12 h; sec. 8). The hold-out
declaration already names Cerny 2003's ceiling a HOLD-OUT to be "evaluated at Cerny's conditions".

**Arm R3.** On the engine state round 2 and arm E leave (E's entry coordinates if E passes its done-tests, else
R2b's state): one new coordinate, a log10 multiplier on BOTH DDP -> MFT branches (`k_ddp_mft`, `k_ddp_mft_hs`;
-2 to +7, start +3), with `k_dpo_ddp` (already a coordinate, on its 10^-8 bound while DDP was a dead end) free to
come back. Rows unchanged.

**Hold-out (scored once, after the fit): Cerny 2003 pot D.** Norfuraneol 72.5 umol + ribose 72.5 umol + cysteine
48.3 umol in a ~0.5 g charge (nominal 145 / 145 / 97 mmol/L: "dissolved to a total amount of 500 mg"), 0.5 M K
phosphate pH 5.00, 95 C, 4 h. Observation operator: the share of MFT formed through `r_nf_mft` (integrated flux) --
an upper bound on the unlabelled share, since in-situ norfuraneol from the labelled ribose also feeds that step
(reported alongside: the fed share of the norfuraneol pool). Pass: <= 0.15 (the measured 0.07 plus the nominal-
concentration and approximation slack). Also reported: the C2 + C3 route's share of MFT (`r_ha_mp_mft`) in Cerny's
ribose + cysteine pot (Cerny 2004: "not relevant"), and Schieberle's MFT at 30 / 60 min (the lag).

**Done-tests:** sec. 4's, plus the hold-out above; leave-one-paper-out on R3's optimum before any acceptance.
**Expected (written now):** the DDP multiplier moves up by >= 10^3 (P 0.6); `k_dpo_ddp` comes off its bound (P 0.55);
Cerny pot D's norfuraneol share <= 0.15 (P 0.45); Schieberle MFT within x3 at every time (P 0.3); `k_nf_mft` falls
(P 0.6). If the fit keeps the DDP route near zero, the evidence and the fit rows disagree about MFT's mechanism, and
the next step is the rows that pin `k_nf_mft` (Hofmann T4 / T10, fed norfuraneol at 145 C), not a new coordinate.

**Dry check (no fit, R2b's optimum):** the hold-out scorer reports pot D's norfuraneol share 1.000 at the R2b
optimum, and still 1.000 with the DDP multiplier at +3 -- because `k_dpo_ddp` sits at 10^-8 there and no DDP forms. With
`k_dpo_ddp` at its shipped value: share 0.98 / 0.35 at DDP x10^3 / x10^5; with `k_dpo_ddp` at 10^-1: 0.65 at x10^5.
The two steps only matter together (a flat start would strand the optimiser), so R3 starts at `k_dpo_ddp` 10^-1
(sec. 4's start) and the DDP multiplier at +4.5.

## 10. Design review of round 3 (step 5) -- done 2026-10-30, before any round-3 fit; sec. 9 superseded where noted

An independent reviewer (one diagnostic with an exact fed-norfuraneol tracer at R2b's optimum) found 2 BLOCKING,
4 SHOULD-FIX, 2 NOTE. Judged under the acceptance protocol adopted the same day (STRATEGY sec. 9: mechanism
signatures lead; done-tests diagnose).

**BLOCKING, and the response:**
1. **The DDP-forming step's pH law forbids a pass.** `r_dpo_ddp` carries the engine's "base" (free-amine) factor:
   x0.29 / 1.43 / 143 / 1401 at pH 4.3 / 5 / 7 / 8 against pH 5 (95 - 100 C). Cerny's pot ends at pH 4.3, Zhai
   2021's 5.5 - 8 ladder (47 DPO + 48 MGO rows) holds `k_dpo_ddp` at 10^-8, and pot D needs ~10^-1. The law is an
   unmeasured assumption (Cerny's route is a proposal; DDP was never detected). -> R3 becomes a contest of the
   step's pH law, rival structures judged by the protocol: "base" (shipped), "none", and a free slope
   (`parameters_sulfur.DDP_FORMATION_PH`; arms R3 / R3n / R3s, and with the entry R3E / R3En / R3Es). The fallback
   in sec. 9 ("the rows that pin k_nf_mft") is withdrawn; the expected pass for R3 as designed is cut to ~0.15.
2. **One coordinate moved a ridge.** Pot D's share depends on the product of `k_dpo_ddp` and the DDP -> MFT rate
   until DDP -> MFT saturates; (-1, +4.5) and (-3, +6.5) give the same share. -> the coordinate is now the route
   PRODUCT (`offset_ddp_route`, start +5.5) with `k_dpo_ddp` (start -1); a second start at the ridge's other end
   (-3, +3.5) for the winner; the H2S budget at the optimum is reported, with B59a's share at 0 as a sensitivity.

**SHOULD-FIX, adopted:**
3. The observation operator is now the exact tracer (the fed norfuraneol as a pool never replenished: its fraction f
   of the pool obeys d ln f / dt = -production / NF; MFT from fed norfuraneol = integral of r_nf_mft x f over all
   MFT formation), scored two-sided in [0.03, 0.15] (measured 0.07). The flux share stays as a check (it was 0.087
   vs the tracer's 0.085 at R3's start). The "fed share" bug (NF[-1] / NF[0]) is gone.
4. Pot D was seen at design (the arm was chosen from it and dry-checked on it): recorded in
   `FIT_HOLDOUT_DECLARATION.md`. Blind companion added: pot C's 2-mercapto-3-pentanone (Cerny: none from ribose;
   in the lane only from norfuraneol): MP3P / MFT made <= 0.1.
5. **Correction of the record (secs. 4 and 9):** B60 did not fit `k_ddp_mft` -- its report lists it frozen. The
   10^-6.54 is B2.3's value, and the constant wandered across waves (B2 -4.62, B2.1 -2.35, B2.2 -5.83, B2.3 -6.54),
   flagged UNIDENTIFIED each time (B38: sigma 800 dex): a numerical leftover, not evidence against Cerny.
6. R3 and R3E both run whatever arm E's verdict (E's thiol rows could fail for the same missing route); E itself is
   judged on its own mechanism signatures (acetic-acid plateau, xylose decay, acetic pH effect), not on the thiols.

**NOTE (declared limitations):** the two DDP -> MFT branches keep B2.3's unidentified HS- / neutral ratio (0.07,
against the lane's own note that HS- is the better nucleophile and FFT's 43); a sensitivity at 43 follows the
winner. Every formation step shares Ea 64.08 and the DDP and norfuraneol decays share 143.5, so the lane's route mix
barely moves with temperature: a pass at 95 C carries to 145 C, where Hofmann T10 measures the C2 + C3 route as the
strongest -- the temperature dependence of the route mix is not representable yet, and is recorded as such.

**Dry check (no fit, R2b's optimum):** at (`k_dpo_ddp` 10^-1, route +5.5) the tracer share is 0.085 ("base") and
0.042 ("none"), pot C MP3P / MFT 0.002 / 0.0005; at the ridge's far end (10^-3, +3.5) 0.81 / 0.58. What decides
between the laws is whether Zhai's pH ladder lets `k_dpo_ddp` sit near 10^-1 -- the question the rival fits answer.

## 11. Draft (2026-10-31, before round 3's result): the MFT route mix vs temperature (plan step B6)

*A draft, written while round 3 runs; it is pre-registered for review only after round 3 decides the DDP step's pH
law, because both arms move the same routes.*

**Why.** The lane gives every formation and assembly step one lumped barrier (64.08 kJ/mol) and the DDP and
norfuraneol decays one (143.5), so its route mix barely moves with temperature (sec. 10, NOTE). The measured route
mix does move: Cerny 2003 (95 C, pH 5) -- MFT keeps the intact skeleton, ~7 % via free norfuraneol, fragment
recombination <= 5 %; Hofmann 1998 T10 (145 C) -- norfuraneol and the C2 + C3 route (hydroxyacetaldehyde +
mercapto-2-propanone) strong; Liu 2023 (LWT, 168 C, pH 5.88, carbon-module labelling) -- ~81 % intact C5 with
norfuraneol named as the intermediate, ~19 % recombined fragments.
**Measured barrier available:** Lau 2003 -- norfuraneol formation from ribose + protein amines, first order,
Ea 81.5 +- 3.2 kJ/mol over 116 - 131 C (whey gel; HPLC peak areas, so a rate, not a yield). A barrier above the DDP
route's makes norfuraneol gain share with temperature, the measured direction.
**Arm (B6):** the norfuraneol-formation barrier at Lau's 81.5 (measured-class, band +-10, a declared ribose-gel ->
aqueous transfer); the DDP route's barrier free (its pH law from round 3); the C2 + C3 route (`r_ha_mp_mft`) with its
own barrier free. **Signatures (declared now):** Cerny 95 C fed-norfuraneol share in [0.03, 0.15]; Liu 168 C
fragment (C2 + C3) share of MFT 0.19 (band 0.10 - 0.30, simulated at Liu's conditions -- volume from the paper, to be
read); Hofmann T10 145 C route ranking; Zhai 140 / 120 C and Yiltirak's temperature ladder directions.
**Open question for the literature (owner downloads on DOIs):** a labelling study at two or more temperatures in one
laboratory -- the search of 2026-10-31 found none (P ~0.1 that one exists).

## 12. Round 3 read so far, and arm R3T (round 3 + temperature-dependent routes) -- design, 2026-10-31, before any R3T fit

**Round 3 so far** (experiment register EX-015, EX-016): with the shipped free-amine law (R3E) the fit switches the
DDP route off; with no pH law (R3En) the data ACCEPT DDP formation (`k_dpo_ddp` 10^-0.61 +- 0.18, reduced chi-square
4.95, the lowest so far) but raise norfuraneol + H2S -> MFT x26 (`k_nf_mft` +1.41), so Cerny 2003's pot D share is
0.994 (measured 0.07). Both fail for one reason: `k_nf_mft` is pinned at 145 C by fed-norfuraneol rows (Hofmann T4 /
T10) while Cerny's 95 C labelling wants the step minor, and one shared barrier (64.08) cannot do both. R3Es (free
slope) is running; R3T takes whichever law the signatures favour (R3En if they tie: fewer coordinates).

**Arm R3T.** R3En's (or R3Es's) state plus ONE coordinate: the barrier of norfuraneol + H2S -> MFT (`ea_k_nf_mft`, kind
Ea, 40 - 200 kJ/mol, start 100). `k_nf_mft`'s reference temperature is 145 C, so its multiplier stays pinned by the
145 C rows and the barrier alone sets its 95 - 100 C value -- the decoupling the conflict needs. Sensitivity R3TnL:
norfuraneol formation (`k_dpo_nf`) at Lau 2003's measured 81.5 kJ/mol (ribose + protein amines, 116 - 131 C; a declared
transfer from a whey gel to water). Signatures and hold-outs as sec. 10, plus the temperature signatures already in
the scorer (Zhai 140 / 120 C for MFT and FFT; Schieberle's lag).

**Dry check (R3En's optimum, the barrier alone moved, no fit):** pot D share 0.994 / 0.991 / 0.971 / 0.877 at
64 / 100 / 140 / 180 kJ/mol. The fed norfuraneol (145 mM) outweighs the DDP made in situ at 95 C, so lowering the
norfuraneol route is not enough: a pass also needs the DDP route to rise at 95 C (the route product is a coordinate;
the DDP route's own barrier is NOT freed here, to keep the arm at one new coordinate).
**Expected (written now):** `ea_k_nf_mft` > 100 kJ/mol (P 0.6); Cerny pot D share in [0.03, 0.15] (P 0.2 -- lowered by
the dry check); Schieberle's MFT lag signature within tolerance (P 0.35); signatures >= 7 / 15 (P 0.3). If the pot D
share stays above ~0.5 with the barrier on its upper half, the next question is the DDP route's own barrier and its
supply at 95 C (Cerny's pots), not more freedom on the norfuraneol side.

## 13. Design review of R3T, and the premortem -- 2026-10-31, before any R3T fit: R3T is NOT run as designed

**Blocking (independent design review).** A grid at R3En's optimum (`ea_k_nf_mft` x `offset_ddp_route`, no refit):
Cerny pot D's share falls into [0.03, 0.15] only where the barrier sits at or beyond the measured sulfur class's upper
end (30 - 133 kJ/mol, median 100.9) AND the DDP route is raised so far that Schieberle 2000's early MFT is over-made by
~1 dex (Ea 160, DDP +3: share 0.038, Schieberle 30 / 60 min +1.30 / +1.17 dex). The conflict that binds is not 95 C vs
145 C but **Cerny 95 C vs Schieberle 100 C** -- the same chemistry (ribose + cysteine, pH 5, phosphate) 5 K apart, which
no Arrhenius term can separate: the DDP route removes Schieberle's measured MFT lag, the norfuraneol route supplies it.
Sec. 12's diagnosis (pinned at 145 C by T4 / T10) missed the binding row. Also (should-fix): the barrier is not
identified on its own (the temperature ladders are all hold-outs); the arm did not nest R3En (it started from START, not
R3En's theta); Lau 2003's 81.5 kJ/mol is an overall ribose -> norfuraneol appearance barrier, not `k_dpo_nf`'s (R3TnL
would be a declared assumption, and moving one DPO branch alone contradicts Lau's temperature-independent plateau);
the Zhai 140 / 120 C MFT signature was read at design (R3En: -0.40 vs measured +0.20, wrong sign) and is recorded as
seen; the pot C MP3P companion cannot fail (the lane's NF -> MP3P : NF -> MFT is 0.003 - 0.075, Whitfield 1999 measures
~5 at 140 C).

**Premortem (same night).** Every joint-fit arm so far stopped on its evaluation budget (11 optimiser iterations), so
EX-015 / EX-016 / EX-017 are inconclusive until converged (STATUS_AND_PLAN sec. 5, B0: E converged, then each rival
from its own optimum under E's fixed measurement model; queued 2026-10-31 ~00:50).

**Decision.** (1) No R3T fit. (2) The question becomes mechanistic: **what makes MFT LATE at 100 C in ribose +
cysteine when the isotope labels say the fed norfuraneol is a minor source?** Rival structures to be pre-registered
after a literature pass (candidates: a slow H2S supply -- cysteine's H2S via its own Strecker degradation, which waits
for dicarbonyls; a ribose-cysteine thiazolidine reservoir delaying both precursors; a two-step intact-C5 route whose
intermediate builds slowly; in-situ norfuraneol reacting before it is released as the free furanone, which fed
norfuraneol cannot). (3) Any later barrier coordinate is bounded to its measured class (<= 133 kJ/mol for sulfur
steps) and starts from the arm it extends. (4) Backlog: tie NF -> MP3P to NF -> MFT at Whitfield 1999's ratio, so
Cerny's "no MP3P from ribose" becomes a real test.

## 14. The 100 C MFT lag: literature pass and dry check T1 -- 2026-10-31, before any converged arm was read

**What the data say (literature pass; dossiers schieberle2000, cerny2003, deroos2005, nedvidek1992, zhai2019, cao2017).**
Schieberle's MFT is not a slow-accumulation curve: MFT = A (t - tau)^n with tau ~ 9 min (the autoclave's come-up) and
n ~ 1.25 on both intervals, then a STOP (720 / 360 min = 1.15) -- formation roughly constant for ~6 h, ending when a
reservoir runs dry (cysteine, 33 mM, is the limiting reagent). An accumulating-norfuraneol route gives the wrong shape
at both ends. Cerny's pot D: 96 % of 2-mercapto-3-pentanone comes from the fed NF but only 7 % of MFT, and pot C (no
fed NF) makes no MP3P -- so fed NF does meet H2S at 95 C; NF -> MFT is slow because the STEP is slow there. de Roos
2005: during boiling NF gives no more MFT than ribose (the NF route needs >= 140 C, closed); the cysteine Amadori
compound, decomposing first order, is the precursor. Nedvidek 1992: 1,4-dideoxyosones form only with an alpha-amino
acid (a reduction coupled to Strecker degradation).

**Ranked mechanisms (agent's probabilities, adopted as the priors of the next design).** (1) An intact-C5 route fed by a
reservoir -- TTCA -> cysteine Amadori compound -> 1-deoxypentosone -> DDP (Strecker-coupled) -> + H2S -> MFT -- rate set
by the reservoir's drain (hours), plateau at cysteine exhaustion; NF -> MFT with a high barrier, minor <= 100 C: 0.45.
(2) Slow H2S supply as a co-limiter: 0.20 (against it as the sole cause: MFT, FFT and MP3P have different shapes in
the same pot). (3) Start conditions / heat-up as the main fix: 0.10. (4) Pot differences (Cerny: 3x the cysteine,
H2S in headspace, pH 4.3 at the end): 0.08. (5) A slowly building two-step intermediate: 0.08 (n ~ 1.25 argues
against). (6) Bound in-situ NF: 0.05. (7) Oxidation / measurement artefacts: 0.04.

**Signatures added (list v2, before any converged arm):** Schieberle MFT plateau 720 / 360 (tol 0.25, group "MFT
supply exhaustion") and early rise 60 / 30 (tol 0.25, sign). At R3En's optimum: plateau -0.20 vs +0.06 (the engine's
MFT falls after 6 h), early rise +0.96 vs +0.49.

**Dry check T1 (no fit; `scripts/generators/diagnose_mft_start_conditions.py`).** Schieberle's pot under four start conditions -- as fitted
(isothermal, all cysteine free); a 13-min linear come-up 20 -> 100 C; TTCA pre-equilibrated at Cao's K (26.9 of 33 mM
cysteine bound); both -- at R3En's optimum and at sec. 13's Cerny-passing point (Ea 160, DDP +3). MFT residuals (dex)
at 30 / 60 / 360 / 720 min:

| point | as fitted | ramp + TTCA |
|---|---|---|
| R3En optimum | -0.65 / -0.17 / +0.16 / -0.10 | -1.15 / -0.37 / +0.15 / -0.10 |
| Cerny-passing | +1.30 / +1.17 / +0.03 / -0.12 | +0.76 / +1.10 / +0.03 / -0.12 |

Mechanisms (1) and (3) predicted the early excess at the Cerny-passing point would fall to <= +0.5 dex: it falls at 30
min (+0.76) but NOT at 60 min (+1.10). So start conditions are worth ~0.5 dex at 30 min (mechanism 3 is real but not
the fix) and **the early DDP supply itself is ~10x too fast in the first hour**: in the engine the intact-C5 route
is supplied by a pool that is full at once, where the data want it fed at the reservoir's drain rate. Both points also
put MFT FALLING after 6 h (720 / 360 err -0.15 to -0.26), against the measured stop.

**Next design (to be pre-registered and reviewed in the day, not run tonight).** Rival R4a: DDP formation fed only by
the cysteine Amadori compound's first-order decomposition (de Roos 2005; Hou 2017's fed cysteine-Amadori pots) with the
TTCA reservoir and the come-up in place; R4b: the same plus an H2S-supply limit; the null: R3En's structure with the
come-up. Decided by the v2 signatures (MFT timing, supply exhaustion), Cerny pot D, and the converged R3En as the
baseline. A designed time-resolved pot-D experiment (fed NF + [13C5]ribose, 0.5 - 6 h at 100 C) separates (1) from
(2): a candidate for the designed time-course protocol.

## Amendments
