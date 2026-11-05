# Case M4: glucose to MFT; protein and other amines as classes

*Plan step M4 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E1 (F23-a, F24-a) and E8 (F25-a,
-b, F06-b) of [ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in sec. 5 (E1: "engine
error SUPPORTED; band-now fix REFUTED"), and the plan item "amine as a variable" (F02-b, F05-a, F19-e). Items 1 - 5 of
sec. 7, written 2026-11-04 before any result of the change exists. No engine was run for this case. Items 6 - 9 are
appended later and never edited into this part. M4 is built on the one network (W3,
[A2_MODULE_CONTRACT](../../guides/A2_MODULE_CONTRACT.md)); sec. 6 says what it needs from it.*

**Seen by the author.** Engine outputs read: ENGINE_VS_LITERATURE Appendix A probes C and E (glucose + Cys MFT
1.5e-14 mM, 100 % C2 + C3; ribose / glucose MFT 4.4e10, FFT 0.47) and B (FFT vs pH); the sec. 5 red-team probe of
`HEXOSE_C5_ENTRY` (ribose / glucose MFT 45 - 78 without amine, 0.79 with 33 mM glycine; a scratchpad probe, not an
artifact); `core_panel_scores.md` rows for the three Hofmann glucose bundles (all "refused") and
`thiamine_cys_glucose_120C_Bolton1994` (predicted 232.9 vs 11.7 ppb). No engine output exists on any v2 E1 or E8 TEST
contrast (both sets are empty). Measured values only for the v1 TEST cells F06#16 (Lee 2022), F14#5 (Hollnagel 2000)
and the v2 E5 TEST cell F06#11 (Thornalley 1999), read in the atlas files.

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml). v2 E1: DEV labs BTBU, CAAS, DFA Garching,
Nestle / Cerny; **no TEST contrast** (UNTESTABLE_on_held_out). v2 E8: DEV labs DFA, Reading, Wageningen; **no TEST
contrast**; Cerny 2007's F13#7, F25#6, F25#7, F25#22 excluded as calibration input. Atlas cells are 0-based indices of
the F files in [pathway_atlas](../pathway_atlas/). Every cell below is v1 DEVELOPMENT unless marked. Two cautions:
- F23#7, F23#8, F24#4, F24#5 and F25#23 (Zhao 2019b, Du 2023) are DEV(E1) but **TEST for E6 and E11 in v2** (M2's
  errors), and their keyed tables (`heldout_du2023_tab2.yml`, `heldout_zhao2019b_tab1_2.yml`) are TEST-only
  ([contrast-benchmark README](../contrast_benchmark/README.md) sec. 11). This case reads only their route pattern
  (intact C5 vs mixed isotopomers); anything they say about the norfuraneol share or the route mix with temperature is
  M2's held-out evidence and is read in item 7 only.
- The Hofmann glucose pH 3 / 7 arms are DEV under both frozen splits, and separately declared Maillard-path hold-out
  bundles (`data/benchmarks/external_validation/maillard_path/mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH{3,7}.json`);
  the pH 5 arm is a trust-loop predictive row, not a B60 fit row (`kinetic_core_b60_fit_targets.json` holds only the
  ribose and xylose T1 rows). See open question 2.

## 1. The claim

Three parts, built and judged separately on the one network.

| part | step(s) | today (B60, read from the code) | M4 |
|---|---|---|---|
| G. hexose entry to MFT | new: Glc + Cys -> cysteine-hexose Amadori compound (D-fructose-L-cysteine) -> 1-DG + Cys released; then `r_odg_dpo` (1-DG -> 1-DPO + C1, [sulfur.py](../../../src/kinetic_core/sulfur.py) ~927) into the shared pentose MFT steps | glucose reaches MFT only through `r_glc_c2c3` (log10 k -8.59 at 145 C) and FFT through `r_glc_fur` (`acid` factor); `HEXOSE_C5_ENTRY.k_per_min_100c = 0` ([parameters_sulfur.py](../../../src/kinetic_core/parameters_sulfur.py) ~623 - 641); cysteine is not an amine on the trunk, so a Glc + Cys pot has no amine route to 1-DG; `engine.unidentified_routes` refuses hexose-only thiol requests (engine.py ~1024 - 1043) | a measured-class Cys-hexose Amadori formation and decay, `r_odg_dpo` on with a declared barrier |
| E. the furfural supply's pH law (E8) | `r_pent_tdp`, `r_arp_tdp`, `r_glc_fur` (`acid` = [H+], a decade per unit; sulfur.py ~1160 - 1175) | as stated | an amine-catalysed 1,2-enolisation term, pH-flat or weakly rising, beside a weaker acid limb |
| A. amine classes | the trunk's `r_schiff` (Glc + Gly), the sulfur lane's Cys-catalysed `r_pent_dpo` / `r_pent_tdp`, B20's `r_glc_lysp` (protein lysine, trunk only) | one amine: MET and PRO folded into `Gly` (engine.py ~2560 - 2570); `LYSP` exists only on the trunk lane; no peptide | classes: free alpha-amino acid, cysteine, peptide, protein-bound lysine; each with its own condensation constant; their competition for each sugar declared |

**Part E depends on M1.** Enabling `CYS_AMADORI_ROUTE` (M1's combination C1) retires `r_pent_dpo` / `r_pent_tdp`
(`inactive_reactions`, sulfur.py ~975 - 978) and carries the pentose supply on de Roos 2005 Fig. 4's measured,
U-shaped Cys-Amadori decay. If C1 ships, E8's pentose part is decided by M1 / M1b, and part E keeps only
`r_arp_tdp` (fed ARP pots) and `r_glc_fur`.

**Earlier waves on these structures** (cited from the record, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| hexose MFT through `k_glc_ha` (C2 + C3), carried by the eight Hofmann T1 level rows | B8 | the route existed only through the rows scored on it | whole-model fit |
| the same with the T1 level rows dropped | B9 (`tasks/data_restructure_plan.md`, 2026-09-03) | glucose and fructose MFT predict zero; `k_glc_ha` at its band floor; B9's Laplace covariance singular in that coordinate | whole-model fit |
| the hexose route reopened | B47, B51.1 | DO NOT SHIP: `k_glc_ha` +3.6 decades bought back one Whitfield row; directional claim SUG-09 flips | whole-model fit |
| hexose -> C5 entry (`HEXOSE_C5_ENTRY`) | sub-model 4 (backlog, unchecked) | never pre-registered or run | -- |
| the entry switched on alone (the "band-now" fix) | red team, ENGINE_VS_LITERATURE sec. 5 | ribose / glucose MFT 45 - 78 with no amine (measured 10.4), 0.79 with 33 mM glycine | read-only probe |
| protein-bound lysine glycation | B20 (SHIPS, trunk only) | Nguyen 2016's constants; answers only with `protein_g_per_l` | fit on printed rates |
| pentose + glycine entry | EX-012 / EX-013 | admitted on audit; norfuraneol sanity bound FAILED (25 - 36 mol %) | fit |
| the furfural supply's pH law | Appendix A probe B | FFT x44 from pH 4.5 to 6 and ~x400 to 6.5 against x6.1 measured (Meynier) | probe |

**Errors stack (CF-107).** The B8 / B9 / B47 / B51.1 outcomes are whole-model fits and weak evidence about the hexose
route. Candidate compensating errors:
- E6 (MFT through free norfuraneol, 99.9 %) and E11 (one 64.08 kJ/mol formation barrier);
- the frozen pot-independent `k_mft_decay` and the Strecker sulfur sent to an inert pool (SULFUR_FLUX_MAP sec. 6);
- MFT destroyed in excess H2S at B9's optimum (backlog, "next sulfur question");
- E2 (the enolisation split frozen in pH) and the trunk's MGO level (B21).

The probe's 45 - 78 is an endpoint ratio that inherits E6 and E11, so it does not reject the entry either. Only the
labelling (sec. 2.1) is a measured share that makes a route inert: it rejects C2 + C3 as the hexose route.

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Du 2023 (F23#8, F24#5; ledger `du2023_table2_camola_isotopomers`) | Cys 1 mmol + 1:1 [13C6]/[12C6]-glucose, 0.2 M phosphate pH 6.5, 150 C, 120 min: MFT M+0 : M+5 = 52 : 48 (57 : 43); no mixed channel printed | hexose MFT is an intact glucose C5; C2 + C3 absent on the label basis | DEV(E1); TEST(E6, E11) |
| Zhao 2019b (F23#7, F24#4, F25#23; `zhao2019b_table2_camola_glucose_labelling`) | GSH + glucose 60 / 60 mM, pH 6.5, 140 C, 5 h: MFT ~100 % glucose C5; FFT ~94 % intact C5 | the same in a peptide-amine pot | DEV(E1); TEST(E6, E11) |
| Liu 2023b (F23#6, F24#6; CF-086) | ribose 168 C: fragments ~6 % after the [M-H]+ correction | fragments stay minor at high T (pentose) | DEV; TEST(E6, E11) |
| Cerny 2003 / 2004, 2007b (F23#0 - #3, F24#0 - #3) | pentose MFT intact C5; free NF 7 % (pot D) | the downstream steps M4 reuses | DEV |
| Bolton 1993 (F27#6; `bolton1993_tab2_mft_34s_split`) | 34S-Cys + thiamine in a glucose bouillon, a_w 0.83, 120 C: Cys supplies ~8 % of MFT sulfur; no MFT without thiamine (< ~2.3 % of the thiamine pots) | an upper bound on Cys + glucose MFT at 120 C in a matrix | DEV |
| Davidek 2005, 2006a (F07#5 - #7) | fed DFG and 1-DG: acetic acid ~70 % from glucose C1 / C2 | 1-DG's own exit (beta-cleavage), the competitor of any 1-DG -> C5 step | DEV |
| none | no positional glucose label (which carbon MFT loses) and no time-resolved hexose label on disk | which C5 (1-DPO vs a 1,4-dideoxy C5) | gap (ENGINE_VS_LITERATURE sec. 3.1, E1) |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Tables 6 / 8 (F24#7, F24#8) | fed C2 + C3 + H2S, 145 C: MFT 0.24 mol % at pH 5; 0.01 / 0.23 / 0.27 at pH 3 / 5 / 7 (rising with pH) | the C2 + C3 route's own pH shape (rival R3) | DEV; B60 fit rows (`hofmann1998_c2c3_recombination_*`) |
| Hofmann 1998 Tables 3 / 4 (F13#2, F23#12, F25#4) | fed NF / furfural + H2S, pH 5 | downstream steps shared with ribose | DEV; fit rows |
| Deng 2021 (F03#11 - #12, F06#6, F07#11) | fed Met-glucose ARP, 100 C: k 0.0048 (pH 5.5), 0.012 /min (pH 7.5); Met returned ~mole for mole | a hexose amino-acid ARP decay and its pH slope (the class analogue for the Cys-hexose ARP) | DEV |
| Martins 2003b (F03#3 - #5, F06#5, F07#10); Davidek 2002, 2005 (F03#0 - #2, F02#11) | fed DFG: amine returned 71 - 96 %; split moves with pH | the hexose ARP steps the claim reuses | DEV; Martins a calibration input |
| Davidek 2006a (F07#5, #9) | fed amine-free 1-DG: acetic 0.20 / 0.56 / 0.80 mol/mol at pH 5 / 6 / >= 6.5; 1-DG gone within 5 min at pH >= 6 | 1-DG's lifetime and its main exit | DEV |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Beck 1988, Nedvidek 1992 (F06#0 - #2, F07#0 - #3) | in-situ OPD on hexose ARPs / glucose + beta-Ala: 1-DG : 3-DG rises x12.5 from pH 4.5 to 7; a secondary-amine ARP goes to 1-DG even at pH 4.5 | the hexose 2,3-enolisation's pH law; amine structure moves the split | DEV (one school) |
| Nedvidek 1992 (F23#16) | the 1,4-dideoxypentosone quinoxaline needs alpha-amino acids | a 1,4-dideoxy C5 made by amino-acid Strecker on a deoxyosone (rival R2) | DEV |
| Gobert 2009 (F06#3, F07#4) | 50 C: air lowers 1-DG x3.4 | O2 on the 1-DG supply (low T only) | DEV |
| Hollnagel 2000 (F14#5) | dry glucose + OPD: deoxyhexosuloses < 1 vs 11 mol % with glycine | amine catalysis of enolisation | **v1 TEST**: listed, not used |
| Thornalley 1999 (F06#11) | BOC-Lys raises 3-DG formation x6 at 37 C | amine catalysis | **v2 TEST(E5)**: listed, not used |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 5 (`hofmann1998_t5_fa_nf_pairs`; dossier sec. 2.5) | Cys 33 + sugar 100 mM, pH 5, 145 C: NF / MFT 2754 (ribose), 7.1 (glucose), 23.9 (rhamnose); furfural / FFT 5.58 / 2.46 / 7.75 | whether glucose MFT passes through the intermediate ribose MFT does (D1) | DEV; ribose NF / furfural are fit rows, the glucose cells are not |
| Hofmann 1998 Table 9 (`hofmann1998_t9_m2p_glucose_cys`) | 1-mercapto-2-propanone from glucose + Cys 11.4 / 59.6 / 26.5 ug at pH 3 / 5 / 7 | a C3 marker in the same pots, maximum at pH 5 | DEV; not in any fit |
| de Roos 2005 Table II (`deroos2005_table2_cysamadori_glucose_aa`) | glucose + Cys + Gly or Pro, 90 C: D-fructose-L-cysteine forms readily; Fru-Cys / Fru-Gly 0.9 - 1.1 (15 % water), 4 - 5 (79 % water) | the Cys-hexose Amadori compound exists in situ; amine competition for glucose | DEV (Givaudan); relative RI areas |
| de Roos 2005 Table I (F23#11; `deroos2005_table1_thiol_yields_cys_carbonyl`) | 1 : 1 Cys + sugar, lactate pH 5, 37 % water, boiling 1.5 h: MFT ribose 43.5, xylose 36, glucose < 0.05 mg/mol Cys | at ~100 C the hexose route is shut although its Amadori compound forms (Table II) | DEV |
| Zhang 2026 (`zhang2026_fig3_gsh_glucose_nonvolatiles`, `zhang2026_tables6_s8_volatiles`) | GSH + glucose 60 / 60 mM, 0.2 M phosphate pH 6.5, 130 - 150 C: Cys-, GSH-, Gly-Amadori by LC-MS/MS (figure-only); MFT only at 90 / 120 min (0.33 - 2.32 ng/g, RF 1) | the Cys-hexose Amadori pool in a peptide pot; MFT's late onset | DEV (BTBU); no atlas cell |
| Hofmann 1999 (F06#4, F07#20) | xylose osones peak at 10 min, glucose at 40 min | the hexose supply's time scale | DEV |
| Brands 2001 (F14#0 - #1, F02#13, F05#0) | 3 % caseinate + glucose, 120 C: net glucose loss identical with and without casein; formic x1.4, acetic x2.6 with casein | protein-bound lysine changes products, not sugar loss | DEV |
| Xia 2022 (F02#4 - #5) | Glc + Gly / diglycine / triglycine: TriG-ARP forms 6.9x (100 C), 5.0x (130 C) faster | peptide amine as its own class | DEV (Ho network) |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Brands 2001 (F05#1, F02#12) | fructose + casein: lysine lost without Heyns or furosine | the protein amine's fate is not the free-amine pattern | DEV |
| Cao 2017 (F04#17, F04#18) | the Cys-bearing pool conserved over 5 h; glycine 16 - 24x Cys drains TTCA to ~0.03 M | amine competition in a cysteine pot | DEV (BTBU) |
| Zhang 2025b (F04#19) | lysine 16x Cys moves TTCA and Cys-Amadori | the same with a free basic amine | DEV |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Hofmann 1998 Table 1 (F23#21, F25#10; `hofmann1998_t1_cys_carbohydrates`) | ribose / xylose / fructose / glucose / rhamnose + Cys, pH 5: MFT 19.8 / 14.3 / 2.5 / 1.9 / 0.8 ug; FFT 12.1 / 9.6 / 3.2 / 2.8 / 0.8 | ribose / glucose MFT 10.4, FFT 4.3 (D4) | DEV; ribose / xylose rows B60 fit rows |
| Hofmann 1998 Table 2 (F23#18, F25#8; dossier sec. 2.2) | pH 3 / 5 / 7: glucose MFT 0.3 / 1.9 / 0.4, FFT 0.7 / 2.8 / 0.6; ribose MFT 55.3 / 19.8 / 2.5, FFT 22.9 / 12.1 / 1.2 | the hexose route's pH shape against the pentose one (D2) | DEV; glucose pH 3 / 7 are Maillard-path hold-out bundles |
| Meynier 1995 Table 3 (F25#9, F13#5; CF-001 conditions) | Cys 50 + ribose 35 mM, 0.2 M pyrophosphate held, 140 C, 1 h, pH 4.5 -> 6.5: FFT 8.5 -> 1.4 (x6.1), furfural 3.9 -> tr (> 39x), MFT 15.2 -> tr | E8: FFT falls far less than furfural (D6) | DEV(E8) |
| Hofmann 1998 Table 2 (F25#8) | ribose FFT x19 from pH 3 to 7, highest at pH 3 | E8 | DEV(E8) |
| Cerny 2007 (F25#7, F13#7) | FFT x2.3 down from pH 4 to 6; zero at pH 7 | E8, third lab | v1 DEV; v2 excluded (calibration input) |
| Mottram 2002 (F23#20, F14#8) | buffered vs unbuffered: MFT x2.2 - 2.6, FFT unchanged; ribose-only furfural unchanged | the branches respond differently (shared with M5) | DEV |
| Blank 1997 (F13#11, F14#10) | HDMF from xylose needs glycine or alanine (< 0.01 vs 2.6 ug/mmol) | amine catalysis at 90 C | DEV; calibration input |
| Tang 2024 (F05#16 - #17) | soy vs pea isolate, wet 85 C: lysine -22 vs -8 %; the order flips dry | protein identity as a class constant | DEV (Wageningen) |
| Zhou 2023c, Feng 2022 (F27#24 - #26) | GSH a weaker MFT source than Cys; the ribose-GSH ARP makes MFT itself | peptide vs Cys class | DEV; Feng a calibration input |
| Bolton 1993 (CF-090) | + 6.3 mM xylose leaves MFT unchanged in the glucose matrix | thiamine dominance in hexose foods | DEV |
| Lee 2022 (F06#16), Kocadagli 2016 (F06#15, F07#18) | amino acid raises osones 4 - 10x; NaCl moves the split | amine catalysis | **v1 TEST**: listed, not used |
| Xu 2007 (`heldout_xu2007_tab1.yml`) | pH 5.6 / 7.0 / 8.0 control arms | E8 | keyed held-out (not in a frozen split): item 7 only |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhang 2026 Table 2 (`zhang2026_table2_rate_constants`) | 19 k at 130 / 140 / 150 C; Cys + glucose -> Cys-Amadori Ea ~101, Cys-Amadori -> alpha-dicarbonyls 85.7 kJ/mol (both already B59b's transferred barriers); units not printed, no CI | barriers for part G | DEV (BTBU) |
| Martins 2005, 2005b (F06#12 - #13, F07#15 - #16) | k4 +0.35, k7 +0.60 dec / pH; Ea 97.1 / 107.3 | the hexose split | DEV; calibration input |
| Brands 2002c, Nguyen 2016, Wolf 1977a (F05#3 - #8; CF-072, CF-096, CF-100) | glycation 0.39 dex / pH; Ea 97 - 138 kJ/mol; lysine return disputed | the protein-amine class constants | DEV |
| engine records | B8, B9, B47, B51.1, B20, EX-012, EX-013, red-team probe | sec. 1 table | engine records |

**Left out, with reasons.**
- The thiamine route (F27, E18), except Bolton's bound: a separate MFT source with its own case.
- Acrylamide amine competition (Claeys 2005, Rydberg 2003; F29#3, #21): E9, the acrylamide lane.
- Pyrazine amine identity (Leahy Lys >> Asn, F19#10; Van Lancker dipeptides): M6.
- Strecker aldehydes from glucose pots: M3. Lipid with glucose (Du 2025, Zhao 2019b Table 1): M5.
- Goncuoglu 2026 PPI extrusion (F29#24, v2 TEST E9) and Hamzalioglu 2026 (F05#10 - #11, v1 TEST): protein amine in acrylamide / milk chemistry; held out.

**Counts by kind (sources).** Labelling 6 (one gap); fed intermediates 6; trapping 5 (two held out); measured
intermediates 9; mass balance 3; one factor within a study 13 (three held out); fitted 4 plus the engine records.

## 3. The rivals

| id | part G: glucose + Cys -> MFT | part E: furfural supply vs pH | part A: amines |
|---|---|---|---|
| R0 | **current engine**: ~0 MFT (refused); FFT through `r_glc_fur` | `acid` = [H+] | one amine, folded |
| R1 | **the claim**: Glc + Cys -> Cys-hexose Amadori -> 1-DG + Cys -> 1-DPO + C1 -> the pentose steps (NF and DDP shares as M2 sets them) | amine-catalysed 1,2-enolisation, flat or weakly rising, plus an acid limb | classes with own constants, competition declared |
| R2 | the same supply, entering **below norfuraneol**: 1-DG + amino acid (Strecker / reduction) -> a 1,4-dideoxy C6 (Huber 1990, CF-099) -> a 1,4-dideoxy C5 (DDP) -> MFT | M1's C1: the supply carries de Roos Fig. 4's measured U-shape (flat pH 3 - 5, rising above 5.5) and the split's pH law; no separate E8 term | -- |
| R3 | **C2 + C3 recombination** (B47's opened `k_glc_ha`) | the supply law is right; the steepness sits downstream (furfural loss, the HS- branch) | one constant per alpha-amino group, whatever the carrier |
| R4 | **the entry alone** on the existing amine-free 1-DG supply (the red team's band-now fix) | -- | -- |

## 4. The decisive data

Directions follow from each structure. "[dry]" marks arithmetic on printed values.

| # | measurement (record) | R1 | R2 | R3 | R4 | separates |
|---|---|---|---|---|---|---|
| D1 | **NF / MFT in the glucose vs ribose pot** (Hofmann 1998 Table 5: 7.1 vs 2754) | glucose's 1-DPO feeds NF and MFT like ribose's, so NF / MFT of the DPO-derived share is the same in both pots; [dry] the measured NF ratio glucose / ribose is 2.5e-4 and the MFT ratio 0.096, 390x apart, so R1 needs most ribose NF to come from a pentose-only route (amine-free `r_pent_thermal_dpo`, `r_pent_caramel`) | glucose MFT bypasses NF: low NF / MFT expected | no NF from fragments: NF / MFT ~ background | as R1 | R1 from R2, scale-free; nearest the intermediate. HexoseC5Entry's own docstring names this falsifier |
| D2 | **glucose MFT pH shape** (Hofmann Table 2: 0.3 / 1.9 / 0.4 at pH 3 / 5 / 7; ribose 55.3 / 19.8 / 2.5) | low at pH 3 (2,3-enolisation minor in acid, P3), falling at pH 7 as beta-cleavage takes 1-DG (Davidek 2006a) | as R1 below 5; above 5 depends on the Strecker step's pH law | rising to pH 7, like fed C2 + C3 (0.01 / 0.23 / 0.27 mol %, F24#8) | the amine-free 1-DG supply's pH law, which no measurement on disk fixes | R3 from R1 / R2 by sign; the pentose shape's opposite sign at pH 3 vs 5 tests that the hexose route is not the pentose route |
| D3 | **mixed isotopomers** (Du 2023, Zhao 2019b) | intact C5 | intact C5 | ~50 % mixed in a 1:1 pool | intact C5 | **rejects R3 by a measured share** (data, not a fit); does not separate R1, R2, R4 |
| D4 | **ribose / glucose MFT 10.4 and FFT 4.3** (Table 1) and their temperature trend: > 870 at ~100 C (de Roos Table I), 10.4 at 145 C | follows the entry's barrier; [dry] if the two labs' pots were comparable, the hexose route's apparent barrier exceeds the pentose route's by at least ~128 kJ/mol (ln(870 / 10.4) over 100 vs 145 C) (they are not: 37 % water lactate vs 0.5 M phosphate) | as R1 | -- | the probe gave 45 - 78 at 145 C | weakly: endpoint levels that E6, E11 and the H2S budget also move. The ratio cancels scale; it is not nearest the mechanism |
| D5 | **the Cys-hexose Amadori pool** (de Roos Table II; Zhang 2026 Fig. 3) with no hexose MFT at ~100 C (de Roos Table I) | the pool exists and its decay to 1-DG is the slow step at 100 C | as R1 | no Amadori step needed | no Cys-hexose Amadori compound: the measured pool is unexplained | R1 / R2 from R4; nearest the supply |
| D6 | **FFT vs furfural vs pH** (Meynier Table 3: FFT x6.1, furfural > 39x, MFT / FFT 1.79 -> < 0.07; Hofmann F25#8 x19 over 3 -> 7) | the supply falls less than a decade per unit; FFT's fall needs H2S or furfural loss at high pH | (E-R2) the supply is flat pH 3 - 5 and rises above 5.5, so the whole FFT fall is downstream or in the split | (E-R3) furfural standing level falls with its own loss; FFT per furfural rises (HS-) | -- | E-R0 from the rest (its supply alone gives ~x400 over 4.5 -> 6.5 in the probe). E-R1, E-R2, E-R3 need the split (M1b's D-rows) and furfural's loss; **not separable on FFT alone** |
| D7 | **amine class at fixed sugar** (Xia 2022 TriG-ARP x5 - 7; de Roos Table II Fru-Cys / Fru-Gly 1 -> 4 - 5 with dilution; Brands 2001 casein leaves glucose loss unchanged) | each class its own constant and dilution law | -- | (A-R3) equal per alpha-amine: Xia's x5 - 7 at equal molarity and de Roos's dilution-dependent ratio contradict it | -- | A-R1 from A-R3, within-study. One lab each |

**The plan's decisive data, checked.**
- "Hofmann ribose / glucose ratios (MFT 10.4, FFT 4.3)": D4 here. The ribose arm is a B60 fit row; the ratio is an
  endpoint that E6, E11 and the H2S budget move, so under sec. 7 item 4 it is secondary to D1 and D2, which come from
  the same pots and cancel scale too.
- "The FFT pH slope": D6. It cannot separate the amine-catalysed term (E-R1) from M1's C1 (E-R2) without the split,
  and both laws act on the same steps M1 changes.

**What the data cannot decide.**
- R1 against R2 rests on D1, one laboratory (DFA), one temperature, with NF and furfural by an unstated method.
- Which glucose carbon is lost: no positional label on disk.
- No fed hexose intermediate to MFT exists (no fed 1-DG + H2S, no fed Cys-hexose Amadori compound).
- Part A: one laboratory per class contrast; no Maillard flavour pot varies protein-bound vs free amine at fixed
  sugar, pH and temperature.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1, D2, D4 Hofmann 1998 Tables 1, 2, 5, 9 | triplicates, "differ by <= 10 %" (<= 0.04 dex per value) for the SIDA thiols; NF and furfural quantitation method not stated, no noise printed | ledger `hofmann1998_t1_cys_carbohydrates` note; `hofmann1998_t5_fa_nf_pairs` flags |
| D2 fed C2 + C3 ladder | as Table 1 (same lab, SIDA) | hofmann1998_reconciliation dossier |
| D3 Du 2023 Table 2 | isotopomer percentages as printed, no SD; two blind readers identical on all 59 rows; [M-H]+ correction not stated (MFT prints only M+0 and M+5, so a leak would have shown as M+4) | `heldout_du2023_tab2.yml` header; du2023 dossier |
| D3 Zhao 2019b Table 2 | triplicate reactions, means; no SD on the isotopomers | zhao2019b dossier |
| D4 de Roos Table I | quantitation not stated; "< 0.05" is a reporting floor; boiling vs SDE time ambiguous | deroos2005 dossier sec. 3.3 |
| D5 de Roos Table II | relative RI peak areas, response factors not given (mol ratios within ~x1.5 by the dossier's estimate) | deroos2005 dossier sec. 3.5 |
| D5 Zhang 2026 | n = 3, mean +- SD; volatiles RF 1 (ratios only); Fig. 3 figure-only, not digitised | ledger `zhang2026_*` |
| D6 Meynier 1995 Table 3 | triplicate; post-heating pH drift <= 0.2 | meynier1995 dossier; F25#9 |
| D6 Hofmann F25#8 | as D1 | -- |
| D7 Xia 2022, de Roos, Brands 2001 | Xia: not checked here; de Roos as D5; Brands digitised +-3 % of axis span | dossiers |
| across laboratories (D4's 100 vs 145 C trend; Bolton's bound) | between-lab spread from H0 (EX-059, pending) | -- |

## 6. What M4 needs from the one network (W3)

- **Species:** Glc, Fru, Cys, the Cys-hexose Amadori compound (new), 1-DG shared by the hexose and sulfur modules,
  1-DPO, DDP, NF, MFT, furfural, FFT; for R2 a 1,4-dideoxy C6. Today 1-DG lives in the trunk and 1-DPO in the sulfur
  lane, which do not compose.
- **Amine classes as species, not a convention:** free alpha-amino acid (with identity factors), cysteine,
  peptide (GSH, diglycine), protein-bound lysine (one key, A2b), each declaring which sugars it competes for, so
  `compose` refuses undeclared competition.
- **One pot with all of them:** Glc + Cys + Gly + protein in one integration (today the lanes refuse to combine).
- **In-situ pH on the hexose Amadori and enolisation steps,** with the Cys-hexose Amadori compound's pKa in the one
  titratable table; and a per-step hook for M5's buffer factor (every Hofmann pot is 0.5 M phosphate).
- **Isotopologue pools:** labelled and unlabelled copies of the sugar carbon (or atom tracking) so D3 is an observable,
  not a hand calculation.
- **The forced-input driver** (A2 contract sec. 2): 1-DG / 3-DG courses (Deng 2021) clamped while the entry is tested.

## Open questions for the build

1. **Order with M1 and M2.** Part G feeds the pentose steps whose NF / DDP shares M2 refits and whose supply M1's C1
   replaces. D1 is read on the network after M2's verdict, or with M2's rivals carried as a factor.
2. **The Hofmann glucose pH 3 / 7 arms.** DEV under both frozen splits, but declared Maillard-path hold-out bundles.
   This case treats them as decisive (plan sec. 7: held-out status is the split files'), scored once in item 6, never
   a fit row. If the older declaration binds, D2 moves to item 7 and R1 vs R2 rests on D1 alone.
3. **The Cys-hexose Amadori step's constants.** Zhang 2026's Ea ~101 is from a GSH pot with units not printed; de Roos
   Table II gives only relative areas. Declare the formation constant's class band before any fit.
4. **E8 after M1.** If C1 ships, does part E survive for `r_glc_fur` and `r_arp_tdp` alone? Decide before building
   the amine-catalysed term.
5. **Part A's constants** come from one laboratory per class (Xia, de Roos, Brands, Tang). Declare them as bands with
   a per-class basis; no class constant is fitted on flavour rows.
6. **W6 request:** the BTBU / CAAS MFT yields per glucose consumed (ENGINE_VS_LITERATURE sec. 3.1, E1) would give D4
   a second laboratory.

Seen by the author: as stated at the top. No engine output on any row named in sec. 4 was computed for this case;
the Appendix A and red-team probe outputs on the Hofmann T1 design were read.
