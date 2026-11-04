# Case M5: buffer, fat and water activity

*Plan step M5 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E5 (F01-b, F03-b, F13-b, F14-a,
F19-c, F23-c, F29-b) and E7 (F30-a, F28-d, F26-c) of
[ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in sec. 5 (E5 "SUPPORTED,
step-specific"; E7 SUPPORTED) and sec. 7 (E5's form; E7 / CF-031); the a_w mismatches F05-a, F16-d, F19-d, F25-d,
which carry no E number). Items 1 - 5 of sec. 7, written 2026-11-04 before any result of the change exists. No engine
was run for this case. Items 6 - 9 are appended later and never edited into this part. Background physics:
[PH_BUFFERS_AND_WATER](../../guides/PH_BUFFERS_AND_WATER.md) secs. 6 - 7.*

**Seen by the author.** No engine output on any held-out row named here. Measured values of held-out sources read
while writing, through the conflict register and ENGINE_VS_LITERATURE sec. 7: the v2 E5 TEST contrasts F02#1
(Potman 1989), F04#2 (Kallen 1971b), F05#23 / F06#11 / F08#2 (Thornalley 1999) as quoted in CF-093 / EX-058; the v2
E7 TEST contrast F25#17 (US 8,920,862, MFT x10.1) as quoted in CF-031; Larralde 1995's ratio (keyed held-out); the
atlas text of F30#11 (Kim 1996b), F30#15 (Zamora 2015) and F30#16 (Hidalgo 1993), measured values only. Yiltirak 2026
engine numbers quoted in CF-068 and SULFUR_FLUX_MAP were read for case M1, not for this one.

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml).
- v2 E5: DEV labs Auburn (Bell), Baynes, Ghent, Leuven, Nestle / Cerny, Reading, Rizzi, SCUT (Guan), Wageningen.
  TEST: F02#1@potman1989, F04#2@kallen1971b, F05#23 / F06#11 / F08#2@thornalley1999. Excluded as calibration input:
  F19#8 (Leahy), F27#9 (Arnold), F28#4 (Zhou 2023).
- v2 E7: DEV labs CAAS, Huazhong (Pan), Reading. TEST: F25#17@us8920862 (viewed: a carried failure, CF-031),
  F30#11@kim1996b, F30#15@zamora2015, F30#16@hidalgo1993. Xu 2011's lard contrast is also **F28#20, a v2 E12 TEST
  contrast**; by CF-031's red-team ruling it must not calibrate an E7 term, so this case reads Xu 2011 in item 7 only.
- Keyed held-out tables outside the frozen splits (README secs. 10 - 11), item 7 only: `heldout_larralde1995.yml`
  (E5), `heldout_zhang2021_table3.yml`, `heldout_zhao2017b_tab1.yml`, `heldout_wang2022_tab2.yml`,
  `heldout_zhao2019b_tab1_2.yml` Table 1, `heldout_zamora2015_tab1.yml` (E7).
- Water activity has no E number and so no v2 entry; v1 governs. Warmbier 1976 (Minnesota / Labuza) and Pereyra
  Gonzales 2010 (Buenos Aires) are v1 TEST labs; Pereyra Gonzales is also a calibration input (the trunk's a_w table).
  Eichner 1972, Bell 1995, Lievonen 2002, Miao 2004, Buera 1993, Modi 2008, Cui 2026, Bassama 2011 and Zhang 2022
  (Hefei) have no atlas cell and sit in no frozen split.

## 1. The claim

Three conditions that every lab constant carries silently, each turned into a term on named steps.

| part | step(s) | today (B60, read from the code) | M5 |
|---|---|---|---|
| B. buffer catalysis | Amadori decay (`k_ama_*`, `r_arp_*`, the Cys-Amadori steps if M1's C1 ships), the 2,3-enolisation branch, `r_schiff` (lumped), TTCA formation (CF-011) | no buffer term anywhere: the buffer is a titrant in the charge balance only ([ph_state.py](../../../src/kinetic_core/ph_state.py)); every constant carries its calibration pot's phosphate (Hofmann 0.5 M, Martins 0.1 M, Cao 0.2 M) | k_s = k_s,u (1 + sum_b kappa_{s,b} c_b) on steps the network resolves; a hyperbolic factor on a lumped step whose uncatalysed part can limit; carboxylates and free-sugar steps on their own bands (ENGINE_VS_LITERATURE sec. 7) |
| F. fat | MFT, FFT (and their disulfides), Cys and H2S against the lipid lane's carbonyls | lanes integrate as a direct sum (`lipid.lane_coupling_verdict`, [lipid.py](../../../src/kinetic_core/lipid.py) ~784 - 830); no oil phase, no thiol partition anywhere in `src/`; the lipid lane's aldehyde-protein `COVALENT_SINK` is disabled by its measured Ea ([parameters_lipid.py](../../../src/kinetic_core/parameters_lipid.py) ~551 - 589) | the plan's oil-water partition of the thiols; this case adds the rival the register names: a lipid-electrophile capture of the MFT pool (CF-031) |
| W. water activity | `k_schiff` (trunk), the sulfur lane's steps | trunk: Pereyra Gonzales's multiplier table on `k_schiff` only (`AW_STEPS`, [trunk_conditions.py](../../../src/kinetic_core/trunk_conditions.py) ~60 - 85); acrylamide: a flat band (CF-104); sulfur and lipid lanes refuse a_w | B54's design: concentration per whole phase, the trunk table re-derived as a residual, a Gordon-Taylor Tg gate; a step-specific term in the sulfur lane only where a sign contrast demands it |

**Earlier waves on these structures** (cited, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| per-step phosphate catalysis | B53 (backlog; PH_BUFFERS_AND_WATER sec. 6) | planned, deferred 2026-10-27 (no blind test), unblocked 10-28 (D5), **never run**. Fit source Davidek 2002's phosphate series; hold-out Mottram & Nobrega 2002 | -- |
| the factor's form | EX-058 (resolve2_pyrbuf, red-teamed) | non-saturating to 0.5 M, order 1.0 - 1.35 (P 0.8): **decided on Potman (E5 TEST)**, with Kallen and Thornalley (E5 TEST) and Larralde (keyed held-out) | arithmetic on data |
| buffer in the Hofmann bundles | B51 data facts (Amendment 44) | 0.5 M phosphate stated on every arm (a titrant fact) | data correction |
| water activity basis and Tg gate | B54 (backlog) | planned, **never run**; needs a pea / soy / gluten isotherm and Tg curve | -- |
| acrylamide a_w | CF-104 red team (2026-11-03) | flat multiplier upheld (P 0.75); band-drop fix applied | review |
| lipid x sulfur | CF-031 cluster 4 (2026-11-01) | a partition operator explains at most 0.1 - 0.4 dex of the 0.5 - 2.2 dex falls (P 0.55) | arithmetic |
| lipid x sulfur | resolve2_sulfur (CF-031) | the pre-registered speciation rule FAILED; partition not excluded | arithmetic |
| aldehyde-protein covalent sink | B8 / K6b | disabled by a measured Ea (15 - 23 kJ/mol) | measured constant |

**Errors stack (CF-107).** No whole-model fit has ever tested a buffer, fat or a_w term, so no refusal needs
discounting. The converse risk is real: every sulfur constant was fitted in 0.5 M phosphate, so adding part B
re-references them, and the joint refit (JR) must carry it. Candidate errors that a buffer term could absorb: E2 (the
split frozen in pH), CF-084 (the pH balance), the Hofmann pH 3 MFT miss (SULFUR_FLUX_MAP sec. 5 names phosphate on the
enolisation as one candidate). For fat: `k_mft_decay` (frozen, pot-independent) can mimic any MFT-selective loss.

## 2. All the evidence, by kind

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Du 2025 (F30#8) | CAMOLA: 2-methylthiophene 100 % glucose carbon without fatty acid, 38 - 48 % with it | lipid carbon enters sulfur heterocycles (F) | DEV(E7) |
| Du 2023 (F30#9) | 2-pentylthiophene only from the 2,4-dienal, ~84 % without glucose carbon | H2S + dienal capture (F) | **v2 TEST(E20)**: listed, not used |
| Zhao 2019b Table 1 / text (F23#7) | label patterns unchanged with chicken fat | fat does not reroute MFT carbon | keyed held-out (E7): item 7 |
| Whitfield 1988, Elmore 2002 (F26#21, #24, #25) | alkylthiophenes / thiapyrans from PUFA + H2S | the sulfur sink into lipid products | DEV(E7) |
| none | no label on buffer or a_w routes | -- | gap |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Davidek 2002 (F03#0 - #2) | fed DFG 0.1 M, 90 C, **pH-stat at 90 C**, water vs 0.1 M phosphate: decay x2.8 - 3.6 at pH 5 - 7, small at 8; reversion up most at pH 5 - 6 | **buffer catalysis at held pH** (D1); B53's fit source | DEV(E5) |
| Zyzak 1995 (F03#7 - #9, F05#22, F06#10) | fed fructoselysine, 37 C, N2 + chelator: t1/2 57.1 / 40.9 / 24.4 / 16.8 d at 10 / 50 / 100 / 200 mM phosphate; HEPES ~12x slower; products unchanged | dose law (D2); step-specific | DEV(E5) |
| Guan 2012 (F01#9, F03#10) | fed proline ARP, 110 - 130 C, pH 5.8: phosphate x3.8 (110 C) / x1.5 (130 C); acetate, citrate +16 / +34 % | carboxylate kappa; kappa falls with T (D4) | DEV(E5) |
| Hofmann 1998 Table 2 dry arm (F24#9) | fed C2 + C3: dry 180 C MFT 1554 vs aqueous 268 ug (T, t, a_w, buffer confounded) | a_w on one fed step | DEV |
| Elmore 2002 (F30#0 - #1, F25#16, F28#19) | Cys 25 + ribose 25 mM, pyrophosphate pH 5.5, 140 C, +- 25 mM methyl linoleate / linolenate: MFT 807 (+-508) -> 58 (+-71) -> tr; FFT 92 -> 87 -> 17; disulfide > 20x down; enals consumed | MFT-selective loss with PUFA (D5, D6) | DEV(E7) |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Esterbauer 1975, Doorn 2002 (F28#17, F30#18 - #19) | thiol + enal Michael adducts, reversible in hours (2-alkenals); thiol >> His > Lys | the capture chemistry exists and reverses | **v1 TEST**: listed, not used |
| Charles-Bernard 2005 (F28#0, #16) | 25 C: carbonyl-competing nucleophile slows thiol loss 4 - 13x | carbonyl capture of thiols (low T) | DEV |
| none at process T | no trapped lipid-thiol adduct at >= 100 C | F | gap |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Cao 2017 (F04#9, F04#17) | 0.2 M phosphate pH 5.5, 90 C: TTCA K flat while pH drifts | the phosphate-catalysed exchange (CF-011) | DEV |
| Zhang 2022 (Hefei; no atlas cell; CF-031) | Xyl + Cys, pH 7.2 unbuffered, 120 C, 90 min, 4 : 1 oil : water, phases separated: water-phase free MFT x0.36 - 4.8; the control's 2.81 ng/mL disulfide vanishes from raw-oil water phases; MFT and disulfide found in the oil | **both species partition into oil** (D7); n, SD not printed | no frozen split |
| Tang 2019 (F02#9 - #10) | GSH + xylose, 80 C, vacuum dehydration: ARP x8 in a tenth of the time | water removal speeds Amadori formation | DEV |
| Kang 2026 (F04#25) | TTCA powder, a_w 0.11 / 0.43 / 0.84, 60 d: retention 86 / 72 / 64 % | a_w on the cysteine reservoir (storage) | DEV; calibration input |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Farmer 1992, Xu 2011 (F30#2, F30#20; CF-088) | lipid aldehydes fall with Cys + pentose while saturated ketones do not: capture, not antioxidant | the reverse flow (lipid lane loses carbonyls) | Farmer DEV(E7); Xu item 7 |
| Potman 1989 (F02#0 - #1) | ARP carries all glucose lost at <= 1 h; phosphate order 1.34, ~x15 at 1.2 M | the dose law | **v1 TEST, v2 TEST(E5)**: listed, not used |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Bell 1997 (F01#7 - #8) | Glc + Gly 0.1 M, pH held 6.86 - 7.05, 25 C: glycine loss 0.00055 / 0.00082 / 0.0019 / 0.0033 M/month at 0.02 / 0.05 / 0.2 / 0.5 M phosphate; citrate not different from zero | **dose law at held pH** (D2); carboxylate (D4) | DEV(E5) |
| Mottram & Nobrega 2002 (F23#20, F13#10, F14#8) | Cys + ribose 50 / 50 mM, 140 C, 30 min: pyrophosphate or phthalate pH 5.6 vs unbuffered (5.6 -> 3.9): MFT x2.2 - 2.6, FFT unchanged; ribose-only furfural 107 vs 103 | **step specificity** (D3); pH-confounded (unbuffered drifts) | DEV(E5) |
| Blank 1997 (F13#11) | HDMF from xylose + Gly, 90 C: water / 0.2 M phosphate / 0.2 M malonate 0.06 / 2.6 / 0.6 ug/mmol | phosphate vs carboxylate on a 2,3-enolisation product | DEV(E5); calibration input |
| Rizzi 2004 (F14#4) | ribose alone, 0 / 0.05 / 0.10 M phosphate vs bis-tris, pH 7.3: A420 0.003 -> 0.045 -> 0.262 | phosphate drives an amine-free sugar step | DEV(E5) |
| Adams 2008 (F19#6) | fed MGO + Ala, pH 7, 130 C: phosphate vs bicarbonate DMP x1.8, pyrrole x0.36 | buffer on pyrazine steps (M6) | DEV(E5) |
| De Vleeschouwer 2006 (F29#6) | acrylamide k_F citrate vs phosphate x0.14 (pH 4), x2.4 (pH 6) | carboxylate vs phosphate by pH | DEV(E5); calibration input |
| Yiltirak 2026 (F30#6; dossier Table S3) | ribose + Cys 25 mM, 0.5 M phosphate pH 5.5, four equivalent cooks: buffer + 10 % canola oil MFT 0.72 vs 6.88 ug/L (x1/9.6, 100 C) while FFT 1.19 vs 1.28; emulsion MFT 3.73 | **fresh-oil MFT-selective fall** (D5) | DEV(E7); the buffer arm is a Maillard-path hold-out bundle, the oil arms are not |
| Whitfield 1988 (F30#3, F25#15) | lecithin: MFT x1/3.0, FFT x1/2.0, methylpyrazine x1/2.1 | moderately selective | DEV(E7) |
| Farmer 1990 (F25#14) | four lipids: FFT x0.62 - 0.72, furfural x0.49 - 0.80 | FFT side | DEV |
| Du 2025 (F30#7) | GSH + glucose 200 / 200 mM, 145 C, +- 10 mM free C18 acids: MFT +19 - 30 % | the small opposite sign (CF-031) | DEV(E7) |
| Xu 2011 (F30#5, F28#20) | fresh / oxidised lard: MFT x1/21 / x1/143, FFT x1/3 / x1/6.4 | oxidation dose | item 7 only (E12 TEST) |
| Hofmann 1998 Table 2 dry arm (F25#11 - #12, F23#18) | dry 180 C / 6 min vs aqueous 145 C / 20 min: ribose FFT x8.0, glucose and rhamnose FFT x0.50; MFT x1.3 / 2.2 / 3.9 | **a_w sign contrast across sugars** (D8); T, t, buffer confounded but shared by the sugars | DEV |
| US 10,076,126 (F25#13) | glycerol vs buffer, 110 C: FFT up to x56 | a low-water solvent on FFT | DEV |
| Granvogl 2006 (F16#8) | silica, 10 % water: Leu / Ile aldehyde x2 - 3; Strecker amine 19 - 27 % | a low-water exit (M3) | DEV |
| Leahy 1989a (F19#11) | a_w 0.32 / 0.58 / 0.75 / 0.84: pyrazine rate 1 : 2.10 : 4.23 : 4.07 | the a_w maximum (D9) | DEV; calibration input (not B18 fit rows) |
| Scalone 2015 (F19#15) | a_w 0.33 -> 0.75: 2,5(6)-DMP -0.3 dex | a_w on pyrazines | DEV |
| Modi 2008 (no atlas cell) | Cys + thiamine + ribose in beef fibre, 130 C: a_w 0.60 / 0.80 / 0.90; peak areas, n = 4 | **the only a_w series for MFT** | no frozen split |
| Bell 1995 (no cell) | at fixed molality, glassy and rubbery polymers give the same rates; reactant concentration dominates | concentration basis (D9) | Auburn, a DEV(E5) lab |
| Eichner 1972 (no cell) | 37 C browning maximum at a_w 0.41 - 0.52 (x38 from 0.84 to 0.41); end points only | a_w maximum (low T) | no frozen split |
| Lievonen 2002, Miao 2004, Buera 1993 (no cells) | (T - Tg) not sufficient; Ea depends on water content; WLF fits | the Tg gate's limits | no frozen split |
| Tang 2024 (F05#17) | dry a_w 0.6 vs wet: dry glycates far more; pea / soy order flips | protein glycation vs regime | DEV |
| De Vleeschouwer 2007 / 2008 (F29#7 - #8) | acrylamide formation flat over a_w 0.34 - 0.99 | concentration basis on a fitted lane | DEV; calibration input |
| Warmbier 1976, Pereyra Gonzales 2010 (F05#18 - #21) | lysine loss maximum 0.5 - 0.6; flat 0.31 - 0.71, x4 slower at 0.98 | the shipped table | **v1 TEST**: listed, not used (Pereyra is already the shipped table) |
| Hartman 1984b (F27#20) | sulfur volatiles peak near a_w 0.72 | a_w on thiamine sulfur | **v2 TEST(E18)**: listed, not used |
| Kallen 1971b (F04#2), Thornalley 1999 | saturating thiazolidine catalysis; convex 10 - 100 mM | the form | **v2 TEST(E5)**: listed, not used |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Bassama 2011 (CF-104) | plantain, a_w 0.97 / 0.90 / 0.43: k_F x8.1, k_E x14; wet-basis mass action fits 0.43 / 0.90, misses 0.90 / 0.97 by x2.6 | concentration basis in a food | no frozen split |
| Martins 2005, Knol 2005 | every constant measured in 0.1 M phosphate | what part B re-references | DEV; calibration inputs |
| engine records | B51 (Amendment 44), EX-058, CF-031 clusters, CF-104 | sec. 1 table | -- |

**Left out, with reasons.** Lipid Strecker and dienal + NH3 pyridines (Zamora, Hidalgo, Kim 1996b; F18, F21): E20.
Lipid oxidation kinetics themselves (Frankel, Schroen): the lipid lane, not its coupling. Matrix binding of
aldehydes to proteins (Snel, Bi): the matrix layer. Thiamine a_w (Hartman): E18.

**Counts by kind (sources).** Labelling 4 (one held out); fed intermediates 5; trapping 3 (two held out);
measured intermediates 4; mass balance 2 (one held out); one factor within a study ~30 (six held out); fitted 2
plus records.

## 3. The rivals

| id | B. buffer | F. fat | W. water activity |
|---|---|---|---|
| R0 | **current engine**: none | direct sum, lipid inert | trunk table on `k_schiff`; refusal elsewhere |
| R1 | **the claim**: linear-with-intercept factor per resolved step, kappa by step and buffer species | **the plan's claim**: oil-water partition of MFT, FFT, their disulfides (fast equilibrium, reactions in water) | B54: whole-phase concentration, residual table, Tg gate |
| R2 | saturating factor (Michaelis, K ~ 17 - 29 mM as Kallen's) on every step | lipid-electrophile capture of the MFT pool, first order in MFT and lipid-carbonyl load, MFT-selective (CF-031's candidate) | concentration only, no residual table |
| R3 | buffer acts only through the in-situ pH it holds (the engine's charge balance already has this) | partition + capture | step-specific terms (the dehydration-limited condensation; a low-water exit), no shared factor |
| R4 | step-specific: Amadori decay and 2,3-enolisation only, not 1,2-enolisation or furfural | supply-side: lipid carbonyls take Cys and H2S (Zhao 2017b; thiophenes), thiols fall through their precursors | -- |

## 4. The decisive data

| # | measurement (record) | where the rivals differ | separates |
|---|---|---|---|
| D1 | **buffer at held pH**: Davidek 2002 pH-stat, water vs 0.1 M phosphate (x2.8 - 3.6 at pH 5 - 7); Bell 1997's pH-held ladder | B-R3 predicts ratio 1 at equal held pH; B-R1, R2, R4 predict > 1 | B-R3 from the rest. Davidek is also B53's fit source: the existence test uses the ratio, not the kappa a fit would take from it |
| D2 | **dose law at held pH**: Bell 1997 0.0019 vs 0.0033 M/month at 0.2 vs 0.5 M (ratio 1.74; 95 % CL in Table 1); Zyzak 1995 k 0.0121 / 0.0169 / 0.0284 / 0.0413 /d at 10 - 200 mM; Guan 2012 0.1 - 0.3 M | B-R2 with K <= 0.1 M saturates: k(0.5) / k(0.2) <= ~1.1 and Zyzak flat above 50 mM; B-R1 rises across the range; [dry] Bell's four rates sit on k0 + kappa c with kappa ~10 /M (least-squares slope 0.0056 M/month per M over a 0.00056 intercept; the backlog's B53 entry reads 10 - 11) | B-R1 from B-R2 on DEV data alone. The red team's form (EX-058) rests on TEST data and is not used to choose |
| D3 | **branch specificity**: Mottram 2002 MFT x2.2 - 2.6 with FFT and ribose-only furfural unchanged; phthalate ~ pyrophosphate at pH 5.6 | B-R4 (2,3-branch and Amadori only) predicts MFT up, FFT flat; a general factor on every sugar step predicts both up; B-R3 is confounded here (the unbuffered pot drifts to 3.9) | B-R4 from a general factor. Scale-free (a ratio of ratios), one lab |
| D4 | **carboxylate vs phosphate**: Guan 2012 citrate / acetate kappa 3.4 / 1.6 vs 28 - 38 /M (110 C); Bell citrate 0 at 25 C | the factor per buffer species | whether foods' citrate / malate act like phosphate; one lab per temperature |
| D5 | **MFT / FFT selectivity of the lipid fall** (Elmore 2002 18:2: MFT x1/14, FFT x0.95; Yiltirak fresh oil: MFT x1/9.6, FFT x0.93 at 100 C, x0.57 at 120 - 130 C; Whitfield x1/3.0 vs x1/2.0) | F-R1 lowers MFT and FFT by 1 / (1 + phi K_i), and K_MFT ~ K_FFT (log K_ow 1.5 - 2.0, an estimate): the MFT / FFT ratio barely moves; F-R2 is MFT-selective; F-R4 lowers both (shared H2S, Cys) | F-R2 / F-R3 from F-R1 and F-R4, scale-free. [dry] Measured MFT / FFT falls x13 in Elmore (8.8 -> 0.67) and x9 (100 C) to x2.2 (130 C) in Yiltirak. Yiltirak's 30 uL aliquot of a non-emulsified oil arm: whether it was water or whole sample decides what partition would show (CF-031 author question) |
| D6 | **lipid unsaturation at equal amount** (Elmore 18:2 vs 18:3; Farmer's four lipids) | F-R1 depends on oil volume, not on unsaturation; F-R2 scales with electrophile load | F-R2 from F-R1. The oxidation-dose series (Xu 2011, Wang 2022, Zhang 2021) are held out |
| D7 | **phase-resolved thiols** (Zhang 2022: MFT and its disulfide in the oil phase; the water-phase pool shifts) | F-R1 requires the oil to hold the missing mass; F-R2 requires the pool to vanish | F-R1 vs F-R2 directly, but n / SD not printed and the free + disulfide rule failed (CF-031); a measured K_oil/water for MFT and FFT is not on disk |
| D8 | **a_w sign across sugars**: Hofmann dry vs aqueous ribose FFT x8.0, glucose FFT x0.50 under one regime change | W-R2 (concentration only) moves both sugars the same way; the opposite sign needs a step-specific term (W-R3) or a sugar-specific route | W-R2 rejected as the sole sulfur-lane law (a sign); T and time are confounded but shared |
| D9 | **the a_w maximum**: Leahy 1989a pyrazine 1 : 2.1 : 4.2 : 4.1 at a_w 0.32 - 0.84; Bell 1995 equal rates at fixed molality; Modi 2008 MFT at 0.60 / 0.80 / 0.90 | W-R2 predicts a monotone rise as water falls until the glass; W-R1's Tg gate and residual table carry the fall below the maximum | W-R1 / R3 from W-R2 below ~0.75. Modi is peak areas without calibration (ratios only) |

**The plan's decisive data, checked.**
- "Potman and Kallen (E5 TEST)": v2 TEST(E5) in `split_v2_per_change.yml`, Potman also v1 TEST. Under sec. 7 they are
  read in item 7 only. **The form the plan row calls "settled" (linear to 0.5 M) was chosen on Potman**; on DEV data
  D2 can choose it (Bell, Zyzak, Guan), which the build must do before item 7.
- "The unseen fat contrasts (Kim 1996, Zamora 2015, Hidalgo 1993)": v2 TEST(E7), item 7 only. They are dienal + amine /
  lipid Strecker / epoxyalkenal browning chemistry (E20, F18), not thiol loss, so they could not decide F-R1 vs F-R2
  even if they were development data.

**What the data cannot decide.** Partition vs capture rests on scale-free selectivity (D5) and one phase-resolved pot
without SDs (D7). No DEV dose series in oil fraction exists. No a_w series for the thiols has calibrated amounts.
Buffer kappa at more than one temperature exists only in Guan (and it falls 0.4 - 0.75 dex from 110 to 130 C).

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1 Davidek 2002 | "All samples analyzed in duplicate", no error bars; first-order k band across three fitting methods (0.5 - 0.8 h induction at pH 5 - 6) | davidek2002 dossier |
| D2 Bell 1997 | Table 1 95 % CL per rate (+-0.00050 at 0.02 M to +-0.0009 at 0.5 M); one sample per time point | bell1997 dossier |
| D2 Zyzak 1995 | t1/2 +- SD, triplicates | zyzak1995 dossier |
| D2, D4 Guan 2012 | figure read errors (+-5 g/L glucose; DFP +-1e-4 mmol/L) | guan2012 dossier |
| D3 Mottram 2002 | triplicate; mean CV ~22 %, no compound > 40 %; RF 1 (ratios only) | mottram2002 dossier |
| D5 Elmore 2002 | four replicates; MFT 807 +- 508 and 58 +- 71 (SD 63 % and 122 % of the mean) | elmore2002 dossier; CF-031 |
| D5 Yiltirak 2026 | triplicate cooks, mean +- SD (Table S3), SIDA with deuterated standards | yiltirak2026 dossier |
| D5 Whitfield 1988 | n = 3, external standards | whitfield1988 dossier |
| D6 Farmer 1992 | n = 4, mean + SD | farmer1992 dossier |
| D7 Zhang 2022 | n and SD not printed; HS-SPME RF 1 | zhang2022 dossier |
| D8 Hofmann 1998 | triplicates within 10 % | hofmann1998 ledger |
| D9 Leahy 1989a | duplicates early, triplicates later | leahy1989a dossier |
| D9 Modi 2008 | peak areas, mean +- SD, n = 4; no calibration | modi2008 dossier |
| across laboratories (D2 across Bell / Zyzak / Guan; D5 across Reading / CAAS) | between-lab spread from H0 (EX-059, pending) | -- |

## 6. What M5 needs from the one network (W3)

- **Buffer as a composition, not a label:** each buffer species with its concentration and pK (already in
  `ph_state` for the charge balance), and a per-step catalytic coefficient table keyed by (step, buffer species) as
  `ConfigKey`s with bands; a hyperbolic form available for lumped steps.
- **Constants re-referenced to buffer-free values** at compose time from each calibration pot's stated buffer, so the
  joint refit sees one convention.
- **Two phases:** an oil phase with a volume fraction, fast partition equilibria per volatile (K_oil/water as declared
  constants), reactions confined to their phase, and observables that say which phase or the whole sample they read
  (the CF-031 tags `mft_species` and `phase`). The headspace operator must use the same partition once.
- **Lane coupling:** the lipid lane's carbonyls (alkenals, dienals) as shared species that react with MFT, Cys and H2S
  (retiring the direct sum), with the lane-coupling guard replaced by declared competition.
- **Water:** moisture and a_w as conditions, a whole-phase concentration basis, a Tg gate (warn within 20 K, refuse
  below), and an a_w hook per step for the step-specific rival.

## Open questions for the build

1. **Choose the buffer form on DEV data first** (D2), then read Potman, Kallen and Thornalley in item 7 with the form
   fixed. The plan row's "settled" form should be re-labelled as provisional.
2. **Which steps get kappa.** The 2,3-branch / Amadori steps exist in the sulfur lane only as `r_pent_dpo`,
   `r_arp_dpo` or, if M1's C1 ships, `r_cama_dpo`. Decide after M1's verdict.
3. **Kappa's temperature dependence.** Guan's fall with T is the only measurement; PH_BUFFERS_AND_WATER sec. 10
   carries kappa T-flat with a x/÷1.5 band. Declare which.
4. **K_oil/water for MFT and FFT** is an estimate (log K_ow 1.5 - 2.0). F-R1 needs a measured value or a declared band.
5. **Yiltirak's aliquot.** Ask the authors whether the 30 uL came from the water phase or a mixed sample (CF-031).
6. **Fat and M1.** F-R2's capture and M1's thiol loss law act on the same pool; build F after M1's verdict, with
   `k_mft_decay`'s status (frozen, banded, free) declared, or the two compensate.
7. **B54's isotherm.** No pea / soy / gluten sorption isotherm or Tg curve is on disk; W-R1 cannot run without one.
8. **"Roos" in the plan row** means Lievonen & Roos 2002 and Miao & Roos 2004 (on disk); neither gives a step rate.

Seen by the author: as stated at the top. No engine output on any row named in sec. 4 was computed for this case.
