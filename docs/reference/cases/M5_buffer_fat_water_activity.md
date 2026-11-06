# Case M5: buffer, fat and water activity

**READY TO FREEZE (items 1 - 5), not frozen** (2026-11-06). All three parts are **development grade on disk**: every
decisive datum on disk is spent or seen, so on disk the data can falsify rivals but credit none. Items 1 - 5 freeze
when the reduced power check of sec. 5.3 has run, its demotions are written in, and the review has answered "What
should block the freeze" (end of file). No engine was run for this case, on any pot. The changelog is at the end.

**The spent-data rule as applied here (coordinator refinement, 2026-11-06; the same rule as M4's FM2).** Spent or
seen data may **falsify** a rival whose prediction there is constant-free (structural: a ratio fixed at 1, a sign, a
bound that holds for every value of the rival's free constants). They never **credit** a rival. The rule applies to
R0 exactly as to the new rivals. Sec. 2.10 lists what is spent and why. Sec. 4 lists, part by part, the
constant-free predictions the spent data falsify, and what is left.

*Plan step M5 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E5 (F01-b, F03-b, F13-b, F14-a,
F19-c, F23-c, F29-b) and E7 (F30-a, F28-d, F26-c) of
[ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, amended in sec. 5 (E5 "SUPPORTED,
step-specific"; E7 SUPPORTED) and sec. 7 (E5's form; E7 / CF-031); the a_w mismatches F05-a, F16-d, F19-d, F25-d,
which carry no E number). It covers M0's group G5 and the M5-owned checks of the
[misfit inventory](../misfit_inventory.md) (sec. 2.11). Items 1 - 5 of sec. 7, written 2026-11-04 and revised
2026-11-06, before any result of the change exists. Items 6 - 9 are appended later and never edited into this part.
Background physics: [PH_BUFFERS_AND_WATER](../../guides/PH_BUFFERS_AND_WATER.md) secs. 6 - 7.*

**Seen by the author.** No engine output on any held-out row named here.
- 2026-11-04 text: measured values of held-out sources read through the conflict register and ENGINE_VS_LITERATURE
  sec. 7: the v2 E5 TEST contrasts F02#1 (Potman 1989), F04#2 (Kallen 1971b), F05#23 / F06#11 / F08#2 (Thornalley
  1999) as quoted in CF-093 / EX-058; the v2 E7 TEST contrast F25#17 (US 8,920,862, MFT x10.1) as quoted in CF-031;
  Larralde 1995's ratio (keyed held-out); the atlas text of F30#11 (Kim 1996b), F30#15 (Zamora 2015) and F30#16
  (Hidalgo 1993), measured values only. Yiltirak 2026 engine numbers quoted in CF-068 and SULFUR_FLUX_MAP were read for
  case M1, not for this one.
- 2026-11-06 revision: the dossiers' printed values for Bell 1997 (Table 1), Zyzak 1995 (Table I), Davidek 2002
  (sec. 4.3 - 4.4 rate table), Elmore 2002 (Table 1 thiol rows and the pot charge), Yiltirak 2026 (Table S3, the
  100 C buffer and buffer + oil rows); CF-011, -031, -093, -094, -104, -120 and PR #88's CF-127 - CF-133; the misfit
  inventory, M0 integrity batch 1 (both halves) and M0 checks batch 2; `trunk_conditions.py` and `engine.py`'s a_w
  refusal (PR #77). The citation sweeps read titles and abstracts of unseen works only. [dry] marks arithmetic on
  printed values, done for the falsify-only units of sec. 4; it is not a result of the change.

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

**Value basis (M1.1).** Every unit in sec. 4 is within one study: a ratio, a within-study contrast or a sign. H0's
between-laboratory term (EX-059) enters only on absolute levels compared across laboratories, so it enters no M5 unit.
Cross-laboratory contrasts (MF-05's Xie vs Jiangnan; D5 across Reading and CAAS) are development evidence only.

## 1. The claim

Three conditions that every lab constant carries silently, each turned into a term on named steps.

| part | step(s) | today (B60, read from the code) | M5 |
|---|---|---|---|
| B. buffer catalysis | Amadori decay (`k_ama_*`, `r_arp_*`, the Cys-Amadori steps if M1's C1 ships), the 2,3-enolisation branch, `r_schiff` (lumped), TTCA formation (CF-011) | no buffer term anywhere: the buffer is a titrant in the charge balance only ([ph_state.py](../../../src/kinetic_core/ph_state.py)); every constant carries its calibration pot's phosphate (Hofmann 0.5 M, Martins 0.1 M, Cao 0.2 M) | k_s = k_s,u (1 + sum_b kappa_{s,b} c_b) on steps the network resolves, **kappa by step and by buffer species** (H2PO4- vs HPO4 2-, carboxylate); a hyperbolic factor on a lumped step whose uncatalysed part can limit; carboxylates and free-sugar steps on their own bands (ENGINE_VS_LITERATURE sec. 7) |
| F. fat | MFT, FFT (and their disulfides), Cys and H2S against the lipid lane's carbonyls and the oil phase | lanes integrate as a direct sum (`lipid.lane_coupling_verdict`, [lipid.py](../../../src/kinetic_core/lipid.py) ~784 - 830); no oil phase, no thiol partition anywhere in `src/`; the lipid lane's aldehyde-protein `COVALENT_SINK` is disabled by its measured Ea ([parameters_lipid.py](../../../src/kinetic_core/parameters_lipid.py) ~551 - 589); the lane charges a lipid carrier as 100 % of the pot (CF-130, pending on PR #88) | the plan's oil-water partition of the thiols, against the rivals the register names: lipid-electrophile capture of the MFT pool (CF-031), oil-phase oxidation (M0 group G9), and both |
| W. water activity | `k_schiff` (trunk), the sulfur lane's steps | trunk: Pereyra Gonzales's multiplier table on `k_schiff` only (`AW_STEPS`, [trunk_conditions.py](../../../src/kinetic_core/trunk_conditions.py) ~60 - 100), **inside the window a_w 0.33 - 1.00** (PR #77): a comparison with an arm below 0.33 is refused, a single pot below it is answered with the multiplier held at the 0.33 value and a glass warning; acrylamide: a flat band (CF-104); sulfur, lipid and pyrazine observables refuse every a_w comparison (`engine.axis_refusal`) | B54's design: concentration per whole phase, the trunk table re-derived as a residual, a Gordon-Taylor Tg gate; a step-specific term in the sulfur lane only where a sign contrast demands it |

**Earlier waves on these structures** (cited, none re-run):

| structure | wave | outcome | kind |
|---|---|---|---|
| per-step phosphate catalysis | B53 (backlog; PH_BUFFERS_AND_WATER sec. 6) | planned, deferred 2026-10-27 (no blind test), unblocked 10-28 (D5), **never run**. Fit source Davidek 2002's phosphate series; hold-out Mottram & Nobrega 2002 | -- |
| the factor's form | EX-058 (resolve2_pyrbuf, red-teamed) | non-saturating to 0.5 M, order 1.0 - 1.35 (P 0.8): **decided on Potman (E5 TEST)**, with Kallen and Thornalley (E5 TEST) and Larralde (keyed held-out). These are now SPENT (misfit inventory MF-26) | arithmetic on data |
| buffer in the Hofmann bundles | B51 data facts (Amendment 44) | 0.5 M phosphate stated on every arm (a titrant fact) | data correction |
| phosphate on ARP decay | CF-120 (M0 CK-55, done) | Deng 2021's 0.1 M arm, extrapolated: weak support at pH 7.5 only | arithmetic |
| buffer in the Schibilsky pots | M0 CK-06 (batch 1) | the pots are unbuffered and pass their pH to the engine: MF-28's buffer candidate (c) is inert there | integrity |
| water activity basis and Tg gate | B54 (backlog) | planned, **never run**; needs a pea / soy / gluten isotherm and Tg curve | -- |
| trunk a_w window | PR #77 (H0 follow-up, owner-approved 2026-11-05) | window 0.33 - 1.00; refusal of comparisons below it | declaration |
| acrylamide a_w | CF-104 red team (2026-11-03) | flat multiplier upheld (P 0.75); band-drop fix applied | review |
| lipid x sulfur | CF-031 cluster 4 (2026-11-01) | a partition operator explains at most 0.1 - 0.4 dex of the 0.5 - 2.2 dex falls (P 0.55) | arithmetic |
| lipid x sulfur | resolve2_sulfur (CF-031) | the pre-registered speciation rule FAILED; partition not excluded | arithmetic |
| aldehyde-protein covalent sink | B8 / K6b | disabled by a measured Ea (15 - 23 kJ/mol) | measured constant |
| lipid matrix pots | M0 integrity batch 1 (MF-42), filed as CF-129 / CF-130 on PR #88 | carrier charged undiluted (Trikusuma x10.4 too much fat; Li SPI 0.258 w/w); Bi 2020's roasted hexanal is a net 3.9x loss of a carried level, now refused | integrity |

**Errors stack (CF-107).** No whole-model fit has ever tested a buffer, fat or a_w term, so no refusal needs
discounting. The converse risk is real: every sulfur constant was fitted in 0.5 M phosphate, so adding part B
re-references them, and the joint refit (JR) must carry it. The same holds for M3's fitted ARP-Phe decay (`k_d`, 0.5 M
phosphate; MF-32) and for M1c's `r_cama_form` if it is fitted on Xie's phosphate pots (MF-05). Candidate errors that a
buffer term could absorb: E2 (the split frozen in pH), CF-084 (the pH balance), the Hofmann pH 3 MFT miss
(SULFUR_FLUX_MAP sec. 5 names phosphate on the enolisation as one candidate). For fat: `k_mft_decay` (frozen,
pot-independent; refuted as too fast at 110 C by M0 CK-23) can mimic any MFT-selective loss, and any oxidation term
shares the vessel's O2 with CF-068. For matrix pots, CF-130's undiluted carrier inflates both the oil volume and the
lipid lane's carbonyl output.

## 2. All the evidence, by kind

The *spent* column says why a datum cannot credit a rival (sec. 2.10 has the rule and the record for each).

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Du 2025 (F30#8) | CAMOLA: 2-methylthiophene 100 % glucose carbon without fatty acid, 38 - 48 % with it | lipid carbon enters sulfur heterocycles (F) | DEV(E7) | CF-031 |
| Du 2023 (F30#9) | 2-pentylthiophene only from the 2,4-dienal, ~84 % without glucose carbon | H2S + dienal capture (F) | **v2 TEST(E20)**: listed, not used | -- |
| Zhao 2019b Table 1 / text (F23#7) | label patterns unchanged with chicken fat | fat does not reroute MFT carbon | keyed held-out (E7): item 7 | -- |
| Whitfield 1988, Elmore 2002 (F26#21, #24, #25) | alkylthiophenes / thiapyrans from PUFA + H2S | the sulfur sink into lipid products | DEV(E7) | CF-031 |
| none | no label on buffer or a_w routes | -- | gap | -- |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Davidek 2002 (F03#0 - #2) | fed DFG 0.1 M, 90 C, **pH-stat at 90 C**, water vs 0.1 M phosphate. Decay ratio k_P / k_W (band = spread of the dossier's fitting methods): **2.79 - 3.11 (pH 5), 3.35 - 3.62 (pH 6), 1.94 - 2.17 (pH 7), 1.24 - 1.27 (pH 8)**; the ratio tracks [H2PO4-], not [HPO4 2-] (dossier sec. 4.4, at an assumed pKa2). Corrected 2026-11-06: the case and plan row read "x2.8 - 3.6 at pH 5 - 7"; that range is pH 5 - 6 | **buffer catalysis at held pH** (D1); the species (part B's kappa by species); B53's fit source | DEV(E5) | atlas F03-b; MF-26 |
| Zyzak 1995 (F03#7 - #9, F05#22, F06#10) | fed fructoselysine, 37 C, pH 7.4, N2 + chelator: t1/2 57.1 / 40.9 / 24.4 / 16.8 d (+- 2.0 / 2.3 / 2.1 / 0.6, triplicates) at 10 / 50 / 100 / 200 mM phosphate; HEPES ~12x slower; products unchanged | dose law on the Amadori decay step (D2) | DEV(E5) | atlas F03-b; CF-093 |
| Guan 2012 (F01#9, F03#10) | fed proline ARP, 110 - 130 C, pH 5.8: phosphate x3.8 (110 C) / x1.5 (130 C); acetate, citrate +16 / +34 % | carboxylate kappa; kappa falls with T (D4) | DEV(E5) | atlas F01-b; CF-094 |
| Hofmann 1998 Table 2 dry arm (F24#9) | fed C2 + C3: dry 180 C MFT 1554 vs aqueous 268 ug (T, t, a_w, buffer confounded) | a_w on one fed step | DEV | atlas F25 |
| Elmore 2002 (F30#0 - #1, F25#16, F28#19) | Cys 25 + ribose 25 mM, pyrophosphate pH 5.5, 140 C, +- 25 mM methyl linoleate / linolenate (an immiscible ester phase, ~0.8 % v/v; dossier), n = 4: MFT 807 (508) -> 58 (71) -> tr; FFT 92 (54) -> 87 (32) -> 17 (5); bis(2-methyl-3-furyl) disulfide 202 -> tr -> nd; C3 - C5 2-alkenals consumed | MFT-selective loss with PUFA (D5); unsaturation at equal ester volume (D6) | DEV(E7) | CF-031; MF-18 |

### 2.3 Trapping

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Esterbauer 1975, Doorn 2002 (F28#17, F30#18 - #19) | thiol + enal Michael adducts, reversible in hours (2-alkenals); thiol >> His > Lys | the capture chemistry exists and reverses | **v1 TEST**: listed, not used | -- |
| Charles-Bernard 2005 (F28#0, #16) | 25 C: carbonyl-competing nucleophile slows thiol loss 4 - 13x | carbonyl capture of thiols (low T) | DEV | seen, not spent |
| none at process T | no trapped lipid-thiol adduct at >= 100 C | F | gap | -- |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Cao 2017 (F04#9, F04#17) | 0.2 M phosphate pH 5.5, 90 C: TTCA K flat while pH drifts | the phosphate-catalysed exchange (CF-011) | DEV | CF-011 |
| Zhang 2022 (Hefei; no atlas cell; CF-031) | Xyl + Cys, pH 7.2 unbuffered, 120 C, 90 min, 4 : 1 oil : water, phases separated: water-phase free MFT x0.36 - 4.8; the control's 2.81 ng/mL disulfide vanishes from raw-oil water phases; MFT and disulfide found in the oil | **both species partition into oil** (D7); n, SD not printed | no frozen split | resolve2_sulfur |
| Tang 2019 (F02#9 - #10) | GSH + xylose, 80 C, vacuum dehydration: ARP x8 in a tenth of the time | water removal speeds Amadori formation | DEV | seen, not spent |
| Kang 2026 (F04#25) | TTCA powder, a_w 0.11 / 0.43 / 0.84, 60 d: retention 86 / 72 / 64 % | a_w on the cysteine reservoir (storage) | DEV; calibration input | calibration |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Farmer 1992, Xu 2011 (F30#2, F30#20; CF-088) | lipid aldehydes fall with Cys + pentose while saturated ketones do not: capture, not antioxidant | the reverse flow (lipid lane loses carbonyls) | Farmer DEV(E7); Xu item 7 | CF-031, CF-088 |
| Potman 1989 (F02#0 - #1) | ARP carries all glucose lost at <= 1 h; phosphate order 1.34, ~x15 at 1.2 M | the dose law | **v1 TEST, v2 TEST(E5)** | EX-058 (SPENT; falsify-only) |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split | spent |
|---|---|---|---|---|
| Bell 1997 (F01#7 - #8) | Glc + Gly 0.1 M, pH held 6.86 - 7.05, 25 C: glycine loss 0.00055 +- 0.00050 / 0.00082 +- 0.00039 / 0.0019 +- 0.0003 / 0.0033 +- 0.0009 M/month (95 % CL; one vial per time point) at 0.02 / 0.05 / 0.2 / 0.5 M phosphate; citrate rates not different from zero | **dose law at held pH** (D2); carboxylate (D4) | DEV(E5) | atlas F01-b; CF-094 |
| Mottram & Nobrega 2002 (F23#20, F13#10, F14#8) | Cys + ribose 50 / 50 mM, 140 C, 30 min: pyrophosphate or phthalate pH 5.6 vs unbuffered (5.6 -> 3.9): MFT x2.2 - 2.6, FFT unchanged; ribose-only furfural 107 vs 103 | **step specificity** (D3); pH-confounded (unbuffered drifts) | DEV(E5) | atlas F13 / F14 / F23; MF-26 |
| Blank 1997 (F13#11) | HDMF from xylose + Gly, 90 C: water / 0.2 M phosphate / 0.2 M malonate 0.06 / 2.6 / 0.6 ug/mmol | phosphate vs carboxylate on a 2,3-enolisation product | DEV(E5); calibration input | atlas F13 |
| Rizzi 2004 (F14#4) | ribose alone, 0 / 0.05 / 0.10 M phosphate vs bis-tris, pH 7.3: A420 0.003 -> 0.045 -> 0.262 | phosphate drives an amine-free sugar step (browning, not a resolved step) | DEV(E5) | atlas F14; CF-094 |
| Adams 2008 (F19#6) | fed MGO + Ala, pH 7, 130 C: phosphate vs bicarbonate DMP x1.8, pyrrole x0.36 | buffer on pyrazine steps (M6) | DEV(E5) | atlas F19 |
| De Vleeschouwer 2006 (F29#6) | acrylamide k_F citrate vs phosphate x0.14 (pH 4), x2.4 (pH 6) | carboxylate vs phosphate by pH | DEV(E5); calibration input | atlas F29 |
| Lin 2021 ("Chang 2021"; MF-37) | Glc + Asn, 180 C, pH set 6.0: 1 % acetic acid vs water, acrylamide 1459 vs 832 ppb at 30 min; heat-up unprinted (M0 batch 1: integrity holds) | a carboxylate arm on a fitted lane; in-situ pH not stated | DEV | MF-37 (CK-04) |
| Yiltirak 2026 (F30#6; dossier Table S3) | ribose + Cys 25 mM, 0.5 M phosphate pH 5.5, four equivalent cooks, n = 3: at 100 C / 4 h buffer MFT 6.88 +- 0.98, FFT 1.28 +- 0.08; buffer + 10 % canola oil MFT 0.72 +- 0.10, FFT 1.19 +- 0.06; emulsion MFT 3.73; FFT x0.57 at 120 - 130 C | **fresh-oil MFT-selective fall** (D5) | DEV(E7); the buffer arm is a Maillard-path hold-out bundle, the oil arms are not | CF-031; MF-18 |
| Whitfield 1988 (F30#3, F25#15) | lecithin: MFT x1/3.0, FFT x1/2.0, methylpyrazine x1/2.1 | moderately selective; a sulfur-free product falls too | DEV(E7) | CF-031; MF-18 |
| Farmer 1990 (F25#14) | four lipids: FFT x0.62 - 0.72, furfural x0.49 - 0.80 | FFT side | DEV | CF-031 |
| Du 2025 (F30#7) | GSH + glucose 200 / 200 mM, 145 C, +- 10 mM free C18 acids: MFT +19 - 30 % | the small opposite sign (CF-031) | DEV(E7) | CF-031 |
| Xu 2011 (F30#5, F28#20) | fresh / oxidised lard: MFT x1/21 / x1/143, FFT x1/3 / x1/6.4 | oxidation dose | item 7 only (E12 TEST) | -- |
| Trikusuma 2020, Li 2026, Bi 2020 (MF-42; CF-129, CF-130 on PR #88) | matrix lipid-aldehyde rows: carriers 3 % w/w PPI at 8 % fat (Trikusuma) and SPI 0.258 w/w (Li) are charged as 100 % of the pot; Bi's roasted hexanal is a net 3.9x loss of the raw level | the oil volume and the carbonyl load in food pots (F); no thiol measured | DEV | MF-42 (CK-05) |
| Hofmann 1998 Table 2 dry arm (F25#11 - #12, F23#18) | dry 180 C / 6 min vs aqueous 145 C / 20 min: ribose FFT x8.0, glucose and rhamnose FFT x0.50; MFT x1.3 / 2.2 / 3.9 | **a_w sign contrast across sugars** (D8); T, t, buffer confounded but shared by the sugars | DEV | atlas F25 |
| US 10,076,126 (F25#13) | glycerol vs buffer, 110 C: FFT up to x56 | a low-water solvent on FFT | DEV | atlas F25; its MFT shape is a B59b hold-out |
| Granvogl 2006 (F16#8) | silica, 10 % water: Leu / Ile aldehyde x2 - 3; Strecker amine 19 - 27 % | a low-water exit (M3) | DEV | atlas F16 |
| Leahy 1989a (F19#11) | a_w 0.32 / 0.58 / 0.75 / 0.84: pyrazine rate 1 : 2.10 : 4.23 : 4.07 | the a_w maximum (D9) | DEV; calibration input (not B18 fit rows); a B54 hold-out | atlas F19 |
| Scalone 2015 (F19#15) | a_w 0.33 -> 0.75: 2,5(6)-DMP -0.3 dex | a_w on pyrazines | DEV | atlas F19 |
| Modi 2008 (no atlas cell) | Cys + thiamine + ribose in beef fibre, 130 C: a_w 0.60 / 0.80 / 0.90; peak areas, n = 4 | **the only a_w series for MFT** | no frozen split | directional panel MOD08-AW-01 (direction read; "not evaluable") |
| Bell 1995 (no cell) | browning in PVP at fixed internal molality: where the three polymers straddle the glass at one a_w, the glassy one is slowest (the author's "7-fold"); browning only, not a step | concentration basis (D9) | Auburn, a DEV(E5) lab | B12 (sets the trunk band's low edge) |
| Eichner 1972 (no cell) | 37 C browning maximum at a_w 0.41 - 0.52 (x38 from 0.84 to 0.41); end points only | a_w maximum (low T) | no frozen split | seen |
| Lievonen 2002, Miao 2004, Buera 1993 (no cells) | (T - Tg) not sufficient; Ea depends on water content; WLF fits | the Tg gate's limits | no frozen split | B54 (motivated the gate) |
| Tang 2024 (F05#17) | dry a_w 0.6 vs wet: dry glycates far more; pea / soy order flips | protein glycation vs regime | DEV | atlas F05 |
| De Vleeschouwer 2007 / 2008 (F29#7 - #8) | acrylamide formation flat over a_w 0.34 - 0.99 | concentration basis on a fitted lane | DEV; calibration input | CF-104; B54 |
| Warmbier 1976, Pereyra Gonzales 2010 (F05#18 - #21) | lysine loss maximum 0.5 - 0.6; flat 0.31 - 0.71, x4 slower at 0.98 | the shipped table | **v1 TEST**: listed, not used (Pereyra is already the shipped table) | calibration |
| Hartman 1984b (F27#20) | sulfur volatiles peak near a_w 0.72 | a_w on thiamine sulfur | **v2 TEST(E18)**: listed, not used | -- |
| Kallen 1971b (F04#2), Thornalley 1999 | saturating thiazolidine catalysis; convex 10 - 100 mM | the form | **v2 TEST(E5)** | EX-058 (SPENT; falsify-only) |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split | spent |
|---|---|---|---|---|
| Bassama 2011 (CF-104; MF-41) | plantain, a_w 0.97 / 0.90 / 0.43: k_F x8.1, k_E x14; wet-basis mass action fits 0.43 / 0.90, misses 0.90 / 0.97 by x2.6 (CI 1.9 - 3.6) | concentration basis in a food | no frozen split | CF-104 |
| Martins 2005, Knol 2005 | every constant measured in 0.1 M phosphate | what part B re-references | DEV; calibration inputs | calibration |
| engine records | B51 (Amendment 44), EX-058, CF-031 clusters, CF-104 | sec. 1 table | -- | -- |

### 2.8 Left out, with reasons

Lipid Strecker and dienal + NH3 pyridines (Zamora, Hidalgo, Kim 1996b; F18, F21): E20. Lipid oxidation kinetics
themselves (Frankel, Schroen): the lipid lane, not its coupling. Matrix binding of aldehydes to proteins (Snel, Bi): the
matrix layer. Thiamine a_w (Hartman): E18. MF-42's mechanism (the lipid lane's under-prediction in matrices): the
lipid lane's own case, after CF-130's carrier fraction is built; this case takes from it only the oil volume and the
carbonyl load (sec. 2.11).

**Counts by kind (sources).** Labelling 4 (one held out); fed intermediates 5; trapping 3 (two held out);
measured intermediates 4; mass balance 2 (one held out); one factor within a study ~32 (six held out); fitted 2
plus records.

### 2.9 The citation sweeps (sec. 7 item 2)

**Run 2026-11-06** with `scripts/tools/openalex.py sweep --direction both` through the wrapper (anonymous, ~1 request
per 1.1 s). Screening used titles and abstracts only, so blinding holds. Each sweep was keyed to its part's decisive
criterion; a title was kept for reading if it named the criterion's terms, otherwise excluded by category.

**Criteria.**
- **B (buffer):** a buffer species or concentration varied at a held or stated in-situ pH on a Maillard step
  (Amadori formation or decay, enolisation, dicarbonyl, furfural, MFT / FFT); a dose ladder of >= 3 concentrations;
  two branches read in one buffered vs unbuffered pot; a carboxylate vs phosphate arm at equal pH; a buffer effect in
  a low-moisture or concentrated matrix.
- **F (fat):** MFT / FFT (or another thiol) with and without a lipid in one pot; lipid unsaturation or oxidation at
  equal amount; a dosed lipid carbonyl without an oil phase; a phase-resolved thiol or a measured oil-water partition
  coefficient for MFT, FFT or the disulfide; a lipid-carbonyl + thiol / H2S adduct at >= 100 C.
- **W (water activity):** a Maillard step rate or a thiol at >= 2 a_w (or moistures) in one study with the
  concentration basis stated; a Tg-resolved step rate; a dry vs wet sign contrast across sugars; a sorption isotherm
  or Tg-moisture curve for pea, soy or gluten protein (B54's need).

**Seeds** (forward and backward each):
- B: Davidek 2002, Bell 1997, Zyzak 1995, Guan 2012, Mottram 2002.
- F: Elmore 2002, Yiltirak 2026, Whitfield 1988, Zhang 2022.
- W: Hofmann 1998 (the quantitative precursor paper), Leahy 1989a, Modi 2008, Bell 1995.

| sweep | works | held | unseen | biomedical / other field | review, book, chapter | computed (excluded by rule) | title outside the criterion (food product, other chemistry) | read by title / abstract | candidates |
|---|---|---|---|---|---|---|---|---|---|
| B | 502 | 79 | 423 | 92 | 8 | 4 | 250 | 69 | 6 (C-B1 - B5, C-W5) |
| F | 409 | 50 | 359 | 25 | 8 | 0 | 163 | 163 | 5 |
| W | 286 | 58 | 228 | 19 | 7 | 0 | 127 | 75 | 4 (C-W1 - W4; C-W2 is C-B2) |

Fourteen distinct candidates in all.

**Read and excluded, by reason.**
- B (63 of 69): Amadori analytics, synthesis or preparation (~25); in vivo or biomedical glycation (~13); reviews and
  chapters (~9, incl. Rizzi 2010's browning chapter, whose 2004 paper is held); food products or matrices (~7, incl. the
  2026 pork AGE and SAPP / STPP glycation papers, where the phosphate salts move the pH and no held-pH contrast is
  stated); computed (3); non-Maillard degradation in buffers (aspartame, rebaudioside: 3); pathway work outside the
  criterion (Strecker from ARPs, retro-aldol labels: passed to M3 / M4).
- F (158 of 163): meat, fat and oil product studies with no defined thiol pot (~110); emulsion and oil physics without
  a Maillard thiol (~15); Maillard x lipid systems with no thiol (glycine + ribose + phospholipid, glucose + lysine +
  linoleate, sesame-oil aroma: ~12); reviews and chapters (~12, incl. Mottram & Elmore 2005, whose data are Reading's
  held papers, and Sagalowicz 2016, a review from the Nestle group whose patent is an E7 TEST source); qualitative
  identification of sulfur volatiles (~8); 2-alkylfuran formation from enals with amino-acid catalysis (Adams 2011:
  passed to the lipid lane's MF-42 / CK-62, not this criterion).
- W (71 of 75): browning kinetics only, product or glass physics without a step rate (~40); water-activity reviews and
  chapters (~15); roasting product aroma (~10); single-a_w or ratio studies from the Minnesota group (Warmbier 1976b:
  one a_w; a v1 TEST laboratory in any case).

**Candidates, with roles fixed before reading** (download needed; the DOIs went to the coordinator with the case's report, not into the repository, which is
public). A candidate that
is obtained enters by amendment before the freeze, read under the role below; after the freeze it enters only by
amendment declared before its data are read. Its laboratory is checked against the splits first.

| id | work (titles / abstracts only) | part | role fixed now | lab check |
|---|---|---|---|---|
| C-B1 | Pepa, dos Santos-Ferreira & Buera 2023, ACS Food Sci. Technol.: organic and inorganic anions on glucose-glycine, 100 C, near-neutral pH; constants correlate with anion basicity | B | D4 (carboxylate vs phosphate by species) and a Bronsted-law sub-rival of R1, **only if** a resolved-step rate (glucose or glycine loss) is printed per anion at equal pH; otherwise a band source | Buenos Aires (Buera): confirm it is not the v1 TEST laboratory of Pereyra Gonzales 2010 |
| C-B2 / C-W2 | Bell, White & Chen 1998, J. Food Sci.: phosphate and citrate concentration in glassy low-moisture solids, pH 7, 25 C; glucose-loss constants **fall** as buffer rises, browning rises | B, W | the concentrated-matrix ceiling (CF-093's third explanation): a sign contrast against the dilute law, decisive for "the dilute factor transfers to low-moisture foods" | Auburn, a DEV(E5) lab |
| C-B3 | Zhang & Serianni 2012, JACS: phosphate-catalysed glucosone degradation, pH 7.5, 37 C, 13C NMR | B (and M4's R5) | a resolved dicarbonyl step with phosphate catalysis: tests R4's "Amadori-side steps only" by a structural sign | -- (spends for M4 too) |
| C-B4 | Mikami & Murata 2015, Food Sci. Technol. Res.: acetate / phosphate / Tris, xylose- and glucose-lysine, 100 C, sugar loss measured | B | D4 band source only (browning is not a step; buffers at stated pH) | -- |
| C-B5 | Dobbs & Bell 2009, Food Res. Int.: tagatose stability in buffers of several compositions | B | the free-sugar factor rho of CF-094 (an amine-free step), band only | Auburn |
| C-F1 | Liu et al. 2023, Foods: nine aliphatic aldehydes (saturated and unsaturated) dosed into glutathione + ribose | F | **decisive for F-R2 vs F-R1 and F-R6**, if MFT / FFT are quantified: a dissolved carbonyl with no oil phase. Capture predicts a thiol fall that rises with electrophilicity (enal > alkanal); partition and oil-phase oxidation predict none | -- |
| C-F2 | Van Ba et al. 2013, Meat Sci.: particular unsaturated fatty acids x pH in meat-like model systems | F | D6 (unsaturation at equal amount), if the thiols and the lipid amount are printed | -- |
| C-F3 | Yang et al. 2015, Flavour Fragr. J.: Cys + glucose +- chicken fat +- oxidised fat, 12 odorants quantified | F | D5 / D6 at equal fat amount (fresh vs oxidised): the oxidation dose | Beijing (BTBU): check against the E7 split |
| C-F4 | Elmore & Mottram 2000, JAFC: (E,E)-2,4-dienals + H2S at 140 C give thiapyrans and alkylthiophenes | F | fed intermediate for F-R4 (the H2S sink): a band, not decisive | Reading, DEV(E7) |
| C-F5 | Yeo et al. 2026, Food Chem.: chicken triglycerides and phospholipids added to a defatted matrix; some S-containing products rise with phospholipid | F | evidence against a uniform lowering (development only; open SDE system) | Reading, DEV(E7) |
| C-W1 | Bell, Touma & White 1998, J. Food Sci.: glycine loss and browning vs a_w and the glass in PVP, 25 C; glycine loss falls on matrix collapse | W | the Tg gate on a step (glycine loss ~ `k_schiff`) vs concentration only (W-R1 vs W-R2), trunk lane, low T | Auburn |
| C-W3 | Sherwin & Labuza 2003, J. Food Sci.: glycerol / sorbitol plasticisation vs reactant dissolution in intermediate-moisture browning | W | concentration vs mobility; **item 7 only** if obtained (Minnesota is a v1 TEST laboratory) | Minnesota (TEST) |
| C-W4 | Bell & Labuza 1994, J. Food Eng.: the low-moisture state's pH and its kinetic implication | W, B | R3-for-foods (the in-situ pH of a dried buffered matrix): a band | Minnesota / Auburn: check |
| C-W5 | Deng et al. 2025, Food Chem.: low-moisture synthesis of the Met-Glc Amadori compound, pH x temperature | W | a_w on Amadori formation / decay: a band | Jiangnan |

**Gaps the sweeps confirm.** No work measures an oil-water (or oil-air) partition coefficient of MFT, FFT or the
disulfide. No work gives a pea, soy or gluten sorption isotherm or Tg curve. No calibrated thiol series over a_w.
These are the decisive measurements proposed in the PR body (sec. 4.4).

### 2.10 Spent data (sec. 7 item 2)

**Spent: these data nominated a rival or a part, and are development evidence only.**

| data | spent by | spent for | consequence |
|---|---|---|---|
| Davidek 2002, Bell 1997, Zyzak 1995, Guan 2012, Mottram 2002, Blank 1997, Rizzi 2004, Adams 2008, De Vleeschouwer 2006 | the atlas cells F01-b, F03-b, F13-b, F14-a, F19-c, F23-c, F29-b that nominated E5; CF-093 / CF-094 tests; M0 MF-26 (Mottram, Davidek) | part B (the claim and every buffer rival) | **falsify-only** (sec. 4.1) |
| Potman 1989, Kallen 1971b, Thornalley 1999, Larralde 1995 | EX-058 (resolve2_pyrbuf) chose the form on them (inventory MF-26: SPENT) | B's form | falsify-only; item 7 still lists them for the side-effect table |
| Deng 2021, Martins 2003 | CF-120 (M0 CK-55, done) | B on ARP decay | development only |
| Cao 2017; the Xie vs Jiangnan CAMA contrast | CF-011; M0 MF-05 | B on Cys-Amadori formation | development only (cross-laboratory) |
| Lin 2021 acetic vs water | M0 MF-37, CK-04 | -- | development only |
| Elmore 2002, Yiltirak 2026, Whitfield 1988, Farmer 1990 / 1992, Du 2025, Zhang 2022, US 8,920,862 | CF-031 cluster 4 and resolve2_sulfur nominated capture (P 0.55) on them; M0 MF-18 | part F (capture nominated; partition bounded) | **falsify-only** (sec. 4.2) |
| Trikusuma 2020, Li 2026, Bi 2020 | M0 MF-42 (CK-05); CF-129 / CF-130 | the lipid lane's matrix rows | not F evidence (no thiol); they set the oil volume and load in food pots |
| Hofmann 1998 dry arm, US 10,076,126, Granvogl 2006, Leahy 1989a, Scalone 2015, Tang 2024, Warmbier 1976, Pereyra Gonzales 2010 | the a_w mismatches F05-a, F16-d, F19-d, F25-d that nominated part W; Pereyra and Bell 1995 also set the shipped table and band (B12) | part W | **falsify-only** (sec. 4.3) |
| De Vleeschouwer 2007 / 2008, Bassama 2011 | CF-104 (acrylamide cluster); M0 MF-41, CK-61 done; B54's "whole-phase concentration explains Pereyra and De Vleeschouwer" | W's concentration basis | falsify-only |
| Lievonen 2002, Miao 2004, Buera 1993 | B54 (they motivated the Tg gate) | W-R1's gate | band and limits only |
| Modi 2008 | the directional panel read its a_w direction (MOD08-AW-01, "not evaluable") | -- | seen; falsify-only |

**Seen but not spent.** Charles-Bernard 2005 and Tang 2019 were read for this case and nominated no rival. Neither
tests a constant-free prediction that separates the surviving rivals (sec. 4), so neither is decisive.

**Reserved, not read for values:** the E7 TEST contrasts F30#11 / #15 / #16 (unseen; outside the thiol-loss
chemistry, so they could not decide F's rivals in any case), Xu 2011 (E12 TEST), the keyed held-out tables, the
B59b hold-outs (US 10,076,126's MFT shape; US 8,920,862's level), Hartman 1984b (E18 TEST).

### 2.11 M0's groups: what an M5 rival predicts for their members on seen data

Directions are structural or [dry]. No engine was run.

**G5 (a step-specific buffer factor: MF-05, -26, -28, -32, -37; inventory sec. 4).**

| member | what a part-B rival predicts there | reading |
|---|---|---|
| MF-26 (no buffer catalysis; Mottram, Davidek, Deng) | every rival but R0 and R3 moves the axis; their sizes are fitted | M5's own misfit; spent; CK-54 (owned by M5) runs inside the prereg, falsify-only |
| MF-05 (CAMA formation, Xie phosphate pots 3 - 6x faster than Jiangnan water pots) | R1 with kappa_P in [5, 40] /M on `r_cama_form` at 0.2 M gives x2 - 9, covering 3 - 6x; R4 (Amadori decay and 2,3-enolisation only) gives x1 | **a structural R1 / R4 split in principle, confounded one-to-one with laboratory and quantification** (CF-011's (b), P 0.3): development only. CK-33 has no seen-data discriminator. Decider: CF-011 buffer arm C |
| MF-28 (Schibilsky furanics pH-flat; candidate (c) HDMF phosphate x4.3) | CK-06 found the pots **unbuffered**: every part-B rival predicts nothing there | candidate (c) is retired for MF-28; MF-28 stays with G2 (the enolisation pH law, CK-21) |
| MF-32 (ARP-Phe fed vs in-pot rates; CF-120) | in 0.5 M phosphate every pot carries the same factor, so B moves none of MF-32's within-laboratory contrasts | B re-references M3's fitted `k_d` at compose time (sec. 6); a JR item, not an M5 unit |
| MF-37 (Glc + Asn 180 C; candidate (b) the phosphate pH law applied to every buffer) | Lin 2021's arms are water vs 1 % acetic acid, not phosphate; Ye 2024 is 0.1 M phosphate. R1 by species replaces the acrylamide lane's phosphate-for-every-buffer law; the misfit's shape (x52 from 10 to 30 min) is a heat-up lead (CK-04), which B cannot move | no B prediction on the shape; the species requirement is recorded for sec. 6 |

**Verdict on G5:** no part-B rival predicts a member beyond its own misfit at a size cited constants fix. The group
gives one structural split (MF-05: R1 vs R4 on `r_cama_form`), confounded by laboratory, and one retirement (MF-28 (c)).

**Fat: MF-18 (lipid inert on the thiols), MF-42 (lipid aldehydes low in matrices), G3 (the thiol sink), G9 (vessel O2).**
- **MF-18** is part F's own misfit (spent). CK-29 (engine-free: does a capture first order in the printed carbonyl
  load fit Whitfield, Elmore and Yiltirak alike?) is owned by M5. It runs inside the prereg as a falsify-only unit: it
  can falsify "first order in the printed load", never credit capture. Its sharpest point is visible in the design
  already: Yiltirak's fall is with **fresh** oil, where a carbonyl load is small unless it forms in situ.
- **CK-37** (the O2 : (Cys + pentose) regression of CF-068) reads Yiltirak, M5's D5, so its Yiltirak part runs here.
  It nominates rival **F-R6**: an oil phase can carry O2 and hydroperoxide, and MFT oxidises faster than FFT. F-R6 shares
  G9 (the vessel and headspace O2 operator) and is listed in sec. 3.
- **G3:** CK-23 refuted the shipped `k_mft_decay` as too fast at 110 C and nominated a cysteine-slowed loss (M1's C2).
  Any F rival acts on the same pool. F's build waits for M1's verdict, with the sink law declared (open question 6).
- **MF-42 and PR #88 (CF-129, CF-130).** No F rival explains MF-42: the matrix pots carry little free thiol, so capture's
  reverse flow (aldehydes lost to thiols, CF-088) is negligible there. After the carrier dilution the misfit is deeper
  (Li -0.86 to -3.15 dex; Trikusuma's decadienal +0.21 -> -0.81) and belongs to the lipid lane. What CF-130 does to part
  F:
  - **F-R1 (partition)** needs the oil volume fraction phi of each food pot. CF-130 says the engine would take phi from
    a carrier charged as 100 % of the pot: x10.4 too much fat in Trikusuma's beverage, x1 / 0.258 in Li's extrudate.
    Partition there would be overstated by up to that factor.
  - **F-R2 (capture) and F-R6 (oil-phase oxidation)** need the carbonyl or oxidant load. In food pots that load is the
    lipid lane's output, which is mis-charged (CF-130) and, once corrected, too low (MF-42). So **no F rival may take its
    load from the lipid lane in any matrix pot** until CF-130's carrier mass fraction is built and MF-42 re-derived. On
    the decisive model pots (Elmore, Yiltirak, the candidates) the load is the printed lipid charge or the measured
    aldehydes. This changes no decisive datum, because none of F's decisive pots is a matrix pot.
  - **CF-129 (Bi 2020).** The roasted hexanal is a net 3.9x loss of a carried level in an open oven, so it is mostly
    volatilisation. Capture predicts a negligible share of that loss at a pea's free-thiol level. The row stays refused
    and is not F evidence.
  - **Lever table.** The W5.4 Th x fat cells need phi from the corrected carrier fraction. Until it exists, F's terms
    may not be evaluated on a food recipe. This is a refusal to declare, not a number.

**Water: MF-41 (acrylamide a_w, Bassama) and PR #77.** CK-61 is done, with no seen-data discriminator left. Bassama
is a falsify-only unit for W-R2 (sec. 4.3). PR #77 sets what R0 answers. Trunk observables are answered inside a_w 0.33
- 1.00. A comparison with an arm below 0.33 is refused. Sulfur, lipid and pyrazine observables refuse every a_w
comparison. So R0 refuses on most of W's units (sec. 4.3), and item 6 must say how a refusing R0 is scored (a freeze
blocker).

## 3. The rivals

| id | B. buffer | F. fat | W. water activity |
|---|---|---|---|
| R0 | **current engine**: none | direct sum, lipid inert | trunk table on `k_schiff` inside a_w 0.33 - 1.00 (PR #77); refusal elsewhere |
| R1 | **the claim**: linear-with-intercept factor per resolved step, kappa by step and buffer species (Davidek's [H2PO4-] reading); hyperbolic on lumped steps (CF-093) | **the plan's claim**: oil-water partition of MFT, FFT and their disulfides (fast equilibrium, reactions in water), phi from the printed oil volume, K_oil/water a declared band | B54: whole-phase concentration, residual table, Tg gate |
| R2 | a saturating factor (Michaelis) on every step; **split by K**: R2a K <= 0.1 M (Kallen's 17 - 29 mM), R2b K 0.3 M, R2c K >= 0.5 M (indistinguishable from linear below 0.2 M) | lipid-electrophile capture of the MFT pool, first order in MFT and lipid-carbonyl load, MFT-selective (CF-031's candidate) | concentration only, no residual table, no gate |
| R3 | buffer acts only through the in-situ pH it holds (the engine's charge balance already has this) | partition + capture | step-specific terms (the dehydration-limited condensation; a low-water exit), no shared factor |
| R4 | step-specific: Amadori decay and 2,3-enolisation only, not 1,2-enolisation, furfural or free-sugar steps | supply-side: lipid carbonyls take Cys and H2S (Zhao 2017b; thiophenes), thiols fall through their precursors (first order in H2S on both thiol routes, the engine's form) | -- |
| R5 | -- | measurement: Yiltirak's 30 uL aliquot from a creaming oil arm under-samples a lipophilic thiol (CF-031 (c), P 0.25). A measurement-model alternative scored beside every F rival on Yiltirak, not a mechanism | -- |
| R6 | -- | **added 2026-11-06 (CK-37, G9):** oil-phase oxidation: thiols partitioned into the oil meet its O2 and hydroperoxides; MFT is lost faster than FFT. First order in the oil's oxidant load and in MFT; zero under N2 with peroxide-free oil | -- |

## 4. The decisive data

### 4.1 Part B (buffer): DEVELOPMENT GRADE on disk

Every buffer datum on disk is spent (sec. 2.10). The units below test constant-free predictions only. "Falsified"
means the rival's prediction fails at every value of its free constants.

| unit | measurement and its noise | R0 / R3 | R2 (Michaelis, any k0, A >= 0) | R1, R4 | reading |
|---|---|---|---|---|---|
| FU-B1 (D1) Davidek 2002 | k_P / k_W at pH-stat 5 / 6 / 7 / 8: 2.79 - 3.11 / 3.35 - 3.62 / 1.94 - 2.17 / 1.24 - 1.27 (method-spread band; analytical CV < 10 %; one reaction run per cell) | predict 1 exactly: **falsified** at pH 5 - 7 (>= 0.29 dex against a <= 0.05 dex band; >= 3 sigma even at a 10 % run-to-run CV); pH 8 is ~1.5 sigma under that CV | predicts > 1 (size fitted) | predict > 1 (size fitted) | R0 and R3 falsified. Credits nothing |
| FU-B2 (D2) Bell 1997 | four rates at held pH (95 % CL) | R3: flat: **falsified** (x6 across the ladder) | [dry] sup k(0.5) / k(0.2) = 1.136 / 1.25 / 1.56 at K 0.05 / 0.1 / 0.3 M; measured 1.737, sigma_ln 0.161: **2.6 / 2.0 / 0.7 sigma**. Four-point non-negative fit, chi2 (2 dof): linear 1.37; K 0.05: 9.84; K 0.1: 4.27; K 0.3: 0.82; K 0.5: 0.16 | linear fits (chi2 1.37) | R2a with K <= 0.05 M falsified on the lumped condensation (p ~0.007); **K 0.1 M only marginally (p ~0.12 by chi2, 2.0 sigma by the ratio)**, so the plan row's "rejects K <= 0.1 M" is too strong on Bell alone. One vial per time point |
| FU-B3 (D2) Zyzak 1995 | fFL decay k at 10 / 50 / 100 / 200 mM, pH 7.4 held by the buffer at every concentration, sigma_k from the triplicates' SD (conservative) | R3: flat: **falsified** | [dry] chi2 (2 dof): linear 2.80; K 0.05: 71; K 0.1: 41; K 0.3: 13.1; K 0.5: 7.9 | linear fits | **R2a and R2b falsified on the Amadori decay step** (K <= 0.3 M, p <= 0.0015; K 0.5 at p ~0.02). 37 C, N2 + chelators: the transfer to process temperature is not shown |
| FU-B4 (D3) Mottram 2002 | MFT x2.2 - 2.6, FFT unchanged, ribose-only furfural 107 vs 103 (CV ~22 %) | the unbuffered arm drifts 5.6 -> 3.9, so R3's prediction needs the engine's own pH laws: **not constant-free** | one shared factor on every sugar step predicts FFT and furfural up | R4 predicts MFT up, FFT flat; R1's split depends on which steps carry kappa (CK-54) | development only. It falsifies at most "one shared factor on every sugar step", and only through the pH confound |
| FU-B5 Rizzi 2004 | ribose-alone A420 x87 from 0 to 0.10 M phosphate, pH 7.3 | -- | -- | R4 (no factor on amine-free steps) predicts ~no change | falsifies R4-as-written **only if** A420 tracks a sugar step ("browning is not the Amadori step", P 0.7): development only |

**What is left on disk:** R1, R2c (K >= 0.5 M) and R4, none credited. The spent TEST set agrees with the
falsifications (Potman's I(0.5) / I(0.1) = 5.7 against <= 1.67 for K <= 0.1 M; EX-058), and is itself falsify-only.
R1 vs R4 needs one pot read on two branches at held pH with a dose ladder. R1 vs R2c needs a point at >= 0.5 M, with
the concentrated-matrix ceiling (Potman Series A; candidate C-B2) kept apart from saturation.

### 4.2 Part F (fat): DEVELOPMENT GRADE on disk

| unit | measurement and its noise | R0 | R1 alone (partition) | R4 alone (supply) | R2, R3, R6 | reading |
|---|---|---|---|---|---|---|
| FU-F1 (D5) Elmore 2002, 18:2 | MFT 807 (508) -> 58 (71), FFT 92 (54) -> 87 (32), n = 4; ester phase phi ~0.008; MFT / FFT 8.8 -> 0.67 (x13.2, 1.12 dex; [dry] SE ~0.31 dex on the mean) | no change: **falsified** (MFT -1.14 dex, ~4 SE) | [dry] the ratio moves by (1 + phi K_FFT) / (1 + phi K_MFT): at most x3.3 (0.52 dex) for any log K in [1.0, 2.5], x1.44 for the case's [1.5, 2.0]; the absolute MFT fall at most x3.6. Measured x13.2: **falsified at 1.9 SE (wide band) to 3.1 SE (narrow)** | both thiols fall together: **falsified** (FFT x0.95 while MFT x1/14) | MFT-selective by construction, sizes fitted | R0, R1 alone and R4 alone falsified. Credits nothing |
| FU-F2 (D6) Elmore 2002, 18:3 vs 18:2 at equal ester amount | FFT 87 (32) -> 17 (5) (x5.1, 0.71 dex; SE ~0.10 dex); MFT 58 -> tr | -- | equal ester volume and near-equal ester polarity: no difference (constant-free): **falsified at ~7 SE (SE of the mean), ~3.5 SD** | predicts a fall only through more carbonyl (not constant-free) | scale with unsaturation / oxidant load | **the cleanest falsifier of partition alone, independent of the K band** |
| FU-F3 (D5) Yiltirak 2026, 100 C / 4 h, fresh canola 10 % | MFT 6.88 +- 0.98 -> 0.72 +- 0.10, FFT 1.28 +- 0.08 -> 1.19 +- 0.06 (n = 3): selectivity 0.95 dex (SE ~0.05 dex) | **falsified** | at phi 0.10 the bound is 0.42 dex for log K in [1.5, 2.0] (**falsified**, ~10 SE), but 0.89 dex when log K_MFT exceeds log K_FFT by 1.0, and 1.21 dex across [1.0, 2.5] (not falsified). **Undecided until K_MFT and K_FFT are measured** | **falsified** (selectivity 0.95 dex vs 0) | R2 needs a carbonyl load that fresh oil may not carry (CK-29's question); R6 does not; R5 (aliquot) explains the fall without either | R0, R4 alone falsified; R1 alone depends on the K split; R2 vs R5 vs R6 open |
| FU-F4 Yiltirak, emulsion vs buffer + oil (same oil volume) | MFT 3.73 vs 0.72 at 100 C | -- | equilibrium partition depends on volume only: equal (constant-free) **if both arms are sampled alike** | -- | the interface's role is not in any rival | **R1 vs R5**: decided only by the aliquot answer (CF-031 author question) |
| FU-F5 Whitfield 1988 lecithin | MFT x1/3.0, FFT x1/2.0, methylpyrazine x1/2.1 (n = 3) | falsified | not decisive (0.18 dex selectivity) | the sulfur-free methylpyrazine also falls: a supply or pH effect | -- | development only |
| D7 Zhang 2022 | phase-resolved, no n / SD | -- | the oil holds both species: consistent | -- | the water-phase pool falls | development only (no SD; the pre-registered rule failed) |

**What is left on disk:** R2, R3 and R6, with R5 as Yiltirak's measurement alternative. None is credited. Capture
and oil-phase oxidation differ on a **dosed carbonyl with no oil phase** (C-F1) and on an **N2 vs air arm with
peroxide-free oil** (designed). Partition + capture (R3) differs from capture alone only by partition's share,
<= 0.1 - 0.4 dex (cluster 4), so R3 needs phase-resolved free + disulfide MFT and a measured K. **No datum on disk,
spent or not, separates R2 from R3.**

### 4.3 Part W (water activity): DEVELOPMENT GRADE on disk

| unit | measurement and its noise | R0 (PR #77) | R2 (concentration only) | R1 (B54), R3 | reading |
|---|---|---|---|---|---|
| FU-W1 Bassama 2011 (acrylamide; MF-41) | k_F ratios 5.10 (0.904 / 0.972) and 1.58 (0.430 / 0.904); wet-basis mass action predicts 1.97 and 1.56; residual x2.60 (CI 1.9 - 3.6) and 1.01 (0.66 - 1.50) | the acrylamide lane answers (flat band, wet basis): the same as R2 here | **falsified at the wet arm** (constant-free ratio; the lower CI edge 1.9 > 1); supported in the two drier arms (spent, no credit) | B54 is not built on the acrylamide lane | R2 alone falsified at high moisture, with CF-104's caveats (a lumped k_F; a_w set at 25 C) |
| FU-W2 Leahy 1989a | pyrazine rate 1 : 2.10 : 4.23 : 4.07 at a_w 0.32 / 0.58 / 0.75 / 0.84 | the pyrazine lane has no a_w term, and 0.32 is below the trunk window: **refused** | the rate must fall as a_w rises (dilution) at any isotherm, unless the matrix is glassy: measured rise x4.2 to 0.75: **falsified if the matrix is rubbery** (its Tg is not on disk) | R1's gate refuses a glassy arm rather than predicting it | R2 alone falsified, conditional on the glass |
| FU-W3 (D8) Hofmann 1998 dry vs aqueous | ribose FFT x8.0 vs glucose FFT x0.50; triplicates within 10 % | sulfur lane: **refused** | same sign for both sugars only if temperature and time act alike on both: **not constant-free** (180 vs 145 C, 6 vs 20 min) | R3 can fit any sign | development only |
| FU-W4 (D9) Modi 2008 | MFT peak areas at a_w 0.60 / 0.80 / 0.90, n = 4, ratios only | sulfur lane: **refused** | the same sign as R1 without an isotherm | same sign | cannot decide |

**What is left on disk:** R1 and R3, plus R0 where it answers. R2 alone is falsified at high moisture and below the
a_w maximum. **R0 refuses on FU-W2 - W4.** A rival cannot beat a refusal on those units; it can only answer where R0
declines, which the acceptance rule (calibrated leave-one-paper-out adequacy) judges, not this comparison. W-R1 cannot
run at all without a protein isotherm and Tg curve (B54's need).

### 4.4 The decisive data that would separate the survivors (none on disk)

| part | survivors | the decisive measurement | record that writes it | proposed decisive-measurements entry |
|---|---|---|---|---|
| B | R1, R2c, R4 | xylose + Gly and glucose + Gly at held in-situ pH (pH-stat), phosphate 0 / 0.05 / 0.2 / 0.5 M and 0.2 M citrate at the same pH, at two pH values (5.5 and 7.0, for the species), reading Amadori decay, 1-DG / 3-DG trapped and, with Cys, MFT and FFT in the same pots | ENGINE_VS_LITERATURE sec. 3.1 E5 row; CF-011 buffer arm C with CF-094's citrate arm; CF-120's decider (ARP decay at 0 / 0.1 / 0.5 M, pH held) | `buffer_step_specificity` (PR body) |
| F | R2, R3, R6 (R5 a measurement model) | Cys + ribose pH 5.5, 120 C: buffer / + fresh oil / + oxidised oil (PV series) / + hexanal and (E,E)-2,4-decadienal dosed without oil; **air and N2 arms with peroxide-free oil**; free + disulfide (or TCEP-reduced) MFT and FFT in **both phases**; K_oil/water of MFT, FFT and the disulfide measured at 25 C and at process T | ENGINE_VS_LITERATURE sec. 3.1 E7 row; CF-031 remedy (designed); the N2 arm and the K measurement are added here (design not yet written) | `lipid_thiol_loss_route` (PR body) |
| W | R0, R1, R3 | acrylamide + Asn in one food matrix at >= 3 moistures spanning 10 - 70 % wb, a_w at temperature estimated (CF-104's decider); a pea / soy / gluten sorption isotherm and Tg-moisture curve (B54); a calibrated MFT / FFT series over >= 3 a_w in one defined matrix (design not yet written) | CF-104; B54 (backlog) | `water_activity_basis` (PR body) |

**The plan's decisive data, checked.**
- "Potman and Kallen (E5 TEST)": spent by EX-058 (inventory MF-26), so falsify-only. Their falsifications agree with
  FU-B2 / B3. They cannot settle the form, and no DEV datum credits it either (sec. 4.1).
- "The unseen fat contrasts (Kim 1996, Zamora 2015, Hidalgo 1993)": v2 TEST(E7), item 7 only. They are dienal + amine /
  lipid Strecker / epoxyalkenal browning chemistry (E20, F18), not thiol loss, so they could not decide F's rivals even
  if they were development data.
- The plan row's evidence column needs three corrections, proposed for the plan's writer:
  - Davidek's x2.8 - 3.6 holds at pH 5 - 6; at pH 7 it is x1.9 - 2.2.
  - Bell's 1.74 rejects K <= 0.1 M at only ~2 sigma. Zyzak's ladder is the stronger falsifier (K <= 0.3 M on the
    Amadori decay step).
  - "The data favour capture" rests on spent data. On disk they falsify partition alone and supply alone; capture,
    oil-phase oxidation and partition + capture stay undecided.

## 5. Tolerances, power and odds

### 5.1 Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| FU-B1 Davidek 2002 | "All samples analyzed in duplicate", no error bars; CV < 10 % above 1.5 mmol/L; one reaction run per cell; the ratio band is the min-max across the dossier's fitting methods (a) - (c) | davidek2002 dossier secs. 4.3 - 4.4 |
| FU-B2 Bell 1997 | Table 1 95 % CL per rate (sigma = CL / 1.96); one sample per time point, the CLs from regression scatter | bell1997 dossier |
| FU-B3 Zyzak 1995 | t1/2 +- SD, triplicates; sigma_k / k = SD / t1/2 (SD, not SE: with the SE even the linear law scores chi2 8.4 on 2 dof, so the SE understates the scatter) | zyzak1995 dossier |
| D2, D4 Guan 2012 | figure read errors (+-5 g/L glucose; DFP +-1e-4 mmol/L) | guan2012 dossier |
| FU-B4 Mottram 2002 | triplicate; mean CV ~22 %, no compound > 40 %; RF 1 (ratios only) | mottram2002 dossier |
| FU-F1, F2 Elmore 2002 | four replicates; mean (SD) per compound; SE of the mean used, SD-based sensitivity stated | elmore2002 dossier; CF-031 |
| FU-F3, F4 Yiltirak 2026 | triplicate cooks, mean +- SD (Table S3), SIDA with deuterated standards | yiltirak2026 dossier |
| FU-F5 Whitfield 1988 | n = 3, external standards | whitfield1988 dossier |
| D6 Farmer 1992 | n = 4, mean + SD | farmer1992 dossier |
| D7 Zhang 2022 | n and SD not printed; HS-SPME RF 1 | zhang2022 dossier |
| FU-W1 Bassama 2011 | 95 % CI of the fitted k_F, k_E (Table 2) | CF-104; bassama2011 dossier |
| FU-W2 Leahy 1989a | duplicates early, triplicates later | leahy1989a dossier |
| FU-W3 Hofmann 1998 | triplicates within 10 % | hofmann1998 ledger |
| FU-W4 Modi 2008 | peak areas, mean +- SD, n = 4; no calibration | modi2008 dossier |
| K_oil/water (F-R1) | **no measurement on disk**: a declared band, log K in [1.0, 2.5] for MFT, FFT and the disulfide (the case's log K_ow 1.5 - 2.0, widened by 0.5 for oil vs octanol). An unsourced width, so it needs the reviewer's agreement | -- |
| across laboratories | H0's between-laboratory term (EX-059) on absolute levels only; no M5 unit uses it | -- |

### 5.2 Bands the comparison needs (declared, unscored sources)

- **kappa_P** [5, 40] /M at pH 5 - 7.4 (DEV: Guan, Davidek, Zyzak, Bell), not T-flat (Guan's 0.4 - 0.75 dex fall
  110 -> 130 C); by species, Davidek's [H2PO4-] reading at an assumed pKa2 6.8 - 7.2.
- **carboxylates** on Amadori steps [0, ~4] /M; free-sugar steps rho x kappa_P, rho [0.1, 1] (CF-094).
- **F-R2's capture constant**: no measurement at process temperature on disk (Charles-Bernard's 25 C slowing is the
  only kinetic anchor). A free constant, so F-R2 is falsify-only everywhere it is fitted.
- **F-R6's oxidant load**: the printed PV of the oil where printed; otherwise not evaluable (refused).
- **W-R1's isotherm and Tg**: not on disk (B54); W-R1 is not evaluable on any food pot until they exist.

### 5.3 Power check: a reduced forced-input mini-model per part (SPEC; not built)

Synthetic data only. No measured value of any unit is read; it reuses M3 phase 1's forced-input layer and the
M4 mini-model's harness. Sizing ~1 day for all three parts (work), < 1 h of single-core compute. It is a power check,
not a comparison, so it does not wait on M1's verdict. It must simulate:

- **B.** One aqueous pot: sugar, amine, Amadori pool, 1-DG / 3-DG, a pentose-to-furfural branch and an MFT branch
  (with Cys). Rivals R1 (by total phosphate and by [H2PO4-]), R2a / b / c, R3 (via the charge balance only), R4.
  Noise per unit as sec. 5.1. Jobs:
  - (a) the falsifications of sec. 4.1 hold at every band corner, and the probability that a TRUE R1 is falsified by
    FU-B2's one-vial design (false falsification) is reported;
  - (b) the designed decider (sec. 4.4) at n = 3 and n = 6 per cell: P(R1 vs R4), P(R1 vs R2c), P(total vs species
    kappa) separated. R2c vs R1 needs a >= 0.5 M arm; the check says whether 1.0 M is needed;
  - (c) the concentrated-matrix ceiling as a separate factor, so a ceiling is not read as saturation.
- **F.** Two phases (water, oil at phi), fast partition per species (MFT, FFT, disulfide, H2S) with K in its band,
  reactions in water (R1), plus capture (R2), both (R3), supply drain (R4) and oil-phase oxidation (R6). The
  measurement model carries the phase tag (water / oil / whole) and the species tag (free / pool), and R5's aliquot
  model on the Yiltirak arms. The thiol sink is declared at M1's form (frozen, banded or C2). Jobs:
  - (a) FU-F1 / F2 / F3 falsifications at every K corner;
  - (b) the E7 designed row with its N2 arm at n = 3 / 6: P(R2 vs R6), P(R2 vs R3), P(R1 vs R5);
  - (c) whether C-F1's dosed-aldehyde design (if its thiols are quantified) separates R2 from R1 / R6 at its printed n.
- **W.** One phase with a moisture input and a declared isotherm band; rivals R1 (with gate), R2, R3 (a step term on
  the sulfur condensation). Jobs:
  - (a) FU-W1 / W2 falsifications under CF-104's caveats;
  - (b) the power of CF-104's >= 3-moisture acrylamide design and of a calibrated thiol a_w series, with the
    isotherm band at its declared width and at a measured width (+-10 % water).

### 5.4 P(the decisive data separate the rivals)

Calibrated by judgement; every number is to be revised by the power check.

| part | falsifications on disk hold (spent, falsify-only) | survivors separated by data on disk | with the sweep candidates read under fixed roles | by the designed decider (sec. 4.4) |
|---|---|---|---|---|
| B | R0 / R3: 0.9; R2a / R2b on the Amadori decay step: 0.8 | **0.05** (nothing unspent) | 0.15 (C-B1 / C-B2 / C-B3 bear on species, ceiling and step set, not on R1 vs R4 at held pH) | ~0.6 |
| F | R0: 0.95; R1 alone: 0.8 (FU-F2 is K-free); R4 alone: 0.7 | **0.05** | 0.25 (R2 vs R6 via C-F1, if quantified: P ~0.4 that it is) ; R2 vs R3: < 0.1 | ~0.55 (R2 vs R6 vs R3, with K measured) |
| W | R2 alone: 0.65 (FU-W1 at high moisture; FU-W2 conditional on the glass) | **0.05** | 0.15 (C-W1 / C-W2: trunk, 25 C) | ~0.5 (trunk and acrylamide); sulfur lane: no design written |

**P(M5 ships any term over R0 on data on disk) ~0.05.** In the late-November slot, M5 can build its network
pieces (sec. 6), record the falsifications above as development evidence, and size the W6 requests. It cannot credit
a buffer, fat or a_w term until a candidate or a designed experiment delivers unspent data.

## 6. What M5 needs from the one network (W3)

- **Buffer as a composition, not a label:** each buffer species with its concentration and pK (already in
  `ph_state` for the charge balance), and a per-step catalytic coefficient table keyed by (step, buffer species) as
  `ConfigKey`s with bands; a hyperbolic form available for lumped steps; the species fraction from the in-situ pH (so
  Davidek's [H2PO4-] reading is expressible). The acrylamide lane's phosphate-for-every-buffer law (MF-37 (b)) is
  replaced by the same table.
- **Constants re-referenced to buffer-free values** at compose time from each calibration pot's stated buffer, so the
  joint refit sees one convention (sulfur constants at 0.5 M; M3's `k_d`; M1c's `r_cama_form` if fitted on phosphate
  pots).
- **Two phases:** an oil phase with a volume fraction, fast partition equilibria per volatile (K_oil/water as declared
  constants), reactions confined to their phase, and observables that say which phase or the whole sample they read
  (the CF-031 tags `mft_species` and `phase`). The headspace operator must use the same partition once. **The oil
  volume of a food pot comes from a carrier mass fraction** (CF-130's schema field and engine change); without it, F's
  terms refuse on matrix pots.
- **Lane coupling:** the lipid lane's carbonyls (alkenals, dienals) and hydroperoxides as shared species that react
  with MFT, Cys and H2S (retiring the direct sum), with the lane-coupling guard replaced by declared competition. Their
  load in food pots stays refused until MF-42 is re-derived on the corrected carriers.
- **O2 in the oil:** an O2 solubility per phase, shared with the vessel operator of G9 (CF-068), so R6 is expressible.
- **Water:** moisture and a_w as conditions, a whole-phase concentration basis, a Tg gate (warn within 20 K, refuse
  below), and an a_w hook per step for the step-specific rival. The PR #77 window and its refusal stay the default
  until a rival is accepted.

## Open questions for the build

1. **Choose the buffer form on DEV data first** (D2): on disk, R2a / R2b are falsified and R1 vs R2c is open. Read
   Potman, Kallen and Thornalley in item 7 with the form fixed. The plan row's "settled" form is provisional.
2. **Which steps get kappa.** The 2,3-branch / Amadori steps exist in the sulfur lane only as `r_pent_dpo`,
   `r_arp_dpo` or, if M1's C1 ships, `r_cama_dpo`. Decide after M1's verdict (CK-54 runs inside the prereg).
3. **Kappa's temperature dependence.** Guan's fall with T is the only measurement; PH_BUFFERS_AND_WATER sec. 10
   carries kappa T-flat with a x/÷1.5 band. Declare which.
4. **Kappa's species.** Davidek's ratio tracks [H2PO4-] at an assumed pKa2; Rizzi's pH 7.3 and Mottram's pH 5.6 do not
   test it. Declare total vs species, and the pKa2 at temperature.
5. **K_oil/water for MFT, FFT and the disulfide** is an estimate. FU-F3's verdict on partition depends on the MFT /
   FFT difference; FU-F2's does not.
6. **Yiltirak's aliquot.** Ask the authors whether the 30 uL came from the water phase or a mixed sample (CF-031).
7. **Fat and M1.** Capture, oil-phase oxidation and M1's thiol loss law act on the same pool. Build F after M1's
   verdict, with `k_mft_decay`'s status (frozen, banded, C2) declared, or they compensate.
8. **B54's isotherm.** No pea / soy / gluten sorption isotherm or Tg curve is on disk; W-R1 cannot run without one.
9. **"Roos" in the plan row** means Lievonen & Roos 2002 and Miao & Roos 2004 (on disk); neither gives a step rate.

## What should block the freeze (for the review and the coordinator)

1. **The reduced power check (sec. 5.3)** must be built (~1 day) and run, and its demotions written in.
2. **The K_oil/water band** (sec. 5.1) is an unsourced width. It needs the reviewer's agreement, or a source.
3. **How a refusing R0 is scored** on W's units (FU-W2 - W4). Proposed: units R0 refuses are scored among the other
   rivals only, and a term that answers there ships through the acceptance rule's adequacy test, not by "beating" R0.
   This is a coordinator / owner call.
4. **PR #88 (CF-129, CF-130) is not merged.** This case cites them as pending. Their substance changes no decisive
   unit, but sec. 2.11's food-pot rule depends on CF-130's carrier fraction.
5. **An independent premortem / adversarial read of items 1 - 5** (the M4 precedent: one round plus one confirmation).
   Points to test:
   - whether FU-F2's "equal ester volume, near-equal polarity" holds;
   - whether Zyzak's 37 C N2 ladder may falsify R2b for process-temperature foods;
   - whether R6 is distinct enough from R2 to be scored (an oxidant vs an electrophile);
   - whether Leahy's matrix was rubbery, and whether Leahy 1989a Table IV, which B54 lists as a hold-out, may serve
     as a falsify-only unit at all (it is seen: quoted in this case and a calibration input).
6. **Candidates:** the owner decides on downloads (sec. 2.9). A candidate read before the freeze is folded in under its
   fixed role; one read after enters by amendment.
7. **Not blockers, but the owner should know:**
   - On disk, M5 credits no term in any part. It records falsifications only: R0 and R3 (buffer), R2 with K <= 0.3 M
     on the Amadori decay step, partition alone and supply alone (fat), and concentration-only at high moisture and
     below the a_w maximum.
   - The deciders are designed experiments, which come last by the owner's ruling (2026-11-05).
   - The network pieces (buffer table, two phases, O2 per phase, moisture basis) are still worth building for M6 and
     W5.

Seen by the author: as stated at the top. No engine output on any row named in sec. 4 was computed for this case.

## Revision log

- **2026-11-04 (v1).** Items 1 - 5 written before any result; no engine run.

### Changelog: 2026-11-06, READY TO FREEZE (items 1 - 5), not frozen

Still before any result of the change; no engine run on any pot.
- **The spent-data rule applied** (coordinator refinement: falsify constant-free predictions, never credit). Secs.
  2.10 and 4 rewritten around falsify-only units; every part marked development grade on disk.
- **M0 folded in.** G5's members (MF-05, -26, -28, -32, -37), MF-18, MF-41, MF-42, and the M5-owned checks (CK-54,
  CK-29, CK-33, and CK-37's Yiltirak part) placed (sec. 2.11). MF-28's buffer candidate retired (the pots are unbuffered,
  CK-06). Rival F-R6 (oil-phase oxidation) added from CK-37 / G9; R5 (Yiltirak's aliquot) made an explicit
  measurement alternative.
- **PR #88 folded in.** CF-130's undiluted carriers mis-state the oil volume (F-R1) and the carbonyl load (F-R2, F-R6)
  in every food pot. No F rival may take its load from the lipid lane on a matrix pot until the carrier fraction is
  built and MF-42 re-derived. CF-129's Bi row is not F evidence. No decisive unit moves.
- **PR #77 folded in.** R0 for W is the trunk table inside a_w 0.33 - 1.00, refusing below. Most W units are R0
  refusals, so a scoring rule is needed (blocker 3).
- **Corrections.**
  - Davidek's x2.8 - 3.6 is pH 5 - 6 (pH 7 x1.9 - 2.2; pH 8 x1.24 - 1.27).
  - The Michaelis bound at K 0.1 M is k(0.5) / k(0.2) <= 1.25, not ~1.1. On that bound Bell rejects K <= 0.1 M at
    ~2 sigma, not decisively.
  - [dry] Zyzak's four-point ladder falsifies K <= 0.3 M on the Amadori decay step.
  - Elmore's 18:3 vs 18:2 FFT contrast at equal ester amount falsifies partition alone without needing K.
  - Bell 1995's row restated: the 11-04 text's "glassy and rubbery give the same rates" contradicts the dossier,
    where the glassy system is the slowest where the polymers straddle the glass.
- **Existing evidence kept, with its status.** Davidek's pH-stat, Bell's ratio, Mottram 2002, the Elmore / Yiltirak
  selectivity and Hofmann dry vs wet are all kept, and all are spent: falsify-only, never credit.
- **Citation sweep run** (three sweeps, 13 seeds, 1,010 unseen works by title / abstract, 14 candidates with roles
  fixed before reading; sec. 2.9). It found no partition coefficient, no protein isotherm and no calibrated thiol
  a_w series.
- **Power check specified** (sec. 5.3), not built. **P(separation)** per part (sec. 5.4).
- **Decisive-measurement entries proposed** in the PR body (`buffer_step_specificity`, `lipid_thiol_loss_route`,
  `water_activity_basis`); the file itself is unchanged.
