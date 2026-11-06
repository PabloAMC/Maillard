# Case M5: buffer, fat and water activity

**READY TO FREEZE (items 1 - 5), not frozen** (2026-11-06). All three parts are **development grade on disk**: every
decisive datum on disk is spent or seen, so on disk the data can falsify rivals but credit none. Items 1 - 5 freeze
when the reduced power check of sec. 5.3 has run, its demotions are written in, and the review has answered "What
should block the freeze" (end of file). No engine was run for this case, on any pot. The changelog is at the end.

**Independent premortem (2026-11-06): "freeze after fixes".** All ten fixes are applied in this revision; still not
frozen. Every falsification is now a posterior-predictive p after Holm (sec. 4.0). Sec. 4.5 declares the S3
consequences for shipped outputs. **The late-November slot is scoped to the power check (sec. 5.3) and the S3
consequence computation (sec. 4.5).** The two-phase and O2-per-phase network builds wait until a decider exists.

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
  refusal (PR #77). The citation sweeps read titles and abstracts of unseen works only. For the premortem fixes:
  Zyzak's air ladder and HEPES labels, Elmore's Table 1 sulfur-free markers and Table 2 aldehydes, Leahy 1989a's a_w
  method, De Vleeschouwer 2008's k_Fref by column, and Guan 2012's phosphate factor at 110 / 130 C (all dossiers),
  and PR #96's S3 / S4 text. [dry] marks arithmetic on
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
| MF-26 (no buffer catalysis; Mottram, Davidek, Deng) | every rival but R0 (which now includes the pH-only reading) moves the axis; their sizes are fitted | M5's own misfit; spent; CK-54 (owned by M5) runs inside the prereg, falsify-only |
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
    lipid lane's output, which is mis-charged (CF-130) and, once corrected, too low (MF-42). The same charge rule
    over-states the model pots too: Yiltirak's 10 % canola would be charged ~10x and Elmore's ester phase (~0.74 % w/v)
    ~130x. So **no F rival may take its load from the lipid lane on ANY pot** until CF-130's carrier mass fraction is
    built and MF-42 re-derived. On every pot the load is the printed lipid charge (and PV where printed) or the
    measured aldehydes. This changes no decisive datum, because every F unit already uses printed charges.
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
| R0 | **current engine**: no buffer term; the buffer acts only through the in-situ pH its charge balance holds (ex-R3, merged 2026-11-06: the engine already does it) | direct sum, lipid inert | trunk table on `k_schiff` inside a_w 0.33 - 1.00 (PR #77); refusal elsewhere |
| R1 | **the claim**: linear-with-intercept factor per resolved step, kappa by step and buffer species (Davidek's [H2PO4-] reading); hyperbolic on lumped steps (CF-093) | **the plan's claim**: oil-water partition of MFT, FFT and their disulfides (fast equilibrium, reactions in water), phi from the printed oil volume, K_oil/water a declared band | B54: whole-phase concentration, residual table, Tg gate |
| R2 | a saturating factor (Michaelis) on every step; **split by K**: R2a K <= 0.1 M (Kallen's 17 - 29 mM), R2b K 0.3 M, R2c K >= 0.5 M (indistinguishable from linear below 0.2 M) | lipid-electrophile capture of the MFT pool, first order in MFT and lipid-carbonyl load, MFT-selective (CF-031's candidate) | concentration only, no residual table, no gate |
| R3 | merged into R0 (id kept so the references stay stable) | partition + capture | step-specific terms (the dehydration-limited condensation; a low-water exit), no shared factor |
| R4 | step-specific: Amadori decay and 2,3-enolisation only, not 1,2-enolisation, furfural or free-sugar steps | supply-side: lipid carbonyls take Cys and H2S (Zhao 2017b; thiophenes), thiols fall through their precursors (first order in H2S on both thiol routes, the engine's form). Tested **as written, at the B60 route structure**; re-test after M1c | **added 2026-11-06:** a solubility-capped concentration basis (the dissolved sugar / amino acid stops rising at its solubility; crystallisation, as of lactose, sets the maximum without a Tg gate) |
| R5 | -- | measurement: Yiltirak's 30 uL aliquot from a creaming oil arm under-samples a lipophilic thiol (CF-031 (c), P 0.25). A measurement-model alternative scored beside every F rival on Yiltirak, not a mechanism | -- |
| R6 | -- | **added 2026-11-06 (CK-37, G9):** oil-phase oxidation: thiols partitioned into the oil meet its O2 and hydroperoxides; MFT is lost faster than FFT. First order in the oil's oxidant load and in MFT; zero under N2 with peroxide-free oil. Products: disulfide, mixed disulfides, over-oxidised sulfur; no thioether (sec. 4.2) | -- |
| R7 | **added 2026-11-06:** saturating at K ~2 mM plus linear (a fast-saturating catalysed path on top of a linear one); explains a ladder intercept far above an inert-buffer rate (Zyzak's HEPES) | -- | -- |
| R8 | **added 2026-11-06:** ionic strength, not buffer species, sets the factor (I co-varies with buffer in Potman, Bell, Guan) | -- | -- |
| R9 | -- | **added 2026-11-06, measurement:** SPME fibre competition in Elmore's two-fibre room-temperature SPME (lipid volatiles crowd the fibre), handled by normalising to sulfur-free markers | -- |

## 4. The decisive data

### 4.0 How a falsification is stated (premortem fix 4)

Every falsification below is stated as a posterior-predictive p, with Holm applied across the datasets on which that
rival was tested (S3 on PR #96; M1.1). [dry] The p is the chi2 (or normal) tail of the rival's best non-negative fit
to the printed values. This is the posterior-predictive p under flat priors and Gaussian noise; item 6 recomputes it
under the full measurement model. A rival is falsified when its smallest Holm-adjusted p is < 0.05. Noise choices
that decide a verdict are disclosed in the row.

### 4.1 Part B (buffer): DEVELOPMENT GRADE on disk

Every buffer datum on disk is spent (sec. 2.10). The units below test constant-free predictions only.

| unit | measurement and its noise | R0 (no buffer term; incl. the pH-only reading, ex-R3) | R2 (Michaelis, any k0, A >= 0) | R1 (linear), R7 (saturating + linear) | reading |
|---|---|---|---|---|---|
| FU-B1 (D1) Davidek 2002 | k_P / k_W at pH-stat 5 / 6 / 7 / 8, low ends of the method band: 2.79 / 3.35 / 1.94 / 1.24. One reaction run per cell, so the run-to-run CV is unknown; **taken as >= 20 %** per k (sigma_ln of the ratio 0.28) | predicts 1: p = 3e-4 / 2e-5 / **0.019** / 0.45 at pH 5 / 6 / 7 / 8 | > 1, size fitted | > 1, size fitted | **R0 falsified at pH 5 - 7 by Davidek alone** (pH 7 is the weakest cell, p 0.019); **not at pH 8**. This is the only unit that tests the pH-only reading at a measured held pH |
| FU-B2 (D2) Bell 1997, glycine loss (lumped condensation), 25 C, pH 6.86 - 7.05 printed | four rates, 95 % CL (sigma = CL / 1.96); one vial per time point | flat: chi2 47.3 on 3 dof, p 3e-10 | sup k(0.5) / k(0.2) = 1.136 / 1.25 / 1.56 at K 0.05 / 0.1 / 0.3 M. Best fit p: K 0.05 **0.007**; K 0.1 0.12; K 0.3 0.66; K 0.5 0.93 | linear p 0.50; R7 (K 2 mM) p 0.28 | K = 0.05 M falsified on Bell; K = 0.1 M not |
| FU-B3 (D2) Zyzak 1995, fFL decay (the Amadori decay step), 37 C, pH 7.4, N2 + chelators | t1/2 +- SD of triplicates; **sigma_k / k = SD / t1/2 (SD, not SE)**. With the SE (÷√3) every two-parameter law fails, linear included (chi2 8.4 on 2 dof here and 21 on the air ladder), so the SE understates the scatter. The SD choice is what keeps R1 alive | flat: p ~1e-85 (not used for the pH-only reading: the ladder's pH is set, not shown held) | K 0.05: 3e-16; K 0.1: **1e-9**; K 0.3: **0.0014**; K 0.5: 0.019 | linear p 0.25; R7 (K 2 mM) p 0.094 | -- |
| FU-B3' Zyzak 1995 air ladder (no chelator), t1/2 10.7 / 9.8 / 5.8 / 4.2 d | as FU-B3 | flat: p ~1e-139 | chi2 (2 dof): K 0.1: **17.6** (p 1.5e-4); K 0.3: 7.27 (p 0.026); K 0.5: 5.60 (p 0.061) | linear chi2 **6.99** (p 0.030); R7 (K 2 mM) p 0.011 | the air ladder is tight at 10 and 100 mM (SD 2 %), so every smooth law is strained |
| HEPES sensitivity (R1 with the HEPES point as the buffer-free rate) | HEPES t1/2 ~200 d (N2), ~56 d (air), figure labels with no SD (sigma taken as 20 % of k) | -- | K 0.5 with HEPES: p 6e-14 (N2) | **linear with HEPES as k_u: p 3e-16 (N2), 5e-61 (air): falsified.** R7 (K 2 mM) with HEPES: p 0.086 (N2), 0.035 (air) | R1's intercept (0.0105 /d) sits 0.48 dex above HEPES (0.0035 /d). R1 survives only with an intercept band (sec. 4.5) or as R7. HEPES itself may not be inert |
| FU-B4 (D3) Mottram 2002 | MFT x2.2 - 2.6, FFT unchanged, ribose-only furfural 107 vs 103 (CV ~22 %) | the unbuffered arm drifts 5.6 -> 3.9: not constant-free | one shared factor on every sugar step predicts FFT and furfural up | R4 predicts MFT up, FFT flat; R1's split depends on the step set (CK-54) | development only |
| FU-B5 Rizzi 2004 | ribose-alone A420 x87 from 0 to 0.10 M phosphate, pH 7.3 | -- | -- | R4 as written predicts ~no change | development only (browning is not a step, P 0.7) |

**Holm across the datasets tested** (smallest adjusted p; falsified if < 0.05):

| rival | datasets | smallest Holm p | verdict |
|---|---|---|---|
| R0 (incl. the pH-only reading) | Davidek pH 7 (CV 20 %), Bell | 6e-10 (Bell); Davidek alone 0.019 | **falsified at pH 5 - 7**; not tested at pH >= 8 |
| R2, K <= 0.05 M | Bell, Zyzak N2, Zyzak air | 9e-16 | falsified |
| R2a, K = 0.1 M | same | 4e-9 (N2), 3e-4 (air); Bell 0.12 | **falsified on Amadori decay at 37 C (both Zyzak ladders)**; K = 0.05 M on Bell's condensation at 25 C |
| R2b, K = 0.3 M | same | 0.004 (N2); air 0.053; Bell 0.66 | **falsified only anoxically at 37 C** |
| R2c, K = 0.5 M | same | 0.057 | not falsified |
| R1 (linear, free intercept) | same | 0.091 (air ladder, raw p 0.030) | not falsified; falsified if its intercept must equal the HEPES rate |
| R7 (saturating at K ~2 mM + linear) | same; + HEPES | 0.031 at K = 2 mM without HEPES (air ladder); 0.07 with HEPES. With K free it has as many constants as a four-point ladder has points | carried (not testable at its best fit on disk) |
| R8 (ionic strength, not species) | none on disk: I co-varies with buffer in Potman, Bell and Guan; only Kallen held I | -- | carried, untested |

**Restated (premortem fix 2):** R2 with K <= 0.1 M is falsified on Amadori decay at 25 - 37 C (both Zyzak ladders;
K = 0.05 M on Bell). K = 0.3 M is falsified only anoxically at 37 C. K = 0.5 M is not falsified. The transfer of
any of these to 90 - 140 C is not shown.

**What is left on disk:** R1, R2c, R4, R7 and R8, none credited. R0 (and with it the pH-only reading) is falsified at
pH 5 - 7 only. The spent TEST set agrees (Potman's I(0.5) / I(0.1) = 5.7 against <= 1.67 for K <= 0.1 M; EX-058), and
is falsify-only.

### 4.2 Part F (fat): DEVELOPMENT GRADE on disk

Measurement models carried beside every rival: **R5** (Yiltirak's 30 uL aliquot from a creaming oil arm) and **R9**
(SPME fibre competition, for Elmore: a two-fibre room-temperature SPME of a pot whose lipid volatiles compete for
fibre sites lowers every analyte's signal, so Elmore's contrasts are normalised by sulfur-free Maillard markers).

| unit | measurement and its noise | R0 | R1 alone (partition) | R4 as written (supply, at the B60 route structure) | R2, R3, R6 | reading |
|---|---|---|---|---|---|---|
| FU-F1 (D5) Elmore 2002, 18:2 | MFT 807 (508) -> 58 (71), FFT 92 (54) -> 87 (32), n = 4; MFT / FFT 8.8 -> 0.67 (1.12 dex) | no change: p 4e-5 | **not used** (the earlier "falsified at any log K in [1.0, 2.5]" was wrong: at the wide band partition reaches 0.52 dex, and the rest is within noise and fibre competition). FU-F2 carries partition alone | both thiols first order in the same H2S / Cys: the ratio stays: p 4e-4 | MFT-selective by construction | R0 and R4 as written falsified |
| FU-F2 (D6) Elmore 2002, 18:3 vs 18:2 at equal ester amount | FFT 87 (32) -> 17 (5): 0.71 dex, SE 0.10. **Fibre-competition normalisation (R9):** the sulfur-free markers fall too (2-methylfuran 170 -> 128, 1-hydroxy-2-propanone 43 -> 28, 2-butanone 24 -> 10). Normalised FFT contrast: 0.59 dex (SE 0.14, by 2-methylfuran), 0.52 (0.15, by hydroxypropanone), 0.33 (0.15, by 2-butanone); by the mean of the three 0.48 dex (SE 0.12) | -- | equal ester volume and near-equal ester polarity: no difference. p 2e-5 / 4e-4 / 0.026 by marker; **4.5e-5 by the mean**. Holm over the two units that test it (FU-F2, FU-F3): 9e-5 | not constant-free | scale with unsaturation or oxidant load | **partition alone falsified after normalisation and Holm.** On the 2-butanone normaliser alone it is marginal (Holm 0.051) |
| FU-F3 (D5) Yiltirak 2026, 100 C / 4 h, fresh canola 10 % | MFT 6.88 +- 0.98 -> 0.72 +- 0.10, FFT 1.28 +- 0.08 -> 1.19 +- 0.06 (n = 3): selectivity 0.95 dex (SE ~0.05) | p ~1e-86 | the bound is 0.42 dex for log K in [1.5, 2.0], 0.89 dex when log K_MFT exceeds log K_FFT by 1.0, 1.21 dex across [1.0, 2.5]: **not falsified at the declared band** | p ~1e-69 | R2 needs a carbonyl load fresh oil may not carry (CK-29); R6 does not; R5 explains the fall without either | R0 and R4 as written falsified; partition alone undecided here |
| FU-F4 Yiltirak, emulsion vs buffer + oil (same oil volume) | MFT 3.73 vs 0.72 at 100 C | -- | equal at equilibrium **if both arms are sampled alike** | -- | the interface's role is in no rival | R1 vs R5 only (CF-031 author question) |
| FU-F6 Elmore 2002 Table 2, the reverse flow | 18:3 pot vs 18:3 + Cys + ribose: (E)-2-butenal 497 -> nd (> 250x), 2-propenal 132 -> nd, (E)-2-pentenal 437 -> tr (> 44x), (E,E)-2,4-heptadienal 3083 (972) -> 297 (126) (x10.4, p 4e-19), (E,Z)-2,4-heptadienal 646 -> 43 (x15) | **the direct sum predicts the lipid aldehydes unchanged by Cys + ribose: falsified in reverse** | -- | consumes the enals as H2S / Cys sinks: consistent | R2 consumes them by capture; R6 does not require it | R0's direct sum is falsified in both directions |
| FU-F5 Whitfield 1988 lecithin | MFT x1/3.0, FFT x1/2.0, methylpyrazine x1/2.1 (n = 3) | falsified | not decisive | the sulfur-free methylpyrazine falls too | -- | development only |
| D7 Zhang 2022 | phase-resolved, no n / SD | -- | the oil holds both species | -- | the water-phase pool falls | development only |

**Holm:** R0 smallest p 9e-86 (Yiltirak; Elmore MFT 4e-5; the reverse flow 8e-19): falsified. Partition alone 9e-5
(FU-F2, normalised): falsified. R4 as written 2e-69 (Yiltirak; Elmore 4e-4): falsified **as written, at the B60 route
structure** (both thiol routes first order in the shared H2S). It is to be re-tested after M1c, whose supply structure
can change which pools the lipid drains.

**R6's product predictions (premortem fix 5).** Oil-phase oxidation loses MFT to bis(2-methyl-3-furyl) disulfide, to
mixed disulfides (with Cys, with FFT) and to over-oxidised sulfur (sulfinic / sulfonic forms, not GC-volatile). The
oil's hydroperoxides fall as the matching hydroxy esters rise. It forms no thiol-enal thioether. It needs O2 or a
printed PV: zero loss under N2 with peroxide-free oil. Capture (R2) instead forms 3-(2-methyl-3-furylthio)alkanals and
other thiol-enal Michael adducts, and consumes the enals. It needs an enal, whether dosed or formed in situ (which needs
O2 when the oil is fresh). Elmore's pool fall (free + disulfide x1/17) is consistent with both, since over-oxidation
and mixed disulfides also remove the pool.

**What is left on disk:** R2, R3 and R6, with R5 and R9 as measurement models. None is credited.

### 4.3 Part W (water activity): DEVELOPMENT GRADE on disk

| unit | measurement and its noise | R0 (PR #77; the acrylamide lane's flat band) | R2 (concentration only) | R1 (B54), R3, R4 (solubility-capped basis) | reading |
|---|---|---|---|---|---|
| FU-W1 Bassama 2011 (acrylamide; MF-41) | k_F ratios 5.10 (0.904 / 0.972) and 1.58 (0.430 / 0.904); 95 % CI of the fitted k | the flat multiplier with wet-basis mass action, its band drawn with **common random numbers across the arms** (so it cancels in the ratio): predicts 1.97 at the wet arm; **residual x2.60 (CI 1.9 - 3.6), p 5e-9: falsified above a_w ~0.9**. De Vleeschouwer 2008's own k_Fref at its wettest column, x0.41 of the 0.92 column, agrees in direction | **basis-independent:** wet basis misses the wet arm x2.60; dry basis misses it more (1 vs 5.10); per-water basis gives x1.52 at the wet arm but misses the dry arm ~x5 (1.58 vs 7.99) | B54 is not built on this lane; R4 (a capped dissolved concentration) predicts the dry arm's flattening, untested | **R0's acrylamide lane and R2 falsified above a_w ~0.9**, with CF-104's caveats (a lumped k_F; a_w set at 25 C) |
| FU-W2 Leahy 1989a (**demoted to development evidence**) | pyrazine rate 1 : 2.10 : 4.23 : 4.07 at a_w 0.32 - 0.84 in nonfat dry milk at 95 C; amorphous lactose can crystallise over this range, and the "ppm/hr" are per initial dry weight from a purged powder (dossier Flag 2), not solution concentrations | refused (no pyrazine a_w term; 0.32 below the trunk window) | the sign test needs a dissolved-concentration basis the units do not give, and crystallisation changes the dissolved sugar: **not constant-free** | R4 fits a maximum by construction | development only |
| FU-W3 (D8) Hofmann 1998 dry vs aqueous | ribose FFT x8.0 vs glucose FFT x0.50 | refused | not constant-free (180 vs 145 C, 6 vs 20 min) | R3 can fit any sign | development only |
| FU-W4 (D9) Modi 2008 | MFT peak areas at a_w 0.60 / 0.80 / 0.90 | refused | same sign as R1 without an isotherm | same sign | cannot decide |

**Stale hold-out tags (for the coordinator's next data PR; no data edited here).** Leahy 1989a Table IV is named as a
B54 hold-out (backlog) but is spent: it is a calibration input, its numbers are quoted in this case, and it is atlas
F19's a_w mismatch. The Yiltirak 2026 buffer arms are a Maillard-path hold-out bundle, yet their values are printed in
CF-031 (Table S3 rows). Both tags should read "seen / spent".

**What is left on disk:** R1, R3 and R4, plus R0 where it answers (the trunk inside 0.33 - 1.00; acrylamide below
a_w ~0.9). R0 refuses on FU-W2 - W4. A rival cannot beat a refusal there; answering where R0 declines is judged by the
acceptance rule's adequacy test. W-R1 cannot run without a protein isotherm and Tg curve (B54's need).

### 4.4 The decisive data that would separate the survivors (none on disk)

| part | survivors | the decisive measurement | record that writes it | proposed decisive-measurements entry |
|---|---|---|---|---|
| B | R1, R2c, R4, R7, R8 | xylose + Gly and glucose + Gly at held in-situ pH (pH-stat), phosphate 0 / 0.005 / 0.05 / 0.2 / 0.5 M, 0.2 M citrate and an NaCl arm at matched ionic strength (R8), an inert-buffer zero (to anchor the intercept against HEPES), at two pH values (5.5 and 7.0, for the species), reading Amadori decay, 1-DG / 3-DG trapped and, with Cys, MFT and FFT in the same pots | ENGINE_VS_LITERATURE sec. 3.1 E5 row; CF-011 buffer arm C with CF-094's citrate arm; CF-120's decider; the 5 mM, NaCl and inert-zero arms are added here (design not yet written) | `buffer_step_specificity` (PR body) |
| F | R2, R3, R6 (R5, R9 measurement models) | Cys + ribose pH 5.5, 120 C, free + disulfide (or TCEP-reduced) MFT and FFT read in **both phases**, with: (a) buffer vs fresh oil under **air and N2**; (b) **an enal dosed without oil** (C-F1's design); (c) **PV-spiked oil under N2 vs the enal dosed into oil under N2**; (d) the products read in both phases: thiol-enal thioether adducts (R2) vs disulfide and oxidised MFT (R6); (e) K_oil/water of MFT, FFT and the disulfide | ENGINE_VS_LITERATURE sec. 3.1 E7 row; CF-031 remedy; arms (a) - (e) added here (design not yet written) | `lipid_thiol_loss_route` (PR body) |
| W | R0, R1, R3, R4 | acrylamide + Asn in one food matrix at >= 3 moistures spanning 10 - 70 % wb, a_w at temperature estimated (CF-104's decider); a pea / soy / gluten sorption isotherm, Tg-moisture curve and sugar solubility (B54; R4); a calibrated MFT / FFT series over >= 3 a_w in one defined matrix (design not yet written) | CF-104; B54 (backlog) | `water_activity_basis` (PR body) |

**What the F arms decide (corrected mapping).**
- Arm (a), the N2 arm with fresh oil, separates **{R2, R6} from {R1, R5}**. Both chemical routes need O2 (R6 directly;
  R2 for enals formed in situ from fresh oil), while partition and the aliquot artefact do not. It does **not**
  separate R2 from R6.
- R2 vs R6 is separated by:
  - (b) an enal without oil: R2 lowers MFT, R6 does not;
  - (c) PV-spiked oil under N2 vs an enal in oil under N2: R6 responds to the first, R2 to the second, with the aldehydes
    read to catch peroxides decomposing into enals;
  - (d) the products: thioethers mean R2, disulfide / oxidised MFT means R6.
- R3 vs R2 needs (e) and both phases.

**The plan's decisive data, checked.**
- "Potman and Kallen (E5 TEST)": spent by EX-058 (inventory MF-26), so falsify-only. Their falsifications agree with
  FU-B2 / B3. They cannot settle the form, and no DEV datum credits it either (sec. 4.1).
- "The unseen fat contrasts (Kim 1996, Zamora 2015, Hidalgo 1993)": v2 TEST(E7), item 7 only. They are dienal + amine /
  lipid Strecker / epoxyalkenal browning chemistry (E20, F18), not thiol loss, so they could not decide F's rivals even
  if they were development data.
- The plan row's evidence column needs these corrections, proposed for the plan's writer:
  - Davidek's x2.8 - 3.6 holds at pH 5 - 6; at pH 7 it is x1.9 - 2.2, and at pH 8 the buffer effect is not resolved.
  - Bell's 1.74 does not reject K = 0.1 M (p 0.12). K <= 0.1 M is falsified on Amadori decay by both Zyzak ladders;
    K = 0.3 M only anoxically.
  - "The data favour capture" rests on spent data. On disk they falsify the direct sum (both directions), partition
    alone and supply as written; capture, oil-phase oxidation and partition + capture stay undecided.

### 4.5 S3 consequences: what the falsifications do to shipped outputs (premortem fix 1)

Under S3 (PR #96, not merged), a falsified incumbent with no rival credited means refusing the outputs that depend on
it: those whose sensitivity to the falsified structure's switch exceeds the lane's within-laboratory term. Where every
carried rival agrees in sign, the direction is shown and the size comes from the union of the rivals' envelopes, or is
refused under S4 when a rival has no band. This section declares the consequences. It does not apply them: the
application is its own S3 case, with its review.

| part | output family | scope (pH, T, composition) | consequence | basis |
|---|---|---|---|---|
| F | MFT, FFT (and the disulfide) in fat-containing recipes with PUFA oil or ester | Cys + pentose, pH ~5.5, 100 - 140 C (Elmore, Yiltirak, Whitfield) | **direction shown** (lipid lowers MFT, MFT more than FFT; every carried rival agrees in sign); **size refused** under S4: capture's constant has no band at process T, and R6's load needs a printed PV | R0 falsified (FU-F1, F3); R2, R3, R6 agree in sign |
| F | the same with free fatty acids or a structured-lipid mesophase | GSH + glucose 145 C (Du 2025: MFT +19 - 30 %); the Nestle mesophase (US 8,920,862, x10 up, a viewed TEST carried failure) | **direction refused** (signs disagree across lipid forms) | CF-031 |
| F | lipid-aldehyde outputs (alkenals, dienals) in recipes with Cys + sugar | pH ~5.5, 140 C (Elmore Table 2) | **direction shown** (lowered, 10 - > 250x in Elmore); size refused (no band) | the reverse flow, FU-F6 |
| F | thiol and lipid-aldehyde outputs in fat-free recipes | all | unchanged | -- |
| B | outputs whose constants came from 0.5 M phosphate pots (the sulfur lane's Hofmann-fitted constants; M3's `k_d`), served in unbuffered or carboxylate-buffered foods | **pH 5 - 7 only**, 90 - 140 C | **widened**: the union of R1 / R2c / R4 / R7 / R8 envelopes; refused where the step set matters and rivals give zero vs a fall at a sensitivity above the within-lab term (FFT and furfural under R1 vs R4) | R0 falsified at pH 5 - 7 (FU-B1) |
| B | the same at **pH >= 8** (the alkaline pyrazine outputs, Leahy's pH 9 borate pots) | pH >= 8 | **unchanged**: R0 is not falsified at pH 8 (Davidek p 0.45 at CV 20 %) | FU-B1 |
| B | outputs from 0.1 M phosphate constants (Martins 2005, Knol 2005) in unbuffered foods | pH 5 - 7 | widened by the same union, scaled to 0.1 M | FU-B1 |
| W | acrylamide formation contrasts crossing a_w ~0.9, and acrylamide levels above a_w ~0.9 (high-moisture plant matrices, 55 - 70 % water) | 140 - 200 C | **contrasts refused; levels widened** (the formation band's low edge to x0.28, CF-104's proposed edge) | R0's acrylamide lane falsified on FU-W1 with common random numbers |
| W | trunk outputs inside a_w 0.33 - 1.00; acrylamide below a_w ~0.9 | -- | unchanged (no constant-free test fails) | -- |
| W | sulfur, lipid and pyrazine a_w comparisons | -- | unchanged: already refused (PR #77) | -- |

**Declared bands for the buffer widening.**
- **kappa(T), bounded.** The phosphate rate factor (1 + kappa c, at Guan's 0.1 M) falls by a factor in
  **[1.5, 2.6] from 110 to 130 C**. Guan 2012's two readings bracket it: x~4 -> x~2.5 in the dossier (1.6) and
  x3.8 -> x1.5 in sec. 2.2's row (2.5). On kappa itself the second reading is a fall of x5.6, so the band is declared
  on the factor, not on kappa. The factor is held at its 90 - 110 C value below 110 C (Davidek, Guan) and carried with
  the full fall above 130 C. Nothing beyond 140 C is evaluated.
- **The zero-buffer intercept, 0.5 dex.** The phosphate ladder's extrapolated intercept (Zyzak N2, 0.0105 /d) sits
  0.48 dex above the HEPES rate (0.0035 /d). The buffer-free k_u of a re-referenced constant therefore carries a 0.5 dex
  band (one-sided, downward). This band is what lets R1 survive the HEPES sensitivity.
- **pH scope.** R0 is falsified at pH 5 - 7, not at pH >= 8. No buffer consequence is applied at pH >= 8, so the
  alkaline pyrazine outputs are unaffected.

**The sensitivity computation and its cost.** S3 requires carried rivals to be runnable. Until W3 runs them in JAX,
they run in the NumPy engine behind switches.
- **Fat:** no run is needed. The size is refused and the direction follows from the rivals' sign agreement, so no
  two-phase or oxidation switch is built for S3.
- **Water:** no run is needed. The refusal and widening are declared over an a_w region.
- **Buffer:** one switch is needed in the NumPy engine: the per-step factor table on the steps each rival names. Work
  about 0.5 day. Each B output in scope is then run under R0 and the five carried buffer rivals (six settings) over S4's
  joint draws, with common random numbers across settings. Its sensitivity to the switch is compared with the lane's
  within-lab term.
- **Cost:** the per-pot cost of a draw is not on record. The first step is therefore one timed setting on the
  affected rows (one core, niced). The full cost (about 6x that) is stated before the full run, capped at 8 of 10
  cores.

**Scope of the late-November slot (coordinator, 2026-11-06).** The slot holds two things: the power check (sec. 5.3),
and this S3 consequence computation (the buffer switch plus the declarations above). The two-phase and O2-per-phase
network builds (sec. 6) are **deferred until a decider exists** (a candidate read under its role, or a W6 result).

## 5. Tolerances, power and odds

### 5.1 Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| FU-B1 Davidek 2002 | "All samples analyzed in duplicate", no error bars; CV < 10 % above 1.5 mmol/L (analytical); one reaction run per cell, so the run-to-run CV is unknown and is **taken as >= 20 %** per k for the falsification; the method band (min-max of the dossier's fitting methods (a) - (c)) is reported beside it | davidek2002 dossier secs. 4.3 - 4.4 |
| FU-B2 Bell 1997 | Table 1 95 % CL per rate (sigma = CL / 1.96); one sample per time point, the CLs from regression scatter | bell1997 dossier |
| FU-B3, B3' Zyzak 1995 | t1/2 +- SD, triplicates (N2 and air ladders); sigma_k / k = SD / t1/2. **The SD, not the SE, is used**: with the SE every two-parameter law fails, linear included (chi2 8.4 N2, 21 air, 2 dof). The choice keeps R1 alive and is disclosed as such. HEPES half-lives are figure labels with no SD (20 % assumed in the sensitivity) | zyzak1995 dossier |
| D2, D4 Guan 2012 | figure read errors (+-5 g/L glucose; DFP +-1e-4 mmol/L) | guan2012 dossier |
| FU-B4 Mottram 2002 | triplicate; mean CV ~22 %, no compound > 40 %; RF 1 (ratios only) | mottram2002 dossier |
| FU-F1, F2, F6 Elmore 2002 | four replicates; mean (SD) per compound; SE of the mean used; two-fibre SPME at room temperature, so contrasts are normalised by sulfur-free markers (R9) and the normaliser is reported | elmore2002 dossier; CF-031 |
| FU-F3, F4 Yiltirak 2026 | triplicate cooks, mean +- SD (Table S3), SIDA with deuterated standards | yiltirak2026 dossier |
| FU-F5 Whitfield 1988 | n = 3, external standards | whitfield1988 dossier |
| D6 Farmer 1992 | n = 4, mean + SD | farmer1992 dossier |
| D7 Zhang 2022 | n and SD not printed; HS-SPME RF 1 | zhang2022 dossier |
| FU-W1 Bassama 2011 | 95 % CI of the fitted k_F, k_E (Table 2) | CF-104; bassama2011 dossier |
| FU-W2 Leahy 1989a (development) | duplicates early, triplicates later; per initial dry weight from a purged powder | leahy1989a dossier |
| FU-W3 Hofmann 1998 | triplicates within 10 % | hofmann1998 ledger |
| FU-W4 Modi 2008 | peak areas, mean +- SD, n = 4; no calibration | modi2008 dossier |
| K_oil/water (F-R1) | **no measurement on disk**: a declared band, log K in [1.0, 2.5] for MFT, FFT and the disulfide (the case's log K_ow 1.5 - 2.0, widened by 0.5 for oil vs octanol). An unsourced width, so it needs the reviewer's agreement. Since fix 4 no falsification depends on it | -- |
| across laboratories | H0's between-laboratory term (EX-059) on absolute levels only; no M5 unit uses it | -- |

### 5.2 Bands the comparison needs (declared, unscored sources)

- **kappa_P** [5, 40] /M at pH 5 - 7.4 (DEV: Guan, Davidek, Zyzak, Bell); its temperature fall bounded as in sec. 4.5
  (the factor falls x1.5 - x2.6 from 110 to 130 C); by species, Davidek's [H2PO4-] reading at an assumed pKa2 6.8 -
  7.2.
- **The zero-buffer intercept**: a 0.5 dex one-sided band (sec. 4.5).
- **carboxylates** on Amadori steps [0, ~4] /M; free-sugar steps rho x kappa_P, rho [0.1, 1] (CF-094).
- **R7's K** in [1, 5] mM (a declared width around the ~2 mM that the HEPES point implies; needs agreement).
- **F-R2's capture constant**: no measurement at process temperature on disk (Charles-Bernard's 25 C slowing is the
  only kinetic anchor). No band, so under S4 any size it sets is refused.
- **F-R6's oxidant load**: the printed PV of the oil where printed; otherwise refused.
- **W-R1's isotherm and Tg, W-R4's solubility**: not on disk (B54); not evaluable on any food pot until they exist.

### 5.3 Power check: a reduced forced-input mini-model per part (SPEC; not built)

Synthetic data only. No measured value of any unit is read; it reuses M3 phase 1's forced-input layer and the
M4 mini-model's harness. Sizing ~1 day for all three parts (work), < 1 h of single-core compute. It is a power check,
not a comparison, so it does not wait on M1's verdict. It must simulate:

- **B.** One aqueous pot: sugar, amine, Amadori pool, 1-DG / 3-DG, a pentose-to-furfural branch and an MFT branch
  (with Cys). Rivals R0 (the charge balance only), R1 (by total phosphate and by [H2PO4-]), R2a / b / c, R4, R7 and R8
  (ionic strength from the pot's salts). Noise per unit as sec. 5.1. Jobs:
  - (a) the falsifications of sec. 4.1 hold at every band corner, with the probability that a TRUE R1 is falsified
    by FU-B2's one-vial design or by Zyzak's SD choice (false falsification) reported;
  - (b) the designed decider (sec. 4.4) at n = 3 and n = 6 per cell: P(R1 vs R4), P(R1 vs R2c), P(R1 vs R7),
    P(R8 vs species), P(total vs species kappa) separated. R2c vs R1 needs a >= 0.5 M arm, and the check says whether
    1.0 M is needed. R7 needs the 5 mM and inert-zero arms;
  - (c) the concentrated-matrix ceiling as a separate factor, so a ceiling is not read as saturation.
- **F.** Two phases (water, oil at phi), fast partition per species (MFT, FFT, disulfide, H2S) with K in its band,
  reactions in water (R1), plus capture (R2), both (R3), supply drain (R4) and oil-phase oxidation (R6). The
  measurement model carries the phase tag (water / oil / whole), the species tag (free / pool), R5's aliquot model on
  the Yiltirak arms and R9's fibre competition on the Elmore arms. The thiol sink is declared at M1's form (frozen,
  banded or C2). Jobs:
  - (a) the FU-F1 / F2 / F3 / F6 falsifications at every K corner and normaliser;
  - (b) the E7 arms (a) - (e) at n = 3 / 6: P({R2, R6} vs {R1, R5}) on the N2 arm; P(R2 vs R6) on arms (b) - (d);
    P(R2 vs R3) with K measured;
  - (c) whether C-F1's dosed-aldehyde design (if its thiols are quantified) separates R2 from R1 / R6 at its printed n.
- **W.** One phase with a moisture input and a declared isotherm band; rivals R0 (acrylamide flat band, CRN), R1 (with
  gate), R2, R3 (a step term on the sulfur condensation), R4 (solubility cap). Jobs:
  - (a) FU-W1 under CF-104's caveats, on every concentration basis;
  - (b) the power of CF-104's >= 3-moisture acrylamide design and of a calibrated thiol a_w series, with the
    isotherm band at its declared width and at a measured width (+-10 % water), and with R4's solubility cap.

### 5.4 P(the decisive data separate the rivals)

Calibrated by judgement; every number is to be revised by the power check.

| part | falsifications on disk hold (Holm p < 0.05; spent, falsify-only) | survivors separated by data on disk | with the sweep candidates read under fixed roles | by the designed decider (sec. 4.4) |
|---|---|---|---|---|
| B | R0 at pH 5 - 7: 0.85 (Davidek's run-to-run CV is assumed); R2 with K <= 0.1 M on Amadori decay at 37 C: 0.85; K = 0.3 M (anoxic only): 0.6 | **0.05** (nothing unspent) | 0.15 (C-B1 / C-B2 / C-B3 bear on species, ceiling and step set, not on R1 vs R4 at held pH) | ~0.5 (more rivals now: R7, R8) |
| F | R0 (both directions): 0.95; partition alone (FU-F2, normalised): 0.75; R4 as written: 0.7 | **0.05** | 0.25 (R2 vs R6 via C-F1, if quantified: P ~0.4 that it is); R2 vs R3: < 0.1 | ~0.55 (arms (a) - (e), with K measured) |
| W | R0's acrylamide lane and R2 above a_w ~0.9: 0.75 (FU-W1; basis-independent) | **0.05** | 0.15 (C-W1 / C-W2: trunk, 25 C) | ~0.5 (trunk and acrylamide); sulfur lane: no design written |

**P(M5 ships any term over R0 on data on disk) ~0.05.** The late-November slot holds the power check and the S3
consequence computation (sec. 4.5), nothing else from M5. No buffer, fat or a_w term is credited until a candidate
or a designed experiment delivers unspent data.

## 6. What M5 needs from the one network (W3)

**Deferred (coordinator, 2026-11-06):** the two-phase and O2-per-phase builds below wait until a decider exists. Only
the buffer switch of sec. 4.5 is built now, in the NumPy engine.

- **Buffer as a composition, not a label:** each buffer species with its concentration and pK (already in
  `ph_state` for the charge balance), and a per-step catalytic coefficient table keyed by (step, buffer species) as
  `ConfigKey`s with bands; a hyperbolic form available for lumped steps; the species fraction from the in-situ pH (so
  Davidek's [H2PO4-] reading is expressible); the ionic strength of the pot (R8). The acrylamide lane's
  phosphate-for-every-buffer law (MF-37 (b)) is replaced by the same table.
- **Constants re-referenced to buffer-free values** at compose time from each calibration pot's stated buffer, with
  the 0.5 dex intercept band, so the joint refit sees one convention (sulfur constants at 0.5 M; M3's `k_d`; M1c's
  `r_cama_form` if fitted on phosphate pots).
- **Two phases (deferred):** an oil phase with a volume fraction, fast partition equilibria per volatile (K_oil/water
  as declared constants), reactions confined to their phase, and observables that say which phase or the whole sample
  they read (the CF-031 tags `mft_species` and `phase`). The headspace operator must use the same partition once.
  **The oil volume and the lipid load of any pot come from printed charges or a carrier mass fraction** (CF-130's
  schema field and engine change); without them, F's terms refuse.
- **Lane coupling (deferred):** the lipid lane's carbonyls (alkenals, dienals) and hydroperoxides as shared species
  that react with MFT, Cys and H2S (retiring the direct sum, which FU-F6 falsifies in reverse too), with the
  lane-coupling guard replaced by declared competition.
- **O2 in the oil (deferred):** an O2 solubility per phase, shared with the vessel operator of G9 (CF-068), so R6 is
  expressible.
- **Water:** moisture and a_w as conditions, a whole-phase concentration basis with a solubility cap (R4), a Tg gate
  (warn within 20 K, refuse below), and an a_w hook per step for the step-specific rival. The PR #77 window and its
  refusal stay the default, with sec. 4.5's acrylamide refusal above a_w ~0.9 added once S3 applies.

## Open questions for the build

1. **Choose the buffer form on DEV data first** (D2). On disk, R2 with K <= 0.1 M is falsified on Amadori decay at
   37 C. K = 0.3 M is falsified only anoxically. R1, R2c, R7 and R8 are open. Read Potman, Kallen and Thornalley in
   item 7 with the form fixed. The plan row's "settled" form is provisional.
2. **Which steps get kappa.** The 2,3-branch / Amadori steps exist in the sulfur lane only as `r_pent_dpo`,
   `r_arp_dpo` or, if M1's C1 ships, `r_cama_dpo`. Decide after M1's verdict (CK-54 runs inside the prereg).
3. **Kappa's temperature dependence:** bounded in sec. 4.5. Whether PH_BUFFERS_AND_WATER sec. 10's T-flat x/÷1.5 band
   is retired in its favour is for the S3 case.
4. **Kappa's species.** Davidek's ratio tracks [H2PO4-] at an assumed pKa2; Rizzi's pH 7.3 and Mottram's pH 5.6 do not
   test it. Declare total vs species, and the pKa2 at temperature.
5. **K_oil/water for MFT, FFT and the disulfide** is an estimate. FU-F3's reading on partition depends on the MFT /
   FFT difference; FU-F2's does not.
6. **Yiltirak's aliquot.** Ask the authors whether the 30 uL came from the water phase or a mixed sample (CF-031).
7. **Fat and M1.** Capture, oil-phase oxidation and M1's thiol loss law act on the same pool. Build F after M1's
   verdict, with `k_mft_decay`'s status (frozen, banded, C2) declared, or they compensate. R4 is re-tested after M1c.
8. **B54's isotherm.** No pea / soy / gluten sorption isotherm, Tg curve or sugar solubility is on disk; W-R1 and W-R4
   cannot run without them.
9. **"Roos" in the plan row** means Lievonen & Roos 2002 and Miao & Roos 2004 (on disk); neither gives a step rate.

## What should block the freeze (for the review and the coordinator)

1. **The reduced power check (sec. 5.3)** must be built (~1 day) and run, and its demotions written in.
2. **Unsourced widths needing agreement:** the K_oil/water band (no falsification depends on it any longer), R7's K
   band, the >= 20 % run-to-run CV assumed for Davidek, the 20 % on Zyzak's HEPES labels, the kappa(T) band.
3. **How a refusing R0 is scored** on W's units (FU-W2 - W4). Proposed: units R0 refuses are scored among the other
   rivals only, and a term that answers there ships through the acceptance rule's adequacy test, not by "beating" R0.
   This is a coordinator / owner call.
4. **PR #88 (CF-129, CF-130) and PR #96 (S1 - S4) are not merged.** This case cites them as pending. Sec. 4.5 depends
   on S3 / S4. Sec. 2.11's load rule depends on CF-130.
5. **The premortem's confirmation round** on these ten fixes (the M4 precedent: one round plus one confirmation).
6. **Candidates:** the owner decides on downloads (sec. 2.9). A candidate read before the freeze is folded in under its
   fixed role; one read after enters by amendment.
7. **Stale hold-out tags** (Leahy 1989a Table IV; the Yiltirak buffer arms) go to the coordinator's next data PR
   (sec. 4.3); no data file is edited here.
8. **Not blockers, but the owner should know:**
   - On disk, M5 credits no term in any part. It records falsifications only:
     - buffer: R0 (with the pH-only reading) at pH 5 - 7; R2 with K <= 0.1 M on Amadori decay;
     - fat: the direct sum, in both directions; partition alone; supply as written;
     - water: the acrylamide lane's flat band and concentration-only above a_w ~0.9.
   - Sec. 4.5's S3 consequences would refuse or widen shipped outputs once S3 applies: fat-thiol sizes, high-moisture
     acrylamide contrasts, and phosphate-referenced sulfur outputs in unbuffered foods at pH 5 - 7.
   - The deciders are designed experiments, which come last by the owner's ruling (2026-11-05).

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
- **PR #88 folded in.** CF-130's undiluted carriers mis-state the oil volume (F-R1) and the carbonyl load (F-R2, F-R6).
  CF-129's Bi row is not F evidence. No decisive unit moves.
- **PR #77 folded in.** R0 for W is the trunk table inside a_w 0.33 - 1.00, refusing below. Most W units are R0
  refusals, so a scoring rule is needed (blocker 3).
- **Corrections.** Davidek's x2.8 - 3.6 is pH 5 - 6 (pH 7 x1.9 - 2.2; pH 8 x1.24 - 1.27). The Michaelis bound at K
  0.1 M is k(0.5) / k(0.2) <= 1.25, not ~1.1. Bell 1995's row restated: the 11-04 text's "glassy and rubbery give the
  same rates" contradicts the dossier, where the glassy system is the slowest where the polymers straddle the glass.
- **Existing evidence kept, with its status.** Davidek's pH-stat, Bell's ratio, Mottram 2002, the Elmore / Yiltirak
  selectivity and Hofmann dry vs wet are all kept, and all are spent: falsify-only, never credit.
- **Citation sweep run** (three sweeps, 13 seeds, 1,010 unseen works by title / abstract, 14 candidates with roles
  fixed before reading; sec. 2.9). It found no partition coefficient, no protein isotherm and no calibrated thiol
  a_w series.
- **Power check specified** (sec. 5.3), not built. **P(separation)** per part (sec. 5.4).
- **Decisive-measurement entries proposed** in the PR body; the file itself is unchanged.

### Changelog: 2026-11-06 (later), independent premortem "freeze after fixes": ten fixes applied; still READY TO FREEZE, not frozen

Still before any result of the change; no engine run. [dry] arithmetic on printed values only.
1. **S3 consequences** (new sec. 4.5): per output family refused / widened / unchanged, scoped by pH, temperature and
   composition. Fat: thiol direction shown, size refused; the reverse flow added. Buffer: phosphate-referenced outputs
   in unbuffered foods are widened at pH 5 - 7 and unchanged at pH >= 8. Water: R0's acrylamide lane is falsified above
   a_w ~0.9 with common random numbers. Declared: a bounded kappa(T) band, a 0.5 dex zero-buffer intercept band, the
   pH scope, and the sensitivity computation (NumPy switches; cost measured on one timed setting first).
2. **R2 rescoped.** K <= 0.1 M is falsified on Amadori decay at 25 - 37 C (both Zyzak ladders; K = 0.05 M on Bell).
   K = 0.3 M is falsified only anoxically. Further changes:
   - Zyzak's air ladder reported (linear 6.99, K 0.3 M 7.27, K 0.1 M 17.6).
   - The SD-vs-SE choice disclosed.
   - R1 scored with the HEPES point as a sensitivity: falsified if its intercept must equal HEPES.
   - Rival R7 (saturating at ~2 mM plus linear) added.
   - "pH held at every concentration" dropped for Zyzak: Davidek alone falsifies the pH-only reading.
3. **FU-W2 (Leahy) demoted** to development evidence (nonfat dry milk, lactose crystallisation, a per-dry-weight
   recovery basis). Leahy 1989a Table IV and the Yiltirak buffer arms are recorded as stale hold-out tags.
4. **FU-F1's partition claim dropped**; FU-F2 carries partition alone. Every falsification is restated as a
   posterior-predictive p after Holm (sec. 4.0), e.g. Bell at K 0.1 M p 0.12 (not falsified), Zyzak at K 0.5 M Holm
   0.057 (not falsified), Davidek at pH 7 with a run-to-run CV of 20 % p 0.019.
5. **R6's product predictions written.** The E7 mapping is corrected: the N2 arm separates {R2, R6} from {R1, R5}, not R2
   from R6. The R2 / R6 separators are added (an enal without oil; PV-spiked oil vs an enal in oil, both under N2;
   thioethers vs disulfide / oxidised MFT in both phases).
6. **R4 labelled** "as written, at the B60 route structure"; re-test after M1c.
7. **Rivals:**
   - R3 (buffer) merged into R0.
   - Added R7 (saturating plus linear) and R8 (ionic strength).
   - Added R9, SPME fibre competition, as a measurement model for Elmore.
   - Added W-R4, a solubility-capped concentration basis.
8. **FU-W1 stated as basis-independent:** wet basis x2.60 at the wet arm; per-water basis x1.52 at the wet arm but ~x5
   at the dry arm.
9. **FU-F2 normalised for fibre competition.** By 2-methylfuran / hydroxypropanone 0.59 / 0.52 dex (~4 / 3.5 SE). By
   the mean of three markers 0.48 dex, Holm p 9e-5: still falsified. By 2-butanone alone it is marginal (Holm 0.051).
10. **The CF-130 rule widened** to ANY pot. Yiltirak's 10 % canola would be charged ~10x and Elmore's ester ~130x.
- **Slot scope:** the late-November slot is the power check plus the S3 consequence computation. The two-phase and
  O2-per-phase builds are deferred until a decider exists.
