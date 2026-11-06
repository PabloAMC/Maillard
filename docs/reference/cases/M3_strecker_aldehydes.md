# Case M3: Strecker aldehydes

*Plan step M3 ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); error E3 (F16-a, F16-b) of
[ENGINE_VS_LITERATURE](../pathway_atlas/ENGINE_VS_LITERATURE.md) sec. 3, with E4 (F16-c), E15 (F17-a) and CF-077 as
neighbours on the same step). Items 1 - 5 of sec. 7, written 2026-11-04 before any result of the change exists. No
engine and no forced-input layer was run for this case. Items 6 - 9 are appended later and never edited into this
part. The H2 tool bake-off (W4) runs on the layer this case defines.*

**Seen by the author.** No engine output on any held-out Strecker row. Read while writing: the plan's sentence that the
shipped engine "predicts no change on all 20 methionine items" of the W1.1 benchmark (the score file itself was not
opened); the MEASURED values in the atlas text of the v1 TEST cells F16#10, F16#13, F16#14 and of the v2 TEST(E3) cells
F16#11, F16#15 - F16#17; the keying-queue entry for Huang 2017; dossier summaries of the same sources (measured values
only). The forced-input results EX-049 and EX-051 (single laboratory, development rows) were read in full: they shaped
the rivals, and the glucose-pot shapes they used are therefore SEEN, not blind (sec. 4).

**Split notation.** "v1" = [split.yml](../contrast_benchmark/split.yml) (`contrast_split_v1`); "v2" =
[split_v2_per_change.yml](../contrast_benchmark/split_v2_per_change.yml) (`contrast_split_v2_per_change`), error E3.
v2 DEV(E3) labs: DFA Garching, Missouri (Feather), Quest. v2 TEST(E3): F16#11@monforte2021 (Porto), F16#14@huang2017
(TUM Becker), F16#15@balagiannis2009, F16#16@martin2001, F16#17@mandin1999 (Reading); F16#12@chan1994b excluded as a
calibration input. v1 TEST labs among the sources below: Bordeaux (F16#10), Munster / Eichner (F16#13), TUM Becker
(F16#14), Iowa State (F31#0 - 1), SayFood (F06#16, F08#9). Atlas cells are 0-based indices of
[F16_strecker_aldehyde.yml](../pathway_atlas/F16_strecker_aldehyde.yml) and
[F17_strecker_acid_o2.yml](../pathway_atlas/F17_strecker_acid_o2.yml) unless named. Every other cell named here is
v1 DEVELOPMENT.

**The TUM Becker contrast, checked.** `F16#14@huang2017` is TEST in v1 (`partition: test`, lab
`tum_weihenstephan_becker`) and in v2 TEST(E3). It is not keyed (keying queue rank 34; ledger
`huang2017_tab1_tab2_leu_ile_methylbutanal_constants`), and its re-read was single-reader (contrast-benchmark README
sec. 9). This case uses it only in item 7.

## 1. The claim

Strecker aldehydes (phenylacetaldehyde PA first; 3- and 2-methylbutanal, methional as identity-scaled extensions)
form by two routes, and the second route's Amadori compound decays with ONE constant in every pot:

- **Route A, dicarbonyl:** donor (MGO, glyoxal, 3-DG, glucosone) + free amino acid -> aldehyde + aminoketone + CO2.
  The aldehyde is O2-independent; the acid is a parallel exit with an O2-dependent channel per donor class (EX-042's
  arm A').
- **Route B, Amadori oxidative:** the amino-acid Amadori compound (ARP-Phe) -> aldehyde (and acid) with O2.
- **The claim proper (the step EX-051 asked for):** ARP-Phe decays first order with one constant k_d, the same in the
  fed ARP-Phe pot and in the glucose + Phe pot. So route B's rate constant, calibrated on the fed ARP-Phe pot, carries
  to the glucose pot with no free parameter.

Then, as phase 2, the layer is fed from the network instead of the measured courses.

| part | step(s) | today | M3 |
|---|---|---|---|
| L. forced-input layer (standalone, not the engine) | `scripts/generators/feasibility_strecker_forced_input.py`: route A forced by hofmann2000 Fig. 3; route B forced by the ARP-Phe course (hofmann2000b Table 3); `strecker.log10_k_arp_decay` in [-2.2, -1.15] | two-route layer INCONCLUSIVE (EX-049); its optima fail the blind fed-ARP shape (EX-051) | k_d shared between the fed-ARP and the glucose pot, then scored |
| N. network (phase 2, on A2) | new species PHE, PHEARP, PHEAL, PHEAC (provisional keys, prereg sec. 9), LEU / MEBUAL3 / MET; donor + AA steps for MGO, GO, 3-DG, glucosone; ARP-AA -> aldehyde (O2) | no Strecker-aldehyde species. `r_go_ak` / `r_mgo_ak` (B18, [network.py](../../../src/kinetic_core/network.py) ~359 - 365) take glycine only and book the aldehyde to `FRAG_C`. `METHIONINE_SHIPPED = False`. Leucine requests are refused | the layer's steps on the trunk's own courses |

**Earlier waves on these structures** (all cited from the record, none re-run):

| structure | earlier wave | outcome | kind of evidence |
|---|---|---|---|
| route A alone (measured donors) | EX-049 arm DC | not feasible (min violation 3.9; 1.3 without Fig. 2's 10-min point); recovery rejects it 6 / 6 when route B is real | fitted shapes, one lab |
| route B alone | EX-049 arm ARP | fails (PAA never overtakes PA) | fitted shapes |
| A + B, independent k_d | EX-049 DC_ARP; EX-051 | INCONCLUSIVE: passes post hoc at band walls (k_d 0.058 - 0.069 /min, O2 ~1); blind fed-ARP shape fails for all five optima (+0.5 dex at 10 min); Fig. 1's own course needs k_d ~0.010 /min (post hoc) | fitted shapes, one lab |
| A + B, one shared k_d (the claim) | EX-051 follow-up ("Pre-register an Amadori-route variant ... no two-route fit before it") | **never pre-registered or run**; no prereg file exists | -- |
| O2 upstream of the donor (trunk) | B29 (`kinetic_core_b29_prereg.md` sec. 7; ship rule) | DO NOT SHIP: the trunk's oxidative share is smaller on the Amadori pot than on the sugar pot, the wrong order; 0.79x vs 2.63x | whole-model fit |
| Amadori-Met -> methional, first order | B22b | DO NOT SHIP: model peaks at 30 min against 120 | whole-model fit |
| donor + Met / Pro on B18's constants | B22, B24, B24b | DO NOT SHIP (B22: Pan 2025 3.6 - 5.6 dex low, Deng 2022 1.5 - 2.9 dex high) | whole-model fit |
| the B18 Strecker step itself | EX-021; prereg sec. 7 | B18's constant is a pyrazine-yield constant, ~2.3 dex below Hofmann's fed MGO + Phe rate | arithmetic on fed rows |
| a per-AA second-order Strecker constant in water | B19 reading logs (`tasks/data_restructure_plan.md`, 2026-09-08 / 09) | none exists on disk; four laboratories print a barrier or a yield | literature |
| the step's pH law | E4 prereg (`kinetic_core_pyrazine_strecker_step_ph_prereg.md`) | expected REFUSE (P 0.85); declared band | prereg, no fit |

**Errors stack (owner, 2026-11-04; CF-107).** The B22, B22b, B24 and B29 refusals are whole-model fits on the trunk.
They are weak evidence against their structures. Candidate compensating errors, named so that items 6 - 7 can register
them:
- B18's pyrazine-yield Strecker constant (EX-021);
- no 3-DG or glucosone donor (F16-b);
- the trunk's MGO level (B21: 7.8 mM against Xia's <= 0.5 mM);
- the Leahy pH law on the step (E4);
- glycine as the only amine (methionine charged as glycine in B22b);
- the trunk's oxidative vs non-oxidative Amadori branching (B29 sec. 7);
- instant aminoketone condensation (E14);
- the missing dicarbonyl-amine adduct sink (CF-095, >= ~1e3x).

The forced-layer outcomes (EX-049, EX-051) are fits too, of one laboratory. Their compensators:
- forcing calibration (quinoxaline scales; GO sigma >= 0.2 dex);
- vessel O2 (reflux vs closed vial, O2 declared in [0, 1]);
- CF-070's reassignment of the Phe course;
- unmeasured donors (1-DG, C4 / C5 osones, DHA, diacetyl);
- donor speciation at 100 vs ~1 mmol/L;
- the 0 - 10 min heat-up;
- ARP-Phe purity (~90 %).

## 2. All the evidence, by kind

"for A / B / k" = bears on route A, route B, or the shared k_d.

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 2000 Table 5 (F17#2; ledger `hofmann2000_table5_labelled_pa_oxidation`) | [13C2]PA in the glucose pot: < 1 % appears as labelled PAA in 30 min | the acid is a parallel exit, not PA oxidation (fixes the layer's acid structure) | DEV |
| Zhou 2024 (F03#15, F19#3, F19#4, F20#10; CF-087, EX-052) | Ala-xylose ARP + [15N]Ala, 120 C: pyrazine 15N : 14N 1.16 - 1.20, so about half the aminoketone N comes from free Ala on ARP-derived dicarbonyls | route A runs in an ARP pot (aminoketone side, pentose ARP) | DEV (Ho network); v2 calibration input |
| Low 2007 (F19#1, F20#2 - 3) | [2-13C]Gly in potato cake: pyrazines <= 7 % labelled; the aminoketone carbons are sugar dicarbonyls | route A's co-product in a food | DEV (Reading) |
| Fang 2009 (F31#0, F31#1) | 13C / 15N NMR: 27 +- 4 % (dry) and 33 +- 4 % (solution) of glycine decarboxylated in melanoidins | total Strecker-type flux, not the route | **v1 TEST** (Iowa State): listed, not used |
| Mundt 2004 (F31#3) | 14C glucose + glycine: whole : decarboxylated glycine 0.289 : 0.662 per glucose | total flux | DEV (Leeds) |
| Blank 1996 (blank1996 dossier) | [13C] labels: the Strecker aldehyde (HCHO / acetaldehyde) is consumed by 1-deoxypentosone into furanones | an aldehyde sink (Gly / Ala aldehydes only) | DEV (Nestle) |
| none | no label traces a Strecker aldehyde to its dicarbonyl donor vs the ARP | the route-A / route-B split | gap |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 2000 Table 4 (F16#0, F16#1; F17#0, F17#1; `strecker_hofmann2000.yml`) | 0.1 M donor + Phe, pH 7, 98 C, 60 min. PA argon / air: MGO 11.1 / 10.2, GO 9.2 / 8.3, 3-DG 8.4 / 7.8, glucosone 8.1 / 6.9 umol/mmol; PAA x4.34 / 5.5 / 2.2 / 1.6 in air | route A per donor; aldehyde O2-independent (x0.85 - 0.93) | DEV; v2: data copy of a calibration input |
| Hofmann 2000 Table 3 (F16#2, F17#5) | fed MGO + Phe 10 - 500 min: PA 6.2 -> 13.5 umol, plateau 120 - 240 min | route A time course; MGO lifetime ~15 - 25 min | DEV |
| Hofmann 2000 Cu arm (F17#3) | air + Cu: PAA x1.33, PA x0.98 | route A acid channel | DEV |
| Hofmann & Schieberle 2000b Table 2 (F17#6; ledger `hofmann2000b_table2_arp_phe_oxygen`, a calibration input via B29) | fed ARP-Phe 0.1 M, 100 C, 120 min, closed vial: PA 0.6 / 5.5 / 13.8, PAA 0.3 / 3.0 / 7.6 (argon / air / air + Cu); PA / PAA 1.8 - 2.0 in every atmosphere | **route B exists; air / argon 9.2** | DEV |
| Hofmann & Schieberle 2000b Fig. 1 circles (`strecker_hofmann2000b_fig1.yml`) | fed ARP-Phe PA, PAA 10 - 500 min, reflux, air | **k_d of the fed pot** (in-sample after EX-051) | DEV |
| Weenen 2001 Table III (F16#5; CF-073) | 15 carbonyls + Phe, pH 3.2, SDE: MGO 21.0 %, 3-deoxyerythrosone 17.2, GO 8.8, 3-deoxyxylosone 2.6, 3-DG 0.6 % | donor ranking (agrees); magnitude x14 - 20 above DFA for MGO (CF-073: per-donor MGO +1.24, GO +0.89, 3-DG -0.20 dex; a donor-generic bias rejected) | DEV (Quest) |
| Weenen 2001 Table II (F17#9) | fed Glc-Phe ARP vs glucose + Phe: x4.2 at pH 7.0, x267 at pH 3.2 | route B or ARP -> donor; cannot separate them | DEV |
| Granvogl 2006 (F16#7) | MGO + Leu / Ile / Phe / Met, 5 mM, 170 C, pH 7: aldehyde 2.25 - 12.2 mol %; Strecker amine <= 1.1 % (3 % Met); Ile 1.68 : Leu 1.16 : Phe 1 : Met 0.31 | identity at fixed donor; no amine exit in water | DEV (DFA network; not cross-lab) |
| Ghiron 1988 (F16#9) | 3-DG + Phe gives PA (unquantified); 3-DG loss x3 faster with Ala | the 3-DG donor in a second lab | DEV (Missouri) |
| Munch 1998 (munch1998 dossier; ingested 11-02, not an atlas cell) | 145 C, pH 7: aldehyde from MGO vs glucose, Leu 4.3 vs 0.4, Phe 1.0 vs 0.3 mol %; Phe + MGO PA plateau from 90 min | route A at a second temperature, same lab | DFA (a v1 DEV lab) |
| Deng 2022 (F17#10) | fed Met-ARP 0.2 M vs Met + Glc, unbuffered, initial pH 7.5, 120 C: methional ARP / sugar 1.37 -> 2.59 -> 1.10 over 30 - 240 min; reversed at 130 C; GO, MGO, diacetyl in the same pots (figure only); pH falls faster in the ARP pot | **the only second-lab fed-ARP vs sugar pair** | DEV (Ho network) |
| Voigt 2010 (F10#1, F12#6 - 7) | fed 1-DG + Ala, 37 / 50 C, pH 7.4: acetaldehyde x1.5 - 1.6 aerated (the deaerated arm carries DTPA) | O2 can raise the aldehyde on a dicarbonyl route too, so the air / argon ratio is not by itself a route marker | DEV (Glomb) |
| Vanlancker 2010 (F19#13, F20#7 - 8) | Lys + X or Lys-X dipeptides + glucose / MGO / GO, 130 C: 3-MB and PA need the free amino acid | route A needs a free alpha-amino acid | DEV (Ghent) |
| Piloty 1979 (piloty1979 dossier; CF-095) | fed MGO + 13 amino acids, 100 C: Phe 33 %, Leu 26 % lost in 3 h, pH-flat; no Strecker product measured | the non-Strecker amine sink (the layer's nu) | not an atlas cell |
| Pripis-Nicolau 2000 (F16#10) | aldehydes at 25 C, pH 3.5, in weeks; existence only | existence at low T | **v1 TEST** (Bordeaux): listed, not used |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann & Schieberle 2000b Figs. 2 - 3 (F17#8; ledger `hofmann2000b_fig2_fig3_arp_phe_dicarbonyls`; **not keyed**) | ARP-Phe 0.5 M, 85 C, in-situ OPD: argon -> 1-deoxyosone the main dicarbonyl; air + Cu 25 mM -> glucosone the main one; Phe liberated ~55 % at 320 min in both | **route B vs ARP -> osone -> route A** (sec. 4, D3) | DEV |
| Nedvidek 1992 (F06#2, F07#3, F15#0) | in-situ OPD, glucose / xylose + beta-Ala, pH 5 - 7: 1,4-dideoxyosones much higher with alpha-amino acids (Strecker on the 1-deoxyosone) | 1-DG as an unmeasured Strecker donor | DEV (Ledl school) |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann 2000 Fig. 3 (`strecker_hofmann2000.yml`) | 3-DG, GO, MGO, glucosone in the glucose + Phe pot, 10 - 500 min (1-DG not measured) | route A's forcing | DEV |
| Hofmann & Schieberle 2000b Table 3 (`strecker_hofmann2000b.yml`) | ARP-Phe in the glucose pot: 42 / 48 / 32 / 21 / 9 umol at 10 / 30 / 60 / 120 / 300 min; decay k_d >= ln(48 / 32) / 30 = 0.0135 /min (0.0126 - 0.0144 within print rounding; a lower bound whatever the formation, since formation only slows the net fall; corrected 2026-11-05 from 0.0062, CF-121) | route B's forcing; **k_d's floor in the glucose pot** | DEV |
| Deng 2021 (F03#11 - 12, F06#6, F07#11, F08#5, F09#7) | fed Met-ARP 0.2 M, 100 C: k 0.0048 /min (pH 5.5), 0.012 (pH 7.5), unbuffered; phosphate ~x1.8; added Met leaves k unchanged (+-15 %); 88 % of reacted ARP returns free Met; glucosone up to 10.6 mM at pH 7.5 | k_d of an amino-acid ARP in a second lab; k_d amine-independent; pH / buffer move k_d | DEV (Ho network) |
| Martins 2003, 2003b (F03#3 - 6) with Martins 2005b (F02#2, F06#13, F07#15) | fed DFG 10 mM, 100 C, pH 6.8: k1 + k2 + k3 = 3.68e-2 /min; in the glucose + glycine pot, k4 + k6 + k7 = 2.47e-2 (x1.5). At pH 5.5, x1.07 (the screen's arithmetic on the same tables, not re-derived here). Different network topologies and charges | **fed vs in-pot ARP decay in one lab** (glycine, not Phe) | DEV (Wageningen); v2 calibration inputs |
| Wei 2022 (F03#16, F04#23) | GX-ARP 5 mM: decay second order, slowed ~15 % by Cys | the order of the decay law | DEV (Ho network) |
| Guan 2012 (F01#9, F03#10) | proline ARP, 110 C, 60 min remaining: phosphate 3.2 %, acetate 46.7, citrate 40.7, water ~51 % | k_d depends on buffer species | DEV (SCUT) |
| Lee 2022 (F06#16, F08#9) with Lee 2024 (3-MB panel not keyed; `heldout_lee2024_fig2.yml` holds only the pyrazine panel) | cake: 3-DG ~11x higher with Leu; 3-MB only in glucose + Leu | donors and aldehyde in one bake | **v1 TEST** (SayFood) / keyed held-out (not in a frozen split): listed, not used |
| Pan 2025 (pan2025 dossier; not an atlas cell) | Met + sugars, 50 mM citrate pH 6.2, 100 - 140 C: methional rate constants (unit unprinted); 3-DG 0.8 and MGO 0.45 ug/mL at one time (120 C, 4 min) | a candidate forced pot for methional; one donor time point | B22 used its rates (refused) |
| Gigl 2021 (F31#16 - 17) | melanoidin binding of 3-MB: 99.4 % free at 60 min; hours to days | the aldehyde sink is slow at <= 98 C over 2 h | DEV (TUM Freising = DFA) |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Hofmann & Schieberle 2000b Table 2 (F17#6) | fed ARP-Phe: argon PA is 0.6 / 5.5 = 11 % of air PA, so >= ~89 % of the fed pot's aldehyde in air needs O2 (arithmetic on printed values) | route B, or an O2-made donor (sec. 4, D3) | DEV |
| Cremer 2000 (F16#13; CF-020) | glucose + Leu, a_w 0.52, 90 C: Leu + Fru-Leu + 3-MB ~100 % to 120 min; 3-MB linear while Fru-Leu is steady | the aldehyde flux follows the ARP pool (either route) | **v1 TEST** (Munster): listed, not used |
| Balagiannis 2009 (F16#15) | liver extract, 120 - 140 C: < 5 % of reacted Leu / Ile become aldehyde; the aldehyde tracks glucose, not free amino acid | supply-limited | v1 DEV; **v2 TEST(E3)**: listed, not used |
| Hofmann 2000 Table 3 Phe column (F16#2; CF-070) | the "(PA + PAA) is 18 % -> 5.8 % of Phe consumed" share uses a Phe series CF-070 confirms is reprinted (P 0.9), probably from the glucose pot (P 0.7) | **voids F16#2's share and the review's "nu = 1 refuted"**: the MGO pot has no Phe measurement | DEV |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Hofmann 2000 Figs. 4 - 5 (F16#3, F17#4) | pH 3 - 9: PA largest at pH 5 for every donor; PAA rises with pH | route A's pH law (E4); declared conditions (60 min, air) | DEV |
| Hofmann 2000 Table 2 (F16#4) | glucose + Phe / Leu / Ala, 30 min: Leu / Phe 2.1 | identity; no route marker | DEV |
| Hofmann & Schieberle 2000b Table 2, glucose pot (F17#7; CF-074) | argon / air / air + Cu: PA 0.4 / 1.4 / 2.6, PAA 0.2 / 1.8 / 3.9; PA / PAA 2.0 -> 0.78 -> 0.67 | the glucose pot's O2 response; PA rows were B29 fit rows (diagnostic only); possible vessel comparison (CF-074) | DEV |
| Weenen 2001 Table II (F16#6; CF-020) | glucose + Phe pH 7.0 vs 3.2: 0.15 vs 0.006 % (buffer co-varies) | supply's pH, not the step's | DEV |
| Granvogl 2006 model A vs B (F16#8) | silica, 10 % water: Leu / Ile aldehyde x2 - 3; Strecker amine 19 % (Phe), 27 % (Met) | water activity (M5) | DEV |
| Chan 1994b (F16#12), Chan 1994 | glucose + Met / Phe / Leu / Pro, pH 6 - 8, 75 - 115 C, N2: apparent Ea 3-MB 80.3, PA 90.0 kJ/mol; PA lost to aldol with 3-MB | supply barriers; an aldol sink | v1 DEV; v2 excluded (calibration input) |
| Martin 2001 (F16#16) | fried potato: six amino acids on one glucose pool: Leu aldehyde x1.24, Ile / Phe / Met /2.4 / 2.7 / 5.2 | competition for one donor pool | v1 DEV; **v2 TEST(E3)**: listed, not used |
| Mandin 1999 (F16#17) | Met +- glucose, pH 6, SDE: no methional from Met alone | route A needs a sugar carbonyl | v1 DEV; **v2 TEST(E3)**: listed, not used |
| Monforte 2021 (F16#11) | Phe + gallic-acid quinone, pH 3.4: pH 7 > 5 > 3.4; metals x4 | a quinone donor, not a sugar route | v1 DEV; **v2 TEST(E3)**: listed, not used |
| Huang 2017 (F16#14; CF-018: huang2016 is the same paper) | glucose / maltose + Leu / Ile, pH 5.2, 90 - 130 C: 3-MB x~100; sugar type matters more than amino acid | supply barrier | **v1 TEST, v2 TEST(E3)**: item 7 only |
| Meynier 1995 Table 4 (`heldout_meynier1995_tab1_2_4.yml`) | methional vs pH 4.5 - 6.5, ribose + Met | methional pH | keyed held-out (not in a frozen split): item 7 only |
| Kocadagli 2021 (F29#23), Low 2006 (F29#22, F19#1) | potato dough / cake: added amino acids raise all aldehydes 3 - 5x, less when Gly + Ala compete; glycine lowers them 30 - 70 % | a shared carbonyl pool in foods | DEV (Reading) |
| Cheng 2020 | methional needs a sugar or ascorbate | route A qualitatively | no atlas cell |
| food-matrix factor studies: Solina 2007, Wang 2025, Alim 2018, de Sousa Fontes 2024, Devaere 2022, Goncuoglu 2026, Zhang 2023, Zhang 2023b, Yu 1995, Munch 1997 | levers in real matrices (sugar, pH, heme, protein, temperature) on headspace aldehyde levels, response factor 1 or FD factors | W5 and item 7; no route marker, supply-confounded | no atlas cell (Alim 2018: BTBU, a v1 DEV lab) |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Piornos 2025 (F02#15) | malt kilning, 7 ARPs measured, donor pool not measured: one shared ARP-decay constant for all ARPs; a direct ARP -> aldehyde route rejected on AIC (672 vs 639) | **against route B, for one k_d**, both from a fit with an unmeasured pool | DEV (Reading) |
| Balagiannis 2009 (F16#15), 2010 | glucose -> pool -> aldehyde; Strecker step declared fast; Ile : Leu 1.64 per mole; 3-MB sink 2.72e-3 /min at 130 C | supply-limited; the sink | 2009: v2 TEST(E3); 2010: no atlas cell |
| Cremer 2000 (F16#13) | Ea 115 - 124 kJ/mol, the supply's | supply barrier | v1 TEST |
| Desclaux 2006 (F06#14, F07#17) | xylose-Gly ARP -> 1-DX 2.79e-2, -> 3-DX 5.50e-4 /min, 100 C | pentose ARP decay | DEV |
| Ge 1997 (F01#0 - 1, F02#6 - 7) | Phe-ARP formation; quarantined (method disputed) | supply only | DEV |
| engine and layer records | EX-021, EX-042, EX-049, EX-051; B18 (shipped), B19, B22, B22b, B24, B24b, B29; E4 prereg; EX-052, EX-057 (CF-077), EX-058 (CF-095) | sec. 1 table | engine records |

**Left out, with reasons.**
- Lipid-carbonyl Strecker (Zamora / Hidalgo 2004 - 2016, Ruan 2025, Yeo 2022, Kim 1996b): E20 / F18, not M3.
- Cysteine Strecker and H2S (Kobayasi 1965, Du 2023, Cao 2017 and others): E19 / M2; the sulfur lane's `r_dpo_ddp`.
- Pyrazine-only labelling and kinetics (Amrani-Hemaimi 1995, Hwang 1995a / b, Leahy 1989, Shibamoto 1977 and others): M6.
- Tressl 1994: methional's label is not reported.
- Thiol-aldehyde adducts (Lienhard 1966, Wondrak 2002; 25 - 37 C): a cysteine-pot sink, M2.
- Elmore 1997, Schutte 1972, Yao 2025: products of the aldehyde, no rate at 98 - 145 C.
- Huang 2016: a duplicate of Huang 2017 (CF-018).
- The remaining ~120 dossiers that mention Strecker were screened (two read-only passes) and mention it in passing, with no data on formation, donors or the ARP.

**Counts by kind (sources).** Isotope labelling 6 (none on the route split); fed intermediates 16 (route B: Hofmann
2000b, Weenen 2001 and Deng 2022 only); trapping 2; measured intermediates and lifetimes 10; mass balance 4 (one of them
voided by CF-070); one factor within a study ~25; fitted time courses 5 plus the engine records. Nine of these sources
are held out (v1 TEST: Fang, Pripis-Nicolau, Lee, Cremer; v2 TEST(E3): Balagiannis 2009, Martin, Mandin, Monforte; Huang
in both), and Meynier Table 4 is keyed held-out. They are listed, never used to choose a rival.

## 3. The rivals

| id | explanation | fed ARP-Phe pot | glucose + Phe pot |
|---|---|---|---|
| R0 | **current engine** | refused (no species) | refused. B18 makes AKG / AKM + `FRAG_C` from GO / MGO + Gly only; the W1.1 methionine items are all "no change" |
| R1 | **the claim**: routes A + B, one shared k_d | PA from route B, PA ~ integral of [ARP]; k_d set by this pot | route B's constant carried in; no free route-B parameter |
| R2 | **route A only, with donors the forcing misses**: the ARP's O2 response is O2 making an osone (glucosone; 1-DG under argon) that then does route A with released Phe | PA needs released Phe x a newly made osone | the residual over the four measured donors comes from 1-DG, C4 / C5 osones or diacetyl, not from route B |
| R3 | **routes A + B with a pot-specific k_d** (the EX-049 optima) | k_d ~0.010 (Fig. 1) | k_d 0.058 - 0.069 (band wall); the fed pot does not transfer (vessel, purity, concentration) |
| R4 | **R1 with k_d moving with pH, buffer and order** (Deng 2021, Guan 2012, Wei 2022) | as R1 at pH 7, 0.5 M phosphate | as R1 for Hofmann; differs where pH drifts or the buffer changes (Deng 2022) |

R0 makes no aldehyde prediction, so its stand-ins on the decisive data are the prereg's two nulls (no change and
typical yield; prereg sec. 6.3). R1 against R3 is the claim's test. R1 against R2 is the route question that both
rest on. R4 matters only beyond the Hofmann pots.

## 4. The decisive data

Directions follow from each structure. "[dry]" marks arithmetic on printed values or on the layer's calibration
identity, not a run.

| # | measurement (record) | R1 | R2 | R3 | R4 | separates |
|---|---|---|---|---|---|---|
| D1 | **early shape of the fed ARP-Phe aldehyde**, PA(10, 30, 60) / PA(120): hofmann2000b Fig. 1 circles | first order from t = 0: concave, no lag | needs released Phe and a new osone at once, so PA rises with ~zero initial slope (a lag) | as R1 in this pot | as R1 | R1 / R3 / R4 from R2. In-sample for k_d (EX-051) but never scored for the lag |
| D2 | **fed vs in-pot ARP decay, measured**: Martins 2003 vs 2005b (glycine, same lab); Deng 2021 (Met-ARP, +- free Met); the Hofmann glucose-pot floor k_d >= 0.0126 - 0.0144 /min (corrected 2026-11-05 from 0.0062, CF-121) | fed = in-pot within noise; k_d at 100 C, pH ~7 in the range those measurements span | -- | in-pot k_d several times the fed one | k_d scales with pH and buffer | R1 from R3 on decay alone. Martins gives x1.5 at pH 6.8. **No measured amino-acid ARP decay at 100 C on disk reaches R3's 0.058 - 0.069 /min** (fastest: Martins DFG 0.037 at 0.1 M phosphate). Hofmann runs 0.5 M phosphate, and E5's linear phosphate law (EX-058, P 0.8) would raise its k_d above 0.1 M values, so D2 bounds the transfer, not the absolute k_d |
| D3 | **the fed ARP's O2 mass balance**: argon PA <= 11 % of air PA (Table 2); trapped osones in argon vs air + Cu (Figs. 2 - 3, unkeyed); Deng 2021's standing glucosone | the air excess comes without a matching free osone | [dry] needs >= ~0.5 - 0.7 mol free C2 - C6 osone per mol ARP within 120 min (4.9 umol/mmol of excess PA over the fed donors' 6.9 - 10.2), more with less free Phe | as R1 | as R1 | R1 from R2, nearest the mechanism. Deng's standing glucosone is ~5 % of the ARP, but a standing level is a floor on flux, not the flux |
| D4 | **PA / PAA fingerprints (scale-free)**: fed ARP 1.8 - 2.0 in argon, air and air + Cu while its trapped donor switches 1-DG -> glucosone; fed donors in argon 2.4 (MGO), 4.2 (GO), 5.6 (3-DG), 6.2 (glucosone); glucose pot 2.0 (argon), 0.78 (air) | the ARP fingerprint is one route's, atmosphere-free | the fingerprint follows the dominant osone, so it should move when the donor switches, unless 1-DG's fingerprint equals glucosone's | as R1 | as R1 | R1 from R2. **1-DG has never been fed with Phe**, so its fingerprint is unknown. The glucose pot's argon 2.0 lies below every measured donor's argon fingerprint, but at the printed precision (0.4 / 0.2, +-0.05) its range is 1.4 - 3.0 |
| D5 | **glucose-pot shapes**: hofmann2000 Fig. 2 (PA = PAA crossover, PAA / PA at 500 min) and hofmann2000b Fig. 1 squares (G5') | [dry] route B's constant falls by about the ratio of the fed ARP's time-integrals, ~x4 from k_d 0.065 to 0.010, so route B's share falls below EX-049's 45 - 65 % at 120 min and part of the misfit that rejected route A alone returns | the residual tracks an unmeasured donor's course (1-DG builds later, Fig. 2 of 2000b) | fits, at band walls (EX-049) | as R1 | R1 from R3. **SEEN and in-sample** (EX-049 used Fig. 2; EX-051 used Fig. 1). No blind development glucose-pot shape remains on disk |
| D6 | **second laboratory**: Deng 2022 methional, ARP / sugar ratio vs time (1.37 -> 2.59 -> 1.10, 120 C) with Deng 2021's k_d and in-pot pH | the ratio's timing from k_d and the measured ARP | the ratio from donor courses only | free | the ratio moves with the pot's pH drift | R1 / R4 from R3 across labs. Needs a methional sink (B22) and a pH course; semi-quantitative (RF 1), donors figure-only. **Phase 2 only** (corrected 2026-11-05): phase 1's layer is forced by Hofmann's Phe courses and has no Met or methional species (MET arrives with part N, sec. 1), so D6 cannot be scored before the network phase |

**The plan's decisive data, checked.**
- "Hofmann 2000 rows": Tables 3 - 4 (and the Table 2 fed ARP) are exact calibrations in every rival (EX-049: "fed-pot
  calibrations exact"), so they cannot discriminate. Only the glucose-pot shapes (D5) discriminate, and they are seen.
- "The held-out TUM Becker contrast": TEST in v1 and v2(E3). It is a sugar pot with no donor or ARP measured, absolute
  concentrations unprinted, and Tables 1 and 2 disagreeing. It sits at pH 5.2 and 90 - 130 C, outside the layer's
  scope (pH 7, 98 - 100 C). Its x~100 is a supply barrier that the trunk's errors also move. Under sec. 7 item 4 it is
  not decisive. By the M1b precedent a held-out contrast is read only in item 7.

**What the data cannot decide.**
- R1 against R2 rests on one laboratory (DFA). D3 needs Figs. 2 - 3 keyed, and those have an air + Cu arm, not plain
  air, at 85 C and 0.5 M. D4 needs a fed 1-DG + Phe pot that does not exist. Piornos's rejection of a direct route is
  an AIC on an unmeasured pool. Voigt shows O2 also raises a dicarbonyl-route aldehyde.
- R1 against R3 rests on D5 (seen, single lab) and on D2 (other amines, other buffers).
- R4 cannot be separated from R1 on Hofmann's pots, which are all 0.5 M phosphate, pH 7.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1, D5 Hofmann 2000b Fig. 1 | two blind readers, agree within 0.011 umol/mmol (~0.004 dex); the 10-min squares 0.126 / 0.127 dex. No replicates printed. Run-to-run for this recipe: Table 2 vs Fig. 1 at 120 min differ 7 - 16 % (0.03 - 0.065 dex) | `strecker_hofmann2000b_fig1.yml` header; CF-075 |
| D5 Hofmann 2000 Fig. 2 | two-reader read errors per row; within-paper near-replicate pairs give PA ~0.007, PAA ~0.06 dex per measurement (pair RMS / sqrt 2); run scatter of the recipe 7 - 18 % | `strecker_hofmann2000_fig2.yml`; prereg sec. 5; CF-074 |
| D1 - D5 forcing | quinoxaline calibration, no response factor printed: per-donor scale N(0, 0.10) dex (3-DG, glucosone), N(0, 0.15) (MGO), >= 0.20 (GO; text vs figure differ 0.23 dex) | prereg sec. 4 and review |
| D2 Martins 2003 / 2005b | 95 % HPD on each constant; data "at least duplicate" | martins2003, martins2005b dossiers |
| D2 Deng 2021 | triplicates, mean +- SD; digitising error +-3 mmol/L; k is the screen's read of Fig. 3 | deng2021 dossier |
| D2 Hofmann 2000b Table 3 | none printed; integers in umol | hofmann2000b dossier |
| D3 / D4 Hofmann 2000b Table 2 | none printed; values to 0.1 umol/mmol (+-0.05: 0.04 - 0.11 dex on the argon cells); keyed read_error per row | `strecker_hofmann2000b.yml` |
| D3 Figs. 2 - 3 | not keyed; no replicates | hofmann2000b dossier |
| D4 Hofmann 2000 Table 4 | "Replicates not stated anywhere; no error bars or SDs printed"; the dossier treats ratios below ~1.3x as directional | hofmann2000 dossier |
| D6 Deng 2022 | triplicates, mean +- SD; response factor 1 (ratios only) | deng2022 dossier |
| across laboratories (D2, D6) | between-lab spread from H0 (EX-059, pending) | -- |

The EX-049 tolerances (T1 0.151, G5 0.154 dex) were calibrated by an oracle on the model at "realistic noise". They
are a power statement, not a measurement's noise, and are not used here.

## Open questions for the build

1. **Pre-register the shared-k_d step first.** EX-051 asked for it; no prereg exists. It must fix:
   - which rows set k_d (Fig. 1 circles alone, or with the Table 2 fed rows, which are a different run, CF-075);
   - O2 for the reflux pots: held at 1 (EX-051 rule 5) or free. A free O2 can absorb the very difference R1 against R3
     tests, so it is a candidate compensator;
   - the 0 - 10 min heat-up (the point that failed DC_ARP in EX-049).
2. **Key hofmann2000b Figs. 2 - 3** with two readers before any D3 score. They are the only trapped-donor data in an
   ARP pot.
3. **D1's lag test needs a declared R2 rate law** (released-Phe kinetics: Phe liberation ~55 % at 320 min at 85 C is
   the only datum; Deng's 88 % return is Met at 100 C).
4. **CF-070.** The glucose pot's [Phe](t) is 2000b Table 3 (P 0.7). Run the sensitivity with Phe held constant, and
   drop F16#2's "18 % -> 5.8 %" share and the review's nu = 1 refutation from any argument.
5. **Phase 2 shares a step with pyrazines.** If the aldehyde comes off `r_go_ak` / `r_mgo_ak`:
   - their constants must rise ~2 dex (EX-021);
   - the E4 pH law (x0.091 at pH 5) contradicts DFA's pH-5 maximum;
   - B18's pH slopes are conditional on CF-077 (EX-057).
   Decide whether the layer keeps its own constants until W3. If shipped pyrazines move, ask the owner first.
6. **Phase 2 inherits the trunk's supply errors:** B29's Amadori branching, E2, the MGO level, and the missing 3-DG /
   glucosone donors. Feeding the layer from the network is then a test of the trunk, not of M3. Score phase 2 on its
   own decisive data.
7. **The plan row's held-out sources.** Cremer (v1 TEST) and Balagiannis 2009 (v2 TEST(E3)) are named as evidence and
   TUM Becker as decisive. All three should move to item 7, and the row should be edited.
8. **H2.** The bake-off's "deliberately wrong structure" can be route A alone (rejected 6 / 6 on recovery in EX-049).
   Its positive control is R1 on synthetic data from the measured forcing.
9. **CF-073.** Sec. 5 of ENGINE_VS_LITERATURE says "3-DG band DFA-only (Quest 0.03x)". Cluster 4 found Weenen's
   absolute 3-DG within 0.2 dex of DFA; the 0.03x is MGO-driven. Use the 3-DG donor as cross-lab supported and the
   MGO / GO magnitudes as conflicted.

Seen by the author: as stated at the top. No engine or layer output on any row named in sec. 4 was computed for this
case; EX-049 / EX-051 outputs on D1 and D5 rows were read.
