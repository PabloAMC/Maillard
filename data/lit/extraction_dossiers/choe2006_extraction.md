# Choe & Min 2006 — EXTRACTION (REVIEW: autoxidation vs photosensitised oxidation of edible oil; hydroperoxide isomer distributions, beta-scission rules, 2-pentylfuran route, and the primary sources behind each)
### A review, so a shorter dossier: what is drawn or tabulated, and which primary paper the review says it comes from. Nothing here is a measurement made by the authors.

**Source on disk:** `data/articles/choe2006.pdf` (owner's download, 2026-09-08). Read from the `pypdf`
text layer (18 pages); Tables 1-3 are clean and re-typed; Figures 3, 4, 7, 8, 9 were read from
renderings of pages 3-8 and are described in words. The second half of the review (factors: metals,
phospholipids, chlorophyll, antioxidants, processing) was grepped, not re-typed; it contains no
mechanism relevant to the three refused routes.

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanisms and Factors for Edible Oil Oxidation" |
| Authors | Eunok Choe (Inha Univ.), David B. Min (Ohio State Univ.) |
| Venue | Comprehensive Reviews in Food Science and Food Safety 2006, 5, 169-186 |
| DOI | 10.1111/j.1541-4337.2006.00009.x |
| Type | narrative review; all numbers are cited from primary sources, chiefly Frankel 1985 |
| Cited by | cao2020, chen2017, yang2024 (all three cite this review for the four oleate hydroperoxides and for the beta-scission rule) |

## 1. Why it matters

Three of the papers in this batch cite this review as their mechanistic authority. It carries (i) the
autoxidation hydroperoxide isomer distributions of oleate, linoleate and linolenate (Table 1, Frankel
1985), (ii) the singlet-oxygen distributions (Table 3, Frankel 1985), (iii) the product slate per fatty
acid methyl ester (Table 2, Frankel 1985), (iv) the generic alkoxyl beta-scission scheme with the A / B
branches and the fates of the radical fragments (Figure 4), (v) the ene-reaction hydroperoxide
formation for oleate and linoleate (Figures 7, 8), and (vi) the 2-pentylfuran route from linoleate via
the singlet-oxygen 10-hydroperoxide -> 3-nonenal -> 4-oxononanal (Figure 9, Min et al. 2003). For the
rule writer this is the place that says which primary paper backs each step; the repo should cite the
primary where possible and this review as the map.

## 2. Methods as they matter to a model

None (review). Provenance of the mechanistic content as stated by the review:

| content | primary source named by the review | on disk? |
|---|---|---|
| Table 1 (autoxidation isomer %), Table 2 (products of FAME autoxidation), Table 3 (1O2 isomer %) | Frankel EN. 1985. "Chemistry of autoxidation: mechanism, products and flavor significance." In Min & Smouse (eds) Flavor chemistry of fats and oils, AOCS, p 1-34 | no |
| C-H bond energies (C11 of linoleate 50 kcal/mol; C8/C14 75; C17/C18 ~100) | Min DB, Boff JM. 2002. Lipid oxidation of edible oil. In Akoh & Min (eds) Food lipids, Dekker, p 335-63 | no |
| O-O homolysis preferred over O-H ("activation energy to cleave the oxygen-oxygen bond is 46 kcal/mol lower than that to cleave the oxygen-hydrogen bond") | Hiatt R, Mill T, Irwin KC, Mayo TR, Gould CW, Castleman JK. 1968. J Org Chem 33:1416-41 | no |
| Figure 4 (beta-scission A/B scheme) | not attributed to a single source (general) | - |
| ene reaction, Figures 7, 8 | Gollnick K. 1978. In Ranby & Rabek (eds) Singlet oxygen, Wiley, p 111-34 | no |
| Figure 9 / 10 (2-pentylfuran, 2-pentenylfuran by 1O2) | Min DB, Callison AL, Lee HO. 2003. J Food Sci 68:1175-8 | no |
| 2-pentylfuran as soybean reversion flavour | Chang SS et al. 1966. Chem Ind 11:1926-7; Smouse TH, Chang SS. 1967. JAOCS 44:509-14; Ho CT, Smagula MS, Chang SS. 1978. JAOCS 55:233-7; Smagula MS, Ho CT, Chang SS. 1979. JAOCS 56:516-9; Chang SS et al. 1983. JAOCS 60:553-7 | no |
| SPME volatile shares in soybean / corn oil at PV 5 | Steenson DFM, Lee JH, Min DB. 2002. J Food Sci 67:71-6 | no |
| pentane, hexanal, propenal, 2,4-decadienal high in canola stored at 60 C | Vaisey-Genser M et al. 1999. Canola Council report | no |
| relative autoxidation rates oleate:linoleate:linolenate 1:40-50:100 (O2 uptake) | Min DB, Bradley GD. 1992. Wiley encyclopedia of food science and technology, p 828-32 | no |
| 1O2 rate constants with stearic/oleic/linoleic/linolenic acid | Vever-Bizet C et al. 1989 (not re-typed here) | no |
| linoleate + 1O2 1450x faster than + 3O2 | Rawls HR, Van Santen PJ. 1970 | no |
| hexanal, pentane, 2,4-decadienal as oxidation indicators | Warner K et al. 1978. JAOCS 55:252-6; Przybylski & Eskin 1995; Choe 1997; Heinonen et al. 1997 | no |

## 3. Tables re-typed

### Table 1 — "Hydroperoxides of fatty acids by autoxidation" (a: Frankel 1985)

| fatty acid | hydroperoxide at | relative amount (%) |
|---|---|---:|
| oleic acid | C8 | 26-28 |
| | C9 | 22-25 |
| | C10 | 22-24 |
| | C11 | 26-28 |
| linoleic acid | C9 | 48-53 |
| | C13 | 48-53 |
| linolenic acid | C9 | 28-35 |
| | C12 | 8-13 |
| | C13 | 10-13 |
| | C16 | 28-35 |

### Table 3 — "Hydroperoxides of fatty acids by singlet oxygen oxidation" (a: Frankel 1985)

| fatty acid | hydroperoxide at | relative amount (%) | type |
|---|---|---:|---|
| oleic acid | C9 | 48 | |
| | C10 | 52 | |
| linoleic acid | C9 | 32 | conjugated |
| | C10 | 17 | nonconjugated |
| | C12 | 17 | nonconjugated |
| | C13 | 34 | conjugated |
| linolenic acid | C9 | 23 | conjugated |
| | C10 | 13 | nonconjugated |
| | C12 | 12 | conjugated |
| | C13 | 14 | conjugated |
| | C15 | 13 | nonconjugated |
| | C16 | 25 | conjugated |

Text: "Production of nonconjugated hydroperoxides is not observed in the autoxidation." "Autoxidation of
linoleic and linolenic acids produces only conjugated products."

### Table 2 — "Secondary oxidation products of fatty acid methyl ester by autoxidation" (a: Frankel 1985; b: footnote marker on octanal, text not printed)

| class | oleic acid | linoleic acid | linolenic acid |
|---|---|---|---|
| aldehydes | octanal (b), nonanal, 2-decenal, decanal | pentanal, hexanal, 2-octenal, 2-nonenal, 2,4-decadienal | propanal, butanal, 2-butenal, 2-pentenal, 2-hexenal, 3,6-nonadienal, decatrienal |
| carboxylic acid (esters) | methyl heptanoate, methyl octanoate, methyl 8-oxooctanoate, methyl 9-oxononanoate, methyl 10-oxodecanoate, methyl 10-oxo-8-decenoate, methyl 11-oxo-9-undecenoate | methyl heptanoate, methyl octanoate, methyl 8-oxooctanoate, methyl 9-oxononanoate, methyl 10-oxodecanoate | methyl heptanoate, methyl octanoate, methyl nonanoate, methyl 9-oxononanoate, methyl 10-oxodecanoate |
| alcohol | 1-heptanol | 1-pentanol, 1-octen-3-ol | - |
| hydrocarbons | heptane, octane | pentane | ethane, pentane |

Note: 2-undecenal is NOT in Frankel's oleate list as re-typed here (Cao 2020 and Chen 2017 both find it
as a major oleate product); 1-hexanol appears nowhere in Table 2.

### Numbers in the text (all cited)

- Odour thresholds (Frankel 1985): hydrocarbons 90-2150 ppm; alkanals 0.04-1; 2-alkenals 0.04-2.5;
  trans,trans-2,4-alkadienals 0.04-0.3 ppm.
- Steenson et al. 2002 (SPME, PV 5): soybean oil hexanal 23.5 % and 2-decenal 34.3 %; corn oil
  2-heptenal 29.5 % and trans-2-octenal 18.1 % (shares of volatiles).
- Frankel 1985 ranking of oxidised-flavour significance: trans,cis-2,4-decadienal > trans,trans-2,4-
  decadienal > trans,cis-2,4-heptadienal > 1-octen-3-ol > butanal > hexanal.
- Photo-oxidised vs autoxidised oleate (Frankel 1985): 1O2 gives more 2-decenal and octane; autoxidation
  gives more octanal and 10-oxodecanoate. 2-Heptenal and 2-butenal "noticeable" in 1O2-oxidised
  linoleate / linolenate, negligible in autoxidised.
- Relative autoxidation rate oleic:linoleic:linolenic 1 : 40-50 : 100 (Min & Bradley 1992); relative
  1O2 rates 1.0 : 1.4 : 1.9 (Vever-Bizet 1989; absolute 5.3, 7.3, 10.0 x 10^4 M-1 s-1, stearic 1.2 x
  10^4); linoleate + 1O2 1450x faster than + 3O2 (Rawls & Van Santen 1970). Temperature has little
  effect on 1O2 oxidation (Ea 0-6 kcal/mol; Yang & Min 1994).
- 1O2 energy 93.6 kJ above 3O2; lifetimes 2 / 17 / 700 µs in water / hexane / CCl4.

### Figure 3 — linoleate autoxidation to the 9- and 13-hydroperoxides (drawn; described)

Linoleic acid (C9=C10, C12=C13 both cis) -> - H• at C11 -> pentadienyl radical delocalised C9-C13 ->
radical localised at C13 (with C9=C10 cis, C11=C12 trans) or at C9 (C10=C11 trans, C12=C13 cis) -> + 3O2,
H• -> 13-hydroperoxide (9Z,11E) or 9-hydroperoxide (10E,12Z). The figure draws only these two.

### Figure 4 — "Mechanisms of hydroperoxide decomposition to form secondary oxidation products" (drawn; described)

Generic R2-CH=CH-CH2-R1 -> (O2, H•) -> R2-CH=CH-CH(OOH)-R1 -> - •OH -> alkoxyl R2-CH=CH-CH(O•)-R1. Two
scissions are drawn on the alkoxyl carbon: **A** breaks the bond to R1 (the saturated side) giving
R2-CH=CH-CHO (a 2-alkenal) + •R1; **B** breaks the bond between the alkoxyl carbon and the vinyl carbon
giving R2-CH=CH• (vinyl radical) + OHC-R1 (an alkanal). Fates: A1: •R1 + •OH -> R1-OH (alcohol); A2: •R1 +
R3H -> R1-H (alkane) + •R3. B1: R2-CH=CH• + •OH -> R2-CH=CH-OH (enol) <-> R2-CH2-CHO (alkanal); B2:
R2-CH=CH• + R3H -> R2-CH=CH2 (1-alkene) + •R3. Text: "The alkoxy radical then undergoes homolytic
beta-scission of the carbon-carbon bond and produces oxo-compounds and saturated or unsaturated alkyl
radicals (Figure 4). After electron rearrangement, the addition of hydroxyl radical, or hydrogen
transfer, the ultimate secondary lipid oxidation products are mostly low-molecular-weight aldehydes,
ketones, alcohols, and short-chain hydrocarbons, as shown in Table 2."

This is the same geometry as Cao 2020's "A-scission" (away from the C=C, 2-alkenal + alkyl radical) and
"B-scission" (toward the C=C, alkanal + vinyl radical -> second alkanal via the enol). Applied to
oleate: 8-OOH A -> 2-undecenal, B -> decanal (via vinyl radical) + methyl 8-oxooctanoate; 9-OOH A ->
2-decenal + methyl octanoate/heptane-type fragments, B -> nonanal + methyl 9-oxononanoate; 10-OOH A ->
methyl 10-oxo-8-decenoate + octyl radical (-> octane, 1-octanol), B -> nonanal; 11-OOH A -> methyl
11-oxo-9-undecenoate + heptyl radical (-> heptane, 1-heptanol), B -> octanal + methyl 10-oxodecanoate.
Every entry of Table 2's oleate column is generated by this reading; the review itself does not print
the oleate mapping sentence by sentence.

### Figure 7 — ene reaction of oleic acid with 1O2 (drawn; described)

1O2 adds to either end of C9=C10 with allylic H transfer: attack at C9 with H from C11 gives the
9-hydroperoxide with C10=C11 (trans); attack at C10 with H from C8 gives the 10-hydroperoxide with
C8=C9. Only these two (Table 3: 48 / 52 %).

### Figure 8 — ene reaction of linoleic acid with 1O2 (drawn; described)

Four products: 9-OOH (10E,12Z conjugated), 10-OOH (8E,12Z nonconjugated), 12-OOH (9Z,13E nonconjugated),
13-OOH (9Z,11E conjugated), each by a six-membered ene transition state drawn on one of the two double
bonds. Text: "When hydroperoxide is formed, double bond migration and trans fatty acid occur, producing
both conjugated and nonconjugated hydroperoxides".

### Figure 9 — "Formation of 2-pentylfuran from linoleic acid by singlet oxygen oxidation" (drawn; described; Min et al. 2003)

Linoleic acid CH3(CH2)4-CH=CH-CH2-CH=CH-(CH2)6-COOH -> 1O2 -> **10-hydroperoxide**
CH3(CH2)4-CH=CH-CH2-CH(OOH)-CH=CH-(CH2)6-COOH (nonconjugated, C8=C9 and C12=C13) -> cleavage of the
C10-C11 ... as drawn the chain breaks to give **3-nonenal** CH3(CH2)4-CH=CH-CH2-CHO -> 1O2 adds across the
C3=C4 double bond of 3-nonenal to a **1,2-dioxetane** -> ring-opened radical CH3(CH2)4-CH(OO•)-CH•-CH2-CHO
-> + 2 RH -> **4-hydroperoxynonanal** CH3(CH2)4-CH(OOH)-CH2-CH2-CHO (+ 2 R•) -> alkoxyl
CH3(CH2)4-CH(O•)-CH2-CH2-CHO -> **4-oxononanal** CH3(CH2)4-C(=O)-CH2-CH2-CHO -> dienol
CH3(CH2)4-C(OH)=CH-CH=CH-OH -> - H2O -> **2-pentylfuran**. Figure 10 is the analogue from linolenic acid
via its 1O2 10-hydroperoxide -> 3,6-nonadienal -> 2-(2-pentenyl)furan.

Text: "2-Pentylfuran and pentenylfuran were reported to be responsible for the beany flavor (Chang and
others 1966, 1983; Smouse and Chang 1967; Ho and others 1978; Smagula and others 1979), and 1O2 was
involved in their formation from linoleic and linolenic acids present in soybean oil, as shown in Figure
9 and 10 (Min and others 2003)."

## 4. Routes and numbers the repository can use

| route | reactant -> product | mechanism as drawn (figure) | numbers (all cited, none measured here) | evidence class |
|---|---|---|---|---|
| OL-ISO-AUTOX | oleate -> 8-/9-/10-/11-OOH | text + Table 1 | 26-28 / 22-25 / 22-24 / 26-28 % (Frankel 1985) | measured_ratio (cited; primary not on disk) |
| LA-ISO-AUTOX | linoleate -> 9-/13-OOH | Figure 3 | 48-53 / 48-53 % (Frankel 1985) | measured_ratio (cited) |
| LN-ISO-AUTOX | linolenate -> 9-/12-/13-/16-OOH | text | 28-35 / 8-13 / 10-13 / 28-35 % | measured_ratio (cited) |
| OL-ISO-1O2 | oleate + 1O2 -> 9-/10-OOH | Figure 7 (ene) | 48 / 52 % | measured_ratio (cited) |
| LA-ISO-1O2 | linoleate + 1O2 -> 9-/10-/12-/13-OOH | Figure 8 (ene) | 32 / 17 / 17 / 34 % | measured_ratio (cited) |
| SCISSION-A/B | allylic alkoxyl -> 2-alkenal + alkyl radical (A) or alkanal + vinyl radical (B); radicals -> alcohol / alkane / enol-alkanal / 1-alkene | Figure 4 | none | mechanism_drawn (generic) |
| OL-PRODUCTS | methyl oleate -> octanal, nonanal, 2-decenal, decanal; methyl 8-oxooctanoate, 9-oxononanoate, 10-oxodecanoate, 10-oxo-8-decenoate, 11-oxo-9-undecenoate; 1-heptanol; heptane, octane | Table 2 (Frankel 1985) | presence only | product slate (cited) |
| LA-PRODUCTS | methyl linoleate -> pentanal, hexanal, 2-octenal, 2-nonenal, 2,4-decadienal; methyl 8-oxooctanoate, 9-oxononanoate, 10-oxodecanoate; 1-pentanol, 1-octen-3-ol; pentane | Table 2 | presence only | product slate (cited) |
| LA-PF-1O2 | **linoleate 10-OOH (1O2) -> 3-nonenal -> 4-hydroperoxynonanal -> 4-oxononanal -> 2-pentylfuran** | Figure 9 (Min et al. 2003) | none | mechanism_drawn (cited primary: Min, Callison & Lee 2003) |
| VOLATILE-SHARES | soybean oil PV 5: hexanal 23.5 %, 2-decenal 34.3 %; corn oil: 2-heptenal 29.5 %, 2-octenal 18.1 % | - | Steenson et al. 2002 | measured_ratio (cited) |
| RATES | oleate:linoleate:linolenate autoxidation 1:40-50:100; 1O2 1.0:1.4:1.9; 1O2/3O2 for linoleate 1450 | - | Min & Bradley 1992; Vever-Bizet 1989; Rawls & Van Santen 1970 | rate ratios (cited); NOT for the hypothesis layer (rules never carry rates) |
| BDE | C11-H 50, C8-H/C14-H 75, C17/C18-H ~100 kcal/mol (linoleate); O-O vs O-H cleavage Ea difference 46 kcal/mol | - | Min & Boff 2002; Hiatt 1968 | cited thermochemistry; not computed here, but not a product measurement either — use only as ordering, not as a number |

## 5. Rule sketches (repository suggestions, not the review's)

Keys: `hexanal`, `nonanal`, `2_pentylfuran`, `1_octen_3_ol`, `e_2_octenal`, `heptanal`, `acrolein`
exist; **octanal, decanal, 2-decenal, 2,4-decadienal, pentane (only `PENTANE` in structures.yml),
1-pentanol, 1-heptanol** do not. Structures: `LOOH_9_ct`/`LOOH_13_ct` (methyl esters) exist; oleate
isomers are lumped (`LOOH_OL`); no 10-/12-HpODE.

**S1. The generic A/B scission pair (Figure 4) is what R18a/R18b already encode for the conjugated
diene.** For the mono-ene (oleate) the same two SMIRKS shapes apply with a single C=C; see cao2020 §5
S1/S2 for the four positive controls. Cite this review's Figure 4 as the generic form and Frankel 1985 as
the primary; cite Cao 2020 Fig. 6 for the oleate-specific drawing.

**S2. Singlet-oxygen ene hydroperoxidation (Figures 7, 8) — a formation rule, if the hypothesis layer
ever needs photo-oxidation.** Pattern: `[CH2:1][CH1:2]=[CH1:3]` + 1O2 -> `[CH1:1]=[CH1:2][CH1:3](OO)`
(allylic H moves, C=C shifts, OOH on the former vinyl carbon). Positive: methyl oleate -> methyl
9-hydroperoxy-10E-octadecenoate and methyl 10-hydroperoxy-8E-octadecenoate; linoleate -> the four of
Table 3. Negative: methyl stearate. Autoxidation must NOT produce the 10-/12-linoleate isomers ("not
observed in the autoxidation"), so a radical-formation rule should be a separate rule with only 9-/13-.

**S3. 2-Pentylfuran via 4-oxononanal (Figure 9; net).** Reactant 3-nonenal (from 10-HpODE); products
2-pentylfuran + H2O (net, after O2 and two H). Steps that could be written separately: (a) 3-nonenal +
O2 -> 4-hydroperoxynonanal; (b) 4-hydroperoxynonanal -> 4-oxononanal (+ H2O); (c) 4-oxononanal ->
2-pentylfuran + H2O (Paal-Knorr-type cyclodehydration of a 1,4-dicarbonyl).
- positive (c): `CCCCCC(=O)CCC=O` (4-oxononanal) -> `CCCCCc1ccco1` (2-pentylfuran)
- positive (a-c net): `CCCCC/C=C\CC=O` ((Z)-3-nonenal) -> `CCCCCc1ccco1`
- negative: nonanal `CCCCCCCCC=O` (no C3=C4); 2-nonenal `CCCCCC/C=C/C=O` (conjugated; would give a 4-hydroperoxy-2-alkenal / HNE, not a furan by this rule — Wanjala 2021 covers that branch); hexanal.
- source chain: Choe & Min 2006 Fig 9 <- Min, Callison & Lee 2003 (J Food Sci 68:1175). Alternative route
  from the same 10-HpODE precursor in Miyazaki 2023 (S4 there) goes to the C8 slate, and Miyazaki's
  2-pentylfuran comes from 13-HpODE by a furyl-hydroperoxide — the two papers disagree on the parent
  isomer of 2-pentylfuran (10- vs 13-HpODE); record both, do not merge.

**S4. Radical-fragment fates (Figure 4 A1/A2/B1/B2) as net rules:** alkyl radical -> alcohol (R-OH) or
alkane (R-H); vinyl radical -> alkanal (via enol) or 1-alkene. These explain 1-heptanol / heptane
(oleate 11-OOH A), 1-octanol / octane (10-OOH A), 1-pentanol / pentane (linoleate 13-OOH). None of them
is aldehyde -> alcohol; the review has no hexanal -> 1-hexanol step.

## 6. Flags

1. **Every number is second-hand** (Frankel 1985 chiefly, a book chapter not on disk). The evidence
   class for Tables 1-3 is "measured_ratio as cited"; the repo should not record them as verified until
   Frankel 1985 (or Frankel's primary papers behind it: Frankel, Neff et al. 1977-1984) is on disk.
2. **Figure 4's B branch produces a vinyl radical**, the step Miyazaki 2023 rejects for HpODE as
   "not chemically facile". The review does not discuss radical stability at all; it treats A and B as
   symmetric options.
3. **2-Undecenal is absent from Table 2's oleate column** although it is a major oleate product in Cao
   2020 and Chen 2017 (and follows from 8-OOH A-scission). Footnote (b) on octanal is not printed in the
   text layer.
4. **Figure 9 assigns 2-pentylfuran to the singlet-oxygen 10-hydroperoxide**; the review offers no
   autoxidative (9-/13-OOH) route to 2-pentylfuran. Miyazaki 2023 (13-HpODE furyl route, thermal) and
   Wanjala 2021 (HNE cyclisation, lysine) are the other two accounts on disk; Yang 2024 assigns it to
   9-HPOD. Four papers, four parent assignments.
5. **The text says autoxidation hydroperoxides "are conjugated dienes"**: this holds for linoleate /
   linolenate, not for oleate (mono-ene, no diene) — the sentence is loose.
6. **Rate ratios and bond energies** are recorded for orientation only; the hypothesis layer does not
   carry rates, and the BDE values are textbook thermochemistry, not measurements from these oils.
7. No DFT in the review; nothing to mark inadmissible.
