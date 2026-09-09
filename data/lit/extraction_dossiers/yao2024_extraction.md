# Yao 2024 — EXTRACTION (neat methyl linoleate, 150 mg in a capped 25 mL vial, 180 C, 0.5-2 h; non-volatiles at 1 h by Ph3P / NaBH4 reduction + SPE + TMS GC-MS; volatiles by HS-SPME-GC-MS; B3LYP/6-311G(d,p) mechanisms)
### The high-temperature methyl-ester companion to Miyazaki 2023: it detects and puts a number on the 9-hydroperoxy-dihydrofuran intermediate of the furyl route (R31) and on the epoxy-hydroperoxides at 180 C, and argues (by DFT, inadmissible) that methyl 9-oxononanoate does not come from the direct scission of the 9-hydroperoxide.

**Source on disk:** `data/articles/Yao2024.pdf` (11 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Yao2024.txt`); Tables 2 and 3 came through clean and are re-typed
below. Table 1 is a table of drawn structures with no text layer and was read from a 110 dpi raster
of p. 705, as were the schemes in Figures 1, 7 and 8 (pp. 707, 712, 713) and the axis labels of
Figures 3-4 (p. 709). No value was read off any bar or curve. The Supporting Information (species
library and coordinates for the DFT calculations only) is not on disk and would not be admissible.

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanisms of the Formation of Nonvolatile and Volatile Oxidation Products from Methyl Linoleic Acid at High Temperatures" |
| Authors | Yunping Yao, Tianliang Wang, Zhiyuan Qiang, Wenqi Du, Changmo Li* (Tianjin University of Science and Technology) |
| Venue | J. Agric. Food Chem. 2024, 72, 704-714. Received June 29 2023, accepted December 11 2023, published December 22 2023 |
| DOI | 10.1021/acs.jafc.3c04405 |
| Naming | Me-LA = methyl linoleate (the paper says "methyl linoleic acid"). 9-/13-/10-/8-hydroperoxides are the methyl esters (Table 1 structures: 9-OOH 10E,12Z; 13-OOH 9Z,11E; 8-OOH 9Z,12Z; 10-OOH 8E,12Z). "Alkyl radical A" = •CH2(CH2)6COOCH3 (the C1-C8 ester radical); "pentane radical" = n-pentyl; "alkoxy radical A" = the 9-alkoxyl of the 10E,12Z-diene; "alkoxy radical B" = the 11-alkoxyl of the 12,13-epoxy-9-ene (from epoxide B); "alkoxy radical C" = the 9-alkoxyl of 9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate (from dihydrofuran C). "Epoxides A" = methyl (Z)-8-(3-(1-hydroperoxyoct-2-en-1-yl)oxiran-2-yl)octanoate (9,10-epoxy-11-hydroperoxy-12Z); "Dihydrofurans C" = methyl 9-hydroperoxy-9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate. Epoxide B's name is not printed (footnote misprint; §5) — by the drawn mechanism it is the 12,13-epoxy-11-hydroperoxy-9-ene. "Peroxide-linked dimers A-D" = ROO-adducts of the 9- and 13-hydroperoxides (Figure 1b). "2-aldehyde-5-pentylfuran" = 5-pentyl-2-furaldehyde. |
| Companions | Miyazaki 2023 (their ref 6; `miyazaki2023_extraction.md`, the repo's R31/R32 anchor; free-acid HpODE isomers at 120 C); Berdeaux 2012 (ref 17; high-temperature Me-LA degradation products); Ding 2022 (ref 29; the group's methyl oleate DFT paper); Li 2019 (ref 28; trilinolein isomerisation DFT) |

## 1. Why it matters

The lipid rules R18a / R18b (schroen2022 anchor, Frankel 1989 slate) and R31 (miyazaki2023 anchor, the
furyl route to 2-pentylfuran) are written as net rules from hydroperoxide to volatiles. Miyazaki
heated purified free-acid hydroperoxides at 120 C for 5 min and proposed, without isolating them, a
2,5-dihydrofuran hydroperoxide (furyl route) and an epoxy-allyl hydroperoxide (hexanal route) as
intermediates. This paper heats the methyl ester itself, neat, at 180 C — the frying regime the rules'
condition strings claim ("25-200 C") — and

- **detects by GC-MS and puts a level on exactly those intermediates**: "dihydrofurans C" (Miyazaki's
  furyl-route hydroperoxide, 0.59 mmol/kg at 1 h) and the epoxy-hydroperoxides A + B (2.47 mmol/kg),
  alongside the four positional hydroperoxides (2.24 mmol/kg total, 9- > 13- > 10- > 8-) and
  "peroxide-linked dimers" (0.83 mmol/kg);
- prints **text ranges (mmol/kg, 0.5-2 h) for eleven volatiles**, including 2,4-decadienal, methyl
  9-oxononanoate, methyl octanoate, methyl 13-oxo-9,11-tridecadienoate and pentanal — the R18a/R18b
  product pairs for the 9- and 13-hydroperoxides — and lists hexanal, 2-heptenal, 2-octenal,
  2-pentylfuran, 2-octylfuran, 1-octen-3-ol, 4,5-epoxy-2-decenal and methyl 10-oxo-8-decenoate as
  detected (levels FIGURE-ONLY);
- draws (Figures 1, 6, 7, 8) a mechanism in which the two "side A" scissions of R18a are ordinary
  beta-scissions, while methyl 9-oxononanoate (R18b's product from the 9-hydroperoxide) is assigned
  mainly to the alkoxyls of the dimers, the epoxide and the dihydrofuran rather than to the
  9-hydroperoxide's own alkoxyl. The argument is DFT and inadmissible as numbers; the DETECTION of the
  intermediates and the measured product levels are not.

Nothing here is a rate. There is one temperature (180 C), one non-volatile time point (1 h), a closed
oxygen-limited vial, and single-internal-standard SPME quantification.

## 2. Methods as they matter to a model

- **Substrate and heating.** Me-LA 99 % (Macklin), **150 mg** neat, in a **25 mL capped glass vial
  with PTFE septum**, silicone oil bath **180 C**, **0.5, 1, 1.5, 2 h**, cooled to room temperature.
  No stirring, no added water, metals or antioxidant. The vial is closed: 25 mL of air holds ~0.21
  mmol O2 (my arithmetic, 21 % of 25 mL at 25 C) against 0.51 mmol Me-LA (MW 294.47), so the system
  is **oxygen-limited** and the headspace also holds the volatiles.
- **Non-volatiles (1 h sample only).** (i) Hydroperoxides and "heterocyclic compounds": 100 µL
  oxidised Me-LA + 100 µL 0.25 mol/L Ph3P in benzene, 10 min (hydroperoxides -> alcohols); silica
  SPE (500 mg): hexane/ether 98:2 (non-polar), 90:10 ("heterocyclic compounds"), then 24 mL ether
  ("hydroxy compounds"); dried; IS 5 µL methyl pentadecanoate 0.862 mg/mL (= 4.31 µg); 500 µL
  derivatised with 100 µL BSTFA 75 C 20 min; GC-MS DB-5MS 30 m x 0.25 mm x 0.25 µm, 70 C (1 min) ->
  170 C at 8 C/min (6 min) -> 280 C at 8 C/min (30 min), EI 70 eV, scan 10-500. (ii) Peroxide-linked
  dimers: 100 µL oxidised Me-LA + 5 mL 2 % NaBH4 in methanol / 0.05 M NaOMe, 30 min, acetic acid,
  dried, hexane, centrifuged, same SPE and TMS derivatisation, same IS and GC-MS. **Consequence:**
  every "hydroperoxide" is measured as the corresponding TMS-hydroxy ester; pre-existing hydroxy
  esters (e.g. from the radical-substitution route the paper itself proposes) are not distinguished
  from reduced hydroperoxides. Identification by m/z of the TMS ethers (Figure 2: 8-OH m/z 239, 10-OH
  m/z 271, dimers m/z 173 and 259, epoxides A + B m/z 199 and 285 "from the same peak", dihydrofuran
  C m/z 259, 4-hydroperoxy-2-nonenal m/z 157, "4-hydroperoxy-5-peroxide-2-decenal" m/z 173,
  5-pentyl-2-furaldehyde m/z 97). Calibration: single IS, response factors not stated (assumed 1).
- **Volatiles.** HS-SPME, DVB/CAR/PDMS 50/30 µm, "the heated samples containing 5 µL of
  2-methyl-3-heptanone (0.816 mg/mL ethanol)" (= **4.08 µg IS**; if added to the whole 150 mg charge
  that is 27 mg/kg = 0.212 mmol/kg), 70 C water bath, 40 min exposure; desorption 250 C 7 min
  splitless; DB-5MS; oven 30 C (5 min) -> 150 C at 4 C/min (1 min) -> 250 C at 8 C/min (6 min); EI
  70 eV, scan 10-500; NIST14 + RI vs C8-C40 alkanes. "According to the content of the internal
  standard, the amounts of oxidation compounds contained in the samples were calculated" —
  **single-IS semi-quantification, response factor 1 implied, reported as mmol/kg** (per kg of the
  Me-LA charge, presumably). Whether the SPME was run on the whole vial or an aliquot is not stated.
- **Replicates.** "All analyses were repeated in triplicate", mean ± SD (Origin 2018). Error bars are
  drawn on Figure 4.
- **DFT.** Gaussian 09, B3LYP/6-311G(d,p), transition states per the group's methyl oleate paper;
  isopropyl hydroperoxide used as the model hydroperoxide for the substitution step. Every barrier
  and every rate constant in the paper (Table 3, Figures 1, 6, 7 and the text) is computed, not
  measured — **INADMISSIBLE as numbers** in this repository; their existence is recorded in §3.

## 3. Tables re-typed

### Table 1. "Molecule Structures Related to the Oxidation Mechanisms" (drawn structures; described)

| name (paper) | structure (read from the drawing) |
|---|---|
| 9-hydroperoxide | methyl 9-hydroperoxy-10E,12Z-octadecadienoate (= repo `LOOH_9_ct`) |
| 13-hydroperoxide | methyl 13-hydroperoxy-9Z,11E-octadecadienoate (= repo `LOOH_13_ct`) |
| 8-hydroperoxide | methyl 8-hydroperoxy-9Z,12Z-octadecadienoate (non-conjugated; no repo key) |
| 10-hydroperoxide | methyl 10-hydroperoxy-8E,12Z-octadecadienoate (non-conjugated; Miyazaki's 10-HpODE as the ester; no repo key) |
| alkyl radical A | H3COOC(CH2)6CH2• (methyl 8-octanoate radical, C1-C8) |
| pentane radical | H3C(CH2)3CH2• |
| alkoxy radical A | 9-oxyl of the 10E,12Z-diene ester |
| alkoxy radical B | 11-oxyl of methyl 12,13-epoxy-9-octadecenoate (the alkoxyl of epoxide B) |
| alkoxy radical C | 9-oxyl of methyl 9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate (the alkoxyl of dihydrofuran C) |

### Table 2. "Content of Nonvolatile Oxidation Products of Methyl Linoleic Acid at 180 C for 1 h"

Unit mmol/kg, mean ± SD, n = 3. Footnote a gives the name of epoxides A; footnotes b and c BOTH give
the name of dihydrofurans C (epoxide B's name is missing; misprint).

| class | compound | content (mmol/kg) |
|---|---|---:|
| hydroperoxides | 9-hydroperoxides | 1.06 ± 0.06 |
| | 13-hydroperoxides | 0.66 ± 0.07 |
| | 10-hydroperoxides | 0.33 ± 0.01 |
| | 8-hydroperoxides | 0.19 ± 0.00 |
| heterocyclic compounds | epoxides A and B (one GC peak) | 2.47 ± 0.31 |
| | dihydrofurans C | 0.59 ± 0.00 |
| peroxide-linked dimers | peroxide-linked dimers | 0.83 ± 0.04 |

Arithmetic (mine): hydroperoxides sum to 2.24 (Figure 8 prints 2.14 — misprint); heterocyclics 3.06;
2.24 : 0.83 : 3.06 = **2.70 : 1 : 3.69**, the abstract's ratio, reproduces; 1.06 : 0.66 : 0.33 : 0.19 =
5.58 : 3.47 : 1.74 : 1 reproduces; 9-/13- = 1.61 reproduces. The text once quotes 13-OOH as 0.67.
Total identified non-volatile oxygenates 6.13 mmol/kg = 0.18 mol % of the Me-LA charge (3396 mmol/kg).

### Table 3. "Reaction Rate Constants Related to the Substitution Reaction of Peroxide Radicals with LA" — DFT, INADMISSIBLE

| pathway | 60 C | 180 C |
|---|---:|---:|
| A (H-abstraction at the bis-allylic C11) | 36.50 s-1 M-1 | 4.82e4 s-1 M-1 |
| B (H-abstraction at C8 / C14, mono-allylic) | 2.24 s-1 M-1 | 7.00e3 s-1 M-1 |

Recorded for existence only; evidence class dft_inadmissible. (The text twice cites "Table 3" for
structures that are in Table 1.)

### Volatile levels printed in the text (mmol/kg; ranges over the 0.5-2 h series of Figure 4)

Class shares (basis not stated — presumably of the summed mmol/kg or peak area): aldehydes 40.77 %,
"alkanes" 19.89 % (the only non-oxygenated products named are the esters methyl octanoate and methyl
heptanoate), alcohols 9.02 %, furans 6.11 %, epoxides 0.46 %, acids 2.50 % (sum 78.75 %; the remainder
is unassigned). Aldehydes with C > 9 = 64.07 % of aldehydes. Alcohols listed as pentanol 1.62 %,
1-octen-3-ol 1.68 %, methyl 8-hydroxyoctanoate 9.02 % (these do not sum to the class's 9.02 %).

| Fig. 4 no. | compound | printed range (mmol/kg) | printed elsewhere |
|---:|---|---|---|
| 15 | methyl 9-oxononanoate | 52.45-68.39 | "large amounts" |
| 16 | 9-methoxy-9-oxononanoic acid (monomethyl azelate; the acid of no. 15) | 9.99-11.42 | |
| 10, 11 | 2,4-decadienal (two isomers, one entry) | 28.65-33.50 | |
| 20 | methyl 13-oxo-9,11-tridecadienoate | 1.68-2.45 | |
| 19 | methyl 11-oxo-9-undecenoate | 10.91-12.34 | |
| 12 | methyl 8-oxooctanoate | 17.22-23.00 | |
| 1 | pentanal | 2.07-2.86 | |
| 14 | methyl 8-hydroxyoctanoate | "29.52-24.99" (printed in that order, i.e. falling) | "34.99 mmol/kg" in §3.3 — inconsistent |
| 2 | pentanol | 8.68-9.70 | |
| 9 | methyl octanoate | 94.80-117.54 | "the most readily formed volatile product" |
| 7 | methyl heptanoate | 5.56-6.14 | |
| 3 | hexanal | not printed | FIGURE-ONLY |
| 4 | 2-heptenal | not printed | FIGURE-ONLY |
| 8 | 2-octenal | not printed | FIGURE-ONLY |
| 5 | 1-octen-3-ol | not printed (1.68 % share) | FIGURE-ONLY |
| 6 | 2-pentylfuran | not printed | FIGURE-ONLY |
| 18 | 2-n-octylfuran | not printed | FIGURE-ONLY |
| 13 | 4,5-epoxy-2-decenal | not printed (epoxides 0.46 %) | FIGURE-ONLY |
| 17 | methyl 10-oxo-8-decenoate | not printed ("exceeded that of 1-octen-3-ol") | FIGURE-ONLY |

Also printed: "The Me-LA concentration was 443 mmol/kg at 180 C for 1 h" (see §5, flag 3); "the
heating process only altered the contents of each volatile oxidation product without affecting the
quantitative relationship among the three volatile oxidation products [methyl octanoate, methyl
8-hydroxyoctanoate, 2,4-decadienal]".

Figure 3 = total ion chromatogram at 180 C (peaks 1-20 as numbered above; methyl octanoate is by far
the tallest peak). Figure 4 = bar chart, "Concentration (mmol/kg)" vs compound, four bars (30, 60,
90, 120 min) per compound, error bars drawn. Figure 2 = mass spectra (identification evidence).
Figure 5 = ESP / HOMO of isopropyl hydroperoxide (DFT).

### Schemes described in words (Figures 1, 6, 7, 8); numbers on the arrows are DFT and inadmissible

**Figure 1a** (free-energy profile, DFT): ROO• + Me-LA -> ROOH + pentadienyl radical via C11-H
(pathway A, TS 17.17 kcal/mol, product -11.94) or via a mono-allylic C-H (pathway B, TS 21.26,
product -1.01). Pathway B is the drawn origin of the 8- and 10-hydroperoxides (O2 at C8 or C10 of
the non-conjugated allyl radical; "an energy gap of 1.13 kcal/mol between the 10-peroxyl radicals
and the 8-peroxyl radicals").

**Figure 1b** (peroxide-linked dimers, four equations): a peroxyl radical ROO• ADDS across the
conjugated diene of a hydroperoxide; O2 then adds to the resulting allyl radical and LH donates H.
(1) 9-OOH + ROO• at C13 -> dimer A (barrier 11.44 kcal/mol, 3.02e7 s-1 M-1); (2) 9-OOH + ROO• at
C10 -> dimer B (13.05 kcal/mol, 5.09e6); (3) and (4) 13-OOH + ROO• at C9 or C12 -> dimers C and D.
Each "dimer" carries the original OOH, a new OOR bridge and a new OOH. Evidence: TMS ions m/z 173 and
259 after NaBH4 reduction (see §5, flag 6).

**Figure 1c** (epoxides): 9-OOH -> -HO• -> 9-alkoxyl -> cyclisation onto C10 (1.26 kcal/mol, 6.76e12
s-1) -> 9,10-epoxide with the allyl radical at C11-C13 -> O2 at C11, LH -> **epoxide A**
(9,10-epoxy-11-OOH-12Z). Same from 13-OOH: 13-alkoxyl onto C12 -> 12,13-epoxide, radical C9-C11, O2
at C11 -> **epoxide B** (12,13-epoxy-11-OOH-9-ene). Both elute as one GC peak (m/z 199, 285); the
authors read A > B from the ion abundances, "consistent with" 9-OOH > 13-OOH.

**Figure 1d** (dihydrofuran C): 9-OOH + ROO• at C13 (as in dimer A) -> the C13-OOR adduct cyclises:
the inner peroxide oxygen attacks C10 with expulsion of RO• (17.15 kcal/mol, 5.06e4 s-1) ->
**dihydrofuran C** = 2,5-dihydrofuran ring O-bridging C10 and C13, C11=C12, pentyl on C13, and the
untouched 9-OOH on the C1-C9 ester arm. This is the SAME molecule (as the methyl ester) as the
furyl-route hydroperoxide Miyazaki draws from the 13-alkoxyl cyclising onto C10 followed by O2 at
C9 (their 24 -> 25 -> 26). The two papers thus give the intermediate two different parents (Yao: the
9-OOH plus a peroxyl radical; Miyazaki: the 13-OOH alone). Evidence: m/z 259 (Figure 2e).

**Figure 6** (hydroxy compounds; DFT on isopropyl hydroperoxide): ROOH -> RO• + HO• (29.30 kcal/mol,
k 8.17e-2 s-1 at 180 C, "cannot be ignored"); the alternative, an alkyl radical abstracting OH from
ROOH (radical substitution, "the -OH of the hydroperoxide approached the alkyl radicals", 16.18
kcal/mol, 1.60e5 s-1 M-1) -> R'OH + RO•, called the predominant route to alcohols and a major source
of alkoxyls. Basis for methyl 8-hydroxyoctanoate and pentanol.

**Figure 7** (ten volatile-forming steps; R1 = butyl, R2 = (CH2)7COOCH3; barriers in kcal/mol and
DFT k):
1. alkoxy radical A (9-oxyl, 10E,12Z) -> C9-C10 scission -> methyl 9-oxononanoate + a C9 dienyl
   radical (24.46; 8.48e2 s-1) — the paper's "difficult" step (conjugation);
2. dimer-B-derived 9-alkoxyl (10-OOR, 13-OOH) -> methyl 9-oxononanoate + C10 radical (2.56; 5.52e11
   s-1) -> fast -> **4-hydroperoxy-2-nonenal** + RO• (detected, m/z 157). The same from dimer C's
   10-OOH position;
3. alkoxy radical C (dihydrofuran C's 9-oxyl) -> C9-C10 scission -> **methyl 9-oxononanoate +
   2-pentyl-2,5-dihydrofuran-5-yl radical** (0.86; 3.63e12 s-1) — the radical is the 2-pentylfuran
   skeleton one H-loss away; the paper does not name 2-pentylfuran here;
4. alkoxy radical A -> C8-C9 scission -> **2,4-decadienal + alkyl radical A** (6.36; 8.11e9 s-1) —
   R18a for LOOH_9;
5. 13-alkoxyl -> C13-C14 scission -> **methyl 13-oxo-9,11-tridecadienoate + pentyl radical** (6.90;
   4.45e9 s-1) — R18a for LOOH_13;
6. the C8 allyl radical from the 10-hydroperoxide (drawn "•(CH)3R1") + ROOH -> **1-octen-3-ol** +
   RO• (20.57; 1.11e3 s-1 M-1), the partner of methyl 10-oxo-8-decenoate — the R32 pair;
7. alkoxy radical B -> C11-C12 scission -> **methyl 11-oxo-9-undecenoate + a 2-pentyloxiranyl
   radical** (7.31; 2.81e9 s-1; the competing C10-C11 scission 18.71, 9.12e3). The figure labels the
   product "methyl 9-oxononanoate"; the text and the drawn fragment say methyl 11-oxo-9-undecenoate
   (label misprint);
8. alkoxy radical C -> C8-C9 scission -> **alkyl radical A + 5-pentyl-2,5-dihydrofuran-2-carbaldehyde**
   (9.11; 3.84e8 s-1) -> detected as **5-pentyl-2-furaldehyde** (m/z 97; "2-aldehyde-5-pentylfuran");
9. alkyl radical A + Me-LA -> **methyl octanoate** + pentadienyl radical (9.54; 2.38e8 s-1 M-1);
10. alkyl radical A + ROOH -> **methyl 8-hydroxyoctanoate** + RO• (16.18; 1.60e5 s-1 M-1).
Also in the text: dimers A and D (9-OOH retained) -> 9-alkoxyl -> C8-C9 scission -> alkyl radical A +
"4-hydroperoxy-5-peroxide-2-decenal" (8.11; 8.11e9 s-1; detected m/z 173); pentyl radical + ROOH ->
pentanol; pentanol -> pentanal and methyl 8-hydroxyoctanoate -> methyl 8-oxooctanoate by oxidation;
methyl 9-oxononanoate -> 9-methoxy-9-oxononanoic acid by aldehyde autoxidation. No mechanism is
drawn or stated for **hexanal, 2-heptenal, 2-octenal, 2-pentylfuran, 2-octylfuran or
4,5-epoxy-2-decenal**.

**Figure 8** (summary): linoleate -> hydroperoxides ("2.14 mmol/kg"; 9:13:10:8 = 5.58:3.47:1.74:1)
-> (i) polyunsaturated aldehydes + alkyl radicals; (ii) peroxide-linked dimers (0.83) -> methyl
9-oxononanoate; (iii) epoxides (2.47) -> monounsaturated aldehydes; (iv) dihydrofurans (0.59) ->
methyl 9-oxononanoate and -> alkyl radicals; alkyl radicals -> methyl esters (9.54) or hydroxy
compounds (16.18). The abstract's "secondary oxidation products mainly came from peroxides (77 %)"
is a DFT-based apportionment, not a measurement.

## 4. Numbers / routes the repository can use

Repo keys: `LOOH_9_ct`, `LOOH_13_ct`, `ME_9_OXONONANOATE`, `ME_OCTANOATE`, `DECADIENAL`,
`ME_13_OXO_TRIDECADIENOATE`, `PENTANE`, `HEXANAL` exist in `data/species/structures.yml`;
`2_pentylfuran`, `hexanal`, `1_octen_3_ol`, `e_2_octenal` in `data/keys/compounds.yml`. No keys for
the 8-/10-hydroperoxides, epoxides A/B, dihydrofuran C, methyl heptanoate, methyl 8-hydroxy-/8-oxo-
octanoate, methyl 10-oxo-8-decenoate, methyl 11-oxo-9-undecenoate, pentanal, pentanol, 2-heptenal,
5-pentylfurfural, 2-octylfuran, 4,5-epoxy-2-decenal (R31's control names `PENTYLFURAN`,
`ME_8_FURYL_OCTANOATE`, `LOOH_10` are used in `reaction_rules.yml` but are not in structures.yml).

| quantity or route | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| 9- / 13- / 10- / 8-hydroperoxide (as Ph3P-reduced TMS-hydroxy esters) | 1.06 ± 0.06 / 0.66 ± 0.07 / 0.33 ± 0.01 / 0.19 ± 0.00 | mmol/kg | neat Me-LA, closed 25 mL vial, 180 C, 1 h, n = 3, single-IS GC-MS | Table 2 | level_only |
| 9-OOH : 13-OOH : 10-OOH : 8-OOH | 5.58 : 3.47 : 1.74 : 1 (9/13 = 1.61) | — | same | text §3.1, Table 2 | within_study_ratio (same derivatisation and IS for all four) |
| epoxides A + B (9,10-epoxy-11-OOH-12-ene + 12,13-epoxy-11-OOH-9-ene, co-eluting) | 2.47 ± 0.31 | mmol/kg | same | Table 2 | level_only |
| dihydrofuran C = methyl 9-hydroperoxy-9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate (R31's intermediate) | 0.59 ± 0.00 | mmol/kg | same | Table 2; m/z 259, Fig. 2e | level_only; the detection itself supports R31 mechanism_drawn |
| peroxide-linked dimers | 0.83 ± 0.04 | mmol/kg | same, NaBH4 route | Table 2 | level_only (identification weak, flag 6) |
| hydroperoxides : dimers : heterocyclics | 2.70 : 1 : 3.69 | — | same | abstract, §3.1 | within_study_ratio (across two work-ups) |
| methyl octanoate / methyl 8-hydroxyoctanoate / methyl 8-oxooctanoate (the three fates of alkyl radical A) | 94.80-117.54 / 24.99-29.52 (or 34.99) / 17.22-23.00 | mmol/kg | 180 C, 0.5-2 h series | text §3.2-3.3 | level_only (SPME, single IS, RF 1) |
| 2,4-decadienal (R18a partner of alkyl radical A from LOOH_9) | 28.65-33.50 | mmol/kg | same | text §3.2 | level_only |
| (methyl octanoate + 8-hydroxy + 8-oxo) / 2,4-decadienal | ~4.1-5.9 | — | same | derived, mine | within_study_ratio (SPME responses differ between an ester and a dienal; order-of-magnitude only) |
| methyl 13-oxo-9,11-tridecadienoate / (pentanol + pentanal) (R18a pair from LOOH_13) | 1.68-2.45 vs 8.68-9.70 + 2.07-2.86 | mmol/kg | same | text §3.2 | level_only; the C5 side exceeds the C13 side ~4-7x (derived) |
| methyl 9-oxononanoate (R18b product from LOOH_9) and its acid | 52.45-68.39 and 9.99-11.42 | mmol/kg | same | text §3.2 | level_only |
| methyl 11-oxo-9-undecenoate (epoxide B scission product) | 10.91-12.34 | mmol/kg | same | text §3.2 | level_only |
| methyl heptanoate (from 8-OOH, side-A analogue) | 5.56-6.14 | mmol/kg | same | text §3.2 | level_only |
| methyl octanoate / methyl heptanoate | ~15-21 (vs 9-OOH / 8-OOH = 5.58) | — | same | derived, mine | within_study_ratio |
| hexanal, 2-heptenal, 2-octenal, 2-pentylfuran, 2-octylfuran, 1-octen-3-ol, 4,5-epoxy-2-decenal, methyl 10-oxo-8-decenoate | detected; levels not printed | mmol/kg (Fig. 4 axis) | same | Figs 3, 4 | figure_only |
| class shares: aldehydes 40.77, "alkanes" (esters) 19.89, alcohols 9.02, furans 6.11, epoxides 0.46, acids 2.50 % | as printed | % (basis unstated) | same | text §3.2 | level_only |
| 5-pentyl-2-furaldehyde, 4-hydroperoxy-2-nonenal, "4-hydroperoxy-5-peroxide-2-decenal" | detected (m/z 97, 157, 173) | — | 180 C | Fig. 2f-h | level_only (presence) |
| Me-LA remaining at 1 h | "443 mmol/kg" | mmol/kg | 180 C, 1 h | text §3.2 | level_only, DOUBTFUL (flag 3) |
| route: LOOH_9 -> 2,4-decadienal + C8 ester radical -> methyl octanoate (H from LH) / methyl 8-hydroxyoctanoate (OH from ROOH) -> methyl 8-oxooctanoate | — | — | drawn; all products measured | Fig. 7(4), (9), (10); Fig. 8 | mechanism_drawn (= R18a side A for LOOH_9, with the alkyl radical's three fates) |
| route: LOOH_13 -> methyl 13-oxo-9,11-tridecadienoate + pentyl radical -> pentanol -> pentanal | — | — | drawn; products measured; pentane not reported | Fig. 7(5); text | mechanism_drawn (= R18a side A for LOOH_13) |
| route: LOOH_9 (+ ROO•) -> dihydrofuran C -> 9-alkoxyl -> methyl 9-oxononanoate + 2-pentyl-dihydrofuranyl radical (-> 2-pentylfuran, my completion) OR -> alkyl radical A + 5-pentylfurfural | — | — | intermediate detected and quantified; both scission partners detected | Fig. 1d, 7(3), 7(8); Fig. 2e, 2h | mechanism_drawn (= R31's product pair via the same intermediate Miyazaki draws from LOOH_13) |
| route: LOOH_9 / LOOH_13 -> alkoxyl cyclises onto the adjacent diene carbon -> epoxy-allyl radical -> O2 -> epoxides A / B; epoxide B -> methyl 11-oxo-9-undecenoate + C7 oxiranyl radical | — | — | epoxides detected (one peak); the undecenoate measured | Fig. 1c, 7(7) | mechanism_drawn (Miyazaki's epoxy-allyl route to hexanal is the same first step) |
| route: 10-OOH -> methyl 10-oxo-8-decenoate + 2-octenyl radical -> 1-octen-3-ol (OH from ROOH) | — | — | both products detected | Fig. 7(6) | mechanism_drawn (= R32) |
| route: methyl 9-oxononanoate mainly from dimer / epoxide / dihydrofuran alkoxyls, not from the 9-alkoxyl's C9-C10 scission | — | — | argued from DFT barriers (24.46 vs 0.86-2.56 kcal/mol) and from the "large amounts" of the aldehyde | §3.4 | mechanism_drawn, ground = dft_inadmissible; net product unchanged |
| all DFT barriers and rate constants (Table 3; Figs 1, 6, 7; text: 21.26, 17.17, 11.44, 13.05, 1.26, 17.15, 29.30 / 8.17e-2 s-1, 16.18 / 1.60e5, 24.46 / 8.48e2, 2.56 / 5.52e11, 0.86 / 3.63e12, 6.36 / 8.11e9, 6.90 / 4.45e9, 20.57 / 1.11e3, 18.71 / 9.12e3, 7.31 / 2.81e9, 8.11 / 8.11e9, 9.11 / 3.84e8, 9.54 / 2.38e8, 1.13 kcal/mol; "77 % from peroxides") | exist | kcal/mol; s-1; s-1 M-1 | B3LYP/6-311G(d,p), gas phase, isopropyl hydroperoxide model for the substitution step | throughout | **dft_inadmissible** |

## 5. Flags

1. **Single temperature (180 C), single non-volatile time point (1 h), no rates.** The 0.5-2 h volatile
   series exists only as bars (Figure 4); the text prints ranges without saying which end is which
   time (methyl 8-hydroxyoctanoate is printed "29.52-24.99", so at least that one falls).
2. **Closed, oxygen-limited vial.** 150 mg Me-LA (0.51 mmol) under ~0.21 mmol O2 (my arithmetic).
   Whether the vial was re-opened between time points is not stated (the four times may be four
   vials). Product ratios at 180 C in air-saturated frying oil may differ; the paper itself argues
   that low dissolved O2 favours the dihydrofuran over the dimers.
3. **"Me-LA concentration was 443 mmol/kg at 180 C for 1 h"** — neat Me-LA is 3396 mmol/kg, so this
   would be 87 % consumption in an oxygen-limited vial whose identified products sum to ~0.3 mol/kg.
   Unexplained; do not use.
4. **SPME quantification is semi-quantitative**: one internal standard (2-methyl-3-heptanone, 4.08
   µg), response factors not stated, headspace partitioning at 70 C not corrected; "mmol/kg" values
   for esters vs dienals vs alcohols are not on a common footing. The volatile numbers (~300 mmol/kg
   summed) exceed the non-volatile primary products (6.1 mmol/kg) by 50x, which is plausible for a
   transient hydroperoxide pool at 180 C but also what an over-responding SPME analyte set would
   produce. Treat every volatile level as order-of-magnitude and every cross-class ratio as weaker
   than that.
5. **Hydroperoxides are measured after Ph3P reduction as TMS-hydroxy esters**, so hydroxy esters
   already present (the paper's own radical-substitution products) are counted as hydroperoxides;
   the 9- > 13- ordering could partly be a 9-OH vs 13-OH ordering.
6. **Dimer identification is weak** (mine): the evidence is TMS ions m/z 173 and 259, which are the
   standard alpha-cleavage fragments of 13-O-TMS and 9-O-TMS C18 esters and would also arise from
   monomeric reduced hydroperoxides; a C36 peroxide-linked dimer is unlikely to survive NaBH4 /
   NaOMe and elute intact on a 280 C DB-5MS run. The 0.83 mmol/kg is a level for "whatever gave
   those ions in that SPE fraction".
7. **Epoxides A and B co-elute** (one peak, 2.47 mmol/kg combined); A > B is an ion-ratio inference.
   Epoxide B's name is missing from Table 2 (footnotes b and c both name dihydrofuran C).
8. **Misprints**: Figure 8 "2.14 mmol/kg" vs Table 2's 2.24; 13-OOH 0.66 (table) vs 0.67 (text); methyl
   8-hydroxyoctanoate 24.99-29.52 vs 34.99; Figure 7(7) product labelled methyl 9-oxononanoate but
   drawn and discussed as methyl 11-oxo-9-undecenoate; "Table 3" cited for Table 1 structures;
   alcohol sub-shares (1.62 + 1.68 + 9.02) exceed the class share (9.02); class shares sum to 78.75 %.
9. **No mechanism for hexanal, 2-heptenal, 2-octenal, 2-pentylfuran, 2-octylfuran or
   4,5-epoxy-2-decenal**, and no printed level for any of them, although hexanal and 2-heptenal are
   visibly among the larger bars in Figure 4 (not read). For R18b's hexanal and R31's 2-pentylfuran
   this paper contributes the intermediate (dihydrofuran C) and the co-products, not the volatile.
10. **Pentane not reported** (as in Miyazaki flag 8): the pentyl radical is accounted for as pentanol
    and pentanal; a DB-5MS run starting at 30 C with SPME would show pentane poorly.
11. **Every barrier and rate constant is DFT** (B3LYP/6-311G(d,p), gas phase, isopropyl hydroperoxide
    as the ROOH model); none may enter a parameter file. The qualitative claim that rests on them
    (methyl 9-oxononanoate is not mainly from the 9-alkoxyl's direct scission) is recorded as a
    mechanism annotation only. Miyazaki 2023 makes the analogous claim for hexanal from 13-HpODE on
    product-pattern grounds; the two papers agree that R18b's "aldehyde-side" products are real but
    reach them through epoxy- / dihydrofuran-hydroperoxides.
12. **The dihydrofuran intermediate has two different parents in the two anchors**: Yao draws it from
    the 9-hydroperoxide plus a peroxyl radical (Figure 1d), Miyazaki from the 13-hydroperoxide's
    alkoxyl. Since 9- and 13-HpODE interconvert on heating (Miyazaki) and both are present here, the
    net rule R31 can stay keyed to either; the positive-control list should mention that Yao's
    parent is LOOH_9.
13. Methyl ester, neat, no water: the ester/acid difference does not touch the scissions; the
    absence of water and of metals does bear on the hydroperoxide decomposition rate but nothing
    here quantifies it.

## 5b. Rule sketches (repository suggestions, not the paper's; SMILES are mine)

Structures used below: dihydrofuran C `CCCCCC1C=CC(O1)C(OO)CCCCCCCC(=O)OC`; epoxide A
`CCCCC/C=C\C(OO)C1OC1CCCCCCCC(=O)OC`; epoxide B `CCCCCC1OC1C(OO)/C=C\CCCCCCCC(=O)OC`; 8-OOH
`CCCCC/C=C\C/C=C\C(OO)CCCCCCC(=O)OC`; 10-OOH `CCCCC/C=C\CC(OO)/C=C/CCCCCCC(=O)OC`; methyl
11-oxo-9-undecenoate `O=C/C=C/CCCCCCCC(=O)OC`; methyl 10-oxo-8-decenoate `O=C/C=C/CCCCCCC(=O)OC`;
methyl 8-hydroxyoctanoate `OCCCCCCCC(=O)OC`; methyl 8-oxooctanoate `O=CCCCCCCC(=O)OC`;
5-pentyl-2-furaldehyde `CCCCCc1ccc(C=O)o1`; 2-heptenal `CCCC/C=C/C=O`; 4-hydroperoxy-2-nonenal
`CCCCCC(OO)/C=C/C=O`; methyl heptanoate `CCCCCCC(=O)OC`; 1-pentanol `CCCCCO`; pentanal `CCCCC=O`.

**Y1. R18a stands; annotate the alkyl radical's three fates.** Reactant LOOH_9_ct -> DECADIENAL +
the C8 ester radical, which the paper measures as ME_OCTANOATE (94.8-117.5), methyl
8-hydroxyoctanoate (25-30) and methyl 8-oxooctanoate (17-23) mmol/kg; likewise LOOH_13_ct ->
ME_13_OXO_TRIDECADIENOATE + pentyl, measured as 1-pentanol (8.7-9.7) and pentanal (2.1-2.9), PENTANE
not reported. If the lane needs the alcohol/aldehyde branch, add a net "alkyl radical + ROOH ->
alcohol" product (Miyazaki S6 says the same for 1-octen-3-ol) rather than an aldehyde reduction.
- positive (unchanged): `CCCCC/C=C\C=C\C(OO)CCCCCCCC(=O)OC` -> `CCCCC/C=C/C=C/C=O` + `CCCCCCCC(=O)OC`
- optional second product set: `OCCCCCCCC(=O)OC` and `O=CCCCCCCC(=O)OC` as fates of the same radical
- negative (unchanged): `CCCCCC=O` -> no fire.

**Y2. R31 anchor extension: the intermediate is detected here.** Add to R31's source anchor: "Yao
2024 detects methyl 9-hydroperoxy-9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate (0.59 mmol/kg, neat
methyl linoleate, 180 C, 1 h) and draws its 9-alkoxyl scission to methyl 9-oxononanoate + the
2-pentyl-dihydrofuranyl radical (Fig. 7(3)); Yao's parent is the 9-hydroperoxide plus a peroxyl
radical, Miyazaki's the 13-hydroperoxide." No SMIRKS change. Do NOT add LOOH_9_ct -> 2-pentylfuran +
methyl 9-oxononanoate as an R31 positive control: R31's SMIRKS makes the furan at the OOH-bearing
carbon and the aldehyde at the far diene carbon, so on LOOH_9_ct it already fires to HEXANAL +
ME_8_FURYL_OCTANOATE (Miyazaki's route LA-9-F, the existing second control). Yao's LOOH_9 route is a
different transformation: the ring oxygen comes from a peroxyl radical added at C13 and the 9-OOH
becomes the aldehyde, giving the SAME product pair as LOOH_13 (2-pentylfuran + methyl
9-oxononanoate). Both can be true. Keep R31 keyed to LOOH_13, carry the 9/13 interconversion
(Miyazaki flag 5), and record Yao's variant in the anchor text only:
- Yao's net, for the record: `CCCCC/C=C\C=C\C(OO)CCCCCCCC(=O)OC` + ROO• -> `CCCCCc1ccco1` +
  `O=CCCCCCCCC(=O)OC` (bimolecular in a peroxyl radical; not a one-reactant SMIRKS).

**Y3. New net rule candidate, "dihydrofuran C, other scission": R31-intermediate -> 5-pentyl-2-
furaldehyde + C8 ester radical (-> methyl octanoate / 8-hydroxyoctanoate).** Reactant LOOH_13_ct (or
LOOH_9_ct per Yao). Change: as R31 but the C8-C9 bond breaks instead of C9-C10; C9 stays on the ring
as CHO; ring aromatises.
- positive: `CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)OC` -> `CCCCCc1ccc(C=O)o1` + `CCCCCCCC(=O)OC`
- negative: oleate 9-OOH `CCCCCCC/C=C/C(OO)CCCCCCCC(=O)OC` -> no fire; `CCCCCC=O` -> no fire.
Evidence: 5-pentylfurfural detected (m/z 97), level not printed; mechanism_drawn. Terminal.

**Y4. Epoxy-hydroperoxide route (net) to methyl 11-oxo-9-undecenoate; the mirror to 2-octenal is my
inference.** Reactant LOOH_13_ct. Change: 13-O bridges C12-C13 (epoxide), O2 lands on C11, C11-C12
breaks, C11 becomes CHO. Products: methyl 11-oxo-9-undecenoate (measured 10.9-12.3 mmol/kg) + a
2-pentyloxiranyl fragment (fate not drawn; 2-heptenal would be the ring-opened aldehyde, not
stated by the paper).
- positive: `CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)OC` -> `O=C/C=C/CCCCCCCC(=O)OC` + C7 fragment
- mirror (mine, not the paper's): LOOH_9_ct -> `CCCCC/C=C/C=O` (2-octenal) + C10 epoxy-ester fragment
- negative: oleate 9-OOH (no conjugated partner carbon for the epoxide) -> no fire; hexanal -> no fire.
Status: mechanism_drawn; the intermediate class (epoxides A + B, 2.47 mmol/kg) is the largest
non-volatile pool measured, which is the one quantitative argument for carrying this route at all.

**Y5. R32 second anchor.** Fig. 7(6) and the text: 10-hydroperoxide -> methyl 10-oxo-8-decenoate +
2-octenyl radical -> 1-octen-3-ol; both products detected here at 180 C in the ester (levels
figure-only). Add "yao2024: Fig. 7(6), ester, 180 C" to R32's anchor; no SMIRKS change (R32 already
uses the ester-agnostic pattern).

**Y6. 8-hydroperoxide side-A analogue (new, low priority).** 8-OOH (0.19 mmol/kg) -> methyl
heptanoate (5.6-6.1 mmol/kg) + a C11 dienal (2,4-undecadienal, not reported). Only the ester is
measured; the paper draws no scheme; mechanism by analogy with R18a. Positive: 8-OOH SMILES above ->
`CCCCCCC(=O)OC` + `CCCCC/C=C/C=C/CC=O`; negative: LOOH_9_ct -> must not give methyl heptanoate.
