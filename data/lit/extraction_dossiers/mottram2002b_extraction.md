# Mottram & Elmore 2002 (ACS Symp. Ser. 826, ch. 5) — EXTRACTION (review chapter; 2,4-alkadienal + H2S -> 2-alkylthiophene and 2-alkyl-2H-thiapyran drawn; one table of ng for the alkylthiophenes and alkylthiapyrans in cooked beef and lamb and in cysteine + ribose + methyl linoleate or linolenate pots, 140 C / 30 min)
### The group's summary of the lipid–Maillard sulfur heterocycles, with the only mass-unit numbers on disk for the 2-alkylthiophene / 2-alkylthiapyran pair from a defined fatty acid.

**Source on disk:** `data/articles/mottram2002b.pdf` (owner's download, 2026-09-08). Born-digital
ACS chapter, clean `pypdf` text layer; Table I checked against a rendering of p. 99. Figures 1-2
are generic structures (described), Figure 3 is the drawn dienal + H2S scheme (rendered, described
atom by atom in §4), Figure 4 is a list of aldehyde structures (transcribed in §3). Not to be
confused with `mottram2002_extraction.md`, which is a different 2002 paper.

## 0. Identity

| field | value |
|---|---|
| Title | "Novel Sulfur Compounds from Lipid-Maillard Interactions in Cooked Meat" |
| Authors | Donald S. Mottram, J. Stephen Elmore (School of Food Biosciences, University of Reading) |
| Venue | Heteroatomic Aroma Compounds, ACS Symposium Series 826, chapter 5, pp. 93-101, American Chemical Society, 2002 (publication date 7 Aug 2002) |
| DOI | 10.1021/bk-2002-0826.ch005 |
| Kind | review of the group's own work with one new table; the meat numbers come from refs 15, 16, 18, 23 (Elmore et al. 1997, 1999, 2000), the thiazoline mechanism from ref 20 (`elmore1997_extraction.md`), the alkylthiophene / thiapyran identifications from ref 21 (`farmer1990_extraction.md`), the pure-dienal pot from ref 22 (van den Ouweland, Demole, Enggist 1989, not on disk) |
| Naming | "2-alkyl-(2H)-thiapyran" = 2-alkyl-2H-thiopyran, a six-membered S ring with two C=C; isomeric with the 2-alkylthiophene carrying one more CH2 in the chain (both C10H16S for pentylthiapyran / hexylthiophene, MW 168) |

## 1. Why it matters

It closes the loop the 1988 and 1990 papers left open: Figure 3 draws, for a generic 2,4-alkadienal,
both H2S routes (1,6-addition -> thiapyran; 1,4-addition -> thiophene), and Table I gives the two
product families in ng from cysteine + ribose + a single fatty acid methyl ester (18:2 n-6 or
18:3 n-3), which is the closest thing in the group's work to a yield from a defined lipid. It also
states the omega-3 / omega-6 mapping the rule writer needs (18:3 -> 2,4-heptadienal -> 2-ethylthiapyran
and 2-propylthiophene; 18:2 -> 2,4-decadienal -> 2-pentylthiapyran and 2-hexylthiophene), and it
records two discrepancies that any rule must carry as flags: in pure dienal + H2S pots the thiapyran
exceeds the isomeric thiophene by up to 100x, while in meat the two are about equal (so the
thiophenes probably have a second source); and the 18:2 model pot in Table I made 2-pentylthiophene,
not 2-pentylthiapyran, the opposite of Farmer 1990's lecithin pot (flag 3).

## 2. Methods as they matter to a model

- **Model systems (Table I, right-hand columns):** cysteine + ribose + fatty acid methyl ester,
  "0.5 mmol each"; methyl linoleate (18:2 n-6) in one pot, methyl alpha-linolenate (18:3 n-3) in the
  other; "heated at 140 C in a sealed container for 30 min"; headspace by SPME (fibre, time,
  temperature not stated); GC-MS. **Volume, buffer, pH, water content and vessel are not stated.**
  Quantities are "ng per 0.5 mmol fatty acid" in the headspace above the mixture (the calibration
  basis is not stated). Replication not stated.
- **Meat (Table I, left columns):** cooked beef and lamb from animals on control, linseed or fish-oil
  diets (refs 15, 16, 23); headspace quantities in ng per 100 g meat.
- **Pure dienal pots (text only):** "2,4-alkadienals were reacted with hydrogen sulfide in aqueous
  solution"; no conditions given; for 2,4-decadienal the cited ref 22 used pH 8. Result stated in
  words only: thiapyrans up to 100 times the isomeric thiophenes.
- **Thiazolines and thiazoles in meat:** 46 alkyl-3-thiazolines in beef, 12 in lamb, higher on the
  PUFA diets; no numbers in this chapter (refs 15, 16, 18).
- Odour: thiazolines / thiazoles with 2-n-alkyl chains "slightly fatty", not detected by GC-O;
  thiapyrans weak, thiophene-like, the small ones garlic-like. No thresholds measured.

## 3. Tables re-typed

### Table I. "Quantities of Alkylthiophenes and Alkylthiapyrans Found in Beef and Lamb Fed on PUFA Rich Diets, Compared with Quantities Produced in Model Reaction Systems." Meat: ng per 100 g meat (headspace). Model: ng per 0.5 mmol fatty acid (headspace). "-" = not detected, tr = trace.

| compound | lamb control | lamb linseed | lamb fish | beef control | beef linseed | beef fish | model 18:2 | model 18:3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **2-alkylthiophenes** | | | | | | | | |
| 2-ethylthiophene | 13 | 14 | 32 | 5 | 6 | 15 | 3 | 99 |
| 2-propylthiophene | 2 | 3 | 4 | tr | tr | tr | 7 | 9 |
| 2-butylthiophene | 2 | 4 | 4 | - | tr | tr | 8 | tr |
| 2-pentylthiophene | 4 | 3 | 4 | tr | tr | tr | **21** | tr |
| 2-hexylthiophene | 4 | 4 | 5 | tr | tr | tr | 9 | tr |
| **2-alkyl-(2H)-thiapyrans** | | | | | | | | |
| 2-methylthiapyran | 3 | 3 | 6 | - | - | - | - | 19 |
| 2-ethylthiapyran | 22 | 59 | 131 | 3 | 6 | 14 | - | **399** |
| 2-propylthiapyran | 2 | 4 | 15 | - | tr | 2 | - | 3 |
| 2-butylthiapyran | tr | tr | 2 | - | - | - | tr | - |
| 2-pentylthiapyran | 5 | 4 | 8 | - | tr | tr | tr | - |

Footnote as printed: "Model systems: cysteine + ribose + fatty acid methyl ester (18:2 n-6 or 18:3
n-3); 0.5 mmol each."

### Figure 4 (transcribed structures). "Some aldehydes formed in the autoxidation of linolenic and linoleic acids."

| parent | aldehydes drawn |
|---|---|
| linolenic 18:3 n-3, CH3-CH2-CH=CH-CH2-CH=CH-CH2-CH=CH-(CH2)7-COOH | 2-butenal CH3-CH=CH-CHO; 3-hexenal CH3-CH2-CH=CH-CH2-CHO; 2-pentenal CH3-CH2-CH=CH-CHO; 2-hexenal CH3-CH2-CH2-CH=CH-CHO; 2,4-heptadienal CH3-CH2-CH=CH-CH=CH-CHO; 2,4,7-decatrienal CH3-CH2-CH=CH-CH2-CH=CH-CH=CH-CHO |
| linoleic 18:2 n-6, CH3(CH2)4-CH=CH-CH2-CH=CH-(CH2)7-COOH | alkanals CH3(CH2)n-CHO, n = 3-5 (pentanal, hexanal, heptanal); 2-alkenals CH3(CH2)n-CH=CH-CHO, n = 3-6 (2-heptenal to 2-decenal); 2,4-decadienal CH3(CH2)4-CH=CH-CH=CH-CHO |

### Figures 1-2 (structures only)

Figure 1: alkyl-3-thiazolines found in cooked beef and lamb, R1 = H, CH3 or C2H5 (4-position), R2 =
CH3 or C2H5 (5-position), 2-CnH2n+1 with n = 3-9 and 15; alkylthiazoles with R1 and R2 = CH3 or C2H5,
n = 4-8 and 15. Figure 2: 2-alkylthiophenes (ring-CH2-CnH2n+1) and 2-alkyl-2H-thiapyrans
(ring-CnH2n+1), n = 1-6 in both series as drawn (six of each found in beef and lamb; Table I lists
five of each).

## 4. Routes and numbers the repository can use

| route or quantity | reactant -> product | mechanism as drawn (Figure 3) | measured numbers, units, conditions | evidence class |
|---|---|---|---|---|
| M02-A dienal -> thiapyran | R-CH=CH-CH=CH-CHO + H2S -> 2-R-2H-thiapyran + H2O | top branch: H2S adds across the diene's far end, S onto the R-bearing carbon (C5 of the dienal, counting CHO as C1) and H onto C4, giving R-CH(SH)-CH2-CH=CH-CHO (drawn with the arrow from the aldehyde O pulling electrons, i.e. a 1,6-conjugate addition); the S lone pair then attacks the aldehyde carbon C1, closing a six-membered ring S-C1(OH)-C2=C3-C4H2-C5(R); loss of water from C1 gives the 2H-thiapyran with R on C2 (the former C5) | 18:3 pot: 2-ethylthiapyran 399, 2-methylthiapyran 19, 2-propylthiapyran 3 ng per 0.5 mmol FAME; 18:2 pot: 2-butyl tr, 2-pentyl tr, others nd. Pure dienal + H2S (text): thiapyran up to 100x the isomeric thiophene | mechanism_drawn + level_only |
| M02-B dienal -> thiophene | R-CH=CH-CH=CH-CHO + H2S -> 2-(R-CH2)-thiophene + H2O | bottom branch: H2S adds S onto C4 and H onto C5 (1,4-conjugate addition), giving R-CH2-CH(SH)-CH=CH-CHO; S attacks C1, closing a five-membered ring S-C1(OH)-C2=C3-C4(CH2R); dehydration gives the thiophene with CH2R on C2 (the former C4) | 18:2 pot: 2-pentylthiophene 21, 2-hexylthiophene 9, 2-butylthiophene 8, 2-propylthiophene 7, 2-ethylthiophene 3 ng per 0.5 mmol; 18:3 pot: 2-ethylthiophene 99, 2-propylthiophene 9 | mechanism_drawn + level_only |
| M02-C chain bookkeeping | a Cn 2,4-alkadienal gives 2-(Cn-5)-thiapyran and 2-(Cn-4)-thiophene | follows from A and B: the thiapyran ring takes C1-C5, the thiophene ring takes C1-C4 | 2,4-heptadienal (from 18:3) -> 2-ethylthiapyran + 2-propylthiophene: 399 and 9 (ratio 44); 2,4-decadienal (from 18:2) -> 2-pentylthiapyran + 2-hexylthiophene: tr and 9 (see flag 3) | within_study_ratio |
| M02-D omega-3 / omega-6 fingerprint | 18:3 -> short-chain (C1-C3) thiapyrans and C2-C3 thiophenes; 18:2 -> C3-C6 thiophenes, few thiapyrans | Figure 4 aldehyde slate feeding A and B | Table I: the 18:3 column has 2-methyl / ethyl / propyl thiapyrans (19 / 399 / 3) and nothing longer; the 18:2 column has 2-propyl to 2-hexylthiophene (7-21) and no C1-C3 thiapyran; meat on linseed / fish oil diets shows the 18:3 fingerprint (lamb fish: 2-ethylthiapyran 131 vs control 22 ng/100 g) | level_only |
| M02-E thiophene second source | meat thiophene:thiapyran about 1:1 while the dienal pot gives up to 1:100 | stated: "a mechanism other than that given in Figure 3 may be required to explain the formation of alkylthiophenes" (candidates named elsewhere in the group's work: H2S on 2-alkylfurans, on 4-oxo-alkanals, on 2-alkenals) | qualitative | level_only |
| M02-F dienal sink | 2,4-alkadienal + H2S -> thiapyran removes a potent aroma dienal and makes a weak odorant | stated | none | level_only |
| M02-G thiazoline formation (review of ref 20) | alpha-hydroxycarbonyl (1-hydroxy-2-butanone, 3-hydroxy-2-butanone) + (NH4)2S + n-alkanal or Strecker aldehyde (2-methylbutanal) -> 3-thiazoline with the alkanal chain at C2 | see `elmore1997_extraction.md` §4 | none in this chapter | mechanism_drawn (elsewhere) |

Mass-unit arithmetic (mine, for scale only): 21 ng 2-pentylthiophene (M 154.3) per 0.5 mmol
methyl linoleate is 0.14 nmol, i.e. 2.7e-7 mol per mol fatty acid in the sampled headspace; 399 ng
2-ethylthiapyran (M 126.2) per 0.5 mmol methyl linolenate is 3.2 nmol, 6.3e-6 mol per mol. These are
headspace amounts recovered by SPME, not totals in the pot, so they are lower bounds of unknown
factor.

## 5. Rule sketches (reactant -> product in words; controls are mine, SMILES mine)

**S1. 2,4-alkadienal + H2S -> 2-alkyl-2H-thiapyran + H2O (M02-A; net of 1,6-addition, ring
closure, dehydration; terminal).** Required substructure: O=CH-CH=CH-CH=CH-C (a conjugated
2,4-dienal with at least one carbon beyond C5).
- positive: 2,4-decadienal `CCCCC/C=C/C=C/C=O` (key `DECADIENAL`) + `S` -> 2-pentyl-2H-thiapyran `CCCCCC1C=CC=CS1`; 2,4-heptadienal `CC/C=C/C=C/C=O` + `S` -> 2-ethyl-2H-thiapyran `CCC1C=CC=CS1` (399 ng in the 18:3 pot).
- negative: (E)-2-nonenal `CCCCCC/C=C/C=O` + `S` -> no fire (one C=C only); hexanal `CCCCCC=O` + `S` -> no fire; 2,4-hexadienoic (sorbic) acid `C/C=C/C=C/C(=O)O` + `S` -> no fire (acid, not aldehyde).

**S2. 2,4-alkadienal + H2S -> 2-(alkyl-CH2)-thiophene + H2O (M02-B; net of 1,4-addition, ring
closure, dehydration; terminal).** Same required substructure as S1; the two rules are the two
regiochemistries of the same addition and should be written as a pair.
- positive: `CCCCC/C=C/C=C/C=O` + `S` -> 2-hexylthiophene `CCCCCCc1cccs1`; 2,4-nonadienal `CCCC/C=C/C=C/C=O` + `S` -> 2-pentylthiophene `CCCCCc1cccs1` (Elmore 1997: 21 % of area); `CC/C=C/C=C/C=O` + `S` -> 2-propylthiophene `CCCc1cccs1`.
- negative: as S1; and 2-pentylfuran `CCCCCc1ccco1` (`PENTYLFURAN`) + `S` -> must NOT fire under this rule (the furan -> thiophene exchange is a different, doubted route: Whitfield 1988's 7:1 versus 1:10 ratio argument).

**Branch note (mine):** in a pure dienal + H2S pot S1 dominates S2 by up to 100:1 (text); in
cysteine + ribose + methyl linoleate at 140 C / 30 min Table I shows the reverse (thiophenes 7-21
ng, thiapyrans trace). Do not encode a branch ratio from this chapter; record both as hypotheses.

## 6. Flags

1. **Review chapter with one new table.** The model-pot conditions (volume, pH, buffer, headspace,
   SPME fibre and calibration, replication) are not given; "ng per 0.5 mmol fatty acid" is a
   headspace-SPME quantity of unstated basis. Treat Table I model columns as level_only.
2. **Meat and model units differ** (ng/100 g meat vs ng per 0.5 mmol FAME); only within-column
   patterns are comparable.
3. **Internal contradiction on 2,4-decadienal.** The chapter says 2-pentylthiapyran and
   2-hexylthiophene were "major products" of decadienal + H2S at pH 8 (ref 22) and that pure dienal
   pots give thiapyran >> thiophene; yet the 18:2 model column has 2-pentylthiapyran "tr" and
   2-hexylthiophene 9 ng, with 2-pentylthiophene (from 2,4-nonadienal) the largest. Farmer 1990's
   egg-PC pot (`farmer1990_extraction.md` Table 2) had 2-pentylthiapyran 28x 2-hexylthiophene by
   relative area. Either the 18:2 pot made little decadienal under these conditions, or SPME
   under-samples the thiapyran, or the ratio is condition-dependent; the chapter does not say.
4. **The 18:3 column is dominated by one compound** (2-ethylthiapyran 399 ng, 2-ethylthiophene 99
   ng); the ratio 4:1 there is far from the 100:1 claimed for pure dienal pots.
5. **"n = 1-6" in Figure 2** versus "C2-C7 thiapyrans and C3-C8 thiophenes" quoted from ref 21 in
   the text; the figure's n counts the chain beyond the ring CH2 for the thiophenes, so the two
   statements are consistent but easy to misread.
6. **No thiazoline numbers** in this chapter, despite its title; those are in refs 15, 16, 18 (not
   on disk) and the mechanism in `elmore1997_extraction.md`.
7. **No kinetic content anywhere:** one temperature, one time, no time course, no replicate, no
   error bar.
