# Miyazaki et al. 2023 — EXTRACTION (six purified linoleic acid hydroperoxide isomers heated neat, 120 C / 5 min; SPME-GC-MS volatiles, LC-MS/MS non-volatiles)
### The only paper in the corpus that decomposes each HpODE positional/geometric isomer separately and prints an isomer-by-product table.

**Source on disk:** `data/articles/Miyazaki2023.pdf` (owner's download, 2026-09-08). Read from the
`pypdf` text layer; Table 1's column alignment is lost in the plain text layer, so page 6 was
re-extracted in layout mode and cross-checked against a rendering of the page. Figures 1-3 are drawn
reaction schemes and are described in words below (no numbers read off them). Figures 4-6 are
chromatograms / product-ion spectra: FIGURE-ONLY. Supplementary material 1-2 was not available.

## 0. Identity

| field | value |
|---|---|
| Title | "Elucidation of decomposition pathways of linoleic acid hydroperoxide isomers by GC-MS and LC-MS/MS" |
| Authors | Ruriko Miyazaki, Shunji Kato, Yurika Otoki, Halida Rahmania, Masayoshi Sakaino, Shigeo Takeuchi, Toshiro Sato, Jun Imagi, Kiyotaka Nakagawa (Tohoku Univ.; J-Oil Mills) |
| Venue | Biosci. Biotechnol. Biochem. 2023, 87(2), 179-190 |
| DOI | 10.1093/bbb/zbac189 |
| Naming | HpODE = hydroperoxy octadecadienoic acid (free acid). 9-EZ = 9-hydroperoxy-10E,12Z; 9-EE = 10E,12E; 13-ZE = 13-hydroperoxy-9Z,11E; 13-EE = 9E,11E; 10-HpODE and 12-HpODE are the non-conjugated singlet-oxygen isomers (Figure 1a). Compound numbers in bold below are the paper's scheme numbers. |
| Scission nomenclature (theirs) | "alpha-scission" = cleavage of the C(n-1)-Cn bond on the CARBOXYL side of the alkoxyl carbon Cn; "beta-scission" = cleavage of Cn-C(n+1) on the METHYL side. This is NOT the repo's R18a/R18b "side A / side B" convention; see §5. |
| Companions | Kato et al. 2018 (npj Sci Food 2), Rahmania et al. 2020 (Sci Rep 10): the six isomers in canola / soybean / rice bran oil; Kato et al. 2022 (npj Sci Food 6): standard preparation |

## 1. Why it matters

The repo's lipid rules R18a/R18b take the linoleate hydroperoxide product slate from Frankel 1989 via
schroen2022 and treat "13-OOH -> hexanal" and "9-OOH -> 2,4-decadienal" as simple beta-scissions.
This paper heats each of the six food-relevant HpODE isomers separately and shows (i) which volatiles
each isomer gives (Table 1), (ii) that hexanal comes abundantly from 9-HpODE as well as from
13-HpODE, (iii) that 2-pentylfuran comes from 9-, 10- and 13-HpODE but NOT from 12-HpODE, (iv) that
10-HpODE is the isomer that gives the C8 set (1-octen-3-ol, 1-octen-3-one, 2-octenal, 2-octen-1-ol),
(v) that 12-HpODE gives almost only 2-heptenal, and (vi) that 9- and 13-HpODE interconvert on heating
(and 10 -> 8, 12 -> 14). It also argues, against the textbook, that hexanal from 13-HpODE does not
come from direct C12-C13 scission (vinyl radical) but from an epoxy-allyl / vinyl-ether-hydroperoxide
route, and that 2-pentylfuran comes from a furyl-hydroperoxide made by alkoxyl cyclisation onto the
gamma-carbon, not from 1-hydroperoxy-1,3-nonadiene. For the rule writer the net transformations are
what matter; the mechanistic re-assignment changes the annotation, not the products (§5).

## 2. Methods as they matter to a model

- **Substrate:** each HpODE isomer standard, purity > 95 % by LC-MS (m/z 50-400), prepared from
  linoleic acid (Nacalai) per Kato et al. 2022. Free acids, not esters.
- **Charge:** 100 µg of one isomer per 2 mL amber screw-cap vial; solvent evaporated under N2; vial
  sealed under ambient air. Neat film, no matrix, no water, no added metal, no antioxidant.
- **Heating:** 120 C for 5 min, then cooled immediately; analysed at once. Authors: "Somewhat lower
  temperature was set to avoid the rapid decomposition of intermediate substances by excessive high
  temperature (~180 C)." No time series, no replicate: **Table 1 is n = 1.**
- **Volatiles:** HS-SPME (fibre type not stated) at 40 C / 20 min, fibre 1 cm above the sample;
  desorbed 250 C / 4.5 min; GC-EI-MS (Shimadzu QP2010 SE), DB-WAX-UI 60 m x 0.25 mm x 0.25 µm; He 25
  cm/s; 30 C (10 min) -> 250 C at 5 C/min -> 250 C (5 min); scan m/z 41-500. Identification: NIST 17
  spectral match only (no retention-index confirmation, no authentic standards stated).
- **Quantification: NONE.** Table 1 reports raw **peak areas** (arbitrary units). No internal standard,
  no calibration, no response factors. Areas are comparable only within one column (one isomer, one
  fibre exposure) and only loosely across columns.
- **Non-volatiles:** LC-MS/MS, 6500 QTRAP; COSMOSIL 2.5C18-MS-II 2.0 x 100 mm, 40 C; A = H2O + 0.1 %
  HOAc, B = MeOH + 0.1 % HOAc, 0.2 mL/min; post-column 0.2 mM sodium acetate in MeOH to force [M+Na]+.
  Sample: decomposed vial dissolved in 1 mL MeOH, diluted 2x with H2O, 10 µL injected. Detected as
  Q1 XIC of m/z 335 ([HpODE+Na]+) and m/z 351 ([LA+3O+Na]+, the "cyclised hydroperoxides"), with
  product-ion scans. Qualitative only.
- **Isomer distribution in oils:** NOT measured here; the paper cites Kato 2018 / Rahmania 2020 for
  the statement that these six isomers are the HpODE found in canola, soybean and rice bran oil.

## 3. Tables re-typed

### Table 1. "Volatile compounds derived from each HpODE isomer (n = 1)." Peak area (arbitrary units); blank = not listed for that isomer.

Compound No. is the paper's scheme number. Two retention times under one compound are two geometric
isomers (see notes). Column placement was confirmed against the rendered page.

| RT (min) | compound | No. | 9-EZ-HpODE | 9-EE-HpODE | 10-HpODE | 12-HpODE | 13-ZE-HpODE | 13-EE-HpODE |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 7.3 | 2-propenal (acrolein) | 45 | | | 14 156 958 | 3 849 489 | | |
| 14.3 | pentanal | 19 | 7 782 157 | 11 024 828 | 15 239 960 | | 9 609 917 | 3 643 309 |
| 19.5 | **hexanal** | 5 | 67 492 032 | 89 039 976 | 28 479 695 | | 139 131 593 | 57 428 632 |
| 24.7 | 2-hexenal | 52 | | | | 3 830 714 | | |
| 25.2 | **2-pentylfuran** | 23 | 2 571 331 | 3 729 530 | 1 683 610 | | 2 547 384 | 3 423 267 |
| 25.9 | 1-pentanol | 21 | 4 340 288 | 6 760 388 | 11 465 119 | | 7 841 460 | 4 216 503 |
| 27.6 | 1-octen-3-one | 38 | | | 22 591 928 | | | |
| 28.3 | **2-heptenal** | 49 | 3 689 599 | 5 097 082 | 3 578 519 | 264 993 230 | 10 441 954 | 10 307 237 |
| 28.6 | 4-nonene | 53 | 5 891 668 | | | | | |
| 30.7 | 2-octenal (2Z, per authors) | 13 | 51 730 240 | 3 790 551 | 6 089 884 | | 23 977 557 | 19 633 041 |
| 31.5 | 2-octenal (2E, per authors) | 13 | 36 356 004 | 126 936 201 | 11 686 096 | | | |
| 30.8 | 3,5-octadien-2-ol | 54 | | | | 5 126 715 | | |
| 31.5 | **1-octen-3-ol** | 43 | | | 144 784 414 | | | |
| 33.0 | formic acid | 55 | 19 693 120 | 26 314 876 | 4 575 490 | | 15 949 697 | 6 454 059 |
| 33.2 | 2,4-heptadienal | 56 | | | | 1 868 437 | | |
| 34.3 | 2-propyltetrahydrofuran | 57 | 4 364 549 | 4 731 185 | | | | |
| 35.7 | 2,4-octadienal | 58 | 2 435 708 | | | | | |
| 36.3 | 2-octen-1-ol | 46 | | | 15 918 439 | | | |
| 37.5 | vinyl hexanoate | 59 | 8 152 529 | 14 161 612 | | | 50 250 033 | 28 890 419 |
| 38.5 | 2,4-nonadienal | 60 | 2 467 057 | 3 096 982 | | | | |
| 39.9 | **2,4-decadienal** (2E,4Z, per authors) | 3 | 9 051 292 | 998 262 | | | 2 727 622 | 2 384 563 |
| 41.0 | 2,4-decadienal (2E,4E, per authors) | 3 | 7 039 238 | 7 810 015 | | | | |
| 40.9 | 1-undecyn-4-ol | 61 | | | 12 632 207 | | | |
| 41.4 | hexanoic acid | 6 | 9 641 678 | 15 214 035 | 7 429 005 | | 29 711 559 | 26 354 349 |
| 44.7 | 4,5-epoxy-2-decenal | 35 | 4 876 361 | 5 157 803 | | | 3 662 424 932 (sic) | 3 982 179 |
| 45.1 | 4,5-epoxy-2-decenal (2nd peak) | 35 | 12 296 791 | 14 144 874 | | | 9 325 191 | 10 776 651 |
| 46.0 | octanoic acid | 1 | 39 195 107 | 31 343 775 | | 48 983 221 | 20 369 657 | 23 497 489 |

Notes on the table: (a) The 13-ZE entry for 4,5-epoxy-2-decenal at 44.7 min is printed as a
ten-digit number "3 662 424 932", 25x larger than anything else in the table and 400x its own second
peak; almost certainly a typesetting error (a dropped space or an extra digit group). Do not use.
(b) Geometric assignments of the 2-octenal and 2,4-decadienal pairs are the authors' inference from
relative abundance in the EZ vs EE columns ("most likely"), not from standards. (c) 3-nonenal was
looked for and not found in any column. (d) No blank/control vial is reported.

### Isomer-specific product summary (from Table 1; areas are within-column relative)

| isomer | dominant volatiles (largest areas first) | absent / notable |
|---|---|---|
| 9-EZ-HpODE | hexanal, 2-octenal (both peaks), octanoic acid, formic acid, 4,5-epoxy-2-decenal, 2,4-decadienal (both peaks), vinyl hexanoate | 2-pentylfuran present (2.6 M); no C8 alcohols/ketone |
| 9-EE-HpODE | 2-octenal (2E peak 127 M), hexanal, octanoic acid, formic acid, 4,5-epoxy-2-decenal, hexanoic acid, vinyl hexanoate, pentanal | 2,4-decadienal (2E,4Z) 40x smaller than in 9-EZ |
| 10-HpODE | **1-octen-3-ol** (145 M), hexanal, 1-octen-3-one, 2-octen-1-ol, pentanal, 2-propenal, 1-undecyn-4-ol, 2-octenal, 1-pentanol | no 2,4-decadienal, no octanoic acid, no 4,5-epoxy-2-decenal |
| 12-HpODE | **2-heptenal** (265 M, ~85 % of its listed area), octanoic acid, 3,5-octadien-2-ol, 2-propenal, 2-hexenal, 2,4-heptadienal | no hexanal, no 2-pentylfuran, no pentanal, no 2-octenal |
| 13-ZE-HpODE | **hexanal** (139 M), vinyl hexanoate, hexanoic acid, 2-octenal (2Z peak), octanoic acid, formic acid, 2-heptenal, pentanal, 4,5-epoxy-2-decenal | 2,4-decadienal present (2.7 M) -> evidence of 13 -> 9 isomerisation |
| 13-EE-HpODE | hexanal, vinyl hexanoate, hexanoic acid, octanoic acid, 2-octenal, 4,5-epoxy-2-decenal, 2-heptenal | similar to 13-ZE at lower area |

### Figure 1b — the four predicted pathway types (drawn; described)

(1) Scission: alkoxyl radical -> alkyl radical + conjugated aldehyde (type i, e.g. 9-HpODE) or
conjugated aldehyde + allyl radical (type ii, e.g. 10-HpODE); the allyl radical is drawn being
further oxidised to a hydroperoxide (iii), or delocalising (iv) and then oxidising (v).
(2) Epoxy-allyl radical: the alkoxyl oxygen closes onto the adjacent alkene carbon giving an epoxide
with the radical on the allylic carbon (vi); further oxidation gives an epoxy-hydroperoxide (vii)
which cleaves (viii) to an epoxy fragment + a 2-alkenal; or the oxirane C-C bond breaks (ix) giving a
vinyl-ether radical, oxidised (x) to a vinyl-ether hydroperoxide that cleaves (xi) to an aldehyde
fragment + an alkoxyl fragment.
(3) Furyl-hydroperoxide: the alkoxyl oxygen attacks the gamma-carbon (xii) closing a five-membered
oxygen ring (2,5-disubstituted 2,3-dihydrofuran-type radical, "furyl radical"); further oxidation
(xiii) gives a furyl-hydroperoxide; scission (xiv) gives an alkylfuran + an aldehyde fragment.
(4) Hydroperoxyl transfer / isomerisation: •OH abstracts H from ROOH giving ROO•, which reversibly
loses O2 to the pentadienyl (or allyl) radical; O2 re-adds at the other terminus (xv: 13 <-> 9; xvi:
12 <-> 14, and 10 <-> 8).

### Figure 2 — 9-HpODE and 13-HpODE (drawn; described; scheme numbers in bold)

- **9-HpODE** (OOH at C9, 10E,12Z) -> alkoxyl at C9 (**7**) -> alpha-scission C8-C9 -> **2,4-decadienal (3)** (C9-C18) + octanoate C8 alkyl radical (**2**) -> **octanoic acid (1)**.
- Alkoxyl **7** -> cyclisation onto C12 (gamma) -> 9,12-oxygen-bridged radical (**8**) -> O2 -> 9,12-furyl-13-hydroperoxy-10-octadecenoic acid (**9**) -> scission C12-C13 -> **hexanal (5)** (C13-C18) + a furan-bearing C1-C12 fragment (**4**, 8-(furan-2-yl)octanoic acid as drawn); hexanal -> **hexanoic acid (6)**.
- Alkoxyl **7** -> epoxy-allyl radical (**10**, 9,10-epoxide, radical over C11-C13) -> O2 at C11 -> 9,10-epoxy-11-hydroperoxide (**11**) -> alkoxyl at C11, scission C11-C12 -> **2-octenal (13)** (C11-C18) + epoxy C1-C10 fragment (**12**).
- **10** -> radical shifted to C13 (**14**) -> O2 -> 9,10-epoxy-13-hydroperoxide (**15**) -> scission C13-C14 -> epoxy-oxo C1-C13 fragment (**16**) + pentyl radical (**17**) -> pentyl hydroperoxide (**18**) -> **pentanal (19)** / **1-pentanol (21)**.
- **13-HpODE** (OOH at C13, 9Z,11E) -> alkoxyl at C13 (**24**) -> beta-scission C13-C14 -> 13-oxo-9,11-tridecadienoic acid (**20**) + pentyl radical (**17**) -> pentanal / 1-pentanol as above.
- Alkoxyl **24** -> cyclisation onto C10 (gamma) -> 10,13-oxygen-bridged radical (**25**, radical at C9) -> O2 -> 9-hydroperoxy-10,13-furyl-11-octadecenoic acid (**26**) -> scission C9-C10 -> **2-pentylfuran (23)** (ring O, C10-C13; pentyl C14-C18) + **9-oxononanoic acid (22)**.
- Alkoxyl **24** -> epoxy-allyl radical (**27**, 12,13-epoxide, radical over C9-C11) -> O2 at C11 -> 11-hydroperoxy-12,13-epoxide (**28**) -> scission -> 11-oxo-9-undecenoic acid (**29**) + epoxy-C5 radical (**30**). **27** -> radical at C9 (**31**) -> O2 -> 9-hydroperoxy-12,13-epoxy-10-octadecenoic acid (**34**) -> scission C8-C9 -> octanoate radical (**2**) + **4,5-epoxy-2-decenal (35)**. **27** -> oxirane C12-C13 bond opens -> vinyl-ether radical (**32**) -> O2 -> 12-((1-hydroperoxyhexyl)oxy)dodeca-9,11-dienoic acid (**33**) -> **hexanal (5)** + 12-oxo-dodecadienoic acid fragment (**36**).
- Left margin: 9-HpODE <-> 9-alkoxyl/peroxyl <-> pentadienyl radical <-> 13-peroxyl <-> 13-HpODE (the isomerisation manifold).

### Figure 3 — 10-HpODE and 12-HpODE (drawn; described)

- **10-HpODE** (OOH at C10, 8E,12Z) -> alkoxyl (**37**) -> beta-scission C10-C11 -> 10-oxo-8-decenoic acid (**39**) + 2-octenyl allyl radical (**40**) <-> 1-octen-3-yl radical (**41**). **41** -> O2 -> 1-octen-3-hydroperoxide (**42**) -> **1-octen-3-ol (43)** (reduction) or **1-octen-3-one (38)**. **40** -> O2 -> 2-octen-1-hydroperoxide (**44**) -> **2-octenal (13)** or **2-octen-1-ol (46)**; **44** -> scission -> **2-propenal (45)** + pentyl radical (**17**) -> pentyl hydroperoxide (**18**) -> pentanal (**19**) / 1-pentanol (**21**). Left margin: 10-HpODE <-> allyl radical <-> **8-HpODE**.
- **12-HpODE** (OOH at C12, 9Z,13E) -> alkoxyl (**47**) -> alpha-scission C11-C12 -> **2-heptenal (49)** (C12-C18) + 9-undecenoic acid allyl radical (**48**) -> rearranged (**50**) -> O2 -> 9-hydroperoxy-10-undecenoic acid (**51**) -> scission -> octanoate radical (**2**) + 2-propenal (**45**); **2** -> octanoic acid (**1**). Left margin: 12-HpODE <-> allyl radical <-> **14-HpODE**.

### LC-MS/MS findings (Figures 4-6; FIGURE-ONLY, qualitative)

- m/z 351 ([LA+3O+Na]+) detected at 8.7-10.7 min from every heated isomer; product ions m/z 191, 195,
  221, 263 from heated 9- and 13-HpODE (O at C9 -> 191/195; O at C13 -> 263; 221 assigned to
  11-hydroperoxy-12,13-epoxy species **28**). Structures of the m/z 351 from 10- and 12-HpODE "could not
  [be] elucidate[d]".
- m/z 335 at ~12 min after heating: 9-HpODE (335 > 195) appears in heated 13-HpODE and 13-HpODE
  (335 > 247) in heated 9-HpODE; ions assigned to 8-HpODE (152, 181, 189) seen in heated 10-HpODE and
  to 14-HpODE (261) in heated 12-HpODE. 8- and 14-HpODE standards were not available; assignments rest
  on Kato 2021 fragmentation rules.

## 4. Routes and numbers the repository can use

Conditions for every row: neat isomer (100 µg), air, 120 C, 5 min, n = 1, peak area units.

| route | reactant -> product | mechanism as drawn (figure) | measured numbers (units, conditions) | evidence class |
|---|---|---|---|---|
| LA-9-A | 9-HpODE -> 2,4-decadienal + octanoic acid (via C8 alkyl radical) | alkoxyl C9, alpha-scission C8-C9 (Fig 2: 7 -> 3 + 2 -> 1) | 2,4-decadienal area 9.05 M + 7.04 M (9-EZ), 1.00 M + 7.81 M (9-EE); octanoic acid 39.2 M / 31.3 M | mechanism_drawn; measured_level (area, n=1) |
| LA-9-F | 9-HpODE -> hexanal + furan-bearing C12 acid | alkoxyl cyclisation onto C12, O2 at C13, scission C12-C13 (Fig 2: 7 -> 8 -> 9 -> 4 + 5) | hexanal area 67.5 M (9-EZ), 89.0 M (9-EE) | mechanism_drawn (proposed); measured_level |
| LA-9-E | 9-HpODE -> 2-octenal + 9,10-epoxy C10 fragment | epoxy-allyl radical, 11-OOH, scission C11-C12 (Fig 2: 7 -> 10 -> 11 -> 12 + 13) | 2-octenal areas 51.7 M + 36.4 M (9-EZ), 3.8 M + 126.9 M (9-EE) | mechanism_drawn; measured_level |
| LA-9-E2 | 9-HpODE -> 4,5-epoxy-2-decenal (via 13-HpODE isomerisation, 34) and -> pentanal/1-pentanol (via 15 -> 17) | Fig 2 | 4,5-epoxy-2-decenal 4.9 M + 12.3 M (9-EZ); pentanal 7.8 M; 1-pentanol 4.3 M | mechanism_drawn; measured_level |
| LA-13-B | 13-HpODE -> pentyl radical (-> pentanal, 1-pentanol) + 13-oxo-9,11-tridecadienoic acid | beta-scission C13-C14 (Fig 2: 24 -> 20 + 17 -> 18 -> 19/21) | pentanal 9.6 M (13-ZE), 3.6 M (13-EE); 1-pentanol 7.8 M / 4.2 M; pentane not reported (SPME/DB-WAX, m/z >= 41: pentane would not be seen reliably) | mechanism_drawn; measured_level |
| LA-13-F | **13-HpODE -> 2-pentylfuran + 9-oxononanoic acid** | alkoxyl cyclisation onto C10, O2 at C9, scission C9-C10 (Fig 2: 24 -> 25 -> 26 -> 22 + 23) | 2-pentylfuran area 2.55 M (13-ZE), 3.42 M (13-EE); also 2.57 M / 3.73 M from 9-EZ / 9-EE and 1.68 M from 10-HpODE; **none from 12-HpODE** | mechanism_drawn (proposed, novel); measured_level |
| LA-13-H | 13-HpODE -> hexanal + 12-oxo-dodecadienoic acid | epoxy-allyl radical 27 -> oxirane C-C cleavage -> vinyl-ether hydroperoxide 33 -> hexanal (Fig 2); authors reject direct C12-C13 alpha-scission ("not chemically facile since an energetically unfavorable vinyl radical is generated") | hexanal area 139.1 M (13-ZE), 57.4 M (13-EE); hexanoic acid 29.7 M / 26.4 M | mechanism_drawn (proposed); measured_level |
| LA-13-E | 13-HpODE -> 4,5-epoxy-2-decenal + octanoic acid | 27 -> 31 -> 34 -> 35 + 2 (Fig 2) | 4,5-epoxy-2-decenal 9.3 M (13-ZE, 45.1 min peak; 44.7 min entry unusable), 3.98 M + 10.8 M (13-EE) | mechanism_drawn; measured_level |
| LA-10-B | **10-HpODE -> 10-oxo-8-decenoic acid + 2-octenyl radical -> 1-octen-3-ol, 1-octen-3-one, 2-octenal, 2-octen-1-ol** | beta-scission C10-C11; allyl radical 40 <-> 41; O2; hydroperoxide reduction (Fig 3) | 1-octen-3-ol 144.8 M; 1-octen-3-one 22.6 M; 2-octen-1-ol 15.9 M; 2-octenal 6.1 M + 11.7 M; 2-propenal 14.2 M; pentanal 15.2 M; 1-pentanol 11.5 M | mechanism_drawn; measured_level |
| LA-12-A | **12-HpODE -> 2-heptenal + 9-undecenoic acid radical (-> octanoic acid + acrolein)** | alpha-scission C11-C12 (Fig 3: 47 -> 48 + 49) | 2-heptenal 265.0 M (12-HpODE) vs 3.6-10.4 M in every other column; octanoic acid 49.0 M; 2-propenal 3.8 M | mechanism_drawn; measured_level; within-study ratio: 2-heptenal/(all listed 12-HpODE areas) ~ 0.80 |
| LA-ISO | 9-HpODE <-> 13-HpODE; 10-HpODE -> 8-HpODE; 12-HpODE -> 14-HpODE | •OH H-abstraction from ROOH, O2 loss/re-addition on the pentadienyl/allyl radical (Fig 1b xv, xvi) | qualitative (LC-MS/MS product-ion chromatograms; GC: 2,4-decadienal 2.7 M from 13-ZE, 2.4 M from 13-EE) | figure_only (LC-MS/MS); measured_level (GC) |
| LA-3NON | 9-/10-HpODE -> 3-nonenal | NOT observed: "3-nonenal was not detected" in any column; authors take this as evidence against 9-HpODE beta-scission and against Hock fragmentation | null result | measured_level (absence) |

Within-study ratios that a later wave could register (same isomer column, same run):
- 13-ZE-HpODE: hexanal : 2-pentylfuran : 2,4-decadienal(2E,4Z) : pentanal = 139.1 : 2.55 : 2.73 : 9.61 M.
- 9-EZ-HpODE: hexanal : 2,4-decadienal(sum) : 2-octenal(sum) : 2-pentylfuran = 67.5 : 16.1 : 88.1 : 2.57 M.
- 10-HpODE: 1-octen-3-ol : 1-octen-3-one : 2-octen-1-ol : 2-octenal(sum) : hexanal = 144.8 : 22.6 : 15.9 : 17.8 : 28.5 M.
These are SPME peak-area ratios with no response-factor correction; treat as order-of-magnitude.

## 5. Rule sketches (repository suggestions, not the paper's)

Compound keys: `hexanal`, `2_pentylfuran`, `1_octen_3_ol`, `e_2_octenal`, `heptanal`, `acrolein`
exist in `data/keys/compounds.yml`; **octanal, decanal, (E)-2-decenal, (E,E)-2,4-decadienal,
2-heptenal, pentanal, 1-pentanol, 1-octen-3-one, pentane do not** (pentane and 2,4-decadienal exist
only as `PENTANE` / `DECADIENAL` in `data/species/structures.yml`). Structures: `LOOH_9_ct`,
`LOOH_9_tt`, `LOOH_13_ct`, `LOOH_13_tt` are METHYL ESTERS; no 10-/12-/8-/14-HpODE entries, no free
linoleic acid entry. Miyazaki's substrates are free acids; the ester/acid difference does not touch
the scissions below.

Convention note: the repo's R18a "side A" (alkane + 2,4-dienal) for LOOH_9 is Miyazaki's
"alpha-scission C8-C9"; R18b "side B" (saturated aldehyde) for LOOH_13 is Miyazaki's "alpha-scission
C12-C13" — the step they call unfavourable. Their measured product slate still contains all of R18a/b's
products, so R18a/R18b stand as NET rules; the mechanism anchor should say "net; Miyazaki 2023
attributes hexanal from 13-HpODE to the vinyl-ether-hydroperoxide route and hexanal from 9-HpODE to
the furyl route".

**S1. Furyl route to 2-pentylfuran (new rule, net).** Reactant: 13-hydroperoxy-9Z,11E-octadecadienoate.
Change: O of the 13-OOH becomes the ring oxygen bonded to C10 and C13; the C9-C10 bond breaks; C9
becomes an aldehyde carbon; the C10-C13 ring aromatises to furan (net loss of H2O and one H). Products:
2-pentylfuran (C10-C18) + 9-oxononanoate (C1-C9).
- positive control: `CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)OC` (LOOH_13_ct) -> `CCCCCc1ccco1` + `O=CCCCCCCCC(=O)OC` (ME_9_OXONONANOATE)
- second positive (free acid): `CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)O` -> `CCCCCc1ccco1` + `O=CCCCCCCCC(=O)O`
- negative control: methyl 9-hydroperoxy-10E-octadecenoate (oleate 9-OOH, no second C=C, cannot close a furan): `CCCCCCC/C=C/C(OO)CCCCCCCC(=O)OC` -> no fire; also stearic acid `CCCCCCCCCCCCCCCCCC(=O)O` -> no fire.
- required substructure: C(OOH)-CH=CH-CH=CH- (hydroperoxide alpha to a conjugated diene, gamma carbon sp2). 12-HpODE must NOT fire (measured: no 2-pentylfuran); 12-HpODE has OOH at C12 flanked by C13=C14 on one side and an isolated C9=C10 two carbons away — the SMIRKS should demand the conjugated diene on the side of the gamma carbon.

**S2. 9-HpODE furyl route to hexanal (new, net).** Reactant: 9-hydroperoxy-10E,12Z-octadecadienoate.
Change: 9-O bridges C9 and C12; C12-C13 breaks; C13 becomes CHO. Products: hexanal (C13-C18) +
8-(furan-2-yl)octanoate.
- positive: `CCCCC/C=C\C=C\C(OO)CCCCCCCC(=O)OC` (LOOH_9_ct) -> `CCCCCC=O` (HEXANAL) + `COC(=O)CCCCCCCc1ccco1`
- negative: oleate 9-OOH `CCCCCCC/C=C/C(OO)CCCCCCCC(=O)OC` (no second C=C at the gamma position, cannot close the ring) -> no fire; stearic acid -> no fire. Note that LOOH_13_ct fed to the mirror-image pattern is S1 (2-pentylfuran), so this rule and S1 are one SMIRKS shape applied at the two ends of the diene; hexanal from 13-HpODE is S3, a different transformation.

**S3. 13-HpODE -> hexanal (existing R18b net product; mechanism re-annotated).** Keep R18b's positive
control (LOOH_13_ct -> HEXANAL + oxo-alkenoate); add anchor to this dossier, route LA-13-H.

**S4. 10-HpODE -> C8 slate (new, net, two products).** Reactant: 10-hydroperoxy-8E,12Z-octadecadienoate.
Change: C10-C11 breaks; C10 becomes CHO (10-oxo-8-decenoate); C11-C18 leaves as the 2-octenyl radical,
which after O2 and reduction is written as 1-octen-3-ol (or 1-octen-3-one, 2-octenal, 2-octen-1-ol).
- positive: `CCCCC/C=C\CC(OO)/C=C/CCCCCCC(=O)O` -> `O=C/C=C/CCCCCCC(=O)O` + `C=CC(O)CCCCC` (1-octen-3-ol; key `1_octen_3_ol`)
- negative: LOOH_13_ct (conjugated 13-OOH cannot give a C8 allyl radical by this scission) -> no fire; stearic acid -> no fire.
- required substructure: CH=CH-CH(OOH)-CH2-CH=CH (non-conjugated, bis-allylic-flanked hydroperoxide).

**S5. 12-HpODE -> 2-heptenal (new, net).** Reactant: 12-hydroperoxy-9Z,13E-octadecadienoate. Change:
C11-C12 breaks; C12 becomes CHO; product 2-heptenal (C12-C18) + 11-carbon acid radical (-> octanoic
acid + acrolein per Fig 3, or 9-undecenoate).
- positive: `CCCC/C=C/C(OO)C/C=C\CCCCCCCC(=O)O` -> `CCCC/C=C/C=O` (2-heptenal; no key) + C11 fragment
- negative: LOOH_9_ct -> must not give 2-heptenal (measured 3.7 M vs 265 M, i.e. only via isomerisation).

**S6. Hydroperoxide -> alcohol (reduction), covers 1-octen-3-ol and 1-pentanol here.** The paper's
alcohols come from R-OOH -> R-OH (Fig 3: 42 -> 43, 44 -> 46, 18 -> 21), NOT from aldehyde reduction. No
paper in this batch draws aldehyde -> alcohol (hexanal -> 1-hexanol). If the repo wants 1-octen-3-ol it
should be S4 (hydroperoxide reduction), not a hexanal-type reduction.

## 6. Flags

1. **n = 1, peak areas, no internal standard, no calibration, NIST-only identification.** Nothing in
   Table 1 is a concentration. Cross-isomer comparisons assume identical SPME uptake, which is not
   checked. Use only as presence/absence and within-column order-of-magnitude ratios.
2. **The 13-ZE 4,5-epoxy-2-decenal entry "3 662 424 932" is a misprint** (see Table 1 note a).
3. **Neat film at 120 C for 5 min, ambient air in a 2 mL sealed vial**: not a food matrix, no water, no
   metals; the relative importance of cyclisation vs scission may differ in oil. Authors chose 120 C
   to keep intermediates alive.
4. **All mechanisms are "predicted pathways"** built before the experiment and supported by product
   presence and by [M+Na]+ m/z 351 product ions; no intermediate (furyl-hydroperoxide, vinyl-ether
   hydroperoxide) was isolated or NMR-confirmed ("this study could not elucidate the structure of m/z
   351 derived from 10-HpODE and 12-HpODE"). Evidence class for every mechanism row is
   mechanism_drawn, not established.
5. **9 <-> 13 isomerisation during heating** means any single-isomer product slate is contaminated
   by its partner's slate (2,4-decadienal and 2-octenal from 13-HpODE; 2-pentylfuran from 9-HpODE).
   A rule set that keys products to isomers should carry the isomerisation step too (LA-ISO).
6. **Vinyl hexanoate** (major from 13-HpODE, 50 M) has no assigned pathway ("its generation pathway
   was not determined").
7. **No DFT or other computed numbers** in the paper; nothing to mark inadmissible.
8. **Pentane** (the repo's R18a product from LOOH_13) is not reported; the DB-WAX/SPME method starting
   at m/z 41 with a 30 C hold would show pentane poorly; its absence in Table 1 is not evidence against
   R18a. Pentanal and 1-pentanol (the oxidised pentyl radical) are reported instead.
9. Substrate purity > 95 % by LC-MS; the remaining < 5 % is unstated and could include the partner
   isomer.
