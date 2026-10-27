# Hofmann & Schieberle 1997 — EXTRACTION
### *"Identification of Potent Aroma Compounds in Thermally Treated Mixtures of Glucose/Cysteine and Rhamnose/Cysteine Using Aroma Extract Dilution Techniques"*
### Glucose or rhamnose 100 mmol/L + cysteine 33 mmol/L, 0.5 M phosphate pH 5.0, 100 mL, autoclave ramped 20 -> 145 °C in 20 min. AEDA, static-headspace olfactometry and seven odour thresholds. **No concentration, no yield, no time course: every number is an FD factor, a headspace volume, an RI or a threshold.** Read 2026-10-27.

**Full read of `data/articles/hofmann1997.pdf`** (9 pp., Acrobat Distiller text layer, clean). Text
from `pdftotext -layout`; Table 1 re-checked against a 300-dpi raster of p. 899 (the bracketed
24/25 row and the column placement of every FD value agree with the text layer). Tables 2-7 came
out of the text layer intact. Figures 1-9 are mass spectra and a synthesis route; Figure 10 (odour
profile spider plots) was not digitised (sensory, no use to the engine).

**Provenance:** **[M]** printed number · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

**Benchmark verdict, up front: NO new benchmark-quality row.** The pots are the same charge as
Hofmann & Schieberle 1998 (JAFC 46, 235; `hofmann1998_reconciliation.md` Table 1: cysteine 3.3 mmol
+ carbohydrate 10.0 mmol in 100 mL, pH 5.0, autoclave), which already *quantifies* FFT and MFT for
exactly these two sugars (glucose FFT 2.8 / MFT 1.9 µg; rhamnose FFT 0.8 / MFT 0.8 µg per 100 mL
pot). This paper is the qualitative companion: which odorants form, not how much.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Thomas Hofmann, Peter Schieberle* |
| Title | "Identification of Potent Aroma Compounds in Thermally Treated Mixtures of Glucose/Cysteine and Rhamnose/Cysteine Using Aroma Extract Dilution Techniques" |
| Venue | *J. Agric. Food Chem.* **1997**, 45 (3), **898-906** |
| DOI | 10.1021/jf960456t (printed as the manuscript code "JF960456T" on p. 906 and "S0021-8561(96)00456-6" on p. 898; DOI string confirmed on Crossref 2026-10-27: JAFC 45(3) 898-906, 1997-03-01) |
| Affiliation | Deutsche Forschungsanstalt für Lebensmittelchemie, Garching |
| Dates | received 24 Jun 1996; revised 18 Dec 1996; accepted 19 Dec 1996; ASAP 15 Feb 1997 |
| Companions | `hofmann1995_extraction.md` (ribose + cysteine, same protocol, AEDA only); `hofmann1998_reconciliation.md` (the quantitative SIDA paper on the same sugar set); `hofmann1996_extraction.md` (thiol oxidative stability) |

---

## 1. THE ONE-PARAGRAPH ANSWER

Heated glucose/cysteine smells sulfury, roasty and caramel-like ("somewhat resembling roast
chicken"). By AEDA its top odorants are **2-furfurylthiol (FD 1024)** and **5-acetyl-2,3-dihydro-
1,4-thiazine (FD 1024)**, then 3-mercapto-2-butanone, 3-mercapto-2-pentanone and Furaneol (FD 512
each) and 2-(1-mercaptoethyl)furan (FD 256, glucose-specific). **MFT is absent from the glucose AEDA
(FD < 1) but is detected by static-headspace olfactometry at FD 8** (0.63 mL headspace) — the text
says MFT "has not yet been reported to be formed from cysteine in the presence of hexoses". Rhamnose
shifts the profile to Furaneol (FD 65536) and 3-hydroxy-6-methyl-2(2H)-pyranone (FD 16384), and
makes 5-methyl-2-furfurylthiol (FD 2048, rhamnose-specific). FFT, the thiazine and 2-acetyl-2-
thiazoline have similar FD factors from ribose, glucose and rhamnose (Table 6). One fed side-
experiment: **5-methylfurfural + H2S (1 mmol each), same conditions, gives 5-methyl-2-furfurylthiol
in 0.4 % yield, rising with H2S** (data not shown). New odorant: 2-propionyl-2-thiazoline,
threshold 0.07 ng/L air.

---

## 2. METHODS [M] (pp. 898-901)

| item | as printed |
|---|---|
| Model reaction | "L-Cysteine (**3.3 mmol**) and either D-glucose (**10 mmol**) or L-rhamnose (**10 mmol**) were dissolved in phosphate buffer (**100 mL; 0.5 M; pH 5.0**) and allowed to react in a laboratory autoclave (model II, **200 mL total volume**; Roth) by **raising the temperature within 20 min from 20 to 145 °C**." [Z: Cys 33 mmol/L, sugar 100 mmol/L, 1 : 3] |
| [!] cysteine charge | Table 1 footnote a says **3 mmol**; Methods and Table 3 footnote say **3.3 mmol**. Take 3.3 (matches Hofmann 1998). |
| [!] thermal history | Abstract and Results say "20 min, 145 °C"; Methods says a **ramp 20 -> 145 °C over 20 min**. No hold time, no cooling time. Hofmann 1998 prints "20 min at 145 °C"; the same ambiguity sits under every pot of this protocol family (see `hofmann1998_reconciliation.md`). |
| Headspace | autoclave 200 mL total, 100 mL liquid: ~100 mL headspace, closed |
| Isolation | diethyl ether extraction + sublimation in vacuo (high-vacuum distillation); split into **acidic (fraction I)** and **neutral/basic (fraction II)** by NaHCO3 |
| AEDA | HRGC-O on FFAP and SE-54, cold on-column; two assessors in duplicate, "data differed by fewer than two FD factors" |
| Scale-up for ID | 25 pots pooled (2.5 L) for identification |
| SHO | static headspace olfactometry: 20 mL maximum headspace, cryofocused; FD = 20 mL / lowest volume perceived |
| Thresholds | HRGC-O in air (Ullrich & Grosch 1987), reference (E)-2-decenal 2.7 ng/L air |
| Quantification | **none** |

---

## 3. TABLES (verbatim)

### Table 1. "Intense Odorants (FD ≥ 4) Formed by Reacting Glucose and Cysteine in Aqueous Solution" (FD in acidic fraction I / neutral-basic fraction II)

| no. | odorant | RI FFAP | RI SE-54 | FD I | FD II | lit. ref. |
|---|---|---|---|---|---|---|
| 1 | 3-mercapto-2-butanone | 1282 | 820 | 512 | <1 | – |
| 2 | ethyl 2-mercaptopropionate | 1300 | 869 | 8 | <1 | – |
| 3 | 3-mercapto-2-pentanone (e) | 1347 | 907 | 512 | <1 | – |
| 4 | mercapto-2-propanone | 1357 | <800 | 8 | <1 | 1 |
| 5 | 2-(1-mercaptoethyl)furan | 1396 | 951 | 256 | <1 | – |
| 6 | 2-furfurylthiol | 1431 | 909 | 1024 | <1 | 2,3,4 |
| 7 | 2,3-diethyl-5-methylpyrazine (f) | 1479 | 1163 | <1 | 32 | – |
| 9 | 2-methyltetrahydrothiophen-3-one | 1512 | 1017 | 8 | <1 | 5,6,7 |
| 10 | 2-acetylthiazol | 1615 | 1020 | <1 | 32 | 1,7,8 |
| 11 | 2-(1-mercaptoethyl)thiophene | 1655 | 1136 | 16 | <1 | – |
| 12 | 2-formylthiophene | 1669 | 1000 | <1 | 128 | 1,2,6,7,9 |
| 13 | 2-thenylmercaptan | 1682 | 1092 | 32 | <1 | 2,3 |
| 14 | 2-propionylthiazol | 1700 | 1107 | <1 | 64 | – |
| 16 | 2-acetyl-2-thiazoline | 1720 | 1111 | <1 | 128 | 8 |
| 17 | 2-acetylthiophene | 1740 | 1090 | <1 | 128 | 1,3,6,7,9 |
| 18 | 5-methyl-2-formylthiophene | 1744 | 1125 | <1 | 256 | 1,2,6,7,9 |
| 20 | unknown (earthy, roasty) | 1784 | – | <1 | 128 | – |
| 22 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 2016 | 1100 | 512 | <1 | 1,4 |
| 24 | 4-hydroxy-2,5-dimethyl-3(2H)-thiophenone | 2040 | 1181 | 128 (24+25 coelute) | <1 | 1 |
| 25 | 2-ethyl-4-hydroxy-5-methyl-3(2H)-furanone (f) | 2050 | 1178 | (bracketed with 24) | <1 | – |
| 26 | 2-mercaptopropionic acid | 2130 | 1057 | 16 | <1 | 1 |
| 27 | 3-hydroxy-4,5-dimethyl-2(5H)-furanone (Sotolon) (f) | 2153 | 1112 | 32 | <1 | (blank) |
| 28 | 5-acetyl-2,3-dihydro-1,4-thiazine | 2177 | 1374 | <1 | 1024 | – |
| 29 | 3-mercaptopropionic acid | 2193 | 1057 | 4 | <1 | 1 |
| 30 | 5-propionyl-2,3-dihydro-1,4-thiazine | 2235 | 1456 | <1 | 8 | – |
| 31 | bis(1-(2′-furyl)-1-ethyl) disulfide | 2422 | 1732 | 8 | 8 | – |
| 32 | (2-furfuryl)-(1-(2′-furyl)-1-ethyl) disulfide | 2449 | 1696 | 8 | 8 | – |
| 33 | bis(2-furfuryl) disulfide | 2465 | 1673 | 8 | 8 | – |

Footnotes (condensed): (a) glucose 10 mmol + cysteine **3 mmol** [!], 20 min, phosphate 100 mL, 0.5 M,
pH 5.0. (e) contained ~25 % of its isomer 2-mercapto-3-pentanone. (f) MS too weak; identified on RI,
odour quality and intensity only. Literature refs: (1) Martin 1988; (2) Shibamoto & Russel 1976; (3)
1977; (4) Yeo & Shibamoto 1991; (5) Kato 1973; (6) Sakaguchi & Shibamoto 1978a; (7) Zhang & Ho 1991;
(8) Sakaguchi & Shibamoto 1978b; (9) Scanlan 1973. The printed title reads "FD = 4"; the glyph is ≥
(Table 3 prints "FD g 4", the font's ≥).

### Table 2. Static headspace olfactometry, glucose/cysteine (vol = lowest headspace volume perceived; FD = 20 mL / vol)

| odorant | RI FFAP | vol (mL) | FD |
|---|---|---|---|
| hydrogen sulfide | <900 | 5.0 | 1 |
| methylmercaptan | <900 | 1.25 | 4 |
| ethylmercaptan | <900 | 0.31 | 16 |
| 3-mercapto-2-butanone | 1282 | 2.5 | 2 |
| **2-methyl-3-furanthiol** | 1300 | **0.63** | **8** |
| 3-mercapto-2-pentanone (d) | 1347 | 1.25 | 4 |
| 2-(1-mercaptoethyl)furan | 1396 | 0.31 | 16 |
| 2-furfurylthiol | 1431 | 0.08 | 64 |
| 2-(1-mercaptoethyl)thiophene (e) | 1655 | 2.5 | 2 |
| 2-thenylthiol | 1682 | 1.25 | 4 |

[!] FD 20/vol does not hold exactly for every row as printed (20/0.63 = 32, not 8; 20/0.31 = 65,
not 16; 20/0.08 = 250, not 64). The printed FDs are consistent with vol × 4 being a 5 mL basis
(text: "In a 5 mL sample of air above the model, 10 odor-active regions were detectable"), i.e.
FD = 5 / vol. Use the volumes, not the FD column, if anything quantitative is ever drawn from this.

### Table 3. "Intense Odorants (FD ≥ 4) Formed by Reacting Rhamnose with Cysteine in Aqueous Solution" (rhamnose 10 mmol + cysteine 3.3 mmol, 100 mL 0.5 M phosphate pH 5.0, 20 min)

| no. | odorant | RI FFAP | RI SE-54 | FD I | FD II | lit. ref. |
|---|---|---|---|---|---|---|
| 1 | 3-mercapto-2-butanone | 1282 | 820 | 32 | <1 | – |
| 3 | 3-mercapto-2-pentanone | 1347 | 907 | 128 | <1 | – |
| 6 | 2-furfurylthiol | 1431 | 909 | 512 | <1 | – |
| 7 | 2,3-diethyl-5-methylpyrazine | 1479 | 1163 | <1 | 64 | – |
| 8 | 5-methyl-2-furfurylthiol | 1500 | 995 | 2048 | <1 | 1,2 |
| 15 | unknown (sulfury) | 1707 | – | 32 | <1 | – |
| 16 | 2-acetyl-2-thiazoline | 1720 | 1111 | <1 | 256 | – |
| 17 | 2-acetylthiophene | 1740 | 1090 | <1 | 32 | 1,2,3 |
| 18 | 5-methyl-2-formylthiophene | 1744 | 1125 | <1 | 256 | – |
| 19 | 5-methyl-2-thenylthiol | 1755 | 1179 | 32 | <1 | – |
| 20 | unknown (earthy, roasty) | 1784 | – | <1 | 64 | – |
| 21 | 2-propionyl-2-thiazoline | 1803 | 1196 | <1 | 32 | – |
| 22 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 2016 | 1100 | 65536 | <1 | 1,2 |
| 23 | 3-hydroxy-6-methyl-2(2H)-pyranone | 2035 | 1047 | 16384 | <1 | 1,2 |
| 25 | 2-ethyl-4-hydroxy-5-methyl-3(2H)-furanone | 2050 | 1178 | 64 | <1 | – |
| 27 | 3-hydroxy-4,5-dimethyl-2(5H)-furanone | 2153 | 1112 | 128 | <1 | – |
| 28 | 5-acetyl-2,3-dihydro-1,4-thiazine | 2177 | 1374 | <1 | 512 | – |
| 30 | 5-propionyl-2,3-dihydro-1,4-thiazine | 2235 | 1456 | <1 | 32 | – |

Refs: (1) Martin 1988; (2) Silwar 1992; (3) Yamaguchi 1979.

### Table 4. Static headspace olfactometry, rhamnose/cysteine

| odorant | RI FFAP | vol (mL) | FD |
|---|---|---|---|
| hydrogen sulfide | <900 | 2.5 | 2 |
| methylmercaptan | <900 | 2.5 | 2 |
| ethylmercaptan | <900 | 0.3 | 16 |
| unknown (sulfury) | 1000 | 5.0 | 1 |
| 2-methyl-3-furanthiol | 1300 | 5.0 | 1 |
| 3-mercapto-2-pentanone (d) | 1347 | 2.5 | 2 |
| dimethyl trisulfide | 1356 | 5.0 | 1 |
| 2-furfurylthiol | 1431 | 0.3 | 16 |
| 5-methyl-2-furfurylthiol | 1500 | 0.04 | 125 |
| 2-thenylmercaptan | 1682 | 5.0 | 1 |
| 5-methyl-2-thenylthiol | 1755 | 1.25 | 4 |

(Same FD = 5 / vol pattern as Table 2.)

### Table 5. "Odor Thresholds of Selected Sulfur-Containing Odorants in the Cysteine/Carbohydrate Mixtures" (HRGC-O, ng/L air)

| odorant | threshold (ng/L air) |
|---|---|
| 2-(1-mercaptoethyl)furan | 0.009 |
| 5-methyl-2-furfurylthiol | 0.006 |
| 2-(1-mercaptoethyl)thiophene | 0.016 |
| 5-methyl-2-thenylthiol | 0.019 |
| 3-hydroxy-6-methyl-2(2H)-pyranone | 0.76 |
| 4-hydroxy-2,5-dimethyl-3(2H)-thiophenone | 1.6 |
| bis(2-(1-mercaptoethyl)furane) disulfide | 0.005 |

### Table 6. "Comparison of the Key Odorants (FD Factors >128 to 65536) Formed by Reacting Cysteine in the Presence of Carbohydrates"

| no. | odorant | glucose | rhamnose | ribose (a) |
|---|---|---|---|---|
| 6 | 2-furfurylthiol | 1024 | 512 | 4096 |
| 28 | 5-acetyl-2,3-dihydro-1,4-thiazine | 1024 | 512 | 1024 |
| 8 | 5-methyl-2-furfurylthiol | <1 | 2048 | <1 |
| 3 | 3-mercapto-2-pentanone | 512 | 128 | 2048 |
| 1 | 3-mercapto-2-butanone | 512 | 32 | 512 |
| 22 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 512 | 65536 | 128 |
| 5 | 2-(1-mercaptoethyl)furan | 256 | <1 | <1 |
| 16 | 2-acetyl-2-thiazoline | 128 | 256 | 256 |
| 23 | 3-hydroxy-6-methyl-2(2H)-pyranone | <1 | 16384 | <1 |
| 27 | 3-hydroxy-4,5-dimethyl-2(5H)-furanone | 16 | 128 | 64 |
| A | 2-methyl-3-furanthiol | <1 | <1 | 1024 |
| B | 4-hydroxy-5-methyl-3(2H)-furanone | <1 | <1 | 256 |
| 24 | 4-hydroxy-2,5-dimethyl-3(2H)-thiophenone | 128 (b) | <1 | <1 |

(a) adapted from Hofmann & Schieberle 1995. (b) coelution with ethylfuraneol.

### Table 7. Odour thresholds of roast-smelling odorants (ng/L air; literature compilation)

| odorant | threshold | identified in |
|---|---|---|
| 2-acetyl-1-pyrroline | 0.02 | bread crust; popcorn; proline/carbohydrates; ornithine/carbohydrates |
| 2-propionyl-1-pyrroline | 0.02 | popcorn |
| 6-acetyl-tetrahydropyridine | 0.06 | bread crust; popcorn |
| 5-acetyl-2,3-dihydro-1,4-thiazine | 0.05 | cysteine/ribose |
| 5-propionyl-2,3-dihydro-1,4-thiazine | 0.10 | cysteine/ribose |
| 2-acetylthiazol | 4.10 | cooked meat |
| 2-propionyl-2-thiazoline | 0.07 | cysteine/rhamnose |
| 2-acetyl-2-thiazoline | 0.05 | roast meat; cysteine/ribose |
| 2-acetylpyrazine | 0.4 | wheat bread crust; popcorn |

### Text-only numbers [M]

| quantity | value | where |
|---|---|---|
| 5-methylfurfural + H2S (1 mmol each) -> 5-methyl-2-furfurylthiol | **0.4 %** yield, "significantly dependent on the H2S concentration" | p. 905, "data not shown"; "same reaction conditions as in the carbohydrate/cysteine model" (pH 5.0, 0.5 M phosphate, 100 mL -> 10 + 10 mmol/L, 145 °C ramp) [Z] |
| 2-acetylfuran -> 2-(1-mercaptoethyl)furan | "to the same extent" as above (≈ 0.4 %) | p. 905, data not shown |
| synthesis yields (MFFT 46 %, MEF 52 %, 5-MTT 49 %, MET 79 %, 2-propionyl-2-thiazoline 5 % from cysteamine + 2-oxobutanal 1 h 80 °C pH 7, HMP 23 %, HDMT 6 %) | — | Experimental; preparative, not kinetic |

---

## 4. OVERLAP WITH THE EXISTING HOFMANN DOSSIERS

| this paper | already in the corpus | relation |
|---|---|---|
| glucose + Cys pot, 33 / 100 mmol/L, pH 5, 145 °C | `hofmann1998_reconciliation.md` Table 1 / Table 2 (glucose FFT 2.8 µg, MFT 1.9 µg per 100 mL at pH 5; pH 3 and 7 also) | **same charge, same autoclave**; 1998 quantifies, 1997 only screens |
| rhamnose + Cys pot | same 1998 tables (rhamnose FFT 0.8, MFT 0.8 µg) | same |
| mercapto-2-propanone FD 8 (glucose) | 1998 Table 9 quantifies mercapto-2-propanone from glucose + Cys | 1998 is the number |
| ribose column of Table 6 | `hofmann1995_extraction.md` | re-print of 1995 FD factors |
| thiol stability in work-up | `hofmann1996_extraction.md` | same caveat on disulfides (below) |

So nothing here is new benchmark material. The only quantitative statements are the fed
5-methylfurfural + H2S yield (0.4 %, data not shown) and the thresholds.

---

## 5. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor (qualitative / ordinal only):**
1. **Hexose MFT is real but minor.** MFT is below AEDA detection in the glucose extract (FD < 1)
   yet perceptible in the headspace (Table 2, 0.63 mL). That is consistent with Hofmann 1998's
   glucose MFT of 1.9 µg per pot and with the engine giving glucose pots a small, non-zero MFT. It
   also corrects the wording in `data/benchmarks/quarantined/README.md` ("MFT does not appear
   among them"): MFT does appear, by SHO, at low potency. The quarantine verdict itself (the
   fabricated Farmer1999 ppb values) is unaffected.
2. **Sugar-ordering of FFT by FD:** ribose 4096 > glucose 1024 > rhamnose 512. The same order as
   Hofmann 1998's µg (12.1 > 2.8 > 0.8), so FD tracks amount here. Ordinal check only.
3. **Rhamnose -> 5-methylfurfural -> 5-methyl-2-furfurylthiol** (FD 2048) and **glucose -> 2-acetyl-
   furan -> 2-(1-mercaptoethyl)furan** (FD 256): two sugar-specific furfurylthiol analogues the engine
   does not carry. The fed 0.4 % yield (5-methylfurfural + H2S, 10 + 10 mmol/L) is the analogue of
   an FFT-from-furfural + H2S step; if the engine has a furfural + H2S -> FFT step, 0.4 % at 145 °C
   ramp, pH 5 is an order-of-magnitude cross-check for the methyl homologue (not the parent).
   Data not shown; single point; no time.
4. **Glucose gives the C3/C4 mercaptoketones strongly** (3-mercapto-2-butanone, 3-mercapto-2-
   pentanone FD 512; mercapto-2-propanone FD 8): consistent with the engine's dicarbonyl + H2S
   chemistry existing for hexoses, but no amounts.

**Cannot anchor:**
- No concentration, no yield, no time series, no temperature or pH series. AEDA FD factors are
  extraction-and-olfactometry potencies, not amounts; they cannot be scored against model output.
- The SHO FD column is internally inconsistent with its stated definition (Table 2 flag above).
- The thermal history is a 20-min ramp to 145 °C with no stated hold; the same unresolved
  ambiguity as the Hofmann 1998 benchmarks.

**Caveats:**
- [!] Cysteine 3 vs 3.3 mmol between Table 1 footnote and Methods.
- [!] Disulfides (nos. 31-33) in a solvent extract after vacuum distillation can be work-up
  artefacts of the thiols (see `hofmann1996_extraction.md`); their FD 8 says nothing about pot
  disulfide content.
- 2,3-diethyl-5-methylpyrazine, Sotolon, ethylfuraneol identified without MS confirmation (f).
