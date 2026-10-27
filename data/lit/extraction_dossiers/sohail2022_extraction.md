# Sohail et al. 2022 — EXTRACTION
### *"Aroma compounds identified in cooked meat: A review"*
### A list of 332 compounds detected by GC–O in cooked meat. Each entry gives a NIST retention index, an odour descriptor, a van Gemert (2011) odour threshold and the meat species with reference numbers. There are no concentrations. The numbered primary references are in Supplementary file 1, which is not in the repo. Read 2026-10-27.

**Full read of `data/articles/sohail2022.pdf`** (14 pages, born-digital text layer). The rows below were
taken from the text layer and checked against 300-dpi rasters of pp. 4, 9 and 10. The rasters were
needed for the **bold reference numbers**, which mark the paper each odour descriptor was copied from;
the text layer does not keep bold.

**Provenance:** **[M]** printed in the review · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Amjad Sohail, Sam Al-Dalali, Jianan Wang, Jianchun Xie\*, Ashbala Shakoor, Sailimuhan Asimi, Haroon Shah, Prasanna Patil |
| Title | "Aroma compounds identified in cooked meat: A review" |
| Venue | *Food Research International* **157** (**2022**) **111385** (Review) |
| DOI | 10.1016/j.foodres.2022.111385 (printed) |
| Affiliations | Beijing Key Lab of Flavor Chemistry, Beijing Technology and Business University; Hefei University of Technology |
| Dates | received 23 Feb 2022, revised 27 Apr 2022, accepted 16 May 2022, online 19 May 2022 |
| Funding | NSFC 31671895; Beijing Municipal NSF 6172004 |
| PDF character | Table 1 (157 lipid-derived compounds, Nos. 1–157), Table 2 (Maillard and lipid–Maillard compounds, Nos. 158–273), Supplementary Table 1 (59 "other source" compounds) and **Supplementary file 1 (the numbered reference list behind every table entry), neither of which is in the repo**, Fig. 2 (FD/OAV of the key odorants in Fan et al. 2018), Fig. 4 (PLS-DA) |

---

## 1. THE ONE-PARAGRAPH ANSWER

The review lists **332 odorants** detected by GC–O in thermally cooked meat (beef, pork, poultry, sheep,
plus alligator, crab and snail) over roughly 40 years. They fall into 157 lipid-derived, 98 Maillard
(43 S-, 29 N-, 26 O-containing), 18 lipid–Maillard and 59 "other" compounds. **The tables hold no
concentrations.** Retention indices are taken from the NIST 11 database, not measured. **Odour
thresholds are all copied from one secondary compilation**, van Gemert (2011), in water (µg/kg) unless
marked. Odour descriptors are copied from one primary paper per row, marked in bold. Species
occurrence is given as reference numbers into an unavailable supplement. The only quantitative meat
data are OAV/FD values: **Fig. 2 (Fan et al. 2018, stewed Beijing-You chicken broth)** and two values
in the text (Gąsior 2021; Du 2020). The PLS-DA is run on how often each compound is reported in the
literature, not on amounts.

---

## 2. TARGET-COMPOUND ROWS [M]

Thresholds are in µg/kg water (van Gemert 2011, as quoted by the review). "Descriptor source" is the
bold reference number, i.e. Supplementary file 1 number and species; "–" means no bold number was
visible. The last column gives the model's current value from `data/species/*.yml`, for comparison
only. **Nothing was changed.**

| No. | compound | RI polar / non-polar | odour | threshold | reported in (Supp. ref. nos.) | descriptor source | model threshold |
|---|---|---|---|---|---|---|---|
| 160 | 2-methyl-3-furanthiol (MFT) | 1323 / 860 | roasted meat, cooked ham | **0.007** | beef 2,11,30,34,40,50–51,63,65,72; chicken 4,15,17; duck 19; sheep 20,26; pork 14 | beef [2] | 0.007 |
| 162 | 2-furfurylthiol (FFT) | 1439 / 915 | meaty, coffee | **0.12–0.13** | beef 11,30,32,40–41; chicken 4,15,17,31; duck 19,60,66; sheep 26; pork 14,32 | beef [11] | 0.01 [!] ×12 |
| 180 | bis(2-methyl-3-furyl) disulfide | 2150 / 1530 | meaty broth | **0.000032** | beef 3,32,51,72; chicken 15; pork 32 | chicken [15] | 0.02 [!] ×600 |
| 176 | 2-methyl-3-(methyldithio)furan | 1668 / 1169 | vitamin-like, meaty | 0.004 | beef 40 | beef [40] | — |
| 174 | 2-methyl-3-(methylthio)furan | 1337 / 938 | meaty | 25–30 | beef 3,40; chicken 15; sheep 7 | chicken [15] | — |
| 175 | methional | 1480 / 894 | cooked potato | **4.5** | beef 1–2,11,23,27–30,32,34,40–41,48,50,54,63,65,71–72; chicken 4,15–17,31,42; crab 38,55; duck 13,19,60,66; pork 32; sheep 20,25,35 | beef [11] | 0.2 [!] ×22 |
| 200 | 2-acetyl-2-thiazoline | 1752 / 1106 | roasty | **0.12** | beef 11,28–29,32,40–41,48,57,71; chicken 4,31; pork 32; sheep 25,35 | – | no record |
| 202 | "2-acetylthiazoline" | 1785 / 1090 | roasty | / | beef 50; chicken 16 | chicken [16] | — |
| 199 | 2-acetylthiazole | 1652 / 1025 | rice | 3 | beef 1–2,32,40,72; chicken 4,15,49; crab 18; duck 19; pork 5,10,14,32,61,69; sheep 26; snail 45 | pork [5] | — |
| 255 | furfural | 1460 / 835 | almond, woody | **9.56** | beef 2; chicken 15; crab 18; pork 5–6,10,39,69 | beef [2] | 3000 [!] ×300 |
| 40 | hexanal | 1078 / 802 | green, fresh grass | **5** | alligator 8,21; beef 1–3,11,22–23,27,29,32,40,43,59,63–65,71; chicken 15–17,24,31,33; crab 12,18,44; duck 13,19,60,66; pork 5–6,10,39,61,67–69; sheep 7,70 | beef [22] | 4.5 (Belitz/Czerny-sourced) |
| 213 | 2,5-dimethylpyrazine | 1323 / 915 | popcorn, roasted | **1.75** | beef 23,71; chicken 15,49; crab 12; pork 10; sheep 25–26; snail 45 | pork [10] | 1800 [!] ×1000 |
| 223 | 2-ethyl-3,5-dimethylpyrazine | 1464 / 1079 | earthy, musty | **0.04** | beef 3,11,23,28–29,40,48; chicken 15; duck 19; sheep 25–26 | beef [40] | 1.0 [!] ×25 |
| 211 | methylpyrazine | 1264 / 827 | roast | 30 | beef 71; chicken 15,42,49; crab 12; pork 10; sheep 70; snail 45 | chicken [49] | — |
| 212 | 2,3-dimethylpyrazine | 1346 / 930 | dry ham, cooked rice, peanut | 800 | beef 2; chicken 49 | beef [2] | — |
| 214 | 2,6-dimethylpyrazine | 1328 / 906 | stale, roasted | 720 | duck 19; sheep 25 | duck [19] | — |
| 217 | trimethylpyrazine | 1406 / 1005 | roasty, green, sweet | 350 | beef 23,29; chicken 15,49; crab 12 | beef [23] | — |
| 222 | 3-ethyl-2,5-dimethylpyrazine | 1447 / 1069 | nutty, roasted sunflower seeds | 8.6 | beef 23; chicken 49; pork 10 | pork [10] | — |
| 225 | tetramethylpyrazine | 1474 / 1072 | roast | 2.52 | chicken 49; sheep 7 | chicken [49] | — |
| 229 | 2,3-diethyl-5-methylpyrazine | 1461 / 1158 | earthy, musty | 0.0031 | beef 1,11,23,28–29,40–41,48; sheep 20 | not checked | — |
| 243 | 2-acetyl-1-pyrroline | 1327 / 926 | roasted, nutty | 0.0053 | beef 3,11,29,32,40,63,65,71; chicken 4,16; crab 38,55; sheep 25,35; pork 14,32 | not checked | 0.1 |

Related sulfur rows, for completeness: H2S 0.01–0.1 (beef 1); methanethiol 0.2 (beef 2,30,40,46,50;
chicken 15; pork 14,69,74); dimethyl disulfide 1.1; dimethyl trisulfide 0.1; 2-furfuryl methyl
disulfide 0.04 (beef 51); bis(2-furylmethyl) disulfide 0.00015–0.0006 **in air** (beef 51; chicken 31).

### 2.1 The only quantitative meat values (Fig. 2, from Fan et al. 2018; stewed Beijing-You chicken broth, SAFE–GC–MS/AEDA–GC–O, recombination and omission) [M]

| compound | FD | OAV |
|---|---|---|
| 2-methyl-3-furanthiol | 2048 | **2686** |
| 3-(methylthio)propanal (methional) | 4096 | **3953** |
| "2-furanthiol" (drawn as 2-furfurylthiol, CH2SH) [!] | 4096 | **543** |
| bis(2-methyl-3-furyl) disulfide | 128 | **810** |
| 2-methyl-3-(methylthio)furan | 1024 | 102 |
| dimethyl disulfide (label; structure drawn as R–S–S–R) | 1024 | 10 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 128 | 1 |
| 2-ethyl-3,5-dimethylpyrazine | 2048 | 1 |
| hexanal (Fig. 2a) | 1024 | 23 |

In the text: 2-acetyl-1-pyrroline, FD 1024 and OAV 1403, roasted goose (Gąsior et al. 2021);
2-ethyl-3,5-dimethylpyrazine, FD 2048 and OAV 15, grilled mutton shashlik (Du et al. 2020).

**[Z] Back-converting to concentrations is not safe from the review alone,** because Fan's own
thresholds are not given. With the review's 0.007 µg/kg, MFT ≈ 2686 × 0.007 ≈ **19 µg/kg broth**. With the
review's 4.5 µg/kg, methional would be ≈ 18 mg/kg, which is implausible, so Fan must have used a
different methional threshold. **Get the concentrations from Fan 2018 itself.**

---

## 3. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Verdict: none of it can serve as a prior or a score. It is useful as a pointer to primary
papers and as a cross-check of thresholds.**

- **Cannot anchor any concentration, yield or rate.** The tables have no amounts, and the only OAVs come
  from one chicken-broth paper that is not yet in the corpus.
- **Thresholds come from a secondary source.** Every value is van Gemert (2011), copied as one
  number or a range with no matrix or method, so none of them is a primary source for
  `odour_threshold_ug_per_kg`. Several disagree with the model by 12–1000×. These are flags to check
  against primary threshold papers, not corrections:
  - 2,5-dimethylpyrazine 1.75 against 1800. Probably a unit slip, mg/kg printed as µg/kg.
  - Furfural 9.56 against 3000. Aqueous furfural thresholds are usually in the mg/kg range, so 9.56 is
    likely a different medium or a unit slip.
  - Bis(2-methyl-3-furyl) disulfide 0.000032 against 0.02. The text says "ca. 10⁻¹³ g/g"; the table
    gives 3.2 × 10⁻¹⁴ g/g.
  - Methional 4.5 against 0.2.
  - FFT 0.12–0.13 against 0.01.
  - 2-Ethyl-3,5-dimethylpyrazine 0.04 against 1.0.
  - 2-Acetyl-1-pyrroline 0.0053 against 0.1.

  The MFT (0.007) and hexanal (5 against 4.5) values agree with the model.
- **Qualitative support.** MFT, FFT, methional and 2-acetyl-2-thiazoline are among the GC–O
  odorants reported most often across cooked beef, chicken, pork and sheep. This fits their status as
  "critical" model targets. Bis(2-methyl-3-furyl) disulfide appears among the high-VIP compounds that
  separate beef and pork in the PLS-DA. That analysis is a literature-frequency statistic and says
  nothing about chemistry.

**Caveats [!].**
1. **Duplicate rows:** No. 200 "2-acetyl-2-thiazoline" and No. 202 "2-acetylthiazoline" are almost
   certainly the same compound, listed with different NIST RIs. Phenylethyl alcohol (No. 249) and
   β-phenylethanol (No. 250) are another duplicate pair.
2. **Mislabels in Fig. 2(b):** "2-furanthiol" is drawn as 2-furfurylthiol. The text repeats "2-furanthiol
   … mostly reported", which almost certainly means FFT.
3. **Mixed media in the threshold column:** some entries are in air (µg/m³), in ethanol or in oil, with no
   flag beyond a trailing word. Table 1's header reads "µg/L or µg/m³" while its footnote says µg/kg.
4. **Lipid–Maillard assignments** (marked \*) are the authors' structural judgement from alkyl-chain
   length, not mechanistic evidence.
5. **The reference numbers in Tables 1 and 2 cannot be resolved without Supplementary file 1.** The
   primary papers that quantify thiols in meat (likely the Grosch-school SIDA work) are hidden there.
   Download the supplement from the article page before using any species/reference mapping.

---

## 4. PRIMARY PAPERS CITED IN THE MAIN TEXT — CORPUS CHECK

Each DOI is as printed in the review's reference list and was checked against `data/keys/papers.yml`
on 2026-10-27.

**Not in the corpus. Candidates for meaty-thiol or key-odorant quantification in real cooked meat:**

| reference | DOI | why |
|---|---|---|
| Fan et al. 2018, *JAFC* 66:10242 (Beijing-You and broiler chicken broths) | 10.1021/acs.jafc.8b03297 | **Top priority.** Its OAVs for MFT 2686, methional 3953, FFT 543 and MFT-SS 810 (Fig. 2) imply quantified concentrations in a stewed broth of stated recipe |
| Guth & Grosch 1994, *JAFC* 42:2862 (stewed beef juice, character-impact odorants) | 10.1021/jf00048a039 | Grosch-school stewed beef. Concentrations for MFT and methional are probably reported (from memory; not stated in the review) |
| Liu et al. 2019, *JAFC* 67:5847 (Beijing roasted duck; GC–O–MS, OAVs, recombination) | 10.1021/acs.jafc.9b01564 | OAVs require concentrations. Duck MFT and FFT are cited in Table 2 |
| Du et al. 2020, *Flavour Fragr. J.* 36:111 (grilled mutton shashlik ± suet) | 10.1002/ffj.3621 | OAVs; pyrazine OAV 15. Sheep MFT and FFT appear in Table 2 |
| Gąsior et al. 2021, *JAFC* 69:5986 (roasted White Kołuda goose) | 10.1021/acs.jafc.1c01475 | OAVs (2-acetyl-1-pyrroline 1403). Goose appears in Table 2 only for 2-phenylethanethiol, so lower priority for MFT and FFT |
| Qi et al. 2022, *Food Chem.* 371:131111 (stewed goat meat + thyme) | 10.1016/j.foodchem.2021.131111 | Key-odorant study with sensory verification; content for thiols unknown |
| Zhao et al. 2017a, *Food Chem.* 226:51 (black-pig pork broth) | 10.1016/j.foodchem.2017.01.011 | Pork broth volatiles; whether thiols are quantified is unknown |

**Not in the corpus. Model systems rather than meat (lower priority):**
- Wang et al. 2020, glucose–GSH–aldehydes: 10.1016/j.foodchem.2019.125482
- Zhang et al. 2021, phospholipid–xylose–cysteine: 10.1016/j.foodres.2020.109918
- Ma et al. 2020, xylose–cysteine + hydrolysed mutton tallow: 10.1002/ejlt.201900212
- Wang et al. 2021, 2-acetylthiazole pathways: 10.1016/j.foodchem.2020.128761
- Zhao et al. 2019b, GSH–glucose ± oxidised fat, *Food Chem.* 270:436 (no DOI printed)

**Already in the corpus:** Hou et al. 2017 (10.1016/j.foodchem.2017.03.133; `hou2017_extraction.md`)
and Mottram 1998 (10.1016/S0308-8146(98)00076-4).

**[!] One registry gap:** Zhao et al. 2019a (cysteine–xylose–glycine, 10.1016/j.foodchem.2018.08.096)
is on disk as `data/articles/zhao2019.pdf`, but its DOI is not in `papers.yml` and there is no dossier for it.
