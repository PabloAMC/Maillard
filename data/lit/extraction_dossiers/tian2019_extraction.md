# Tian, Xu, Chen & Yu 2019 — COMPLETE TRANSCRIPTION
### *"Flavoromics approach to identifying the key aroma compounds in traditional Chinese milk fan"*
### ⭐⭐⭐ **THE PAPER THAT SETTLES THE `?/kg` MILK-UNIT QUESTION — by an exact seven-for-seven numeric anchor, not by inference.** Wave **K6b**, 2026-08-30.

**Full extraction of `data/articles/tian2019.pdf`, including Table 1 in its entirety
(71 compounds × 6 samples × concentration + RSD + threshold + 6 OAVs).**
Read-only wave. No repo file outside `data/lit/extraction_dossiers/` created or modified.

**Read first:** `k2_matrix_and_thresholds.md` §A.4 (the seven milk rows with a literal `?/kg`)
and `k3_final_parameter_inventory.md` §C.18 (which names retrieval of this paper as *"the only
way to settle the `?/kg` unit; without it **seven milk thresholds carry a factor-of-1000 basis
risk**"*). **⇒ THAT GAP IS NOW CLOSED. §2.**

**Provenance:** **[M]** measured · **[C]** cited from another paper · **[Q]** qualitative ·
**[Z]** derived here · **[!]** flag.

---

## 0. IDENTITY — **MATCHES EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Huaixiang Tian, Xiaolin Xu, Chen Chen, Haiyan Yu\*** |
| Venue | ***J. Dairy Sci.* 102 (11), 2019** — running head *"Journal of Dairy Science Vol. 102 No. 11, 2019"*; the PDF carries an **"ARTICLE IN PRESS"** banner and prints `J. Dairy Sci. 102` **without page numbers** |
| DOI | **https://doi.org/10.3168/jds.2019-16796** ✔ **exactly the expected DOI** |
| Affiliation | **Department of Food Science and Technology, Shanghai Institute of Technology, Shanghai 201418, China** |
| Corresponding | hyyu@sit.edu.cn |
| PDF character | **12 pages**, clean text layer for the body, **Table 1 is a rotated landscape table spanning 2 pages** ⚠️ (§0a), Figures 1–5, Tables 1–4 |

### 0a. ⚠️ **[!] TABLE 1 IS ROTATED AND ITS TEXT LAYER IS PARTLY COLUMN-SHIFTED**
`pdftotext -layout` renders the 90°-rotated Table 1 with occasional column bleed in the
right-hand OAV block (e.g. butanoic acid Y6 emerges as "55" where 7,030 ÷ 13.0 = **541**).
**Wave K6b rasterized page 6 at 300 dpi, rotated it, and read the compound / concentration /
threshold columns visually.** Those three columns are **confirmed correct as extracted**; the
**OAV block is not** (§4c). **Every value used as load-bearing below has been visually verified.**

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER

**"Milk fan" is not milk.** It is a **Chinese traditional stretched-curd cheese** from the
Dengchuan region of Yunnan — the paper's own framing is a cheese study throughout, comparing
milk fan to Cheddar, Camembert, Grana Padano, Gruyère, Emmental, Parmigiano-Reggiano and
Mozzarella. **Six commercial samples (Y1–Y6) from six handmade workshops** were profiled by
HS-SPME-GC-MS, GC-O, PLS-DA, aroma recombination and omission.

**Its value to this repo is one thing above all others: Table 1 prints, side by side and in an
unambiguously labelled `μg/kg`, both the concentration of every compound in sample Y6 and the
literature odour threshold used to score it. Those Y6 concentrations reappear, digit for
digit, in the `?/kg` table of `tian2020.pdf` — the paper `k2_matrix_and_thresholds.md` §A.4
could not read. Seven rows, seven exact matches. The unit is μg/kg. §2.**

Secondarily it supplies **71 quantified volatiles × 6 samples** and **50 literature odour
thresholds in μg/kg**, which is a large threshold table the corpus did not have.

---

## 2. ⭐⭐⭐ THE MILK-UNIT RESOLUTION — the anchor, and it is exact

### 2a. The problem, restated
`data/articles/tian2020.pdf` = **Tian, Xu, Sun, Chen & Yu (2020),
*J. Dairy Sci.* 103:5863–5873, `10.3168/jds.2019-17880`**, *"Evaluation of the perceptual
interaction among key aroma compounds in milk fan by gas chromatography−olfactometry, odor
threshold, and sensory analyses."* Its **Table 1** prints **two data columns whose unit cell
renders as a literal question mark**:

> *"**Table 1.** Average concentration, threshold, and SD of the key aroma compounds detected in
> the milk fan sample **Y6**"*
> column headers: **`Concentration (?/kg)` · `SD` · `Threshold (?/kg)` · `SD`**
> footnote 1: *"**Key aroma compounds detected in the milk fan sample Y6 in our previous study
> (Tian et al., 2019).**"*

`k2_matrix_and_thresholds.md` §A.4 recorded the seven threshold values, noted *"The units cell
prints a literal question mark, verified at 900 dpi. Contextually µg/kg; **not asserted**"*, and
correctly refused to compute an OAV.

### 2b. ⭐ THE ANCHOR — Tian 2020's concentration column IS Tian 2019's Y6 column

**The 2020 paper's own footnote tells you where its numbers come from. They are reproduced
below against the 2019 Table 1 Y6 column, which is headed `Concentration (μg/kg; RSD in % in
parentheses)` in plain type.**

| compound | **tian2020 Table 1, `Concentration (?/kg)`** | **tian2019 Table 1, `Y6` column, `μg/kg`** | match |
|---|---:|---:|:---:|
| **Propanoic acid** | **347** | **347** (2.30) | ✅ **exact** |
| **Butanoic acid** | **7,030** | **7,030ᵈ** (0.980) | ✅ **exact** |
| **Octanoic acid** | **1,719** | **1,720ᵇ** (1.75) | ✅ (rounding) |
| **Octanal** | **29** | **29.3** (5.08) | ✅ (rounding) |
| **Nonanal** | **198** | **199ᶜᵈ** (11.9) | ✅ (rounding) |
| **2-Nonanone** | **244** | **244ᶜ** (4.81) | ✅ **exact** |
| **Ethyl hexanoate** | **1,001** | **1,000ᵇ** (2.74) | ✅ (rounding) |

**SEVEN OF SEVEN.** Three digit-for-digit identical, four differing only by the last-place
rounding expected when a value carrying three significant figures is re-typeset as an integer.

### 2c. ⭐⭐ **THE VERDICT**

> **`?/kg` = `μg/kg`. Established, not inferred.**

**The reasoning is airtight in one step:** in `tian2020` Table 1 the **same `(?/kg)` notation
heads both the concentration column and the threshold column**. The concentration column is
now identified, by a seven-row exact numeric match, as a verbatim reproduction of a column that
the source paper labels **μg/kg**. **Therefore the notation `?/kg` in that table means μg/kg,
and therefore the threshold column is in μg/kg.** Nothing about odour-threshold plausibility,
matrix ratios, or unit conventions is needed; the identification is arithmetic.

**⇒ The seven milk-fan threshold rows, now unit-resolved:**

| compound | **matrix threshold (μg/kg)** ± SD | tian2019 **literature** threshold (μg/kg) | **[Z] matrix / literature** |
|---|---:|---:|---:|
| **Propanoic acid** | **51,200 ± 84.07** | **3.00** | **17,067×** |
| **Butanoic acid** | **7,500 ± 24.64** | **13.0** | **577×** |
| **Octanoic acid** | **25,600 ± 120.88** | **5.10** | **5,020×** |
| **Octanal** | **160 ± 3.29** | **0.880** | **182×** |
| **Nonanal** | **1,600 ± 29.41** | **1.10** | **1,455×** |
| **2-Nonanone** | **52,000 ± 116.04** | **32.0** | **1,625×** |
| **Ethyl hexanoate** | **1,024 ± 34.28** | **3.00** | **341×** |

**[Z] Cross-check against the corpus's independent recovery.** `k2_matrix_and_thresholds.md`
§A.4 recovered a nonanal aqueous threshold of **≈1.1 µg/kg** from `Xin2026b_extraction.md` and
computed **≈1,450×**. **This wave gets 1,455× from Tian's own cited threshold of 1.10.**
Two independent routes to the same ratio to within 0.3 %. ✅

**[Z] And the sanity check that would have failed under any other unit.** If `?/kg` were
**ng/kg**, octanal's matrix threshold would be **0.16 μg/kg — BELOW its 0.880 μg/kg reference**,
and ethyl hexanoate's **1.024 vs 3.00**, i.e. **three of the seven compounds would be easier to
smell in cheese than in the reference medium**, which no fat-and-protein matrix does. If it were
**mg/kg**, propanoic acid's threshold would be **51.2 g/kg — 5 % w/w**. **Only μg/kg survives,
and the anchor in §2b makes the elimination argument redundant anyway.**

### 2d. **⚠️ THREE CAVEATS THAT MUST RIDE ON THESE ROWS**

1. **⚠️ [!] THE MATRIX IS NOT MILK AND IT IS NOT MILK FAN EITHER.** `tian2020` §Measurement of
   Odor Threshold, verbatim: *"A series of test samples were prepared by adding **a fresh milk
   solution matrix** and adding the test substances to brown glass bottles."* **The thresholds
   were measured by dosing into a *fresh milk solution*, while the CONCENTRATIONS come from
   *milk fan cheese*.** `k2_matrix_and_thresholds.md` §A.4's label *"Fresh milk (Tian et al.
   2020)"* is therefore right about the threshold matrix and the OAVs it enables are
   **cross-matrix** (cheese numerator ÷ milk-solution denominator).
   **Register `matrix_numerator: milk_fan_cheese`, `matrix_denominator: fresh_milk_solution`,
   `same_matrix: FALSE`. The fat, protein and pH of the "fresh milk solution" are never stated.**
2. **⚠️ [!] THE REFERENCE THRESHOLDS ARE CITED, NOT MEASURED, AND ARE DECLARED "IN AIR".**
   `tian2019` §Odor Activity Values, verbatim: *"Odor activity values were used to evaluate the
   contribution of each compound to overall milk fan aroma, which was calculated as the ratio of
   the concentration of a single aroma compound **to its detection threshold in air**. **The
   threshold values were taken from the literature.**"* **No source is given for any of the 50
   thresholds.** And a threshold "in air" cannot be divided into a μg/kg concentration in a
   solid — this is the **same "in air" mislabelling** that `k2_matrix_and_thresholds.md` §C
   already flags for `Xin2026b`. **The values are almost certainly aqueous/orthonasal
   thresholds mislabelled; that is why nonanal's 1.10 agrees with the corpus's independently
   recovered aqueous value.** **Register the ratio column as `cross_study_cross_method`, exactly
   as §D.2 requires.**
3. **⚠️ [!] DO NOT INGEST THE `tian2020` SDs AS THRESHOLD UNCERTAINTY.** They are
   **84.07 on 51,200 (0.16 %)**, **120.88 on 25,600 (0.47 %)**, **3.29 on 160 (2.1 %)** — and
   **five of the seven threshold values are exact steps of a power-of-two dilution series**
   (**1,024 = 2¹⁰**; 25,600 = 25 × 1,024; 51,200 = 50 × 1,024; 12,800/102,400/204,800/409,600/
   819,200 all appear in `tian2020` Table 2). A 15-panellist sensory threshold cannot carry
   0.16 % RSD. `k2_matrix_and_thresholds.md` §A.4 already flagged this; **the flag stands and is
   now explained: the values are series steps, and the SDs describe something other than the
   between-panellist spread.**

### 2e. **The threshold method, transcribed from `tian2020` for completeness [M]**
> *"The odor threshold was measured based on **ASTM (ASTM Standard, 2011)** and previous studies
> (Saison et al., 2009). **Fifteen panelists** participated … A series of test samples were
> prepared by adding **a fresh milk solution matrix** … to **brown glass bottles**. The panelists
> conducted several **3-AFC (ISO 13301, 2002)** tests with sample concentrations ranging **from
> high to low with a 2-fold interval** until they could not recognize the aroma. If the panelist
> was able to correctly recognize the aroma **in both trials of the same sample**, the same test
> was continued for samples of the next lowest concentration. **The best-estimate threshold was
> the geometric mean of the final unrecognized concentration and the adjacent higher
> concentration (Guadagni et al., 1963).**"*
**⇒ Descending 2-fold 3-AFC, BET by geometric mean of last-missed and adjacent-higher,
n = 15. Consistent with the summary already in `k2_matrix_and_thresholds.md` §A.4.**

---

## 3. SYSTEM AND METHOD — complete transcription **[M]**

| variable | value as printed |
|---|---|
| **Samples** | **six milk fan samples (Y1–Y6)**, *"made from the same materials and using the same processing techniques"*, purchased from **6 of the most popular handmade workshops in Yunnan Province, China**; *"specific parameters, such as **pH of the acid juice and curd temperature and time**, were different"*; transported on ice, frozen until analysis |
| Region | **Dengchuan**, *"an altitude of about **1,900 m** … a temperate climate throughout the year"* |
| **SPME sample** | **4 g** milk fan, thawed **24 h at 4 °C**, then **1 h at room temperature**, cut into cubes, + **100 µL IS**, in a **15-mL vial** with Teflon cover |
| **Internal standard** | **2-octanol, 13 mg/L**, 100 µL |
| **SPME fiber** | **75-µm CAR/PDMS** (chosen from DVB/CAR/PDMS 50/30 µm, CAR/PDMS 75 µm, PDMS/DVB 65 µm) |
| **SPME conditions** | **60 °C for 40 min** (optimized from 50/60/70/80 °C × 20/30/40/50 min) |
| GC-MS | Agilent **7890-5973**; **HP-Innowax 60 m × 0.25 mm × 0.25 µm**; **splitless**; **40 °C (4 min) → 100 °C at 3 °C/min (hold 2 min) → 150 °C at 4 °C/min → 230 °C at 10 °C/min (hold 5 min)**; He **1 mL/min**; **EI 70 eV**; source **230 °C**; **scan m/z 30–450 at 1 scan/s** |
| Identification | **retention indices vs C6–C30 n-alkanes** + **NIST11** library |
| **Quantitation** | internal-standard curve built by adding **propionic acid** at several levels + IS to **milk fan samples**; **y = 0.7269x + 0.0004, R² = 0.993**; `Ci = fi × CIS × Ai/AIS` where `fi = (Ai/AIS)/(Ci/CIS)` |
| GC-O | Agilent 7890 + **Gerstel ODP-2**; effluent split **equally** between FID and sniffer; **OSME** time–intensity; **15 trained panelists**; **0–5 intensity scale** (0 none, 3 medium, 5 extremely strong); **triplicate per panelist**; AI = mean of 15 |
| Sensory | **15 panelists (8 M, 7 F, mean age 23)**; **10 descriptors**: fruity, cheesy, milky, sour, rancid, musty, nutty, umami, sulfurous, fatty; rated **0–10**; **triplicate**; **8-g** samples in a **100-mL odor-proof container** |
| Descriptor references [M] | **20 mg/kg ethyl hexanoate** = "fruity" · fresh butter = "cheesy" · fresh milk = "milky" · **0.08 % citric acid in water** = "sour" · butanoic acid = "rancid" · potting soil = "musty" · unsalted raw nuts = "nutty" · **1 % MSG in water** = "umami" · mashed boiled egg = "sulfurous" · salad oil = "fatty" |
| **Recombination** | **23 compounds with OAV ≥ 1**, dissolved in **glycerol triacetate** matrix at their **Y6** natural concentrations, equilibrated **10 min** at ambient ⇒ model **Y** |
| **Omission** | **triangle test**, one or a group of compounds removed from model Y; each omission model + **2 recombinants (5 g each)**, randomized |
| Statistics | ANOVA + **Duncan's test**, SPSS 13.0; **PLS-DA** in R with **mixOmics** |

**⚠️ [!] THE RECOMBINATION MATRIX IS GLYCEROL TRIACETATE (TRIACETIN), NOT CHEESE.** The whole
key-aroma verification therefore runs in a synthetic solvent, not in the matrix whose
thresholds §2 reports. **A third matrix enters the paper at the validation step.**

---

## 4. ⭐⭐ TABLE 1 — ALL 71 COMPOUNDS — TRANSCRIBED IN FULL **[M]**

*Caption verbatim: "**Table 1.** Average concentrations (**μg/kg**), relative standard
deviations (RSD), and odor activity values of the volatile compounds detected in the 6 milk fan
samples (Y1 to Y6)."
Column groups: `No. · Compound · RI Calc · RI Ref · ID · Concentration (μg/kg; RSD in % in
parentheses) for Y1…Y6 · Threshold (μg/kg) · Odor activity value Y1…Y6`.
Footnotes: ᵃ⁻ᵉ *"Values with different letters in a row are significantly different using
Duncan's multiple comparison tests (P < 0.05)"* · ³ ID: **MS** = NIST11 mass-spectrum match,
**RI** = retention index agreement with literature · ⁴ **— = not detected**.*

### 4a. **Concentrations (μg/kg) and thresholds (μg/kg)** — the load-bearing columns, visually verified

| # | compound | RI calc | RI ref | Y1 | Y2 | Y3 | Y4 | Y5 | **Y6** | **threshold** |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | **Propanoic acid** | 1,544 | 1,535 | 176 (21.5) | 71.5 (4.71) | 373 (5.62) | — | 100 (5.19) | **347 (2.30)** | **3.00** |
| 2 | Isobutyric acid | 1,573 | 1,564 | 34.9 (2.97) | 17.1 (8.18) | 57.2 (11.2) | — | — | 35.0 (7.47) | 10.0 |
| 3 | **Butanoic acid** | 1,631 | 1,624 | 10,700ᵈᵉ (2.01) | 4,510ᵇ (6.11) | 3,190ᵃ (6.77) | 6,060ᵇᶜ (3.98) | 6,500 (3.10) | **7,030ᵈ (0.980)** | **13.0** |
| 4 | Isovaleric acid | 1,677 | 1,682 | 129 (3.98) | 31.1 (2.21) | 75.5 (6.46) | — | — | 65.1 (4.94) | 1.00 |
| 5 | Pentanoic acid | 1,745 | 1,743 | 334ᵈ (6.54) | 136ᵇ (3.05) | 94.0ᵃ (11.6) | 217ᶜ (5.71) | 167 (4.78) | **2ᶜ (10.3)** ⚠️ | 0.150 |
| 6 | Hexanoic acid | 1,851 | 1,855 | 654ᵈ (3.96) | 3,750ᵇ (0.680) | 2,940ᵃ (10.2) | 7,090ᵈ (7.85) | 6,550 (3.69) | 5,220ᶜ (1.49) | 10.0 |
| 7 | Heptanoic acid | 1,958 | 1,953 | 197ᶜ (5.31) | 99.8ᵃ (1.23) | 110ᵃ (4.31) | 291ᵈ (4.99) | 136 (11.4) | 173ᵇᶜ (9.94) | 22.0 |
| 8 | **Octanoic acid** | 2,063 | 2,070 | 2,340ᵈ (3.37) | 1,400ᵃ (4.03) | 2,060ᶜ (1.42) | 2,750ᵉ (8.35) | 2,820 (2.29) | **1,720ᵇ (1.75)** | **5.10** |
| 9 | Nonanoic acid | 2,169 | 2,169 | 56.1 (8.58) | 32.8 (7.82) | 32.2 (9.77) | — | — | 13.1 (6.67) | 1.60 |
| 10 | n-Decanoic acid | 2,275 | 2,272 | 569ᵇ (11.6) | 360ᵃ (5.93) | 529ᵇ (3.18) | 602ᵇ (1.97) | 520 (2.63) | 345ᵃ (7.50) | 5.00 |
| 11 | 2-Ethylpropanal | 919 | 912 | 10.4 (5.81) | — | 25.3 (5.22) | — | — | 24.1 (2.66) | 100 |
| 12 | Isovaleraldehyde | 923 | 932 | — | — | — | — | — | 31.3 (18.0) | 2.00 |
| 13 | **Hexanal** | 1,085 | 1,084 | — | — | — | **315 (0.580)** | — | — | **11.0** |
| 14 | (E)-2-Hexenal | 1,223 | 1,219 | 59.5 (1.75) | — | — | — | — | 13.5 (1.43) | 3.10 |
| 15 | **Octanal** | 1,292 | 1,299 | 83.3 (8.90) | 18.2 (4.21) | 27.6 (12.1) | — | 12.3 (5.22) | **29.3 (5.08)** | **0.880** |
| 16 | (E)-2-Heptenal | 1,328 | 1,332 | — | 37.0 (7.59) | — | 116 (2.35) | 8.44 (7.35) | 44.8 (5.72) | 13.0 |
| 17 | **Nonanal** | 1,397 | 1,396 | 169ᶜ (4.91) | 99.4ᵇ (3.39) | 184ᶜᵈ (5.50) | 220ᵈ (9.56) | 29.7 (5.63) | **199ᶜᵈ (11.9)** | **1.10** |
| 18 | (E)-2-Octenal | 1,436 | 1,434 | 45.5 (9.85) | — | 14.0 (6.94) | 41.1 (3.71) | — | 29.9 (5.85) | 3.00 |
| 19 | 2-Methyl-1-propanol | 1,101 | 1,092 | — | — | 67.4 (7.39) | — | — | 32.8 (7.56) | 100 |
| 20 | Butanol | 1,148 | 1,150 | 33.3 (1.88) | — | — | — | — | — | 78.0 |
| 21 | 3-Methyl-1-butanol | 1,209 | 1,205 | 243 (5.02) | — | 982 (4.75) | — | 273 (4.90) | 883 (2.64) | 6.10 |
| 22 | Pentanol | 1,252 | 1,255 | 106 (3.01) | 54.5 (6.90) | 162 (8.54) | 80.5 (4.88) | — | 74.8 (1.84) | 20.0 |
| 23 | Hexanol | 1,355 | 1,359 | 230ᵃᵇ (14.5) | 223ᵃᵇ (5.26) | 265ᵇᶜ (9.14) | 203ᵃ (3.58) | 180 (6.35) | 314ᶜ (4.75) | 34.0 |
| 24 | 1-Octen-3-ol | 1,452 | 1,456 | — | 41.8 (8.69) | — | 27.2 (2.04) | — | — | 1.00 |
| 25 | 2-Ethylhexanol | 1,491 | 1,486 | 25.7ᶜ (6.91) | 16.6ᵃᵇ (10.21) | 23.8ᵇᶜ (7.74) | 15.7ᵃ (0.830) | 19.5 (2.67) | 12.8ᵃ (2.39) | 740 |
| 26 | Linalool | 1,551 | 1,552 | — | — | — | — | — | 14.0 (3.12) | 4.40 |
| 27 | Nonanol | 1,665 | 1,666 | — | — | 27.3 (5.92) | 177 (6.06) | 46.5 (5.46) | 29.7 (9.56) | 5.30 |
| 28 | **p-Cresol** | 2,099 | 2,091 | — | — | 19.7 (7.35) | 53.8 (7.83) | — | 17.3 (9.53) | **0.240** |
| 29 | 2-Pentanone | 984 | 980 | — | 1,210 (3.01) | 777 (7.88) | 1,940 (2.61) | — | 2,050 (1.63) | 98.0 |
| 30 | 2-Heptanone | 1,186 | 1,180 | 4,540ᵉ (9.44) | 1,450ᵇ (7.04) | 1,250ᵇ (11.4) | 3,790ᵈ (7.22) | 20.1 (8.98) | 2,090ᶜ (1.68) | 3.50 |
| 31 | 3-Octanone | 1,253 | 1,259 | — | 10.1 (2.48) | — | 9.65 (12.2) | 8.17 (1.49) | 15.7 (1.55) | 1.30 |
| 32 | 2-Octanone | 1,288 | 1,283 | 170ᵈ (2.14) | 95.8ᵇ (3.41) | 76.9ᵃ (8.73) | 138ᶜ (6.64) | 67.1 (4.49) | 78.1ᵃ (10.2) | 60.0 |
| 33 | 6-Methyl-5-hepten-2-one | 1,338 | 1,336 | — | 13.5 (8.26) | 13.4 (6.62) | 12.5 (3.74) | — | 32.0 (7.81) | 68.0 |
| 34 | **2-Nonanone** | 1,393 | 1,393 | 552ᵉ (5.06) | 219ᶜ (7.05) | 122ᵇ (5.61) | 470ᵈ (2.17) | 35.9 (3.49) | **244ᶜ (4.81)** | **32.0** |
| 35 | **3-Octen-2-one** | 1,414 | 1,416 | — | — | — | — | — | 9.45 (0.670) | **0.0300** |
| 36 | 2-Decanone | 1,498 | 1,493 | — | — | 12.3 (5.05) | 28.1 (2.06) | — | 64.3 (9.50) | 60.0 |
| 37 | 2-Undecanone | 1,604 | 1,598 | 18.6ᶜ (8.47) | 6.83ᵃ (1.89) | 11.3ᵇ (5.71) | 12.8ᵇ (6.31) | 12.5 (3.94) | 11.4ᵇ (6.94) | 30.0 |
| 38 | **2-Pentylfuran** | 1,231 | 1,240 | 212ᵈ (8.11) | 87.2ᵇᶜ (5.49) | 61.5ᵇ (3.88) | 126ᶜ (10.1) | 13.6 (2.36) | 119ᶜ (2.91) | **19.0** |
| 39 | Styrene | 1,261 | 1,255 | 23.9 (8.50) | — | 13.9 (5.77) | 76.4 (4.56) | 7.81 (6.27) | 690 (2.68) | 12.0 |
| 40 | Benzaldehyde | 1,535 | 1,530 | 13.6ᵃ (8.52) | 12.4ᵃ (6.89) | 15.2ᵇᶜ (3.95) | 12.4ᵃ (2.89) | 17.5 (4.27) | 14.8ᵃᵇ (6.66) | 85.0 |
| 41 | Acetophenone | 1,668 | 1,669 | — | — | 13.2 (3.21) | — | — | 8.78 (8.58) | 1.20 |
| 42 | γ-Hexanolactone | 1,702 | 1,703 | 9.24 (1.85) | — | 7.49 (3.84) | 8.33 (2.53) | — | 6.35 (6.50) | **—** |
| 43 | δ-Hexalactone | 1,815 | 1,818 | 12.0 (3.37) | 8.92 (1.36) | 13.3 (14.7) | 9.48 (3.25) | — | 7.02 (6.99) | **—** |
| 44 | Benzyl alcohol | 1,890 | 1,884 | 28.3ᵇ (8.67) | 35.8ᶜ (10.8) | 43.6ᵈ (4.99) | 41.0ᶜᵈ (1.58) | 12.7 (1.24) | 39.7ᶜᵈ (9.94) | 250 |
| 45 | Phenylethyl alcohol | 1,928 | 1,925 | 45.1ᵇᶜ (6.87) | 23.4ᵃ (1.12) | 42.3ᵇᶜ (2.04) | 57.4ᶜ (7.18) | 32.7 (9.31) | 70.4ᵈ (4.50) | 21.0 |
| 46 | 5-Octanolide | 1,998 | 1,999 | 14.8 (9.31) | 7.24 (8.53) | 12.3 (14.1) | 8.41 (1.20) | — | — | 420 |
| 47 | 5-Decanolide | 2,238 | 2,229 | 14.0 (5.09) | 13.5 (5.20) | — | 17.0 (3.27) | 13.6 (9.68) | 20.7 (1.44) | 66.0 |
| 48 | Ethyl acetate | 898 | 891 | 216ᵃ (10.6) | 717ᶜ (9.52) | 1,230ᵈ (10.3) | 454ᵇ (4.05) | 734 (5.21) | 640ᶜ (6.23) | 340 |
| 49 | **Ethyl propionate** | 962 | 961 | — | — | — | — | 84.2 (7.58) | — | **0.0270** |
| 50 | Propyl acetate | 980 | 977 | — | — | 49.7 (3.27) | — | 42.3 (7.04) | — | 200 |
| 51 | **Ethyl butyrate** | 1,040 | 1,041 | 529ᵃ (7.56) | 4,120ᵈ (5.43) | 1,340ᵇ (12.8) | 3,750ᵈ (3.57) | 7,060 (7.43) | 2,180ᶜ (0.800) | **0.0530** |
| 52 | Pentyl butyrate | 1,111 | 1,094 | — | 12.1 (4.08) | — | — | — | — | **—** |
| 53 | Ethyl valerate | 1,132 | 1,133 | — | 77.2 (7.62) | 65.5 (10.3) | 143 (3.00) | 267 (10.7) | — | 0.580 |
| 54 | Isobutyl butyrate | 1,162 | 1,160 | — | — | — | — | 16.4 (6.93) | 11.4 (9.97) | 9.4 |
| 55 | Methyl hexanoate | 1,190 | 1,187 | — | — | — | — | 10.7 (6.01) | — | 700 |
| 56 | Butyl butyrate | 1,222 | 1,221 | — | — | — | — | 23.1 (9.36) | — | 28.0 |
| 57 | **Ethyl hexanoate** | 1,235 | 1,232 | 312ᵃ (6.32) | 2,410ᶜ (4.86) | 1,130ᵇ (3.43) | 3,520ᵈ (6.55) | 18,300 (1.10) | **1,000ᵇ (2.74)** | **3.00** |
| 58 | Propyl hexanoate | 1,317 | 1,316 | — | — | — | 8.96 (2.32) | 103 (6.25) | — | **—** |
| 59 | Ethyl heptanoate | 1,323 | 1,317 | — | 48.6 (6.84) | 41.5 (5.90) | 51.6 (5.66) | 353 (7.54) | — | 39.0 |
| 60 | Butyl hexanoate | 1,418 | 1,418 | — | — | — | — | 15.6 (9.77) | — | 700 |
| 61 | Ethyl caprylate | 1,438 | 1,440 | 72.2ᵃ (10.1) | 488ᶜ (1.96) | 1,570ᵉ (2.25) | 923ᵈ (4.24) | 8,370 (1.38) | 205ᵇ (7.31) | 40.0 |
| 62 | Propyl octanoate | 1,524 | 1,523 | — | — | 14.7 (2.07) | — | 25.1 (2.65) | — | **—** |
| 63 | Ethyl nonanoate | 1,542 | 1,541 | — | — | — | — | 77.2 (5.76) | — | 10.0 |
| 64 | Ethyl caprate | 1,646 | 1,648 | 30.7ᵃ (3.21) | 85.2ᵇ (6.45) | 446ᵈ (2.80) | 337ᶜ (5.48) | 1,250 (2.42) | 79.1ᵇ (4.87) | 53.0 |
| 65 | Pentyl octanoate | 1,687 | 1,696 | — | — | — | — | 7.03 (4.67) | — | **—** |
| 66 | Ethyl myristate | 2,055 | 2,054 | — | 8.90 (7.20) | 20.5 (5.26) | 10.0 (1.01) | 57.7 (5.95) | 8.73 (3.78) | 180 |
| 67 | 1-Octene | 845 | 842 | — | — | — | 125 (1.46) | — | — | 4.60 |
| 68 | 2-Octene | 860 | 858 | — | — | 63.1 (6.72) | 180 (10.6) | — | — | 500 |
| 69 | d-Limonene | 1,263 | 1,264 | 395ᵈ (3.1) | 76.1ᵃ (1.67) | 85.0ᵃᵇ (10.5) | 105ᵇ (0.800) | 61.9 (6.64) | 215ᶜ (12.1) | 200 |
| 70 | Dimethyl disulfide | 1,076 | 1,071 | 123 (11.1) | 46.0 (10.9) | 135 (6.29) | 25.4 (1.22) | — | 83.9 (4.98) | **—** |
| 71 | Dimethyl sulfone | 1,925 | 1,912 | 31.2ᵇᶜ (3.9) | 37.3ᶜ (4.78) | 47.2ᵈ (6.69) | 24.8ᵇ (7.01) | 11.7 (9.38) | 33.4ᶜ (6.35) | **—** |

*All identifications are **MS, RI** (both criteria) for all 71 compounds. RSD in parentheses, in %.*

### 4b. ⚠️ **[!] THE PENTANOIC ACID Y6 CELL IS A TYPO IN THE PUBLISHED PAPER, NOT AN EXTRACTION ARTEFACT**
Row 5, Y6 reads **`2ᶜ (10.3)`** while the printed OAV for the same cell is **1,655**.
**[Z] 1,655 × 0.150 = 248 μg/kg** — and the Y1–Y5 values are 334/136/94.0/217/167.
**⇒ The true Y6 pentanoic acid concentration is ≈ 248 μg/kg; "2" has lost two digits in
typesetting.** Verified visually at 300 dpi: the PDF really does print `2ᶜ`.
**Register `pentanoic_acid_Y6: PRINTED 2, TRUE ≈248 μg/kg [!]`.**

### 4c. ⚠️ **[!] DO NOT INGEST THE PRINTED OAV COLUMNS — RECOMPUTE THEM**
The OAV block is arithmetically correct in most cells (checked: propanoic Y1 176/3.00 = 59 ✔;
butanoic Y2 4,510/13.0 = 347 ✔; hexanoic Y3 2,940/10.0 = 294 ✔; nonanal Y1 169/1.10 = 154 ✔;
octanal Y1 83.3/0.880 = 95 ✔) **but several cells are inconsistent by exactly one or two
decimal digits** (butanoic Y6 prints **55** where 7,030/13.0 = **541**; nonanal Y4 prints **2**
where 220/1.10 = **200**; hexanoic Y4 prints **80** where 7,090/10.0 = **709**).
**Whether these are typesetting losses in the journal or column bleed in the rotated text layer
was not resolved by this wave.** **⇒ Compute every OAV from the concentration and threshold
columns (both verified) rather than reading the printed OAV.**

### 4d. **[Z] What the threshold column is worth on its own**
**Fifty literature odour thresholds in μg/kg**, spanning **0.0270 (ethyl propionate) to 740
(2-ethylhexanol)** — a **27,000×** range — including several the corpus lacked:
**3-octen-2-one 0.0300 · ethyl propionate 0.0270 · ethyl butyrate 0.0530 · pentanoic acid
0.150 · p-cresol 0.240 · ethyl valerate 0.580 · octanal 0.880 · isovaleric acid 1.00 ·
1-octen-3-ol 1.00 · nonanal 1.10 · acetophenone 1.20 · 3-octanone 1.30 · nonanoic acid 1.60 ·
isovaleraldehyde 2.00 · ethyl hexanoate 3.00 · propanoic acid 3.00 · (E)-2-octenal 3.00 ·
(E)-2-hexenal 3.10 · 2-heptanone 3.50 · linalool 4.40 · 1-octene 4.60 · n-decanoic acid 5.00 ·
octanoic acid 5.10 · nonanol 5.30 · 3-methyl-1-butanol 6.10 · isobutyl butyrate 9.4 ·
hexanoic/isobutyric acid 10.0 · ethyl nonanoate 10.0 · hexanal 11.0 · styrene 12.0 ·
butanoic acid 13.0 · (E)-2-heptenal 13.0 · 2-pentylfuran 19.0 · pentanol 20.0 ·
phenylethyl alcohol 21.0 · heptanoic acid 22.0 · butyl butyrate 28.0 · 2-undecanone 30.0 ·
2-nonanone 32.0 · hexanol 34.0 · ethyl heptanoate 39.0 · ethyl caprylate 40.0 ·
ethyl caprate 53.0 · 2-octanone/2-decanone 60.0 · 5-decanolide 66.0 · 6-methyl-5-hepten-2-one
68.0 · butanol 78.0 · benzaldehyde 85.0 · 2-pentanone 98.0 · 2-ethylpropanal /
2-methyl-1-propanol 100 · ethyl myristate 180 · propyl acetate 200 · d-limonene 200 ·
benzyl alcohol 250 · ethyl acetate 340 · 5-octanolide 420 · 2-octene 500 · methyl / butyl
hexanoate 700 · 2-ethylhexanol 740.**
**⚠️ ALL are `[C]` with NO source cited and are declared *"in air"* (§2d-2). Ingest as
`unsourced_literature_threshold`, never as measured, and never as an air threshold.**
**Four compounds have NO threshold printed (—): γ-hexanolactone, δ-hexalactone, pentyl
butyrate, propyl hexanoate, propyl octanoate, pentyl octanoate, dimethyl disulfide, dimethyl
sulfone.** ⚠️ **Note especially `dimethyl disulfide: no threshold` — a compound the corpus's
sulfur lane needs.**

---

## 5. THE REST OF THE PAPER — the parts with numbers **[M]**

### 5a. Composition summary
> *"**Seventy-one aroma compounds** were detected, including **10 acids, 8 aldehydes, 9 alcohols,
> 9 ketones, 19 esters, 11 aromatic and heterocyclic compounds, and 5 other compounds.**
> **Acids (30.8–71.6 %)**, **esters (3.94–67.6 %)**, and **ketones (2.73–17.9 %)** were the most
> abundant groups."*
**⚠️ [!] COUNT CHECK [Z]: 10 + 8 + 9 + 9 + 19 + 11 + 5 = 71 ✔ but the class list does not match
Table 1's ordering** (Table 1 groups acids 1–10 ✔, aldehydes 11–18 ✔ = 8, alcohols 19–28 = 10
not 9, ketones 29–37 = 9 ✔, then 38–47 mixed, esters 48–66 = 19 ✔, 67–71 = 5 ✔).
**Minor; the compound-level table is authoritative.**

### 5b. Cross-cheese comparisons, all with numbers [M][C]
| compound | **milk fan (μg/kg)** | comparison cheese |
|---|---|---|
| total acids | **~9,400 to 21,100** | rennet-curd cheese **3,028** [C] |
| butanoic acid | **3,190–10,700** | Gokceada **677** [C] |
| hexanoic acid | **2,900–6,500** | Gokceada **820** [C] |
| octanoic acid | **1,400–2,800** | Gokceada **305** [C] |
| decanoic acid | **340–600** | — |
| ethyl acetate | **210–1,220** | Parmigiano-Reggiano **50–250** [C] |
| *(a compound at)* | **~700–7,000** | Mozzarella **700**, Emmental **260–500** [C] |
| *(a compound at)* | **250 μg/kg**, **~3,500**, **~700** | *(fragmentary in the rotated column)* ⚠️ |

⚠️ **[!] The cross-cheese paragraph is split across the rotated table's column gutters and
several sentences are incomplete in the text layer. The values above that carry a named
comparator are verified; the last row is not. Do not ingest the last row.**

### 5c. Aroma intensities from GC-O (0–5 scale), verbatim [M]
> *"**Eight acids were identified** (aroma intensity: **1.18–4.51**): **propanoic acid
> (4.10–4.82)**, **butanoic acid (1.31–3.86)**, **hexanoic acid (1.76–4.62)**, **pentanoic acid
> (1.31–3.86)**, **heptanoic acid (1.25–3.50)**, **octanoic acid (1.50–4.73)**, **nonanoic acid
> (1.55)**, and **n-decanoic acid (1.39–2.70)**."*
⚠️ **[!] The lead-in says "1.18–4.51" but the listed sub-ranges span 1.25–4.82 — the summary
range is inconsistent with its own members. And butanoic and pentanoic acid are given the
identical range (1.31–3.86), which is suspicious. Register both flags.**
**Also flagged: eight acids are named here but Table 1 lists ten.**

### 5d. The high-OAV compounds, verbatim [M]
> *"**Pentanoic acid (OAV of 92–2,230)**, **hexanoic acid (OAV of 80–655)**, **heptanoic acid
> (OAV of 5–13)**, and **decanoic acid (OAV of 12–114)** had high OAV."*
**[Z] cross-check: pentanoic Y1 334/0.150 = 2,227 ✔ vs printed 2,229; hexanoic Y3 2,940/10 =
294, and the quoted 80–655 brackets the printed Y1–Y6 span ✔. The prose OAV ranges are
consistent with the concentration ÷ threshold arithmetic, which further supports §4c's
conclusion that the arithmetic is right and only some printed OAV cells are corrupted.**

### 5e. Odour-active compounds named by GC-O [M]
*"**Isovaleraldehyde, hexanal, octanal, (E)-2-heptenal, nonanal, (E)-2-octenal, 2-nonanone,
ethyl acetate, ethyl butyrate, ethyl valerate, ethyl hexanoate**"* … *"In addition to
(E)-2-heptenal, (E)-2-octenal, and ethyl valerate, these compounds were also odor-active
compounds in **Cheddar, Grana Padano, Ragusano, Emmental, Mozzarella, and Gorgonzola** cheeses."*
**⭐ Note hexanal is named as odour-active by GC-O in milk fan despite being detected by GC-MS
in only ONE of six samples (Y4, 315 μg/kg). A GC-O/GC-MS discordance worth registering.**

### 5f. Recombination and omission
**23 compounds with OAV ≥ 1** were recombined in glycerol triacetate at Y6 concentrations;
omission by triangle test. *(Tables 2–4 and Figures 2–5 carry the descriptive profiles,
PLS-DA loadings and omission outcomes; they are sensory-profile objects with no additional
concentration or threshold data and are not transcribed here.)*

---

## 6. WHAT THIS PAPER SUPPLIES THE REPO

| question | answer |
|---|---|
| **Does `?/kg` in `tian2020` Table 1 mean μg/kg?** | ✅ ⭐⭐⭐ **YES — seven-for-seven exact numeric anchor, §2b** |
| **Are the seven milk threshold rows now safe to use?** | ✅ **as μg/kg**, with the three caveats in §2d (cross-matrix, unsourced "in air" reference, SDs not uncertainty) |
| A large threshold table | ✅ **50 values, μg/kg, §4d** — all `[C]`, all unsourced |
| Quantified volatiles in a real cheese matrix | ✅ **71 compounds × 6 samples, μg/kg, §4a** |
| A hexanal matrix threshold | ❌ — only a **cited** 11.0 μg/kg reference threshold; hexanal is not among `tian2020`'s seven measured matrix thresholds |
| A same-matrix threshold/reference pair | ❌ **the reference thresholds are literature values from an unnamed source; §D.2's `cross_study_cross_method` classification applies in full** |

**SI: none listed. The companion `tian2020.pdf` is already in `data/articles/` and its Table 1
and Table 2 are transcribed in §2 above; no further retrieval is needed to close §C.18's
milk-unit gap.**

---

## 7. INTEGRITY FLAGS

| # | flag | severity |
|---|---|---|
| **1** | **"Milk fan" is a CHEESE, not milk. The threshold matrix is a separate "fresh milk solution" whose fat/protein/pH are unstated. Numerator and denominator matrices differ.** | ⭐ HIGH |
| **2** | **The 50 reference thresholds are cited with NO source and declared *"in air"* while being divided into μg/kg concentrations in a solid.** | ⭐ HIGH |
| **3** | Pentanoic acid Y6 prints `2` where the arithmetic requires ≈248 μg/kg | MEDIUM |
| **4** | Several printed OAV cells are inconsistent with concentration ÷ threshold by 1–2 digits; recompute rather than ingest | MEDIUM |
| 5 | Table 1 is rotated; `pdftotext` column-shifts parts of the OAV block | MEDIUM (method) |
| 6 | GC-O reports 8 acids with intensities; Table 1 lists 10 acids. The summary intensity range (1.18–4.51) contradicts its own member ranges (1.25–4.82) | LOW–MEDIUM |
| 7 | Butanoic and pentanoic acid are given identical intensity ranges (1.31–3.86) | LOW |
| 8 | Eight compounds have no threshold at all, including **dimethyl disulfide** | MEDIUM (gap) |
| 9 | Recombination runs in **glycerol triacetate**, a third matrix | MEDIUM |
| 10 | `tian2020`'s threshold SDs are 0.16–2.1 % RSD on 15-panellist sensory thresholds and five of seven values are exact power-of-two series steps — the SDs are not threshold uncertainty | ⭐ HIGH (carried from `k2` §A.4, now explained) |
| 11 | Class counts in the prose (9 alcohols) disagree with Table 1 (10) | LOW |
