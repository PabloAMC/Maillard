# Xin et al. 2026b — Wave K3 extraction 2026-08-28

**Xin, R.; Cui, S.; Liu, X.; Nie, C.; Huang, X.; Qin, L. "Multiprotein blends-regulated
high-moisture extruded seafood analogs: texture, structure, and flavor profile."
*Food Research International* **233** (2026) 119010.**
DOI 10.1016/j.foodres.2026.119010.
School of Food Science and Technology, State Key Laboratory of Marine Food Processing &
Safety Control, National Engineering Research Center of Seafood, **Dalian Polytechnic
University**, Dalian 116034, China.
Received 24 November 2025; revised 16 March 2026; accepted 16 March 2026; online 18 March 2026.
Funding: Dalian 2023JB11SN007; NSFC 32472290, 32572615.

**Source / how it was read:** `data/articles/Xin2026b.pdf`, **15 pages**, PDF text layer
extracted in full (layout mode) to `k3_xin2026b.txt` — **the entire body text, Table 1 and
Table 2 came from the text layer, not from renders.** Figures 1, 3, 4 were re-read from
**pdftoppm renders at 400–600 dpi** into `k3img/` (`fig1_top-05.png`, `fig1_mid-05.png`,
`fig3_AB-07.png`, `fig3_all-07.png`, `fig4B-10.png`, `fig4D-10.png`, `fig4E_left-10.png`,
`fig4E_inset-10.png`). Every number below carries a page anchor; figure-derived numbers are
tagged **[fig]** with the dpi.

---

## 0. THE HEADLINE — FOUR FINDINGS

**1. ⚠️ JACKPOT: PARTIAL — and it is worth more than a full table would have been.**
The odour-threshold *table* is **Table S5, online supplementary, NOT in this PDF**. BUT:
**one threshold is printed verbatim in the body text with units** (1-octen-3-ol,
**1.5 μg/kg**, p. 9), and **four thresholds are recoverable exactly by back-solving the
paper's own printed OAV ranges against its own printed concentrations** — three of them
**exact to 4–5 significant figures at BOTH ends of the range**. **No source citation is
given for any threshold, anywhere in the paper.** And the recovered values are
**unambiguously the classical AQUEOUS (water) thresholds**, while §2.11 declares the
denominator to be an **"air threshold value (μg/kg)"**. See §4.

**2. ⚠️ THE CROSS-PAPER REPRODUCIBILITY BOMB.** This paper's **ENPI** formulation
(**70 SPI + 20 NPI + 10 WG**) is *bit-for-bit* the companion paper's base blend
(**NPI:SPI:WG = 20:70:10**), run on **the same extruder at the same seven barrel
temperatures, the same 240 rpm, the same 2 kg/h, the same 70 °C cooling die, the same ~65 %
moisture**, extracted with **the same 2 g / 20 μL of 40 μg/mL cyclohexanone** protocol, by
**the same six authors**. So the companion's **no-carbohydrate control** and this paper's
**ENPI** are nominally the same material. Their reported volatile concentrations differ by
**~10× (hexanal), ~17× (2-pentylfuran), ~23× (nonanal)**, with **non-overlapping ranges**
across the whole sample sets for 2-pentylfuran and nonanal. See §6 — this is a direct,
same-lab measurement of how far these semiquantitative μg/kg numbers travel: **they don't.**

**3. A 12.4× swing in total aldehydes from swapping 20 % of the protein.** ERPI 2398.22 vs
ENPI 193.95 μg/kg (p. 9). Hexanal 1264.62 (EPPI) vs 74.86 (ENPI) = **16.9×**; nonanal 195.04
(ERPI) vs 12.44 (ENPI) = **15.7×**. **Protein source, at constant process and zero added
carbohydrate, moves lipid-oxidation aldehydes by more than an order of magnitude.**

**4. Furans dominate the volatile mass and are almost entirely one compound.** Total furans
**≈18 900 μg/kg in ERPI [fig]**, of which **2-pentylfuran is 17 033.09 (90 %)** (p. 13).
Same SPME-partition-artefact hazard as the companion, one order of magnitude smaller.

---

## 1. SYSTEM DEFINITION — VERBATIM

### 1a. The formulations (Table 1, p. 2 — transcribed verbatim, MEASURED design)

> **"Table 1. Formulations of mixed protein pastes before extrusion and corresponding
> extrudates."**

| No. | Mixed protein pastes before extrusion | Extruded samples | Composition (%, w/w, dry basis) |
|---|---|---|---|
| 1 | SPIP | **ESPI** | **70 SPI + 20 SPI + 10 WG** |
| 2 | NPIP | **ENPI** | **70 SPI + 20 NPI + 10 WG** |
| 3 | RPIP | **ERPI** | **70 SPI + 20 RPI + 10 WG** |
| 4 | YPP | **EYP** | **70 SPI + 20 YP + 10 WG** |
| 5 | MPIP | **EMPI** | **70 SPI + 20 MPI + 10 WG** |
| 6 | PPIP | **EPPI** | **70 SPI + 20 PPI + 10 WG** |

> Note (verbatim): "Abbreviations: SPI, soybean protein isolate; WG, wheat gluten; RPI, rice
> protein isolate; NPI, peanut protein isolate; PPI, pea protein isolate; MPI, mung bean
> protein isolate; YP, yeast protein. The suffix "P" denotes mixed protein pastes before
> extrusion, whereas the prefix "E" indicates samples obtained after high-moisture extrusion."

⚠️ **Row 1 as printed reads "70 SPI + 20 SPI + 10 WG"** — SPI listed twice. The control is
therefore **90 % SPI + 10 % WG**. Not a typo to fix silently: it is the *reference* arm and
its composition is only legible by arithmetic.
⚠️ **The substitution is SUBSTITUTIVE, not additive** — the alternative protein *replaces*
20 % SPI. (Contrast the companion paper, where the carbohydrate was **added on top of** the
protein blend.) So the comparison here is at **constant total protein**, which is cleaner.
⚠️ **No carbohydrate is added in this study at all.** Every Maillard-relevant sugar is
whatever the isolates carry endogenously (Table S1, online, not in this PDF).

**Protein suppliers, verbatim p. 2:**
> "SPI was obtained from Linyi Shansong Biological Products Co., Ltd. (Shangdong Province,
> China); RPI was obtained from Wuxi Jinnong Biotechnology Co., Ltd. (Jiangsu Province,
> China); MPI was obtained from Harbin Hada starch Co., Ltd. (Heilongjiang Province, China);
> NPI was obtained from Qingdao Changshou Food Co. Ltd. (Shangdong Province, China); PPI was
> obtained from Shandong Yuwang Ecological Food Industry Co., Ltd. (Shangdong Province,
> China); YP was obtained from Angel Yeast Co., Ltd. (Hubei Province, China); WG was obtained
> from Henan Lanchi Biotechnology Co., Ltd. (Henan Province, China)."

⚠️ **The proximate composition of the seven proteins is Table S1 — online, NOT in this PDF.**
This is the single most damaging omission for the repo: **without it there is no protein,
fat, or carbohydrate content for any arm**, so the lipid-oxidation and Maillard substrate
loading cannot be normalised. (The companion paper prints SPI 86.5 % protein / 2.4 % carb /
4.2 % ash / 0.1 % fat, NPI 77.1 / 10.5 / 0.9 / 5.3, WG 72.9 / 15.7 / 0.91 / 0.5 — **the same
SPI, NPI and WG suppliers**, so those three can be carried across, but RPI, MPI, PPI and YP
are unknown.) Methods used: **GB 5009.4-2016 (ash), GB 5009.3-2016 (moisture), GB 5009.5-2016
(protein), 5009.6-2016 (fat), GB 28050-2011 (carbohydrate)** (p. 2).

### 1b. The extrusion conditions (§2.3, p. 3 — verbatim, MEASURED/SET)

> "Based on preliminary experiments in conjunction with the extrusion conditions reported in
> previous studies (Tian et al., 2023; Zhao et al., 2024), the formulation of HME protein raw
> materials (Table 1) and the corresponding extrusion parameters were established for this
> study. **A corotating intermeshing twin-screw extruder** (Hunan Fumaco Food Engineering
> Technology Co., Ltd., Changsha, China) **with a screw length of 704 mm, a screw diameter of
> 22 mm (Length/Diameter = 32:1) and a cooling die** was used for all extrusion experiments.
> **Nine zones were created in the barrel: the feeding zone (T1), which serves as a gate for
> solid feeding without temperature control; seven temperature-controlled zones (T2-T7); a
> cooling zone kept at 70 °C by a tap water line. From the feeding zone to the die, the
> temperatures in the seven zones were adjusted to 60 °C, 80 °C, 110 °C, 135 °C, 140 °C,
> 140 °C, and 140 °C, in that order. 240 rpm and 2 kg/h were the screw speed and feed rate
> settings. The final product has a moisture content of around 65 %.** After the extrusion
> process reached a steady state, the extrudates had been collected at the end of the cooling
> die. After cooling to ambient temperature, **the extrudates were kept for further
> examination at −20 °C. The samples were allowed to equilibrate at room temperature for at
> least two hours before being analyzed after being thawed overnight in a refrigerator set at
> 4 °C.**"

**Consolidated process block:**

| parameter | value | status |
|---|---|---|
| extruder | corotating intermeshing **twin-screw**, Hunan Fumaco | set |
| screw | **704 mm × 22 mm, L/D = 32:1** | set |
| barrel zones | **9 total**: T1 feed (uncontrolled) + **T2–T7 seven controlled** + cooling zone | set |
| zone temperatures (feed→die) | **60 / 80 / 110 / 135 / 140 / 140 / 140 °C** | set |
| cooling die | **70 °C**, tap-water line | set |
| screw speed | **240 rpm** | set |
| feed rate | **2 kg/h** | set |
| **final product moisture** | **≈65 %** | stated |
| post-extrusion | collected at cooling-die exit, cooled to ambient, **stored −20 °C**; thawed overnight at 4 °C then **≥2 h at RT** before analysis | set |
| replication | **n = 3** (inferred from Fig. 4A PLS-DA: ESPI-1/2/3 … EPPI-1/2/3, p. 10 [fig]) | inferred |

⚠️ **The counting of zones is internally inconsistent as printed:** "Nine zones … the feeding
zone (T1) … seven temperature-controlled zones (T2-T7) … a cooling zone". T2–T7 is **six**
labels, not seven; and 1 + 7 + 1 = 9 only if T2–T7 is read as seven. The seven listed
temperatures resolve it (T2…T8 in effect). **Record the seven temperatures, discard the
zone-label arithmetic.**
⚠️ **The water-injection rate is NOT given** (the companion gives 3.7 kg/h, which with 2 kg/h
dry feed yields 64.9 % — consistent with "around 65 %"). **This wave's inference, not the
paper's:** water injection was almost certainly **3.7 kg/h**, as in the companion.
⚠️ **No die dimensions, no screw profile, no specific mechanical energy, no die pressure, no
residence time.**

### 1c. The reference seafood samples (§2.1, p. 2 — verbatim)

> "**Seven commercial grilled fish reference samples** (including **eel (*Astroconger
> Myriaster*), cod (*Gadus morhua*), Squid (*Teuthida*), monkfish (*Lophius americanus*),
> octopu (*Octopodidae*), needlefish (*Belonidae*), Black Scraper (*Thamnaconus modestus*)**)
> were purchased from **Dalian Shuiyifang Trading Co., Ltd.**"

They are used **only as texture and colour comparators** (Fig. 1E PCA, p. 5; Fig. S2 colour
PCA, online). **No volatile data were collected on the seafood references** — they are not a
flavour benchmark. ("octopu" is as printed.)

---

## 2. THE VOLATILE METHOD — verbatim, with the OAV definition

**§2.11, p. 3 (verbatim):**
> "With some adjustments, the technique described by Chen et al. (2023) was used to extract
> volatile chemicals using **headspace solid-phase microextraction (SPME)**. A **20 mL
> headspace vial** was filled with **20 μL (40 μg/mL) of a cyclohexanone internal standard
> solution** and a **2 g sample that had been chopped**. A **5977 A-7890B GC–MS system**
> (Agilent) fitted with a **non-polar HP-5MS capillary column (30 m × 0.25 mm × 0.25 μm)**
> was used … **splitless mode**; carrier gas **high purity helium**; **injector 250 °C**;
> temperature program **45 °C for 3 min, [in]creased to 250 °C at a rate of 5 °C/min, and
> stayed at 250 °C for 10 min**. … **ion source EI; ionization voltage 70 eV; full scan;
> ion fragment scanning range 35–500 m/z; quadrupole 150 °C; ion source 230 °C; transfer line
> 280 °C.**
> The **retention index (RI) was calculated using n-alkane (C7–C30) external references**. RI,
> retention time (RT), and mass spectra were compared with **GC–MS NEST libraries** in order
> to identify volatiles."

⇒ IS loading: 20 μL × 40 μg/mL = **0.8 μg into 2 g = 400 μg/kg nominal** (identical to the
companion). Data processing: **MS-DIAL 3.96**; stats **SPSS 23.0**, ANOVA + **Tukey HSD**,
p < 0.05; **MetaboAnalyst 5.0** for heatmaps/PCA/PLS-DA (§2.15, p. 4).

⚠️ **THE SPME FIBRE COATING IS NEVER STATED.** Nor the incubation temperature, incubation
time, or extraction time. The companion specifies **DVB/CAR/PDMS**; this paper specifies
nothing. **The single most response-factor-determining choice in the method is absent.**
⚠️ **"GC–MS NEST libraries"** — almost certainly **NIST**, misspelled. The library *version*
is not given (companion: NIST14 + Wiley11).
⚠️ **C7–C30** alkane ladder here vs **C6–C30** in the companion. Minor, but it means the two
papers' RI scales are not identical at the light end.
⚠️ **Quantification basis:** §2.11 never states the semi-quantification formula. §3.9 (p. 9)
calls it "**qualitative and relative quantitative analysis**" and the text says "**relative
content**" throughout; Fig. 4E's axis is **μg/kg**. **Single internal standard, no
compound-specific response factors ⇒ these are peak-area ratios wearing μg/kg clothing.**

### ⚠️ THE OAV DEFINITION, verbatim, §2.11, p. 3

> "The contribution of volatile chemicals to the sample's flavor profile was evaluated using
> the **Odor Activity Value (OAV)**. The OAV for every constituent was determined by
> **dividing the volatile compound's concentration in the sample (μg/kg) by its air threshold
> value (μg/kg)**."

**This is verbatim the same defective sentence as the companion paper** (`Xin2026_extraction.md`
§2), reproduced in a different journal. Three defects:
1. **"air threshold value (μg/kg)"** — air thresholds are given in μg/m³ or ng/L of air, never
   in μg/kg. A μg/kg-of-air threshold is a unit essentially never used.
2. **Dividing a matrix concentration by an air threshold is a category error** — numerator per
   kg of a 65 %-moisture protein solid, denominator per kg of air, **no partition coefficient
   anywhere**.
3. **They did not use air thresholds.** §4 proves it four ways, and the paper itself prints an
   aqueous threshold in the body text (1-octen-3-ol, 1.5 μg/kg, p. 9).

**⇒ SECOND INDEPENDENT INSTANCE, same group, different journal, of the exact practice the
mission exists to correct.** Two papers is a habit, not a slip.

---

## 3. ⚠️⚠️ JACKPOT CHECK — EXPLICIT ANSWER

**(a) Is there an odour-threshold list / OAV table with threshold values and stated basis?**

**NO complete table. PARTIAL, and better-evidenced than the companion.** Specifically:

| what exists | where | in this PDF? |
|---|---|---|
| Full OAV table for all 82 volatiles | **Table S5** | ❌ **online supplementary, NOT in PDF** |
| Full concentration table for all 82 volatiles | **Table S4** | ❌ **online supplementary, NOT in PDF** |
| Proximate composition of the 7 proteins | **Table S1** | ❌ online, not in PDF |
| Sensory scoring criteria | **Table S2** | ❌ online, not in PDF |
| LF-NMR relaxation/signal table | **Table S3** | ❌ online, not in PDF |
| Taste-compound concentrations | **Table S6** | ❌ online, not in PDF |
| **TAV values of 32 taste compounds** | **Table S7** | ❌ online, not in PDF |
| **ONE threshold printed verbatim in body text** | **p. 9** | ✅ **YES** |
| **Four thresholds recoverable exactly from printed OAV ranges** | **p. 9** | ✅ **YES (§4)** |
| **The stated basis of the thresholds** | **§2.11, p. 3: "air threshold value (μg/kg)"** | ✅ **YES — and it is WRONG (§4)** |

**(b) Is a source citation given for each threshold?**

**NO. Not for any threshold. Not one.** The 1.5 μg/kg for 1-octen-3-ol (p. 9) is printed
**bare, with no reference**. §2.11 (p. 3) cites nothing for the threshold set. §3.9 (p. 9)
says only "the OAV of each volatile compound was determined based on **its flavor threshold**"
— no source. **The provenance chain for every threshold in this paper is broken at the first
link.** If Table S5 is ever retrieved, checking whether *it* carries citations is the single
highest-value follow-up question.

**Verbatim, the one printed threshold (p. 9):**
> "1-Octen-3-ol was the most abundant alcohol and was found in all samples. Its content in
> ERPI (974.28 μg/kg) was significantly higher than in other samples. The lowest content was
> found in EYP (138.24 μg/kg), followed by ENPI (197.56 μg/kg) and ESPI (198.79 μg/kg). The
> entire flavor character of the PBMAs is significantly influenced by 1-octen-3-ol **because
> of its exceptionally low odor threshold (1.5 μg/kg).**"

**1.5 μg/kg is the classical 1-octen-3-ol threshold IN WATER.** The air threshold of
1-octen-3-ol is of order 10⁻² μg/L of air. **The paper prints an aqueous threshold in the body
while its Methods calls the denominator an air threshold. The contradiction is internal to the
paper and needs no external evidence.**

---

## 4. THE RECOVERED THRESHOLDS — arithmetic, both ends where available

The paper prints, p. 9 (verbatim):
> "It can be observed that among these 24 compounds, **2-Pentylfuran (OAV: 990.04-2936.74),
> 1-Octen-3-ol (OAV: 92.16-649.52), 2-Heptanone (OAV: 18.55-21.88), and Hexanal (OAV:
> 14.97-252.92)** not only showed OAV values >10, indicating significant contributions to the
> overall flavor, but also exhibited relatively high content and marked differences among
> samples, suggesting that they may be key volatile compounds in HME-PBMAs."

Cross-multiplying against the printed concentrations (**this wave's arithmetic — DERIVED, not
printed by the paper**):

| compound | printed OAV range (p. 9) | printed conc. min / max (μg/kg) | **implied threshold (μg/kg)** | check |
|---|---|---|---|---|
| **hexanal** | 14.97 – 252.92 | **74.86** (ENPI, p. 9) / **1264.62** (EPPI, p. 9) | **5.000** | 74.86/14.97 = **5.0007**; 1264.62/252.92 = **5.0001** — **BOTH ENDS EXACT** |
| **1-octen-3-ol** | 92.16 – 649.52 | **138.24** (EYP, p. 9) / **974.28** (ERPI, p. 9) | **1.500** | 138.24/92.16 = **1.50000**; 974.28/649.52 = **1.50000** — **BOTH ENDS EXACT, and 1.5 is PRINTED (p. 9)** |
| **2-pentylfuran** | 990.04 – 2936.74 | **5742.21** (ENPI, p. 13 & abstract) / **17033.09** (ERPI, p. 13) | **5.800** | 5742.21/990.04 = **5.8000**; 17033.09/2936.74 = **5.8000** — **BOTH ENDS EXACT** |
| **2-heptanone** | 18.55 – 21.88 | max **3063.41** (EMPI, p. 9); min not printed | **140.0** | 3063.41/21.88 = **140.01** — **one end only**; implies min conc ≈ **2597 μg/kg** |

**All four are the classical AQUEOUS (water) odour thresholds**: hexanal ≈5 ppb,
1-octen-3-ol ≈1.5 ppb, 2-pentylfuran ≈5.8 ppb, 2-heptanone ≈140 ppb in water. **None is an
air threshold.**

**⇒ CONFIRMED, by four independent back-solves plus one printed value: Xin et al. computed
`matrix concentration (μg/kg of a 65 %-moisture protein extrudate) ÷ aqueous odour threshold
(μg/kg of water)`, with no partition step, and labelled the denominator an air threshold.**

**⚠️ THE 5.800 CROSS-CHECK IS THE STRONGEST SINGLE PIECE OF EVIDENCE IN THIS DOSSIER.**
The companion paper (`Xin2026_extraction.md` §5) yielded **5.800 μg/kg for 2-pentylfuran**
from *entirely different* concentrations (50 230.59 / 100 259.10) and *entirely different*
OAVs (8 660.45 / 17 286.05). This paper yields **5.800** from 5 742.21 / 17 033.09 and
990.04 / 2 936.74. **Four independent numbers, two papers, two journals, one threshold, agreeing
to four significant figures.** The threshold set is a fixed, reused, uncited internal table —
and the repo now has four of its entries with certainty.

**Thresholds NOT recoverable** (no OAV range printed alongside a concentration): nonanal,
heptanal, octanal, decanal, pentanal, (E)-2-octenal, (E,E)-2,4-decadienal, hexanol,
1-heptanol, 2-butylfuran, 2-octanone, 2-nonanone, 2-decanone, 2-undecanone, 3-octanone,
3,5-octadien-2-one, limonene, 2-pentylpyridine, 2-ethyl-5-methylpyrazine, acetylthiazole
(20 of the 24 OAV > 1 compounds). **Those 20 are in Table S5.**

---

## 5. EVERY QUANTIFIED VOLATILE NUMBER IN THE PDF — transcribed, not summarised

⚠️ **There is no volatile *table* in this PDF.** Table S4 (82 compounds × 6 samples) and
Table S5 (OAVs) are online. What follows is **every quantified volatile datum the printed
paper contains** — body text (p. 9, p. 13, abstract p. 1) plus Fig. 4E digitised at 600 dpi.

### 5a. Compound inventory (p. 9, verbatim, MEASURED)

> "**82 volatile substances were found in the samples, including 12 alcohols, 12 aldehydes,
> 17 ketones, 1 ester, 1 phenol, 6 furan compounds, 2 acids, 19 alkenes and alkanes,
> 3 aromatic compounds, and 9 other compounds (Table S4).**"

12+12+17+1+1+6+2+19+3+9 = **82 ✓**

> "Among the 82 compounds, **24 volatiles had OAVs >1, including 3 alcohols, 8 aldehydes,
> 7 ketones, 2 furans, 1 alkene/alkane, and 3 other compounds (Table S5). Eleven of these
> compounds had OAVs >10, namely 1-Octen-3-ol, (E, E)-2,4-Decadienal, Hexanol, Heptanal,
> Nonanal, Octanal, Hexanal, 2-Heptanone, 2-Octanone, 2-Butylfuran, and 2-Pentylfuran.**"

3+8+7+2+1+3 = **24 ✓**; the OAV>10 list has **11 ✓**.

**The 24 OAV>1 compounds, read off Fig. 4D (Sankey source nodes, p. 10, [fig, 600 dpi]) —
this list is NOT printed in the text and is recovered here:**
- **alcohols (3):** 1-Heptanol · 1-Octen-3-ol · Hexanol
- **aldehydes (8):** (E,E)-2,4-Decadienal · (E)-2-Octenal · Heptanal · Decanal · Hexanal ·
  Nonanal · Pentanal · Octanal
- **ketones (7):** 2-Heptanone · 2-Decanone · 2-Nonanone · 2-Undecanone · 2-Octanone ·
  3,5-Octadien-2-one · 3-Octanone
- **furans (2):** 2-Butyl-Furan · 2-Pentyl-Furan
- **alkene (1):** Limonene
- **other (3):** 2-Pentyl-pyridine · **2-Ethyl-5-methyl-Pyrazine** · **Acetyl-thiazole**

Class counts match the text exactly (3/8/7/2/1/3). ⚠️ **The three "other" compounds are the
paper's only Maillard N/S-heterocycles, and the text never names them.** A pyrazine and a
thiazole are OAV > 1 in every sample and are discussed nowhere in the body.

**The 6 furans (p. 11, verbatim):** "(E)-2-(2-pentenyl) furan, 2-propylfuran, 2-butylfuran,
2-pentylfuran, 2-ethylfuran, and 2-octylfuran."

### 5b. CLASS TOTALS, μg/kg — text values + Fig. 4E digitisation

Bar order in Fig. 4E is **ESPI (blue) · ENPI (orange) · ERPI (grey) · EYP (yellow) ·
EMPI (light blue) · EPPI (green)**. Axis 0–20 000 μg/kg. Letters as printed on the figure.

| class | ESPI | ENPI | ERPI | EYP | EMPI | EPPI |
|---|---:|---:|---:|---:|---:|---:|
| **Alcohols** | ≈450 **[fig]** `d` | **701.72** (p. 9) `b` | **1252.64** (p. 9) `a` | ≈420 **[fig]** `d` | **715.23** (p. 9) `b` | ≈540 **[fig]** `c` |
| **Aldehydes** | ≈1140 **[fig]** `c` | **193.95** (p. 9) `d` | **2398.22** (p. 9) `a` | **1868.31** (p. 9) `b` | **1746.01** (p. 9) `b` | **1871.78** (p. 9) `b` |
| **Ketones** | ≈3930 **[fig]** `c` | **3079.28** (p. 9) `f` | **4127.31** (p. 9) `b` | **3648.19** (p. 9) `d` | **4398.58** (p. 9) `a` | ≈3360 **[fig]** `e` |
| **Furans** | ≈6750 **[fig]** `e` | ≈6600 **[fig]** `e` | **≈18 880 [fig]** `a` | ≈14 960 **[fig]** `b` | ≈8 120 **[fig]** `c` | ≈7 940 **[fig]** `d` |
| **Esters** | ≈0 **[fig]** `a` | ≈0 `a` | ≈0 `a` | ≈0 `a` | ≈0 `a` | ≈0 `a` |
| **Phenols** | ≈20 **[fig]** `b` | ≈3 **[fig]** `b` | ≈83 **[fig]** `a` | ≈22 **[fig]** `b` | ≈25 **[fig]** `b` | ≈62 **[fig]** `a` |
| **Alkanes & Alkenes** | ≈172 **[fig]** `e` | ≈165 **[fig]** `e` | ≈300 **[fig]** `a` | ≈185 **[fig]** `d` | ≈240 **[fig]** `c` | ≈262 **[fig]** `b` |
| **Aromatic compounds** | ≈48 **[fig]** `e` | **≈470 [fig]** `a` | ≈180 **[fig]** `d` | ≈275 **[fig]** `c` | ≈252 **[fig]** `c` | ≈350 **[fig]** `b` |
| **Acids** | ≈0 **[fig]** | ≈0 | ≈25 **[fig]** `a` | ≈5 **[fig]** `b` | ≈0 | ≈0 |
| **Others** | ≈297 **[fig]** `b` | ≈180 **[fig]** `e` | **≈650 [fig]** `a` | ≈272 **[fig]** `c` | ≈305 **[fig]** `b` | ≈230 **[fig]** `d` |

**[fig] values digitised from `k3img/fig4E_left-10.png` and `k3img/fig4E_inset-10.png`
(600 dpi crops of p. 10). Accuracy ±3 % on the main panel, ±5 % on the 0–600 inset.**

⚠️ **ENPI is the LOWEST arm in six classes but the HIGHEST in aromatic compounds (≈470 μg/kg,
group `a`, ~2.8× the next-lowest).** The text never mentions this. **The one class ENPI leads
is the one the paper's "NPI binds off-flavours" story cannot explain**, since aromatic
compounds would bind at least as well as the aliphatic markers.

### 5c. Class-total text quotes, verbatim (p. 9)

> "Alcohols are created when peroxides are broken down, and **12 alcohols in all were found
> across all samples. The alcohol content was greatest in RPI (1252.64 μg/kg), followed by
> EPPI (1871.78 μg/kg), EMPI (715.23 μg/kg), and ENPI (701.72 μg/kg).**"

⚠️ **ARITHMETICALLY IMPOSSIBLE AS PRINTED: 1871.78 > 1252.64, so "greatest … followed by" is
false.** **1871.78 μg/kg is the EPPI *aldehyde* total** (quoted two paragraphs later) —
copy-pasted into the alcohol sentence. **Fig. 4E puts EPPI alcohols at ≈540 μg/kg, group `c`**,
between EMPI/ENPI (`b`, ~700–715) and ESPI/EYP (`d`, ~420–450) — which is **also inconsistent
with the letters**, since `c` should sit *below* `b`. Reading the letters as authoritative:
**ERPI (a) 1252.64 > ENPI (b) 701.72 ≈ EMPI (b) 715.23 > EPPI (c) ≈540 > ESPI (d) ≈450 ≈
EYP (d) ≈420.** **DO NOT INGEST 1871.78 as an alcohol value.** ("RPI" should read "ERPI".)

> "Among the samples, **ERPI exhibited the highest aldehyde content (2398.22 μg/kg), followed
> by EPPI (1871.78 μg/kg), EYP (1868.31 μg/kg), and EMPI (1746.01 μg/kg). Notably, the
> aldehyde content in ENPI was only 193.95 μg/kg, substantially lower than in the other
> samples.**"

> "**Ketone levels were relatively high across all samples, with the highest content observed
> in EMPI (4398.58 μg/kg), followed by ERPI (4127.31 μg/kg) and EYP (3648.19 μg/kg). The
> lowest ketone content was found in ENPI (3079.28 μg/kg).**"

⚠️ **The ketone sentence skips ESPI, which Fig. 4E ranks `c` — i.e. ABOVE EYP (`d`).** The
text's "followed by … EYP" is therefore wrong in rank order. ESPI ketones ≈3930 μg/kg **[fig]**.

### 5d. INDIVIDUAL COMPOUNDS, μg/kg — every printed value

**All MEASURED (semi-quantitative, relative to cyclohexanone).**

| compound | ESPI | ENPI | ERPI | EYP | EMPI | EPPI | anchor |
|---|---:|---:|---:|---:|---:|---:|---|
| **1-octen-3-ol** | **198.79** | **197.56** | **974.28** (max) | **138.24** (min) | — | — | p. 9; ERPI also p. 13 |
| **hexanal** | — | **74.86** (min) | — | — | — | **1264.62** (max) | p. 9; ENPI also abstract p. 1 |
| **nonanal** | — | **12.44** (min) | **195.04** (max) | **119.01** | — | **110.78** | p. 9; ENPI also abstract p. 1; ERPI also p. 13 |
| **2-heptanone** | **2969.74** | **2778.1** | **3004.15** | — | **3063.41** (max) | — | p. 9 (ESPI/ERPI/EMPI); **abstract p. 1** (ENPI) |
| **2-pentylfuran** | — | **5742.21** (min) | **17033.09** (max) | — | — | — | **abstract p. 1** (ENPI); **p. 13** (both) |

**That is the complete set of individually quantified volatiles in the printed paper: five
compounds, 17 numbers.** Everything else is in Table S4.

**Derived, this wave (DERIVED, not printed):**
- hexanal **EPPI / ENPI = 1264.62 / 74.86 = 16.9×**
- nonanal **ERPI / ENPI = 195.04 / 12.44 = 15.7×**
- total aldehydes **ERPI / ENPI = 2398.22 / 193.95 = 12.4×**
- 2-pentylfuran **ERPI / ENPI = 17033.09 / 5742.21 = 2.97×**
- 2-heptanone spread **EMPI / ENPI = 3063.41 / 2778.1 = 1.10×** — **nearly flat**
- 1-octen-3-ol **ERPI / EYP = 974.28 / 138.24 = 7.05×**
- **2-pentylfuran is 90.2 % of ERPI's furan class** (17033.09 / ≈18 880 **[fig]**) and
  **87.0 % of ENPI's** (5742.21 / ≈6 600 **[fig]**)
- **2-heptanone is 90.2 % of ENPI's ketone class** (2778.1 / 3079.28) and
  **69.6 % of EMPI's** (3063.41 / 4398.58)

⚠️ **The 2-heptanone spread (1.10×) is the flattest of any key compound, while aldehydes span
12.4×.** Protein source moves lipid-*aldehyde* chemistry hard and methyl-*ketone* chemistry
barely at all. The paper does not notice. **This is a clean, usable dissociation.**

### 5e. Abstract-level summary of ENPI (p. 1, verbatim)

> "Volatile compound analysis indicated that **ENPI contained substantially lower levels of
> key off-flavor compounds, including hexanal (74.86 μg/kg), nonanal (12.44 μg/kg),
> 2-heptanone (2778.1 μg/kg), and 2-pentylfuran (5742.21 μg/kg)**, compared with other
> formulations."

⚠️ **"2-heptanone (2778.1 μg/kg)" is presented as "substantially lower" but is only 9.3 %
below the maximum (3063.41).** The abstract's framing overstates the effect for that compound.

### 5f. Conclusion-level ERPI values (p. 13, verbatim)

> "The ERPI blend, likely limited by the relatively low water absorption capacity of RPI,
> exhibited higher hardness and relatively higher levels of several volatile off-flavor
> compounds (**e.g., 1-octen-3-ol: 974.28 μg/kg, nonanal: 195.04 μg/kg, 2-pentylfuran:
> 17033.09 μg/kg**)."

### 5g. VIP scores — Fig. 4B (p. 10), digitised at 600 dpi **[fig]**

Axis 1.4 → 2.4. Twenty compounds have VIP > 1 in the figure. **Values ±0.02.**

| compound | VIP **[fig]** | | compound | VIP **[fig]** |
|---|---:|---|---|---:|
| 2-Undecanone | **2.40** | | 1,3-Octadiene | 1.56 |
| **2-Acetyl-thiazole** | **2.32** | | Indole | 1.56 |
| 2-Octyl-Furan | **2.26** | | Propiophenone | 1.53 |
| 3-Octanone | 1.96 | | Hexanal | 1.49 |
| 2,6-Dimethyl-Undecane | 1.94 | | Pentanal | 1.45 |
| **2-Methyl-pyrazine** | **1.91** | | Butyldihydro-2(3H)-Furanone | 1.44 |
| Hexanol | 1.86 | | Nonanal | 1.43 |
| (E,E)-2,4-Decadienal | 1.67 | | 5-Methyl-Tridecane | 1.41 |
| 2-Methyl-3-Octanol | 1.64 | | **Dimethyl phthalate** | 1.38 |
| Benzyl alcohol | 1.61 | | | |
| Nonanoic acid | 1.56 | | | |

**The text's VIP list (p. 9, verbatim):**
> "**Ketones (2-Undecanone, 3-Octanone, Propiophenone, 5-Butyldihydro-2(3H)-furanone),
> alcohols (Hexanol, 2-Methyl-3-octanol, Benzyl alcohol), aldehydes ((E, E)-2,4-Decadienal,
> Pentanal, Hexanal, Nonanal), and alkene/alkane compounds (2,6-Dimethylundecane,
> 1,3-Octadiene, 5-Methyltridecane) are the main compounds with VIP scores >1**"

⚠️ **The text names 14 of the 20. The six it omits are 2-acetyl-thiazole (VIP 2.32, rank 2),
2-octyl-furan (2.26, rank 3), 2-methyl-pyrazine (1.91, rank 6), nonanoic acid, indole, and
dimethyl phthalate.** **The two highest-ranked Maillard markers in the entire dataset — a
thiazole at rank 2 and a pyrazine at rank 6 — are silently dropped from the narrative.**
⚠️ **Dimethyl phthalate (VIP 1.38) is a plasticiser, not a food volatile** — a contamination
marker discussed nowhere.

**PLS-DA variance (p. 9, verbatim):** "With a **cumulative contribution rate of 56.6 %**, the
first two main components explained **29.1 % and 27.5 %** of the total variance."
(Fig. 4A axis labels confirm: **PC 1 (29.1 %)**, **PC 2 (27.5 %)** **[fig]**.)
⚠️ 29.1 + 27.5 = **56.6 ✓**. **No R²/Q² and no permutation test are reported** — the PLS-DA is
unvalidated. (Companion: 27.8 % + 27.5 % = 55.3 %, likewise unvalidated.)

**Fig. 4C** is a circular heatmap of all 82 volatiles — **z-scored "Expression" −2…+3, no
concentrations recoverable**. **Fig. 4D** Sankey line widths are relative content, **not
quantified**.

### 5h. E-nose (§3.8, p. 9)

Ten sensors; response read at a **unified 90 s** after stabilisation; **3 g homogenised
sample in a 20 mL headspace vial, incubated 40 °C for 30 min**; injection **200 mL/min**,
pre-run 5 s, measurement 120 s, sampling 1.0 s, wash 120 s (§2.10, p. 3).
Radar axes (Fig. 3E **[fig]**): **W1C, W5S, W3C, W6S, W5C, W1S, W1W, W2S, W2W, W3S**, scale
**0–3.0**, all responses ≤ ~2.3.
**PCA (Fig. 3F): PC1 92.9 %, PC2 3.9 %, cumulative 96.8 %** (p. 9, confirmed on axis labels
**[fig]**).
> "**EMPI and ERPI exhibited notably higher responses in W5S, W1W, and W2W sensors compared to
> other formulations, indicating increased nitrogen oxides and sulfur-containing compounds in
> these samples. ESPI and EPPI exhibited elevated response values on the W1S sensor …
> ENPI demonstrated the lowest response values across most sensors.**"

⚠️ **The sulfur claim is unsupported by the GC-MS**: the only S-containing volatile in the
whole 82-compound set is **acetyl-thiazole** (Fig. 4D), and no S-compound concentration is
reported. **E-nose sensor semantics ≠ identified compounds. Directional at best.**

---

## 6. ⚠️⚠️ THE CROSS-PAPER REPRODUCIBILITY TEST — the most important thing in this dossier

**This paper's ENPI and the companion's control are the same nominal material.**

| | **Xin 2026b, ENPI** | **Xin 2026 (Food Hydrocolloids 182, 113124), Control** |
|---|---|---|
| protein blend | **70 SPI + 20 NPI + 10 WG** (Table 1, p. 2) | **NPI:SPI:WG = 20:70:10** (companion p. 2) — **identical** |
| added carbohydrate | none | **none** (the no-carbohydrate control) |
| extruder | Hunan Fumaco twin-screw, 704 mm × 22 mm, L/D 32:1, cooling die | **identical** |
| zone temps | 60/80/110/135/140/140/140 °C | **identical** |
| screw / feed / moisture / die | 240 rpm, 2 kg/h, ≈65 %, 70 °C | **identical** |
| SPME | 2 g + 20 μL of 40 μg/mL cyclohexanone, 20 mL vial | **identical** |
| GC-MS | 5977A-7890B, 45 °C 3 min → 5 °C/min → 250 °C 10 min, EI 70 eV, 35–500 m/z | **identical** |
| authors | Xin, Cui, Liu, Nie, Huang, Qin | **Xin, Liu, Cui, Nie, Huang, Qin — same six** |

**The reported concentrations (DERIVED comparison, this wave):**

| compound | **Xin 2026b ENPI** | **Xin 2026 control** | **discrepancy** | ranges overlap? |
|---|---:|---:|---:|---|
| **hexanal** | **74.86** (p. 9) | ≈750 **[fig]** (companion Fig. 7B) | **10.0×** | 2026b: 74.86–1264.62; 2026: 453–1099 → **partial** |
| **nonanal** | **12.44** (p. 9) | ≈290 **[fig]** (companion Fig. 7B) | **23.3×** | 2026b: 12.44–195.04; 2026: ≈292–484 → **NO OVERLAP** |
| **2-pentylfuran** | **5742.21** (p. 13) | ≈97 000 **[fig]** (companion Fig. 7B) | **16.9×** | 2026b: 5742–17 033; 2026: 50 231–100 259 → **NO OVERLAP** |
| **2-heptanone** | **2778.1** (abstract p. 1) | ≈2 700 **[fig]** (companion Fig. 7B) | **1.03×** | overlapping |

**Geometric-midpoint comparison across the FULL sample sets (this wave's arithmetic):**
- 2-pentylfuran: 2026b √(5742.21 × 17033.09) = **9 890**; 2026 √(50 230.59 × 100 259.10) =
  **70 970** → **7.2× apart, zero overlap.**
- hexanal: 2026b √(74.86 × 1264.62) = **308**; 2026 √(453 × 1099) = **706** → **2.3× apart.**

**⇒ WHAT THIS ESTABLISHES.** Two papers, one lab, one extruder, one protocol, one internal
standard, six shared authors, the same nominal formulation — and the reported μg/kg for the
repo's two most load-bearing markers (hexanal, 2-pentylfuran) differ by **one order of
magnitude**, with **non-overlapping ranges across the entire sample sets**. **No physical
difference in the described systems can account for it.** The most likely cause is the
**unspecified SPME fibre / extraction conditions in this paper versus DVB/CAR/PDMS in the
companion**, i.e. the numbers are dominated by the extraction partition, not by the sample.

**This is the cleanest available empirical demonstration that HS-SPME single-IS "μg/kg" values
are not absolute concentrations and must never be ingested as such.** It is stronger than the
companion's 100 ppm 2-pentylfuran sanity check, because it is **internal to the source
literature and requires no external prior.**

⚠️ **Corollary: 2-heptanone agrees to 3 %.** So the discrepancy is compound-specific, tracking
volatility/hydrophobicity — exactly the signature of a fibre/partition change, and further
evidence against a real chemical difference.

---

## 7. TEXTURE, STRUCTURE, WATER, COLOUR — the non-volatile numbers

### 7a. Texture, Fig. 1 (p. 5) — MEASURED, ≥6 replicates per sample (§2.4, p. 3)

Method verbatim (§2.4, p. 3): TPA with **P/50 probe**, **trigger 5 g**, **pre-test 2 mm/s,
test 1 mm/s, post-test 2 mm/s**, **3 s gap between two compressions**, **50 % compression**.
Fibrous degree: **HDP/BS probe**, transverse (FV) and longitudinal (FL) shear to **75 % of
initial thickness at 1 mm/s**, **fibrous degree = FV/FL**, **≥6 measurements per sample**.

| sample | hardness (g) | springiness | chewiness | **fibrous degree** |
|---|---:|---:|---:|---:|
| **ESPI** | ≈23 200 **[fig]** `c` | ≈0.98 **[fig]** `a` | ≈19 700 **[fig]** `d` | ≈1.05 **[fig]** `c` |
| **ENPI** | **≈17 000 [fig]** `e` (lowest) | ≈0.98 **[fig]** `a` | **≈13 700 [fig]** `e` (lowest) | **1.59** (p. 5, text) `a` (highest) |
| **ERPI** | ≈27 100 **[fig]** `b` | ≈0.96 **[fig]** `a` | ≈24 100 **[fig]** `b` | ≈1.06 **[fig]** `c` |
| **EYP** | **≈35 500 [fig]** `a` (highest) | ≈0.96 **[fig]** `a` | **≈33 100 [fig]** `a` (highest) | ≈1.04 **[fig]** `c` |
| **EMPI** | ≈23 900 **[fig]** `c` | ≈0.98 **[fig]** `a` | ≈21 500 **[fig]** `c` | **1.27** (p. 5, text) `b` |
| **EPPI** | ≈21 800 **[fig]** `d` | ≈0.92 **[fig]** `b` | ≈18 700 **[fig]** `d` | **1.24** (p. 5, text) `b` |

**[fig] = digitised from `k3img/fig1_top-05.png`, 600 dpi crop of p. 5. Hardness/chewiness
±3 %; springiness ±0.02; fibrous degree ±0.03.** Only **1.59 / 1.27 / 1.24** are printed
(p. 5, verbatim): "**ENPI exhibited the highest fibrous degree (1.59), followed by EMPI (1.27)
and EPPI (1.24).**"

**Derived (this wave):** **EYP / ENPI hardness = 2.09×**, **chewiness = 2.42×**,
**fibrous degree 1.59 / 1.04 = 1.53×**. Springiness is essentially flat (0.92–0.98) and carries
no information.

**PCA on texture vs the 7 grilled-fish references (Fig. 1E, p. 5, verbatim):**
> "The top two main components explained **94.4 % and 3.5 %** of the total variance … The PCA
> score plot indicated that **ESPI, EPPI, and EMPI clustered closely and partially overlapped
> with grilled cod fillet** … **ENPI displayed distinct textural characteristics … and
> partially overlapped with grilled eel, grilled squid, and grilled black scraper.**"

⚠️ **PC1 = 94.4 % means the PCA is essentially a one-dimensional ranking by hardness/chewiness**
(the two variables with by far the largest numeric range). "Overlaps with grilled eel" reduces
to "has a similar hardness". **Directional only; the multivariate framing is decorative.**
⚠️ Fig. 1E confirms **EYP sits at the extreme left of PC1, alone**, and **ENPI at the far
right with Squid and Black Scraper** **[fig]**. **No seafood texture values are printed
anywhere.**

### 7b. Protein secondary structure, Fig. 3A (p. 7) — FTIR amide I, deconvolution + Gaussian fit

Method (§2.7, p. 3): **Spectrum II PerkinElmer**, freeze-dried sample **1 mg** in **100 mg
KBr**, **4000–400 cm⁻¹, 64 scans, 4 cm⁻¹ resolution**, **Omnic 9.2 + Peakfit 4.12**, baseline
correction / smoothing / Fourier deconvolution / **Gaussian peak fitting**.

**FITTED (deconvolution output), not directly measured.** Read at 600 dpi from
`k3img/fig3_AB-07.png`; **the percentages are printed on the bars, so these are exact, not
estimated**:

| sample | **α-helix** | **β-sheet** | **β-turn** | **random coil** | Σ |
|---|---:|---:|---:|---:|---:|
| **ESPI** | **19 %** `a` | **20 %** `c` (lowest) | **40 %** `a` (highest) | 21 % `bc` | 100 ✓ |
| **ENPI** | 16 % `b` (lowest) | **27 %** `b` | 35 % `b` | 22 % `ab` | 100 ✓ |
| **ERPI** | 18 % `ab` | 22 % `c` | 36 % `b` | **24 %** `a` | 100 ✓ |
| **EYP** | 18 % `ab` | 22 % `c` | 37 % `ab` | 23 % `ab` | 100 ✓ |
| **EMPI** | 17 % `ab` | **32 %** `a` (highest) | **29 %** `c` (lowest) | 22 % `ab` | 100 ✓ |
| **EPPI** | 18 % `ab` | **31 %** `a` | 31 % `c` | 20 % `c` | 100 ✓ |

Text confirmation (p. 7–8, verbatim): "**ESPI has the most α-helix (19 %) and β-turn (40 %),
and the lowest β-sheet content (20 %)** … Among the six samples, **EMPI demonstrated the
highest β-sheet content (32 %), followed by EPPI (31 %) and ENPI (27 %)**."

⚠️ **THE ABSTRACT OVERSELLS ENPI.** Abstract (p. 1): "**The NPI-containing extrudate (ENPI)
exhibited a relatively high β-sheet content (27 %) and the highest fibrous degree (1.59)**".
**ENPI is THIRD on β-sheet (27 %) behind EMPI (32 %) and EPPI (31 %) — yet it has by far the
highest fibrous degree (1.59 vs 1.27 and 1.24).** **The paper's own headline mechanism
(β-turn → β-sheet conversion drives fibrousness) is contradicted by its own ranking:**
β-sheet ranks EMPI > EPPI > ENPI while fibrous degree ranks ENPI ≫ EMPI ≈ EPPI. §3.12 (p. 13)
nonetheless reports "**sensory attributes related to organization and mouthfeel showed strong
positive correlations with fibrous degree and β-sheet content**" — the correlation matrix is
**Fig. S13, online, with no coefficients printed**. **Record the β-sheet→fibrousness claim as
UNSUPPORTED by the printed data.**

### 7c. LF-NMR water distribution (§3.6, p. 8) — MEASURED, CPMG

| quantity | value | sample | anchor |
|---|---:|---|---|
| **T₂₁** | **19.34 ms** | ERPI | p. 8 |
| **T₂₁** | **16.83 ms** | EYP | p. 8 |
| **T₂₂** | **652.08 ms** | EMPI | p. 8 |
| **T₂₂** | **1367.45 ms** | ENPI (longest) | p. 8 |
| **M₂b** | **8.72 %** | EYP (highest bound-water fraction) | p. 8 |

> "**a T₂₂ peak was detected only in these two samples (EMPI: 652.08 ms; ENPI: 1367.45 ms)**"
Fig. 3B confirms: intensity axis 0–800 with a break at ~40/190; **T₂b ≈0.3–1.5 ms (intensity
≤40), T₂₁ ≈10–40 ms (intensity ≈600–660), T₂₂ ≈500–1500 ms (intensity ≈5, EMPI/ENPI only)**
**[fig, 400 dpi]**. **Full table is Table S3, online — the other four samples' T₂₁ and all
M-fractions are NOT in this PDF.**

### 7d. Water absorption capacity, Fig. 3D (p. 7–8) — MEASURED

Method (§2.9, p. 3): dried **40 °C / 24 h**, rehydrated in **40 mL deionised water, 50 mL
tube, 50 °C water bath, 16 h**, drained **5 min** at 26 °C / 50 ± 5 % RH.
WAC (%) = (W2 − W1)/W1 × 100.

| sample | WAC (%) | letter |
|---|---:|---|
| **ENPI** | **565.07** (p. 8, highest) | `a` |
| **EMPI** | **546.39** (p. 8) | `b` |
| **ESPI** | **545.51** (p. 8) | `b` |
| **EPPI** | ≈433 **[fig, 400 dpi]** | `c` |
| **ERPI** | ≈410 **[fig, 400 dpi]** | `d` |
| **EYP** | **359.12** (p. 8, lowest) | `e` |

> "**EYP showed the lowest water absorption rate (359.12 %), approximately 36.5 % lower than
> that of ENPI.**" ✓ check: 1 − 359.12/565.07 = **36.4 %** — consistent.

### 7e. Water loss (§3.7, p. 8 + Fig. 3C)

> "**During the initial hour, samples exhibited approximately 3 % weight loss due to water
> evaporation, which increased to 5 % cumulative loss by the fourth hour.**"
Conditions (§2.9, p. 3): pieces **20 × 20 × 10 mm**, **ten pieces per sample**, stored **26 °C,
50 ± 5 % RH**, weighed every **10 min for 60 min** then every **30 min for 180 min**.

⚠️ **THREE PROBLEMS WITH FIG. 3C.**
1. **The x-axis is labelled "Time (s)" running 0–240, but the protocol is 240 MINUTES.** Unit
   error in the figure.
2. **The y-axis is labelled "Rate of water loss (%)" but runs 90 → 100 DESCENDING** — it is
   **mass retention**, not loss.
3. **The plotted endpoint is ≈92.2–93.3 % at t = 240 [fig, 400 dpi], i.e. 6.7–7.8 % loss —
   not the 5 % the text states.** **The text and the figure disagree by ~50 %.**
   Ranking from the figure **[fig]**: **EMPI lowest curve (≈92.2 %, most loss) < ENPI ≈92.4 <
   ESPI/ERPI ≈92.8 < EPPI ≈93.0 < EYP ≈93.3 (least loss)** — consistent with the text's
   "EMPI and ENPI … significantly higher moisture loss rates … EYP exhibited the slowest water
   loss". **Ingest the RANKING; do not ingest 3 % / 5 %.**
   Eq. (2) as printed is **"Water loss rate (%) = Wb/Wa × 100 %"** — that is the *retention*
   ratio, not the loss. **The equation as printed is wrong for its own label**, which explains
   the axis.

### 7f. Colour (Table 2, p. 5) — MEASURED, transcribed verbatim

**UltraScan Pro HunterLab; 20 × 20 × 10 mm cubes; calibrated against a standard white plate;
8 measurements per sample; outliers > 2 SD excluded (§2.5, p. 3). ΔE is vs the WHITE PLATE,
not vs a control** (Eq. 1, p. 3; confirmed p. 5: "ESPI and ENPI exhibited the smallest color
differences … indicating their color was closest to the standard white plate").

| Sample | L* | a* | b* | ΔE |
|---|---|---|---|---|
| **ESPI** | **56.36 ± 0.21** `a` | 0.19 ± 0.03 `c` | 8.55 ± 0.08 `bc` | **43.99 ± 0.22** `c` |
| **ENPI** | **56.20 ± 0.61** `a` | 0.17 ± 0.03 `c` | 8.68 ± 1.28 `bc` | **44.19 ± 0.60** `c` |
| **ERPI** | **53.09 ± 0.62** `c` | 0.84 ± 0.09 `ab` | 8.14 ± 0.38 `c` | **47.13 ± 0.57** `a` |
| **EYP** | 53.36 ± 0.30 `cd` | **0.91 ± 0.06** `a` | 9.43 ± 0.36 `ab` | **47.12 ± 0.24** `a` |
| **EMPI** | 55.02 ± 0.48 `b` | **0.06 ± 0.00** `d` | **9.72 ± 0.14** `a` | 45.54 ± 0.50 `b` |
| **EPPI** | 53.89 ± 0.08 `d` | 0.79 ± 0.06 `b` | 9.04 ± 0.29 `abc` | 46.52 ± 0.09 `a` |

> Note (verbatim): "Comparisons were carried out between values of the same column; values with
> different letter(s) indicate a significant difference at p ≤ 0.05."

⚠️ **THE LETTERING IN THE L* COLUMN IS INTERNALLY INCONSISTENT AS PRINTED.**
ERPI = 53.09 `c`, EYP = 53.36 `cd`, **EPPI = 53.89 `d`**, EMPI = 55.02 `b`. A `d` group at
53.89 that excludes 53.09 (`c`) while 53.36 is `cd` is a coherent chain — **but EPPI at 53.89
being significantly different from ERPI at 53.09 while EYP at 53.36 differs from neither, with
SDs of 0.08–0.62, is arithmetically implausible under Tukey HSD.** **Treat the L* letters as
unreliable.**
⚠️ **ΔE is measured against a white plate, so it is NOT a browning index relative to a control**
— it is dominated by L*. **DO NOT use these ΔE values as Maillard-browning extents.** (Contrast
the companion, where ΔE has the same white-plate basis and was likewise over-read.)
✅ **The usable colour signal is a*:** EYP 0.91 ≈ ERPI 0.84 ≈ EPPI 0.79 **≫** ESPI 0.19 ≈
ENPI 0.17 **>** EMPI 0.06. **Redness tracks the off-flavour ranking** (EYP/ERPI highest on
both), which is at least consistent with more Maillard/oxidation in those arms.

### 7g. Microstructure (Fig. 2, p. 6) — QUALITATIVE ONLY

SEM **S8020 Hitachi**, freeze-dried, gold-sputtered, **×40, ×200, ×2000** (§2.6, p. 3).
**No pore sizes, no pore counts, no fibre widths are quantified anywhere.** The claims
("EYP denser", "ENPI pores uniformly arranged … distinct layered fibrous network",
"ERPI/EYP/ESPI larger and fewer pores") are **eyeball descriptions**. **NOT INGESTABLE as
numbers.**

### 7h. Rheology (§3.1, p. 4 + Fig. S1) — NO NUMBERS IN THE PDF

**DHR-2 TA Instruments, 40 mm parallel plate, pastes at 25 % (w/w) moisture, mixed 600 r/min
for 1 h, held 4 °C 12 h, frequency sweep at 25 °C over 0.1–100 rad/s in the LVR** (§2.2, p. 2).
Only an ordinal claim is printed: "**NPIP exhibits similar and higher storage modulus than
PPIP, followed by MPIP, whereas YPP and RPIP displays relatively lower storage modulus.**"
⚠️ **Fig. S1 is online. No G′ or G″ value appears in this PDF.** Also ⚠️ **the rheology is run
at 25 % moisture and 25 °C — nowhere near the 65 % / 140 °C extrusion state**, so it is a
material-screening proxy, not a process measurement.

---

## 8. NON-VOLATILE TASTE COMPOUNDS (§3.10, p. 11–13) — MEASURED, HILIC-MS/MS

Method (§2.13, p. 3): **AB Sciex 5500 QTrap + Shimadzu LC-30 CE**, **InfinityLab Poroshell 120
HILIC-Z, 2.1 × 150 mm, 2.7 μm** + guard **2.1 × 5 mm**; **MRM**; curtain gas 35 psi;
**TEM 600 °C**; GS1 = GS2 = 60 psi; positive mode **EP 10 V, IS 5500 V**; negative mode
**EP −10 V, IS −4500 V**. Method from the group's own **Xin et al. 2023** (validated
HILIC-MS/MS). **This is genuine absolute quantification with authentic standards — unlike the
volatiles.**

**Inventory (p. 11):** "**A total of 32 key taste-active compounds were detected across all
samples, primarily comprising 7 organic acids, 18 free amino acids, 6 nucleotide metabolites,
and 1 organic base (Table S6).**" 7+18+6+1 = **32 ✓**

**PLS-DA (p. 11):** "**With a cumulative contribution rate of 60.2 %, the first two main
components explained 38.5 % and 21.7 % of the variance.**" (38.5 + 21.7 = 60.2 ✓)
**15 VIP > 1 (Fig. 5C, p. 12), named in text (p. 11):** 3 organic acids (**cis-aconitic acid,
pyroglutamic acid, pyruvic acid**), 1 organic base (**betaine**), and 11 FAAs (**Arg, Gly, Trp,
Thr, His, Asp, Ser, Pro, Tyr, Lys, Phe**). 3+1+11 = **15 ✓**

**Every printed taste number (units μg/g — note: NOT μg/kg):**

| quantity | value | sample | anchor |
|---|---:|---|---|
| **total FAAs** | **3347.27 μg/g** (highest) | EPPI | p. 11 |
| total FAAs | **2632.81 μg/g** | EYP | p. 11 |
| total FAAs | **1940.65 μg/g** | ERPI | p. 11 |
| **total FAAs** | **1613.74 μg/g** (lowest) | ESPI | p. 11 |
| **Asp** | **119.45 μg/g** | EMPI | p. 13 |
| **Asp** | **98.22 μg/g** | ENPI | p. 13 |
| **Arg** | **982.23 μg/g** (highest FAA overall) | EPPI | p. 13 |
| **Arg** | **903.10 μg/g** | EYP | p. 13 |
| **Arg** | **273.33 μg/g** | ENPI | p. 13 |
| **lactic acid** | **290.43 – 747.88 μg/g** (range across samples) | — | p. 12 |
| **malic acid** | **65.10 – 271.57 μg/g** (range across samples) | — | p. 12 |

⚠️ **ENPI and EMPI total FAAs are NOT printed** (only EPPI, EYP, ERPI, ESPI). ⚠️ **The
predominant FAAs (p. 13): "Arg, Ser, Lys, His, Gly, and Glu"; highest-TAV FAAs: "Asp, Arg, His,
Glu, and Lys" — but no TAV VALUES are printed; they are Table S7, online.**
⚠️ **No taste thresholds are printed either** — the TAV denominators are as unrecoverable as
the OAV denominators, and here there is no printed TAV range to back-solve from.

**Directly relevant to the repo:** **total FAAs span 1613.74 (ESPI) → 3347.27 (EPPI) μg/g,
a 2.07× range in Maillard amine substrate at constant total protein and constant process.**
That is a genuine, absolutely-quantified precursor-loading benchmark — **the best-quantified
number in the entire paper** — and the paper itself makes the link (p. 13): "**FAAs participate
in the Maillard reaction with reducing sugars.**" ⚠️ **But there is no reducing-sugar
measurement, so the pairing is one-sided.**

---

## 9. SENSORY (§2.14 p. 4, §3.11 p. 13) — NO NUMBERS AT ALL

**10 trained evaluators, aged 20–30, 5 F / 5 M; samples cut to 30 × 30 × 10 mm on disposable
plates; six attributes — colour, appearance, organization, mouthfeel, taste, flavor — scored
1–10; criteria in Table S2 (online). Ethics: institutional approval not required; written
informed consent obtained.**

⚠️ **Every sensory score is in Fig. S12, ONLINE. Not one sensory number appears in this PDF.**
The §3.11 discussion is entirely ordinal: "**EYP received the lowest scores for most
attributes** … **ENPI exhibited overall superior sensory characteristics** … **EMPI and EPPI
exhibited intermediate sensory acceptance** … **ESPI scored relatively high in color,
appearance, and taste but showed lower flavor acceptance** … **ERPI exhibited generally lower
sensory scores**."
⚠️ **n = 10, and unlike the companion (n = 12) the authors do NOT flag the panel size as a
limitation.**

**⇒ SENSORY: DIRECTIONAL RANKING ONLY — ENPI > {EMPI, EPPI, ESPI} > ERPI > EYP.**

---

## 10. MOLECULAR DOCKING (§2.12 p. 3, §3.9 p. 11) — COMPUTATIONAL, NO NUMBERS PRINTED

**LibDock module, Discovery Studio 2019**; proteins **SPI, WG, NPI, RPI, YP, PPI, MPI** from
**PDB**, waters removed, H added, protonation optimised, sites from **Site Finder**; ligands
**nonanal, 2-pentylfuran, hexanol, 1-octen-3-ol, 2-heptanone, heptanal** from **PubChem**,
geometry-optimised with Prepare Ligands; **default LibDock parameters**.
> "**Compared with the other four alternative plant proteins, NPI exhibited higher LibDock
> Scores for these off-flavor compounds (Fig. S4), suggesting a relatively stronger potential
> binding capacity, which may promote their retention or reduce their release.**"

⚠️ **ALL LibDock Scores are in Figs. S4–S11, ONLINE. Not one score is printed.**
⚠️ **"the other four alternative plant proteins" — there are FIVE alternatives (RPI, YP, MPI,
PPI, NPI ⇒ four others besides NPI). As written it is ambiguous but arithmetically it must
mean RPI/YP/MPI/PPI.**
⚠️ ⚠️ **THE DOCKING ARGUMENT IS DIRECTIONALLY BACKWARD FOR THE PAPER'S OWN CONCLUSION.** The
paper argues NPI binds off-flavours *more strongly*, which would **retain** them in the matrix
— **and then reports ENPI as having the LOWEST measured headspace levels.** Stronger binding
lowers *headspace* concentration but **raises** the amount actually present in the food; on an
HS-SPME measurement the two are confounded. **So the docking result, if true, is an alternative
explanation for ENPI's low numbers that has nothing to do with less off-flavour being formed.**
The paper half-notices this ("may promote their retention or reduce their release") and then
proceeds to treat ENPI as genuinely lower in the abstract and conclusions. **A crisp,
recordable confound: HS-SPME cannot distinguish "less formed" from "more bound".**

---

## 11. ⚠️ EVERY DEFECT AND UNREADABLE ITEM, CONSOLIDATED

### Unreadable / absent from this PDF (all supplementary, online only)
| item | content | impact |
|---|---|---|
| **Table S1** | proximate composition of all 7 proteins | **fatal for substrate normalisation** |
| **Table S2** | sensory scoring criteria | minor |
| **Table S3** | LF-NMR relaxation times + signal fractions, all 6 samples | moderate |
| **Table S4** | **82 volatiles × 6 samples, μg/kg** | **major — 492 cells** |
| **Table S5** | **OAVs of the 24 OAV>1 compounds + THE THRESHOLD LIST** | **the jackpot, online** |
| **Table S6** | 32 taste compounds × 6 samples | major |
| **Table S7** | **TAVs of 32 taste compounds + taste thresholds** | major |
| **Fig. S1** | G′/G″ frequency sweeps | moderate |
| **Fig. S2** | colour PCA vs 7 grilled fish | minor |
| **Fig. S3** | macro images of seafood references | minor |
| **Figs. S4–S11** | **all LibDock Scores** | moderate |
| **Fig. S12** | **all sensory scores** | moderate |
| **Fig. S13** | structure/texture/sensory correlation matrix | moderate |

**Nothing in the PDF itself was unreadable.** The text layer was complete and clean; all five
figures rendered legibly at 400–600 dpi. Fig. 4C (circular heatmap, z-scored) and Fig. 4D
(Sankey widths) are **inherently non-quantitative** rather than unreadable.

### Internal errors and inconsistencies found
1. **Table 1 row 1 reads "70 SPI + 20 SPI + 10 WG"** (SPI twice) — control is 90 SPI + 10 WG. p. 2
2. **Alcohol sentence is arithmetically impossible**: "greatest in RPI (1252.64), followed by
   EPPI (1871.78)". **1871.78 is EPPI's ALDEHYDE total, copy-pasted.** p. 9
3. **Ketone rank order skips ESPI**, which Fig. 4E ranks `c`, above EYP `d`. p. 9
4. **Barrel-zone arithmetic**: "seven temperature-controlled zones (T2-T7)" — T2–T7 is six. p. 3
5. **Eq. (2) "Water loss rate (%) = Wb/Wa × 100 %"** is the retention ratio, not the loss. p. 3
6. **Fig. 3C x-axis "Time (s)" 0–240** where the protocol is 240 **minutes**. p. 7
7. **Fig. 3C endpoint ≈92.2–93.3 % ⇒ 6.7–7.8 % loss, vs the text's "5 %".** p. 7 vs p. 8
8. **Abstract calls ENPI's β-sheet (27 %) "relatively high"** while EMPI (32 %) and EPPI (31 %)
   are higher; **the β-sheet ranking does not match the fibrous-degree ranking**, contradicting
   the paper's stated mechanism. p. 1 vs p. 7–8
9. **Abstract lists 2-heptanone (2778.1) as "substantially lower"** when it is 9.3 % below max. p. 1
10. **§2.11 "GC–MS NEST libraries"** — NIST, misspelled; **no library version given**. p. 3
11. **SPME fibre coating, incubation T/time and extraction time are never stated.** p. 3
12. **"air threshold value (μg/kg)"** — wrong unit, wrong basis, wrong operation (§4). p. 3
13. **No threshold source citation anywhere.** p. 3, p. 9
14. **Text's VIP list omits 6 of the 20 VIP>1 compounds**, including the #2 and #6 ranked ones
    (**2-acetyl-thiazole 2.32, 2-methyl-pyrazine 1.91**) — the only Maillard markers. p. 9 vs Fig. 4B
15. **Dimethyl phthalate (a plasticiser) has VIP 1.38** and is never discussed. Fig. 4B, p. 10
16. **E-nose sulfur claim** has no GC-MS support (only acetyl-thiazole in 82 compounds). p. 9
17. **PLS-DA has no R², no Q², no permutation test** (volatiles or taste). p. 9, p. 11
18. **PCA on texture: PC1 = 94.4 %** — effectively a hardness ranking dressed as multivariate. p. 5
19. **ΔE is vs a white plate, not vs a control** — not a browning index. Table 2, p. 5
20. **Table 2 L* significance letters are implausible** given the SDs (§7f). p. 5
21. **Docking argument is directionally confounded with the HS-SPME readout** (§10). p. 11
22. **Corresponding-author email typo**: "xuhuisos@dlpu.**deu**.cn". p. 2
23. **"octopu (Octopodidae)"** — typo, p. 2. **"creased to 250 °C"** for "increased", p. 3.
24. ⚠️ **REFERENCE-LIST CONTAMINATION — directly relevant to the repo's citation-integrity work:**
    - **Lv et al. (2011)**, a soymilk-blanching paper, is given **doi:10.3389/fnut.2022.973677**
      — which is **Liu et al. (2022)**'s DOI, listed correctly elsewhere on the same page.
      **A 2011 paper carrying a 2022 DOI.** p. 14
    - **Kobayashi et al. (1995)** appears **twice** in the reference list. p. 14
    - **Webb et al. (2023)** *Foods* **12**, 1586 and **Wen et al. (2022)** *Foods* **12**, 1586
      share the identical malformed DOI **10.3390/foods1208158** (13 digits, invalid). p. 14–15
    - **Yuliarti et al. (2025)**, *Int. J. Biol. Macromol.* 330, 147991, carries
      **doi:10.1016/B978-0-323-89842-3.00009-9** — a **book-chapter** DOI belonging to
      **Yuliarti et al. (2023)**. p. 15
    **⇒ At least 5 defective citations in a ~60-item list (~8 %), in a March-2026 Elsevier
    paper. Worth logging against the repo's own 30–45 % contamination finding as a base rate
    for the current literature.**

---

## 12. WHAT THIS UNLOCKS

1. **A same-lab, same-process, same-protocol REPRODUCIBILITY MEASUREMENT** on the repo's two
   most load-bearing markers: **hexanal ~10×, 2-pentylfuran ~17×, nonanal ~23×, with
   non-overlapping ranges** (§6). **This is the empirical justification for never ingesting
   HS-SPME μg/kg as absolute.**
2. **Four odour thresholds recovered exactly** (hexanal **5.000**, 1-octen-3-ol **1.500**,
   2-pentylfuran **5.800**, 2-heptanone **140.0** μg/kg), **three of them double-ended**, plus
   **one printed verbatim in the body text** (§3–4). **2-pentylfuran = 5.800 now confirmed from
   two independent papers and four independent numbers** — the group's threshold table is fixed
   and reusable, and **uncited**.
3. **A second independent instance of the "air threshold" laundering sentence**, in a different
   journal, with the paper itself printing an aqueous value one page later. **Two papers makes
   it a documented practice, not an anomaly.**
4. **A protein-source × volatile-class factorial at constant process and zero added
   carbohydrate**: total aldehydes span **12.4×** (194 → 2398 μg/kg) while **2-heptanone spans
   only 1.10×**. Protein source drives lipid-aldehyde chemistry hard and methyl-ketone
   chemistry barely at all. **The paper does not notice this dissociation.**
5. **An absolutely-quantified Maillard amine-substrate loading** (total FAAs
   **1613.74 → 3347.27 μg/g, 2.07×**), by validated HILIC-MS/MS with standards — the only
   trustworthy absolute number set in the paper.
6. **A complete, reusable HME process block** identical to the companion's (§1b), letting the
   two papers be treated as one experimental platform.
7. **The complete 24-compound OAV>1 list recovered from Fig. 4D**, which the text never prints
   — including three Maillard N/S heterocycles the paper silently drops.
8. **A recorded HS-SPME confound**: the docking result implies NPI *retains* off-flavours, which
   would lower headspace without lowering formation (§10). **"Less formed" vs "more bound" is
   not separable from these data**, and the paper's conclusion assumes the former.
9. **A citation-defect base rate** (~8 % of references) for a March-2026 Elsevier food-science
   paper (§11.24).

---

## USABILITY VERDICT

| tier | content | why |
|---|---|---|
| ✅ **BENCHMARK-GRADE** | **Process block** (§1b): 7 zone temps, 240 rpm, 2 kg/h, 65 % moisture, 70 °C cooling die, L/D 32:1 | fully specified, matches the companion exactly, two-paper corroboration |
| ✅ **BENCHMARK-GRADE** | **Formulations** (§1a): 70 SPI + 20 X + 10 WG, X ∈ {SPI, NPI, RPI, YP, MPI, PPI} | substitutive at constant total protein; clean design |
| ✅ **BENCHMARK-GRADE** | **Four recovered thresholds** (§4), as `used_by_source` records | three exact at both ends; one printed verbatim; 2-pentylfuran corroborated across two papers |
| ✅ **BENCHMARK-GRADE** | **Total FAAs 1613.74–3347.27 μg/g** (§8) | validated HILIC-MS/MS with authentic standards |
| ✅ **BENCHMARK-GRADE** | **Secondary structure percentages** (§7b) | printed on the bars; all six sum to 100 % |
| ✅ **BENCHMARK-GRADE** | **Colour a\* values** (Table 2, §7f) | instrumental, with SDs and n = 8 |
| ⚠️ **DIRECTIONAL ONLY** | **All volatile μg/kg** (§5) | semi-quantitative, single IS, no response factors, **unspecified fibre**, refuted as absolute by §6 |
| ⚠️ **DIRECTIONAL ONLY** | **Texture values** (§7a) | figure-digitised ±3 %; only fibrous degree is printed |
| ⚠️ **DIRECTIONAL ONLY** | **WAC, LF-NMR, water loss** (§7c–e) | partial printing; Fig. 3C contradicts its own text |
| ⚠️ **DIRECTIONAL ONLY** | **Sensory ranking** (§9) | zero numbers in the PDF; n = 10, unflagged |
| ⚠️ **DIRECTIONAL ONLY** | **E-nose** (§5h) | sensor semantics unsupported by GC-MS |
| ❌ **DO NOT INGEST** | **Every OAV in this paper** | matrix concentration ÷ aqueous threshold, no partition step, mislabelled as air |
| ❌ **DO NOT INGEST** | **ΔE as a browning index** (§7f) | measured against a white plate, dominated by L* |
| ❌ **DO NOT INGEST** | **"alcohols EPPI 1871.78 μg/kg"** (§5c) | arithmetically impossible; it is the aldehyde value |
| ❌ **DO NOT INGEST** | **"3 % / 5 %" water loss** (§7e) | contradicted by the paper's own figure (6.7–7.8 %) |
| ❌ **DO NOT INGEST** | **"β-sheet drives fibrous degree"** (§7b) | the paper's own rankings contradict it; correlations are in an online figure with no coefficients |
| ❌ **DO NOT INGEST** | **PCA/PLS-DA as evidence of separation** | PC1 = 94.4 % (texture) is a hardness proxy; no R²/Q²/permutation anywhere |
| ❌ **DO NOT INGEST** | **Anything from Fig. 2 (SEM) or Fig. 4C/4D** | qualitative / z-scored / relative widths |

## INGESTION RECOMMENDATION

**(a) INGEST the process + formulation block** as a **second, corroborating** record of the
same HME platform already recorded from the companion, keyed so the two papers share it:
`zones_C: [60,80,110,135,140,140,140]`, `cooling_die_C: 70`, `LD: 32`, `screw_mm: 704×22`,
`rpm: 240`, `dry_feed_kg_h: 2`, `moisture_pct: 65`, `water_kg_h: 3.7 (INFERRED from companion,
not printed here)`, `n: 3 (inferred from Fig. 4A)`, `design: substitutive, 20 % of SPI replaced`,
`added_carbohydrate: none`.

**(b) INGEST the four recovered thresholds as `used_by_source` records, NOT as repo thresholds**
— `hexanal 5.000`, `1-octen-3-ol 1.500`, `2-pentylfuran 5.800`, `2-heptanone 140.0` μg/kg,
tagged `basis_claimed: air`, `basis_actual: water`, `source_citation: NONE GIVEN`,
`recovery: back-solved, both ends exact (2-heptanone: one end)`. **Cross-link the 5.800 to the
companion's independent recovery of the same value.**

**(c) INGEST §6 (the cross-paper reproducibility test) AT TOP PRIORITY as a METHOD-LIMIT
RECORD**, not as a chemistry record: `same_lab_same_process_discrepancy: {hexanal: 10.0x,
2-pentylfuran: 16.9x, nonanal: 23.3x, 2-heptanone: 1.03x}`,
`ranges_overlap: {2-pentylfuran: false, nonanal: false}`,
`most_likely_cause: unspecified SPME fibre vs DVB/CAR/PDMS`. **This is the paper's single most
valuable contribution to the repo and it is a negative result about measurement, not a
finding about food.**

**(d) INGEST the protein-source volatile ordering as DIRECTIONAL ONLY**, as within-paper
ratios, never as ppb: `aldehyde_total: ERPI 2398 > EPPI 1872 ≈ EYP 1868 > EMPI 1746 > ESPI
≈1140 [fig] ≫ ENPI 194`; `hexanal: EPPI 1264.62 ≫ ENPI 74.86 (16.9x)`;
`2-heptanone: flat, 1.10x span`. Tag `quantification: semiquantitative_relative`,
`matrix: HME_65pct_moisture_plant_protein_no_added_sugar`, `fibre: UNSPECIFIED`.

**(e) INGEST total FAAs (1613.74–3347.27 μg/g, 2.07× span) as a genuine amine-substrate
loading benchmark**, flagged `no_paired_reducing_sugar_measurement`.

**(f) DO NOT INGEST any μg/kg volatile number as an absolute concentration, and DO NOT INGEST
any OAV.** Both for the companion's reasons and, now, for §6's stronger internal reason.

**(g) RECORD as named hazards:** the repeated **"odor threshold in air (μg/kg)"** sentence
(now two papers, two journals, one group); the **HS-SPME "less formed vs more bound" confound**
(§10); the **1871.78 copy-paste** (§5c); the **Fig. 3C unit and magnitude errors** (§7e); the
**abstract's β-sheet overstatement** (§7b); and the **five defective references** (§11.24).

## THE ONE RETRIEVAL WORTH REQUESTING NEXT

**The Supporting Information of `10.1016/j.foodres.2026.119010` — specifically Table S5.**
It contains **the OAVs of the 24 OAV>1 compounds AND the threshold list**, i.e. the remaining
**20 thresholds** this dossier could not recover, and — critically — **whether that table
carries source citations**. Retrieving it would settle the provenance question for the whole
threshold set that this group reuses across at least two papers. **Table S4 (82 × 6 = 492
quantified cells) and Table S1 (proximate composition of all 7 proteins, without which none of
the volatile data can be normalised to fat or carbohydrate) are the joint runners-up, and S1 is
arguably the higher scientific priority even though S5 is the higher mission priority.**
