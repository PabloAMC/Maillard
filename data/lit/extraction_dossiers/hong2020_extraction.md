# Hong & Kim 2020 — COMPLETE TRANSCRIPTION
### The corpus's FIRST true same-panel, same-method, paired matrix-vs-water orthonasal threshold set.

**Full extraction of every number in `data/articles/hong2020.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

Reads against, and materially revises, `k2_matrix_and_thresholds.md` §A, §D.1, §D.2(i).

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **Hyun Hee Hong, Mina K. Kim** |
| Title | **"Comparison of orthonasal thresholds of key volatile flavor compounds responsible for traditional doenjang flavor in two matrices: Water-based and soybean-based model system"** |
| Venue | *Journal of Sensory Studies* **2020;e12567** (article number; no volume/page) |
| Type | **SHORT CONTRIBUTION** — 5 pages |
| DOI | **10.1111/joss.12567** |
| Dates | Received 20 May 2019; Revised 21 Nov 2019; Accepted 26 Jan 2020 |
| Affiliation | Dept. of Food Science and Human Nutrition **and Fermented Food Research Center, Jeonbuk National University**, Jeonju-si, Jeollabukdo, South Korea |
| ⚠️ affiliation inconsistency | Methods p. 2 says testing was done "**in the sensory laboratory of the Department of Food Science and Human Nutrition at Chonbuk National University**" — Chonbuk was renamed Jeonbuk in 2016. Same institution, two names in one paper. Cosmetic. |
| Funding | NRF Korea, Ministry of Science ICT & Future Planning, **NRF-2017R1D1A3B03032520** |
| PDF character | **Clean digital text layer, fully machine-readable. Not one cell was unreadable.** Table 1 and Table 2 re-read at 400 dpi to confirm the ± column. |

**Provenance codes (same convention as `k1_kinetic_parameters.md`):**
- **[M]** measured / reported directly by the authors in a table
- **[C]** cited by them from another paper
- **[Z]** derived by this wave from the paper's own numbers — **never printed by the authors**
- **[D]** digitised — **not used in this dossier; there are no data figures in this paper**

---

## 1. THE ONE-PARAGRAPH ANSWER

**Every one of the 20 BET values (10 compounds × 2 matrices) is printed, legible and
complete, together with its dispersion statistic, a Z statistic, and a log P.** This is
the first dataset in the corpus where a matrix threshold and its water reference were
produced by **the same panel, the same day, the same 7-step 3-AFC ASTM E679-91 protocol
and the same BET rule** — so the ratios are `same_study_same_method`, unlike **every**
ratio in `k2_matrix_and_thresholds.md` §A.8. **Nine of the ten compounds are harder to
smell in the soybean matrix (29× to 2 035×, geometric mean 277×); one — ethyl-4-methyl-
pentanoate — is 6.5× EASIER, an inversion no other paper in the corpus reports.** Two
things the data say and the paper does not: **(i) the 1-σ band on the log-ratio collapses
from the 27–41× of the cross-study K2 estimates to 4.7×**, i.e. most of K2's spread was
method noise, not chemistry; and **(ii) hydrophobicity does not predict the shift at all
(r = −0.05 between log P and the log matrix shift)** — which rules out both reversible
hydrophobic binding and air/matrix partition as the dominant mechanism, and leaves the
paper's own proposal, **background-odour masking**, standing. Two defects: the ± column's
scale is **never stated** and cannot be linear (§4), and the printed Z statistic is
**arithmetically not the formula the paper writes down** (§5).

---

## 2. SYSTEM COMPOSITION AND METHOD — applies to every number below

Stated on p. 2 ("Materials and Method", §§2.1–2.4).

| variable | value as printed |
|---|---|
| Compounds (10) | **2,3-dimethyl pyrazine, ethyl-4-methylpentanoate, 2-pentyl furan, 4-vinyl phenol, hexanal, 2-methyl butanal, 3-methyl butanal, butyric acid, 4-ethyl phenol, dimethyl disulfide** |
| Compound selection basis | **[C]** Jo, Cho, Song, Shin & Kim (2011) — previously reported key doenjang volatiles |
| Source | **all from Sigma-Aldrich (St. Louis, MO)**; purities **not stated** |
| Stock solvent | **70% ethanol**, homogenised on a vortex mixer (VM-10, DAIHAN Scientific) |
| **Water matrix** | **distilled water**, **12 mL** per cup |
| **Soybean matrix** | **12 g** per cup |
| Soybean model preparation | **2.5 kg soybeans** (local grocery, Jeonju) **soaked in 3 L distilled water, 20 h, room temperature**; **autoclaved 40 min at 121 °C, 0.1 MPa** (LAC-5060S, Lab Tech); cooled to room temperature; **ground in a food processor (NINJA BL682KR) 1 min at 24 000 rpm** |
| Soybean model basis | "standard method of manufacturing soybean paste" **[C]** Park, Hwang, Jung & Lee (2002) |
| Dosing, water | aliquot of stock spiked into water; **"test samples at each row were prepared as soon as the aliquots of stock solution were spiked"** — i.e. **NO equilibration hold** |
| Dosing, soybean | spiked one compound at a time, bagged (polyethylene, CSE Co.), **stomached 30 s at 8 strokes/second** (BagMixer 400 cc, Interscience), **vacuum-sealed** (Freshield Elite IV), **stored 4 °C for 24 h** |
| ⚠️ **asymmetric equilibration** | **water = 0 h at ~22 °C; soybean = 24 h at 4 °C.** This is a protocol asymmetry, not a matrix property. See §6.1. |
| Cups | **70 mL (70 mm × 30 mm) white opaque plastic**, three-digit random code, lidded |
| Blank | **same amounts of water (12 mL) and soybean (12 g)** — i.e. **matrix-matched blanks** |
| Method | **7-series of 3-Ascending Forced Choice (7-3AFC)**, **ASTM E679-91** (ASTM 1997) |
| Panel | **n = 20** (stated in the Table 2 note); recruited by campus flyers + word of mouth among friends, staff and faculty of the department |
| Panel training | instructed on **proper sniffing technique** before testing; **no screening or performance criterion reported** |
| Palate/nose cleansing | **rested 1 min between each set of three**; **instructed to sniff their sleeves** between cups **[C]** Leksrisompong et al. (2010) |
| Replication | **each compound tested individually in each matrix. Number of replicate sessions per panellist is NOT stated** (contrast Leksrisompong 2010, which states duplicate on different days). |
| Individual BET rule | **geometric mean of the last concentration with an incorrect answer and the first with a correct answer; if the subject answered "not sure" on a correct choice, that concentration was multiplied by 1.41** before the geometric mean **[C]** Lawless, Hartono & Hernandez (2000) |
| Group BET rule | **geometric mean of the individual BETs** |
| Statistic | **Z = (BET2 − BET1)/(√[SD1² + SD2²])** with **BET1 = water, BET2 = soybean**; **"This study used a pairwise comparison using standard deviations, unlike the delta method that used the SE values."** |
| Serving temperature | **not stated** (room temperature implied) |
| Soybean matrix pH, protein, fat, moisture | **NOT STATED — none of the four.** See §6.2. |

---

## 3. TABLE 1 — the dosing series. **Anchor: Table 1, p. 2.** All cells legible. [M]

Title as printed: *"The concentration of volatile compounds of stock solution"*.
Column headers as printed: `Concentration (mg/kg) of stock solution prepared in 70% ethanol`
and `Concentration ranges covered for 7-AFC threshold tests (μg/kg)`, once per matrix.

| compound | water: stock **mg/kg** | water: series range **µg/kg** | soybean: stock **mg/kg** | soybean: series range **µg/kg** | span, water [Z] | span, soybean [Z] |
|---|---:|---|---:|---|---:|---:|
| 2,3-dimethyl pyrazine | **10,009** | **100–72,900** | **40,036** | **1,800 – 145,800** | 729× | ⚠️ **81×** |
| Ethyl-4-methylpentanoate | **16,839** | **3–2,187** | **574** | **0.07–54.0** | 729× | 771× |
| 2-pentyl furan | **17,188** | **0.67–486** | **17,188** | **1,300 – 947,000** | 725× | 728× |
| 4-vinyl phenol | **2,078** | **9.0–6,561** | **34,633** | **200–145,800** | 729× | 729× |
| Hexanal | **409** | **2.0–1,458** | **13,077** | **20–14,580** | 729× | 729× |
| 3-methyl butanal | **3,819** | **0.1–72.9** | **3,055** | **4.0–2,916** | 729× | 729× |
| Butyric acid | **3,871** | ⚠️ **0.075–6.075** | **9,544** | **20–14,580** | ⚠️ **81×** | 729× |
| 2-methyl butanal | **42,750** | **0.07–51.0** | **21,375** | **20–14,580** | 729× | 729× |
| 4-ethyl phenol | **1,039** | **9.0–6,561** | **24,750** | **200–145,800** | 729× | 729× |
| Dimethyl disulfide | **5,035** | **0.1–72.9** | **10,070** | ⚠️ **400–29,160** | 729× | ⚠️ **72.9×** |

**⚠️ THREE ROWS ARE INCONSISTENT WITH A 7-STEP SERIES.** A 7-concentration ascending
series at the implied dilution factor of 3 spans **3⁶ = 729×**, and **17 of the 20 ranges
do exactly that** (725–771× after rounding). Three do not:
- **butyric acid, water: 0.075–6.075 µg/kg = 81× = 3⁴, a FIVE-step series.**
- **2,3-dimethylpyrazine, soybean: 1,800–145,800 = 81× = 3⁴, a FIVE-step series.**
- **dimethyl disulfide, soybean: 400–29,160 = 72.9×.** Note **360 × 3⁴ = 29,160** exactly;
  the printed lower bound of 400 is the only value in the whole table that is not an exact
  member of its own geometric series. **Most parsimonious reading: a typo for 360.** Not
  asserted.

Consequence: for those three rows the effective series is shorter or mis-stated, and the
BET for a panellist who never detected (or always detected) is bounded by a narrower
window. **All 20 group BETs nonetheless fall inside their own printed series ranges [Z]**
— checked cell by cell — so no BET is an extrapolation beyond the series.

---

## 4. TABLE 2 — THE PAIRED THRESHOLD SET. **Anchor: Table 2, p. 3.** All cells legible. [M]

Title as printed: *"The BET values in water and soybean models"*.
Note as printed: *"Numbers represent Group BET ± SD of individual threshold value for each
compound (n = 20). Z-value = (BET2 − BET1)/(√[SD1² + SD2²]), while, BET1 represent best
estimate threshold value in water-based model system and BET2 represents best estimate
threshold value in soybean-based model system. Aroma descriptions were collected from
previously reported studies."*

**Units as printed: µg/kg in both matrices.** (Note: the soybean matrix is dosed **by mass**
— 12 g — so µg/kg is dimensionally correct there; the water matrix is dosed by **volume**
— 12 mL — and µg/kg is being used for what is really µg/L. At water density the two
coincide to <0.5 %; not a defect, but record the basis.)

| compound | **water BET µg/kg** [M] | ± | **soybean BET µg/kg** [M] | ± | **soybean/water ratio** [Z] | **Δlog10** [Z] | Z printed [M] | log P [C] | aroma description [C] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 2,3-dimethyl pyrazine | **401.0** | 0.7 | **19,600.8** | 0.6 | **48.9×** | +1.689 | 19,173 | 0.54 | Nutty-roast ᵃ |
| **Ethyl-4-methyl pentanoate** | **9.7** | 0.9 | **1.5** | 0.5 | **0.155× (⇒ 6.47× EASIER)** | **−0.811** | −8.1 | 2.10 | Citrus, pineapple ᵇ |
| 2-pentyl furan | **8.6** | 0.7 | **11,690.2** | 0.6 | **1,359×** | +3.133 | 11,681 | 3.70 | Licorice, beany ᶜ |
| 4-vinyl phenol | **166.0** | 0.6 | **25,139.9** | 0.6 | **151×** | +2.180 | 24,942 | 2.40 | Phenolic, medicinal, spicy ᵈ |
| **Hexanal** | **12.6** | 0.8 | **1,669.6** | 0.7 | **132.5×** | +2.122 | 1,657 | 1.77 | Green/sweet ᵉ |
| 3-methyl butanal | **0.5** | 0.6 | **131.6** | 0.6 | **263×** | +2.420 | 131 | 1.79 | Dark chocolate/malty ᵉ |
| Butyric acid | **0.9** | 0.9 | **1,831.2** | 0.8 | **2,035×** | +3.308 | 1,830 | 1.00 | Cheesey/rancid ᵉ |
| 2-methyl butanal | **2.1** | 0.6 | **548.9** | 0.7 | **261×** | +2.417 | 546 | 1.76 | Dark chocolate/malty ᵉ |
| 4-ethyl phenol | **23.2** | 0.4 | **675.1** | 0.6 | **29.1×** | +1.464 | 643 | 2.58 | Leather, phenol, spice, stable, smoke, creosote ᵈ |
| Dimethyl disulfide | **2.7** | 0.7 | **4,737.9** | 0.4 | **1,755×** | +3.244 | 4,734 | 1.22 | Brothy ᶠ |

Footnote sources as printed: ᵃ Akiyama et al. (2005); ᵇ Kishimoto, Wanikama, Kono & Shibata
(2006); ᶜ Lam & Proctor (2003); ᵈ Lentz (2018); ᵉ Carunchia-Whetstine, Cadwallader & Drake
(2005); ᶠ Whitson, Miracle & Drake (2010).

### 4.1 ⚠️ THE ± COLUMN'S SCALE IS NEVER STATED AND CANNOT BE LINEAR

The note calls the ± column "**SD of individual threshold value**". Read as a **linear** SD
in µg/kg this is impossible:

| compound / matrix | BET | printed ± | **implied RSD if linear** |
|---|---:|---:|---:|
| 4-vinyl phenol, soybean | 25,139.9 | 0.6 | **0.0024 %** |
| 2,3-dimethylpyrazine, soybean | 19,600.8 | 0.6 | **0.0031 %** |
| dimethyl disulfide, soybean | 4,737.9 | 0.4 | **0.0084 %** |
| 3-methyl butanal, water | 0.5 | 0.6 | **120 %** |

A 20-panellist sensory threshold cannot have a 0.0024 % relative SD, and the four-decade
swing in implied RSD across rows is itself disqualifying. **The only self-consistent
reading is that the ± values are SDs of the log-transformed individual BETs** — 0.4 to 0.9
**decades**, which is squarely the normal range for a 20-person orthonasal panel (Vega 1994
and Brewer 1995 in `k2_matrix_and_thresholds.md` §A.1–A.2 show 15×–37× individual ranges,
i.e. 1.2–1.6 decades total spread). **The base of the logarithm is not stated.** This
dossier assumes **log₁₀** — that is the base implied by the "log transformed then geometric
mean" convention of Leksrisompong et al. (2010), which this paper cites for its statistics.
**Record as `dispersion_scale: log10_assumed_not_stated`. If a downstream wave needs the
uncertainty, this is a hard blocker until the base is confirmed.**

Under the log₁₀ reading, the 20 dispersions are:

| ± value | count | as a fold-range (10^±SD) |
|---:|---:|---|
| 0.4 | 1 | 2.5× |
| 0.5 | 1 | 3.2× |
| 0.6 | 7 | 4.0× |
| 0.7 | 6 | 5.0× |
| 0.8 | 2 | 6.3× |
| 0.9 | 3 | 7.9× |

**Mean 0.65 decades ⇒ a typical individual-panellist spread of ±4.5× about the group BET.**
That is a directly usable per-compound sensory-threshold uncertainty prior — **the first one
in the corpus that is per-compound and internally consistent** (contrast Tian 2020's
0.16–0.47 % RSDs, named laundering hazard #13 in `k2_matrix_and_thresholds.md`).

### 4.2 ⚠️ THE PRINTED Z-VALUES ARE NOT THE FORMULA THE PAPER WRITES DOWN

The note defines `Z = (BET2 − BET1)/(√[SD1² + SD2²])`. Evaluate it:

| compound | BET2 − BET1 | √(SD1²+SD2²) | formula ⇒ Z | **Z printed** | printed ÷ (BET2−BET1) |
|---|---:|---:|---:|---:|---:|
| hexanal | 1,657.0 | 1.063 | 1,558.8 | **1,657** | **1.0000** |
| 2-pentyl furan | 11,681.6 | 0.922 | 12,671 | **11,681** | **0.99999** |
| butyric acid | 1,830.3 | 1.204 | 1,520.2 | **1,830** | **0.99984** |
| 3-methyl butanal | 131.1 | 0.849 | 154.4 | **131** | **0.99924** |
| dimethyl disulfide | 4,735.2 | 0.806 | 5,874 | **4,734** | **0.99975** |
| 2-methyl butanal | 546.8 | 0.922 | 592.9 | **546** | **0.99854** |
| ethyl-4-methylpentanoate | −8.2 | 1.030 | −7.96 | **−8.1** | **0.98780** |
| 4-ethyl phenol | 651.9 | 0.721 | 904.2 | **643** | **0.98635** |
| 2,3-dimethylpyrazine | 19,199.8 | 0.922 | 20,824 | **19,173** | **0.99860** |
| 4-vinyl phenol | 24,973.9 | 0.849 | 29,417 | **24,942** | **0.99872** |

**Every printed Z equals the raw difference `BET2 − BET1` to within 1.4 %.** The denominator
was effectively 1.00, not √(SD1²+SD2²). **The printed Z column is therefore a difference in
µg/kg, not a test statistic** — which is also why its values reach 24,942, a number no
standard normal deviate can take. **Do not ingest the Z column, and do not use it to decide
significance.**

**The correct statistics, computed by this wave [Z].** Because the BETs are log-normal and
the SDs are (§4.1) log-scale, the test must be on Δlog₁₀:

| compound | Δlog₁₀ | **Z using SDs** (as the paper says it did) | **Z using SEs = SD/√20** (the delta-method convention it says it departed from) |
|---|---:|---:|---:|
| 2,3-dimethyl pyrazine | +1.689 | +1.83 | **+8.19** |
| **ethyl-4-methylpentanoate** | −0.811 | **−0.79 (NOT significant)** | **−3.52 (significant)** |
| 2-pentyl furan | +3.133 | +3.40 | +15.20 |
| 4-vinyl phenol | +2.180 | +2.57 | +11.49 |
| hexanal | +2.122 | +2.00 | +8.93 |
| 3-methyl butanal | +2.420 | +2.85 | +12.76 |
| butyric acid | +3.308 | +2.75 | +12.29 |
| 2-methyl butanal | +2.417 | +2.62 | +11.73 |
| 4-ethyl phenol | +1.464 | +2.03 | +9.08 |
| dimethyl disulfide | +3.244 | +4.02 | +18.00 |

**⚠️ Load-bearing consequence.** The paper's headline exception — *"ethyl-4-methylpentanoate
had higher BET in water-based model (9.7 µg/kg), while BET in the soybean-based model was
1.5 µg/kg (p < .05)"* (p. 4) — **is significant on the SE convention (Z = −3.52) and NOT
significant on the SD convention the paper explicitly says it used (Z = −0.79)**. Under the
SD convention, **2,3-dimethylpyrazine (Z = +1.83) also drops below |1.96|**. The paper
prints **no p-value anywhere** and the "(p < .05)" claim is unsupported by any computable
statistic in the paper.
**⇒ The inversion is real in the point estimates but its significance is convention-
dependent. Flag it; do not treat "ethyl-4-methylpentanoate is easier in soy" as established.**

---

## 5. THE COMPARISONS THE PAPER MAKES TO PRIOR LITERATURE — all [C], all in water

Quoted from pp. 3–4. These are the paper's own reality checks on its water arm.

| compound | **this study, water [M]** | prior value [C] | source as cited | **ratio [Z]** |
|---|---:|---:|---|---:|
| 2,3-dimethyl pyrazine | **401.0 µg/kg** | 2,500 µg/kg | Guadagni, Buttery & Turnbaugh (1972) | **6.2× lower here** |
| butyric acid | **0.9 µg/kg** | 0.23 µg/kg | Kim, Drake & Drake (2011) | 3.9× higher here |
| butyric acid | **0.9 µg/kg** | 1,274 µg/kg | Karagul-Yuceer et al. (2003) | **1,416× lower here** |
| 2-methyl butanal | **2.1 µg/kg** | 0.21 µg/kg | Kim et al. (2011) | 10× higher here |
| 3-methyl butanal | **0.5 µg/kg** | 160 µg/kg | Kim et al. (2011) | 320× lower here |
| 3-methyl butanal | **0.5 µg/kg** | 0.2 µg/kg | Guadagni et al. (1972) | 2.5× higher here |

**⚠️ The paper calls all of these "within the range of variability of threshold values."
They span 1,416×.** Two independent water values for **the same compound, butyric acid,
differ by 5,540×** (0.23 vs 1,274) and two for **3-methyl butanal by 800×** (0.2 vs 160).
This is a **within-water, cross-study** spread and it is *larger* than most of this paper's
own matrix shifts. **That is the single most important calibration fact in this dossier and
it belongs on every threshold record: the cross-study reproducibility of an aqueous
orthonasal threshold is worse than three orders of magnitude, so any matrix/water ratio
that crosses studies is dominated by that noise.** It is the quantitative justification for
`k2_matrix_and_thresholds.md` §D.2(i) and it is why Hong's *within*-study ratios matter.

**⚠️ Citation defect.** The reference list gives *"Guadagni, D. G., Buttery, R. G., &
Turnbaugh, J. G. (1972). Odour thresholds and similarity ratings of some potato Chip
components. **Journal of Agricultural and Food Chemistry, 23**, 1435–1444."* The correct
venue is ***Journal of the Science of Food and Agriculture* 23, 1435–1444**. Wrong journal,
right volume/pages. The same paper is the primary source of four of the six aqueous anchors
in `k2_matrix_and_thresholds.md` §A.8 (there reached second-hand through Vega 1994) — so
this is now the **second** corpus paper that mis-cites it.
**⚠️ Year mismatch.** The text cites "**Karagul-Yuceer et al. (2003)**" twice (pp. 3, 4);
the reference list has only **Karagul-Yuceer, Drake & Cadwallader (2004)**, *J. Sens. Stud.*
19, 1–13. One of the two is wrong; not resolvable from the PDF.

---

## 6. WHAT THE PAIRED SET ACTUALLY SHOWS — the analysis the paper does not do [Z]

### 6.1 The distribution of the shift, and how it compares to K2's cross-study estimates

Over the **nine elevated compounds** (excluding the inverted ester):

| statistic | **Hong 2020, soybean/water, same panel** | K2 §D.1, beef/water | K2 §D.1, gelatin/water |
|---|---:|---:|---:|
| geometric mean | **277×** | 906× | 33× |
| minimum | **29.1×** (4-ethylphenol) | 77× | 3.4× |
| maximum | **2,035×** (butyric acid) | 6,714× | 914× |
| **max / min** | **70×** | 87× | 269× |
| **log₁₀ SD** | **0.668 decades** | 0.71 decades | 0.80 decades |
| **1-σ band width** | **4.7×** | **27×** | **41×** |

**⚠️ Note the two different things this table says.** The **max/min spread (70×) is as bad
as K2's** — so `k2_matrix_and_thresholds.md` §D.4's central conclusion, *no general matrix
correction factor*, **survives intact and is now confirmed on same-method data**. But the
**1-σ band collapses from 27–41× to 4.7×**, because a log-SD of 0.668 decades over nine
compounds is far tighter than the same statistic computed across labs. **The distribution
is tighter than K2 could see, but the tails are just as long, and it is the tails that
break a single factor.** Including the inverted ester, the ten-compound geometric mean is
**131×** with a log₁₀ SD of **1.206 decades (16× band)** — the inversion alone nearly
triples the dispersion.

### 6.2 ⚠️ HYDROPHOBICITY DOES NOT PREDICT THE SHIFT. r = −0.05.

Regressing the printed log P against the log₁₀ matrix shift over all ten compounds:

**r = −0.052; slope = −0.070 decades per log P unit.**

There is **no relationship whatsoever**, and the sign of the (null) slope is the *opposite*
of what binding or partition predicts. Concretely:

| compound | log P | matrix shift |
|---|---:|---:|
| 2-pentyl furan | **3.70** (most hydrophobic) | 1,359× |
| **ethyl-4-methylpentanoate** | **2.10** | **0.155× — inverted** |
| **butyric acid** | **1.00** | **2,035× — the largest shift, from a nearly hydrophilic acid** |
| 2,3-dimethyl pyrazine | **0.54** (least hydrophobic) | 48.9× |
| 4-ethyl phenol | 2.58 | **29.1× — the smallest elevated shift, from a hydrophobic phenol** |

**This is the strongest single-study refutation of the partition/binding account of matrix
threshold shifts anywhere in the corpus.** `k2_matrix_and_thresholds.md` §B.4 capped
reversible hydrophobic binding at 1.9–7.6× at 100 g/L and §B.5 showed perception is *less*
sensitive to binding than headspace is. Hong now adds: **even the rank order is wrong.**
A mechanism that scales with hydrophobicity cannot produce a set whose extremes are a
log P 1.00 acid (2,035×) and a log P 2.58 phenol (29.1×).

**What the paper itself proposes, and what the data are consistent with:** p. 4 —
*"Typically, odor thresholds increase in complex food matrices due to background aromas
present in the food model system."* **Background masking is compound-specific (it depends
on the overlap between the target's odour quality and the matrix's own volatile profile),
which is exactly the pattern observed and exactly what no partition model can reproduce.**
Autoclaved ground soybean is itself rich in **hexanal and 2-pentylfuran** (the two canonical
beany markers) — and those two compounds show shifts of 132× and 1,359×, at the top of the
distribution, while **ethyl-4-methylpentanoate, a fruity ester with no soybean counterpart,
is the one compound that gets EASIER.** The masking hypothesis predicts the inversion; the
binding hypothesis cannot. **Not asserted as proven — the paper measured no background
concentrations — but it is the reading the data support and it should be recorded as a
named candidate mechanism.**

### 6.3 The hexanal number, in context with everything the corpus has

| matrix | protein / composition | hexanal matrix threshold | aqueous reference | **ratio** | pairing quality |
|---|---|---:|---:|---:|---|
| 3 % gelatin, 22 °C (Vega 1994) | 30 g/L collagen hydrolysate | 58 ppb | 4.5 ppb (Guadagni 1972) | **12.9×** | cross-study, cross-method |
| **soybean paste, ~22 °C (Hong 2020)** | **not stated** | **1,669.6 µg/kg** | **12.6 µg/kg, same panel** | **132.5×** | **SAME STUDY, SAME METHOD** ✅ |
| cooked ground beef 1:1, 45 °C (Brewer 1995) | ~100 g/L | 5.87 ppm | 4.5 ppb (Guadagni 1972) | **1,304×** | cross-study + `dose_added_pre_cook` |

**Hong's 132.5× is the only defensible hexanal matrix factor in the corpus**, and it lands
an order of magnitude above gelatin and an order below Brewer's added-dose figure — exactly
where `k2_matrix_and_thresholds.md` §D.2(ii) predicted a properly-controlled number would
fall once Brewer's pre-cook thermal loss was removed. **This is an independent
corroboration of the Brewer reclassification decision (K2 owner item #2).**

### 6.4 What the soybean matrix does NOT come with

**No composition is reported.** Not pH, not protein, not fat, not moisture, not solids. A
soaked-and-autoclaved whole-soybean grind is roughly **36–40 % protein and 18–20 % fat on
dry basis**, and after a 20 h soak and a wet grind the wet paste is perhaps **60–65 %
moisture**, i.e. **on the order of 130–150 g protein/kg and 65–75 g fat/kg** — but **none of
that is in the paper and none of it is asserted here**. Without it, the 132.5× hexanal
factor **cannot be scaled to another protein loading**, and cannot be decomposed into
protein and lipid terms. `matrix_composition: not_stated` is a mandatory field on every one
of these ten records.

---

## 7. NAMED LAUNDERING HAZARDS — this paper's contributions to the "342/200 list"

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| H-1 | **The Z-value column** (up to 24,942) presented as a test statistic from `Z = (BET2 − BET1)/√(SD1²+SD2²)` | Every entry equals **BET2 − BET1** to within 1.4 %; the denominator is 1.00, not the stated root. It is a **difference in µg/kg**, not a Z. | Table 2 + its own note, p. 3 |
| H-2 | **"BET in the soybean-based model was 1.5 µg/kg (p < .05)"** for ethyl-4-methylpentanoate | **No p-value is computed anywhere in the paper.** On the SD convention the paper says it used, the correct Z is **−0.79 — not significant.** It reaches significance only on the SE convention the paper explicitly says it did *not* use. | p. 4 |
| H-3 | **"these values are within the range of individual variability"** covering water butyric acid at 0.9 vs 0.23 and 1,274 µg/kg | The two cited prior values differ from each other by **5,540×**. Calling a 1,416× disagreement "within the range of variability" launders a three-order-of-magnitude reproducibility failure into a footnote. | pp. 3–4 |
| H-4 | **"2,3-dimethyl pyrazine ... has two methyl groups and amino group chemical structure, and this methyl group increases the water-solubility. Thus, this methyl group affects a higher BET value in water-based model system."** | **Pyrazine has no amino group** — it is a 1,4-diazine. And **methyl substitution decreases, not increases, water solubility.** Both halves of the mechanistic sentence are chemically wrong. Do not propagate the reasoning; the number (401.0 µg/kg) is unaffected. | p. 3 |
| H-5 | **"the ± SD of individual threshold value"** with no stated scale | Cannot be linear (0.0024 % RSD on one row, 120 % on another). Must be a **log SD**, base not stated. Ingesting it as a linear µg/kg SD would understate threshold uncertainty by **four orders of magnitude** on the soybean rows. | Table 2 note, p. 3 |
| H-6 | **"Guadagni ... (1972). ... *Journal of Agricultural and Food Chemistry*, 23, 1435–1444"** | Correct venue is ***J. Sci. Food Agric.*** 23, 1435–1444. Second corpus paper to mis-cite the source of four of K2's six aqueous anchors. | Reference list, p. 4 |
| H-7 | Table 1 dosing ranges for **butyric acid/water (81×)**, **2,3-DMP/soybean (81×)**, **DMDS/soybean (72.9×)** | Inconsistent with the stated **7-step** series (should be 729×). Two are 5-step; the DMDS lower bound is the only value in the table that is not an exact member of its own geometric series (**360**, not 400, would fit). | Table 1, p. 2 |

**Plus three unresolved gaps, recorded but not resolvable from the PDF:**
**(a)** the base of the logarithm behind the ± column; **(b)** the number of replicate
sessions per panellist (never stated; contrast Leksrisompong 2010's "in duplicate on
different days"); **(c)** the entire composition of the soybean matrix.

---

## 8. CONSOLIDATED NEW-PARAMETER TABLE

**Class M = measured and printed; Z = derived by this wave; C = cited by the paper.**
Conditions common to all: **orthonasal, 3-AFC ASTM E679-91, 7-step ascending, factor-3
dilution, n = 20 untrained-but-instructed panellists, matrix-matched blank, group BET =
geometric mean of individual BETs, room temperature, Korean population, Jeonbuk/Chonbuk
National University, 2019.** Matrix A = **distilled water, 12 mL, dosed and served
immediately**. Matrix B = **soaked (20 h) + autoclaved (121 °C, 0.1 MPa, 40 min) + ground
(24 000 rpm, 1 min) whole soybean paste, 12 g, stomached 30 s, vacuum-sealed, 24 h at 4 °C**.

| # | parameter | value | units | matrix / conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1 | BET 2,3-dimethylpyrazine | **401.0** | µg/kg | A | M | T2 p.3 |
| 2 | BET 2,3-dimethylpyrazine | **19,600.8** | µg/kg | B | M | T2 p.3 |
| 3 | BET ethyl-4-methylpentanoate | **9.7** | µg/kg | A | M | T2 p.3 |
| 4 | BET ethyl-4-methylpentanoate | **1.5** | µg/kg | B | M | T2 p.3 |
| 5 | BET 2-pentylfuran | **8.6** | µg/kg | A | M | T2 p.3 |
| 6 | BET 2-pentylfuran | **11,690.2** | µg/kg | B | M | T2 p.3 |
| 7 | BET 4-vinylphenol | **166.0** | µg/kg | A | M | T2 p.3 |
| 8 | BET 4-vinylphenol | **25,139.9** | µg/kg | B | M | T2 p.3 |
| 9 | **BET hexanal** | **12.6** | µg/kg | A | M | T2 p.3 |
| 10 | **BET hexanal** | **1,669.6** | µg/kg | B | M | T2 p.3 |
| 11 | BET 3-methylbutanal | **0.5** | µg/kg | A | M | T2 p.3 |
| 12 | BET 3-methylbutanal | **131.6** | µg/kg | B | M | T2 p.3 |
| 13 | BET butyric acid | **0.9** | µg/kg | A | M | T2 p.3 |
| 14 | BET butyric acid | **1,831.2** | µg/kg | B | M | T2 p.3 |
| 15 | BET 2-methylbutanal | **2.1** | µg/kg | A | M | T2 p.3 |
| 16 | BET 2-methylbutanal | **548.9** | µg/kg | B | M | T2 p.3 |
| 17 | BET 4-ethylphenol | **23.2** | µg/kg | A | M | T2 p.3 |
| 18 | BET 4-ethylphenol | **675.1** | µg/kg | B | M | T2 p.3 |
| 19 | BET dimethyl disulfide | **2.7** | µg/kg | A | M | T2 p.3 |
| 20 | BET dimethyl disulfide | **4,737.9** | µg/kg | B | M | T2 p.3 |
| 21 | **matrix shift, hexanal** | **132.5** | ×, dimensionless | B/A, same panel | Z | T2 p.3 |
| 22 | matrix shift, 2-pentylfuran | **1,359** | × | B/A | Z | T2 p.3 |
| 23 | matrix shift, butyric acid | **2,035** | × | B/A — largest | Z | T2 p.3 |
| 24 | matrix shift, 4-ethylphenol | **29.1** | × | B/A — smallest elevated | Z | T2 p.3 |
| 25 | matrix shift, ethyl-4-methylpentanoate | **0.155** | × (INVERTED) | B/A | Z | T2 p.3 |
| 26 | **geometric-mean matrix shift, 9 elevated compounds** | **277** | × | B/A, same panel | Z | §6.1 |
| 27 | **log₁₀ SD of the matrix shift, 9 compounds** | **0.668** | decades (⇒ 4.7× 1-σ band) | " | Z | §6.1 |
| 28 | geometric-mean matrix shift, all 10 | **131** | × | " | Z | §6.1 |
| 29 | log₁₀ SD, all 10 | **1.206** | decades (⇒ 16.1× band) | " | Z | §6.1 |
| 30 | **corr(log P, log matrix shift)** | **−0.052** | r, n = 10 | " | Z | §6.2 |
| 31 | slope of log-shift on log P | **−0.070** | decades per log P unit | " | Z | §6.2 |
| 32 | **per-compound panel dispersion** | **0.4–0.9 (mean 0.65)** | decades (⚠️ base assumed log₁₀, NOT STATED) | individual BETs, n = 20 | M(scale Z) | T2 note |
| 33–42 | log P, 10 compounds | 0.54 / 2.10 / 3.70 / 2.40 / 1.77 / 1.79 / 1.00 / 1.76 / 2.58 / 1.22 | — | cited, source not given | C | T2 p.3 |
| 43 | prior aqueous BET, 2,3-dimethylpyrazine | **2,500** | µg/kg | water | C | Guadagni 1972 via p.3 |
| 44 | prior aqueous BET, butyric acid | **0.23** | µg/kg | water | C | Kim et al. 2011 via p.3 |
| 45 | prior aqueous BET, butyric acid | **1,274** | µg/kg | water | C | Karagul-Yuceer 2003/4 via p.4 |
| 46 | prior aqueous BET, 2-methylbutanal | **0.21** | µg/kg | water | C | Kim et al. 2011 via p.4 |
| 47 | prior aqueous BET, 3-methylbutanal | **160** | µg/kg | water | C | Kim et al. 2011 via p.4 |
| 48 | prior aqueous BET, 3-methylbutanal | **0.2** | µg/kg | water | C | Guadagni 1972 via p.4 |
| 49 | **cross-study aqueous reproducibility, butyric acid** | **5,540** | × spread between two cited water values | — | Z | §5 |
| 50 | cross-study aqueous reproducibility, 3-methylbutanal | **800** | × | — | Z | §5 |

**Mandatory provenance fields on all 20 threshold records** (extending
`k2_matrix_and_thresholds.md` §D.4.2): `criterion: 3AFC_ASTM_E679-91`,
`thermal_step_after_dosing: false`, `concentration_verified: false`,
`cross_study_cross_method: FALSE` ← **the first records in the corpus that can carry this**,
`equilibration_asymmetric: true (water 0 h @22 °C vs soybean 24 h @4 °C)`,
`matrix_composition: not_stated`, `dispersion_scale: log10_assumed_not_stated`,
`printed_Z_statistic: DO_NOT_INGEST`.

---

## 9. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **This is a proposal only.** Hong 2020 is a **new source**. Under the repo's
> `docs/reference/FIT_HOLDOUT_DECLARATION.md` regime, **a declaration amendment must be
> approved before any wave fits to any row below.** This dossier does not edit the
> declaration and no wave may treat this section as authorisation.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **Table 2, water column (10 BETs)** | 10 | **HOLD-OUT** | These are aqueous thresholds obtained by an independent panel and population. They are the natural check on whatever aqueous threshold list the repo already ships. Fitting them would destroy their only value. |
| **Table 2, soybean column (10 BETs)** | 10 | **HOLD-OUT — highest-value hold-out in the wave** | This is the corpus's only same-panel matrix threshold set. If the repo ever ships a matrix correction, this is the set that adjudicates it. **Do not fit.** |
| **The 10 paired ratios (§4, col. 6)** | 10 | **HOLD-OUT, and the recommended primary acceptance target for any matrix-correction work** | Suggested criterion: a matrix model is acceptable only if it reproduces **≥7/10 ratios within 5×** (the measured 1-σ band) **and reproduces the sign on all 10**, including the ethyl-ester inversion. Nothing in the repo today can do the last part. |
| **The log P ⊥ shift null result (§6.2)** | 1 | **HOLD-OUT — a falsification test, not a parameter** | Any shipped matrix term that is a monotone function of log P is refuted by r = −0.05 on this set. Use it as a guard, not a fit target. |
| **Per-compound panel dispersion (± column)** | 20 | **NOT FITTABLE, NOT INGESTIBLE, pending clarification** | The log base is not stated (§4.1). Blocked until confirmed. If confirmed as log₁₀, propose role = **PRIOR** on sensory-threshold uncertainty (mean 0.65 decades). |
| **Printed Z column** | 10 | **REJECT — do not ingest under any role** | Arithmetically not the stated statistic (§4.2, H-1). |
| **Table 1 dosing ranges** | 20 | **METADATA only** | Useful for bounding BETs against their series; three rows internally inconsistent (H-7). |
| **The six cited prior aqueous values (§5)** | 6 | **REJECT as parameters; RETAIN as a calibration fact** | Second-hand, methods unstated, mutually inconsistent by up to 5,540×. Their value is the §5 reproducibility argument, not the numbers. |

**Recommended companion action for the orchestrator (not taken here):** this paper's
existence makes the K2 owner-decision item #2 (**reclassify Brewer 1995 as
`dose_added_pre_cook`**) considerably stronger — see §6.3 — and it supplies the first
`cross_study_cross_method: false` rows the threshold lane has ever had.

---

## 10. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Kim, M. K., Drake, S. L. & Drake, M. A. (2011)**, *J. Sensory Studies* **26**(4),
   278–290 — orthonasal thresholds of volatile compounds in **three matrices (water, pH 5.5
   buffer, and oil)** by the same protocol lineage. **This would be a second same-panel
   multi-matrix paired set, and it includes an oil arm** — directly comparable to
   Leksrisompong 2010 (this wave) and to the Guadagni paraffin-oil values held second-hand
   in `k2_matrix_and_thresholds.md` §A.3. **Highest-yield retrieval in wave K4b.**
2. **Karagul-Yuceer, Y., Drake, M. A. & Cadwallader, K. R. (2004)**, *J. Sensory Studies*
   **19**, 1–13 — **BETs of 13 volatile compounds in skim milk**, the source of the 1,274
   µg/kg butyric acid value. Would settle the §5 reproducibility question and add a
   third matrix (skim milk) on the same protocol lineage.
3. **Jo, Y. J., Cho, I. H., Song, C. K., Shin, H. W. & Kim, Y. S. (2011)**, *J. Food Sci.*
   **76**(3), C368–C379 — the volatile profile of doenjang from which the ten compounds were
   chosen. **This is the only route to the soybean matrix's own background concentrations**,
   which §6.2 identifies as the load-bearing missing variable behind the masking hypothesis.
4. **Guadagni, Buttery & Turnbaugh (1972)**, *J. Sci. Food Agric.* **23**, 1435–1444 —
   already a K2 runner-up retrieval; this paper is the second to depend on it and the second
   to mis-cite it.
