# Ma, Cao, Li, Zhang, Guo & Li 2026 — EXTRACTION (D-xylose with L-cysteine or L-methionine dosed into high-moisture extruded soy protein: MFT and FFT by GC×GC-SCD peak area, water thresholds by 3-AFC, binary hexanal masking, sensory)

**Source on disk:** `data/articles/Ma2026.pdf` (the publisher's PDF, 16 pages, with a text layer). Read and
checked by eye on 2026-10-09: every number below was read from the page image and cross-checked against
`pdftotext -layout`. Fig. 3 (sensory, p. 7), Fig. 4 (MFT and FFT bars, p. 10) and Fig. 5 (threshold curves,
p. 12) have no data table: Fig. 4 was re-rendered at 200 dpi and its bar heights are "read from graph,
approx."; the threshold and R labels printed inside Fig. 5 were re-rendered at 500 dpi and read as text.
The supplement (Table S1 sensory criteria, Table S2 e-nose sensors, Table S3 GC-O, Fig. S1 single-compound
thresholds) is not on disk. Not to be confused with `ma2024_extraction.md` (a different paper). Written
for the plant-based meat-flavour use case; it is the most product-shaped paper on disk for the sulfur lane.

| field | value |
|---|---|
| Title | "Mechanistic insights into off-flavor masking of textured soy protein using a Maillard reaction system with D-xylose and sulfur-containing amino acids" |
| Authors | Jian Ma, Rui Cao, Xuejie Li, Wentao Zhang, Zengwang Guo, Jian Li (Beijing Technology and Business University; Tobacco Research Institute CAAS, Qingdao; Northeast Agricultural University, Harbin) |
| Venue | Food Research International 2026, 233, 119042; received 14 Oct 2025, accepted 19 Mar 2026 |
| DOI | 10.1016/j.foodres.2026.119042 |

## 1. Methods

**Formulation (§2.1, §2.3, p. 2-3).** Protein base: low-temperature defatted soybean meal (50.10 % protein),
soy protein concentrate (69.00 %) and soy protein isolate (90.50 %), all Yuwang (Shandong), blended
SM:SPC:SPI = 4:3:3 w/w/w. Laccase (10,000 U/g) 0.2 % w/w of total dry mix. Two precursor systems, D-xylose
+ L-cysteine (D-LC) or D-xylose + L-methionine (D-LM), xylose 99.00 %, Cys and Met 99.00 %:
- ratios D-LC or D-LM = 4:1, 3:2, 1:1, 2:3, 1:4 (read as xylose : amino acid from §3.5.3, "the ratio of
  D-xylose to L-cysteine was 3:2"); **whether the ratio is by mass or by mole is not printed**;
- doses 0.2, 0.4, 0.6, 0.8, 1.0 % w/w of total dry mix;
- **which dose the ratio series was run at, and which ratio the dose series, is not printed**;
- controls: protein base without precursors (no control bar appears in Fig. 4).
Upper limit 1.0 %: above it D-LC gave melt instability and lost fibres; D-LM above 1.0 % was extrudable but
"overly intense onion/garlic-like" and rejected; below 0.2 % changes were minor (§2.3, p. 2).

Derived here, assuming the ratio is by mass: molar xylose:amino acid is 3.23, 1.21, 0.81, 0.54, 0.20
(Cys; MW 150.13 and 121.16) and 3.98, 1.49, 0.99, 0.66, 0.25 (Met; 149.21) for 4:1 to 1:4. At the chosen
D-LC optimum (3:2, 0.8 % of dry mix) in a 60 % moisture melt: xylose 0.48 % of dry mix → 1.92 g/kg melt →
21.3 mM in the melt water; cysteine 0.32 % → 1.28 g/kg → 17.6 mM. At the D-LM optimum (1:4, 0.6 %):
xylose 5.3 mM, methionine 21.4 mM.

**Extrusion (§2.3, p. 3).** AHT36-32D co-rotating twin screw (Shandong Arrow), screw diameter 36 mm, L/D 32;
dry feed 10 kg/h; water injected to **60 % moisture (wet basis)**; barrel zones 1-8 at 30, 70, 110, 130,
**170**, 130, 60, 40 °C; screw 360 r/min; cooling slit die 1000 mm long, 70 × 6 mm cross-section, two
sections at 60 and 40 °C. Steady state judged by torque and appearance. **Residence time, melt temperature,
die pressure and SME are not printed.** §3.3 (p. 8) mentions "the hot extrusion process at 150 °C", which
does not match the 170 °C peak zone. Extrudates were frozen (−20 °C) for texture, or **freeze-dried, ground
and sieved (80 mesh) before every flavour analysis**.

**Volatiles (§2.8, p. 3-4).** HS-SPME on 1.000 g of the freeze-dried powder in a 20 mL vial, 15 min at
60 °C, fibre 40 min at 60 °C (fibre type not printed), desorbed 5 min at 250 °C. GC-O-MS on DB-WAX
60 m (identification only, Table S3, not on disk). GC×GC-SCD: DB-WAX 30 m × DB-5 2 m, modulation 5 s,
SCD 250 °C / plasma 800 °C. Identification by retention-time match to single standards plus NIST 17
(match > 700) and RI. **Quantification: "the peak area of 2-methyl-3-furanthiol and 2-furfurylthiol"
(§3.5.3, p. 9). No internal standard, no SIDA, no calibration curve, no units.** Hexanal was not measured
in any extrudate.

**Thresholds and masking (§2.6, p. 3).** 3-AFC in water, start 1 µmol/L, 1:3 (v/v) serial dilution,
sigmoid fit, threshold at 50 % detection (ASTM E1432). Binary mixtures hexanal + FFT and hexanal + MFT at
4:1, 3:2, 1:1, 2:3, 1:4 (basis not printed). Theoretical threshold from P(AB) = P(A) + P(B) − P(A)P(B);
R = experimental / theoretical; R < 0.5 synergy, 0.5-1 additive, R > 1 masking (Wang et al. 2024).

**Sensory (§2.6, p. 3).** 19 trained assessors (9 male, 10 female), ISO 8589 room, criteria in Table S1
(not on disk); attributes on radar charts with a 0-9 axis (Fig. 3): meaty, sweet, roasted, onion and
garlic, pure(ity), beany, sensory acceptance. "The average of three experiments for each indicator was
used". The direction of the beany axis (whether high means more or less beany) is defined only in Table
S1; the text reads a beany score of 7-9 as "effectively masked the beany odor". Statistics: triplicates,
mean ± SD, ANOVA and Pearson in SPSS 22, P < 0.05.

**Also in the paper, not used here:** e-nose PCA and PLSR, texture, OR51E2 docking, MD and a DFT ligand
optimisation (§2.9). By owner policy no DFT or docking output enters the repository.

## 2. Findings that matter

### 2.1 MFT and FFT, GC×GC-SCD peak area (Fig. 4c-f, p. 10; read from graph, approx.; y-axis unlabelled, scale 0 to 1.2 × 10⁹)

| series | level | MFT | FFT |
|---|---|---|---|
| D-LC, ratio (Fig. 4c) | 4:1 | 5.2e8 (a) | 4.5e8 (a) |
| | 3:2 | 8.7e8 (b) | 5.85e8 (b) |
| | 1:1 | 9.1e8 (b) | 6.45e8 (b) |
| | 2:3 | 9.3e8 (b) | 6.5e8 (b) |
| | 1:4 | 8.95e8 (b) | 6.2e8 (b) |
| D-LM, ratio (Fig. 4d) | 4:1 | 4.1e8 (d) | 2.75e8 (d) |
| | 3:2 | 5.8e8 (c) | 3.0e8 (cd) |
| | 1:1 | 6.2e8 (bc) | 3.5e8 (c) |
| | 2:3 | 6.65e8 (b) | 4.2e8 (b) |
| | 1:4 | 7.65e8 (a) | 5.4e8 (a) |
| D-LC, dose (Fig. 4e) | 0.2 % | 5.0e8 (e) | 2.3e8 (e) |
| | 0.4 % | 5.9e8 (d) | 4.0e8 (d) |
| | 0.6 % | 7.1e8 (c) | 5.8e8 (c) |
| | 0.8 % | 7.95e8 (b) | 6.55e8 (b) |
| | 1.0 % | 9.0e8 (a) | 7.3e8 (a) |
| D-LM, dose (Fig. 4f) | 0.2 % | 3.0e8 (d) | 3.85e8 (c) |
| | 0.4 % | 4.2e8 (c) | 4.15e8 (c) |
| | 0.6 % | 6.55e8 (b) | 4.85e8 (b) |
| | 0.8 % | 7.55e8 (a) | 5.25e8 (b) |
| | 1.0 % | 8.0e8 (a) | 6.3e8 (a) |

Letters as printed. What the bars say: in D-LC, MFT and FFT plateau from 3:2 on (xylose no longer
limiting below about 1.2 mol xylose per mol Cys, if the ratio is by mass); in D-LM both keep rising with
methionine share; both rise monotonically with dose in both systems, about 1.8x (MFT) and 3.2x (FFT) from
0.2 to 1.0 % in D-LC. **The abstract's claim that the highest MFT and FFT were reached "at a 3:2 ratio and
0.8% addition level" is contradicted by Fig. 4e**, where 1.0 % is higher than 0.8 % with a different
letter; 0.8 % is the sensory/texture optimum, not the yield maximum. Methionine, with no thiol of its own,
gives MFT 1.05-1.7x below cysteine at matched dose (Fig. 4e vs 4f); the paper does not discuss the route.

### 2.2 Sulfur volatiles identified (Table 1, p. 11; §3.5.2, p. 9)

68 compounds in total by retention-time match: 8 thiols, 15 thioesters, 9 thioethers, 10 thiazoles,
17 thiophenes, 1 thiane, 4 sulfhydryls, 4 heterocyclics (sum checks); 55 in D-LC, 38 in D-LM, 25 common.
Present: MFT, FFT, bis(2-methyl-3-furyl) disulfide, methyl furfuryl disulfide, difurfuryl sulfide,
dimethyl disulfide, thiazole, 2-acetylthiazole, 3-mercapto-2-pentanone, 2-mercapto-3-butanol. Not
quantified. The list also holds 2-chlorothiazole, 2-bromothiazole, 4-methyldibenzothiophene and ethyl
methanesulfonate, which are implausible Maillard products in a soy extrudate: single-standard RT matching
on SCD is weak identification, and the list should be read as tentative. GC-O (Table S3, not on disk)
found 6 odour-active new compounds (2,6-dimethylpyrazine, 2-ethyl-3,5-dimethylpyrazine, 2-acetylthiazole,
furfuryl alcohol in both; DMDS and furfural in D-LM only).

### 2.3 Water thresholds and binary masking (p. 12; Fig. 5d-m labels)

Single compounds in water (§3.6.1): hexanal **0.0040 mg/L**, FFT **0.0073 mg/L**, MFT **0.0014 mg/L**. The
authors note these are "generally higher than those reported in previous studies".

| mixture (hexanal : thiol) | FFT: threshold_exp / threshold_the (mg/L) | R | MFT: threshold_exp / threshold_the (mg/L) | R |
|---|---|---|---|---|
| 4:1 | 0.562 / 0.081 | 6.918 | 0.650 / 0.042 | 15.596 |
| 3:2 | 0.564 / 0.049 | 11.350 | 0.571 / 0.042 | 13.583 |
| 1:1 | 0.608 / 0.064 | 9.397 | 0.703 / "0.0.036" (sic) | 19.588 |
| 2:3 | 0.895 / 0.070 | 12.76 | 1.396 / 0.038 | 36.644 |
| 1:4 | 0.416 / 0.066 | 6.295 | 0.136 / 0.031 | 44.771 |

Checks (derived here). R = exp/the reproduces to within rounding in nine panels; in Fig. 5m (MFT 1:4)
0.136 / 0.031 = 4.4, not 44.771, so one of the three printed numbers is wrong by 10x. Every mixture
threshold (0.14-1.40 mg/L) is above the top of the dilution series as described (1 µmol/L hexanal =
0.100 mg/L; 1 µmol/L MFT = 0.114 mg/L), and the Fig. 5 x-axis "Log (concentration, mg/L)" runs 0.8-4.0,
i.e. 6-10,000 mg/L if taken literally. The concentration basis of the mixture thresholds cannot be
reconstructed from what is printed. Direction only: R > 1 at every ratio, i.e. the thiols raised the
detection threshold of the hexanal-thiol mixture above independent-probability addition.

### 2.4 Sensory (Fig. 3, p. 7; §3.3, p. 8, values as printed in the text)

- D-LC: meaty and roasted-nut dominant. Sensory acceptance highest at 0.8 % (7.6 points). Purity 7.1 at
  4:1 (meaty "insufficient", slight bitterness, sweet 4 points); purity fell from 7.9 to 7.2 as Cys rose
  from 1:1 to 1:4 ("pungent sulfur odor"). Acceptance vs control "generally minimal (6-8 points)".
- D-LM: sweet plus onion/garlic (3-5 points). At 1:4: sweet 7.8, meaty 4.7, roasted nut 4.3, onion and
  garlic 5. Acceptance highest at 0.6 % (7.5), falling to 6.7 beyond; attributed to DMDS/DMTS from
  methanethiol.
- "Soybean flavor score" 7-9 with precursors, read by the authors as masking. No control values on the
  radar plots; no significance tests on sensory attributes are shown.
- PLSR: meaty, roasted-nut and onion/garlic load with the W1W (sulfide) e-nose sensor; 82 % of X and
  44 % of Y variance.

## 3. What it means for the model

**What can be checked, and what cannot.** No number in this paper is an absolute concentration, so no
benchmark row in ppb can be built from it. It gives **within-study ratios and orderings** of MFT and FFT
across 20 formulations on one instrument, which is the form the engine's comparative verbs are scored on.
Three checks a run of the engine could face, once residence time is known or bounded (it is not printed;
a 36 mm, L/D 32 extruder at 10 kg/h is typically on the order of one to a few minutes, an assumption that
must be declared):

1. **Dose response.** MFT and FFT rise monotonically from 0.2 to 1.0 % D-LC, by about 1.8x and 3.2x
   (Fig. 4e). A sign and rough-slope test.
2. **Cysteine saturation.** At fixed dose, both plateau once xylose:Cys drops to 3:2 (molar ≈ 1.2 if by
   mass) (Fig. 4c): sugar becomes limiting. A shape test of the sugar/thiol stoichiometry.
3. **Cysteine vs methionine.** Matched-dose MFT from methionine is 1.05-1.7x below cysteine (Fig. 4e vs
   4f; the two dose series may sit at different ratios, which is not printed). The engine carries
   methionine as a Strecker substrate to methional (`species.py` B22, `parameters_methionine.py`); if it
   cannot make MFT from Met + xylose at all, that is a measured miss in direction.

Caveats that keep it from being a scoring row: peak areas from HS-SPME at 60 °C on freeze-dried powder
(thiol losses on drying, oxidation to the disulfide, fibre competition all unquantified); the melt
(60 % water, 170 °C peak, laccase, soy fibre) is far outside every pot on the current panels; the ratio
basis and the fixed level of each series are not printed. An MFT/FFT ratio from SCD peak areas would need
an equimolar-sulfur-response assumption and equal SPME recovery, neither tested here.

**Live values for orientation (not a comparison: different system).** The nearest engine prediction is the
held-out aqueous xylose + cysteine pot `mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5` in
`results/validation/core_prediction_uncertainty.json`: MFT predicted point 546 ppb (p5-p95 316-750,
measured 143, outside the interval), FFT predicted point 191 ppb (p5-p95 28-1722, measured 96, inside).
Ma has neither a pH nor a time to match it.

**Thresholds.** The engine's water thresholds (`WATER_THRESHOLDS`, `src/kinetic_core/matrix_oav.py`) carry
hexanal 4.5 µg/L (Guadagni via Vega 1994) and MFT 0.005 µg/L, FFT 0.006 µg/L (Zhou 2023 SI Table S2,
provenance uncited). Ma's hexanal 4.0 µg/L agrees within 1.1x; its MFT 1.4 µg/L and FFT 7.3 µg/L are 280x
and 1200x above the engine's values (derived here). The authors concede theirs are high and blame the
glass/water set-up. They should not be added as threshold records for the thiols; the hexanal value is a
fourth consistent water value at most.

**Masking.** The engine has no odour-interaction term and the matrix layer refuses every matrix-corrected
threshold. Ma's R values are the only binary-mixture data on disk for a meaty thiol against the beany
off-note, but with the concentration basis unrecoverable (section 2.3) they support only the direction
(meaty thiols suppress hexanal detection) and cannot set a size.

**For a formulator.** The actionable recipe: xylose + cysteine at 0.8 % of dry mix, about 3:2, in a
60 % moisture SM:SPC:SPI 4:3:3 HME melt with 0.2 % laccase, peak barrel 170 °C, gives meaty and
roasted-nut notes and the best acceptance (7.6/9) the panel recorded; more cysteine raises the thiols
further but weakens fibres (Cys reduces disulfides) and above 1.0 % the melt fails; extra cysteine beyond
3:2 adds sulfur pungency and no MFT. Methionine keeps texture intact but steers to sweet and onion/garlic
(DMDS) and caps at 0.6 %. Precursors added to the melt make the meat aroma in situ; the paper shows they
change the e-nose and sensory picture strongly, but it does not show how much hexanal remains.

## What it does not give

- Any absolute concentration (MFT, FFT or anything else): peak areas only, no internal standard, no SIDA.
- Any hexanal (or other off-note) measurement in the extrudates, with or without precursors.
- Control-sample values for MFT and FFT (Fig. 4 has no control bar).
- Residence time, melt temperature, pressure or SME; pH of the melt or product.
- The basis (mass or molar) of the precursor ratios and of the binary-mixture ratios; the fixed dose of
  the ratio series and the fixed ratio of the dose series.
- The sensory scale definitions (Table S1), the GC-O table (S3) and the single-compound threshold curves
  (Fig. S1): supplement not on disk.
- Any kinetics: one process condition, one time point.
