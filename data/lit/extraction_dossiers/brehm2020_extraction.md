# Brehm, Frank, Ranner & Hofmann 2020 — EXTRACTION (thiamine-derived pyrimidinylmethyl thioethers: pH, temperature and time series in buffered thiamine, thiamine/cysteine and thiamine/cysteine/ribose model systems)

**Source on disk:** `data/articles/brehm2020.pdf` (the publisher's PDF, 9 pages, journal pp. 6181-6189;
PDF page n is journal page 6180 + n) and `data/articles/brehm2020_supplementary.pdf` (the Supporting
Information, 7 pages). Both files were read and checked by eye on 2026-10-09: every table value below was
read from the rendered page image and cross-checked against `pdftotext -layout` output (both PDFs have a
text layer, and the two agreed on every number quoted here). Figure 4 has no tabulated counterpart
anywhere (see section 2.4); its values are read from a 300 dpi render and are labelled as such. Written
for the thiamine route in the sulfur lane (`src/kinetic_core/sulfur.py`, "THE THIAMINE ROUTE").

| field | value |
|---|---|
| Title | "Quantitative Determination of Thiamine-Derived Taste Enhancers in Aqueous Model Systems, Natural Deep Eutectic Solvents, and Thermally Processed Foods" |
| Authors | Laura Brehm, Oliver Frank, Josef Ranner, Thomas Hofmann (Chair of Food Chemistry and Molecular and Sensory Science, TU Munich, Freising) |
| Venue | Journal of Agricultural and Food Chemistry 2020, 68(22), 6181-6189 (received 20 Mar 2020, published 1 May 2020); funded by Lucta S.A. |
| DOI | 10.1021/acs.jafc.0c01849 |

## 1. Methods

**Aqueous model systems** (Materials and Methods, PDF p. 3 / journal p. 6183). Thiamine hydrochloride
1 mmol (systems I-III), plus L-cysteine 1 mmol (II and III), plus D-ribose 1 mmol (III), dissolved in
aqueous KH2PO4 buffer (0.1 M, 10 mL), i.e. 0.1 M of each precursor. Heated "in a closed vessel" for 30,
60, 90, 120, 180 or 240 min at 80, 120 or 180 °C at pH 3.0, 6.5 or 9.0. Stopped by cooling and dilution.
Not stated: the vessel, the heating device, the heat-up and cool-down times, whether the pH was re-measured
after heating, the number of replicates (an RSD is printed for each cell, so at least duplicates).

The three series actually reported (Results, PDF pp. 5-6):
- **pH series**: pH 3.0 / 6.5 / 9.0, all at 120 °C, 120 min (SI Table S3).
- **temperature series**: 80 / 120 / 180 °C, all at pH 6.5, 120 min (SI Table S4).
- **time series**: 30-240 min at 120 °C, pH 6.5 (Figure 4 only).

**NADES systems IV-VI** (thiamine / + cysteine / + cysteine + ribose, 1 mmol each in 5 g of NADES 1-7),
120 °C, 120 min (Table 1, PDF p. 2; SI Table S5). Each table includes an aqueous "Buffer" reference row.

**Analysis.** LC-MS/MS (QTrap 5500, Luna PFP column, ESI+, MRM; transitions in SI Table S1) with two
ethyl-pyrimidine stable analogues (15, 16) as internal standards; reference compounds quantified by qNMR.
Validation in a pork matrix, main Table 2 (PDF p. 4): recovery 78-109 %; RSD intraday 2-16 %, interday
2-13 %. Thiamine itself (compound 1): recovery 98 %, LOD 0.007 nmol/L, LOQ 0.022 nmol/L, RSD 3 / 3 %.

**What is measured.** Compound 1 is intact thiamine. Compounds 3-7 are NOT the free sulfur nucleophiles:
each is a (4-amino-2-methylpyrimidin-5-yl)methyl thioether, formed (per the authors, ref. 1-2) by S_N1
substitution of a SECOND, intact thiamine's thiazole by the nucleophile. Free HMP, free MFT and free
3-mercapto-2-pentanone were not measured anywhere in the paper.

**Units.** Every model-system value is printed as "formation rate [µmol/mmol thiamine] (RSD%)". It is not a
rate: it is an amount after the stated heating, per mmol thiamine charged. For compound 1 it is therefore
the RESIDUAL thiamine per mmol charged, so 1000 µmol/mmol would be no loss.

## 2. Findings that matter

### 2.1 Compound identities (Figure 1, PDF p. 2; confirmed by the structures and by the text, PDF p. 5)

| no. | name as printed | trapped nucleophile | model species |
|---|---|---|---|
| 1 | thiamine | — | `THI` |
| 3 | S-((4-amino-2-methylpyrimidin-5-yl)methyl)-L-cysteine | cysteine | — |
| 4 | 3-(((4-amino-2-methylpyrimidin-5-yl)methyl)thio)-5-hydroxypentan-2-one | 5-hydroxy-3-mercapto-2-pentanone | `HMP` |
| 5 | 2-methyl-5-(((2-methylfuran-3-yl)thio)methyl)pyrimidin-4-amine | 2-methyl-3-furanthiol | `MFT` |
| 6 | 3-(((4-amino-2-methylpyrimidin-5-yl)methyl)thio)pentan-2-one | 3-mercapto-2-pentanone | `MP2P` |
| 7 | 5-(((furan-2-ylmethyl)thio)methyl)-2-methylpyrimidin-4-amine | 2-furfurylthiol | — |

Figure 1 draws 4 as S on C3 of a pentan-2-one carrying a 5-OH; 5 as S on C3 of 2-methylfuran; 6 as S on
C3 of pentan-2-one with no OH. The text (PDF p. 5, journal p. 6185) states it directly: "The precursors
for compounds 4, 5, and 6, 5-hydroxy-3-mercapto-2-pentanone, MFT, and 3-mercapto-2-pentanone,
respectively, are conversion products of thiamine." Compound 7 (FFT adduct) was not found in any
pH-series system and only "minor amounts" in system III at 180 °C (no number printed).

### 2.2 Temperature series: SI Table S4 (SI p. 5), pH 6.5, 120 min, µmol/mmol thiamine (RSD %)

| system | T | 1 (residual thiamine) | 3 | 4 (HMP adduct) | 5 (MFT adduct) | 6 (MP2P adduct) |
|---|---|---|---|---|---|---|
| I | 80 °C | 537.1 (6.9) | n.det. | 2.5 (20.3) | 0.4 (10.1) | n.det. |
| II | 80 °C | 803.0 (1.3) | 27.6 (7.3) | 1.3 (30.3) | 0.1 (20.8) | n.det. |
| III | 80 °C | 904.0 (0.0) | 16.1 (15.1) | 0.7 (39.9) | 0.1 (69.1) | n.det. |
| I | 120 °C | 348.5 (6.6) | n.det. | 42.1 (14.8) | 4.6 (13.4) | 0.1 (17.5) |
| II | 120 °C | 354.3 (1.5) | 196.2 (5.2) | 21.2 (3.6) | 0.8 (1.3) | n.det. |
| III | 120 °C | 412.0 (7.2) | 157.3 (9.0) | 27.8 (7.6) | 0.9 (8.8) | 0.1 (24.1) |
| I | 180 °C | 34.7 (23.4) | n.det. | 16.0 (7.5) | 1.0 (9.0) | 0.1 (22.0) |
| II | 180 °C | 42.6 (42.8) | 2.8 (24.3) | 4.0 (10.4) | 1.7 (3.8) | 0.1 (19.9) |
| III | 180 °C | 45.1 (48.5) | 2.5 (92.9) | 4.4 (7.5) | 1.9 (1.8) | 0.1 (13.3) |

"n.det." is defined under SI Table S5 as "not detectable". The text (PDF p. 6) quotes compound 5 in
system I as 0.4 / 4.6 / 1.9 µmol/mmol at 80 / 120 / 180 °C; the table's system I value at 180 °C is
**1.0**, and 1.9 is the system III value. The table is taken here; the text's 1.9 for system I is a
misquote.

### 2.3 pH series: SI Table S3 (SI p. 4), 120 °C, 120 min, µmol/mmol thiamine (RSD %)

| system | pH | 1 (residual thiamine) | 3 | 4 (HMP adduct) | 5 (MFT adduct) | 6 (MP2P adduct) |
|---|---|---|---|---|---|---|
| I | 3.0 | 760.0 (11.0) | n.det. | 4.3 (13.5) | 0.3 (18.5) | n.det. |
| II | 3.0 | 680.1 (4.0) | 39.8 (6.7) | 1.7 (19.3) | 0.1 (24.0) | n.det. |
| III | 3.0 | 392.4 (1.0) | 13.2 (8.0) | 1.1 (7.0) | 0.1 (11.2) | n.det. |
| I | 6.5 | 348.5 (6.6) | n.det. | 42.1 (14.8) | 4.6 (13.4) | 0.1 (17.5) |
| II | 6.5 | 354.3 (1.5) | 196.2 (5.2) | 21.2 (3.6) | 0.8 (1.3) | n.det. |
| III | 6.5 | 312.0 (7.2) | 157.3 (9.0) | 27.8 (7.6) | 0.9 (8.8) | 0.1 (24.1) |
| I | 9.0 | 106.7 (4.1) | n.det. | 0.8 (6.2) | 0.2 (5.7) | n.det. |
| II | 9.0 | 45.5 (32.3) | 17.5 (22.2) | 0.4 (22.2) | 0.1 (17.7) | n.det. |
| III | 9.0 | 56.6 (9.3) | 14.6 (6.7) | 0.9 (2.8) | 0.2 (9.4) | n.det. |

**Internal inconsistencies, all read by eye:**
- The pH 6.5 rows of Table S3 and the 120 °C rows of Table S4 are the same condition and agree to the
  digit in every cell (including the RSDs) except residual thiamine in system III: **312.0 (7.2)** in
  Table S3, **412.0 (7.2)** in Table S4. Identical RSD makes it one measurement and one typo; which is
  right cannot be told from the paper. Both are carried in section 3.
- The text (PDF p. 5) gives compound 4 in system I at pH 9.0 as 0.9 µmol/mmol; Table S3 has 0.8 for
  I-pH9.0 (0.9 is III-pH9.0).
- The text (PDF p. 5) gives compound 6 in system III at pH 6.5 as "<0.02"; Table S3 has 0.1 (24.1).
- The aqueous reference rows of SI Table S5 (SI p. 6), nominally the same recipe (120 °C, 120 min,
  pH 6.5, per PDF p. 6: "All yields were compared to aqueous buffered systems at optimized parameters"),
  differ from Tables S3/S4:

| row (Table S5) | 1 | 3 | 4 | 5 | 6 | same condition in S3/S4 |
|---|---|---|---|---|---|---|
| Buffer-IV (thiamine only) | 189.5 (16.4) | n.det. | 35.9 (15.8) | 4.1 (3.7) | 0.08 (0.9) | system I: 348.5 / n.det. / 42.1 / 4.6 / 0.1 |
| Buffer-VI (thiamine + cysteine + ribose) | 140.0 (16.4) | 175.6 (6.1) | 10.1 (5.4) | 0.5 (7.5) | n.det. | system III: 312.0 or 412.0 / 157.3 / 27.8 / 0.9 / 0.1 |

  So a nominal replicate of system I lost 81 % of its thiamine rather than 65 %, and of system III 86 %
  rather than 59-69 %. This between-batch spread is far larger than the printed RSDs, and is the honest
  error bar on any single residual-thiamine value here.

### 2.4 Time series at 120 °C, pH 6.5: Figure 4 (PDF p. 6), read from graph, approx.

The caption says "compounds 1, 3, 4, 5, and 6" and "Quantitative data are disclosed in the Supporting
Information". Neither is true: **the plotted panels contain no compound 1 and no compound 6**, and the SI
has no time-series table (SI Tables S1-S5 are MS parameters, foods, pH, temperature, NADES). **There is no
time course of residual thiamine anywhere in the paper or SI.** What is plotted (µmol/mmol thiamine,
read from graph, approx., about ±5 % of each axis span):

| t (min) | I: 4 | I: 5 | II: 3 | II: 4 | II: 5 | III: 3 | III: 4 | III: 5 |
|---|---|---|---|---|---|---|---|---|
| 30 | 8 | 0.5 | 68 | 3.6 | 0.1 | 74 | 4.2 | 0.1 |
| 60 | 21 | 1.4 | 90 | 6.1 | 0.3 | 139 | 9.2 | 0.2 |
| 90 | 36 | 2.4 | 155 | 11 | 0.3 | 165 | 14 | 0.4 |
| 120 | 49 | 3.9 | 222 | 17.3 † | 0.5 | 185 | 16 | 0.4 |
| 180 | 60 | 5.3 | 255 | 20.2 † | 0.7 | 192 | 21 | 0.6 |
| 240 | 67 | 6.0 | 280 | 21.1 † | 0.9 | 220 | 24 | 0.7 |

Axes: panel I, 4 on left (0-80), 5 on right (0-7); panels II and III, 3 on left (0-350 / 0-300), 4 and 5
on right (0-25 / 0-30). † printed in the text (PDF p. 6): "in system II, the concentration of compound 4
increases from 17.3 at 120 min to 20.2 at 180 min to 21.1 at 240 min", which matches the graph and fixes
the axis assignment. The 120 min points are a third batch again: system I compound 4 ≈ 49 in Figure 4
against 42.1 (Table S3/S4) and 35.9 (Table S5); system III compound 4 ≈ 16 against 27.8 and 10.1.

Shape: every adduct rises roughly linearly to 90-120 min and then flattens ("after 120 min, the curve was
flattening", PDF p. 6); none falls by 240 min at 120 °C. Compound 5 in panel I is sigmoid-ish (lag over
the first 30 min, approx.), consistent with a consecutive thiamine → MFT → adduct sequence, though the
error bars do not exclude a straight line.

### 2.5 Heated foods: SI Table S2 (SI p. 3), nmol/kg (RSD %), meat rows only

| sample | 1 | 4 | 5 | 7 |
|---|---|---|---|---|
| pork pan (5 min, max. heat) | 4572 (0.9) | n.det. | 0.4 (8.4) | n.det. |
| roasted pork (180 °C oven, 120 min) | 4363 (7.9) | n.det. | 0.2 (19.2) | 8.7 (29.5) |
| beef pan | 825 (10.3) | n.det. | 0.4 (21.6) | 17.7 (10.9) |
| roasted beef | 207.0 (2.6) | n.det. | 0.2 (19.0) | 18.3 (5.3) |
| roasted beef sauce | 224.2 (1.7) | n.det. | 0.6 (10.9) | 5.6 (38.0) |
| chicken pan | 283.2 (0.6) | n.det. | 0.4 (15.0) | 11.6 (44.0) |

Compound 6 was detected in no food ("could not be detected in any of the food samples", PDF p. 8).
Compound 4 only in egg yolk (2.0), egg white (2.3) and Brazil nut (5.2) nmol/kg. No unheated control,
so none of these is a conversion.

## 3. What it means for the model

### 3.1 Crude first-order constant for total thiamine loss (derived here)

Assumptions, all of them strong: (i) compound 1 is the residual of 1000 µmol charged per mmol, so
f = value / 1000; (ii) loss is first order in thiamine with one constant; (iii) the full 120 min is
isothermal at the set temperature (heat-up is not reported, so the effective time is shorter, most of all
at 180 °C); (iv) one time point per temperature, no zero-time measurement of thiamine after dissolution.
k = −ln(f) / 120 min.

| system | T | value (Table S4) | f | −ln f | k (min⁻¹) | k (s⁻¹) |
|---|---|---|---|---|---|---|
| I | 80 °C | 537.1 | 0.5371 | 0.622 | 5.18e-3 | 8.6e-5 |
| I | 120 °C | 348.5 | 0.3485 | 1.054 | 8.78e-3 | 1.46e-4 |
| I | 180 °C | 34.7 | 0.0347 | 3.361 | 2.80e-2 | 4.67e-4 |
| II | 80 °C | 803.0 | 0.8030 | 0.219 | 1.83e-3 | 3.0e-5 |
| II | 120 °C | 354.3 | 0.3543 | 1.038 | 8.65e-3 | 1.44e-4 |
| II | 180 °C | 42.6 | 0.0426 | 3.156 | 2.63e-2 | 4.38e-4 |
| III | 80 °C | 904.0 | 0.9040 | 0.101 | 8.41e-4 | 1.4e-5 |
| III | 120 °C | 412.0 (S4) / 312.0 (S3) | 0.412 / 0.312 | 0.887 / 1.165 | 7.39e-3 / 9.71e-3 | 1.23e-4 / 1.62e-4 |
| III | 180 °C | 45.1 | 0.0451 | 3.099 | 2.58e-2 | 4.30e-4 |

Nominal-replicate spread at 120 °C (derived here, same formula): system I 8.78e-3 (Table S4) vs
1.39e-2 min⁻¹ (Buffer-IV, Table S5, f = 0.1895); system III 7.4e-3 to 9.7e-3 vs 1.64e-2 min⁻¹
(Buffer-VI, f = 0.140). A factor of about 1.6-2 between batches.

pH series at 120 °C (derived here): k = 2.29e-3 (I, pH 3.0), 8.78e-3 (I, pH 6.5), 1.87e-2 min⁻¹ (I,
pH 9.0); II: 3.21e-3 / 8.65e-3 / 2.58e-2; III: 7.80e-3 / 9.71e-3 (S3 value) / 2.39e-2. Thiamine loss
rises monotonically with pH (the authors say so, PDF p. 5), while every adduct peaks at pH 6.5.

### 3.2 Apparent activation energy (derived here)

Least-squares fit of ln k against 1/T with T = 353.15, 393.15, 453.15 K, R = 8.314 J mol⁻¹ K⁻¹:

| system | Ea, 3-point fit (kJ/mol) | residuals of ln k (80 / 120 / 180 °C) | 2-point 80→120 | 2-point 120→180 |
|---|---|---|---|---|
| I (thiamine alone) | **22.6** | +0.09 / −0.17 / +0.08 | 15.2 | 28.6 |
| II (+ cysteine) | **35.3** | −0.12 / +0.22 / −0.10 | 44.8 | 27.5 |
| III (+ cysteine + ribose), S4 412.0 | **45.2** | −0.21 / +0.40 / −0.18 | 62.7 | 30.9 |
| III, S3 312.0 | 45.0 | −0.31 / +0.58 / −0.27 | 70.6 | 24.2 |

How far to trust these: not far. The 180 °C point almost certainly spent a substantial part of its
120 min heating up (unreported), which lowers its effective k and flattens the slope; the 120 °C batch
spread alone (factor ~1.6-2 in k) moves a two-point Ea between 80 and 120 °C by roughly ±10-15 kJ/mol;
and system I's 46 % loss at 80 °C against II's 20 % and III's 10 % is not explained in the paper (any
reaction of thiamine with cysteine would make II lose MORE, not less). The curvature (the 80→120 slope
steeper than 120→180 in II and III) is the signature of a heat-up-limited high-temperature point.
My reading: these data say the total-thiamine-loss Ea in 0.1 M phosphate at pH 6.5 is **probably below
~60 kJ/mol over 80-120 °C (confidence ~60 %)**, and they cannot pin it closer than 20-70 kJ/mol.

What it is an Ea OF: total loss of intact thiamine by every route (thiazole release, pyrimidine
substitution, hydrolysis, the adduct-forming S_N1 itself). It is NOT the Ea of the model's `r_thi_hmp`
step (thiamine → HMP), which is one branch of that loss. The sulfur lane books `k_thi_hmp` to the
`thiol_assembly` formation-Ea route, prior centre 100 kJ/mol, band (55, 145) kJ/mol
(`parameters_sulfur.FORMATION_EA_BOUNDS_BY_ROUTE`, `FORMATION_EA_PRIOR_CENTRE`; the centre is from Chan &
Reineccius 1994, a class prior). That band is the search space, not the value the engine runs: B10's route
split was re-merged ("barriers not identified"), so the live engine gives `k_thi_hmp` the frozen
`lumped_formation_Ea_kJ_mol` = 64.08 kJ/mol from `results/validation/kinetic_core_b9_fit_report.json`,
with its rate anchored at 145 °C and not sampled by the envelope. All three apparent values here sit below
the band's floor, and below the live 64. This is a
flag, not a correction: a branch can carry a higher Ea than the total, but the total cannot be slower
than any one branch, so if `k_thi_hmp` at 80-120 °C exceeds the total loss rate derived above, the model
is wrong.

### 3.3 Bearing on "thiamine adds MFT at 100 °C/20 min but little at 140 °C/5 min"

(a) **How much thiamine converts at the two sweep conditions** (derived here, extrapolating the 3-point
fits of 3.2 to 100 and 140 °C; same assumptions): k·t at 100 °C/20 min vs 140 °C/5 min is 0.143 vs 0.073
(system I, i.e. 13 % vs 7 % of thiamine lost), 0.078 vs 0.059 (II), 0.047 vs 0.049 (III). The ratio
k·t(100 °C, 20 min) / k·t(140 °C, 5 min) = 4 · exp[−Ea/R · (1/373.15 − 1/413.15)] is 1.98 at 22.6, 1.33 at
35.3, 0.98 at 45.2, 0.72 at 55, 0.18 at 100 and 0.04 at 145 kJ/mol. So with this paper's low apparent Ea
the two sweep points convert about the same thiamine (100/20 up to twice as much); with the lane's prior
centre (100 kJ/mol) 140/5 converts about 5.6 times more; at the live engine value (64.08 kJ/mol) the ratio
is 0.54, so 140/5 converts about 1.8 times more. If the sweep's contrast is driven by the
thiamine branch itself, it rests on a low effective Ea, which this paper weakly supports; if the lane runs
`k_thi_hmp` at 64 kJ/mol or above, the contrast must instead come from the competing sugar route outgrowing the
thiamine route at 140 °C, and this paper says nothing about that.

(b) **Yield of the MFT adduct per thiamine consumed** (derived here: compound 5 / (1000 − compound 1),
Table S4): system I 0.0009 (80 °C), 0.0071 (120 °C), 0.0010 (180 °C); system II 0.0005 / 0.0012 /
0.0018; system III 0.0010 / 0.0015 / 0.0020. Same for the HMP adduct (compound 4): I 0.0054 / 0.0646 /
0.0166; II 0.0066 / 0.0328 / 0.0042; III 0.0073 / 0.0473 / 0.0046. In thiamine alone the partition of lost
thiamine into the HMP and MFT adducts peaks sharply at 120 °C; at 80 °C much of the loss goes elsewhere
(or into HMP not yet trapped). With cysteine and ribose present (II, III, the meat-like case), the MFT
adduct per thiamine lost is small (0.1-0.2 %) and rises weakly with temperature, while the HMP adduct
collapses at 180 °C. The authors' reading (PDF p. 6): at 180 °C "either thiamine degrades to other
compounds than 5-hydroxy-3-mercapto-2-pentanone or MFT or the taste modulators 3, 4, 5, and 6 are degrading
themselves".

(c) **The caveat that dominates (b).** Compound 5 needs both free MFT and a second intact thiamine. At
180 °C only 3.5-4.5 % of thiamine survives, so the adduct can fall because the trap is gone, not because
less MFT formed. The adducts are a lower bound on the nucleophile flux and a convolution of nucleophile
formation with residual thiamine; they cannot be converted to free MFT without a trapping-rate model the
paper does not give. Adding cysteine (II) cut the MFT adduct at 120 °C from 4.6 to 0.8 µmol/mmol, which
the authors attribute to cysteine competing as nucleophile (compound 3 at 196.2). In a meat-like pot,
cysteine is present, so system II/III numbers, not system I, are the relevant ones.

(d) **Qualitative support.** The paper supports two things the route assumes: thiamine alone at pH 6.5
does produce HMP and MFT (as their adducts) with HMP about 9-fold ahead of MFT at 120 °C (42.1 vs 4.6,
Table S4), and 3-mercapto-2-pentanone (via 6) at about 1/400 of HMP (0.1 vs 42.1). If adduct ratios
tracked free-thiol ratios (unproven: the three thiols' trapping rates may differ), the branch ratio
k_hmp_mft : k_hmp_mp2p would be about 46:1 in system I at 120 °C (4.6 / 0.1, derived here; 0.1 is at the
table's resolution limit so this is ±50 % or worse).

## What it does not give

- **No time course of thiamine (compound 1).** Figure 4's caption lists it; the panels do not plot it and
  the SI has no time-series table. Only single 120 min points exist, so first-order behaviour cannot be
  tested and every k above is a one-point estimate.
- **No heat-up / cool-down profile, vessel or heating device**, so the effective isothermal time at 180 °C
  (and to a lesser degree 120 °C) is unknown.
- **No free HMP, MFT or 3-mercapto-2-pentanone.** Only the pyrimidinylmethyl thioethers 4, 5, 6, which
  require a second intact thiamine. No step rate constants, no per-step Ea, no kinetic fit of any kind.
- **No 100 °C or 140 °C points, no time series other than at 120 °C, no time series at other pH.**
- **No pH re-measurement after heating**, no buffer capacity check at 180 °C (0.1 M KH2PO4 against
  0.1 M thiamine hydrochloride).
- **No replicate count** behind the RSDs, and the between-batch spread (Table S3/S4 vs Table S5 vs
  Figure 4, section 2.3-2.4) is several-fold larger than the printed RSDs.
- **Compound 7 at 180 °C in system III**: described as "minor amounts", no number.
- **Foods**: no unheated controls, no time/temperature series; the meat rows give no conversion.
- **Minor labelling clash**: SI Table S1 names 13 as the pyrimidinylmethanol and 14 as the aminomethyl
  pyrimidine; Figure 1 and main Table 2 have them the other way round. Not used here.
