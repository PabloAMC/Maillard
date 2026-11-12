# Why laboratories disagree: method and inter-laboratory papers, read against the seven sign-conflict cells

Reader's note on files. `data/articles/Zhou2023.pdf` is NOT the zhou2023 source. It is Zhu et al. 2023, Food Chem 400:133854 (kafirin and sulfur odorants in baijiu); `data/keys/pdf_dossier_aliases.yml` agrees. The conflict-cell source zhou2023 (Zhou, ..., Ho, JAFC 71:2472, 2023) is `Zhou2023b.pdf`, and I read its methods there. I also used the Zhu kafirin paper, because it measures protein binding of 2-methyl-3-furanthiol (MFT) and bis(2-methyl-3-furyl) disulfide (BMFD), the two analytes in the pH cells. Dex conversions use RSD/2.3 throughout.

---

## Part A. The five method papers

### A1. Roberts, Pollien & Milo 2000, J Agric Food Chem 48:2430 (SPME method development; Nestle)
**Read in full**, with Figures 6 and 8 read from the page images.
- **What and how:** HS-SPME (PDMS, PDMS/DVB, CAR/PDMS and others) on coffee brew and on aqueous flavour solutions. Relative peak-area quantitation is compared with isotope dilution (SPME-IDA) during coffee storage.
- **Numbers:**
  - **Repeatability:** CV about 7% (range 2-15%), so about 0.03 dex (0.01-0.06). (p. 2433)
  - **Linear range:** usually under 1 ppm. Table 1: polar compounds (2,3-pentanedione, pyrazine, furfural, 2-acetylfuran) gave 56-60% of full area at 50% strength instead of 50%, which is +0.05 to +0.08 dex of compression.
  - **Competition (Fig. 6, PDMS, 30 min):** adding 2-isobutyl-3-methoxypyrazine had little effect up to 1 ppm. From my reading of the figure, at about 100 ppm guaiacol fell from about 4.0e5 to about 2.7e5 (about -0.17 dex) and 3-methylbutanal from about 1.6e5 to about 0.4e5 (about -0.6 dex). At 1000 ppm the low-affinity compounds went to nearly zero (more than 1 dex). The authors' conclusion: if one sample of a pair has a large load of a high-K(fibre-air) compound, the two "could not be adequately compared unless isotopic standards were used". That is a lever-dependent bias whenever the lever raises total volatile load.
  - **pH-dependent bias (the key mechanism):** during storage the coffee pH fell from 5.25 to 4.85. Relative SPME showed acetic acid rising, while IDA showed no change. The cause is a larger undissociated fraction partitioning into the headspace. The same mechanism applies to any weak acid, including a thiol whose pKa is near the lever's pH range. The figure gives no number for the size of the acetic-acid artefact; I read it only qualitatively.
  - **Tray standing time (Fig. 10):** low-volatility compounds and alkylpyrazines rise slightly with time in the autosampler tray, and N-methylpyrrole falls. Samples run in lever order therefore carry an order-correlated bias.
- **Bearing:** MFT and BMFD x pH (a pH-dependent partition bias in non-IDA SPME), and hexanal x sugar (competition). See Part C.

### A2. Aparicio-Ruiz et al. 2022, Food Control 135:108756 (5-laboratory peer validation, harmonised SPME-GC-MS, virgin olive oil)
**Read in full (text and all tables).**
- **What and how:** five laboratories ran the same 15 oils in duplicate. Fibre DVB/CAR/PDMS, 40 C, 40 min. Internal standard (IS) 4-methyl-2-pentanol. Three quantitation modes: QM1 (calibration ratioed to the IS), QM2 (external calibration, no IS), QM3 (IS-slope ratio).
- **Hexanal numbers:**
  - **Repeatability RSD (Table 4):** 5.5 / 4.8 / 4.5% (QM1 / QM2 / QM3), about 0.02 dex.
  - **Reproducibility RSD across laboratories (Table 5):** 23.0% QM1 (0.10 dex), 25.8% QM2 (0.11 dex), 53.9% QM3 (0.23 dex).
  - **Means over the 18 compounds:** 30.9 / 48.0 / 55.4%, i.e. 0.13 / 0.21 / 0.24 dex.
  - **Recovery of a 5 mg/kg spike (Table 6):** 106 ± 11% (QM1), 150 ± 67% (QM2), 266 ± 221% (QM3). In the low calibration range: 119 ± 42%, 114 ± 53%, 188 ± 142%.
  - **Example of the absolute scatter:** sample S2 hexanal ranged 4.33-7.47 mg/kg across laboratories, a factor of 1.7 (0.24 dex).
- **Competition in real samples (section 3.5):** IS area RSD averaged 9.66% across samples within a laboratory. Two lampante oils with heavy defects (high volatile load) had significantly lower IS areas. Without them the RSD was 7.15%. The authors attribute recovery deviations to competition that differs between the calibration mixture and the sample, which is a sample-dependent bias. Its size in real oils is modest: the depressed IS areas account for a few percent, roughly 0.01-0.03 dex.
- **Bearing:**
  - Absolute between-laboratory scatter is about 0.10 dex with an IS and 0.11-0.23 dex without one, even under a harmonised protocol on identical samples. That is the same order as the pooled between-laboratory SD (0.31 dex), but it is a scatter in absolute values, and a constant per-laboratory part cancels in a ratio.
  - The sample-dependent (competition) part is about 0.03 dex or less in olive oil. It cannot by itself produce 0.5-1 dex sign flips.
  - Hexanal x sugar cells: a quantitative bound (Part C).

### A3. Casadei et al. 2021, Food Control 123:107823 (3-laboratory peer validation, same protocol with FID, same 15 oils)
Read: abstract, methods sections 2.3-2.7, all of the results (3.1-3.7) and Tables 5-9. Introduction and references skimmed.
- **Hexanal numbers:**
  - **Repeatability (Table 5):** 7.1 ± 2.3 / 6.9 / 6.3% (QM1 / QM2 / QM3), about 0.03 dex.
  - **Reproducibility (Table 6, QM1):** 28.0%, about 0.12 dex. Mean over the 18 compounds 38.8% (0.17 dex). Ethyl propanoate reached 122% near its detection limit.
  - **Recovery (Table 7):** 91% (laboratory range 69-104) with QM1; 140% (93-181) with QM2. The QM2 range between laboratories is log10(181/93) = 0.29 dex.
  - **Manual integration by four operators on one chromatogram:** up to 7% RSD (0.03 dex), higher for small peaks among crowded ones.
- **FID vs MS on the same oils:** for hexanal in S2, FID read 5.14-6.96 and MS read 4.33-7.47 mg/kg, so detector choice adds little beyond the laboratory spread.
- **Bearing:**
  - Without an IS, a laboratory's hexanal response can be off by 0.1-0.3 dex. That is constant within a run, so it cancels in a ratio unless it is matrix-dependent.
  - The authors say competition "may be also different for the analyte in the calibration mixture and in a given sample". This is a mechanism statement; the paper does not measure how much the competition varies from sample to sample.

### A4. Cordero et al. 2019, Anal Chim Acta (accepted manuscript; MHS-SPME on cocoa nibs, cocoa mass and powder, parallel MS and FID; Torino and Firmenich)
Read: abstract, introduction, methods 2.1-2.7, results 3.2-3.4, Tables 2-3 and conclusions. Section 3.1 (cocoa background) skimmed.
- **What and how:** multiple headspace extraction (4 extractions) gives total analyte amounts in solid matrices. External calibration on MS. FID relative response factors cross-check MS within ±20% (acetoin +24%).
- **Precision (Table 2):** RSD 1.3-10.1% (benzaldehyde 2.4% MS / 5.1% FID; trimethylpyrazine (TMP) 3.1% / 4.2%), i.e. 0.006-0.044 dex.
- **Matrix-dependent release (Table 3, Figs. 3-4):**
  - The decay constant β (the fraction left after each extraction) depends on the physical state of the matrix, not on cocoa origin or roast. Benzaldehyde: nibs 0.37, mass 0.57, powder 0.55. TMP: nibs 0.32 ± 0.07, mass 0.58 ± 0.03, powder 0.50.
  - **Consequence for single-shot HS-SPME:** the first-extraction area is proportional to amount × (1 − β). For the same amount, mass gives 0.43/0.63 of the nibs response for benzaldehyde (about −0.17 dex) and 0.42/0.68 for TMP (about −0.21 dex).
  - The authors state that normalised indicators (percent area) and absolute amounts "lead to apparently different conclusions" (Fig. 6A vs 6B).
- **Bearing:**
  - A lever that changes the physical matrix (fat dispersion, solids, a dried versus a liquid sample) induces about 0.2 dex of bias in single-shot HS-SPME, and a ratio does not cancel it.
  - This is relevant where one side of a contrast changes matrix state: browning or melanoidin load with sugar identity, or the freeze-dried versus liquid samples in Chen 2023.
  - It cannot produce the observed 0.5-1.3 dex flips alone.

### A5. Eckert et al. 2023, Metabolomics 19:62 (5-laboratory ring trial, PDMS tubes and TD-GC-MS, clonal Tanacetum vulgare, jasmonic acid treatment)
Read: abstract, methods, results and discussion in full. Figure 5 read from the page image. Supplementary tables not available.
- **What and how:** clonal plants were grown in five donor laboratories under one protocol, sampled before and after jasmonic acid (JA). Tubes were exchanged with four recipient laboratories for analysis.
- **Numbers:**
  - **Internal standard:** 1-bromodecane was detected in only 45% of samples (167/370). Its quantity agreed poorly between recipient laboratories: intraclass correlation near 0 for the mono-chemotype, below 0.5 for the mixed one.
  - **Common detections:** only 8 of 70 VOCs were detected by all recipient laboratories. Peaks per sample ranged from 11-13 (L4) to 33-34 (L2).
  - **Within-study treatment effects were not reproducible (Fig. 5):** the median log2 fold change for the same compound ranged from about −2 to about +7 across donor laboratories, i.e. about −0.6 to +2.1 dex, with sign flips. Examples: eucalyptol about +7 log2 in L2 and about 0 elsewhere; β-thujone about −2 in L1 and L5. The LMM test for donor laboratory × compound gave χ² = 52 and 61 (p < 0.01 and < 0.001).
  - The authors attribute this mainly to pre-analytical realisation of the lever: JA poorly dissolved in some laboratories, plant size differing by laboratory, and a dbRDA effect of plant height (F = 6.8, p < 0.001). Analytical differences come second.
- **Bearing:** this is the general precedent for every cell. A within-study ratio can flip sign between laboratories by more than 1 dex when the lever is applied differently, even under a shared protocol. The measurement step is not needed to explain it. See the MFT and BMFD x pH cells, where the "pH" lever was realised very differently.

### A6. (Extra, used for thiol binding) Zhu et al. 2023, Food Chem 400:133854, the file `Zhou2023.pdf`
Read: abstract, methods 2.1-2.4, results 3.2 and 3.4, Tables 1-2.
- **What:** HS-SPME-GC×GC-SCD of sulfur odorants in 53% ethanol, with and without 200 mg/L kafirin.
- **Numbers (Table 1):** headspace loss on adding protein was BMFD −64% (−0.45 dex), MFT −55% (−0.34 dex), methyl 2-methyl-3-furyl disulfide −37% (−0.20 dex) and dimethyl trisulfide −22% (−0.11 dex).
- **Binding constants (Table 2):** MFT K fell from 1093 to 848 to 404 M⁻¹ at 25, 45 and 65 C; BMFD from 1934 to 1214 M⁻¹ over the same range.
- **Bearing:**
  - Protein or peptide binding biases free-headspace MFT and BMFD by 0.3-0.45 dex. It depends on the sampling temperature and on how much free protein or peptide is present, which a lever that consumes protein (browning) would change.
  - It does not bear on the two pH cells: both xu2007 and zhou2023 are protein-free model systems (ribose or ARP with cysteine).

---

## Part B. Methods of the conflict-cell sources (methods sections only)

| source | system | pH control | heat | extraction / quantitation | noise |
|---|---|---|---|---|---|
| xu2007 (Xu, Gao et al. JSFA 88:328) | ribose + cysteine | **0.2 M phosphate/pyrophosphate buffer**. Final pH (Table 1): 5.6→5.5, 7.0→6.2 (control) or 6.0 (SC-CO2), 8.0→6.8 or 6.3 | 140 C, with and without SC-CO2 at 40 MPa | DVB/CAR/PDMS, 60 C, 20 min; area vs tridecane, response factor 1 ("approximate quantities") | mean CV < 20%, up to 47% (0.09-0.20 dex) |
| zhou2023 = Zhou2023b.pdf (Zhou ... Ho, JAFC 71:2472) | xylose-alanine ARP 20 mM + cysteine 20 mM | **not buffered:** initial pH set with 3 M NaOH. **Final pH (Fig. 2, ARP+Cys): about 3.3 from initial 6, about 3.4 from initial 7, about 5.0 from initial 8** | 120 C, 1 h | CAR/PDMS/DVB, 60 C, 30 min, NaCl; calibration curves against 1,2-dichlorobenzene IS (Table S1; matrix of the standards not stated) | MFT SD about 2.5% (0.01 dex) |
| chen2023 (Lycium seed meal) | defatted seed-meal hydrolysate (1.48% fat) + cysteine + sugar (0.3 g per 1.0 g hydrolysate) | **not buffered:** adjusted to 7.5; xylose gave the largest pH drop | 120 C, 120 min | (a) CAR/PDMS 50 C, 30 min, o-dichlorobenzene IS, ng/g (SI Table S1). (b) GC-IMS on 0.05 g **freeze-dried** product in 2 mL water, raw peak area, no IS (SI Table S2) | benzaldehyde SD 10-17% (0.04-0.07 dex) |
| zhang2024b (chicken carcass hydrolysate) | hydrolysate (1.25% fat on dry basis), **no added cysteine**, sugar dose 0.5-3.5% | **natural pH 6.62, no adjustment**; final pH 5.01 (xylose 3.5%) to 5.67 (glucose 3.5%) | 100 C, 60 min | CAR/PDMS 60 C, 30 min; **GC-FID area, no IS**; also relative peak area (%) | n = 4 |
| theng2024 (yeast hydrolysate) | nutritional-yeast hydrolysate + sugar | adjusted to 6.5, "no external pH control" | 100 C, 2 h | DVB/CAR/PDMS 60 C, 30 min, 2-methyl-3-heptanone IS curve ("relative concentration, ppm") | benzaldehyde SD 7-25% |
| hwang1995 (Hwang, Hartman, Ho JAFC 43:179) | glucose + ¹⁵N-glycine (the control) ± a second amino acid | adjusted to pH 7 (HCl or NaOH), then freeze-dried and brought to 12-14% moisture | 180 C, 1 h, low-moisture solid | purge-and-trap (Tenax/Carbotrap, 80 C, 1 h), d8-toluene IS, µg/g glucose | single values, no SD given |
| park2022 (mealworm) | mealworm protein 2 g + xylose 0.1 g + water 2 mL, ± 0.1 g amino acid | **no pH adjustment reported** | 125 C, 60 min | Likens-Nickerson SDE for 2 h, benzyl acetate IS, peak area relative to IS | n = 3 |

**Metadata errors in the conflict table that matter:**
1. zhou2023 and chen2023lycium are recorded `buffered = True`, but both were only pH-adjusted. Zhou's buffering is what decides the pH cells. hwang1995's `buffered = True` is also only an adjustment.
2. The matrix labels in the TMP cell appear swapped. hwang1995 is a glucose/glycine model system but is recorded as "food or protein matrix". park2022 uses mealworm protein but is recorded as "model system".
3. Effect values against the source tables:
   - **Reproduced:** Zhou per initial-pH unit 6→7 (MFT log10(1588.6/697.0) = +0.36; BMFD log10(102.6/59.7) = +0.24); Chen hexanal (GC-MS +0.46 / +0.33 / +0.34 for xylose vs fructose / arabinose / glucose; GC-IMS +0.07 / +0.11 / +0.46); Chen benzaldehyde log10(87.6/211.3) = −0.38.
   - **Xu's rows do not all reproduce under one normalisation.**
     - MFT, control: log10(2.63/16.62) = −0.80 dex for initial pH 5.6→7.0, i.e. −0.57 per initial-pH unit or −1.14 per final-pH unit. The recorded −0.96 matches neither.
     - MFT, SC-CO2: −0.12 dex, i.e. −0.27 per final-pH unit. This matches the recorded −0.27.
     - BMFD, control: −0.51 dex over 0.7 final-pH units, i.e. −0.72 per unit (recorded −0.71, per final pH).
     - BMFD, SC-CO2: +0.12 dex over 1.4 initial-pH units, i.e. +0.09 per unit (recorded +0.09, per initial pH).
     - So the Xu rows mix final- and initial-pH normalisation. Worth auditing.
   - **TMP orientation:** the cell's signs are oriented arginine→none. Read none→arginine, hwang1995 is log10(73.43/133.57) = −0.26 (Table 1) and park2022 is log10(0.250/0.042) = +0.77 (Table, compound 8). The sign conflict is real either way.

---

## Part C. The seven cells: best-supported explanation and confidence

**1. MFT x pH (xu2007 −0.96 / −0.27 vs zhou2023 +0.36).**
- **Explanation: the "pH" lever was realised differently in the two laboratories.**
  - Xu held pH with 0.2 M phosphate, so the realised contrast was 5.5 vs 6.2 (control) or 5.5 vs 6.0 (SC-CO2).
  - Zhou was unbuffered. Its initial-pH 6 and 7 runs both finished at about pH 3.3-3.4 (Fig. 2). Its "6→7" contrast is therefore a difference in early alkalinity (how much NaOH), at an operating pH near 3.3, not a sustained pH difference.
  - Zhou's own 7→8 step, where the final pH really differs (about 3.4 vs 5.0), gives MFT 1588.6→525.6 µg/L, **−0.48 dex: the same sign as Xu**. BMFD gives −0.31.
  - Eckert 2023 shows this failure mode, a lever realised differently, producing sign flips of up to 2 dex under a shared protocol.
- **A secondary bias pushing Xu negative.** Neither laboratory used isotope dilution. In Xu the sampled vial is at the final pH, which rises with the lever. Higher pH means more thiolate and less MFT in the headspace (Roberts 2000, acetic-acid mechanism). In Zhou both vials sat at about 3.3, so there was no such bias. Its size depends on MFT's thiol pKa, which none of my papers gives, so I cannot size it. Protein binding (Zhu 2023) does not apply: neither system contains protein.
- **What would confirm:** re-key Zhou's effect per realised (final or time-averaged) pH, or drop its 6→7 step from the "pH 6-7.5" band; check that its 7→8 step pools with Xu's sign. A run of Zhou's system in 0.2 M phosphate at pH 6 and 7 should then give a negative MFT slope. To size the pH-partition bias, spike MFT at a fixed amount into the heated Xu matrix at pH 5.5 and 6.8 and compare SPME areas. A difference above about 0.1 dex supports it; a flat response refutes it.
- **What would refute:** a buffered Zhou-type run giving a positive slope from pH 6 to 7.
- **Confidence: medium-high.** The realised-pH difference is a fact in the papers. That it is the whole story is less certain, because the two systems also differ in sugar (ribose vs ARP), temperature (140 vs 120 C) and SC-CO2.

**2. BMFD x pH (xu2007 −0.71 / +0.09 vs zhou2023 +0.24).**
- **Explanation:** the same realised-pH problem as cell 1. Zhou's 7→8 step gives −0.31 dex, the same sign as Xu's control.
- Xu's own two rows already disagree in sign (control −0.71, SC-CO2 +0.09). The within-laboratory spread spans the between-laboratory conflict, so treatment and pressure interact with pH. Xu's CV runs up to 47% (0.2 dex) for minor compounds.
- **What would confirm:** same tests as cell 1.
- **Confidence: medium.**

**3. Benzaldehyde x sugar, xylose vs glucose (chen2023 −0.38 vs theng2024 +0.04 / +0.18 / +0.20).**
- **Explanation: neither laboratory's contrast is significant by its own test, and the analyte is barely a Maillard product in Theng.**
  - Theng: benzaldehyde is already 0.31-0.37 ppm in the unheated yeast extract and hydrolysate. The glucose values (0.32-0.35) and xylose values (0.38-0.51) all carry the same letter 'a', with SDs up to 25% (0.11 dex).
  - Chen: Table S1 also gives xylose 87.6 ± 14.7 and glucose 211.3 ± 21.2 the same letter 'a'.
  - Chen's identification is questionable. Table S1 lists hexanal at RI 872 and benzaldehyde at RI 873 on a 5MS column. By standard literature values (my own knowledge, not from these papers), hexanal is about 800 and benzaldehyde about 960 on that phase. Two compounds about 160 RI units apart cannot both elute near 872, so at least one assignment or RI is wrong.
  - Chen's glucose sample also has an outlying 2,4-dimethylbenzaldehyde (386 vs about 5-28 elsewhere), consistent with peak misassignment.
- **What would confirm:** Chen's raw MS spectrum and RI for the peak at 873. Refuted if a replicate xylose-vs-glucose run in the Lycium system reproduces −0.38 with an authentic-standard RI match.
- **Recommendation:** treat this cell as two noisy near-null contrasts, not a sign conflict. Down-weight Chen's benzaldehyde and hexanal identifications.
- **Confidence: medium-high** that this is not a real disagreement.

**4-6. Hexanal x sugar: xylose vs arabinose, vs fructose, vs glucose (chen2023 positive vs zhang2024b negative).**
- **Ruled out as the main cause: measurement precision and constant bias.** Repeatability is about 0.02-0.03 dex. Absolute between-laboratory scatter is 0.10-0.12 dex with an IS and up to 0.23 dex without (Aparicio, Casadei). A constant per-laboratory hexanal bias (recovery 91-150%) cancels in a ratio. None of this reaches the 0.5-1.3 dex gaps.
- **Unlikely: fibre competition in Zhang (FID, no IS).** Competition would need the xylose samples to carry more total volatile load. From Zhang's text (area and relative-area pairs), the xylose sample's total FID area is about 62.5e6/0.0291 ≈ 2.1e9, against about 585.8e6/0.0893 ≈ 6.6e9 for Fru-05. The xylose sample has less total load, not more. The competition sizes in olive oil (Aparicio, a few percent) and in Roberts below 10 ppm of competitor (≤ 0.1 dex) are also too small.
- **Supported: (i) a substrate and condition interaction.**
  - Zhang's hexanal is a pre-existing off-odourant of the raw meat hydrolysate. With pentose dose it falls to below detection above 2.5% xylose, and Zhang reads this as consumption in Maillard chemistry.
  - Chen's system differs on every axis that controls hexanal's fate: defatted plant meal, added cysteine, 120 C for 2 h, initial pH 7.5 unbuffered.
  - The pooled model treats sugar identity as transferable. Eckert's ring trial shows treatment effects need not transfer.
- **Supported: (ii) Chen's hexanal values are internally unstable.**
  - Its two methods disagree on the same contrast: xylose vs fructose is +0.46 by GC-MS but +0.07 by GC-IMS (a 0.39 dex gap); xylose vs glucose is +0.34 vs +0.46.
  - Galactose has no hexanal by GC-MS but one of the higher GC-IMS values (653).
  - The GC-MS RI (872) does not fit hexanal.
  - The GC-IMS sample is a freeze-dried product. Cordero shows that matrix physical state alone shifts single-shot headspace response by about 0.2 dex, and freeze-drying adds volatile losses that vary with the matrix.
- **What would confirm:**
  - Hexanal in Chen's unreacted hydrolysate, and the RI and spectrum check.
  - Zhang's heated sugar-free control against the pentose samples (consumption predicts a dose-dependent fall below the control).
  - A spiked-hexanal recovery in the xylose and fructose products of each laboratory. Equal recoveries within about 0.1 dex rule out a lever-dependent measurement bias.
- **Confidence: low-medium.** The evidence argues against a measurement artefact on Zhang's side and casts doubt on Chen's identification. It cannot settle whether the true chemistry differs between systems.

**7. TMP x amino acid, none vs arginine (hwang1995 vs park2022).**
- **Explanation: "none" means different things, and only Hwang neutralised arginine's basicity.**
  - Hwang's control is glucose plus glycine. Arginine is a second amino acid that competes with glycine; Hwang reports it lowered glycine's reactivity. The mixture was adjusted to pH 7 and heated as a 12-14% moisture solid at 180 C.
  - Park's control is mealworm protein plus xylose with no free amino acid. Adding 0.1 g arginine, a strongly basic amino acid, with no pH adjustment reported, adds free amino nitrogen and likely raises the pH.
  - Zhou 2023b shows pyrazine formation is strongly pH-dependent: none at pH 6-7 from the ARP alone, present at pH 8.
- **Method differences:** purge-and-trap with no replicate SD (Hwang) against 2 h SDE (Park). They are unlikely to explain a 1 dex sign difference.
- **Metadata:** the matrix labels are swapped in the cell (Part B).
- **What would confirm:** the pH of Park's mixture with and without arginine, or a pH-fixed rerun. Refuted if Park's arginine effect persists at fixed pH.
- **Confidence: medium** that this is a definitional or condition difference rather than measurement; **low** for the pH mechanism specifically, since Park's pH is not reported.

---

## Overall
1. **Measurement precision is much smaller than the effects.** Single-laboratory repeatability is 0.01-0.03 dex (Roberts, Aparicio, Casadei, Cordero). Harmonised-protocol between-laboratory scatter in absolute values is 0.10-0.24 dex (Aparicio, Casadei). Neither comes close to the 0.5-1.3 dex sign flips.
2. **Documented lever-dependent biases are about 0.1-0.45 dex:**
   - competition: ≤ 0.03 dex in real oils; 0.17-0.6 dex only at about 100 ppm of competitor (Roberts, read from the figure);
   - matrix physical state: about 0.2 dex (Cordero);
   - protein binding of MFT and BMFD: 0.34-0.45 dex (Zhu);
   - pH-dependent partition of weak acids: qualitative only (Roberts).
3. **The largest documented driver is how the lever is realised** (Eckert: up to about 2 dex with sign flips). In the cells it shows up as an unbuffered "pH" lever whose final pH does not differ (Zhou), and a "none" control that means glycine in one laboratory and no free amino acid in the other (TMP).
4. **Two cells (benzaldehyde, and partly hexanal) rest on contrasts that are non-significant or carry questionable identifications in the source** (Chen's RIs).
