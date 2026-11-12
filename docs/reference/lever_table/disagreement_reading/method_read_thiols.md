# Method read: why laboratories disagree (thiols, extraction, SPME) — 2026-11-12

Scope: seven assigned method papers + the methods sections (and the one or two table rows that set each conflict cell) of the conflict-cell sources.
Dex conversions: log10(1 + RSD) ~ RSD/2.3; recovery R -> log10(R).

Reading status. The task asked for the Read tool. I used it (page images + text) for hofmann2001b (all 5 pages), tominaga2005 (all 5 pages) and hwang1995 p. 3 (Table 1). For the other files I read a `pdftotext` dump in full (muller2007, sun2017, mateo-vivaracho2006, engel1999, chaintreau2001) and used `pdftotext -layout` to check the table rows that set each cell. So figures with no printed numbers were read from captions and text only. For the conflict-cell sources I read the methods sections and the relevant table rows, not the whole paper.

File note: `Zhou2023.pdf` on disk is NOT the conflict-cell paper. It is Zhu et al. on kafirin and baijiu (see `data/keys/pdf_dossier_aliases.yml`). The cell's zhou2023 is `Zhou2023b.pdf` (Zhou et al., JAFC 71:2472, 10.1021/acs.jafc.2c08360) plus `zhou2023b_supporting.pdf`. I read those.

---

## Part A. Assigned method papers

### A1. Hofmann, Czerny, Calligaris & Schieberle 2001, JAFC 49:2382 — coffee melanoidins bind thiols (read in full)
- **What and how.** A coffee-brew odorant recombinate in phosphate buffer at pH 5.4, with and without coffee melanoidins at 125 mg/10 mL, the same as in a brew. Static headspace GC at 30 °C. Carbonyls were quantified by an external-reference SIDA: the labelled standards sat in a separate buffer, because a standard added to the sample would be bound too. Thiols (FFT, MBT, MMBF) were compared only as headspace before vs after.
- **Numbers.**
  - Thiol loss is fast (Fig. 2): about 50 % of FFT is gone from the headspace after 20 min (-0.30 dex) and it is "nearly absent" at 30 min (≥ 1 dex), at 30 °C.
  - Low-MW melanoidins (fraction IV, 1500–3000 Da) bind most (Fig. 4). Even the least active fraction (II) cut MBT by about 75 % (-0.6 dex).
  - The binding is largely irreversible: dithioerythritol regenerated "only a certain portion" of FFT, and less of MBT and MMBF (Fig. 5).
  - Short carbonyls are NOT bound (Table 3, 1 h / 3 h, control vs melanoidin, µg/L air):
    - acetaldehyde 34.6/25.9 vs 35.8/27.4
    - methylpropanal 13.0/17.2 vs 13.5/15.7
    - 3-methylbutanal 13.8/12.9 vs 13.4/12.6
    - 2,3-butanedione 3.8/5.4 vs 3.8/5.7
    - 2,3-pentanedione 3.7/5.1 vs 3.7/5.5
    - All shifts are ≤ 0.04 dex.
  - Mottram 1996 is cited: heating disulfides with egg albumin gave a >100-fold (>2 dex) fall in disulfides, which were reduced to thiols. No loss occurred in water or maltodextrin.
- **Bears on.** MFT × pH and bis(2-methyl-3-furyl) disulfide (BMFD) × pH. It is direct evidence that a thiol's measured level depends on the browned matrix it sits in. Browning changes with the lever (pH), so this bias does NOT cancel in a ratio. It also gives quantitative evidence AGAINST melanoidin binding of short aldehydes, which bears on the hexanal and benzaldehyde cells by extrapolation only: hexanal and benzaldehyde were not tested.

### A2. Müller & Hofmann 2007, JAFC 55:4095 — FFT loss to phenol conjugates (read in full, text)
- **What and how.** Coffee brew held at 30 °C. FFT was measured by SIDA ([2H]2-FFT) after trapping with 4-vinylpyridine; phenols and conjugates by SIDA LC-MS/MS.
- **Numbers.**
  - Free FFT falls fast (Table 3): 0.17 → 0.10 → 0.07 → 0.02 µmol/L at 0, 10, 20 and 60 min. That is -0.23, -0.39 and -0.93 dex. FFT is undetectable at 90 min (Fig. 2).
  - The trap is hydroxyhydroquinone (238 → 8.2 µmol/L in 60 min). Conjugate 9 formed in 10 min accounts for about 60 % of the initial FFT.
  - The loss rate scales with the trap: conjugates rise up to 10-fold from light to medium roast (Table 4). With FFT added, conjugate 9 was 46 times higher at 20 min.
- **Bears on.** Absolute thiol levels: a lab-to-lab difference of about 1 dex is reachable just from the delay between the end of the reaction and extraction (minutes to an hour). A ratio cancels this only if the delay and the trap concentration are the same in both arms. The trap (oxidisable phenols, browning) is lever-dependent. Coffee-specific; no direct bearing on any single cell.

### A3. Sun et al. 2017, Food Anal. Methods 10.1007/s12161-017-1034-8 — free vs total FFT (read in full, text)
- **What and how.** HS-SPME (DVB/PDMS, 40 °C) on coffee brew. External calibration in a "prefabricated" coffee matrix, stripped of FFT with cysteine and vacuum distillation. 1,2-dichlorobenzene injection standard. Total FFT was released with 16 g/L cysteine.
- **Numbers.**
  - **Free vs total.** Free FFT is 11.3, 15.3 and 21.0 µg/L; total (cysteine-released) is 130.3, 114.9 and 138.2 µg/L. Free/total = 0.09–0.15, i.e. **-0.8 to -1.0 dex**.
  - **Matrix calibration.** The buffer calibration slope (0.017) is 48.6 % of the coffee-matrix slope (0.035), a **0.31 dex** matrix effect. The spike-in-brew slope is 89.6 % of the matrix slope (-0.05 dex).
  - **Recovery and precision.**
    - Recovery 86.8–106.2 % (spread 0.09 dex) (Table 3).
    - CV 7.1 %, n = 5 (0.03 dex) (Table 2).
    - Intra-day RSD 2.0–4.9 %; inter-day RSD 6.5–9.7 % (0.03–0.04 dex).
  - Cysteine release plateaued at 4–8 g/L cysteine (Fig. 2).
- **Bears on.** The ~1 dex between-lab gap in absolute thiol levels. One method's "free headspace thiol" and another's "total" (cysteine, SIDA with label equilibrated in the matrix, derivatisation) differ by about 0.8–1.0 dex in coffee. The paper notes that a labelled standard spiked into the matrix is itself partly bound, so SIDA reports total, not free. Buffer vs matrix calibration adds about 0.3 dex. Both biases scale with the binding matrix, so neither is guaranteed to cancel in a ratio across a browning lever.

### A4. Tominaga & Dubourdieu 2006 (file tominaga2005), JAFC 54:29 — MFT and FFT in wine (read in full)
- **What and how.** pHMB conjugation in the wine at pH 7, then a Dowex column, cysteamine reverse elution, and DCM. GC-MS SIM. Internal standard 4-methoxy-2-methyl-2-mercaptobutane (not isotopic).
- **Numbers.**
  - Repeatability CV: MFT 3.0 % (0.013 dex), FFT 4.2 % (0.018 dex) (Table 2).
  - **Matrix-independent spike recovery** across 2 white and 2 red wines (Table 3). MFT +25 ng/L recovered 23–26; +100 recovered 104–110 (spread ≤ 0.03 dex). FFT 23–28 and 98–108.
  - Overall recovery is 75–85 % for the thiols and > 85 % for the internal standard. Forward elution gave "very poor and not reproducible" MFT recovery (Fig. 4).
  - MFT in a DCM extract at -10 °C declined "clearly" after 24 h. Fig. 5 shows roughly 10–40 % loss by 48–72 h (-0.05 to -0.2 dex). FFT was stable for 48 h.
  - The calibration R² for MFT is only 0.967.
- **Bears on.** MFT × pH and BMFD × pH, as a contrast case. A thiol-selective, protected extraction can make thiol recovery matrix-independent to ≤ 0.03 dex; headspace methods do not demonstrate this. MFT is the least stable analyte, so the timing from extraction to injection matters (about 0.05–0.2 dex per day even when cold).

### A5. Mateo-Vivaracho, Ferreira & Cacho 2006, J. Chromatogr. A 1121:1 — on-fibre PFBBr derivatisation HS-SPME (read in full, text)
- **What and how.** PDMS–DVB fibre loaded with PFBBr and tributylamine, then sample headspace at 55 °C, then GC-NCI-MS, on wine.
- **Numbers.**
  - **MFT could not be quantified:** "the signal was independent of the concentration".
  - **Oxidation in the autosampler queue:** over 14 sequential vials at room temperature, the last sample's signal was 20–50 % lower than the first (-0.10 to -0.30 dex). It was fixed only by N2 bubbling plus 4 °C storage; BHA had no effect and EDTA a weak one.
  - **Matrix effect on slope** (Table 4): FFT synthetic wine 492 vs real wines 711–720 (0.16 dex). 3-Mercaptohexyl acetate (MHA) 96 vs 122–132 (0.10–0.14 dex).
  - Repeatability RSD 10–21 % (0.04–0.08 dex) (Table 3). Spike recovery within ±20 % across 5 wines (Table 5).
  - Extraction did not reach equilibrium and rose with time and temperature (Fig. 3).
- **Bears on.** MFT × pH and BMFD × pH (mechanism, not those labs' data):
  - Sequential-queue oxidation of 0.1–0.3 dex would be lever-aligned if the samples were run in lever order.
  - Matrix-dependent SPME slopes of 0.1–0.16 dex give the order of magnitude for a fibre/matrix bias that moves with matrix composition.
  - Non-equilibrium SPME makes the response sensitive to small differences in time and temperature.

### A6. Engel, Bahr & Schieberle 1999, Eur Food Res Technol 209:237 — SAFE (read in full, text)
- **What and how.** Yields of SAFE vs high-vacuum transfer (HVT) from a solvent, and from a 50 % fat matrix.
- **Numbers.**
  - Solvent (Table 2): SAFE recovers 100 % of the seven odorants.
  - **50 % fat (Table 3), SAFE:**
    - 3-methylbutanoic acid 37 % (-0.43 dex)
    - phenylacetaldehyde 26 % (-0.59 dex)
    - sotolon 4.3 %
    - 2-phenylethanol 12.8 %
    - (E,E)-2,4-decadienal 4.3 % (-1.37 dex)
    - damascenone 4.2 %
    - vanillin 0.5 %
  - HVT is lower still. SD < 3–5 %.
  - n-Alkanes: HVT yield collapses from C14 (48 %); SAFE from C16 (59 %) (Table 1).
- **Bears on.** Fat content changes recovery by 0.4–2 dex, and compound-specifically. This cancels in a within-study ratio only if the lever leaves the lipid phase unchanged. For sugar-identity contrasts it does (same hydrolysate), so: **no direct bearing** on the seven cells. It does explain why absolute aldehyde levels differ between a defatted seed meal (Chen) and a chicken hydrolysate (Zhang).

### A7. Chaintreau 2001, Flavour Fragr. J. 16:136 — SDE review (read in full, text; figure values read from captions and text)
- **What and how.** A review of atmospheric SDE (A-SDE) and vacuum SDE (V-SDE): recovery, artefacts and theory.
- **Numbers.**
  - Godefroot micro-SDE: mean RSD per peak 5.8 % (0.025 dex).
  - Cinnamon (Table 4): benzaldehyde RSD 17.6 % in one sample vs 1.7 % in another (0.07 vs 0.007 dex). Precision is matrix-dependent.
  - Filbertone CV: SDE 4.7 %, SFE 25.5 %, LC/GC 9.1 %.
  - Recovery peaks at 30–45 min and then declines (Fig. 6, benzaldehyde among the analytes). Too much N2 flow cuts recovery (Fig. 7, benzaldehyde).
  - Lipid matrices "greatly decrease" recovery and need longer distillation.
  - **Artefacts:**
    - Oxidation under air (cod lipids → 4(Z)-heptenal; 2,4-decadienals in nut oil).
    - Ester hydrolysis (linalyl acetate 50 % in A-SDE).
    - **Maillard reaction during A-SDE.** A glucose–leucine mixture reacts in A-SDE but not in static V-SDE. Furfural appears in honey only by A-SDE. 2-Methylpyrazine was *synthesised* in an A-SDE at the Shibamoto yield.
    - Mottram's protein–disulfide case: in A-SDE, disulfides came out as thiols.
- **Bears on.** Trimethylpyrazine (TMP) × amino acid, directly. Park 2022 used A-SDE (Likens–Nickerson, 2 h) on the reacted mixture, which still contains residual sugar and the added free amino acid. Pyrazine formation during extraction is then an artefact that scales with the lever (free amino acid). It also bears on oxidation-sensitive aldehydes (benzaldehyde, hexanal) for any lab that used SDE; none of the hexanal or benzaldehyde cells did (all used HS-SPME or GC-IMS).

Between-laboratory SDs. **None of the seven assigned papers reports an interlaboratory study.** The closest figure is Sun's inter-day RSD of 6.5–9.7 %, which is within-lab.

---

## Part B. Conflict-cell sources: methods and the rows that set each cell (methods read; rows checked in the tables)

Sign convention check (TMP cell). As I read the tables:
- Hwang 1995 Table 1: TMP is 133.57 (glycine only) vs 73.43 µg/g glucose (glycine + Arg), i.e. none→Arg = **-0.26**.
- Park 2022 Table 2: TMP is 0.042 (MP-RF) vs 0.250 relative peak area (MPArg-RF), i.e. none→Arg = **+0.77**.

The cell records +0.26 and -0.77. Both signs are flipped relative to my reading, so the lever table appears to use log10(none/Arg). The sign disagreement stands either way, but please confirm the convention. The recorded matrix labels also look swapped. Hwang is a glucose/glycine/Arg + wheat-starch model, dried to 12–14 % moisture and heated at 180 °C. Park is mealworm protein + xylose in water at 125 °C, with no pH adjustment.

MFT and BMFD rows:
- **xu2007** (Table 1, ng/mL headspace-equivalent).
  - MFT: control 16.62 / 2.63 / ND (<0.07) at initial pH 5.6 / 7.0 / 8.0; SC-CO2 34.42 / 25.84 / 7.06.
  - BMFD: control 3.94 / 1.23 / ND; SC-CO2 7.87 / 10.39 / 1.4.
  - **Final pH:** control 5.5 / 6.2 / 6.8; SC-CO2 5.5 / 6.0 / 6.3.
  - The two Xu rows per cell are evidently the control and SC-CO2 arms.
- **zhou2023 (Zhou2023b)** (Table 1, µg/L, ARP + Cys). MFT 697 / 1589 / 526 at pH 6 / 7 / 8, i.e. +0.36 for 6→7 and -0.48 for 7→8. BMFD 59.7 / 102.6 / 50.1. Non-monotonic, with the optimum at pH 7.

Methods:
- **xu2007** (J Sci Food Agric 88:328). Ribose + Cys, 0.05 M each, in 0.2 M phosphate, 140 °C for 1 h, with and without SC-CO2. HS-SPME on DVB/CAR/PDMS, 10 min equilibration + 20 min extraction at 60 °C, of 1 mL undiluted reaction mixture (plus the absorption buffers, summed). Tridecane internal standard with **response factor 1, TIC areas**. The authors call the results "relative ... not absolute concentrations". CV < 20 %, up to 47 %. Samples were stored at ≤ -18 °C "until use". A420 (browning) depends on pH and arm. At pH 7, SC-CO2 browned less than the control.
- **zhou2023b** (JAFC 71:2472). Purified Xyl–Ala ARP 20 mM + Cys 20 mM, pH adjusted with NaOH (no buffer), 120 °C for 1 h. HS-SPME on CAR/PDMS/DVB, 3 mL sample + 2 mL saturated NaCl, 60 °C for 30 min. 1,2-dichlorobenzene internal standard. **External calibration curves prepared "in water"** (SI Table S1; MFT y = 0.001x + 0.0029). This is not isotope dilution and not matrix-matched, so the lever table's "absolute" overstates it: it is calibrated, but it does not correct binding or matrix partition. Final pH dropped (Fig. 2; values only plotted). The ARP route is described as a way to make "light-colored" products.
- **chen2023** (Foods 12:4346).
  - Reaction: Lycium seed meal (oil removed by SC-CO2) hydrolysate 1.0 g + sugar 0.3 g + Cys 0.15 g at 10 %, pH 7.5, 120 °C for 2 h in a 20 mL ampoule.
  - **The supernatant was freeze-dried before the volatile analysis.**
  - GC-MS: HS-SPME CAR/PDMS 75 µm, 50 °C for 30 min, o-dichlorobenzene internal standard, ng/g.
  - GC-IMS: 0.05 g powder in 2 mL water, 60 °C for 15 min.
  - Rows (SI Table S1, GC-MS): hexanal X 549.5, F 190.5, A 257.7, G 251.6; benzaldehyde X 87.6, G 211.3.
  - Rows (SI Table S2, GC-IMS): hexanal X 720.2, F 607.5, A 562.4, G 246.9.
  - **Chen's two recorded rows per hexanal cell are GC-MS and GC-IMS on the same samples.**
- **zhang2024b = zhang2024sugars** (Foods 13:991). Chicken carcass hydrolysate (1.25 % fat on a dry basis), natural pH 6.62, unbuffered, 100 °C for 1 h in 100 mL bottles holding 40 mL (60 mL air headspace). HS-SPME CAR/PDMS (85 µm), 5 mL sample, 60 °C for 30 min. **FID peak area with no internal standard; relative peak area (RPA) %.** Hexanal is already present before heating. It falls with pentoses and is "below the detected level" above 2.5 % xylose (censored). It rises with hexoses.
- **theng2024** (IJFST 59:6325). Yeast hydrolysate, pH 6.5, 100 °C for 2 h, no Cys. HS-SPME DVB/CAR/PDMS, 1 g sample + 1 g 30 % NaCl, 60 °C for 30 min. 2-methyl-3-heptanone internal standard, response factor 1. In this study benzaldehyde is the one aldehyde that does not fall with added sugar.
- **hwang1995**. Glucose / [15N]glycine / test amino acid at equal moles (2.66 µmol each, per the text), plus wheat starch, pH 7, freeze-dried to 12–14 % moisture, 180 °C for 1 h. Purge-and-trap at 80 °C for 1 h onto Tenax/Carbotrap. Deuterated toluene spiked onto the tube. Yields in µg/g glucose.
- **park2022**. Mealworm protein 2 g + xylose 0.1 g + amino acid 0.1 g + 2 mL water, 125 °C for 60 min. **Likens–Nickerson SDE, 2 h, 5 g in 200 mL water.** Benzyl acetate internal standard in the solvent flask. Relative peak area. The text ranks pyrazine totals Thr > Arg > Ser > ... > Cys > MP-RF: **every amino-acid arm exceeds the no-amino-acid base.**

---

## Part C. The seven cells: best-supported explanation and confidence

**1. MFT × pH (Xu − vs Zhou +)**
- **Best-supported: lever-dependent loss in the browned matrix, plus a different chemistry and pH range.**
  - **Browning-linked binding.** Within Xu, the sign depends on the arm. Control: 16.62 → 2.63, -0.80 dex for 5.6 → 7.0, and more browning. SC-CO2: 34.42 → 25.84, -0.12 dex, and less browning at pH 7. Hofmann 2001 shows browned-matrix thiol binding of ≥ 0.3–1 dex within 20–30 min at 30 °C, irreversible. Xu's 30 min of headspace sampling at 60 °C on the undiluted browned mixture gives that loss room to act. Browning rises with pH (Xu Fig. 2), so the loss is lever-aligned, in the negative direction. Zhou's ARP system is described as light-coloured.
  - **The pH axis is not the same.** Xu's final pH spans only 5.5 → 6.2 (control) and 5.5 → 6.0 (SC-CO2), and the pH 8 control point is censored at ND. Zhou is non-monotonic with the optimum at pH 7 (+0.36 then -0.48). Xu's 5.6 starting point lies below Zhou's range, and the precursors differ (free ribose vs pre-formed ARP).
  - **A further possible headspace bias (my inference; not in the papers read).** Heteroaromatic thiols ionise near neutral pH, so headspace SPME at the final pH would under-read MFT more at higher pH. Neither lab corrects for this. It would bias both labs negative, so it can inflate Xu's slope but cannot explain Zhou's positive one.
- **To confirm:** spike MFT (ideally labelled) into Xu-type reacted mixtures at each pH after the reaction and run Xu's HS-SPME. A recovery that falls with pH and browning by ≥ 0.3 dex confirms the bias. A recovery flat within ~0.1 dex refutes it. Alternatively, re-measure by a thiol-selective method (pHMB or SIDA), whose matrix independence is ≤ 0.03 dex (Tominaga Table 3).
- **Confidence:** medium that measurement-side binding is part of it; low on how it splits against the real chemistry.

**2. BMFD × pH (Xu − vs Zhou +)**
- **Best-supported: the same as cell 1.** Xu control 3.94 → 1.23 (-0.51 dex) vs Xu SC-CO2 7.87 → 10.39 (+0.12). The less-browned Xu arm already agrees in sign with Zhou (+0.24). Disulfide–protein/matrix interchange (Mottram, cited in Hofmann 2001 and Chaintreau 2001: >2 dex losses with protein) is a matrix-dependent sink. The cell rests in effect on Xu's control arm alone, and its pH 8 point is ND.
- **To confirm or refute:** as cell 1, with BMFD spiked after the reaction.
- **Confidence:** medium.

**3. Benzaldehyde × sugar identity (Chen − vs Theng +)**
- **Best-supported: a substantive matrix and chemistry difference, with small effects close to noise.** Chen has cysteine, pH 7.5, 120 °C, and a freeze-dried, reconstituted product. Theng has no cysteine, pH 6.5, 100 °C, and liquid samples. Theng's effects (+0.04 to +0.20) sit within a few within-lab SDs: SPME CVs of 10–20 % are 0.04–0.08 dex (Mateo-Vivaracho, Xu).
- **What the papers rule against.** Melanoidin binding of short aldehydes is ≤ 0.04 dex (Hofmann Table 3). Benzaldehyde oxidation and SDE time-dependence (Chaintreau) are not in play here, since neither lab used SDE.
- **Freeze-drying.** Chen's freeze-drying step is a plausible lever-dependent loss (product solids differ by sugar), but none of my papers quantifies it.
- **To confirm:** split each Chen arm before and after freeze-drying and compare; a benzaldehyde loss that differs by sugar by ≥ 0.2 dex supports it.
- **Confidence:** low.

**4–6. Hexanal × sugar identity (X vs A, X vs F, X vs G; Chen + vs Zhang −)**
- **Best-supported: different hexanal histories, plus method-dependent lever bias.**
  - **Direct evidence of method-dependent lever bias.** Chen's own two methods on the same samples disagree on the lever effect: X/A +0.33 (GC-MS) vs +0.11 (GC-IMS); X/F +0.46 vs +0.07; X/G +0.34 vs +0.46. That is 0.12–0.39 dex within one lab. This matches the pooled extra disagreement for mixed methods (~0.17 dex) and is of the order of the between-lab SD (0.45 for sugar identity).
  - **Zhang's quantitation.** Zhang has no internal standard (raw FID area, CAR/PDMS). Its xylose arm is censored (ND above 2.5 % xylose), which exaggerates the negative effect.
  - **Different hexanal histories (substantive).** In Zhang, hexanal pre-exists in a lipid-bearing (1.25 % fat) hydrolysate heated with 60 mL of air. It falls with pentoses, consistent with consumption or antioxidant MRPs in the faster-reacting pentose arms. Chen's matrix is defatted.
  - **Fibre competition.** It is plausible on Carboxen fibres, since pentose arms produce far more volatiles, but **none of my papers quantifies fibre competition**. The only SPME matrix-slope effects I have are 0.10–0.31 dex (Sun; Mateo-Vivaracho). That is smaller than Zhang's -0.45 to -0.84, so a bias alone is unlikely to flip the sign.
- **To confirm:** spike hexanal after the reaction into each Zhang sugar arm and re-run Zhang's HS-SPME-FID. A recovery that falls with xylose by ≥ 0.3 dex points to a measurement bias. A recovery flat within 0.1 dex means real consumption. For Chen, compare before vs after freeze-drying.
- **Confidence:** low–medium that the sign difference is mostly real chemistry (different hexanal sources); medium that method choice contributes 0.1–0.4 dex.

**7. TMP × amino acid identity (none vs Arg; Hwang vs Park)**
- **Best-supported: the two "none" references are not the same contrast, plus an A-SDE artefact in Park.**
  - **Reference mismatch.** Hwang's "none" still contains equimolar glycine, and Arg there competes with it (Hwang's abstract: Arg gives the lowest pyrazine yield, and Arg depresses glycine's reactivity). Park's "none" has no free amino acid: protein + xylose only. Every one of Park's 11 amino-acid arms raises pyrazines over the base, i.e. amino nitrogen is limiting there.
  - **Lever-aligned A-SDE artefact.** Park extracted by 2 h of atmospheric SDE on a mixture still holding residual xylose and the added free Arg. Chaintreau documents that A-SDE runs the Maillard reaction and can make methylpyrazine at preparative yield. Extra free amino acid then means extra pyrazine formed during extraction.
- **To confirm:**
  - Re-extract Park's arms by SAFE, V-SDE or headspace at ≤ 40 °C. If TMP(Arg)/TMP(none) falls well below 6× (0.77 dex), the artefact is confirmed.
  - SDE an unheated Park Arg mixture; any TMP found is pure artefact.
  - For the reference mismatch: Park's TMP for Gly vs Arg arms is a matched contrast to Hwang's.
- **Confidence:** medium–high that the cell is a definition mismatch rather than a lab disagreement; medium for the SDE artefact.

**Cross-cutting: absolute thiol levels (~1 dex between labs).**
- **Best-supported: a free vs total difference plus calibration matrix.** Free vs total FFT differs by 0.8–1.0 dex (Sun). Buffer vs matrix calibration adds 0.31 dex (Sun). The delay to extraction costs 0.2–0.9 dex in 10–60 min (Müller; Hofmann). Zhou's calibration in water and Xu's response factor of 1 are both free-headspace measures with different response factors.
- **Confidence:** high that these mechanisms can produce about 1 dex. They are lever-dependent whenever the lever changes browning or phenol/melanoidin load.
